#!/usr/bin/env julia
# SeedRef 2.1 - seed reads and a reference plastome for every read library.
#
# Dependencies: CodecZlib (already required by the assembler) and the Julia
# standard library only (Downloads, Printf, Dates, SHA). No BLAST+, HTTP.jl or
# ProgressMeter.
#
# Folder mode (prototype for benchmarking; files are named after the library):
#   julia --project=. --threads auto SeedRef.jl --input ./input --output ./seedref \
#         --config seedref_config.txt --template config_v16.txt [--results ./results]
# With --template, a ready assembler config is written per library:
#   seedref/<library>_config.txt  ->  julia --project=. --threads 64 AssemblerComplete_v16.jl seedref/<library>_config.txt
#
# Paired-end libraries only.
# writes <output>/<library>_seed.fasta (single best read), <library>_seeds_all.fasta,
# <library>_reference.gb, <library>_reference.fasta and seedref_summary.tsv.
#
# Per library:
#   1. Marker k-mers come from an automatic panel of RefSeq plastomes spanning land
#      plants (liverworts, mosses, hornworts, lycophytes, gymnosperms, angiosperms),
#      downloaded once and cached, plus optional user FASTA files (<gene>.fasta).
#   2. One streaming pass over the first N reads of R1 finds the best reads per
#      marker gene (deterministic tie-break by read order).
#   3. Reference: one NCBI BLAST URL API request with those reads, restricted to
#      RefSeq complete plastomes; highest summed best bit score wins. A local
#      <library>.gb beside the reads, or Reference GenBank, avoids the network.
#   4. Seed refinement: the same marker genes are cut from the chosen reference and
#      the reads are scanned again, so the final seed matches the closest relative.
#   Genome Range = reference length +/- margin.
#
# Single-sample mode (reads an assembler config): julia SeedRef.jl config_v16.txt

module SeedRef

using Printf
using Dates
using Downloads
using CodecZlib
using SHA

const VERSION_STRING = "2.3.3"

# ---------------------------------------------------------------------------
# Settings and shared key = value parsing
# ---------------------------------------------------------------------------
Base.@kwdef mutable struct SeedSettings
    marker_dir::String = ""                              # Optional user FASTA: <gene>.fasta
    markers::Vector{String} = ["rbcL", "psaA", "psbB", "atpB"]  # LSC genes; priority order.
    panel::Bool = true                                   # Automatic RefSeq marker panel.
    refine::Bool = true                                  # Second scan with the reference's genes.
    kmer::Int = 25
    max_reads::Int = 2_000_000                           # First N records of R1.
    min_hits::Int = 3                                    # K-mer positions per read.
    reads_per_marker::Int = 2                            # BLAST queries per marker.
    batch_size::Int = 20_000
    reference_mode::String = "auto"                      # auto | local | none
    reference_genbank::String = ""                       # local: one GenBank for all samples
    genome_margin::Int = 10_000
    blast_database::String = "nt"
    blast_entrez_query::String = "srcdb_refseq[prop] AND (chloroplast[Title] OR plastid[Title]) AND complete[Title] AND genome[Title]"
    blast_hitlist::Int = 100
    blast_timeout::Float64 = 3600.0
    tool::String = "AssemblerCompleteV16"
    reference_cache::String = ""                         # Shared GenBank cache (batch).
end

"Ordered `key => value` pairs; `#` starts a comment, duplicate keys are errors."
function read_key_values(path::AbstractString)
    pairs = Pair{String,String}[]
    seen = Set{String}()
    for (lineno, raw) in enumerate(eachline(path))
        line = strip(first(split(raw, '#'; limit=2)))
        isempty(line) && continue
        occursin('=', line) || error("$path line $lineno: expected key = value")
        key, value = strip.(split(line, '='; limit=2))
        key in seen && error("$path: duplicate key '$key'")
        push!(seen, key)
        push!(pairs, String(key) => String(value))
    end
    return pairs
end

resolve_path(base::AbstractString, value::AbstractString) =
    isempty(value) ? "" : normpath(joinpath(base, expanduser(value)))

const SEED_KEYS = Set(["Marker folder", "Seed markers", "K-mer size", "SeedIR reads",
    "Marker panel", "Refine seeds with reference",
    "Min k-mer hits for seed", "Seed reads per marker", "Reference mode", "Reference GenBank",
    "Genome range margin", "BLAST database", "BLAST Entrez query", "BLAST hitlist size",
    "BLAST timeout"])

"Build SeedSettings from config pairs; returns (settings, keys consumed)."
function seed_settings(pairs, base::AbstractString)
    s = SeedSettings()
    used = Set{String}()
    for (key, value) in pairs
        key in SEED_KEYS || continue
        push!(used, key)
        if key == "Marker folder"
            s.marker_dir = resolve_path(base, value)
        elseif key == "Marker panel"
            lowercase(value) in ("auto", "none") || error("Marker panel must be auto or none")
            s.panel = lowercase(value) == "auto"
        elseif key == "Refine seeds with reference"
            s.refine = parse(Bool, lowercase(value))
        elseif key == "Seed markers"
            s.markers = String.(filter(!isempty, strip.(split(value, ','))))
        elseif key == "K-mer size"
            s.kmer = parse(Int, value)
        elseif key == "SeedIR reads"
            s.max_reads = parse(Int, value)
        elseif key == "Min k-mer hits for seed"
            s.min_hits = parse(Int, value)
        elseif key == "Seed reads per marker"
            s.reads_per_marker = parse(Int, value)
        elseif key == "Reference mode"
            s.reference_mode = lowercase(value)
        elseif key == "Reference GenBank"
            s.reference_genbank = resolve_path(base, value)
        elseif key == "Genome range margin"
            s.genome_margin = parse(Int, value)
        elseif key == "BLAST database"
            s.blast_database = value
        elseif key == "BLAST Entrez query"
            s.blast_entrez_query = value
        elseif key == "BLAST hitlist size"
            s.blast_hitlist = parse(Int, value)
        elseif key == "BLAST timeout"
            s.blast_timeout = parse(Float64, value)
        end
    end
    if isempty(s.marker_dir) && isdir(joinpath(base, "markers"))
        s.marker_dir = joinpath(base, "markers")
    end
    validate_settings(s)
    return s, used
end

function validate_settings(s::SeedSettings)
    1 <= s.kmer <= 32 || error("K-mer size must be in 1:32")
    s.max_reads > 0 || error("SeedIR reads must be positive")
    s.min_hits > 0 || error("Min k-mer hits for seed must be positive")
    s.reads_per_marker > 0 || error("Seed reads per marker must be positive")
    s.reference_mode in ("auto", "local", "none") || error("Reference mode must be auto, local or none")
    s.reference_mode == "local" && isempty(s.reference_genbank) &&
        error("Reference mode = local requires Reference GenBank (or use auto with <sample>.gb beside the reads)")
    s.genome_margin >= 0 || error("Genome range margin cannot be negative")
    isempty(s.markers) && error("Seed markers is empty")
    s.panel || !isempty(s.marker_dir) || error("Marker panel = none requires a Marker folder")
    return s
end

# ---------------------------------------------------------------------------
# Canonical 2-bit k-mers (identical encoding convention to the assembler)
# ---------------------------------------------------------------------------
const ENC = let t = fill(0xff, 256)
    for (i, b) in enumerate(codeunits("ACGT"))
        t[Int(b)+1] = UInt8(i-1)
        t[Int(b)+0x21] = UInt8(i-1)          # lowercase
    end
    t
end

"Fill `out` with canonical k-mers of every ACGT-only window (windows with other bytes are skipped)."
function canonical_kmers!(out::Vector{UInt64}, seq::AbstractVector{UInt8}, k::Int)
    empty!(out)
    mask = k == 32 ? typemax(UInt64) : (UInt64(1) << (2k)) - UInt64(1)
    shift = 2*(k-1)
    f = UInt64(0); r = UInt64(0); run = 0
    @inbounds for b in seq
        x = ENC[Int(b)+1]
        if x == 0xff
            f = UInt64(0); r = UInt64(0); run = 0
            continue
        end
        f = ((f << 2) | UInt64(x)) & mask
        r = (r >> 2) | (UInt64(x ⊻ 0x03) << shift)
        run += 1
        run >= k && push!(out, min(f, r))
    end
    return out
end

# ---------------------------------------------------------------------------
# Marker k-mer database. A k-mer keeps EVERY (marker, order) group containing it,
# so conserved k-mers no longer go to whichever order happened to be read last.
# ---------------------------------------------------------------------------
struct MarkerDB
    k::Int
    genes::Vector{String}
    orders::Vector{String}
    group_gene::Vector{Int}          # group -> gene index
    group_order::Vector{Int}         # group -> order index
    group_total::Vector{Int}         # distinct canonical k-mers per group
    index::Dict{UInt64,Int32}        # canonical k-mer -> entry
    gene_mask::Vector{UInt32}        # entry -> bitmask of genes
    group_offsets::Vector{Int32}     # CSR entry -> groups
    group_ids::Vector{Int32}
end

"FASTA records as (first header token, uppercase bytes)."
function read_fasta_records(path::AbstractString)
    records = Tuple{String,Vector{UInt8}}[]
    name = ""; buf = IOBuffer()
    flush_record() = (isempty(name) || push!(records, (name, take!(buf))))
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, ">")
            flush_record()
            tokens = split(line[2:end])
            isempty(tokens) && error("$path: empty FASTA header")
            name = String(tokens[1])
        else
            isempty(name) && error("$path: sequence before the first header")
            write(buf, uppercase(line))
        end
    end
    flush_record()
    return records
end

"(gene, group, sequence) marker records; group = order (user FASTA) or organism (panel)."
const MarkerRecord = Tuple{String,String,Vector{UInt8}}

function find_marker_file(dir::AbstractString, gene::AbstractString)
    isempty(dir) && return nothing
    for ext in (".fasta", ".fa", ".fas", ".fna")
        path = joinpath(dir, gene*ext)
        isfile(path) && return path
    end
    return nothing
end

"User marker FASTA files for the configured genes, plus the automatic panel."
function marker_records(s::SeedSettings; io::IO=stdout)
    records = MarkerRecord[]
    for gene in s.markers
        path = find_marker_file(s.marker_dir, gene)
        path === nothing && continue
        for (group, seq) in read_fasta_records(path)
            push!(records, (gene, group, seq))
        end
    end
    s.panel && append!(records, panel_records(s; io=io))
    isempty(records) && error("No marker sequences: panel download failed and no <gene>.fasta in Marker folder")
    return records
end

"K-mer database over marker records; genes keep the priority order of `gene_order`."
function build_marker_db(records::Vector{MarkerRecord}, k::Int, gene_order::Vector{String})
    present = Set(first.(records))
    genes = [g for g in gene_order if g in present]
    append!(genes, sort!(collect(setdiff(present, Set(genes)))))
    length(genes) <= 32 || error("At most 32 marker genes are supported")
    gene_index = Dict(g => i for (i, g) in enumerate(genes))
    orders = String[]; order_index = Dict{String,Int}()
    group_of = Dict{Tuple{Int,Int},Int}()
    group_gene = Int[]; group_order = Int[]
    members = Dict{UInt64,Vector{Int32}}()
    kmers = UInt64[]
    for (gene, group, seq) in records
        g = gene_index[gene]
        o = get!(order_index, group) do
            push!(orders, group); length(orders)
        end
        grp = get!(group_of, (g, o)) do
            push!(group_gene, g); push!(group_order, o); length(group_gene)
        end
        for key in canonical_kmers!(kmers, seq, k)
            list = get!(Vector{Int32}, members, key)
            Int32(grp) in list || push!(list, Int32(grp))
        end
    end
    isempty(members) && error("Marker sequences contain no valid $(k)-mers")
    keys_sorted = sort!(collect(keys(members)))          # Deterministic entry numbering.
    index = Dict{UInt64,Int32}(); sizehint!(index, length(keys_sorted))
    gene_mask = zeros(UInt32, length(keys_sorted))
    offsets = Int32[0]; ids = Int32[]
    group_total = zeros(Int, length(group_gene))
    for (e, key) in enumerate(keys_sorted)
        index[key] = Int32(e)
        for grp in sort!(members[key])
            gene_mask[e] |= UInt32(1) << (group_gene[grp]-1)
            group_total[grp] += 1
            push!(ids, grp)
        end
        push!(offsets, Int32(length(ids)))
    end
    return MarkerDB(k, genes, orders, group_gene, group_order, group_total,
                    index, gene_mask, offsets, ids)
end

build_marker_db(s::SeedSettings; io::IO=stdout) = build_marker_db(marker_records(s; io=io), s.kmer, s.markers)

# ---------------------------------------------------------------------------
# Streaming read scan
# ---------------------------------------------------------------------------
mutable struct TopReads                  # Best-first; ties -> earlier read wins.
    score::Vector{Int}
    index::Vector{Int}
    id::Vector{String}
    seq::Vector{Vector{UInt8}}
end
TopReads() = TopReads(Int[], Int[], String[], Vector{UInt8}[])

@inline better(s1, i1, s2, i2) = s1 > s2 || (s1 == s2 && i1 < i2)

function offer!(t::TopReads, cap::Int, score::Int, index::Int, id::String, seq::Vector{UInt8})
    score > 0 || return t
    n = length(t.score)
    n == cap && !better(score, index, t.score[end], t.index[end]) && return t
    pos = n+1
    while pos > 1 && better(score, index, t.score[pos-1], t.index[pos-1])
        pos -= 1
    end
    insert!(t.score, pos, score); insert!(t.index, pos, index)
    insert!(t.id, pos, id); insert!(t.seq, pos, seq)
    if length(t.score) > cap
        pop!(t.score); pop!(t.index); pop!(t.id); pop!(t.seq)
    end
    return t
end

struct ScanState                         # Owned by one logical chunk, never by threadid().
    matched::BitVector
    top::Vector{TopReads}
    kmers::Vector{UInt64}
    scores::Vector{Int}
end
ScanState(db::MarkerDB) = ScanState(falses(length(db.gene_mask)), [TopReads() for _ in db.genes],
                                    UInt64[], zeros(Int, length(db.genes)))

function open_reads(path::AbstractString)
    raw = open(path, "r")
    endswith(lowercase(path), ".gz") || return raw
    try
        return GzipDecompressorStream(raw; bufsize=1 << 20)
    catch
        close(raw); rethrow()
    end
end

function next_fastq!(io::IO, path::AbstractString, record::Int)
    eof(io) && return nothing
    header = readline(io)
    isempty(header) && eof(io) && return nothing
    seq = readline(io); plus = readline(io); qual = readline(io)
    startswith(header, "@") || error("$path: invalid FASTQ header at record $record")
    startswith(plus, "+") || error("$path: missing '+' line at record $record")
    ncodeunits(seq) == ncodeunits(qual) || error("$path: sequence/quality length mismatch at record $record")
    tokens = split(header[2:end])
    id = isempty(tokens) ? "read$record" : String(tokens[1])
    return id, Vector{UInt8}(codeunits(seq))
end

struct ScanResult
    reads::Int
    top::Vector{TopReads}                # Per marker (db.genes order).
    group_matched::Vector{Int}           # Distinct matched k-mers per (marker, order) group.
end

function scan_reads(path::AbstractString, db::MarkerDB, s::SeedSettings; io::IO=stdout)
    nthreads = Threads.nthreads()
    states = [ScanState(db) for _ in 1:nthreads]
    stream = open_reads(path)
    total = 0
    live = io isa Base.TTY
    try
        while total < s.max_reads
            ids = String[]; seqs = Vector{UInt8}[]
            while length(seqs) < s.batch_size && total + length(seqs) < s.max_reads
                rec = next_fastq!(stream, path, total+length(seqs)+1)
                rec === nothing && break
                push!(ids, rec[1]); push!(seqs, rec[2])
            end
            isempty(seqs) && break
            offset = total
            chunks = min(nthreads, length(seqs))
            Threads.@threads :dynamic for chunk in 1:chunks
                st = states[chunk]
                first = fld((chunk-1)*length(seqs), chunks)+1
                last = fld(chunk*length(seqs), chunks)
                for i in first:last
                    fill!(st.scores, 0)
                    for key in canonical_kmers!(st.kmers, seqs[i], db.k)
                        e = get(db.index, key, Int32(0))
                        e == 0 && continue
                        st.matched[e] = true
                        mask = db.gene_mask[e]
                        while mask != 0
                            st.scores[trailing_zeros(mask)+1] += 1
                            mask &= mask - UInt32(1)
                        end
                    end
                    for g in eachindex(st.scores)
                        offer!(st.top[g], CANDIDATES_PER_MARKER, st.scores[g], offset+i, ids[i], seqs[i])
                    end
                end
            end
            total += length(seqs)
            if live
                print(io, "\r  scanned ", total, " reads"); flush(io)
            elseif total % (50*s.batch_size) < length(seqs)
                println(io, "  scanned ", total, " reads")
            end
            length(seqs) < s.batch_size && break
        end
    finally
        close(stream)
        live && println(io)
    end
    total > 0 || error("No reads in $path")
    matched = falses(length(db.gene_mask))
    top = [TopReads() for _ in db.genes]
    for st in states
        matched .|= st.matched
        for g in eachindex(top), j in eachindex(st.top[g].score)
            t = st.top[g]
            offer!(top[g], CANDIDATES_PER_MARKER, t.score[j], t.index[j], t.id[j], t.seq[j])
        end
    end
    group_matched = zeros(Int, length(db.group_gene))
    for e in findall(matched), j in db.group_offsets[e]+1:db.group_offsets[e+1]
        group_matched[db.group_ids[j]] += 1
    end
    return ScanResult(total, top, group_matched)
end

# ---------------------------------------------------------------------------
# Identification (reported, and cross-checked against the reference lineage)
# ---------------------------------------------------------------------------
"Rows (marker, order, matched, total, percent), best first per marker; plus the top order."
function identification(db::MarkerDB, scan::ScanResult)
    rows = Tuple{String,String,Int,Int,Float64}[]
    summed = zeros(Float64, length(db.orders))
    for grp in eachindex(db.group_gene)
        total = db.group_total[grp]; hit = scan.group_matched[grp]
        pct = total == 0 ? 0.0 : 100*hit/total
        push!(rows, (db.genes[db.group_gene[grp]], db.orders[db.group_order[grp]], hit, total, pct))
        summed[db.group_order[grp]] += pct
    end
    gene_rank = Dict(g => i for (i, g) in enumerate(db.genes))
    sort!(rows; by=r -> (gene_rank[r[1]], -r[5], -r[3], r[2]))
    top_order = ""
    if any(>(0), scan.group_matched)
        best = sort(collect(eachindex(db.orders)); by=o -> (-summed[o], db.orders[o]))[1]
        top_order = db.orders[best]
    end
    return rows, top_order
end

function write_identification(path::AbstractString, rows)
    open(path, "w") do io
        println(io, "marker\tgroup\tmatched_kmers\treference_kmers\tpercent")
        for (gene, order, hit, total, pct) in rows
            @printf(io, "%s\t%s\t%d\t%d\t%.3f\n", gene, order, hit, total, pct)
        end
    end
end

# ---------------------------------------------------------------------------
# Seeds
# ---------------------------------------------------------------------------
"Reads with >= Min k-mer hits, marker priority order: (marker, rank, read id, hits, sequence)."
const CANDIDATES_PER_MARKER = 50        # best-scoring reads kept per marker before validation
const MIN_SEED_LENGTH = 60              # shortest usable seed after cutting at an adapter

# Adapter starts (Illumina TruSeq, Nextera, MGI/DNBSEQ read 1 and 2). The data are NOT
# trimmed: only a seed candidate is cut before its adapter (degraded herbarium DNA gives
# short inserts, so most reads run into the adapter).
const ADAPTER_PREFIXES = ("AGATCGGAAGAGC", "CTGTCTCTTATACACATCT", "AAGTCGGAGGCCAAGCGG", "AAGTCGGATCGTAGCCATG")

"Length of the genomic part: up to a full adapter start, or a >= 8-base adapter start at the read end."
function genomic_length(seq::Vector{UInt8})
    t = String(copy(seq)); cut = length(seq)
    for a in ADAPTER_PREFIXES
        r = findfirst(a, t)
        r === nothing || (cut = min(cut, first(r) - 1))
        for m in length(a)-1:-1:8                      # partial adapter at the 3' end
            if length(t) >= m && endswith(t, a[1:m])
                cut = min(cut, length(t) - m); break
            end
        end
    end
    return cut
end

"Canonical 31-mers at both ends of a read (fewer if an end window contains N)."
function end_kmers(seq::Vector{UInt8})
    buf = UInt64[]; out = UInt64[]
    length(seq) >= 31 || return out
    append!(out, canonical_kmers!(buf, seq[1:31], 31))
    append!(out, canonical_kmers!(buf, seq[end-30:end], 31))
    return out
end

"Number of reads (among the first s.max_reads) containing each given canonical 31-mer (threaded)."
function count_in_reads(path::AbstractString, targets::Vector{UInt64}, s::SeedSettings)
    index = Dict(k => i for (i, k) in enumerate(unique(targets)))
    nthreads = Threads.nthreads()
    counts = [zeros(Int, length(index)) for _ in 1:nthreads]
    bufs = [UInt64[] for _ in 1:nthreads]
    stream = open_reads(path); total = 0
    try
        while total < s.max_reads
            seqs = Vector{UInt8}[]
            while length(seqs) < s.batch_size && total + length(seqs) < s.max_reads
                rec = next_fastq!(stream, path, total+length(seqs)+1)
                rec === nothing && break
                push!(seqs, rec[2])
            end
            isempty(seqs) && break
            chunks = min(nthreads, length(seqs))
            Threads.@threads :dynamic for chunk in 1:chunks
                cnt = counts[chunk]; buf = bufs[chunk]
                for i in fld((chunk-1)*length(seqs), chunks)+1:fld(chunk*length(seqs), chunks)
                    # Each k-mer counts at most once per read: support means reads, and a
                    # repeat inside one read must not look like several supporting reads.
                    for key in unique!(sort!(canonical_kmers!(buf, seqs[i], 31)))
                        j = get(index, key, 0)
                        j == 0 || (cnt[j] += 1)
                    end
                end
            end
            total += length(seqs)
            length(seqs) < s.batch_size && break
        end
    finally
        close(stream)
    end
    summed = reduce(+, counts)
    return Dict(k => summed[i] for (k, i) in index)
end

"31-mers of a read in order, with validity (windows containing N are invalid)."
function ordered_kmers(seq::Vector{UInt8})
    n = length(seq) - 30
    keys = zeros(UInt64, max(n, 0)); ok = falses(max(n, 0)); buf = UInt64[]
    for i in 1:n
        k = canonical_kmers!(buf, seq[i:i+30], 31)
        length(k) == 1 && (keys[i] = k[1]; ok[i] = true)
    end
    return keys, ok
end

"K-mer positions of `seq` that belong to marker `gene` in `db` (the scan's score)."
function marker_hits(seq::Vector{UInt8}, gene::AbstractString, db::MarkerDB)
    g = findfirst(==(gene), db.genes)
    g === nothing && return 0
    bit = UInt32(1) << (g-1); buf = UInt64[]; h = 0
    for key in canonical_kmers!(buf, seq, db.k)
        e = get(db.index, key, Int32(0))
        e != 0 && (db.gene_mask[e] & bit) != 0 && (h += 1)
    end
    return h
end

"""
Seed candidates are cut before any adapter, then shaved to their longest stretch whose
31-mers all occur in >= 2 reads (so at least one read besides the candidate itself) and
never far above typical (adapter, repeat). Sequencing errors near read ends are thereby
removed instead of disqualifying the read. Stretches shorter than 60 bp are dropped; up
to reads_per_marker seeds per marker are kept, re-ranked.
"""
function validate_seeds(seeds, r1::AbstractString, db::MarkerDB, s::SeedSettings; io::IO=stdout)
    isempty(seeds) && return seeds
    seeds = [(g, r, id, h, seq[1:genomic_length(seq)]) for (g, r, id, h, seq) in seeds]
    kmers = [ordered_kmers(x[5]) for x in seeds]
    counts = count_in_reads(r1, reduce(vcat, [k[1][k[2]] for k in kmers]; init=UInt64[]), s)
    all_counts = sort([counts[x] for (k, ok) in kmers for (x, v) in zip(k, ok) if v])
    typical = isempty(all_counts) ? 0 : all_counts[(length(all_counts)+1)÷2]
    kept = Tuple{String,String,Int,Vector{UInt8},Int}[]      # gene, id, hits on stretch, stretch, order
    rejected = Dict{String,Int}(); shaved = 0
    for (i, (gene, _, id, hits, seq)) in enumerate(seeds)
        keys, ok = kmers[i]
        good = [ok[j] && 2 <= counts[keys[j]] <= max(20*typical, 2) for j in eachindex(keys)]
        # Among supported stretches >= MIN_SEED_LENGTH, keep the one that still matches the
        # marker best (hits are re-counted on the stretch; the read's original score no
        # longer applies once it has been shortened).
        best = nothing; j = 1
        while j <= length(good)
            if good[j]
                e = j
                while e < length(good) && good[e+1]
                    e += 1
                end
                if e + 30 - j + 1 >= MIN_SEED_LENGTH
                    stretch = seq[j:e+30]
                    h = marker_hits(stretch, gene, db)
                    (best === nothing || (h, length(stretch)) > (best[1], length(best[2]))) && (best = (h, stretch))
                end
                j = e + 1
            else
                j += 1
            end
        end
        if best === nothing || best[1] < s.min_hits
            rejected[gene] = get(rejected, gene, 0) + 1
            continue
        end
        length(best[2]) < length(seq) && (shaved += 1)
        push!(kept, (gene, id, best[1], best[2], i))
    end
    out = eltype(seeds)[]
    for gene in unique(first.(kept))
        list = sort([k for k in kept if k[1] == gene]; by=k -> (-k[3], k[5]))
        for (r, (g, id, h, stretch, _)) in enumerate(first(list, s.reads_per_marker))
            push!(out, (g, r, id, h, stretch))
        end
    end
    for (gene, n) in sort(collect(rejected))
        println(io, "  seed check: ", n, " ", gene, " candidate(s) without a supported stretch >= ", MIN_SEED_LENGTH,
                " bp that still has >= ", s.min_hits, " marker hits")
    end
    shaved > 0 && println(io, "  seed check: ", shaved, " kept seed(s) shortened to their read-supported part")
    return out
end

function selected_seed_reads(db::MarkerDB, scan::ScanResult, s::SeedSettings)
    out = Tuple{String,Int,String,Int,Vector{UInt8}}[]
    for (g, gene) in enumerate(db.genes), j in eachindex(scan.top[g].score)
        t = scan.top[g]
        t.score[j] >= s.min_hits && push!(out, (gene, j, t.id[j], t.score[j], t.seq[j]))
    end
    return out
end

"Assembler seed FASTA: the best read of each marker, in marker priority order."
function write_seed_fasta(path::AbstractString, seeds)
    written = 0
    open(path, "w") do io
        for (gene, rank, id, hits, seq) in seeds
            rank == 1 || continue
            println(io, ">", gene, "_seed hits=", hits, " read=", id)
            println(io, String(copy(seq)))
            written += 1
        end
    end
    return written
end

function write_seed_table(path::AbstractString, seeds)
    open(path, "w") do io
        println(io, "marker\trank\tread_id\tkmer_hits\tlength\tused_as_seed")
        for (gene, rank, id, hits, seq) in seeds
            println(io, join((gene, rank, id, hits, length(seq), rank == 1), '\t'))
        end
    end
end

# ---------------------------------------------------------------------------
# NCBI access with Downloads (stdlib). BLAST URL API usage rules: at most one
# request every 10 s, poll a RID at most once a minute, identify the tool.
# ---------------------------------------------------------------------------
const BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
const EFETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
const ESEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
const LAST_BLAST_REQUEST = Ref(0.0)
const LAST_EUTILS_REQUEST = Ref(0.0)     # E-utilities without API key: <= 3 requests/s.

function url_encode(text::AbstractString)
    io = IOBuffer()
    for b in codeunits(text)
        if UInt8('A') <= b <= UInt8('Z') || UInt8('a') <= b <= UInt8('z') ||
           UInt8('0') <= b <= UInt8('9') || b in (UInt8('-'), UInt8('_'), UInt8('.'), UInt8('~'))
            write(io, b)
        else
            print(io, '%', uppercase(string(b; base=16, pad=2)))
        end
    end
    return String(take!(io))
end
form(params) = join((url_encode(k)*"="*url_encode(v) for (k, v) in params), '&')

function http_text(url::AbstractString; body::Union{Nothing,String}=nothing, timeout::Real=180)
    out = IOBuffer()
    response = if body === nothing
        Downloads.request(url; output=out, timeout=timeout)
    else
        Downloads.request(url; method="POST", input=IOBuffer(body), output=out, timeout=timeout,
                          headers=["Content-Type" => "application/x-www-form-urlencoded"])
    end
    text = String(take!(out))
    response.status == 200 || error("HTTP $(response.status) from $(first(split(url, '?')))")
    return text
end

function blast_text(params; body::Bool=false)
    wait = LAST_BLAST_REQUEST[] + 10.0 - time()
    wait > 0 && sleep(wait)
    try
        return body ? http_text(BLAST_URL; body=form(params)) : http_text(BLAST_URL*"?"*form(params))
    finally
        LAST_BLAST_REQUEST[] = time()
    end
end

identity_params(s::SeedSettings) =
    ["TOOL" => s.tool]

"Submit, poll and fetch one BLAST search; returns the Tabular text (empty if no hits)."
function remote_blast(query_fasta::String, s::SeedSettings; io::IO=stdout)
    put = ["CMD" => "Put", "PROGRAM" => "blastn", "DATABASE" => s.blast_database,
           "QUERY" => query_fasta, "HITLIST_SIZE" => string(s.blast_hitlist)]
    isempty(s.blast_entrez_query) || push!(put, "ENTREZ_QUERY" => s.blast_entrez_query)
    text = blast_text(vcat(put, identity_params(s)); body=true)
    rid = match(r"RID = (\S+)", text)
    rid === nothing && error("BLAST submission returned no RID")
    rid = rid.captures[1]
    rtoe = match(r"RTOE = (\d+)", text)
    first_wait = rtoe === nothing ? 30 : clamp(parse(Int, rtoe.captures[1]), 10, 300)
    println(io, "  BLAST RID ", rid, "; first status check in ", first_wait, " s")
    started = time()
    sleep(first_wait)
    while true
        info = blast_text(vcat(["CMD" => "Get", "FORMAT_OBJECT" => "SearchInfo", "RID" => rid], identity_params(s)))
        status = match(r"Status=(\w+)", info)
        status === nothing && error("BLAST status for $rid could not be read")
        status = status.captures[1]
        if status == "READY"
            occursin("ThereAreHits=yes", info) || return ""
            break
        elseif status in ("FAILED", "UNKNOWN")
            error("BLAST search $rid ended with status $status")
        end
        time()-started > s.blast_timeout && error("BLAST search $rid exceeded $(s.blast_timeout) s")
        println(io, "  BLAST ", rid, " waiting (", round(Int, time()-started), " s)")
        sleep(60)                                   # Once a minute per RID.
    end
    return blast_text(vcat(["CMD" => "Get", "RID" => rid, "FORMAT_TYPE" => "Tabular", "RESULTS_FILE" => "on",
                            "FORMAT_OBJECT" => "Alignment", "ALIGNMENTS" => string(s.blast_hitlist),
                            "DESCRIPTIONS" => string(s.blast_hitlist)], identity_params(s)))
end

function clean_accession(field::AbstractString)
    parts = split(field, '|'; keepempty=false)
    if length(parts) >= 2 && parts[1] in ("ref", "gb", "emb", "dbj", "gi")
        parts[1] == "gi" && length(parts) >= 4 && return String(parts[4])
        return String(parts[2])
    end
    return String(field)
end

"Parse BLAST Tabular output via its '# Fields:' line (standard 12-column order as fallback)."
function parse_blast_tabular(text::AbstractString)
    fields = String[]
    hits = Tuple{String,String,Float64,Float64}[]        # query, subject, identity, bitscore
    for raw in eachline(IOBuffer(text))
        line = strip(replace(raw, r"<[^>]*>" => ""))      # tolerate <PRE>-wrapped output
        isempty(line) && continue
        if startswith(line, "#")
            m = match(r"^#\s*Fields:\s*(.*)$", line)
            m === nothing || (fields = String.(strip.(split(m.captures[1], ','))))
            continue
        end
        parts = split(line, '\t')
        length(parts) >= 12 || continue
        col(names, default) = something(findfirst(in(names), fields), default)
        q = col(("query acc.ver", "query id", "query acc."), 1)
        sbj = col(("subject acc.ver", "subject id", "subject acc.", "subject ids"), 2)
        pid = col(("% identity",), 3)
        bits = col(("bit score",), 12)
        push!(hits, (String(parts[q]), clean_accession(parts[sbj]),
                     parse(Float64, parts[pid]), parse(Float64, parts[bits])))
    end
    if isempty(hits) && occursin("<html", lowercase(text))
        error("BLAST returned HTML instead of Tabular output")
    end
    return hits
end

"Summed best bit score per query; ties: more queries hit, higher mean identity, accession."
function rank_accessions(hits)
    best = Dict{String,Dict{String,Tuple{Float64,Float64}}}()
    for (query, subject, pid, bits) in hits
        per = get!(Dict{String,Tuple{Float64,Float64}}, best, subject)
        old = get(per, query, (-Inf, 0.0))
        bits > old[1] && (per[query] = (bits, pid))
    end
    ranked = [(acc, sum(first, values(per)), length(per), sum(last, values(per))/length(per))
              for (acc, per) in best]
    sort!(ranked; by=x -> (-x[2], -x[3], -x[4], x[1]))
    return ranked
end

function eutils_text(url::AbstractString, params)
    wait = LAST_EUTILS_REQUEST[] + 0.4 - time()
    wait > 0 && sleep(wait)
    full = vcat(params, ["tool" => "AssemblerCompleteV16"])
    try
        return http_text(url*"?"*form(full))
    finally
        LAST_EUTILS_REQUEST[] = time()
    end
end

"First nuccore UID for an Entrez term, or nothing."
function esearch_first(term::AbstractString)
    text = eutils_text(ESEARCH_URL, ["db" => "nuccore", "term" => term, "retmax" => "1"])
    m = match(r"<Id>(\d+)</Id>", text)
    return m === nothing ? nothing : String(m.captures[1])
end

function efetch_genbank(accession::AbstractString, s::SeedSettings)
    text = eutils_text(EFETCH_URL, ["db" => "nuccore", "id" => accession, "rettype" => "gbwithparts", "retmode" => "text"])
    startswith(lstrip(text), "LOCUS") || error("efetch for $accession did not return a GenBank record")
    return text
end

# ---------------------------------------------------------------------------
# GenBank header (length, topology, organism) - features are left for annotation.
# ---------------------------------------------------------------------------
struct ReferenceInfo
    path::String
    accession::String
    definition::String
    organism::String
    lineage::String
    length::Int
    topology::String
    records::Int
end

function read_genbank_header(path::AbstractString)
    len = 0; topology = ""; accession = ""; definition = ""; organism = ""
    lineage = IOBuffer(); records = 0; field = ""
    for line in eachline(path)
        if startswith(line, "LOCUS")
            records += 1
            records > 1 && continue
            m = match(r"\s(\d+)\s+bp", line)
            m === nothing || (len = parse(Int, m.captures[1]))
            topology = occursin(r"\bcircular\b"i, line) ? "circular" : "linear"
            field = ""
            continue
        end
        records == 1 || continue
        if startswith(line, "FEATURES") || startswith(line, "ORIGIN")
            field = ""
        elseif startswith(line, "DEFINITION")
            definition = strip(line[11:end]); field = "definition"
        elseif startswith(line, "VERSION")
            tokens = split(line)
            length(tokens) >= 2 && (accession = String(tokens[2])); field = ""
        elseif startswith(line, "ACCESSION") && isempty(accession)
            tokens = split(line)
            length(tokens) >= 2 && (accession = String(tokens[2])); field = ""
        elseif startswith(line, "  ORGANISM")
            organism = strip(line[11:end]); field = "lineage"
        elseif startswith(line, "            ") && field == "definition"
            definition *= " "*strip(line)
        elseif startswith(line, "            ") && field == "lineage"
            print(lineage, strip(line), " ")
        elseif !startswith(line, " ")
            field = ""
        end
    end
    len > 0 || error("$path: no sequence length in the LOCUS line")
    return ReferenceInfo(String(path), accession, String(definition), String(organism),
                         String(strip(String(take!(lineage)))), len, topology, records)
end

genome_range(ref::ReferenceInfo, margin::Int) = (max(1000, ref.length-margin), ref.length+margin)

"Sequence and features (key, first /gene, segments) of the first record in a GenBank file."
function read_genbank_record(path::AbstractString)
    seq = IOBuffer()
    features = Tuple{String,String,Vector{UnitRange{Int}}}[]
    state = :header
    key = ""; gene = ""; loc = IOBuffer(); inloc = false
    function flush_feature!()
        isempty(key) && return
        location = String(take!(loc))
        if !occursin(':', location)                       # remote locations unsupported
            segments = UnitRange{Int}[parse(Int, m.captures[1]):parse(Int, m.captures[2])
                                      for m in eachmatch(r"<?(\d+)\.\.>?(\d+)", location)]
            push!(features, (key, gene, segments))
        end
        key = ""; gene = ""
    end
    for line in eachline(path)
        if state == :header
            startswith(line, "FEATURES") && (state = :features)
        elseif state == :features
            if startswith(line, "ORIGIN")
                flush_feature!(); state = :origin; continue
            end
            bytes = codeunits(line)
            length(bytes) >= 22 || continue
            fkey = strip(String(bytes[6:min(20, end)]))
            body = strip(String(bytes[22:end]))
            if !isempty(fkey)
                flush_feature!()
                key = fkey; write(loc, body); inloc = true
            elseif startswith(body, "/")
                inloc = false
                m = match(r"^/gene=\"?([^\"]*)\"?", body)
                m !== nothing && isempty(gene) && (gene = strip(m.captures[1]))
            elseif inloc
                write(loc, body)
            end
        elseif state == :origin
            startswith(line, "//") && break
            for b in codeunits(line)
                UInt8('a') <= b <= UInt8('z') && (b -= 0x20)
                UInt8('A') <= b <= UInt8('Z') && write(seq, b)
            end
        end
    end
    return take!(seq), features
end

"Sequences of the wanted genes (CDS preferred over gene/rRNA features; first copy only)."
function genbank_genes(path::AbstractString, wanted::Vector{String})
    seq, features = read_genbank_record(path)
    lower = Dict(lowercase(g) => g for g in wanted)
    found = Dict{String,Vector{UInt8}}()
    for pass in ("CDS", "gene", "rRNA")
        for (key, gene, segments) in features
            key == pass || continue
            name = get(lower, lowercase(gene), "")
            (isempty(name) || haskey(found, name) || isempty(segments)) && continue
            all(r -> 1 <= first(r) <= last(r) <= length(seq), segments) || continue
            found[name] = reduce(vcat, [seq[r] for r in segments])
        end
    end
    return found, seq
end

function write_fasta_from_genbank(path::AbstractString, genbank::AbstractString, header::AbstractString)
    seq, _ = read_genbank_record(genbank)
    open(path, "w") do io
        println(io, ">", header)
        for i in 1:80:length(seq)
            println(io, String(seq[i:min(i+79, end)]))
        end
    end
    return length(seq)
end

# ---------------------------------------------------------------------------
# Automatic marker panel: one RefSeq plastome per major land-plant lineage.
# Organism names (not accession numbers) are resolved with esearch; a name that
# returns nothing is skipped with a warning. Downloaded once into the cache.
# ---------------------------------------------------------------------------
const PANEL = [                  # Names verified with esearch on 2026-09-26.
    ("liverwort_Marchantia", ["Marchantia polymorpha"]),
    ("liverwort_Apopellia", ["Apopellia endiviifolia", "Pellia endiviifolia"]),
    ("liverwort_Ptilidium", ["Ptilidium pulcherrimum"]),
    ("liverwort_Aneura", ["Aneura mirabilis"]),
    ("moss_Physcomitrium", ["Physcomitrium patens", "Physcomitrella patens"]),
    ("moss_Sphagnum", ["Sphagnum palustre"]),
    ("moss_Tetraphis", ["Tetraphis pellucida"]),
    ("hornwort_Anthoceros", ["Anthoceros angustus", "Anthoceros formosae"]),
    ("lycophyte_Huperzia", ["Huperzia lucidula"]),
    ("fern_Equisetum", ["Equisetum arvense"]),
    ("fern_Adiantum", ["Adiantum capillus-veneris"]),
    ("gymnosperm_Pinus", ["Pinus thunbergii"]),
    ("gymnosperm_Cycas", ["Cycas taitungensis"]),
    ("angiosperm_Amborella", ["Amborella trichopoda"]),
    ("angiosperm_Liriodendron", ["Liriodendron tulipifera"]),
    ("angiosperm_Oryza", ["Oryza sativa"]),
    ("angiosperm_Arabidopsis", ["Arabidopsis thaliana"]),
    ("angiosperm_Nicotiana", ["Nicotiana tabacum"]),
]

function panel_records(s::SeedSettings; io::IO=stdout)
    cache = joinpath(isempty(s.reference_cache) ? joinpath(pwd(), "seedref_cache") : s.reference_cache, "panel")
    mkpath(cache)
    records = MarkerRecord[]
    for (label, names) in PANEL
        path = joinpath(cache, label*".gb")
        if !isfile(path)
            try
                organisms = join(("$(n)[Organism]" for n in names), " OR ")
                uid = esearch_first("($organisms) AND " * s.blast_entrez_query)
                if uid === nothing
                    println(io, "  panel: no RefSeq plastome for ", join(names, " / "), "; skipped")
                    continue
                end
                write(path, efetch_genbank(uid, s))
            catch err
                println(io, "  panel: ", label, " not downloaded (", first(split(sprint(showerror, err), '\n')), ")")
                continue
            end
        end
        genes, _ = genbank_genes(path, s.markers)
        for gene in s.markers
            haskey(genes, gene) && push!(records, (gene, label, genes[gene]))
        end
    end
    groups = length(unique(r[2] for r in records))
    println(io, "Marker panel: ", groups, " plastomes, ", length(records), " gene sequences (cache ", cache, ")")
    return records
end

# ---------------------------------------------------------------------------
# Per-sample preparation
# ---------------------------------------------------------------------------
struct SamplePrep
    seed_fasta::String               # Best read of every marker gene (assembler, several seeds).
    seed_single::String              # Single best read of the highest-priority marker.
    seed_marker::String
    seed_hits::Int
    seeds_written::Int
    reads_scanned::Int
    top_group::String                # Best-matching panel plastome / marker group.
    reference::Union{Nothing,ReferenceInfo}
    reference_status::String         # blast | cached | local | reused | none | failed: ...
    genome_range::Union{Nothing,Tuple{Int,Int}}
    refined::Bool                    # Final seeds come from the reference's own genes.
end

function write_blast_table(path::AbstractString, ranked)
    open(path, "w") do io
        println(io, "accession\tsummed_best_bitscore\tqueries_hit\tmean_identity")
        for (acc, bits, n, pid) in ranked
            @printf(io, "%s\t%.1f\t%d\t%.2f\n", acc, bits, n, pid)
        end
    end
end

"Reference for one sample. Never throws for network/BLAST problems; returns a status."
function obtain_reference(sample::AbstractString, seeds, s::SeedSettings, refdir::AbstractString;
                          local_genbank::AbstractString="", io::IO=stdout)
    s.reference_mode == "none" && return nothing, "none"
    target = joinpath(refdir, "reference.gb")
    if isfile(target)
        return read_genbank_header(target), "reused"
    end
    for (path, status) in ((local_genbank, "local"), (s.reference_genbank, "local"))
        if !isempty(path) && isfile(path)
            cp(path, target; force=true)
            return read_genbank_header(target), status
        end
    end
    s.reference_mode == "local" && return nothing, "failed: Reference GenBank not found"
    isempty(seeds) && return nothing, "failed: no seed reads for BLAST"
    try
        query = join((">$(gene)_$(rank)\n$(String(copy(seq)))" for (gene, rank, _, _, seq) in seeds), '\n')
        println(io, "  BLAST: ", length(seeds), " seed reads, one request")
        ranked = rank_accessions(parse_blast_tabular(remote_blast(query, s; io=io)))
        write_blast_table(joinpath(refdir, "blast_hits.tsv"), ranked)
        isempty(ranked) && return nothing, "failed: no BLAST hits within the Entrez query"
        accession = ranked[1][1]
        length(ranked) > 1 && ranked[2][2] == ranked[1][2] &&
            println(io, "  note: ", accession, " tied with ", ranked[2][1], "; chosen by queries hit/identity/accession")
        cached = isempty(s.reference_cache) ? "" : joinpath(s.reference_cache, accession*".gb")
        if !isempty(cached) && isfile(cached)
            cp(cached, target; force=true)
            return read_genbank_header(target), "cached"
        end
        text = efetch_genbank(accession, s)
        write(target, text)
        isempty(cached) || (mkpath(s.reference_cache); write(cached, text))
        return read_genbank_header(target), "blast"
    catch err
        msg = sprint(showerror, err)
        return nothing, "failed: "*first(replace(msg, r"\s+" => " "), 200)
    end
end

function write_single_seed(path::AbstractString, seeds)
    for (gene, rank, id, hits, seq) in seeds
        rank == 1 || continue
        open(path, "w") do io
            println(io, ">", gene, "_seed hits=", hits, " read=", id)
            println(io, String(copy(seq)))
        end
        return gene, hits
    end
    return "", 0
end

"Seeds, identification and reference for one sample, written under `outdir`."
function prepare_sample(sample::AbstractString, r1::AbstractString, outdir::AbstractString,
                        db::MarkerDB, s::SeedSettings; local_genbank::AbstractString="", io::IO=stdout)
    seeddir = joinpath(outdir, "seed"); refdir = joinpath(outdir, "reference")
    mkpath(seeddir); mkpath(refdir)
    println(io, "  seed scan: first ", s.max_reads, " reads of ", basename(r1), " (k=", db.k, ")")
    scan = scan_reads(r1, db, s; io=io)
    rows, top_group = identification(db, scan)
    write_identification(joinpath(seeddir, "identification.tsv"), rows)
    seeds = validate_seeds(selected_seed_reads(db, scan, s), r1, db, s; io=io)
    isempty(seeds) && error("No valid seed: no candidate with >= $(MIN_SEED_LENGTH) bp supported by other reads (try a larger SeedIR reads)")
    isempty(top_group) || println(io, "  closest panel group: ", top_group)
    reference, status = obtain_reference(sample, seeds, s, refdir; local_genbank=local_genbank, io=io)
    final = seeds; refined = false
    if reference !== nothing && s.refine
        genes, _ = genbank_genes(reference.path, s.markers)
        if !isempty(genes)
            records = MarkerRecord[(g, reference.accession, genes[g]) for g in s.markers if haskey(genes, g)]
            db2 = build_marker_db(records, db.k, s.markers)
            println(io, "  refinement scan with ", length(records), " genes of ", reference.accession)
            scan2 = scan_reads(r1, db2, s; io=io)
            seeds2 = validate_seeds(selected_seed_reads(db2, scan2, s), r1, db2, s; io=io)
            if !isempty(seeds2)
                final = seeds2; refined = true
            else
                println(io, "  refinement found no read above Min k-mer hits; first-pass seeds kept")
            end
        else
            println(io, "  reference has none of the marker genes annotated; first-pass seeds kept")
        end
    end
    open(joinpath(seeddir, "seed_reads.tsv"), "w") do tsv
        println(tsv, "pass\tmarker\trank\tread_id\tkmer_hits\tlength\tused_as_seed")
        for (pass, list, used) in (("panel", seeds, !refined), ("reference", refined ? final : [], refined))
            for (gene, rank, id, hits, seq) in list
                println(tsv, join((pass, gene, rank, id, hits, length(seq), used && rank == 1), '\t'))
            end
        end
    end
    seed_fasta = joinpath(seeddir, "seed.fasta")
    written = write_seed_fasta(seed_fasta, final)
    seed_single = joinpath(seeddir, "seed_single.fasta")
    marker, hits = write_single_seed(seed_single, final)
    for (gene, rank, id, h, _) in final
        rank == 1 && println(io, "  seed ", gene, ": ", h, " k-mer hits (", id, ")", refined ? " [reference genes]" : " [panel]")
    end
    range = nothing
    if reference !== nothing
        range = genome_range(reference, s.genome_margin)
        println(io, "  reference ", reference.accession, " (", reference.organism, "), ",
                reference.length, " bp ", reference.topology, " [", status, "] -> Genome Range ",
                range[1], "-", range[2])
        reference.records > 1 && println(io, "  note: reference file has ", reference.records,
                                         " records; the first one was used")
        reference.topology == "circular" || println(io, "  note: reference LOCUS is not circular")
    else
        println(io, "  reference: ", status)
    end
    return SamplePrep(seed_fasta, seed_single, marker, hits, written, scan.reads, top_group,
                      reference, status, range, refined)
end

# ---------------------------------------------------------------------------
# Sample discovery: paired FASTQ(.gz) in the input folder and its direct subfolders.
# R1/R2 tokens: _R1/_R2 (optionally _001) or _1/_2, also with '.' or '-'.
# ---------------------------------------------------------------------------
struct Sample
    name::String
    r1::String
    r2::String
    genbank::String     # <sample>.gb/.gbk/.gbff beside the reads, or ""
end

const FASTQ_RE = r"\.(fastq|fq)(\.gz)?$"i
const MATE_PATTERNS = (r"^(.+?)[._-]R([12])(_\d{3})?$", r"^(.+?)[._-]([12])$")

sanitize(name::AbstractString) = replace(name, r"[^A-Za-z0-9_.-]" => "_")

function mate_of(stem::AbstractString)
    for pattern in MATE_PATTERNS
        m = match(pattern, stem)
        m === nothing || return String(m.captures[1]), parse(Int, m.captures[2])
    end
    return String(stem), 0
end

function local_genbank(dir::AbstractString, prefix::AbstractString)
    for ext in (".gb", ".gbk", ".gbff")
        path = joinpath(dir, prefix*ext)
        isfile(path) && return path
    end
    return ""
end

function discover_samples(input::AbstractString)
    dirs = [input; [joinpath(input, d) for d in sort(readdir(input)) if isdir(joinpath(input, d))]]
    groups = Dict{Tuple{String,String},Dict{Int,Vector{String}}}()
    for dir in dirs, file in sort(readdir(dir))
        path = joinpath(dir, file)
        isfile(path) && occursin(FASTQ_RE, file) || continue
        prefix, mate = mate_of(replace(file, FASTQ_RE => ""))
        push!(get!(Vector{String}, get!(Dict{Int,Vector{String}}, groups, (dir, prefix)), mate), path)
    end
    samples = Sample[]
    for ((dir, prefix), mates) in sort(collect(groups); by=first)
        for (mate, files) in mates
            length(files) > 1 && error("Several FASTQ files map to $(prefix) mate $(mate): $(join(basename.(files), ", "))")
        end
        # Paired-end only: every library needs exactly one R1 and one R2 file.
        (haskey(mates, 1) && haskey(mates, 2) && !haskey(mates, 0)) ||
            error("Library '$(prefix)' in $(dir) is not a complete R1/R2 pair: " *
                  join(basename.(reduce(vcat, values(mates))), ", "))
        r1 = only(mates[1]); r2 = only(mates[2])
        push!(samples, Sample(sanitize(prefix), r1, r2, local_genbank(dir, prefix)))
    end
    # Same library name in two subfolders: prefix with the folder name.
    counts = Dict{String,Int}()
    for x in samples
        counts[x.name] = get(counts, x.name, 0)+1
    end
    samples = [counts[x.name] == 1 ? x :
               Sample(sanitize(basename(dirname(x.r1))*"_"*x.name), x.r1, x.r2, x.genbank) for x in samples]
    length(unique(x.name for x in samples)) == length(samples) || error("Duplicate sample names after renaming")
    return sort!(samples; by=x -> x.name)
end


# ---------------------------------------------------------------------------
# Folder mode: every library in a folder; outputs named after the library.
# ---------------------------------------------------------------------------
const PER_LIBRARY_KEYS = ("Project name", "Forward reads", "Reverse reads", "Single/Paired",
                          "Seed Input", "Output path", "Genome Range", "Plastid panel")
const WIDE_GENOME_RANGE = "100000-180000"     # land-plant plastomes when no reference is available

"Template config with the per-library keys replaced (comments and all other settings kept)."
function write_assembler_config(path, template, lib, r1, r2, seed, results, range; panel="")
    open(path, "w") do io
        for line in eachline(template)
            body = strip(first(split(line, '#'; limit=2)))
            key = occursin('=', body) ? strip(first(split(body, '='; limit=2))) : ""
            key in PER_LIBRARY_KEYS && continue
            println(io, line)
        end
        println(io, "\n# --- written by SeedRef ", VERSION_STRING, " for ", lib, " ---")
        println(io, "Project name = ", lib)
        println(io, "Single/Paired = PE")
        println(io, "Forward reads = ", r1)
        println(io, "Reverse reads = ", r2)
        println(io, "Seed Input = ", seed)
        println(io, "Output path = ", joinpath(results, lib))
        println(io, "Genome Range = ", range)
        isempty(panel) || println(io, "Plastid panel = ", panel)
    end
end

function run_folder(input::AbstractString, output::AbstractString, s::SeedSettings; only_names=String[],
                    template::AbstractString="", results::AbstractString="")
    samples = discover_samples(input)
    isempty(only_names) || (samples = [x for x in samples if x.name in only_names])
    isempty(samples) && error("No FASTQ libraries found in $input")
    mkpath(output)
    isempty(s.reference_cache) && (s.reference_cache = joinpath(output, "references"))
    @printf("SeedRef %s | %d threads | %d libraries in %s\n", VERSION_STRING, Threads.nthreads(), length(samples), input)
    db = build_marker_db(s)
    println("Marker k-mers: ", length(db.index), " (", join(db.genes, ", "), "; k=", db.k, ")")
    cols = ["library", "r1", "r2", "seed", "seeds_all", "seed_marker", "seed_hits",
            "seed_source", "closest_panel_group", "reference_status", "reference_accession",
            "reference_organism", "reference_length", "reference_topology", "genome_range",
            "reference_genbank", "reference_fasta", "read_length", "assembler_config", "error"]
    rows = Vector{Vector{String}}()
    summary = joinpath(output, "seedref_summary.tsv")
    for (i, x) in enumerate(samples)
        println("\n[", i, "/", length(samples), "] ", x.name)
        row = Dict(c => "" for c in cols)
        row["library"] = x.name
        row["r1"] = x.r1; row["r2"] = x.r2
        try
            detail = joinpath(output, x.name)
            prep = prepare_sample(x.name, x.r1, detail, db, s; local_genbank=x.genbank)
            single = joinpath(output, x.name*"_seed.fasta")
            multi = joinpath(output, x.name*"_seeds_all.fasta")
            cp(prep.seed_single, single; force=true); cp(prep.seed_fasta, multi; force=true)
            row["seed"] = single; row["seeds_all"] = multi
            row["seed_marker"] = prep.seed_marker; row["seed_hits"] = string(prep.seed_hits)
            row["seed_source"] = prep.refined ? "reference_genes" : "panel"
            row["closest_panel_group"] = prep.top_group
            row["reference_status"] = prep.reference_status
            row["read_length"] = string(first_read_length(x.r1))
            if prep.reference !== nothing
                ref = prep.reference
                gb = joinpath(output, x.name*"_reference.gb")
                fa = joinpath(output, x.name*"_reference.fasta")
                cp(ref.path, gb; force=true)
                write_fasta_from_genbank(fa, gb, "$(ref.accession) $(ref.organism)")
                row["reference_accession"] = ref.accession; row["reference_organism"] = ref.organism
                row["reference_length"] = string(ref.length); row["reference_topology"] = ref.topology
                row["genome_range"] = "$(prep.genome_range[1])-$(prep.genome_range[2])"
                row["reference_genbank"] = gb; row["reference_fasta"] = fa
            end
            if !isempty(template)
                cfg = joinpath(output, x.name*"_config.txt")
                write_assembler_config(cfg, template, x.name, x.r1, x.r2, single, results,
                                       isempty(row["genome_range"]) ? WIDE_GENOME_RANGE : row["genome_range"];
                                       panel=s.panel ? joinpath(s.reference_cache, "panel") : "")
                row["assembler_config"] = cfg
            end
        catch err
            err isa InterruptException && rethrow()
            msg = first(split(sprint(showerror, err), '\n'))
            println("  FAILED: ", msg)
            row["error"] = replace(msg, r"[\t\r\n]" => " ")
        end
        push!(rows, [row[c] for c in cols])
        open(summary, "w") do io                    # Rewritten after every library.
            println(io, join(cols, '\t'))
            for r in rows
                println(io, join(r, '\t'))
            end
        end
    end
    println("\nSummary: ", summary)
    if !isempty(template)
        println("Assembler: julia --project=. --threads 64 AssemblerComplete_v16.jl ", joinpath(output, "<library>_config.txt"))
    end
end

function first_read_length(path::AbstractString)
    stream = open_reads(path)
    try
        rec = next_fastq!(stream, path, 1)
        return rec === nothing ? 0 : length(rec[2])
    finally
        close(stream)
    end
end

function main(args=ARGS)
    if isempty(args) || args[1] in ("-h", "--help")
        println("Folder mode:  julia --project=. --threads auto SeedRef.jl --input DIR --output DIR [--config FILE] [--template ASSEMBLER_CONFIG] [--results DIR] [--samples a,b]")
        println("Single mode:  julia --project=. --threads auto SeedRef.jl config_v16.txt")
        return
    end
    if startswith(args[1], "--")
        opts = Dict{String,String}()
        i = 1
        while i <= length(args)
            startswith(args[i], "--") && i < length(args) || error("Options are --name value pairs: $(args[i])")
            opts[args[i][3:end]] = args[i+1]; i += 2
        end
        haskey(opts, "input") && haskey(opts, "output") || error("--input and --output are required")
        if haskey(opts, "config")
            cfg = abspath(opts["config"])
            s, _ = seed_settings(read_key_values(cfg), dirname(cfg))
        else
            s, _ = seed_settings(Pair{String,String}[], pwd())
        end
        names = haskey(opts, "samples") ? String.(strip.(split(opts["samples"], ','))) : String[]
        template = haskey(opts, "template") ? abspath(opts["template"]) : ""
        isempty(template) || isfile(template) || error("Template config not found: $template")
        results = abspath(get(opts, "results", "results"))
        run_folder(abspath(opts["input"]), abspath(opts["output"]), s; only_names=names,
                   template=template, results=results)
        return
    end
    path = abspath(args[1]); base = dirname(path)
    pairs = read_key_values(path)
    s, _ = seed_settings(pairs, base)
    isempty(s.reference_cache) && (s.reference_cache = joinpath(base, "seedref_cache"))
    kv = Dict(pairs)
    haskey(kv, "Forward reads") || error("Config has no Forward reads")
    r1 = resolve_path(base, kv["Forward reads"])
    name = get(kv, "Project name", splitext(basename(r1))[1])
    prep = prepare_sample(name, r1, base, build_marker_db(s), s)
    if haskey(kv, "Seed Input")
        target = resolve_path(base, kv["Seed Input"])
        target == prep.seed_fasta || cp(prep.seed_fasta, target; force=true)
        println("Seed FASTA: ", target)
    end
    prep.genome_range === nothing ||
        println("Set in the config: Genome Range = ", prep.genome_range[1], "-", prep.genome_range[2])
end

end # module SeedRef

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    SeedRef.main()
end
