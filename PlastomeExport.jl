#!/usr/bin/env julia
# PlastomeExport.jl - gene tables, GenBank submission files and genome maps from AssemblerComplete results
# (Base Julia only).
#
#   julia PlastomeExport.jl gene-table --out DIR [--categories FILE] GENES_TSV [GENES_TSV ...]
#   julia PlastomeExport.jl genbank --template METADATA_TSV GENES_TSV [...]
#   julia PlastomeExport.jl genbank --metadata METADATA_TSV --out DIR [--products FILE] GENES_TSV [...]
#   julia PlastomeExport.jl map --out DIR [--metadata METADATA_TSV] GENES_TSV [...]
#   julia PlastomeExport.jl ir-plot --svg FILE [--ref ACC,...] [--metadata METADATA_TSV] GENES_TSV [...]
#
# GENES_TSV is <prefix>_plastome_genes.tsv written by AssemblerComplete (annotation of the
# ordered plastome); <prefix>_plastome_LSC_IR_SSC_IR.fasta.regions.tsv must sit next to it.
# See the docstrings of gene_table, write_template, genbank, plot_genome_map and ir_plot for the outputs.

module PlastomeExport

using Printf

const VERSION_STRING = "0.4.4"

# A cis intron longer than this (along the strand) is not expected in land-plant plastomes
# (the longest group II introns, e.g. trnK-UUU, are ~2.5 kb); a longer or strand-switching
# junction is treated as trans-splicing.
const MAX_CIS_INTRON = 5000

struct Segment
    first::Int        # genomic coordinates of the ordered FASTA; first > last crosses the origin
    last::Int
    strand::Char
end

struct Model
    record::String
    gene::String
    kind::String      # CDS, tRNA, rRNA
    strand::Char
    segs::Vector{Segment}   # transcript order
    parts::Vector{Int}
    p::Float64
    pseudo::Bool
    flags::Vector{String}
end

struct Regions
    length::Int
    names::Vector{String}                  # LSC, IR1, SSC, IR2 (file order)
    spans::Vector{Tuple{Int,Int}}
end

struct Category
    category::String
    subgroup::String
end

# ------------------------------------------------------------------ input

function read_tsv_rows(path::AbstractString)
    rows = Vector{Dict{String,String}}()
    header = String[]
    for line in eachline(path)
        (isempty(strip(line)) || startswith(line, "#")) && continue
        x = split(line, '\t')
        if isempty(header)
            header = String.(x)
            continue
        end
        # trailing empty fields may be dropped by editors (e.g. empty flags); longer rows are errors
        length(x) <= length(header) || error("$path: $(length(x)) columns, header has $(length(header)): $line")
        push!(rows, Dict(header[i] => (i <= length(x) ? String(x[i]) : "") for i in eachindex(header)))
    end
    return rows
end

function parse_segments(s::AbstractString)
    segs = Segment[]
    for t in split(s, ';')
        m = match(r"^(\d+)\.\.(\d+):([+-])$", t)
        m === nothing && error("bad segment: $t")
        push!(segs, Segment(parse(Int, m[1]), parse(Int, m[2]), m[3][1]))
    end
    return segs
end

"Gene models of one genes.tsv in file order."
function read_models(path::AbstractString)
    models = Model[]
    for r in read_tsv_rows(path)
        flags = isempty(r["flags"]) ? String[] : String.(split(r["flags"], ','))
        parts = isempty(r["parts"]) ? Int[] : parse.(Int, split(r["parts"], ','))
        push!(models, Model(r["id"], r["gene"], r["kind"], r["strand"][1], parse_segments(r["segments"]),
                            parts, parse(Float64, r["p"]), r["pseudo"] == "1", flags))
    end
    return models
end

"Regions per record id from a regions.tsv of the assembler."
function read_regions(path::AbstractString)
    out = Dict{String,Regions}()
    acc = Dict{String,Vector{Tuple{String,Int,Int}}}()
    order = String[]
    for r in read_tsv_rows(path)
        id = r["id"]
        haskey(acc, id) || (acc[id] = Tuple{String,Int,Int}[]; push!(order, id))
        push!(acc[id], (r["region"], parse(Int, r["start"]), parse(Int, r["end"])))
    end
    for id in order
        v = acc[id]
        out[id] = Regions(maximum(t[3] for t in v), [t[1] for t in v], [(t[2], t[3]) for t in v])
    end
    return out
end

"Category table: rows (match, name, category, subgroup) in file order."
function read_categories(path::AbstractString)
    isfile(path) || error("categories file not found: $path")
    rows = read_tsv_rows(path)
    for r in rows
        r["match"] in ("exact", "kind", "prefix") || error("$path: unknown match type $(r["match"])")
    end
    return rows
end

# ------------------------------------------------------------------ geometry

seglen(s::Segment, L::Int) = mod(s.last - s.first, L) + 1

"Linear pieces of a segment (two when it crosses the origin)."
pieces(s::Segment, L::Int) = s.first <= s.last ? [(s.first, s.last)] : [(s.first, L), (1, s.last)]

overlap(a1, b1, a2, b2) = max(0, min(b1, b2) - max(a1, a2) + 1)

"Region names overlapped by the model, in region order (e.g. LSC,IR1 for rps12)."
function model_regions(m::Model, R::Regions)
    hit = falses(length(R.names))
    for s in m.segs, (a, b) in pieces(s, R.length), (i, (ra, rb)) in enumerate(R.spans)
        overlap(a, b, ra, rb) > 0 && (hit[i] = true)
    end
    return join(R.names[hit], ",")
end

"Distance along the strand from the 3' end of s1 to the 5' end of s2 (bases between them)."
function strand_gap(s1::Segment, s2::Segment, L::Int)
    s1.strand == '+' ? mod(s2.first - s1.last - 1, L) : mod(s1.first - s2.last - 1, L)
end

"""
Junction types between consecutive segments: :cis or :trans. A junction is trans when the
strand changes or the gap along the strand exceeds MAX_CIS_INTRON. A model flagged
trans_spliced by the annotator with no such junction gets its widest junction marked trans.
"""
function junctions(m::Model, L::Int)
    n = length(m.segs)
    n <= 1 && return Symbol[], Int[]
    gaps = [strand_gap(m.segs[k], m.segs[k+1], L) for k in 1:n-1]
    types = [(m.segs[k].strand != m.segs[k+1].strand || gaps[k] > MAX_CIS_INTRON) ? :trans : :cis for k in 1:n-1]
    if "trans_spliced" in m.flags && !(:trans in types)
        types[argmax(gaps)] = :trans
    end
    return types, gaps
end

incomplete(m::Model) = any(f -> startswith(f, "missing_part_"), m.flags)

# ------------------------------------------------------------------ categories and order

function categorize(gene::AbstractString, kind::AbstractString, cats)
    for r in cats
        r["match"] == "exact" && r["name"] == gene && return Category(r["category"], r["subgroup"])
    end
    for r in cats
        r["match"] == "kind" && r["name"] == kind && return Category(r["category"], r["subgroup"])
    end
    best = nothing; bestlen = 0
    for r in cats
        if r["match"] == "prefix" && startswith(gene, r["name"]) && length(r["name"]) > bestlen
            best = r; bestlen = length(r["name"])
        end
    end
    best === nothing && return Category("Unclassified", "Unclassified")
    return Category(best["category"], best["subgroup"])
end

"Category and subgroup order as first listed in the categories file; Unclassified last."
function category_order(cats)
    corder = String[]; sorder = Tuple{String,String}[]
    for r in cats
        r["category"] in corder || push!(corder, r["category"])
        (r["category"], r["subgroup"]) in sorder || push!(sorder, (r["category"], r["subgroup"]))
    end
    push!(corder, "Unclassified"); push!(sorder, ("Unclassified", "Unclassified"))
    return corder, sorder
end

"Natural order of gene names: rpl2 before rpl14, rrn4.5 before rrn5 before rrn16."
function gene_key(g::AbstractString)
    m = match(r"^([A-Za-z]+)(\d+(?:\.\d+)?)?(.*)$", g)
    m === nothing && return (g, -1.0, "")
    return (String(m[1]), m[2] === nothing ? -1.0 : parse(Float64, m[2]), String(m[3]))
end

# ------------------------------------------------------------------ per-gene summary

struct GeneSummary
    gene::String
    kind::String
    cat::Category
    functional::Vector{Model}     # non-pseudo models, by start
    pseudo::Vector{Model}
end

"Cis introns and trans-splicing of a gene: from its non-pseudo copy with most segments (complete copies first)."
function gene_structure(g::GeneSummary, L::Int)
    cands = [m for m in g.functional if !incomplete(m)]
    isempty(cands) && (cands = g.functional)
    isempty(cands) && return 0, false
    m = cands[argmax([length(x.segs) for x in cands])]
    types, _ = junctions(m, L)
    return count(==(:cis), types), :trans in types
end

function summarize(models::Vector{Model}, cats)
    bygene = Dict{String,Vector{Model}}()
    for m in models
        push!(get!(bygene, m.gene, Model[]), m)
    end
    out = GeneSummary[]
    for (gene, ms) in bygene
        sort!(ms; by=m -> minimum(s.first for s in m.segs))
        kind = ms[1].kind
        push!(out, GeneSummary(gene, kind, categorize(gene, kind, cats),
                               [m for m in ms if !m.pseudo], [m for m in ms if m.pseudo]))
    end
    return out
end

"Gene label for publication tables: psigene for pseudogene-only, intron and copy marks."
function gene_label(g::GeneSummary, L::Int)
    labels = String[]
    if !isempty(g.functional)
        nintr, trans = gene_structure(g, L)
        n = length(g.functional)
        push!(labels, g.gene * repeat("*", nintr) * (trans ? "#" : "") * (n > 1 ? "(x$n)" : ""))
    end
    if !isempty(g.pseudo)
        n = length(g.pseudo)
        push!(labels, "ψ" * g.gene * (n > 1 ? "(x$n)" : ""))
    end
    return join(labels, ", ")
end

# ------------------------------------------------------------------ gene-table

const LEGEND = [
    "# * one cis intron, ** two cis introns; # trans-spliced; (xN) N copies; ψ pseudogene (annotator call).",
    "# Copies and pseudogenes as annotated by PlastomeAnnotator; flags per model are in genes_long.tsv.",
]

safe_name(s::AbstractString) = replace(s, r"[^A-Za-z0-9._-]" => "_")

"""
Gene tables by function from the annotation of ordered plastomes.

Input: one or more <prefix>_plastome_genes.tsv; regions from
<prefix>_plastome_LSC_IR_SSC_IR.fasta.regions.tsv. Categories from a TSV
(gene_categories.tsv next to this file by default).

Outputs in `out`:
  genes_long.tsv          one row per gene model: category, region(s), segments, exon and
                          intron counts, trans-splicing, pseudo, annotator flags
  gene_table_<record>.tsv publication layout per record: category, subgroup, genes with marks
  gene_matrix.tsv         gene x record: functional copies, "+ψN" for pseudogene copies
  intron_genes.tsv        exon and intron lengths of genes with introns or trans-splicing
                          (first complete non-pseudo copy; "trans" for a trans junction)
  gene_summary.tsv        per record: length, region lengths, gene counts

Counts: a gene is counted once per name (tRNAs with their anticodon); "total" counts all
non-pseudo copies. Pseudogene models are reported separately and never counted as genes.
"""
function gene_table(files::Vector{String}, out::AbstractString; categories::AbstractString)
    cats = read_categories(categories)
    corder, sorder = category_order(cats)
    mkpath(out)
    records = String[]
    data = Dict{String,Tuple{Vector{GeneSummary},Regions,Vector{Model}}}()
    for f in files
        if !isfile(f)
            @warn "genes table not found; skipped" file=f
            continue
        end
        endswith(f, "_plastome_genes.tsv") || error("expected <prefix>_plastome_genes.tsv: $f")
        prefix = f[1:end-length("_plastome_genes.tsv")]
        rfile = prefix * "_plastome_LSC_IR_SSC_IR.fasta.regions.tsv"
        isfile(rfile) || error("regions file not found: $rfile")
        regions = read_regions(rfile)
        models = read_models(f)
        dups = duplicate_records(prefix * "_plastome_LSC_IR_SSC_IR.fasta")
        ids = unique(m.record for m in models)
        for id in ids
            if haskey(dups, id)
                @warn "record has the same sequence as an earlier one; skipped" record=id same_as=dups[id]
                continue
            end
            haskey(regions, id) || error("$id: no regions in $rfile")
            id in records && error("record $id given twice")
            ms = [m for m in models if m.record == id]
            push!(records, id)
            data[id] = (summarize(ms, cats), regions[id], ms)
        end
    end
    isempty(records) && error("no annotated records")

    unclassified = sort!(unique([g.gene for id in records for g in data[id][1] if g.cat.category == "Unclassified"]))
    isempty(unclassified) ||
        @warn "genes without a category (add them to the categories file)" genes=join(unclassified, ", ")

    sortkey(g) = (something(findfirst(==((g.cat.category, g.cat.subgroup)), sorder), length(sorder) + 1),
                  gene_key(g.gene))

    # genes_long.tsv
    open(joinpath(out, "genes_long.tsv"), "w") do io
        println(io, join(("record", "category", "subgroup", "gene", "kind", "copy", "strand", "regions",
                          "segments", "parts", "length", "exons", "cis_introns", "trans_spliced", "pseudo",
                          "incomplete", "p", "flags"), '\t'))
        for id in records
            gs, R, _ = data[id]
            for g in sort(gs; by=sortkey)
                for (ci, m) in enumerate(vcat(g.functional, g.pseudo))
                    types, _ = junctions(m, R.length)
                    segs = join(["$(s.first)..$(s.last):$(s.strand)" for s in m.segs], ";")
                    println(io, join((id, g.cat.category, g.cat.subgroup, g.gene, g.kind, ci, m.strand,
                                      model_regions(m, R), segs, join(m.parts, ","),
                                      sum(seglen(s, R.length) for s in m.segs), length(m.segs),
                                      count(==(:cis), types), (:trans in types) ? 1 : 0, m.pseudo ? 1 : 0,
                                      incomplete(m) ? 1 : 0, @sprintf("%.3f", m.p), join(m.flags, ",")), '\t'))
                end
            end
        end
    end

    # gene_table_<record>.tsv
    for id in records
        gs, R, _ = data[id]
        open(joinpath(out, "gene_table_$(safe_name(id)).tsv"), "w") do io
            foreach(l -> println(io, l), LEGEND)
            println(io, "category\tsubgroup\tgenes")
            for (c, s) in sorder
                sel = sort([g for g in gs if g.cat.category == c && g.cat.subgroup == s]; by=g -> gene_key(g.gene))
                isempty(sel) && continue
                println(io, join((c, s, join([gene_label(g, R.length) for g in sel], ", ")), '\t'))
            end
        end
    end

    # gene_matrix.tsv
    allgenes = Dict{String,GeneSummary}()
    for id in records, g in data[id][1]
        haskey(allgenes, g.gene) || (allgenes[g.gene] = g)
    end
    open(joinpath(out, "gene_matrix.tsv"), "w") do io
        println(io, join(vcat(["category", "subgroup", "gene", "kind"], records), '\t'))
        for g in sort(collect(values(allgenes)); by=sortkey)
            cells = String[]
            for id in records
                k = findfirst(x -> x.gene == g.gene, data[id][1])
                if k === nothing
                    push!(cells, "0")
                else
                    x = data[id][1][k]
                    nf, np = length(x.functional), length(x.pseudo)
                    push!(cells, np == 0 ? string(nf) : nf == 0 ? "ψ$np" : "$nf+ψ$np")
                end
            end
            println(io, join(vcat([g.cat.category, g.cat.subgroup, g.gene, g.kind], cells), '\t'))
        end
    end

    # intron_genes.tsv
    rows = Vector{Vector{String}}()
    maxseg = 1
    for id in records
        gs, R, _ = data[id]
        for g in sort(gs; by=sortkey)
            cands = [m for m in g.functional if !incomplete(m) && length(m.segs) > 1]
            isempty(cands) && continue
            m = cands[1]
            types, gaps = junctions(m, R.length)
            cols = String[]
            for k in eachindex(m.segs)
                push!(cols, string(seglen(m.segs[k], R.length)))
                k < length(m.segs) && push!(cols, types[k] == :trans ? "trans" : string(gaps[k]))
            end
            same = all(x -> [seglen(s, R.length) for s in x.segs] == [seglen(s, R.length) for s in m.segs] &&
                            junctions(x, R.length)[2] == gaps, cands)
            maxseg = max(maxseg, length(m.segs))
            push!(rows, vcat([id, g.gene, g.kind, model_regions(m, R), string(length(cands)), same ? "1" : "0"], cols))
        end
    end
    open(joinpath(out, "intron_genes.tsv"), "w") do io
        hdr = ["record", "gene", "kind", "regions", "complete_copies", "copies_same_structure"]
        roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII"]
        for k in 1:maxseg
            push!(hdr, "exon_" * (k <= length(roman) ? roman[k] : string(k)))
            k < maxseg && push!(hdr, "intron_" * (k <= length(roman) ? roman[k] : string(k)))
        end
        println(io, join(hdr, '\t'))
        for r in rows
            println(io, join(vcat(r, fill("", length(hdr) - length(r))), '\t'))
        end
    end

    # gene_summary.tsv
    open(joinpath(out, "gene_summary.tsv"), "w") do io
        println(io, join(("record", "length", "LSC", "IR", "SSC", "genes_unique", "genes_total", "CDS_unique",
                          "tRNA_unique", "rRNA_unique", "genes_multicopy", "genes_with_introns",
                          "trans_spliced", "incomplete_models", "pseudo_models", "unclassified"), '\t'))
        for id in records
            gs, R, ms = data[id]
            rl(name) = (k = findfirst(==(name), R.names); k === nothing ? "." : string(R.spans[k][2] - R.spans[k][1] + 1))
            fun = [g for g in gs if !isempty(g.functional)]
            st = [gene_structure(g, R.length) for g in fun]
            println(io, join((id, R.length, rl("LSC"), rl("IR1"), rl("SSC"), length(fun),
                              sum(length(g.functional) for g in fun; init=0),
                              count(g -> g.kind == "CDS", fun), count(g -> g.kind == "tRNA", fun),
                              count(g -> g.kind == "rRNA", fun), count(g -> length(g.functional) > 1, fun),
                              count(t -> t[1] > 0, st), count(t -> t[2], st),
                              count(m -> !m.pseudo && incomplete(m), ms), count(m -> m.pseudo, ms),
                              count(g -> g.cat.category == "Unclassified", gs)), '\t'))
        end
    end
    println("Records: ", length(records), "; tables written to ", out)
    return records
end

# ------------------------------------------------------------------ GenBank submission

# Closures that end the assembler's search (conclusive); other closures need force=1.
const CONCLUSIVE = Set(["validated", "ir_collapse", "read_bridge", "graph_path"])

# Columns of metadata.tsv that are not source modifiers.
const META_INFO = Set(["record", "seqid", "force", "closure", "length", "N_bases", "N_runs", "ready", "editing"])

const STARTS_11 = Set(["ATG", "GTG", "TTG", "CTG", "ATT", "ATC", "ATA"])   # translation table 11
const STOPS = Set(["TAA", "TAG", "TGA"])
const CU_STOPS = Set(["CAA", "CAG", "CGA"])      # become stops by C-to-U editing

const TRNA_AA = Dict("A" => "Ala", "C" => "Cys", "D" => "Asp", "E" => "Glu", "F" => "Phe", "G" => "Gly",
                     "H" => "His", "I" => "Ile", "K" => "Lys", "L" => "Leu", "M" => "Met", "fM" => "fMet",
                     "N" => "Asn", "P" => "Pro", "Q" => "Gln", "R" => "Arg", "S" => "Ser", "T" => "Thr",
                     "V" => "Val", "W" => "Trp", "Y" => "Tyr")
const RRNA_PRODUCT = Dict("rrn16" => "16S ribosomal RNA", "rrn23" => "23S ribosomal RNA",
                          "rrn5" => "5S ribosomal RNA", "rrn4.5" => "4.5S ribosomal RNA")

function read_fasta(path::AbstractString)
    recs = Tuple{String,String,String}[]          # (id, header, sequence)
    header = ""; buf = IOBuffer()
    flush!() = isempty(header) || push!(recs, (String(first(split(header))), header, uppercase(String(take!(buf)))))
    for line in eachline(path)
        if startswith(line, ">")
            flush!(); header = strip(line[2:end])
        else
            print(buf, strip(line))
        end
    end
    flush!()
    return recs
end

"""
Records of an ordered FASTA whose sequence equals that of an earlier record (the same molecule
written twice, e.g. two assembled SSC isomers ordered to the same orientation): id => first id.
"""
function duplicate_records(fasta::AbstractString)
    dup = Dict{String,String}(); seen = Dict{String,String}()
    isfile(fasta) || return dup
    for (id, _, seq) in read_fasta(fasta)
        haskey(seen, seq) ? (dup[id] = seen[seq]) : (seen[seq] = id)
    end
    return dup
end

"Value of key=value in a FASTA header, or empty."
function header_tag(header::AbstractString, key::AbstractString)
    m = match(Regex("(?:^|\\s)" * key * "=(\\S+)"), header)
    return m === nothing ? "" : String(m[1])
end

revcomp(s::AbstractString) = reverse(map(c -> c == 'A' ? 'T' : c == 'C' ? 'G' : c == 'G' ? 'C' : c == 'T' ? 'A' : 'N', s))

"Sequence of one segment along its strand (circular coordinates)."
function segment_seq(seq::AbstractString, s::Segment)
    L = length(seq)
    fw = s.first <= s.last ? seq[s.first:s.last] : seq[s.first:L] * seq[1:s.last]
    return s.strand == '+' ? fw : revcomp(fw)
end

spliced_seq(seq::AbstractString, m::Model) = join(segment_seq(seq, s) for s in m.segs)

"""
Translation checks of a spliced CDS with table 11 -> (status, reasons). status :ok, :editing
(passes only when C-to-U editing is assumed: ACG start, internal stops, CAA/CAG/CGA terminal
codon) or :fail (length not a multiple of 3, no terminal stop, editing needed but not on).
Codons with N are not counted as stops.
"""
function cds_check(s::AbstractString, editing::Bool)
    reasons = String[]
    length(s) % 3 == 0 || return :fail, ["length_not_multiple_of_3"]
    length(s) >= 6 || return :fail, ["too_short"]
    codons = [s[i:i+2] for i in 1:3:length(s)-2]
    need_edit = false
    if !(codons[1] in STARTS_11)
        if codons[1] == "ACG"
            push!(reasons, "start_ACG"); need_edit = true
        else
            return :fail, ["start_$(codons[1])"]
        end
    end
    internal = count(c -> c in STOPS, codons[1:end-1])
    if internal > 0
        push!(reasons, "internal_stops=$internal"); need_edit = true
    end
    if !(codons[end] in STOPS)
        if codons[end] in CU_STOPS
            push!(reasons, "terminal_$(codons[end])"); need_edit = true
        else
            return :fail, vcat(reasons, ["no_terminal_stop"])
        end
    end
    need_edit || return :ok, reasons
    return editing ? (:editing, reasons) : (:fail, vcat(reasons, ["editing_not_supported_for_record"]))
end

"5-column interval lines of a segment, splitting at the origin (minus strand: start > stop)."
function tbl_intervals(s::Segment, L::Int)
    if s.first <= s.last
        return s.strand == '+' ? [(s.first, s.last)] : [(s.last, s.first)]
    end
    return s.strand == '+' ? [(s.first, L), (1, s.last)] : [(s.last, 1), (L, s.first)]
end

"""
Gene feature segments: one span per group of cis-joined segments (a trans-spliced gene has one
span per RNA molecule, e.g. 5'-rps12 and 3'-rps12).
"""
function gene_spans(m::Model, L::Int)
    types, _ = junctions(m, L)
    groups = Vector{Vector{Segment}}([[m.segs[1]]])
    for k in 2:length(m.segs)
        types[k-1] == :trans ? push!(groups, [m.segs[k]]) : push!(groups[end], m.segs[k])
    end
    spans = Segment[]
    for g in groups
        st = g[1].strand
        push!(spans, st == '+' ? Segment(g[1].first, g[end].last, st) : Segment(g[end].first, g[1].last, st))
    end
    return spans
end

function trna_product(gene::AbstractString)
    m = match(r"^trn(fM|[A-Z])", gene)
    m === nothing && return ""
    aa = get(TRNA_AA, m[1], "")
    return isempty(aa) ? "" : "tRNA-" * aa
end

"CDS product: products table, then rps/rpl rules; empty when unknown."
function cds_product(gene::AbstractString, products::Dict{String,String})
    haskey(products, gene) && return products[gene]
    m = match(r"^rps(\d+)$", gene); m === nothing || return "ribosomal protein S$(m[1])"
    m = match(r"^rpl(\d+)$", gene); m === nothing || return "ribosomal protein L$(m[1])"
    return ""
end

function read_products(path::AbstractString)
    isfile(path) || error("products file not found: $path")
    return Dict(r["gene"] => r["product"] for r in read_tsv_rows(path))
end

"Editing tolerance per record from the editing column of genes.tsv (missing column: nothing)."
function read_editing(path::AbstractString)
    out = Dict{String,Union{Bool,Nothing}}()
    for r in read_tsv_rows(path)
        out[r["id"]] = haskey(r, "editing") ? r["editing"] == "1" : nothing
    end
    return out
end

"N runs as (start, stop)."
function n_runs(seq::AbstractString)
    runs = Tuple{Int,Int}[]
    i = 1; n = length(seq)
    while i <= n
        if seq[i] == 'N'
            j = i
            while j < n && seq[j+1] == 'N'
                j += 1
            end
            push!(runs, (i, j)); i = j + 1
        else
            i += 1
        end
    end
    return runs
end

struct SubmissionRecord
    id::String
    header::String
    seq::String
    closure::String
    regions::Regions
    models::Vector{Model}
    editing::Union{Bool,Nothing}
end

"Ordered records with annotation from genes tables (FASTA, regions and genes of each prefix)."
function load_submission_records(files::Vector{String})
    out = SubmissionRecord[]
    for f in files
        endswith(f, "_plastome_genes.tsv") || error("expected <prefix>_plastome_genes.tsv: $f")
        isfile(f) || (@warn "genes table not found; skipped" file=f; continue)
        prefix = f[1:end-length("_plastome_genes.tsv")]
        fasta = prefix * "_plastome_LSC_IR_SSC_IR.fasta"
        regions = read_regions(fasta * ".regions.tsv")
        models = read_models(f)
        editing = read_editing(f)
        dups = duplicate_records(fasta)
        for (id, header, seq) in read_fasta(fasta)
            if haskey(dups, id)
                @warn "record has the same sequence as an earlier one; skipped" record=id same_as=dups[id]
                continue
            end
            haskey(regions, id) || error("$id: no regions in $fasta.regions.tsv")
            regions[id].length == length(seq) || error("$id: regions length $(regions[id].length) != sequence $(length(seq))")
            push!(out, SubmissionRecord(id, header, seq, header_tag(header, "closure"), regions[id],
                                        [m for m in models if m.record == id], get(editing, id, nothing)))
        end
    end
    return out
end

const TEMPLATE_MODIFIERS = ["organism", "isolate", "specimen_voucher", "geo_loc_name", "collection_date",
                            "collected_by", "location"]

"""
Metadata template: one row per ordered record with information columns (closure, length, N
bases, ready = conclusive closure without N) and empty source-modifier columns. location
defaults to chloroplast. Rows without organism are not exported; force=1 exports a record that
is not ready.
"""
function write_template(files::Vector{String}, path::AbstractString)
    isfile(path) && error("$path exists; not overwritten")
    recs = load_submission_records(files)
    open(path, "w") do io
        println(io, join(vcat(["record", "seqid", "closure", "length", "N_bases", "N_runs", "editing", "ready", "force"],
                              TEMPLATE_MODIFIERS), '\t'))
        for r in recs
            nb = count(==('N'), r.seq)
            ready = r.closure in CONCLUSIVE && nb == 0
            ed = r.editing === nothing ? "." : r.editing ? "1" : "0"
            println(io, join(vcat([r.id, r.id, isempty(r.closure) ? "." : r.closure, string(length(r.seq)), string(nb),
                                   string(length(n_runs(r.seq))), ed, ready ? "1" : "0", "0"],
                                  [m == "location" ? "chloroplast" : "" for m in TEMPLATE_MODIFIERS]), '\t'))
        end
    end
    println("Template with ", length(recs), " records written to ", path)
end

function tbl_feature(io, intervals, key, quals)
    for (k, (a, b)) in enumerate(intervals)
        k == 1 ? println(io, a, '\t', b, '\t', key) : println(io, a, '\t', b)
    end
    for (q, v) in quals
        println(io, "\t\t\t", q, '\t', v)
    end
end

"""
GenBank submission files from annotated ordered plastomes.

Records exported: rows of the metadata table with organism; a record that is not ready
(closure not conclusive, or N in the sequence) only with force=1. Outputs in `out`:
  submission.fsa           sequences, table2asn source modifiers in the deflines
                           ([organism=...] [location=...] [gcode=11], and [topology=circular]
                           [completeness=complete] for records without assembly gaps)
  submission.tbl           5-column feature table: gene, CDS, tRNA, rRNA, repeat_region (IRs),
                           assembly_gap (N runs >= gap_min)
  source_modifiers.tsv     the same source modifiers as a table for the Submission Portal
  report.tsv               every gene model with the action taken and the reason
  report_records.tsv       per record: exported or skipped, counts of each action

Gene models:
  CDS passing translation table 11 checks              -> gene + CDS
  CDS passing only with C-to-U editing, record editing=1 -> gene + CDS /exception="RNA editing"
  CDS failing the checks (also ycf15/ycf68 kept by the annotator) -> gene /pseudo, no CDS
  annotator pseudo                                      -> gene /pseudo
  model with a missing part (missing_part_k)            -> gene /pseudo + note, to be reviewed
  CDS overlapping an assembly gap                       -> gene only, to be reviewed
  trans-spliced: gene with one span per molecule, gene and CDS /trans_splicing
Every feature carries its /gene explicitly (genes nested in introns, e.g. matK in trnK-UUU).
"""
function genbank(files::Vector{String}, out::AbstractString; metadata::AbstractString,
                 products::AbstractString, gap_min::Int=10)
    meta = Dict{String,Dict{String,String}}()
    for r in read_tsv_rows(metadata)
        haskey(meta, r["record"]) && error("$metadata: record $(r["record"]) listed twice")
        meta[r["record"]] = r
    end
    prod = read_products(products)
    recs = load_submission_records(files)
    mkpath(out)
    nexp = 0
    seqids = Set{String}()
    open(joinpath(out, "submission.fsa"), "w") do fsa
    open(joinpath(out, "submission.tbl"), "w") do tbl
    open(joinpath(out, "source_modifiers.tsv"), "w") do src
    open(joinpath(out, "report.tsv"), "w") do rep
    open(joinpath(out, "report_records.tsv"), "w") do rrep
        println(rep, join(("record", "seqid", "gene", "kind", "copy", "segments", "action", "reasons", "annotator_flags"), '\t'))
        println(rrep, join(("record", "seqid", "status", "closure", "length", "assembly_gaps", "ambiguous_N",
                            "genes", "CDS", "CDS_RNA_editing", "pseudo", "review"), '\t'))
        srccols = String[]
        for m in values(meta), k in keys(m)
            (k in META_INFO || k in srccols || isempty(m[k])) || push!(srccols, k)
        end
        sort!(srccols; by=k -> (k != "organism", k))
        println(src, join(vcat(["Sequence_ID"], srccols), '\t'))
        for r in recs
            m = get(meta, r.id, nothing)
            if m === nothing || isempty(get(m, "organism", ""))
                println(rrep, join((r.id, ".", "skipped: no metadata/organism", r.closure, length(r.seq), ".", ".", ".", ".", ".", ".", "."), '\t'))
                continue
            end
            seqid = isempty(get(m, "seqid", "")) ? r.id : m["seqid"]
            occursin(r"^[A-Za-z0-9_.:*#-]+$", seqid) || error("$seqid: sequence ID with unsupported characters")
            seqid in seqids && error("sequence ID $seqid used twice")
            runs = n_runs(r.seq)
            gaps = [x for x in runs if x[2] - x[1] + 1 >= gap_min]
            ambiguous = sum((x[2] - x[1] + 1 for x in runs if x[2] - x[1] + 1 < gap_min); init=0)
            ready = r.closure in CONCLUSIVE && isempty(runs)
            if !ready && get(m, "force", "0") != "1"
                println(rrep, join((r.id, seqid, "skipped: not ready (closure $(r.closure), N runs $(length(runs))); force=1 to export",
                                    r.closure, length(r.seq), length(gaps), ambiguous, ".", ".", ".", ".", "."), '\t'))
                continue
            end
            if r.editing === nothing
                @warn "genes table has no editing column (assembler < 16.47.1): CDS needing editing are exported as pseudo" record=r.id
            end
            editing = r.editing === true
            push!(seqids, seqid)
            L = length(r.seq)
            complete = isempty(gaps)

            # defline and source table
            mods = [(k, m[k]) for k in srccols if !isempty(get(m, k, ""))]
            any(t -> t[1] == "location", mods) || push!(mods, ("location", "chloroplast"))
            defl = join(["[$k=$v]" for (k, v) in mods], " ") * " [gcode=11]" *
                   (complete ? " [topology=circular] [completeness=complete]" : "")
            println(fsa, ">", seqid, " ", defl)
            for i in 1:80:L
                println(fsa, r.seq[i:min(i + 79, L)])
            end
            println(src, join(vcat([seqid], [get(m, k, "") for k in srccols]), '\t'))

            # features
            println(tbl, ">Feature ", seqid)
            gapmask = falses(L)
            for (a, b) in gaps
                gapmask[a:b] .= true
            end
            ncds = nedit = npseudo = nreview = 0
            ms = sort(r.models; by=x -> minimum(s.first for s in x.segs))
            seen = Dict{String,Int}()
            for x in ms
                seen[x.gene] = get(seen, x.gene, 0) + 1
                types, _ = junctions(x, L)
                trans = :trans in types
                segtxt = join(["$(s.first)..$(s.last):$(s.strand)" for s in x.segs], ";")
                ongap = any(s -> any(p -> any(gapmask[p[1]:p[2]]), pieces(s, L)), x.segs)
                action = ""; reasons = String[]
                cds_quals = Tuple{String,String}[]
                if x.pseudo
                    action = "pseudo_annotator"
                elseif incomplete(x)
                    action = "pseudo_incomplete_review"; push!(reasons, join(filter(f -> startswith(f, "missing_part_"), x.flags), ","))
                elseif ongap
                    action = "gene_only_gap_review"
                elseif x.kind == "CDS"
                    st, why = cds_check(spliced_seq(r.seq, x), editing)
                    append!(reasons, why)
                    product = cds_product(x.gene, prod)
                    if st == :fail
                        action = "pseudo_orf_check"
                    elseif isempty(product)
                        action = "gene_only_no_product_review"
                    else
                        action = st == :editing ? "CDS_RNA_editing" : "CDS"
                        push!(cds_quals, ("gene", x.gene), ("product", product), ("transl_table", "11"))
                        st == :editing && push!(cds_quals, ("exception", "RNA editing"))
                        trans && push!(cds_quals, ("trans_splicing", ""))
                    end
                else
                    product = x.kind == "tRNA" ? trna_product(x.gene) : get(RRNA_PRODUCT, x.gene, "")
                    if isempty(product)
                        action = "gene_only_no_product_review"
                    else
                        action = x.kind
                        push!(cds_quals, ("gene", x.gene), ("product", product))
                        trans && push!(cds_quals, ("trans_splicing", ""))
                    end
                end
                ispseudo = startswith(action, "pseudo")
                gq = Tuple{String,String}[("gene", x.gene)]
                ispseudo && push!(gq, ("pseudo", ""))
                action == "pseudo_incomplete_review" && push!(gq, ("note", "incomplete gene model"))
                trans && push!(gq, ("trans_splicing", ""))
                tbl_feature(tbl, vcat([tbl_intervals(s, L) for s in gene_spans(x, L)]...), "gene", gq)
                if !isempty(cds_quals)
                    key = x.kind == "CDS" ? "CDS" : x.kind
                    tbl_feature(tbl, vcat([tbl_intervals(s, L) for s in x.segs]...), key, cds_quals)
                end
                ncds += action in ("CDS", "CDS_RNA_editing")
                nedit += action == "CDS_RNA_editing"
                npseudo += ispseudo
                nreview += endswith(action, "_review")
                println(rep, join((r.id, seqid, x.gene, x.kind, seen[x.gene], segtxt, action, join(reasons, ","),
                                   join(x.flags, ",")), '\t'))
            end
            for (name, label) in (("IR1", "IRa"), ("IR2", "IRb"))
                k = findfirst(==(name), r.regions.names)
                k === nothing && continue
                a, b = r.regions.spans[k]
                tbl_feature(tbl, [(a, b)], "repeat_region",
                            [("rpt_type", "inverted"), ("note", "inverted repeat $label")])
            end
            linkage = r.closure == "ref_scaffold" ? "align genus" : r.closure == "inferred_join" ? "unspecified" : "paired-ends"
            for (a, b) in gaps
                tbl_feature(tbl, [(a, b)], "assembly_gap",
                            [("estimated_length", "unknown"), ("gap_type", "within scaffold"), ("linkage_evidence", linkage)])
            end
            println(rrep, join((r.id, seqid, complete ? "exported" : "exported (force=1, assembly gaps)", r.closure, L,
                                length(gaps), ambiguous, length(ms), ncds, nedit, npseudo, nreview), '\t'))
            nexp += 1
        end
    end; end; end; end; end
    println("Exported ", nexp, " of ", length(recs), " records to ", out, " (see report_records.tsv)")
    return nexp
end

# ------------------------------------------------------------------ circular genome map (SVG)

const MAP_DEFAULT_COLOR = "#8a8a8a"

xmlesc(s) = replace(String(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
fmt_bp(x::Integer) = replace(string(x), r"(?<=\d)(?=(\d{3})+$)" => ",")

struct MapGeom
    cx::Float64
    cy::Float64
    L::Int
    rot::Float64      # added to every angle (radians): puts the SSC at 6 o'clock
end

"Angle (radians, clockwise from 12 o'clock) of the start of base x; coordinates run counterclockwise (OGDRAW)."
map_angle(g::MapGeom, x::Real) = -2π * (x - 1) / g.L + g.rot
polar(g::MapGeom, t::Real, r::Real) = (g.cx + r * sin(t), g.cy - r * cos(t))

"SVG path of an annulus sector from angle t1 to t2 (clockwise) between radii r1 < r2."
function sector_path(g::MapGeom, t1::Real, t2::Real, r1::Real, r2::Real)
    t2 = max(t2, t1 + 0.0015)                 # very short features stay visible
    large = t2 - t1 > π ? 1 : 0
    x1, y1 = polar(g, t1, r2); x2, y2 = polar(g, t2, r2)
    x3, y3 = polar(g, t2, r1); x4, y4 = polar(g, t1, r1)
    return @sprintf("M%.2f,%.2f A%.2f,%.2f 0 %d 1 %.2f,%.2f L%.2f,%.2f A%.2f,%.2f 0 %d 0 %.2f,%.2f Z",
                    x1, y1, r2, r2, large, x2, y2, x3, y3, r1, r1, large, x4, y4)
end

"Open arc at radius r from t1 to t2 (clockwise)."
function arc_line(g::MapGeom, t1::Real, t2::Real, r::Real)
    t2 < t1 && (t2 += 2π)
    x1, y1 = polar(g, t1, r); x2, y2 = polar(g, t2, r)
    return @sprintf("M%.2f,%.2f A%.2f,%.2f 0 %d 1 %.2f,%.2f", x1, y1, r, r, t2 - t1 > π ? 1 : 0, x2, y2)
end

"Angular interval(s) of a segment (split at the origin)."
seg_angles(g::MapGeom, s::Segment) = [(map_angle(g, b + 1), map_angle(g, a)) for (a, b) in pieces(s, g.L)]   # clockwise from..to

"""
Label angles pushed apart until neighbours are at least `delta` radians apart (circular);
labels keep their order. Returns angles in the input order.
"""
function spread_angles(theta::Vector{Float64}, delta::Float64; iters::Int=3000)
    n = length(theta)
    n <= 1 && return copy(theta)
    order = sortperm(theta)
    a = theta[order]
    delta = min(delta, 2π / n)
    for _ in 1:iters
        moved = false
        for i in 1:n
            j = i == n ? 1 : i + 1
            gap = i == n ? a[1] + 2π - a[n] : a[j] - a[i]
            if gap < delta - 1e-9
                d = (delta - gap) / 2
                a[i] -= d; a[j] += d
                moved = true
            end
        end
        moved || break
    end
    out = similar(theta)
    for k in 1:n
        out[order[k]] = a[k]
    end
    return out
end

"Radial text at angle t: reads outward (outward=true) or inward from radius r, never upside down."
function radial_text(io, g::MapGeom, t::Real, r::Real, s::AbstractString; outward::Bool, size::Real=9,
                     style::AbstractString="italic", fill::AbstractString="#222")
    x, y = polar(g, t, r)
    right = mod(t, 2π) < π
    rot = rad2deg(t) + (right ? -90 : 90)
    anchor = outward == right ? "start" : "end"
    println(io, @sprintf("""<text x="%.2f" y="%.2f" transform="rotate(%.2f %.2f %.2f)" font-size="%g" font-style="%s" fill="%s" text-anchor="%s" dominant-baseline="central">%s</text>""",
                         x, y, rot, x, y, size, style, fill, anchor, xmlesc(s)))
end

"""
Sliding GC fraction on the circular sequence: a window of `window` bases centred every `step`
bases (centres 1, 1+step, ...), wrapping around the origin; N and other symbols are ignored.
Returns (centre, GC fraction); windows without A/C/G/T are left out.
"""
function gc_sliding(seq::AbstractString, window::Int, step::Int)
    L = length(seq)
    (window >= 1 && step >= 1) || error("GC window and step must be positive")
    w = min(window, L)
    cg = zeros(Int, 2L + 1); cat = zeros(Int, 2L + 1)       # prefix sums over the doubled sequence
    for i in 1:2L
        c = seq[mod1(i, L)]
        cg[i+1] = cg[i] + (c == 'G' || c == 'C')
        cat[i+1] = cat[i] + (c == 'A' || c == 'T')
    end
    out = Tuple{Int,Float64}[]
    for centre in 1:step:L
        a = mod1(centre - w ÷ 2, L)                      # window a .. a+w-1 on the doubled sequence
        gc = cg[a+w] - cg[a]; at = cat[a+w] - cat[a]
        gc + at > 0 && push!(out, (centre, gc / (gc + at)))
    end
    return out
end

function gc_total(seq::AbstractString)
    gc = count(c -> c == 'G' || c == 'C', seq); at = count(c -> c == 'A' || c == 'T', seq)
    return gc / max(1, gc + at)
end

"""
Circular map of one annotated ordered plastome as SVG (OGDRAW-like layout):
  genes on the outer ring coloured by functional subgroup (colour column of the categories
  file); genes on the + strand outside the ring (transcribed counterclockwise), - strand inside
  (clockwise); exons as blocks, cis introns as thin arcs; pseudogenes pale, labelled ψ;
  quadripartite structure (LSC, IRa, SSC, IRb from regions.tsv; IRa follows the LSC) on the inner ring with lengths;
  GC content as deviation from the genome mean in a sliding window (`window` bases, every `step`
  bases; dark above, light below the mean);
  organism, length and GC in the centre; legend of the subgroups present on the right.
Coordinates increase counterclockwise, as in OGDRAW (LSC at the top, IRa on the left, SSC at the
bottom, IRb on the right); the map is rotated so that the middle of the SSC is at 6 o'clock
(position 1, the IRb/LSC junction of the ordered FASTA, is marked by its IR boundary line).
"""
function plot_genome_map(r::SubmissionRecord, file::AbstractString, cats; title::AbstractString="",
                         window::Int=500, step::Int=50)
    W, H = 1400, 1000
    L = length(r.seq)
    k = findfirst(==("SSC"), r.regions.names)
    rot = 0.0
    if k !== nothing
        a, b = r.regions.spans[k]
        rot = π + 2π * ((a + b + 1) / 2 - 1) / L          # SSC centre to 6 o'clock
    end
    g = MapGeom(500.0, 500.0, L, rot)
    R, h = 330.0, 20.0                        # gene ring: + strand R..R+h, - strand R-h..R
    rreg1, rreg2 = 222.0, 232.0               # region ring
    rgc, amp = 150.0, 30.0                    # GC baseline and maximal bar length
    colors = Dict{String,String}()
    for row in cats
        c = get(row, "color", "")
        isempty(c) || haskey(colors, row["subgroup"]) || (colors[row["subgroup"]] = c)
    end
    _, sorder = category_order(cats)
    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$W" height="$H" viewBox="0 0 $W $H" font-family="Helvetica, Arial, sans-serif">""")
    println(io, """<rect width="$W" height="$H" fill="white"/>""")

    # IR boundary lines (under everything else)
    R2 = r.regions
    for (name, (a, b)) in zip(R2.names, R2.spans)
        startswith(name, "IR") || continue
        for x in (a, b + 1)
            t = map_angle(g, x)
            x1, y1 = polar(g, t, rreg1); x2, y2 = polar(g, t, R + h)
            println(io, @sprintf("""<line x1="%.2f" y1="%.2f" x2="%.2f" y2="%.2f" stroke="#bbbbbb" stroke-width="0.8" stroke-dasharray="3,3"/>""", x1, y1, x2, y2))
        end
    end

    # GC content
    wins = gc_sliding(r.seq, window, step)
    gmean = gc_total(r.seq)
    maxdev = maximum(abs(w[2] - gmean) for w in wins; init=0.0)
    maxdev = max(maxdev, 1e-6)
    for (c, f) in wins                         # bar of one step centred on the window centre
        d = (f - gmean) / maxdev * amp
        abs(d) < 0.05 && continue
        r1, r2 = d >= 0 ? (rgc, rgc + d) : (rgc + d, rgc)
        println(io, """<path d="$(sector_path(g, map_angle(g, c + step / 2), map_angle(g, c - step / 2), r1, r2))" fill="$(d >= 0 ? "#4d4d4d" : "#b3b3b3")"/>""")
    end
    println(io, @sprintf("""<circle cx="%.2f" cy="%.2f" r="%.2f" fill="none" stroke="#666" stroke-width="0.6"/>""", g.cx, g.cy, rgc))

    # quadripartite structure
    for (name, (a, b)) in zip(R2.names, R2.spans)
        label = name == "IR1" ? "IRa" : name == "IR2" ? "IRb" : name
        fill = startswith(name, "IR") ? "#6b6b6b" : name == "LSC" ? "#d9d9d9" : "#a8a8a8"
        println(io, """<path d="$(sector_path(g, map_angle(g, b + 1), map_angle(g, a), rreg1, rreg2))" fill="$fill" stroke="white" stroke-width="0.5"/>""")
        # name and length tangential to the ring, inside it; turned 180 degrees in the lower half
        t = mod(map_angle(g, (a + b + 1) / 2), 2π)
        flip = π / 2 < t < 3π / 2
        rot = rad2deg(t) - (flip ? 180 : 0)
        rname, rlen = flip ? (rreg1 - 20, rreg1 - 9) : (rreg1 - 9, rreg1 - 20)
        for (rr, txt, attrs) in ((rname, label, "font-size=\"11\" font-weight=\"bold\" fill=\"#333\""),
                                 (rlen, fmt_bp(b - a + 1) * " bp", "font-size=\"9\" fill=\"#555\""))
            x, y = polar(g, t, rr)
            println(io, @sprintf("""<text x="%.2f" y="%.2f" transform="rotate(%.2f %.2f %.2f)" text-anchor="middle" dominant-baseline="central" %s>%s</text>""",
                                 x, y, rot, x, y, attrs, txt))
        end
    end

    # genes
    ms = sort(r.models; by=m -> minimum(s.first for s in m.segs))
    labels_in = Tuple{Float64,String,Bool}[]     # (gene angle, text, pseudo) for + strand
    labels_out = Tuple{Float64,String,Bool}[]
    present = Set{String}()
    seen_labels = Set{Tuple{Float64,String}}()
    pseudo_present = false
    println(io, """<circle cx="$(g.cx)" cy="$(g.cy)" r="$R" fill="none" stroke="#444" stroke-width="0.8"/>""")
    for m in ms
        cat = categorize(m.gene, m.kind, cats)
        push!(present, cat.subgroup)
        col = get(colors, cat.subgroup, MAP_DEFAULT_COLOR)
        ring(st) = st == '+' ? (R, R + h) : (R - h, R)     # + outside, - inside; per segment (trans-spliced parts differ)
        op = m.pseudo ? "0.35" : "1"
        m.pseudo && (pseudo_present = true)
        types, _ = junctions(m, g.L)
        for (k, s) in enumerate(m.segs)
            r1, r2 = ring(s.strand)
            for (t1, t2) in seg_angles(g, s)
                println(io, """<path d="$(sector_path(g, t1, t2, r1, r2))" fill="$col" fill-opacity="$op" stroke="#333" stroke-width="0.3"/>""")
            end
            if k < length(m.segs) && types[k] == :cis
                s2 = m.segs[k+1]
                t1, t2 = s.strand == '+' ? (map_angle(g, s2.first), map_angle(g, s.last + 1)) : (map_angle(g, s.first), map_angle(g, s2.last + 1))
                println(io, """<path d="$(arc_line(g, t1, t2, (r1 + r2) / 2))" fill="none" stroke="#333" stroke-width="0.7"/>""")
            end
        end
        for sp in gene_spans(m, g.L)
            len = mod(sp.last - sp.first, g.L) + 1
            t = mod(map_angle(g, sp.first + len / 2), 2π)
            text = (m.pseudo ? "ψ" : "") * m.gene
            key = (round(t; digits=4), text)          # shared 5' exon of trans-spliced copies: one label
            key in seen_labels && continue
            push!(seen_labels, key)
            push!(sp.strand == '+' ? labels_out : labels_in, (t, text, m.pseudo))
        end
    end

    # labels with leader lines
    for (set, outward) in ((labels_out, true), (labels_in, false))
        isempty(set) && continue
        redge = outward ? R + h : R - h
        rl = outward ? redge + 18 : redge - 18
        ts = spread_angles([x[1] for x in set], 11.0 / rl)
        sgn = outward ? 1 : -1
        for (k, (t0, text, ps)) in enumerate(set)
            x1, y1 = polar(g, t0, redge); x2, y2 = polar(g, t0, redge + sgn * 6)
            x3, y3 = polar(g, ts[k], redge + sgn * 15)
            println(io, @sprintf("""<polyline points="%.2f,%.2f %.2f,%.2f %.2f,%.2f" fill="none" stroke="#999" stroke-width="0.5"/>""", x1, y1, x2, y2, x3, y3))
            radial_text(io, g, ts[k], rl, text; outward=outward, fill=ps ? "#777" : "#222")
        end
    end

    # centre
    gcpct = @sprintf("%.1f", 100 * gmean)
    ttl = isempty(title) ? r.id : title
    println(io, """<text x="$(g.cx)" y="$(g.cy - 22)" font-size="18" font-style="$(isempty(title) ? "normal" : "italic")" text-anchor="middle">$(xmlesc(ttl))</text>""")
    println(io, """<text x="$(g.cx)" y="$(g.cy + 2)" font-size="13" text-anchor="middle" fill="#333">plastid genome</text>""")
    println(io, """<text x="$(g.cx)" y="$(g.cy + 22)" font-size="13" text-anchor="middle" fill="#333">$(fmt_bp(g.L)) bp, GC $gcpct%</text>""")

    # legend
    x0, y = 1000, 70
    println(io, """<text x="$x0" y="$y" font-size="13" font-weight="bold">Gene functions</text>""")
    y += 10
    for (c, s) in sorder
        s in present || continue
        y += 18
        println(io, """<rect x="$x0" y="$(y - 10)" width="14" height="12" fill="$(get(colors, s, MAP_DEFAULT_COLOR))" stroke="#333" stroke-width="0.3"/>""")
        println(io, """<text x="$(x0 + 22)" y="$y" font-size="11">$(xmlesc(s))</text>""")
    end
    if pseudo_present
        y += 18
        println(io, """<rect x="$x0" y="$(y - 10)" width="14" height="12" fill="#8a8a8a" fill-opacity="0.35" stroke="#333" stroke-width="0.3"/>""")
        println(io, """<text x="$(x0 + 22)" y="$y" font-size="11">pseudogene (ψ)</text>""")
    end
    y += 30
    for (fill, text) in (("#4d4d4d", "GC above genome mean"), ("#b3b3b3", "GC below genome mean"))
        println(io, """<rect x="$x0" y="$(y - 10)" width="14" height="12" fill="$fill"/>""")
        println(io, """<text x="$(x0 + 22)" y="$y" font-size="11">$text</text>""")
        y += 18
    end
    println(io, """<text x="$x0" y="$y" font-size="10" fill="#555">GC in $(fmt_bp(window)) bp window, step $(fmt_bp(step)) bp; mean $gcpct%</text>""")
    y += 26
    println(io, """<text x="$x0" y="$y" font-size="10" fill="#555">Genes outside the circle: + strand (counterclockwise);</text>""")
    println(io, """<text x="$x0" y="$(y + 13)" font-size="10" fill="#555">inside: &#8722; strand (clockwise).</text>""")
    println(io, "</svg>")
    write(file, take!(io))
    return file
end

"""
Circular maps (SVG) of annotated ordered plastomes: one map_<record>.svg per record in `out`.
Titles: organism from a metadata table (record, organism; the genbank template works), else
the record id.
"""
function genome_maps(files::Vector{String}, out::AbstractString; categories::AbstractString,
                     metadata::AbstractString="", window::Int=500, step::Int=50)
    cats = read_categories(categories)
    titles = Dict{String,String}()
    if !isempty(metadata)
        for m in read_tsv_rows(metadata)
            isempty(get(m, "organism", "")) || (titles[m["record"]] = m["organism"])
        end
    end
    mkpath(out)
    n = 0
    for r in load_submission_records(files)
        isempty(r.models) && (@warn "no gene models; map skipped" record=r.id; continue)
        f = joinpath(out, "map_$(safe_name(r.id)).svg")
        plot_genome_map(r, f, cats; title=get(titles, r.id, ""), window=window, step=step)
        n += 1
    end
    println("Maps: ", n, " written to ", out)
    return n
end

# ------------------------------------------------------------------ IR junction plot (SVG)

const ANNOTATOR = Ref{Any}(nothing)

"Loads PlastomeAnnotator.jl from `dir` once (its PlastomeIR module draws the junction plot)."
function load_annotator_module(dir::AbstractString)
    ANNOTATOR[] === nothing || return ANNOTATOR[]
    path = joinpath(dir, "PlastomeAnnotator.jl")
    isfile(path) || error("PlastomeAnnotator.jl not found in $dir")
    Base.include(@__MODULE__, path)
    ANNOTATOR[] = Base.invokelatest(getfield, @__MODULE__, :PlastomeAnnotator)
    return ANNOTATOR[]
end

# IR names. The annotator follows IRscope (LSC-IRb-SSC-IRa; junctions JLB, JSB, JSA, JLA); our
# outputs name the IR that follows the LSC IRa (LSC-IRa-SSC-IRb; junctions JLA, JSA, JSB, JLB).
# The annotator's first copy (its "IRb", fields irb / pos_irb) is therefore our IRa.
copy_name(c) = c === :IRb ? "IRa" : c === :IRa ? "IRb" : string(c)

"Renames IR copies and junctions in the text of a plot_ir_junctions SVG to LSC-IRa-SSC-IRb."
function relabel_ir_svg(svg::AbstractString)
    swap = Dict("IRb" => "IRa", "IRa" => "IRb", "JLB" => "JLA", "JLA" => "JLB", "JSB" => "JSA", "JSA" => "JSB")
    txt = read(svg, String)
    txt = replace(txt, r"(<text[^>]*>)([^<]*)(</text>)" => function (m)
        mm = match(r"^(<text[^>]*>)([^<]*)(</text>)$", m)
        body = replace(mm.captures[2], r"\b(IRa|IRb|JLA|JLB|JSA|JSB)\b" => k -> swap[k])
        return mm.captures[1] * body * mm.captures[3]
    end)
    write(svg, txt)
    return svg
end

"""
IRscope-style junction plot of assembled plastomes (and optional dataset references) with
plot_ir_junctions of PlastomeAnnotator. One row per record, in the order of the genes
tables, then the references (`refs`, accessions of the annotator dataset). Per record:
sequence of the ordered FASTA, IR from the annotator's detector (find_ir_status, which also
gives the differences between the IR copies), genes from <prefix>_plastome_genes.tsv
(pseudogenes named ψgene). Row names: organism from the metadata table, else the record id.
Records with the same sequence as an earlier record are skipped.

The IR boundaries drawn come from find_ir_status; they are compared with the assembler's
regions.tsv (IR1 = IRa, IR2 = IRb) in <svg without .svg>_ir_check.tsv, and every difference
is also printed as a warning. The differences between the IR copies found by the annotator
(SNPs, indels; positions in each copy) are listed in <svg without .svg>_ir_diffs.tsv, with
in_assembler_IR = yes when both positions lie inside the assembler's IRs.
"""
function ir_plot(files::Vector{String}, svg::AbstractString; annotator::AbstractString, kwargs...)
    A = load_annotator_module(annotator)
    # the annotator was defined in a newer world: run the body there
    return Base.invokelatest(ir_plot_body, A, files, svg; annotator=annotator, kwargs...)
end

function ir_plot_body(A, files::Vector{String}, svg::AbstractString; annotator::AbstractString,
                      dataset::AbstractString="", metadata::AbstractString="", refs::Vector{String}=String[],
                      min_len::Int=500, title::AbstractString="IR junctions")
    IR = A.PlastomeIR
    isempty(dataset) && (dataset = joinpath(annotator, "dataset"))
    names = Dict{String,String}()
    if !isempty(metadata)
        for m in read_tsv_rows(metadata)
            isempty(get(m, "organism", "")) || (names[m["record"]] = m["organism"])
        end
    end
    tracks = IR.IRTrack[]
    fmt(x) = x === nothing ? "." : "$(x[1])..$(x[2])"
    base = replace(svg, r"\.svg$" => "")
    check = open(base * "_ir_check.tsv", "w")
    diffs = open(base * "_ir_diffs.tsv", "w")
    println(diffs, join(("record", "kind", "copy", "pos_IRa", "pos_IRb", "length", "bases", "in_assembler_IR"), '\t'))
    println(check, join(("record", "status", "asm_IRa", "annot_IRa", "asm_IRb", "annot_IRb",
                         "asm_IR_length", "annot_IRa_length", "annot_IRb_length", "identity", "agree"), '\t'))
    for r in load_submission_records(files)
        L = length(r.seq)
        feats = IR.Feature[IR.Feature((m.pseudo ? "ψ" : "") * m.gene, m.kind,
                                      [(s.first, s.last, s.strand == '+' ? Int8(1) : Int8(-1)) for s in m.segs])
                           for m in r.models]
        st, ir, d = IR.find_ir_status(r.seq; min_len=min_len)
        kB = findfirst(==("IR1"), r.regions.names); kA = findfirst(==("IR2"), r.regions.names)
        aB = kB === nothing ? nothing : r.regions.spans[kB]
        aA = kA === nothing ? nothing : r.regions.spans[kA]
        asm_len = aB === nothing ? "." : string(aB[2] - aB[1] + 1)
        name = get(names, r.id, r.id)
        if st != :found
            @warn "annotator IR detector found no IR; row drawn as 'no IR detected'" record=r.id status=st
            push!(tracks, IR.IRTrack(name, L, nothing, feats, nothing))
            println(check, join((r.id, string(st), fmt(aB), ".", fmt(aA), ".", asm_len, ".", ".", ".", "no"), '\t'))
            continue
        end
        agree = aB !== nothing && aA !== nothing && ir.irb == aB && ir.ira == aA
        if !agree
            @warn "IR boundaries differ: assembler regions.tsv vs annotator detector" record=r.id asm_IRa=fmt(aB) annot_IRa=fmt(ir.irb) asm_IRb=fmt(aA) annot_IRb=fmt(ir.ira)
        end
        println(check, join((r.id, "found", fmt(aB), fmt(ir.irb), fmt(aA), fmt(ir.ira), asm_len,
                             ir.len_irb, ir.len_ira, @sprintf("%.6f", ir.identity), agree ? "yes" : "no"), '\t'))
        push!(tracks, IR.IRTrack(name, ir, feats, d))
        inasm(p, sp) = sp !== nothing && sp[1] <= p <= sp[2]
        for x in something(d, IR.IRDiff[])
            println(diffs, join((r.id, x.kind, copy_name(x.copy), x.pos_irb, x.pos_ira, x.len, x.bases,
                                 (inasm(x.pos_irb, aB) && inasm(x.pos_ira, aA)) ? "yes" : "no"), '\t'))
        end
    end
    close(check); close(diffs)
    if !isempty(refs)
        want = [String(first(split(a, '.'))) for a in refs]
        rfeats = IR.read_features_tsv(joinpath(dataset, "labels.tsv"); accession=want)
        org = Dict{String,String}()
        for row in read_tsv_rows(joinpath(dataset, "records.tsv"))
            org[String(first(split(row["accession"], '.')))] = get(row, "organism", row["accession"])
        end
        seqs = Dict{String,String}()
        IR.foreach_fasta(joinpath(dataset, "genomes.fasta.gz")) do id, seq
            k = String(first(split(id, '.')))
            k in want && (seqs[k] = String(seq))
        end
        for a in want
            haskey(seqs, a) || (@warn "reference not in genomes.fasta.gz; skipped" accession=a; continue)
            st, ir, d = IR.find_ir_status(seqs[a]; min_len=min_len)
            fs = get(rfeats, a, IR.Feature[])
            push!(tracks, st == :found ? IR.IRTrack(get(org, a, a), ir, fs, d) :
                                         IR.IRTrack(get(org, a, a), length(seqs[a]), nothing, fs, nothing))
        end
    end
    isempty(tracks) && error("nothing to plot")
    IR.plot_ir_junctions(tracks, svg; title=title)
    relabel_ir_svg(svg)
    println("IR junction plot: ", length(tracks), " rows written to ", svg)
    return svg
end

# ------------------------------------------------------------------ CLI

const USAGE = """
PlastomeExport.jl $VERSION_STRING
  julia PlastomeExport.jl gene-table --out DIR [--categories FILE] GENES_TSV [GENES_TSV ...]
  julia PlastomeExport.jl genbank --template METADATA_TSV GENES_TSV [GENES_TSV ...]
  julia PlastomeExport.jl genbank --metadata METADATA_TSV --out DIR [--products FILE] [--gap-min 10] GENES_TSV [...]
  julia PlastomeExport.jl map --out DIR [--metadata METADATA_TSV] [--categories FILE] [--window 500] [--step 50] GENES_TSV [...]
  julia PlastomeExport.jl ir-plot --svg FILE [--annotator DIR] [--dataset DIR] [--metadata METADATA_TSV]
                          [--ref ACC,ACC] [--min-len 500] [--title TEXT] GENES_TSV [...]
"""

function main(argv::Vector{String})
    (isempty(argv) || argv[1] in ("-h", "--help")) && (print(USAGE); return)
    cmd = argv[1]
    if cmd == "gene-table"
        out = ""; categories = joinpath(@__DIR__, "gene_categories.tsv"); files = String[]
        i = 2
        while i <= length(argv)
            a = argv[i]
            if a == "--out" && i < length(argv)
                out = argv[i+1]; i += 2
            elseif a == "--categories" && i < length(argv)
                categories = argv[i+1]; i += 2
            elseif startswith(a, "--")
                error("unknown option $a\n$USAGE")
            else
                push!(files, a); i += 1
            end
        end
        isempty(out) && error("--out is required\n$USAGE")
        isempty(files) && error("no genes tables given\n$USAGE")
        gene_table(files, out; categories=categories)
    elseif cmd == "genbank"
        opt = Dict("out" => "", "template" => "", "metadata" => "", "gap-min" => "10",
                   "products" => joinpath(@__DIR__, "gene_products.tsv"))
        files = String[]
        i = 2
        while i <= length(argv)
            a = argv[i]
            if startswith(a, "--")
                (haskey(opt, a[3:end]) && i < length(argv)) || error("unknown option $a\n$USAGE")
                opt[a[3:end]] = argv[i+1]; i += 2
            else
                push!(files, a); i += 1
            end
        end
        isempty(files) && error("no genes tables given\n$USAGE")
        if !isempty(opt["template"])
            write_template(files, opt["template"])
        else
            (isempty(opt["metadata"]) || isempty(opt["out"])) && error("--metadata and --out are required\n$USAGE")
            genbank(files, opt["out"]; metadata=opt["metadata"], products=opt["products"],
                    gap_min=parse(Int, opt["gap-min"]))
        end
    elseif cmd == "map"
        opt = Dict("out" => "", "metadata" => "", "window" => "500", "step" => "50",
                   "categories" => joinpath(@__DIR__, "gene_categories.tsv"))
        files = String[]
        i = 2
        while i <= length(argv)
            a = argv[i]
            if startswith(a, "--")
                (haskey(opt, a[3:end]) && i < length(argv)) || error("unknown option $a\n$USAGE")
                opt[a[3:end]] = argv[i+1]; i += 2
            else
                push!(files, a); i += 1
            end
        end
        isempty(files) && error("no genes tables given\n$USAGE")
        isempty(opt["out"]) && error("--out is required\n$USAGE")
        genome_maps(files, opt["out"]; categories=opt["categories"], metadata=opt["metadata"],
                    window=parse(Int, opt["window"]), step=parse(Int, opt["step"]))
    elseif cmd == "ir-plot"
        opt = Dict("svg" => "", "annotator" => joinpath(homedir(), "Assembler", "refseq_plastid"), "dataset" => "",
                   "metadata" => "", "ref" => "", "min-len" => "500", "title" => "IR junctions")
        files = String[]
        i = 2
        while i <= length(argv)
            a = argv[i]
            if startswith(a, "--")
                (haskey(opt, a[3:end]) && i < length(argv)) || error("unknown option $a\n$USAGE")
                opt[a[3:end]] = argv[i+1]; i += 2
            else
                push!(files, a); i += 1
            end
        end
        isempty(opt["svg"]) && error("--svg is required\n$USAGE")
        (isempty(files) && isempty(opt["ref"])) && error("no genes tables or references given\n$USAGE")
        ir_plot(files, opt["svg"]; annotator=opt["annotator"], dataset=opt["dataset"], metadata=opt["metadata"],
                refs=isempty(opt["ref"]) ? String[] : String.(split(opt["ref"], ',')),
                min_len=parse(Int, opt["min-len"]), title=opt["title"])
    else
        error("unknown command $cmd\n$USAGE")
    end
end

end # module PlastomeExport

if abspath(PROGRAM_FILE) == @__FILE__
    PlastomeExport.main(ARGS)
end
