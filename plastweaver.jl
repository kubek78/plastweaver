#!/usr/bin/env julia
# plastweaver.jl - paired-end reads to an assembled, ordered and annotated plastome, with gene
# tables, figures and GenBank submission files, one folder per library (Base Julia + the tools
# beside this file: SeedRef.jl, AssemblerComplete_v16.jl, PlastomeExport.jl,
# refseq_plastid/PlastomeAnnotator.jl).
#
#   julia plastweaver.jl run --input DIR [--samples a,b] [--results DIR] [--threads N]
#                            [--template FILE] [--seedref-config FILE] [--resume]
#   julia plastweaver.jl export  RESULTS/LIBRARY [RESULTS/LIBRARY ...]
#   julia plastweaver.jl genbank RESULTS/LIBRARY [RESULTS/LIBRARY ...]
#
# run      seeds and reference for every library in DIR (SeedRef, one pass, shared panel and
#          reference cache in RESULTS/_seedref), then per library the assembler (ordering and
#          annotation included) and the export step.
# export   gene tables, genome map, IR junction plot and the GenBank metadata template again
#          from an existing assembly (e.g. after a new PlastomeExport.jl).
# genbank  submission files, after the metadata table of the library has been filled in
#          (instructions in RESULTS/LIBRARY/genbank/README.txt).
#
# Layout of RESULTS/LIBRARY:
#   seed/        seeds, SeedRef reference, the assembler config used (LIBRARY_config.txt)
#   assembly/    all assembler outputs (LIBRARY_assembly.fasta, LIBRARY_run.txt, graph, ...)
#   LIBRARY_plastome.fasta / .gff3 / _regions.tsv   ordered plastome(s) and annotation
#   genes/       gene tables by function, intron table, gene summary
#   figures/     circular map(s), IR junction plot (+ IR boundary check and IR copy differences)
#   genbank/     metadata.tsv to fill in, README.txt; submission files after `genbank`
#   plastweaver.log

module PlastWeaver

using Printf
using Dates

const VERSION_STRING = "0.1.3"
const HERE = @__DIR__

include(joinpath(HERE, "PlastomeExport.jl"))

# ------------------------------------------------------------------ settings

Base.@kwdef mutable struct Settings
    results::String = abspath("results")
    threads::String = "auto"
    template::String = ""
    seedref_config::String = ""
    annotation_dir::String = joinpath(HERE, "refseq_plastid")
    dataset::String = joinpath(HERE, "refseq_plastid", "dataset")
    panel::String = joinpath(HERE, "seedref", "references", "panel")
    categories::String = joinpath(HERE, "gene_categories.tsv")
    products::String = joinpath(HERE, "gene_products.tsv")
    resume::Bool = false
end

"Assembler template: config_template.txt beside this file, else config_v16.txt."
function default_template()
    for f in ("config_template.txt", "config_v16.txt")
        p = joinpath(HERE, f)
        isfile(p) && return p
    end
    error("no assembler template (config_template.txt or config_v16.txt) beside plastweaver.jl; use --template")
end

julia(script::AbstractString, args::Vector{String}, threads::AbstractString) =
    `$(Base.julia_cmd()) --project=$HERE --threads=$threads $(joinpath(HERE, script)) $args`

function logline(lib_dir::AbstractString, msg::AbstractString)
    line = Dates.format(now(), "yyyy-mm-dd HH:MM:SS") * "  " * msg
    println(line)
    mkpath(lib_dir)
    open(io -> println(io, line), joinpath(lib_dir, "plastweaver.log"), "a")
end

"key => value of a key = value config (comments after #)."
function config_values(path::AbstractString)
    kv = Dict{String,String}()
    for raw in eachline(path)
        line = strip(first(split(raw, '#'; limit=2)))
        occursin('=', line) || continue
        k, v = strip.(split(line, '='; limit=2))
        kv[String(k)] = String(v)
    end
    return kv
end

"Copy of a config with the given keys replaced (old lines dropped, new values appended)."
function write_config(src::AbstractString, dst::AbstractString, set::Vector{Pair{String,String}})
    keys_set = Set(first.(set))
    open(dst, "w") do io
        for line in eachline(src)
            body = strip(first(split(line, '#'; limit=2)))
            key = occursin('=', body) ? strip(first(split(body, '='; limit=2))) : ""
            key in keys_set && continue
            println(io, line)
        end
        println(io, "\n# --- written by plastweaver ", VERSION_STRING, " ---")
        for (k, v) in set
            println(io, k, " = ", v)
        end
    end
    return dst
end

function read_summary(path::AbstractString)
    rows = Dict{String,Dict{String,String}}()
    isfile(path) || return rows
    lines = readlines(path)
    isempty(lines) && return rows
    hdr = split(lines[1], '\t')
    for l in lines[2:end]
        x = split(l, '\t')
        r = Dict(String(hdr[i]) => (i <= length(x) ? String(x[i]) : "") for i in eachindex(hdr))
        rows[r["library"]] = r
    end
    return rows
end

# ------------------------------------------------------------------ steps

# Files of the reference dataset needed for annotation, nearest reference and orientation
# (downloaded once by the user and unpacked into refseq_plastid/; see README.md).
const DATASET_FILES = ("genomes.fasta.gz", "records.tsv", "labels.tsv", "features.tsv", "ref_sketches.bin",
                       "refindex_all", "models_final")

function check_dataset(st::Settings)
    isfile(joinpath(st.annotation_dir, "PlastomeAnnotator.jl")) ||
        error("PlastomeAnnotator.jl not found in $(st.annotation_dir)")
    missing_files = [f for f in DATASET_FILES if !ispath(joinpath(st.dataset, f))]
    isempty(missing_files) ||
        error("reference dataset incomplete in $(st.dataset) (missing: $(join(missing_files, ", "))); " *
              "download the reference dataset (README.md, Installation) and unpack it into $(dirname(st.dataset))")
    return true
end

"Seeds, references and assembler configs for the libraries in `input` (one SeedRef run)."
function run_seedref(input::AbstractString, samples::Vector{String}, st::Settings)
    out = joinpath(st.results, "_seedref")
    args = ["--input", abspath(input), "--output", out, "--template", st.template, "--results", st.results]
    isempty(st.seedref_config) || append!(args, ["--config", st.seedref_config])
    isempty(samples) || append!(args, ["--samples", join(samples, ',')])
    println("SeedRef: ", join(args, ' '))
    run(julia("SeedRef.jl", args, st.threads))
    return read_summary(joinpath(out, "seedref_summary.tsv")), out
end

"seed/ of one library: seeds, reference, SeedRef details and the assembler config used."
function prepare_library(lib::AbstractString, row::Dict{String,String}, seedref_out::AbstractString, st::Settings)
    L = joinpath(st.results, lib)
    S = joinpath(L, "seed")
    mkpath(S)
    seed = joinpath(S, lib * "_seed.fasta")
    cp(row["seed"], seed; force=true)
    # the assembler takes LIBRARY_seeds_all.fasta beside the seed as extra rescue seeds
    isempty(row["seeds_all"]) || cp(row["seeds_all"], joinpath(S, lib * "_seeds_all.fasta"); force=true)
    for k in ("reference_genbank", "reference_fasta")
        isempty(get(row, k, "")) || cp(row[k], joinpath(S, basename(row[k])); force=true)
    end
    detail = joinpath(seedref_out, lib)
    isdir(detail) && cp(detail, joinpath(S, "details"); force=true)
    set = Pair{String,String}["Seed Input" => seed, "Output path" => joinpath(L, "assembly")]
    isdir(st.panel) && push!(set, "Plastid panel" => st.panel)        # rpoB orientation; path of this installation
    if isfile(joinpath(st.annotation_dir, "PlastomeAnnotator.jl")) && isdir(st.dataset)
        append!(set, ["Reference dataset" => st.dataset, "Annotation dir" => st.annotation_dir,
                      "Annotation dataset" => st.dataset])
    else
        logline(L, "WARNING: annotator or dataset not found ($(st.annotation_dir)); no reference orientation or annotation")
    end
    return write_config(row["assembler_config"], joinpath(S, lib * "_config.txt"), set)
end

function run_assembler(lib::AbstractString, config::AbstractString, st::Settings)
    L = joinpath(st.results, lib)
    log = joinpath(L, lib * "_assembler.log")
    logline(L, "assembler: $(basename(config)) (log $(basename(log)))")
    ok = open(log, "w") do io            # one handle: stdout and stderr interleaved in order
        success(pipeline(julia("AssemblerComplete_v16.jl", [config], st.threads); stdout=io, stderr=io))
    end
    logline(L, ok ? "assembler: done" : "assembler: FAILED (see $(basename(log)))")
    return ok
end

const GENBANK_README = """
GenBank submission for this library
===================================

1. Fill in metadata.tsv (tab-separated; a spreadsheet saved as tab-delimited text is fine).
   One row per ordered plastome of the library. Information columns (do not edit): record,
   closure, length, N_bases, N_runs, editing, ready.
     organism          scientific name; REQUIRED - a row without organism is not exported
     seqid             sequence identifier in the submission (default: the record id)
     isolate / specimen_voucher   at least one is usually expected by GenBank
     geo_loc_name      "Country: region, locality" (if the Submission Portal rejects this
                       column name, rename it to country)
     collection_date   e.g. 2025-06-14 or 14-Jun-2025
     collected_by      collector name(s)
     location          chloroplast (default) or plastid
     force             1 exports a record that is not ready (closure not conclusive, or N in
                       the sequence; such records are written without complete/circular)
   Any other column you add is written as a source modifier with that name.

2. Run:   julia plastweaver.jl genbank RESULTS/LIBRARY

3. Check report_records.tsv (exported / skipped) and report.tsv (one row per gene model;
   actions ending in _review need a decision before submission: missing product, model
   overlapping a gap, incomplete model).

4. Submit submission.fsa with submission.tbl (5-column feature table) and
   source_modifiers.tsv through the GenBank Submission Portal, or validate them first with
   NCBI table2asn (same base name, so the .tbl is found automatically).
"""

"Ordered plastome copies, gene tables, figures and the GenBank metadata template of one library."
function export_library(L::AbstractString, st::Settings)
    L = abspath(L)
    lib = basename(L)
    A = joinpath(L, "assembly")
    prefix = joinpath(A, lib)
    genes = prefix * "_plastome_genes.tsv"
    fasta = prefix * "_plastome_LSC_IR_SSC_IR.fasta"
    status = Dict{String,String}()
    if !isfile(genes) || !isfile(fasta)
        logline(L, "export: no annotated ordered plastome in assembly/ (linear result or annotation off); skipped")
        status["export"] = "no ordered plastome"
        return status
    end
    for (src, dst) in ((fasta, lib * "_plastome.fasta"), (prefix * "_plastome_LSC_IR_SSC_IR.gff3", lib * "_plastome.gff3"),
                       (fasta * ".regions.tsv", lib * "_plastome_regions.tsv"))
        isfile(src) && cp(src, joinpath(L, dst); force=true)
    end
    files = [genes]
    steps = (
        ("genes", () -> PlastomeExport.gene_table(files, joinpath(L, "genes"); categories=st.categories)),
        ("map", () -> PlastomeExport.genome_maps(files, joinpath(L, "figures"); categories=st.categories)),
        ("ir_plot", () -> PlastomeExport.ir_plot(files, joinpath(L, "figures", lib * "_ir_junctions.svg");
                                                 annotator=st.annotation_dir, dataset=st.dataset)),
        ("genbank_template", () -> begin
            G = joinpath(L, "genbank"); mkpath(G)
            meta = joinpath(G, "metadata.tsv")
            if isfile(meta)
                logline(L, "export: genbank/metadata.tsv exists; kept (delete it to write a new template)")
            else
                PlastomeExport.write_template(files, meta)
            end
            write(joinpath(G, "README.txt"), replace(GENBANK_README, "RESULTS/LIBRARY" => L))
        end),
    )
    for (name, f) in steps
        try
            f()
            status[name] = "ok"
        catch err
            msg = first(split(sprint(showerror, err), '\n'))
            status[name] = "failed: " * msg
            logline(L, "export $name FAILED: $msg")
        end
    end
    logline(L, "export: " * join(["$k=$(status[k])" for k in sort(collect(keys(status)))], ", "))
    log_readiness(L)
    return status
end

"GenBank submission files from the filled metadata table of one library."
function genbank_library(L::AbstractString, st::Settings)
    L = abspath(L)
    lib = basename(L)
    genes = joinpath(L, "assembly", lib * "_plastome_genes.tsv")
    meta = joinpath(L, "genbank", "metadata.tsv")
    isfile(genes) || error("$lib: no annotated plastome ($genes)")
    isfile(meta) || error("$lib: no genbank/metadata.tsv; run `julia plastweaver.jl export $L` first")
    n = PlastomeExport.genbank([genes], joinpath(L, "genbank"); metadata=meta, products=st.products)
    logline(L, "genbank: $n record(s) exported to genbank/ (see report_records.tsv)")
    return n
end

"""
GenBank readiness of the ordered records of a library, from genbank/metadata.tsv (written by
export): record => (closure, ready). ready = closure conclusive for the assembler (validated,
ir_collapse, read_bridge, graph_path) and no N in the sequence; other records are exported to
GenBank only with force=1 and then without [completeness=complete] / [topology=circular].
"""
function readiness(L::AbstractString)
    meta = joinpath(L, "genbank", "metadata.tsv")
    out = Dict{String,Tuple{String,Bool}}()
    isfile(meta) || return out
    lines = readlines(meta)
    isempty(lines) && return out
    hdr = split(lines[1], '\t')
    ir = findfirst(==("record"), hdr); ic = findfirst(==("closure"), hdr); iy = findfirst(==("ready"), hdr)
    (ir === nothing || ic === nothing || iy === nothing) && return out
    for l in lines[2:end]
        x = split(l, '\t')
        length(x) >= max(ir, ic, iy) || continue
        out[String(x[ir])] = (String(x[ic]), x[iy] == "1")
    end
    return out
end

"Logs every ordered record that is not ready for GenBank as complete and circular."
function log_readiness(L::AbstractString)
    for (id, (cl, ok)) in sort(collect(readiness(L)); by=first)
        ok && continue
        logline(L, "NOTE: $id (closure $cl) has an inconclusive closure or N in the sequence: " *
                   "GenBank export needs force=1 and is written without complete/circular")
    end
end

"Primary result of a library: first ordered record and its closure, from the ordering table."
function primary_result(L::AbstractString, lib::AbstractString)
    p = joinpath(L, "assembly", lib * "_plastome_ordering.tsv")
    isfile(p) || return "", "", 0
    rows = [split(l, '\t') for l in readlines(p)[2:end]]
    ordered = [r for r in rows if length(r) >= 8 && r[2] == "ordered"]
    isempty(ordered) && return "", "", 0
    return String(ordered[1][1]), String(ordered[1][8]), length(ordered)
end

function run_all(input::AbstractString, samples::Vector{String}, st::Settings)
    mkpath(st.results)
    rows, seedref_out = run_seedref(input, samples, st)
    libs = sort(collect(keys(rows)))
    isempty(samples) || (libs = [l for l in libs if l in samples])
    summary = joinpath(st.results, "plastweaver_summary.tsv")
    out = Vector{Vector{String}}()
    for (i, lib) in enumerate(libs)
        row = rows[lib]
        L = joinpath(st.results, lib)
        println("\n[", i, "/", length(libs), "] ", lib)
        res = Dict("seed" => "", "assembly" => "", "primary" => "", "closure" => "", "ready" => "", "ordered" => "", "export" => "")
        if !isempty(get(row, "error", "")) || isempty(get(row, "assembler_config", ""))
            res["seed"] = "failed: " * get(row, "error", "no assembler config")
            logline(L, "SeedRef " * res["seed"])
        else
            res["seed"] = get(row, "seed_source", "ok") * " " * get(row, "reference_accession", "")
            config = prepare_library(lib, row, seedref_out, st)
            done = isfile(joinpath(L, "assembly", lib * "_assembly.fasta"))
            if st.resume && done
                logline(L, "assembler: existing result kept (--resume)")
                res["assembly"] = "kept"
            else
                res["assembly"] = run_assembler(lib, config, st) ? "ok" : "failed"
            end
            if res["assembly"] != "failed"
                id, cl, n = primary_result(L, lib)
                res["primary"] = id; res["closure"] = cl; res["ordered"] = string(n)
                ex = export_library(L, st)
                res["export"] = all(v -> v == "ok", values(ex)) ? "ok" : join(["$k=$v" for (k, v) in ex if v != "ok"], "; ")
                rd = get(readiness(L), id, nothing)
                res["ready"] = rd === nothing ? "." : rd[2] ? "yes" : "no (force=1)"
            end
        end
        push!(out, [lib, res["seed"], res["assembly"], res["primary"], res["closure"], res["ready"], res["ordered"], res["export"]])
        open(summary, "w") do io                                  # rewritten after every library
            println(io, join(("library", "seed", "assembly", "primary_record", "closure", "genbank_ready", "ordered_records", "export"), '\t'))
            foreach(r -> println(io, join(r, '\t')), out)
        end
    end
    println("\nSummary: ", summary)
    return summary
end

# ------------------------------------------------------------------ CLI

const USAGE = """
plastweaver $VERSION_STRING
  julia plastweaver.jl run --input DIR [--samples a,b] [--results DIR] [--threads N]
                           [--template FILE] [--seedref-config FILE] [--resume]
  julia plastweaver.jl export  RESULTS/LIBRARY [...]
  julia plastweaver.jl genbank RESULTS/LIBRARY [...]
Options for all commands: --annotation-dir DIR, --dataset DIR, --panel DIR, --categories FILE, --products FILE
The reference dataset (refseq_plastid/dataset) is downloaded once; see README.md
"""

function main(argv::Vector{String})
    (isempty(argv) || argv[1] in ("-h", "--help")) && (print(USAGE); return)
    cmd = argv[1]
    st = Settings()
    input = ""; samples = String[]; dirs = String[]
    i = 2
    while i <= length(argv)
        a = argv[i]
        if a == "--resume"
            st.resume = true; i += 1; continue
        end
        if startswith(a, "--")
            i < length(argv) || error("missing value for $a\n$USAGE")
            v = argv[i+1]; i += 2
            a == "--input" ? (input = v) :
            a == "--samples" ? (samples = String.(strip.(split(v, ',')))) :
            a == "--results" ? (st.results = abspath(v)) :
            a == "--threads" ? (st.threads = v) :
            a == "--template" ? (st.template = abspath(v)) :
            a == "--seedref-config" ? (st.seedref_config = abspath(v)) :
            a == "--annotation-dir" ? (st.annotation_dir = abspath(v)) :
            a == "--dataset" ? (st.dataset = abspath(v)) :
            a == "--panel" ? (st.panel = abspath(v)) :
            a == "--categories" ? (st.categories = abspath(v)) :
            a == "--products" ? (st.products = abspath(v)) :
            error("unknown option $a\n$USAGE")
        else
            push!(dirs, a); i += 1
        end
    end
    if cmd == "run"
        isempty(input) && error("--input is required\n$USAGE")
        isdir(input) || error("input folder not found: $input")
        isempty(st.template) && (st.template = default_template())
        isfile(st.template) || error("template not found: $(st.template)")
        check_dataset(st)
        run_all(input, samples, st)
    elseif cmd == "export"
        isempty(dirs) && error("no library folder given\n$USAGE")
        foreach(L -> export_library(L, st), dirs)
    elseif cmd == "genbank"
        isempty(dirs) && error("no library folder given\n$USAGE")
        foreach(L -> genbank_library(L, st), dirs)
    else
        error("unknown command $cmd\n$USAGE")
    end
end

end # module PlastWeaver

if abspath(PROGRAM_FILE) == @__FILE__
    PlastWeaver.main(ARGS)
end
