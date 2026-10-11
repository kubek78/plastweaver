#!/usr/bin/env julia
#=
PlastomeAnnotator.jl - Chloë-style plastome annotation in one file (Base Julia only).

  annotate FASTA files (production index and models by default):
    julia --threads auto PlastomeAnnotator.jl annotate --fasta genome.fa [--fasta ...] [--out DIR]
                                     [--index DIR] [--models DIR] [--exclude-genus G] [--near-ref 0.01]
  evaluate on a held-out split:
    julia --threads auto PlastomeAnnotator.jl annotate --eval validation [--exclude genus|accession]
  build the reference index, training candidates, existence models:
    julia --threads auto PlastomeAnnotator.jl build-index    [--splits all] ...
    julia --threads auto PlastomeAnnotator.jl build-training [--index DIR] [--out training_all] ...
    julia --threads auto PlastomeAnnotator.jl train          [--training training_all] [--out models_final] [--final 1]
  IR scan of the dataset and IR diagnostics:
    julia --threads auto PlastomeAnnotator.jl ir-scan [--dataset DIR] | ir-scan --plot ACC,... --svg F
    julia PlastomeAnnotator.jl ir-check ACC,ACC [genomes.fasta.gz]
    julia PlastomeAnnotator.jl ir-diag  ACC,ACC [--dataset DIR]

Library use (AssemblerComplete: Annotation dir = the directory of this file):
    include("PlastomeAnnotator.jl")
    a = PlastomeAnnotator.load_annotator(dataset_dir)
    models, ir, editing, nearest = PlastomeAnnotator.annotate_record(a, id, sequence)
    PlastomeAnnotator.write_gff(io, id, length(sequence), models; header=true)

Contents: modules PlastomeIR (inverted repeats), PlastomeRef (reference index, MinHash),
PlastomeAnnot (alignment and transfer), PlastomeFeatures (candidate features),
PlastomeTrees (gradient-boosted trees), PlastomeProtein (protein search), the annotator
(gene models, ORF checks, evaluation), and the subcommands BuildIndex, BuildTraining,
TrainExistence, IRScan, IRCheck, IRDiag.
=#
module PlastomeAnnotator

# .gz input: CodecZlib (pure Julia, any operating system) when the active environment provides
# it - plastweaver, the assembler and `julia --project=...` runs; otherwise the gzip command.
const CODECZLIB = try
    Base.require(Base.PkgId(Base.UUID("944b1d66-785c-5afd-91f1-9de20f533193"), "CodecZlib"))
catch
    nothing
end

"Opens a plain or gzip-compressed file for reading (decompressed stream for .gz)."
function open_maybe_gz(path::AbstractString)
    endswith(path, ".gz") || return open(path)
    CODECZLIB === nothing && return open(`gzip -dc $path`)
    return Base.invokelatest(CODECZLIB.GzipDecompressorStream, open(path))
end

# ======================================================================
# PlastomeIR
# ======================================================================
"""
    PlastomeIR

Inverted repeat (IR) detection, IR copy comparison, region assignment, canonical
rotation and IRscope-style junction plots (SVG) for plastid genomes. Base Julia only.

Coordinates are 1-based and circular; a region `(first, last)` with `last < first`
wraps across the origin. Genome layout follows the usual convention
LSC - IRb - SSC - IRa, with junctions JLB (LSC/IRb), JSB (IRb/SSC), JSA (SSC/IRa)
and JLA (IRa/LSC); a junction is stored as the first base of the region to its right.
"""
module PlastomeIR

export IRResult, IRDiff, Feature, IRTrack,
       find_ir, find_ir_status, find_direct_repeat, ir_alignment, ir_differences, junctions,
       region_of, region_vector,
       rotate_canonical, map_position, needs_ssc_flip,
       read_fasta, foreach_fasta, read_features_tsv, parse_segments,
       plot_ir_junctions

# ------------------------------------------------------------------ sequence tables

const INVALID = typemax(UInt64)

const ENC = let t = fill(0xff, 256)
    for (c, v) in zip("ACGT", 0x00:0x03)
        t[Int(c) + 1] = v
    end
    t
end

const COMPL = let t = fill(UInt8('N'), 256)
    for (a, b) in zip("ACGTN", "TGCAN")
        t[Int(a) + 1] = UInt8(b)
    end
    t
end

isacgt(b::UInt8) = ENC[Int(b) + 1] != 0xff
compl(b::UInt8) = COMPL[Int(b) + 1]
seqbytes(seq::AbstractString) = Vector{UInt8}(uppercase(String(seq)))

"Circular inclusion of position p in region (a, b) on a genome of length L."
incirc(p::Integer, r::Tuple{Int,Int}) = r[1] <= r[2] ? (r[1] <= p <= r[2]) : (p >= r[1] || p <= r[2])

# ------------------------------------------------------------------ types

"""
    IRResult

`irb`, `ira`: (first, last) base of each IR copy; `len_*`: lengths of the IR copies,
`lsc`, `ssc`: lengths of the single-copy regions; `gc`: GC fraction per region;
`drift`: (min, max) diagonal offset j - i of the k-mer anchors in the alignment of IRb
with the reverse complement of IRa (sets the alignment band); `identity`: matching /
(matching + mismatched) columns among aligned ACGT positions, so a single long indel
does not lower it; `aligned_frac`: aligned (non-gap) columns / length of the shorter
copy; both NaN when not computed; `n_unknown`: aligned columns involving a non-ACGT
base, -1 when not computed (no differences with n_unknown > 0 does not mean identical
copies).
"""
struct IRResult
    L::Int
    irb::Tuple{Int,Int}
    ira::Tuple{Int,Int}
    len_irb::Int
    len_ira::Int
    lsc::Int
    ssc::Int
    gc::NamedTuple{(:LSC, :IRb, :SSC, :IRa),NTuple{4,Float64}}
    drift::Tuple{Int,Int}
    identity::Float64
    aligned_frac::Float64
    n_unknown::Int
end

"""
    IRDiff

One difference between the IR copies, from the alignment of IRb with the reverse
complement of IRa. `kind` is `:snp` or `:indel`; `copy` is `:both` for SNPs and the
copy carrying the extra bases (`:IRb` or `:IRa`) for indels; `pos_irb` / `pos_ira` are
genome positions of the difference (or its flank) in each copy; `bases` is given in
IRb orientation (`"A>G"` for SNPs: IRb base, then the IRa base on the IRb strand).
"""
struct IRDiff
    kind::Symbol
    copy::Symbol
    pos_irb::Int
    pos_ira::Int
    len::Int
    bases::String
end

"Annotated feature; `segs` are (first, last, strand) in genome coordinates."
struct Feature
    name::String
    kind::String
    segs::Vector{Tuple{Int,Int,Int8}}
end

"""
One genome row of a junction plot; `ir === nothing` draws a 'no IR detected' row,
`diffs === nothing` leaves out the IR copy comparison.
"""
struct IRTrack
    name::String
    L::Int
    ir::Union{IRResult,Nothing}
    features::Vector{Feature}
    diffs::Union{Vector{IRDiff},Nothing}
end

IRTrack(name::AbstractString, ir::IRResult, features::Vector{Feature}=Feature[],
        diffs::Union{Vector{IRDiff},Nothing}=nothing) =
    IRTrack(String(name), ir.L, ir, features, diffs)

# ------------------------------------------------------------------ IR detection

function kmer_codes(s::Vector{UInt8}, k::Int)
    L = length(s)
    codes = fill(INVALID, L)
    mask = k == 32 ? typemax(UInt64) : (UInt64(1) << (2k)) - 1
    v = UInt64(0)
    valid = 0
    for p in 1:(L + k - 1)
        b = s[mod1(p, L)]
        c = ENC[Int(b) + 1]
        if c == 0xff
            valid = 0
            v = UInt64(0)
        else
            v = ((v << 2) | UInt64(c)) & mask
            valid += 1
        end
        st = p - k + 1
        if st >= 1 && st <= L && valid >= k
            codes[st] = v
        end
    end
    return codes
end

function rc_code(v::UInt64, k::Int)
    r = UInt64(0)
    for _ in 1:k
        r = (r << 2) | (UInt64(3) - (v & UInt64(3)))
        v >>= 2
    end
    return r
end

"Ungapped X-drop extension: x moves by dx, y by dy, base x compared with complement(y)."
function ungapped(s::Vector{UInt8}, x::Int, y::Int, dx::Int, dy::Int, xdrop::Int, limit::Int)
    L = length(s)
    score = 0
    bestsc = 0
    bestt = 0
    t = 0
    while t < limit
        t += 1
        a = s[mod1(x + dx * t, L)]
        b = compl(s[mod1(y + dy * t, L)])
        score += (a == b && isacgt(a)) ? 1 : -4
        if score > bestsc
            bestsc = score
            bestt = t
        end
        score < bestsc - xdrop && break
    end
    return bestt
end

"""
    gapped_extend(s, x, y, dx, dy, xdrop, limit) -> (tx, ty)

Extends a pairing of base x with the complement of base y outwards (x by dx, y by dy);
returns the number of bases added on each side (tx + ty <= limit). After each ungapped
stretch it tries to resume past an indel (gx or gy bases skipped, up to max_gap) or a
short substitution block (gx = gy <= max_sub_block), starting up to max_back bases
before the stop, because chance matches can carry the stretch past the indel; the
resumed stretch must begin with min_after exact matches.
"""
function gapped_extend(s::Vector{UInt8}, x::Int, y::Int, dx::Int, dy::Int, xdrop::Int, limit::Int;
                       max_gap::Int=30, min_after::Int=10, max_sub_block::Int=3, max_back::Int=3)
    L = length(s)
    tx = 0
    ty = 0
    while true
        t = ungapped(s, x + dx * tx, y + dy * ty, dx, dy, xdrop, (limit - tx - ty) ÷ 2)
        tx += t
        ty += t
        found = false
        nx = ny = 0
        for g in 1:max_gap
            for back in 0:max_back
                for (gx, gy) in ((g, 0), (0, g), (g, g))
                    (gx == gy && g > max_sub_block) && continue
                    bx, by = tx - back, ty - back
                    bx + gx + by + gy + 2min_after > limit && continue
                    ok = true
                    for r in 1:min_after
                        a = s[mod1(x + dx * (bx + gx + r), L)]
                        b = compl(s[mod1(y + dy * (by + gy + r), L)])
                        if a != b || !isacgt(a)
                            ok = false
                            break
                        end
                    end
                    if ok
                        found = true
                        nx, ny = bx + gx, by + gy
                        break
                    end
                end
                found && break
            end
            found && break
        end
        found || return (tx, ty)
        tx, ty = nx, ny
    end
end

function gc_fraction(s::Vector{UInt8}, first::Int, len::Int)
    len <= 0 && return NaN
    L = length(s)
    gc = 0
    acgt = 0
    for t in 0:(len - 1)
        b = s[mod1(first + t, L)]
        if isacgt(b)
            acgt += 1
            gc += (b == UInt8('G') || b == UInt8('C'))
        end
    end
    return acgt == 0 ? NaN : gc / acgt
end

"Positions of each forward k-mer."
function kmer_occurrences(codes::Vector{UInt64})
    occ = Dict{UInt64,Vector{Int}}()
    sizehint!(occ, length(codes))
    for i in eachindex(codes)
        v = codes[i]
        v == INVALID && continue
        push!(get!(() -> Int[], occ, v), i)
    end
    return occ
end

"At least 8 distinct 3-mers in the k-mer (skips low-complexity k-mers such as (AT)n)."
function complex_enough(v::UInt64, k::Int)
    seen = falses(64)
    n = 0
    for t in 0:(k - 3)
        c = Int((v >> (2t)) & 0x3f) + 1
        if !seen[c]
            seen[c] = true
            n += 1
        end
    end
    return n >= 8
end

"Signed circular offset in (-L/2, L/2]."
circoff(d::Int, L::Int) = mod(d + L ÷ 2, L) - L ÷ 2

"""
Anchor pairs (x, y) (k-mer starts in copy A and copy B) of the outermost
constant-diagonal runs of >= min_seg k-mer matches within `G` bases of the end
(a0, b0) that reach beyond it in either copy. K-mers of any copy number are used
(low-complexity ones are skipped), so a sequence duplicated next to one IR copy (a
tandem copy of the IR end) is found. A run is kept when it has at least one anchor per
10 nt of the indel it implies. `outer`: the end is the A start / B end (A extends
downwards, B upwards); otherwise the A end / B start.
"""
function boundary_blocks(codes::Vector{UInt64}, L::Int, k::Int, a0::Int, b0::Int, outer::Bool,
                         G::Int, gap_max::Int, min_seg::Int)
    idx = Dict{UInt64,Vector{Int}}()
    for t in -G:G
        y = mod1(b0 + t, L)
        v = codes[y]
        (v == INVALID || !complex_enough(v, k)) && continue
        push!(get!(() -> Int[], idx, rc_code(v, k)), y)
    end
    diag = Dict{Int,Vector{NTuple{4,Int}}}()
    for t in -G:G
        x = mod1(a0 + t, L)
        v = codes[x]
        v == INVALID && continue
        ys = get(idx, v, nothing)
        ys === nothing && continue
        for y in ys
            eA, eB = outer ? (-circoff(x - a0, L), circoff(y + k - 1 - b0, L)) :
                             (circoff(x + k - 1 - a0, L), -circoff(y - b0, L))
            push!(get!(() -> NTuple{4,Int}[], diag, mod(x + y, L)), (eA, eB, x, y))
        end
    end
    res = Tuple{Int,Int}[]
    for pts in values(diag)
        sort!(pts)
        n = 1
        for q in (length(pts) - 1):-1:1        # outermost run along the diagonal
            pts[q + 1][1] - pts[q][1] > gap_max && break
            n += 1
        end
        eA, eB, x, y = pts[end]
        if n >= max(min_seg, abs(eA - eB) / 10) && max(eA, eB) > 0
            push!(res, (x, y))
        end
    end
    return sort!(res)
end

"""
Leading gap columns when aligning W bases from x0 (step dx) with the complement of W
bases from y0 (step dy); 0 means the two copies' terminal bases pair up.
"""
function terminal_gap(s::Vector{UInt8}, x0::Int, y0::Int, dx::Int, dy::Int; W::Int=400)
    L = length(s)
    sa = [s[mod1(x0 + dx * t, L)] for t in 0:(W - 1)]
    sb = [compl(s[mod1(y0 + dy * t, L)]) for t in 0:(W - 1)]
    ops, _ = align_affine(sa, sb, 60, 6, 1, 1)
    g = 0
    for op in ops
        op == 0x01 && break
        g += 1
    end
    return g
end

"""
    find_ir_status(seq; k=31, min_len=500, max_shift=2000, gap_max=300, chain_gap=3000,
                   ext_shift=10000, ext_gap=10000, min_dup_dist=1000, seg_jump=30,
                   min_seg=20, xdrop=15, min_identity=0.8, min_aligned=0.5,
                   validate=true) -> (status, ir, diffs)

`status` is :found, :not_detected, :low_identity or :inconsistent; `ir` is an
IRResult (also returned for :low_identity), `diffs` the IR copy differences
(`nothing` when not computed).

K-mers whose reverse complement occurs elsewhere give base pairings p + q = S (mod L);
k-mers occurring once, or twice with copies >= min_dup_dist apart (an IR overlapping a
direct duplication, e.g. three rRNA operons), are used; close duplicates (tandem
repeats near IR ends) are not, as they create spurious pairings. The densest window of
S (width 2*max_shift; indels accumulated along an IR shift S) selects the anchors. The
map p -> S - p has two fixed points, S/2 and S/2 + L/2, and each IR copy lies in one
of the two halves they define, so copies separated by a very short SSC are not merged.
The anchors of one half are split into segments of constant pairing sum (jumps <=
seg_jump, gaps <= gap_max; an indel changes the sum, a stray inverted match elsewhere
in the genome has another sum); segments with fewer than min_seg anchors are dropped.
Segments are chained when collinear (ascending in copy A, partners descending in copy
B, gaps <= chain_gap in both copies, and >= 1 anchor per 10 nt of the implied indel in
the smaller segment), so an IR interrupted by indels or divergent stretches is taken
whole, and the chain with most anchors gives the two copies (copy B = partners of the
chained anchors). Each end of the pair is extended with `gapped_extend` from the
extreme anchor of either copy and, followed up to three times, from duplicated blocks
found beyond the current end (`boundary_blocks`, e.g. a tandem copy of the IR end next
to one IR copy); the candidate whose terminal bases pair in an alignment of the last
400 bases of both copies (no terminal gap) is kept, then the longer one. The copy
preceded by the longer single-copy region is IRb. With `validate`, the copies are
aligned (`ir_alignment`) and a pair with substitution identity below `min_identity`
or with less than `min_aligned` of the shorter copy aligned is reported as
:low_identity.

Two passes: the first with max_shift / chain_gap; when it finds an IR, a second pass
with ext_shift / ext_gap (indels of several kb between IR blocks, e.g. a 7.6-kb
insertion between rrn16 and rrn23 in one copy) replaces it only if it is :found, both
of its copies contain the first-pass copies, and it is longer. The second pass thus
only extends a found IR. ext_shift <= max_shift and ext_gap <= chain_gap turn it off.
"""
function find_ir_status(seq::AbstractString; max_shift::Int=2000, chain_gap::Int=3000,
                        ext_shift::Int=10000, ext_gap::Int=10000, kwargs...)
    r1 = ir_pass(seq; max_shift = max_shift, chain_gap = chain_gap, kwargs...)
    (r1[1] != :found || (ext_shift <= max_shift && ext_gap <= chain_gap)) && return r1
    r2 = ir_pass(seq; max_shift = ext_shift, chain_gap = ext_gap, kwargs...)
    r2[1] == :found || return r1
    a, b = r2[2], r1[2]
    cov(x, y) = incirc(y[1], x) && incirc(y[2], x)
    contains = (cov(a.irb, b.irb) && cov(a.ira, b.ira)) || (cov(a.irb, b.ira) && cov(a.ira, b.irb))
    return contains && a.len_irb + a.len_ira > b.len_irb + b.len_ira ? r2 : r1
end

"Single IR detection pass of `find_ir_status` (see there)."
function ir_pass(seq::AbstractString; k::Int=31, min_len::Int=500, max_shift::Int=2000,
                        gap_max::Int=300, chain_gap::Int=3000, min_dup_dist::Int=1000,
                        seg_jump::Int=30, min_seg::Int=20, xdrop::Int=15, min_identity::Float64=0.8, min_aligned::Float64=0.5,
                        validate::Bool=true)
    1 <= k <= 32 || throw(ArgumentError("k must be in 1:32"))
    s = seqbytes(seq)
    L = length(s)
    none = (:not_detected, nothing, nothing)
    L > 4k || return none
    codes = kmer_codes(s, k)

    occ = kmer_occurrences(codes)
    usable(ps) = length(ps) == 1 ||
                 (length(ps) == 2 && min(abs(ps[2] - ps[1]), L - abs(ps[2] - ps[1])) >= min_dup_dist)
    pa = Int[]
    pj = Int[]
    psum = Int[]
    for (v, ps) in occ
        usable(ps) || continue
        r = rc_code(v, k)
        r == v && continue
        qs = get(occ, r, nothing)
        (qs === nothing || !usable(qs)) && continue
        for i in ps, j in qs
            push!(pa, i)
            push!(pj, j)
            push!(psum, mod(i + j + k - 1, L))  # base i + t pairs with base j + k - 1 - t
        end
    end
    n = length(psum)
    n == 0 && return none

    # densest circular window of the pairing constant
    ord = sortperm(collect(zip(psum, pa, pj)))   # deterministic order of ties
    sv = psum[ord]
    extval(t) = t <= n ? sv[t] : sv[t - n] + L
    W = 2max_shift
    best = 0
    bestlo = 1
    hi = 1
    for lo in 1:n
        hi = max(hi, lo)
        while hi + 1 <= lo + n - 1 && extval(hi + 1) - sv[lo] <= W
            hi += 1
        end
        if hi - lo + 1 > best
            best = hi - lo + 1
            bestlo = lo
        end
    end
    partner = Dict{Int,Int}()
    for t in 0:(best - 1)
        q = ord[mod1(bestlo + t, n)]
        partner[pa[q]] = pj[q]
    end
    selvals = sort!([extval(bestlo + t) for t in 0:(best - 1)])
    S0 = mod(selvals[(length(selvals) + 1) ÷ 2], L)

    # split anchors by the mirror centres: copy A lies in the half [c, c + L/2)
    c = S0 / 2
    uA(p) = mod(p - c, L)
    uB(p) = mod(p - c - L / 2, L)
    pairsum(x) = mod(x + partner[x] + k - 1, L)
    side1 = Tuple{Float64,Int}[]
    for x in keys(partner)
        u = mod(x - c, L)
        u < L / 2 && push!(side1, (u, x))
    end
    (isempty(side1) || length(side1) == length(partner)) && return none
    sort!(side1)
    # segments: anchors with a constant pairing sum (jumps <= seg_jump) and gaps <= gap_max;
    # segments with fewer than min_seg anchors (chance matches) are dropped
    segs = Vector{Vector{Tuple{Float64,Int}}}()
    cur = [side1[1]]
    for q in 2:length(side1)
        p = side1[q]
        ds = abs(pairsum(p[2]) - pairsum(cur[end][2]))
        ds = min(ds, L - ds)
        if p[1] - cur[end][1] <= gap_max && ds <= seg_jump
            push!(cur, p)
        else
            push!(segs, cur)
            cur = [p]
        end
    end
    push!(segs, cur)
    filter!(g -> length(g) >= min_seg, segs)
    isempty(segs) && return none
    # chain collinear segments (A ascending, B descending, gaps <= chain_gap in both
    # copies); the smaller of two chained segments needs >= 1 anchor per 10 nt of the
    # indel implied by the two gaps
    info = [(g[1][1], g[end][1], minimum(uB(partner[x]) for (_, x) in g),
             maximum(uB(partner[x]) for (_, x) in g), length(g)) for g in segs]
    nseg = length(segs)                      # segs are ordered by A start
    chain_sc = zeros(Int, nseg)
    chain_prev = zeros(Int, nseg)
    for q in 1:nseg
        a0, _, _, b1, nq = info[q]
        chain_sc[q] = nq
        for q0 in 1:(q - 1)
            _, c1, d0, _, n0 = info[q0]
            gA, gB = a0 - c1, d0 - b1
            if gA > -k && gB > -k && gA <= chain_gap && gB <= chain_gap &&
               min(nq, n0) >= max(min_seg, abs(gA - gB) / 10) && chain_sc[q0] + nq > chain_sc[q]
                chain_sc[q], chain_prev[q] = chain_sc[q0] + nq, q0
            end
        end
    end
    q = argmax(chain_sc)
    A = Tuple{Float64,Int}[]
    while q > 0
        append!(A, segs[q])
        q = chain_prev[q]
    end
    sort!(A)
    Apts = [x for (_, x) in A]
    ub = [uB(partner[x]) for x in Apts]
    (A[end][1] - A[1][1] + k < min_len || maximum(ub) - minimum(ub) + k < min_len) && return none

    # ends: candidate anchor pairs are the extreme anchor of either copy and, followed up
    # to 3 times, duplicated blocks beyond the current end (boundary_blocks); the
    # candidate whose terminal bases pair (no terminal gap in an alignment of the last
    # 400 bases of both copies) is kept, then the longer one
    function ext_out(i, j)                   # A k-mer at i pairs with B k-mer at j
        lim = mod(i - (j + k - 1) - 1, L)    # bases between B's end and A's start
        tx, ty = gapped_extend(s, i, j + k - 1, -1, +1, xdrop, lim)
        as_, be_ = mod1(i - tx, L), mod1(j + k - 1 + ty, L)
        return (terminal_gap(s, as_, be_, +1, -1), -(uB(be_) - uA(as_)), as_, be_)
    end
    function ext_in(i, j)
        lim = mod(j - (i + k - 1) - 1, L)    # bases between A's end and B's start
        tx, ty = gapped_extend(s, i + k - 1, j, +1, -1, xdrop, lim)
        ae_, bs_ = mod1(i + k - 1 + tx, L), mod1(j - ty, L)
        return (terminal_gap(s, ae_, bs_, -1, +1), -(uA(ae_) - uB(bs_)), ae_, bs_)
    end
    outs = [ext_out(i, partner[i]) for i in unique([A[1][2], argmax(x -> uB(partner[x] + k - 1), Apts)])]
    ins = [ext_in(i, partner[i]) for i in unique([A[end][2], argmin(x -> uB(partner[x]), Apts)])]
    for _ in 1:3
        cur_o = minimum(outs)
        for (x, y) in boundary_blocks(codes, L, k, cur_o[3], cur_o[4], true, chain_gap, gap_max, min_seg)
            push!(outs, ext_out(x, y))
        end
        minimum(outs) == cur_o && break
    end
    for _ in 1:3
        cur_i = minimum(ins)
        for (x, y) in boundary_blocks(codes, L, k, cur_i[3], cur_i[4], false, chain_gap, gap_max, min_seg)
            push!(ins, ext_in(x, y))
        end
        minimum(ins) == cur_i && break
    end
    _, _, a_s, b_e = minimum(outs)
    _, _, a_e, b_s = minimum(ins)

    lenA = mod(a_e - a_s, L) + 1
    lenB = mod(b_e - b_s, L) + 1
    gap_ab = mod(b_s - a_e - 1, L)          # after copy A, before copy B
    gap_ba = mod(a_s - b_e - 1, L)          # after copy B, before copy A
    if gap_ba >= gap_ab                      # LSC precedes A: A is IRb
        irb, ira, lirb, lira, lsc, ssc = (a_s, a_e), (b_s, b_e), lenA, lenB, gap_ba, gap_ab
    else
        irb, ira, lirb, lira, lsc, ssc = (b_s, b_e), (a_s, a_e), lenB, lenA, gap_ab, gap_ba
    end
    lsc + ssc + lirb + lira == L || return (:inconsistent, nothing, nothing)   # overlapping copies

    # diagonal offsets of the anchors in the IRb x rc(IRa) alignment
    base = irb[1] + ira[2]
    dlo, dhi = min(0, lira - lirb), max(0, lira - lirb)
    for x in Apts
        off = mod(base - (x + partner[x] + k - 1) + L ÷ 2, L) - L ÷ 2
        dlo, dhi = min(dlo, off), max(dhi, off)
    end
    gc = (LSC = gc_fraction(s, mod1(ira[2] + 1, L), lsc),
          IRb = gc_fraction(s, irb[1], lirb),
          SSC = gc_fraction(s, mod1(irb[2] + 1, L), ssc),
          IRa = gc_fraction(s, ira[1], lira))
    ir = IRResult(L, irb, ira, lirb, lira, lsc, ssc, gc, (dlo, dhi), NaN, NaN, -1)
    validate || return (:found, ir, nothing)
    al = ir_alignment(s, ir)
    ir = IRResult(L, irb, ira, lirb, lira, lsc, ssc, gc, (dlo, dhi), al.identity, al.aligned_frac, al.n_unknown)
    ok = al.identity >= min_identity && al.aligned_frac >= min_aligned
    return (ok ? :found : :low_identity, ir, al.diffs)
end

"""
    find_ir(seq; kwargs...) -> IRResult or nothing

`find_ir_status` reduced to the IR of a :found pair (see there for the method).
"""
function find_ir(seq::AbstractString; kwargs...)
    st, ir, _ = find_ir_status(seq; kwargs...)
    return st == :found ? ir : nothing
end

"""
    find_direct_repeat(seq; k=31, min_len=1000, max_shift=2000, gap_max=300,
                       min_dup_dist=1000) -> NamedTuple or nothing

Largest repeat with both copies on the same strand: k-mers occurring exactly twice
(copies >= min_dup_dist apart) give shifts d = j - i; the densest shift window (width
2*max_shift) is taken and its longest run of first-copy positions (gaps <= gap_max)
is the repeat. Returns (len, start1, end1, start2, end2), lengths from anchors only
(no boundary extension), or `nothing` below min_len.
"""
function find_direct_repeat(seq::AbstractString; k::Int=31, min_len::Int=1000, max_shift::Int=2000,
                            gap_max::Int=300, min_dup_dist::Int=1000)
    s = seqbytes(seq)
    L = length(s)
    L > 4k || return nothing
    occ = kmer_occurrences(kmer_codes(s, k))
    p1 = Int[]
    dv = Int[]
    for ps in values(occ)
        length(ps) == 2 || continue
        dd = ps[2] - ps[1]
        min(dd, L - dd) >= min_dup_dist || continue
        push!(p1, ps[1])
        push!(dv, dd)
    end
    isempty(dv) && return nothing
    ord = sortperm(collect(zip(dv, p1)))
    sd = dv[ord]
    n = length(sd)
    best = 0
    bestlo = 1
    hi = 1
    for lo in 1:n
        hi = max(hi, lo)
        while hi + 1 <= n && sd[hi + 1] - sd[lo] <= 2max_shift
            hi += 1
        end
        if hi - lo + 1 > best
            best = hi - lo + 1
            bestlo = lo
        end
    end
    pts = sort!([(p1[ord[t]], dv[ord[t]]) for t in bestlo:(bestlo + best - 1)])
    bspan = -1
    bs = be = st = 1
    for q in 2:(length(pts) + 1)
        if q > length(pts) || pts[q][1] - pts[q - 1][1] > gap_max
            sp = pts[q - 1][1] - pts[st][1]
            if sp > bspan
                bspan = sp
                bs, be = st, q - 1
            end
            st = q
        end
    end
    run = pts[bs:be]
    len1 = run[end][1] - run[1][1] + k
    len1 < min_len && return nothing
    return (len = len1, start1 = run[1][1], end1 = mod1(run[end][1] + k - 1, L),
            start2 = mod1(run[1][1] + run[1][2], L), end2 = mod1(run[end][1] + run[end][2] + k - 1, L))
end

# ------------------------------------------------------------------ IR copy differences

"Banded Gotoh alignment of a with b; returns (ops, touched_band_edge). ops: 1 diag, 2 a-only, 3 b-only."
function align_affine(a::Vector{UInt8}, b::Vector{UInt8}, w::Int, go::Int, ge::Int, mismatch::Int)
    n, m = length(a), length(b)
    Wd = 2w + 1
    INF = typemax(Int) ÷ 4
    # tb[d, i + 1], d = j - i + w + 1: bits 0-1 predecessor state of M, 2-3 of X, 4-5 of Y
    # (0 = M, 1 = X: a base against a gap, 2 = Y: b base against a gap)
    tb = zeros(UInt8, Wd, n + 1)
    pM, pX, pY = fill(INF, Wd), fill(INF, Wd), fill(INF, Wd)
    cM, cX, cY = fill(INF, Wd), fill(INF, Wd), fill(INF, Wd)
    pM[w + 1] = 0
    for j in 1:min(m, w)
        pY[j + w + 1] = go + ge * j
        tb[j + w + 1, 1] = (j == 1 ? 0x00 : 0x02) << 4
    end
    argmin3(x, y, z) = x <= y ? (x <= z ? (x, 0x00) : (z, 0x02)) : (y <= z ? (y, 0x01) : (z, 0x02))
    for i in 1:n
        fill!(cM, INF); fill!(cX, INF); fill!(cY, INF)
        ai = a[i]
        for j in max(0, i - w):min(m, i + w)
            d = j - i + w + 1
            code = 0x00
            if j >= 1
                y = b[j]
                sub = (ai == y || !isacgt(ai) || !isacgt(y)) ? 0 : mismatch
                v, k = argmin3(pM[d], pX[d], pY[d])
                if v < INF
                    cM[d] = v + sub
                    code |= k
                end
            end
            if d + 1 <= Wd
                v, k = argmin3(pM[d + 1] + go + ge, pX[d + 1] + ge, pY[d + 1] + go + ge)
                if v < INF
                    cX[d] = v
                    code |= k << 2
                end
            end
            if j >= 1 && d - 1 >= 1
                v, k = argmin3(cM[d - 1] + go + ge, cX[d - 1] + go + ge, cY[d - 1] + ge)
                if v < INF
                    cY[d] = v
                    code |= k << 4
                end
            end
            tb[d, i + 1] = code
        end
        pM, cM = cM, pM
        pX, cX = cX, pX
        pY, cY = cY, pY
    end
    dfin = m - n + w + 1
    _, state = argmin3(pM[dfin], pX[dfin], pY[dfin])
    ops = UInt8[]
    edge = false
    i, j = n, m
    while i > 0 || j > 0
        dd = j - i + w + 1
        (dd == 1 || dd == Wd) && (edge = true)
        code = tb[dd, i + 1]
        if state == 0x00
            push!(ops, 0x01)
            state = code & 0x03
            i -= 1
            j -= 1
        elseif state == 0x01
            push!(ops, 0x02)
            state = (code >> 2) & 0x03
            i -= 1
        else
            push!(ops, 0x03)
            state = (code >> 4) & 0x03
            j -= 1
        end
        (i < 0 || j < 0) && error("align_affine: traceback left the matrix")
    end
    reverse!(ops)
    return ops, edge
end

"""
    ir_alignment(seq, ir; band_extra=100, max_band_extra=1000, gap_open=6, gap_extend=1,
                 mismatch=1) -> (diffs, identity, aligned_frac, n_gap, n_unknown, band, band_limited)

Banded global alignment with affine gap costs (Gotoh; a gap of length l costs
gap_open + l*gap_extend) of IRb with the reverse complement of IRa. The band covers
the anchor drift `ir.drift` (compensating indels shift the alignment diagonal even when
the copies have equal length) plus `band_extra`; when the path touches the band edge,
`band_extra` is doubled up to `max_band_extra` (`band_limited` reports a path still at
the edge). To bound memory the half-width is capped at `max_band`, but never below the
anchor drift / length difference + `band_extra`, so a multi-kb indel stays inside the
band. `identity` is
matches / (matches + mismatches) over aligned ACGT columns,
`aligned_frac` the aligned columns over the shorter copy, `n_gap` the gap columns;
columns with a non-ACGT base are counted in `n_unknown`, not as matches.
"""
function ir_alignment(s::Vector{UInt8}, ir::IRResult; band_extra::Int=100, max_band_extra::Int=1000,
                      max_band::Int=4000, gap_open::Int=6, gap_extend::Int=1, mismatch::Int=1)
    L = ir.L
    length(s) == L || throw(ArgumentError("sequence length differs from ir.L"))
    a = [s[mod1(ir.irb[1] + t, L)] for t in 0:(ir.len_irb - 1)]
    b = [compl(s[mod1(ir.ira[2] - t, L)]) for t in 0:(ir.len_ira - 1)]
    n, m = length(a), length(b)
    dlo, dhi = ir.drift
    extra = band_extra
    local ops, edge, w
    while true
        base = max(abs(dlo), abs(dhi), abs(n - m))
        w = min(base + extra, max(max_band, base + band_extra))
        ops, edge = align_affine(a, b, w, gap_open, gap_extend, mismatch)
        (!edge || extra >= max_band_extra) && break
        extra = min(2extra, max_band_extra)
    end

    diffs = IRDiff[]
    i = j = 0
    q = 1
    nmatch = nmis = ngap = nunk = 0
    while q <= length(ops)
        op = ops[q]
        if op == 0x01
            i += 1
            j += 1
            x, y = a[i], b[j]
            if !isacgt(x) || !isacgt(y)
                nunk += 1
            elseif x == y
                nmatch += 1
            else
                nmis += 1
                push!(diffs, IRDiff(:snp, :both, mod1(ir.irb[1] + i - 1, L), mod1(ir.ira[2] - j + 1, L),
                                    1, string(Char(x), '>', Char(y))))
            end
            q += 1
        else
            r = q
            while r <= length(ops) && ops[r] == op
                r += 1
            end
            len = r - q
            ngap += len
            if op == 0x02                    # bases in IRb only
                push!(diffs, IRDiff(:indel, :IRb, mod1(ir.irb[1] + i, L), mod1(ir.ira[2] - j + 1, L),
                                    len, String(Char.(a[(i + 1):(i + len)]))))
                i += len
            else                             # bases in IRa only
                push!(diffs, IRDiff(:indel, :IRa, mod1(ir.irb[1] + i - 1, L), mod1(ir.ira[2] - j, L),
                                    len, String(Char.(b[(j + 1):(j + len)]))))
                j += len
            end
            q = r
        end
    end
    sub = nmatch + nmis
    return (diffs = diffs, identity = sub == 0 ? NaN : nmatch / sub,
            aligned_frac = min(n, m) == 0 ? NaN : (sub + nunk) / min(n, m), n_gap = ngap,
            n_unknown = nunk, band = w, band_limited = edge)
end

ir_alignment(seq::AbstractString, ir::IRResult; kwargs...) = ir_alignment(seqbytes(seq), ir; kwargs...)

"IR copy differences (see `ir_alignment`)."
ir_differences(seq::AbstractString, ir::IRResult; kwargs...) = ir_alignment(seq, ir; kwargs...).diffs

# ------------------------------------------------------------------ regions and layout

"Junctions as the first base of the region to their right."
junctions(ir::IRResult) = (JLB = ir.irb[1], JSB = mod1(ir.irb[2] + 1, ir.L),
                           JSA = ir.ira[1], JLA = mod1(ir.ira[2] + 1, ir.L))

"Region of position p: :IRb, :IRa, :SSC or :LSC."
function region_of(ir::IRResult, p::Integer)
    p = mod1(p, ir.L)
    incirc(p, ir.irb) && return :IRb
    incirc(p, ir.ira) && return :IRa
    ir.ssc > 0 && incirc(p, (mod1(ir.irb[2] + 1, ir.L), mod1(ir.ira[1] - 1, ir.L))) && return :SSC
    return :LSC
end

"Per-base region codes: 0 LSC, 1 IRb, 2 SSC, 3 IRa."
function region_vector(ir::IRResult)
    v = zeros(UInt8, ir.L)
    for (r, code) in ((ir.irb, 0x01), (ir.ira, 0x03))
        len = mod(r[2] - r[1], ir.L) + 1
        for t in 0:(len - 1)
            v[mod1(r[1] + t, ir.L)] = code
        end
    end
    for t in 0:(ir.ssc - 1)
        v[mod1(ir.irb[2] + 1 + t, ir.L)] = 0x02
    end
    return v
end

"""
    needs_ssc_flip(ir, p) -> Bool or nothing

`p`: a position inside the gene that should open the SSC (e.g. the ndhF midpoint).
True when that position lies closer to JSA than to JSB, i.e. the SSC is in the other
orientation; `nothing` when p is not in the SSC.
"""
function needs_ssc_flip(ir::IRResult, p::Integer)
    region_of(ir, p) == :SSC || return nothing
    d_jsb = mod(p - ir.irb[2], ir.L)
    d_jsa = mod(ir.ira[1] - p, ir.L)
    return d_jsa < d_jsb
end

"""
    map_position(ir, p; flip_ssc=false) -> (position, strand_flipped)

Position of base p in the layout produced by `rotate_canonical` (LSC starting at 1).
"""
function map_position(ir::IRResult, p::Integer; flip_ssc::Bool=false)
    start = mod1(ir.ira[2] + 1, ir.L)
    np = mod1(p - start + 1, ir.L)
    if flip_ssc && region_of(ir, p) == :SSC
        s1 = ir.lsc + ir.len_irb + 1
        s2 = s1 + ir.ssc - 1
        return (s1 + s2 - np, true)
    end
    return (np, false)
end

"""
    rotate_canonical(seq, ir; flip_ssc=false) -> (sequence, IRResult)

Rotates the genome to LSC - IRb - SSC - IRa with the LSC starting at position 1;
`flip_ssc` reverse-complements the SSC (the other SSC orientation isomer).
Whole-genome orientation (strand) is left as given.
"""
function rotate_canonical(seq::AbstractString, ir::IRResult; flip_ssc::Bool=false)
    s = seqbytes(seq)
    L = ir.L
    start = mod1(ir.ira[2] + 1, L)
    r = vcat(s[start:end], s[1:(start - 1)])
    s1 = ir.lsc + ir.len_irb + 1
    s2 = s1 + ir.ssc - 1
    if flip_ssc && ir.ssc > 0
        r[s1:s2] = [compl(x) for x in reverse(r[s1:s2])]
    end
    irb = (ir.lsc + 1, ir.lsc + ir.len_irb)
    ira = (s2 + 1, s2 + ir.len_ira)
    return String(r), IRResult(L, irb, ira, ir.len_irb, ir.len_ira, ir.lsc, ir.ssc, ir.gc,
                               ir.drift, ir.identity, ir.aligned_frac, ir.n_unknown)
end

# ------------------------------------------------------------------ input helpers

"""
    foreach_fasta(f, path)

Calls f(id, sequence) for each record of a FASTA file (plain or .gz, see open_maybe_gz)
without holding the whole file in memory.
"""
function foreach_fasta(f, path::AbstractString)
    io = parentmodule(@__MODULE__).open_maybe_gz(path)
    try
        foreach_fasta(f, io)
    finally
        close(io)
    end
    return nothing
end

"""
    foreach_fasta(f, io::IO)

Same on an open stream.
"""
function foreach_fasta(f, io::IO)
    name = ""
    buf = IOBuffer()
    for line in eachline(io)
        if startswith(line, '>')
            isempty(name) || f(name, String(take!(buf)))
            name = String(first(split(line[2:end])))
        else
            print(buf, strip(line))
        end
    end
    isempty(name) || f(name, String(take!(buf)))
    return nothing
end

"Reads a FASTA file (plain or .gz) into id => sequence pairs."
function read_fasta(path::AbstractString)
    io = parentmodule(@__MODULE__).open_maybe_gz(path)
    out = Pair{String,String}[]
    name = ""
    buf = IOBuffer()
    try
        for line in eachline(io)
            if startswith(line, '>')
                isempty(name) || push!(out, name => String(take!(buf)))
                name = String(first(split(line[2:end])))
            else
                print(buf, strip(line))
            end
        end
        isempty(name) || push!(out, name => String(take!(buf)))
    finally
        close(io)
    end
    return out
end

"Parses 'a..b:+;c..d:-' (also ',' separators and 1/-1 strands)."
function parse_segments(s::AbstractString)
    segs = Tuple{Int,Int,Int8}[]
    for m in eachmatch(r"(\d+)(?:\.\.(\d+))?:([+-]|-?1)", s)
        a = parse(Int, m.captures[1])
        b = m.captures[2] === nothing ? a : parse(Int, m.captures[2])
        push!(segs, (a, b, startswith(m.captures[3], "-") ? Int8(-1) : Int8(1)))
    end
    return segs
end

"""
    read_features_tsv(path; accession=nothing, skip_masked=false, skip_status=("fail",))
        -> Dict{String,Vector{Feature}}

Reads labels.tsv (columns accession, feature, gene, segments, optional mask and status);
accession versions are ignored when matching. By default only features whose QC status
is in `skip_status` are left out: a training mask (e.g. a rare name) is not a reason to
hide a feature on a genome map.
"""
function read_features_tsv(path::AbstractString; accession=nothing, skip_masked::Bool=false,
                           skip_status=("fail",))
    key(a) = first(split(a, '.'))
    want = accession === nothing ? nothing : Set(key.(accession isa AbstractString ? [accession] : accession))
    out = Dict{String,Vector{Feature}}()
    open(path) do io
        hdr = split(readline(io), '\t')
        col(n) = something(findfirst(==(n), hdr), 0)
        ia, ift, ig, isg, im, ist = col("accession"), col("feature"), col("gene"), col("segments"),
                                    col("mask"), col("status")
        (ia == 0 || ift == 0 || ig == 0 || isg == 0) && error("$path: needs accession, feature, gene, segments")
        for line in eachline(io)
            f = split(line, '\t')
            k = key(f[ia])
            want !== nothing && !(k in want) && continue
            skip_masked && im > 0 && f[im] == "1" && continue
            ist > 0 && f[ist] in skip_status && continue
            push!(get!(out, k, Feature[]), Feature(String(f[ig]), String(f[ift]), parse_segments(f[isg])))
        end
    end
    return out
end

# ------------------------------------------------------------------ junction plot (SVG)

const FAMILY_COLORS = [
    ("trn", "#7b4fa0"), ("rrn", "#d1495b"), ("psa", "#2a6fbb"), ("psb", "#3a9a4a"),
    ("pet", "#8fc46a"), ("atp", "#e8862e"), ("ndh", "#e9b949"), ("rps", "#7fb7d9"),
    ("rpl", "#b59bcf"), ("rpo", "#ef8f8a"), ("ycf", "#9a5b2e"), ("chl", "#4fb3a9"),
    ("rbc", "#c9b400"), ("ccs", "#9e9ac8"), ("mat", "#5d8aa8"), ("clp", "#c46ba0"),
    ("acc", "#8c8c3c"), ("inf", "#6c9a8b"),
]
const OTHER_COLOR = "#8a8a8a"

function gene_family(name::AbstractString)
    n = lowercase(name)
    for (p, _) in FAMILY_COLORS
        startswith(n, p) && return p
    end
    return "other"
end

function gene_color(name::AbstractString)
    fam = gene_family(name)
    for (p, c) in FAMILY_COLORS
        p == fam && return c
    end
    return OTHER_COLOR
end

xml(s) = replace(String(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;", "\"" => "&quot;")
fmtbp(x::Integer) = replace(string(x), r"(?<=\d)(?=(\d{3})+$)" => ",")
textw(s::AbstractString, size::Real) = 0.58 * size * length(s)

# signed distance (bases, real) of genome point x from boundary point c on a circle of length L
cdist(x::Real, c::Real, L::Int) = mod(x - c + L / 2, L) - L / 2

"Pieces of features around a boundary: (name, strand, da, db), da/db in bases from it."
function junction_pieces(features::Vector{Feature}, jb::Int, L::Int, R::Real)
    c = jb - 0.5
    out = Tuple{String,Int8,Float64,Float64}[]
    for f in features, (a, b, st) in f.segs
        da = cdist(a - 0.5, c, L)
        db = da + (mod(b - a, L) + 1)
        (db > -R && da < R) && push!(out, (f.name, st, da, db))
    end
    return out
end

"Radius (bases) at which at least one feature piece reaches into the window."
function needed_radius(features::Vector{Feature}, jb::Int, L::Int)
    c = jb - 0.5
    best = Inf
    for f in features, (a, b, _) in f.segs
        da = cdist(a - 0.5, c, L)
        db = da + (mod(b - a, L) + 1)
        d = da < 0 < db ? 0.0 : min(abs(da), abs(db))
        best = min(best, d)
    end
    return best
end

niceup(x) = (p = 10.0^floor(log10(max(x, 1))); ceil(x / p * 2) / 2 * p)

"""
    plot_ir_junctions(tracks, file; width=1200, track_height=110, radius=nothing,
                      min_radius=200, max_radius=5000, title="IR junctions") -> file

IRscope-style junction plot written as SVG. Each track is a bar LSC | IRb | SSC | IRa | LSC
(segment widths fixed, lengths printed); around each junction (JLB, JSB, JSA, JLA)
features within ±radius bases are drawn to a common local scale, + strand above and
- strand below the bar. A feature crossing a junction gets the number of its bases on
each side; when none crosses, the distance to the nearest feature is given. IR copy
differences within the window are marked on the bar (circle SNP, triangle indel),
and their totals are printed in the middle of each IR. `radius` may be a 4-tuple of
bases (JLB, JSB, JSA, JLA); by default the smallest radius showing a feature in every
track is used, clamped to [min_radius, max_radius].
"""
function plot_ir_junctions(tracks::AbstractVector{IRTrack}, file::AbstractString;
                           width::Int=1200, track_height::Int=110, radius=nothing,
                           min_radius::Real=200, max_radius::Real=5000,
                           title::AbstractString="IR junctions")
    left, right, top, bottom = 210, 30, 78, 70
    P = width - left - right
    fr = (0.14, 0.24, 0.26, 0.24, 0.12)
    xs = cumsum([left; [f * P for f in fr]])          # 6 segment edges
    xJ = xs[2:5]                                      # JLB, JSB, JSA, JLA
    half = 0.085P
    nT = length(tracks)
    height = top + nT * track_height + bottom
    jnames = ("JLB", "JSB", "JSA", "JLA")

    jbs(t) = (j = junctions(t.ir); (j.JLB, j.JSB, j.JSA, j.JLA))
    R = if radius === nothing
        Rs = zeros(4)
        for t in tracks
            t.ir === nothing && continue
            jb = jbs(t)
            for q in 1:4
                Rs[q] = max(Rs[q], needed_radius(t.features, jb[q], t.L) + 60)
            end
        end
        Tuple(isfinite(Rs[q]) ? clamp(niceup(Rs[q]), min_radius, max_radius) : Float64(max_radius) for q in 1:4)
    else
        Tuple(Float64(r) for r in radius)
    end

    io = IOBuffer()
    println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" """,
            """viewBox="0 0 $width $height" font-family="Helvetica, Arial, sans-serif">""")
    println(io, """<rect width="$width" height="$height" fill="white"/>""")
    println(io, """<text x="$(width / 2)" y="24" font-size="15" font-weight="bold" text-anchor="middle">$(xml(title))</text>""")

    ytop = top - 8
    ybot = top + nT * track_height
    for q in 1:4
        println(io, """<line x1="$(xJ[q])" y1="$ytop" x2="$(xJ[q])" y2="$ybot" stroke="#444" stroke-width="0.8" stroke-dasharray="3,3"/>""")
        println(io, """<text x="$(xJ[q])" y="$(top - 30)" font-size="12" font-weight="bold" font-style="italic" text-anchor="middle">$(jnames[q])</text>""")
        println(io, """<text x="$(xJ[q])" y="$(top - 16)" font-size="9" fill="#555" text-anchor="middle">±$(fmtbp(round(Int, R[q]))) bp</text>""")
    end

    used = Set{String}()
    for (ti, t) in enumerate(tracks)
        yc = top + (ti - 0.5) * track_height
        println(io, """<text x="$(left - 12)" y="$(yc - 2)" font-size="11" font-weight="bold" font-style="italic" text-anchor="end">$(xml(t.name))</text>""")
        println(io, """<text x="$(left - 12)" y="$(yc + 11)" font-size="9.5" fill="#444" text-anchor="end">$(fmtbp(t.L)) bp</text>""")
        if t.ir === nothing
            println(io, """<rect x="$left" y="$(yc - 7)" width="$P" height="14" fill="#e6e6e6" stroke="#999" stroke-width="0.5"/>""")
            println(io, """<text x="$(left + P / 2)" y="$(yc + 3.5)" font-size="10" text-anchor="middle" fill="#555">no IR detected</text>""")
            continue
        end
        draw_track!(io, t, yc, xs, xJ, half, R, used)
    end

    # legend
    ly = height - bottom + 30
    lx = left
    for (p, c) in FAMILY_COLORS
        p in used || continue
        println(io, """<rect x="$lx" y="$(ly - 8)" width="12" height="9" fill="$c"/>""")
        println(io, """<text x="$(lx + 16)" y="$ly" font-size="9.5" font-style="italic">$p</text>""")
        lx += 52
    end
    if "other" in used
        println(io, """<rect x="$lx" y="$(ly - 8)" width="12" height="9" fill="$OTHER_COLOR"/>""")
        println(io, """<text x="$(lx + 16)" y="$ly" font-size="9.5">other</text>""")
        lx += 60
    end
    println(io, """<circle cx="$(lx + 10)" cy="$(ly - 3.5)" r="3" fill="#d62728" stroke="#222" stroke-width="0.4"/>""")
    println(io, """<text x="$(lx + 17)" y="$ly" font-size="9.5">SNP</text>""")
    lx += 50
    println(io, """<path d="M $(lx + 6) $(ly) L $(lx + 10) $(ly - 7) L $(lx + 14) $(ly) Z" fill="#2ca02c" stroke="#222" stroke-width="0.4"/>""")
    println(io, """<text x="$(lx + 18)" y="$ly" font-size="9.5">indel (IR copies)</text>""")
    println(io, "</svg>")
    write(file, String(take!(io)))
    return file
end

function draw_track!(io::IO, t::IRTrack, yc::Real, xs, xJ, half::Real, R, used::Set{String})
    ir = t.ir
    L = t.L
    bh = 14                                            # bar height
    segcol = ("#cfe3f5", "#f7c873", "#d9ecc9", "#f7c873", "#cfe3f5")
    seglab = ("LSC", "IRb", "SSC", "IRa", "LSC")
    seglen = (ir.lsc, ir.len_irb, ir.ssc, ir.len_ira, ir.lsc)
    for q in 1:5
        println(io, """<rect x="$(xs[q])" y="$(yc - bh / 2)" width="$(xs[q + 1] - xs[q])" height="$bh" fill="$(segcol[q])" stroke="#777" stroke-width="0.5"/>""")
        lab = q == 5 ? "LSC" : "$(seglab[q]) $(fmtbp(seglen[q])) bp"
        println(io, """<text x="$((xs[q] + xs[q + 1]) / 2)" y="$(yc + 3.5)" font-size="9" font-weight="bold" text-anchor="middle" fill="#333">$(xml(lab))</text>""")
    end

    # IR copy differences: totals in the middle of each IR, markers near junctions
    diffs = t.diffs === nothing ? IRDiff[] : t.diffs
    if t.diffs !== nothing
        nsnp = count(d -> d.kind == :snp, diffs)
        nind = count(d -> d.kind == :indel, diffs)
        nunk = ir.n_unknown
        msg = if isempty(diffs)
            nunk > 0 ? "no differences at comparable positions, $(nunk) unknown" :
            nunk == 0 ? "identical copies" : "no differences found"
        else
            "$(nsnp) SNP, $(nind) indel" * (nunk > 0 ? ", $(nunk) unknown" : "")
        end
        for q in (2, 4)
            println(io, """<text x="$((xs[q] + xs[q + 1]) / 2)" y="$(yc + bh / 2 + 13)" font-size="8.5" fill="#555" text-anchor="middle">$(msg)</text>""")
        end
    end
    j = junctions(ir)
    jb = (j.JLB, j.JSB, j.JSA, j.JLA)

    for q in 1:4
        scale = half / R[q]                            # px per base
        xj = xJ[q]
        # differences on the bar
        for d in diffs
            p = q <= 2 ? d.pos_irb : d.pos_ira
            off = cdist(p, jb[q] - 0.5, L)
            abs(off) > R[q] && continue
            x = xj + off * scale
            if d.kind == :snp
                println(io, """<circle cx="$x" cy="$yc" r="2.6" fill="#d62728" stroke="#222" stroke-width="0.4"/>""")
            else
                println(io, """<path d="M $(x - 3.2) $(yc + 3) L $x $(yc - 3.5) L $(x + 3.2) $(yc + 3) Z" fill="#2ca02c" stroke="#222" stroke-width="0.4"/>""")
            end
        end

        pieces = junction_pieces(t.features, jb[q], L, R[q])
        crossing = [p for p in pieces if p[3] < 0 < p[4]]
        # Vertical layout per side, as distance from the bar edge: feature box 3..13,
        # junction numbers (baseline 17 above / 24 below), names in two staggered rows
        # (baselines 29, 40 above / 36, 47 below).
        rows_up = [Tuple{Float64,Float64}[], Tuple{Float64,Float64}[]]
        rows_dn = [Tuple{Float64,Float64}[], Tuple{Float64,Float64}[]]
        for (name, st, da, db) in sort(pieces; by = p -> abs((p[3] + p[4]) / 2))
            push!(used, gene_family(name))
            x1 = max(xj + da * scale, xj - half)
            x2 = min(xj + db * scale, xj + half)
            x2 - x1 < 0.6 && (x2 = x1 + 0.6)
            up = st >= 0                               # + strand above the bar
            yb = up ? yc - bh / 2 - 13 : yc + bh / 2 + 3
            println(io, """<rect x="$x1" y="$yb" width="$(x2 - x1)" height="10" fill="$(gene_color(name))" stroke="#333" stroke-width="0.4"/>""")
            tw = textw(name, 9)
            xc = clamp((x1 + x2) / 2, xj - half + tw / 2, xj + half - tw / 2)
            span = (xc - tw / 2 - 2, xc + tw / 2 + 2)
            rows = up ? rows_up : rows_dn
            r = findfirst(row -> all(sp -> span[2] < sp[1] || span[1] > sp[2], row), rows)
            r === nothing && continue                  # no room: box drawn without a name
            push!(rows[r], span)
            ylab = up ? yc - bh / 2 - 18 - 11r : yc + bh / 2 + 25 + 11r
            println(io, """<text x="$xc" y="$ylab" font-size="9" font-style="italic" text-anchor="middle">$(xml(name))</text>""")
        end

        if !isempty(crossing)
            for (name, st, da, db) in crossing
                ynum = st >= 0 ? yc - bh / 2 - 17 : yc + bh / 2 + 24
                lbp, rbp = round(Int, -da), round(Int, db)
                println(io, """<text x="$(xj - 2)" y="$ynum" font-size="8.5" text-anchor="end">$(fmtbp(lbp)) bp</text>""")
                println(io, """<text x="$(xj + 2)" y="$ynum" font-size="8.5" text-anchor="start">$(fmtbp(rbp)) bp</text>""")
            end
        elseif !isempty(pieces)
            # nearest feature: distance from the junction
            name, st, da, db = first(sort(pieces; by = p -> min(abs(p[3]), abs(p[4]))))
            dist = da >= 0 ? da : -db
            edge = da >= 0 ? xj + min(da * scale, half) : xj + max(db * scale, -half)
            up = st >= 0
            yl = up ? yc - bh / 2 - 14 : yc + bh / 2 + 14
            println(io, """<line x1="$xj" y1="$yl" x2="$edge" y2="$yl" stroke="#222" stroke-width="0.7"/>""")
            println(io, """<line x1="$xj" y1="$(yl - 2.5)" x2="$xj" y2="$(yl + 2.5)" stroke="#222" stroke-width="0.7"/>""")
            println(io, """<line x1="$edge" y1="$(yl - 2.5)" x2="$edge" y2="$(yl + 2.5)" stroke="#222" stroke-width="0.7"/>""")
            ynum = up ? yl - 3 : yl + 9
            println(io, """<text x="$((xj + edge) / 2)" y="$ynum" font-size="8.5" text-anchor="middle">$(fmtbp(round(Int, dist))) bp</text>""")
        else
            println(io, """<text x="$xj" y="$(yc - bh / 2 - 18)" font-size="8.5" fill="#777" text-anchor="middle">no feature within ±$(fmtbp(round(Int, R[q]))) bp</text>""")
        end
    end
end

end # module

# ======================================================================
# PlastomeRef
# ======================================================================
"""
    PlastomeRef

Reference index for the plastome annotator: reference genomes with their labelled
features (from labels.tsv), MinHash sketches, and selection of references for a query.
Base Julia only.

Reference features are stored per part, numbered along the transcript as in Chloë
(`gene/CDS/1`, `gene/intron/2`, `gene/CDS/3`); an intron part is added between two
consecutive segments on the same strand less than `max_intron` apart (cis-splicing),
not across a trans-splicing junction (rps12). Copies of a gene in one reference
(e.g. the two IR copies) are numbered 1, 2, ... in order of their first segment; a
segment shared by two copies (the rps12 5' exon) is stored once.

Reference selection (`select_references`): MinHash (Mash) distances rank the
references; the list holds the `n_near` nearest (at most `max_per_genus` per genus,
the query's genus excluded when given) plus the `per_group` nearest of every group,
because at deep divergence (e.g. ferns against bryophytes) Mash distances saturate and
stop ranking groups correctly. The annotator re-ranks this shortlist by alignment
coverage.
"""
module PlastomeRef

export RefIndex, RefRecord, RefPart, sketch, mash_distance, build_ref_parts,
       write_index, load_index, read_sequence, select_references

const ENC = let t = fill(0xff, 256)
    for (c, v) in zip("ACGT", 0x00:0x03)
        t[Int(c) + 1] = v
    end
    t
end

# ------------------------------------------------------------------ MinHash

@inline function splitmix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

"""
    sketch(seq; k=21, s=1000) -> Vector{UInt64}

Bottom-s MinHash sketch (sorted) of the canonical k-mers of a circular sequence;
k-mers with a non-ACGT base are skipped.
"""
function sketch(seq::AbstractString; k::Int=21, s::Int=1000)
    1 <= k <= 32 || throw(ArgumentError("k must be in 1:32"))
    b = codeunits(uppercase(String(seq)))
    L = length(b)
    mask = k == 32 ? typemax(UInt64) : (UInt64(1) << (2k)) - 1
    shift = 2 * (k - 1)
    v = r = UInt64(0)
    valid = 0
    hs = Vector{UInt64}(undef, 0)
    sizehint!(hs, L)
    for p in 1:(L + k - 1)
        c = ENC[Int(b[mod1(p, L)]) + 1]
        if c == 0xff
            valid = 0
            v = r = UInt64(0)
            continue
        end
        v = ((v << 2) | UInt64(c)) & mask
        r = (r >> 2) | (UInt64(3 - c) << shift)
        valid += 1
        st = p - k + 1
        if valid >= k && 1 <= st <= L
            push!(hs, splitmix64(min(v, r)))
        end
    end
    unique!(sort!(hs))
    return length(hs) > s ? hs[1:s] : hs
end

"""
    mash_distance(a, b; k=21, s=1000) -> Float64

Mash distance from two sorted bottom-s sketches: Jaccard estimate j over the s smallest
hashes of the union, D = -log(2j / (1 + j)) / k; 1.0 when no hash is shared.
"""
function mash_distance(a::AbstractVector{UInt64}, b::AbstractVector{UInt64}; k::Int=21, s::Int=1000)
    i = j = 1
    seen = shared = 0
    while seen < s && i <= length(a) && j <= length(b)
        if a[i] == b[j]
            shared += 1
            i += 1
            j += 1
        elseif a[i] < b[j]
            i += 1
        else
            j += 1
        end
        seen += 1
    end
    (seen == 0 || shared == 0) && return 1.0
    jac = shared / seen
    return -log(2jac / (1 + jac)) / k
end

# ------------------------------------------------------------------ reference parts

"One part of a reference feature (exon, intron, tRNA, rRNA); coordinates 1-based, first <= last."
struct RefPart
    gene::String
    copy::Int
    kind::String       # CDS, tRNA, rRNA, intron
    part::Int          # position along the transcript (1, 2, 3, ...)
    strand::Char
    first::Int
    last::Int
end

path(p::RefPart) = string(p.gene, '/', p.kind, '/', p.part)

const SEG_RE = r"(\d+)\.\.(\d+):([+-])"

"""
    build_ref_parts(features; max_intron=4000) -> Vector{RefPart}

`features`: (gene, type, segments) of one reference, segments as "a..b:+;c..d:-" in
transcript order. Segments of one feature become parts 1, 3, 5, ...; an intron part
(2, 4, ...) is added between consecutive segments on the same strand separated by
1..max_intron bases.
"""
function build_ref_parts(features; max_intron::Int=4000)
    parts = RefPart[]
    ncopy = Dict{String,Int}()
    order = sortperm([first_start(f[3]) for f in features])
    for oi in order
        gene, kind, segs = features[oi]
        segv = [(parse(Int, m[1]), parse(Int, m[2]), m[3][1]) for m in eachmatch(SEG_RE, segs)]
        isempty(segv) && continue
        c = ncopy[gene] = get(ncopy, gene, 0) + 1
        n = 0
        for (q, (a, b, st)) in enumerate(segv)
            if q > 1
                pa, pb, pst = segv[q - 1]
                n += 1
                if pst == st
                    # 5' -> 3' gap between consecutive exons on the strand
                    g0, g1 = st == '+' ? (pb + 1, a - 1) : (b + 1, pa - 1)
                    if 1 <= g1 - g0 + 1 <= max_intron
                        push!(parts, RefPart(gene, c, "intron", n, st, g0, g1))
                    end
                end
            end
            n += 1
            # a segment shared by two copies (rps12 5' exon) is kept once
            dup = any(p -> p.gene == gene && p.kind == kind && p.part == n && p.strand == st &&
                           p.first == min(a, b) && p.last == max(a, b), parts)
            dup || push!(parts, RefPart(gene, c, kind, n, st, min(a, b), max(a, b)))
        end
    end
    return parts
end

first_start(segs::AbstractString) = (m = match(SEG_RE, segs); m === nothing ? typemax(Int) : parse(Int, m[1]))

# ------------------------------------------------------------------ index

struct RefRecord
    accession::String
    group::String
    genus::String
    split::String
    structure::String
    length::Int
    offset::Int        # byte offset of the sequence line in refs.fasta
end

struct RefIndex
    dir::String
    k::Int
    s::Int
    records::Vector{RefRecord}
    sketches::Vector{Vector{UInt64}}
    parts::Dict{String,Vector{RefPart}}           # accession -> parts
    templates::Dict{String,Tuple{Int,Float64}}    # path -> (n references, median length)
end

"""
    write_index(dir, records, sketches, parts; k, s)

Writes refs.tsv, sketches.bin, parts.tsv and templates.tsv (references carrying each
part path and its median length) into `dir`; refs.fasta (one line per sequence, the
byte offset of each sequence line in `records`) is written by the caller.
"""
function write_index(dir::AbstractString, records::Vector{RefRecord},
                     sketches::Vector{Vector{UInt64}}, parts::Dict{String,Vector{RefPart}};
                     k::Int, s::Int)
    mkpath(dir)
    open(joinpath(dir, "refs.tsv"), "w") do io
        println(io, join(["accession", "group", "genus", "split", "structure", "length", "offset",
                          "n_parts"], '\t'))
        for r in records
            println(io, join((r.accession, r.group, r.genus, r.split, r.structure, r.length, r.offset,
                              length(get(parts, r.accession, RefPart[]))), '\t'))
        end
    end
    open(joinpath(dir, "sketches.bin"), "w") do io
        write(io, Int64(length(sketches)), Int64(k), Int64(s))
        for sk in sketches
            write(io, Int64(length(sk)))
            write(io, sk)
        end
    end
    lens = Dict{String,Vector{Int}}()
    nref = Dict{String,Set{String}}()
    open(joinpath(dir, "parts.tsv"), "w") do io
        println(io, join(["accession", "gene", "copy", "kind", "part", "strand", "first", "last"], '\t'))
        for r in records, p in get(parts, r.accession, RefPart[])
            println(io, join((r.accession, p.gene, p.copy, p.kind, p.part, p.strand, p.first, p.last), '\t'))
            push!(get!(() -> Int[], lens, path(p)), p.last - p.first + 1)
            push!(get!(() -> Set{String}(), nref, path(p)), r.accession)
        end
    end
    open(joinpath(dir, "templates.tsv"), "w") do io
        println(io, "path\tn_refs\tmedian_length")
        for pth in sort!(collect(keys(lens)))
            v = sort(lens[pth])
            m = isodd(length(v)) ? v[(end + 1) ÷ 2] : (v[end ÷ 2] + v[end ÷ 2 + 1]) / 2
            println(io, pth, '\t', length(nref[pth]), '\t', m)
        end
    end
    return nothing
end

function load_index(dir::AbstractString)
    records = RefRecord[]
    open(joinpath(dir, "refs.tsv")) do io
        readline(io)
        for ln in eachline(io)
            f = split(ln, '\t')
            push!(records, RefRecord(f[1], f[2], f[3], f[4], f[5], parse(Int, f[6]), parse(Int, f[7])))
        end
    end
    sketches = Vector{Vector{UInt64}}()
    k = s = 0
    open(joinpath(dir, "sketches.bin")) do io
        n = read(io, Int64)
        k = Int(read(io, Int64))
        s = Int(read(io, Int64))
        for _ in 1:n
            m = read(io, Int64)
            push!(sketches, read!(io, Vector{UInt64}(undef, m)))
        end
    end
    length(sketches) == length(records) || error("sketches.bin and refs.tsv differ in length")
    parts = Dict{String,Vector{RefPart}}()
    open(joinpath(dir, "parts.tsv")) do io
        readline(io)
        for ln in eachline(io)
            f = split(ln, '\t')
            push!(get!(() -> RefPart[], parts, f[1]),
                  RefPart(f[2], parse(Int, f[3]), f[4], parse(Int, f[5]), f[6][1], parse(Int, f[7]), parse(Int, f[8])))
        end
    end
    templates = Dict{String,Tuple{Int,Float64}}()
    open(joinpath(dir, "templates.tsv")) do io
        readline(io)
        for ln in eachline(io)
            f = split(ln, '\t')
            templates[f[1]] = (parse(Int, f[2]), parse(Float64, f[3]))
        end
    end
    return RefIndex(String(dir), k, s, records, sketches, parts, templates)
end

"Sequence of reference `i` (read from refs.fasta at its offset)."
function read_sequence(idx::RefIndex, i::Integer)
    open(joinpath(idx.dir, "refs.fasta")) do io
        seek(io, idx.records[i].offset)
        return readline(io)
    end
end

"""
    select_references(idx, qsketch; n_near=24, per_group=2, max_per_genus=2,
                      exclude_genus="", exclude_accession="") -> Vector{Tuple{Int,Float64}}

Shortlist of (reference index, Mash distance), nearest first: the `n_near` nearest
references with at most `max_per_genus` per genus, plus the `per_group` nearest of each
group not yet represented that often. References of `exclude_genus` (and
`exclude_accession`) are never chosen.
"""
function select_references(idx::RefIndex, qsketch::Vector{UInt64}; n_near::Int=24, per_group::Int=2,
                           max_per_genus::Int=2, exclude_genus::AbstractString="",
                           exclude_accession::AbstractString="")
    n = length(idx.records)
    d = Vector{Float64}(undef, n)
    Threads.@threads for i in 1:n
        d[i] = mash_distance(qsketch, idx.sketches[i]; k=idx.k, s=idx.s)
    end
    ord = sortperm(d)
    chosen = Tuple{Int,Float64}[]
    taken = falses(n)
    pergenus = Dict{String,Int}()
    pergroup = Dict{String,Int}()
    ok(r) = !(r.genus == exclude_genus && !isempty(exclude_genus)) &&
            !(r.accession == exclude_accession && !isempty(exclude_accession))
    for i in ord
        length(chosen) >= n_near && break
        r = idx.records[i]
        ok(r) || continue
        get(pergenus, r.genus, 0) >= max_per_genus && continue
        push!(chosen, (i, d[i]))
        taken[i] = true
        pergenus[r.genus] = get(pergenus, r.genus, 0) + 1
        pergroup[r.group] = get(pergroup, r.group, 0) + 1
    end
    for i in ord
        r = idx.records[i]
        (taken[i] || !ok(r)) && continue
        get(pergroup, r.group, 0) >= per_group && continue
        get(pergenus, r.genus, 0) >= max_per_genus && continue
        push!(chosen, (i, d[i]))
        taken[i] = true
        pergenus[r.genus] = get(pergenus, r.genus, 0) + 1
        pergroup[r.group] = get(pergroup, r.group, 0) + 1
    end
    sort!(chosen; by=last)
    return chosen
end

end # module

# ======================================================================
# PlastomeAnnot
# ======================================================================
"""
    PlastomeAnnot

Annotation transfer for plastomes (baseline, step 2): alignment of reference genomes to
a query, transfer of reference feature parts through the aligned blocks, feature
stacks, candidate parts and their ends. Base Julia only. Coordinates are 0-based
inside the module; `Candidate` reports 1-based `first`/`last`.

Alignment (`align`): exact matches seeded by k-mers (k = 12) of the query, extended to
maximal exact matches (MEMs). A MEM is significant when it is at least as long as the
shortest match with P(chance) < 0.1 for the two sequence lengths (Chloë's
matchLengthThreshold), or when a second MEM >= k lies on a nearby diagonal (within
`pair_diag`) within `pair_dist` bases (pairs of short matches are rarely random).
Significant MEMs are chained when consecutive in the reference with gaps <= max_gap in
both sequences differing by <= max_skew; gaps are filled with MEMs of the local
significance threshold (colinear subset), and blocks separated by equal gaps <=
merge_gap (substitutions only) are merged, as in Chloë.

Transfer (`transfer`): every reference part overlapping a block gives a piece on the
query (strand flipped for blocks on the reverse complement), with the projected 5' and
3' ends of the whole part and Chloë's weights: coverage x ((length - offset)/length)^2.

Candidates (`candidates`): pieces of one part path and strand are clustered; the window
of the template median length with the largest coverage-weighted depth gives the
candidate, depth = window sum / (median length x summed coverage of the shortlisted
references carrying that path); ends = weighted mode of the projected ends.

`refine_cds`: baseline rules for single-exon CDS (start codon near the projected 5'
end, ATG preferred over ACG/GTG; 3' end at the first in-frame stop if the ORF length is
within 20% of the projected length).
"""
module PlastomeAnnot

export Piece, Candidate, kmer_index, align, coverage, transfer, candidates, refine_cds,
       match_threshold, revcomp_bytes, transfer_query, mode_ends

const ENC = let t = fill(0xff, 256)
    for (c, v) in zip("ACGT", 0x00:0x03)
        t[Int(c) + 1] = v
    end
    t
end
const COMPL = let t = fill(UInt8('N'), 256)
    for (a, b) in zip("ACGTN", "TGCAN")
        t[Int(a) + 1] = UInt8(b)
    end
    t
end
@inline isacgt(b::UInt8) = ENC[Int(b) + 1] != 0xff
revcomp_bytes(s::Vector{UInt8}) = [COMPL[Int(s[i]) + 1] for i in length(s):-1:1]
@inline at(s::Vector{UInt8}, p::Int) = s[mod(p, length(s)) + 1]

function prob_match(m::Int, n::Int, k::Int)
    (m < k || n < k) && return 0.0
    return 1 - ((1 - 4.0^-k)^(m - k + 1))^(n - k + 1)
end

"Shortest exact-match length with P(chance match) < 0.1 between sequences of lengths m, n."
function match_threshold(m::Int, n::Int)
    for k in 1:25
        prob_match(m, n, k) < 0.1 && return max(k, 2)
    end
    return 26
end

"k-mer (k <= 15) -> 0-based start positions on a circular sequence."
function kmer_index(q::Vector{UInt8}, k::Int=12)
    1 <= k <= 15 || throw(ArgumentError("k must be in 1:15"))
    L = length(q)
    idx = Dict{UInt32,Vector{Int32}}()
    mask = (UInt32(1) << (2k)) - UInt32(1)
    v = UInt32(0)
    valid = 0
    for p in 0:(L + k - 2)
        c = ENC[Int(at(q, p)) + 1]
        if c == 0xff
            valid = 0
            v = UInt32(0)
            continue
        end
        v = ((v << 2) | UInt32(c)) & mask
        valid += 1
        st = p - k + 1
        if valid >= k && 0 <= st < L
            push!(get!(Vector{Int32}, idx, v), Int32(st))
        end
    end
    return idx
end

"""
MEMs (r, q, len) with r in [r0, r1) and q in [q0, q1) (unrolled coordinates, taken
modulo the sequence lengths), length >= minlen, seeded by the k-mers of `qidx`.
"""
function find_mems(r::Vector{UInt8}, q::Vector{UInt8}, qidx::Dict{UInt32,Vector{Int32}}, k::Int,
                   r0::Int, r1::Int, q0::Int, q1::Int, minlen::Int)
    Lr, Lq = length(r), length(q)
    mask = (UInt32(1) << (2k)) - UInt32(1)
    done = Dict{Int,Int}()
    out = NTuple{3,Int}[]
    v = UInt32(0)
    valid = 0
    for p in r0:(r1 + k - 2)
        c = ENC[Int(at(r, p)) + 1]
        if c == 0xff
            valid = 0
            v = UInt32(0)
            continue
        end
        v = ((v << 2) | UInt32(c)) & mask
        valid += 1
        x = p - k + 1
        (valid < k || x < r0 || x + k > r1) && continue
        ys = get(qidx, v, nothing)
        ys === nothing && continue
        for y in ys
            yy = q0 + mod(Int(y) - q0, Lq)
            yy + k > q1 && continue
            d = mod(yy - x, Lq)
            get(done, d, -1) > x && continue
            a, b = x, yy
            while a > r0 && b > q0
                rc = at(r, a - 1)
                (rc == at(q, b - 1) && isacgt(rc)) || break
                a -= 1
                b -= 1
            end
            e, f = x + k, yy + k
            while e < r1 && f < q1
                rc = at(r, e)
                (rc == at(q, f) && isacgt(rc)) || break
                e += 1
                f += 1
            end
            done[d] = e
            e - a >= minlen && push!(out, (a, b, e - a))
        end
    end
    return out
end

"Longest chain (by matched bases) of MEMs increasing in both sequences."
function colinear(mems::Vector{NTuple{3,Int}})
    isempty(mems) && return mems
    ms = sort(mems)
    n = length(ms)
    best = [m[3] for m in ms]
    prev = zeros(Int, n)
    for i in 1:n
        r, q, l = ms[i]
        for j in 1:(i - 1)
            rj, qj, lj = ms[j]
            if rj + lj <= r && qj + lj <= q && best[j] + l > best[i]
                best[i] = best[j] + l
                prev[i] = j
            end
        end
    end
    i = argmax(best)
    ch = NTuple{3,Int}[]
    while i > 0
        push!(ch, ms[i])
        i = prev[i]
    end
    return reverse!(ch)
end

"""
    align(r, q, qidx; k=12, max_gap=3000, max_skew=1000, merge_gap=400, min_fill=20,
          pair_dist=150, pair_diag=8) -> Vector{NTuple{3,Int}}

Blocks (r_start, q_start, length) of the forward reference on one strand of the query
(`q`, indexed by `qidx`); 0-based, q_start in [0, length(q)).
"""
function align(r::Vector{UInt8}, q::Vector{UInt8}, qidx::Dict{UInt32,Vector{Int32}}; k::Int=12,
               max_gap::Int=3000, max_skew::Int=1000, merge_gap::Int=400, min_fill::Int=20,
               pair_dist::Int=150, pair_diag::Int=8)
    Lr, Lq = length(r), length(q)
    thr = match_threshold(Lr, Lq)
    raw = sort!(find_mems(r, q, qidx, k, 0, Lr, 0, 2Lq, k))
    keep = NTuple{3,Int}[]
    for i in eachindex(raw)
        a, b, l = raw[i]
        ok = l >= thr
        j = i - 1
        while !ok && j >= 1 && a - raw[j][1] <= pair_dist
            if abs((b - a) - (raw[j][2] - raw[j][1])) <= pair_diag && raw[j][1] + raw[j][3] <= a
                ok = true
            end
            j -= 1
        end
        j = i + 1
        while !ok && j <= length(raw) && raw[j][1] - a <= pair_dist
            if abs((b - a) - (raw[j][2] - raw[j][1])) <= pair_diag && a + l <= raw[j][1]
                ok = true
            end
            j += 1
        end
        ok && push!(keep, (a, b, l))
    end
    mems = sort!(unique!([(a, mod(b, Lq), l) for (a, b, l) in keep]))

    chains = Vector{Vector{NTuple{3,Int}}}()
    for m in mems
        a, b, l = m
        bestc, bestg = 0, typemax(Int)
        for (ci, ch) in enumerate(chains)
            ra, qb, ll = ch[end]
            gr = a - (ra + ll)
            gq = mod(b - (qb + ll), Lq)
            gq > Lq ÷ 2 && (gq -= Lq)
            if -k < gr <= max_gap && -k < gq <= max_gap && abs(gr - gq) <= max_skew && gr < bestg
                bestc, bestg = ci, gr
            end
        end
        bestc == 0 ? push!(chains, [m]) : push!(chains[bestc], m)
    end

    blocks = NTuple{3,Int}[]
    for ch in chains
        filled = NTuple{3,Int}[ch[1]]
        for m in ch[2:end]
            pa, pb, pl = filled[end]
            a, b, l = m
            g = mod(b - (pb + pl), Lq)
            b_un = pb + pl + (g <= Lq ÷ 2 ? g : g - Lq)
            gr, gq = a - (pa + pl), b_un - (pb + pl)
            if gr >= min_fill && gq >= min_fill && !(gr == gq && gr <= merge_gap)
                t = max(match_threshold(gr, gq), 8)
                if t >= k
                    for s in colinear(find_mems(r, q, qidx, k, pa + pl, a, pb + pl, b_un, t))
                        push!(filled, s)
                    end
                end
            end
            push!(filled, (a, b_un, l))
        end
        merged = [collect(filled[1])]
        for (a, b, l) in filled[2:end]
            pa, pb, pl = merged[end]
            gr, gq = a - (pa + pl), b - (pb + pl)
            if gr == gq && 0 <= gr <= merge_gap
                merged[end][3] = a + l - pa
            elseif gr < 0 || gq < 0
                cut = max(-gr, -gq)
                l - cut > 0 && push!(merged, [a + cut, b + cut, l - cut])
            else
                push!(merged, [a, b, l])
            end
        end
        for m in merged
            push!(blocks, (m[1], mod(m[2], Lq), m[3]))
        end
    end
    return blocks
end

"Fraction of the query covered by forward blocks and reverse-complement blocks."
function coverage(bf::Vector{NTuple{3,Int}}, br::Vector{NTuple{3,Int}}, Lq::Int)
    cov = falses(Lq)
    for (_, b, l) in bf, t in 0:(l - 1)
        cov[mod(b + t, Lq) + 1] = true
    end
    for (_, b, l) in br
        bf0 = Lq - 1 - mod(b + l - 1, Lq)          # forward start of the rc block
        for t in 0:(l - 1)
            cov[mod(bf0 + t, Lq) + 1] = true
        end
    end
    return count(cov) / Lq
end

struct Piece
    path::String
    strand::Char
    qstart::Int        # 0-based, query forward coordinates
    qlen::Int
    q5::Int            # projected 5' end of the whole part
    w5::Float64
    q3::Int
    w3::Float64
    ref::Int
end

"""
    transfer(parts, blocks, orient, Lr, Lq, cov, ref) -> Vector{Piece}

`parts`: (gene, kind, part, strand, first, last) 0-based reference coordinates;
`orient` 1 for blocks on the query forward strand, -1 for blocks on its reverse
complement (q coordinates of the reverse complement).
"""
function transfer(parts, blocks::Vector{NTuple{3,Int}}, orient::Int, Lr::Int, Lq::Int, cov::Float64, ref::Int)
    out = Piece[]
    for (gene, kind, part, st, first, last) in parts
        plen = last - first + 1
        pth = string(gene, '/', kind, '/', part)
        for (rs, qs, ln) in blocks, shift in (0, Lr)
            f, l = first + shift, last + shift
            o1, o2 = max(f, rs), min(l, rs + ln - 1)
            o1 > o2 && continue
            m(x) = orient == 1 ? mod(qs + (x - rs), Lq) : Lq - 1 - mod(qs + (x - rs), Lq)
            qa, qb, sq = orient == 1 ? (m(o1), m(o2), st) : (m(o2), m(o1), st == '+' ? '-' : '+')
            r5, r3 = st == '+' ? (f, l) : (l, f)
            off5 = st == '+' ? o1 - f : l - o2
            off3 = st == '+' ? l - o2 : o1 - f
            w5 = cov * ((plen - off5) / plen)^2
            w3 = cov * ((plen - off3) / plen)^2
            push!(out, Piece(pth, sq, qa, mod(qb - qa, Lq) + 1, m(r5), w5, m(r3), w3, ref))
        end
    end
    return out
end

struct Candidate
    gene::String
    kind::String
    part::Int
    strand::Char
    first::Int         # 1-based
    last::Int
    depth::Float64
    nrefs::Int
    w5frac::Float64    # weight of the chosen 5' end / all 5' votes
    w3frac::Float64
    nvote5::Int        # distinct projected 5' ends
    nvote3::Int
    refs::Vector{Int}  # supporting references
    near::Bool         # ends taken from the nearest reference (candidates(near_ref))
    alt_first::Int     # ends from the weighted mode of all references (= first/last unless near)
    alt_last::Int
end

Candidate(gene, kind, part, strand, first, last, depth, nrefs, w5frac, w3frac, nvote5, nvote3, refs) =
    Candidate(gene, kind, part, strand, first, last, depth, nrefs, w5frac, w3frac, nvote5, nvote3, refs, false,
              first, last)

Candidate(c::Candidate, first::Int, last::Int) =
    Candidate(c.gene, c.kind, c.part, c.strand, first, last, c.depth, c.nrefs, c.w5frac, c.w3frac,
              c.nvote5, c.nvote3, c.refs, c.near, c.alt_first, c.alt_last)

"The candidate with the weighted-mode ends of all references instead of the nearest reference's."
mode_ends(c::Candidate) =
    Candidate(c.gene, c.kind, c.part, c.strand, c.alt_first, c.alt_last, c.depth, c.nrefs, c.w5frac, c.w3frac,
              c.nvote5, c.nvote3, c.refs, false, c.alt_first, c.alt_last)

"Weighted mode of positions -> (mode, weight fraction at the mode, distinct positions)."
function weighted_mode(pos::Vector{Int}, w::Vector{Float64})
    acc = Dict{Int,Float64}()
    for (p, x) in zip(pos, w)
        acc[p] = get(acc, p, 0.0) + x
    end
    best, bw = 0, -1.0
    for p in sort!(collect(keys(acc)))
        if acc[p] > bw
            best, bw = p, acc[p]
        end
    end
    tot = sum(values(acc))
    return best, (tot > 0 ? bw / tot : 0.0), length(acc)
end

"""
    candidates(pieces, templates, path_weight, refw, Lq; min_depth=0.25, near_ref=0) -> Vector{Candidate}

`templates`: path -> (n_refs, median length); `path_weight`: path -> summed coverage of
the shortlisted references carrying the path; `refw`: reference -> coverage.
`near_ref` > 0: when that reference projects both ends of the part fully (its end lies
inside an aligned block), the candidate takes its ends instead of the weighted mode
(near = true); the vote statistics (w5frac, nvote5, ...) stay those of all references.
"""
function candidates(pieces::Vector{Piece}, templates::Dict{String,Tuple{Int,Float64}},
                    path_weight::Dict{String,Float64}, refw::Dict{Int,Float64}, Lq::Int;
                    min_depth::Float64=0.25, near_ref::Int=0)
    by = Dict{Tuple{String,Char},Vector{Piece}}()
    for p in pieces
        push!(get!(Vector{Piece}, by, (p.path, p.strand)), p)
    end
    out = Candidate[]
    for ((pth, sq), ps) in by
        m = max(1, round(Int, get(templates, pth, (0, 100.0))[2]))
        nr = max(1e-9, get(path_weight, pth, 1.0))
        sort!(ps; by=p -> p.qstart)
        clusters = Vector{Vector{Piece}}()
        cur = Piece[]
        cend = -1
        for p in ps
            if !isempty(cur) && p.qstart > cend + m
                push!(clusters, cur)
                cur = Piece[]
            end
            push!(cur, p)
            pe = p.qstart + p.qlen - 1
            cend = length(cur) > 1 ? max(cend, pe) : pe
        end
        isempty(cur) || push!(clusters, cur)
        gene, kind, part = rsplit(pth, '/'; limit=3)
        for cl in clusters
            c0 = minimum(p.qstart for p in cl)
            c1 = maximum(p.qstart + p.qlen - 1 for p in cl)
            span = c1 - c0 + 1
            cov = zeros(span + m)
            for p in cl
                w = get(refw, p.ref, 1.0)
                for t in 0:(p.qlen - 1)
                    cov[p.qstart - c0 + t + 1] += w
                end
            end
            s = sum(@view cov[1:m])
            best, bpos = s, 0
            for i in 1:(span - 1)
                s += cov[i + m] - cov[i]
                if s > best
                    best, bpos = s, i
                end
            end
            depth = best / (m * nr)
            depth < min_depth && continue
            ws, we = c0 + bpos, c0 + bpos + m - 1
            inw = [p for p in cl if p.qstart <= we && p.qstart + p.qlen - 1 >= ws]
            e5, f5, n5 = weighted_mode([p.q5 for p in inw], [p.w5 for p in inw])
            e3, f3, n3 = weighted_mode([p.q3 for p in inw], [p.w3 for p in inw])
            af, al = sq == '+' ? (e5, e3) : (e3, e5)
            near = false
            if near_ref > 0
                wfull = 0.999 * get(refw, near_ref, Inf)
                a5 = [p for p in inw if p.ref == near_ref && p.w5 >= wfull]
                a3 = [p for p in inw if p.ref == near_ref && p.w3 >= wfull]
                if !isempty(a5) && !isempty(a3)
                    e5 = weighted_mode([p.q5 for p in a5], [p.w5 for p in a5])[1]
                    e3 = weighted_mode([p.q3 for p in a3], [p.w3 for p in a3])[1]
                    near = true
                end
            end
            f, l = sq == '+' ? (e5, e3) : (e3, e5)
            rs = sort!(unique(p.ref for p in inw))
            push!(out, Candidate(String(gene), String(kind), parse(Int, part), sq, mod(f, Lq) + 1, mod(l, Lq) + 1,
                                 depth, length(rs), f5, f3, n5, n3, rs, near, mod(af, Lq) + 1, mod(al, Lq) + 1))
        end
    end
    return out
end

"""
    transfer_query(q, refs; min_depth=0.25, templates) -> (candidates, coverage by reference)

`q`: query bytes (uppercase); `refs`: vector of (reference id, sequence bytes, parts) with
parts (gene, kind, part, strand, first, last) in 1-based reference coordinates. Aligns
every reference to both query strands, transfers its parts and builds candidates (depth
normalised by the summed coverage of the copies of each path in the references).
"""
function transfer_query(q::Vector{UInt8}, refs, templates::Dict{String,Tuple{Int,Float64}}; min_depth::Float64=0.25,
                        near_ref::Int=0)
    Lq = length(q)
    qrc = revcomp_bytes(q)
    idxf, idxr = kmer_index(q), kmer_index(qrc)
    pieces = Piece[]
    refw = Dict{Int,Float64}()
    pw = Dict{String,Float64}()
    for (ri, r, parts) in refs
        bf = align(r, q, idxf)
        br = align(r, qrc, idxr)
        c = coverage(bf, br, Lq)
        refw[ri] = c
        rp = [(p[1], p[2], p[3], p[4], p[5] - 1, p[6] - 1) for p in parts]
        append!(pieces, transfer(rp, bf, 1, length(r), Lq, c, ri))
        append!(pieces, transfer(rp, br, -1, length(r), Lq, c, ri))
        for p in rp                       # every copy counts (IR genes have two)
            pth = string(p[1], '/', p[2], '/', p[3])
            pw[pth] = get(pw, pth, 0.0) + c
        end
    end
    return candidates(pieces, templates, pw, refw, Lq; min_depth=min_depth, near_ref=near_ref), refw, qrc
end

const STOPS = ("TAA", "TAG", "TGA")
const STARTS = ("ATG", "ACG", "GTG")

"""
    refine_cds(c, q, qrc, single_exon_genes; window=30) -> Candidate

Baseline rules for single-exon CDS: 5' end to the nearest in-frame start codon within
+/- window of the projection (first ATG by distance; otherwise the nearest ACG/GTG);
3' end to the first in-frame stop from the 5' end if the ORF length is within 20% of
the projected length.
"""
function refine_cds(c::Candidate, q::Vector{UInt8}, qrc::Vector{UInt8}, single_exon_genes::Set{String};
                    window::Int=30)
    (c.kind == "CDS" && c.gene in single_exon_genes) || return c
    Lq = length(q)
    fwd = c.strand == '+'
    s = fwd ? q : qrc
    a = fwd ? c.first - 1 : Lq - c.last           # 0-based 5' end on the strand
    b = fwd ? c.last - 1 : Lq - c.first           # 0-based 3' end
    plen = mod(b - a, Lq) + 1
    cod(p) = String([at(s, p), at(s, p + 1), at(s, p + 2)])
    best = nothing
    for d in sort(collect(-window:3:window); by=abs)
        x = cod(a + d)
        if x == "ATG"
            best = a + d
            break
        end
        (x in STARTS && best === nothing) && (best = a + d)
    end
    best === nothing || (a = best)
    e = nothing
    for t in 0:3:(floor(Int, plen * 1.2) + 2)
        if cod(a + t) in STOPS
            e = a + t + 2
            break
        end
    end
    (e !== nothing && 0.8plen <= e - a + 1 <= 1.2plen) && (b = e)
    f, l = fwd ? (mod(a, Lq) + 1, mod(b, Lq) + 1) : (mod(Lq - 1 - b, Lq) + 1, mod(Lq - 1 - a, Lq) + 1)
    return Candidate(c, f, l)
end

end # module

# ======================================================================
# PlastomeFeatures
# ======================================================================
"""
    PlastomeFeatures

Candidate features shared by build_training.jl (training table) and annotate.jl
(inference): the same code computes them in both places. Base Julia only.
"""
module PlastomeFeatures

export FEATURES, candidate_features, overlap, ilen, cdist, anticodon_dna, strand_seq,
       region_code, trna_aa, name_conflict

const FEATURES = [
    "depth", "nrefs", "nrefs_frac", "rel_depth", "rank_in_path", "n_cand_path",
    "w5frac", "w3frac", "nvote5", "nvote3", "supp_cov", "max_cov", "supp_mash_min",
    "len_ratio", "template_frac", "part", "last_part",
    "start_codon", "stop_ok", "inframe_stops", "len_mod3", "anticodon_found",
    "gc", "region", "dist_junction"]

segs1(a, b, L) = a <= b ? ((a, b),) : ((a, L), (1, b))
function overlap(a0, a1, b0, b1, L)
    s = 0
    for (x0, x1) in segs1(a0, a1, L), (y0, y1) in segs1(b0, b1, L)
        s += max(0, min(x1, y1) - max(x0, y0) + 1)
    end
    return s
end
ilen(a, b, L) = mod(b - a, L) + 1
cdist(a, b, L) = mod(b - a + L ÷ 2, L) - L ÷ 2

const COMPL = Dict('A' => 'T', 'C' => 'G', 'G' => 'C', 'T' => 'A', 'U' => 'A')

"Expected anticodon (DNA, 5'->3') from a tRNA name such as trnG-UCC; nothing if absent."
function anticodon_dna(gene::AbstractString)
    m = match(r"^trn[A-Za-z]+-([ACGU]{3})$", gene)
    m === nothing ? nothing : replace(String(m[1]), 'U' => 'T')
end

"Bases of a 1-based circular interval on the given strand (5'->3')."
function strand_seq(q::Vector{UInt8}, qrc::Vector{UInt8}, first::Int, last::Int, strand::Char)
    L = length(q)
    n = ilen(first, last, L)
    if strand == '+'
        return String([q[mod1(first + t, L)] for t in 0:(n - 1)])
    else
        a = L - last + 1                         # position of the 5' end on the reverse complement
        return String([qrc[mod1(a + t, L)] for t in 0:(n - 1)])
    end
end

"Region code of a 1-based position: 0 LSC, 1 IRb, 2 SSC, 3 IRa (NaN without an IR)."
function region_code(ir, p, L)
    ir === nothing && return NaN
    b0, b1, a0, a1 = ir
    inr(x, y) = x <= y ? (x <= p <= y) : (p >= x || p <= y)
    inr(b0, b1) && return 1.0
    inr(a0, a1) && return 3.0
    inr(mod1(b1 + 1, L), mod1(a0 - 1, L)) && return 2.0
    return 0.0
end

"""
    candidate_features(c, cands, q, qrc, L, refw, mash, templates, nref_total, parts_by_path, ir)
        -> Vector{Float64} (order of FEATURES)
"""
function candidate_features(c, cands, q, qrc, L, refw, mash, templates, nref_total, parts_by_path, ir)
    pth = "$(c.gene)/$(c.kind)/$(c.part)"
    same_path = [x for x in cands if x.gene == c.gene && x.kind == c.kind && x.part == c.part]
    maxd = maximum(x.depth for x in same_path)
    rank = count(x -> x.depth > c.depth, same_path) + 1
    nsl = get(parts_by_path, pth, 0)
    clen = ilen(c.first, c.last, L)
    tn, tm = get(templates, pth, (0, NaN))
    nextpart = get(templates, "$(c.gene)/$(c.kind)/$(c.part + 2)", (0, 0.0))[1]
    lastp = nextpart < 0.5 * max(tn, 1)
    s = strand_seq(q, qrc, c.first, c.last, c.strand)
    gc = count(ch -> ch == 'G' || ch == 'C', s) / max(1, length(s))
    startc = stopok = instops = mod3 = acf = NaN
    if c.kind == "CDS"
        if c.part == 1 && length(s) >= 3
            st = s[1:3]
            startc = st == "ATG" ? 1.0 : st == "GTG" ? 2.0 : st == "ACG" ? 3.0 : 0.0
        end
        if lastp && length(s) >= 3
            stopok = s[end-2:end] in ("TAA", "TAG", "TGA") ? 1.0 : 0.0
        end
        if c.part == 1 && lastp
            mod3 = length(s) % 3 == 0 ? 1.0 : 0.0
            instops = Float64(count(i -> s[i:i+2] in ("TAA", "TAG", "TGA"), 1:3:(length(s) - 5)))
        end
    elseif c.kind == "tRNA" && c.part == 1 && lastp
        ac = anticodon_dna(c.gene)
        if ac !== nothing && length(s) >= 45
            acf = occursin(ac, s[28:45]) ? 1.0 : 0.0
        end
    end
    reg = region_code(ir, c.first, L)
    dj = NaN
    if ir !== nothing
        b0, b1, a0, a1 = ir
        dj = Float64(minimum(abs(cdist(p, j, L)) for p in (c.first, c.last) for j in (b0, b1, a0, a1)))
    end
    supp = [refw[r] for r in c.refs if haskey(refw, r)]
    feats = [c.depth, c.nrefs, nsl > 0 ? c.nrefs / nsl : NaN, maxd > 0 ? c.depth / maxd : NaN, rank,
             length(same_path), c.w5frac, c.w3frac, c.nvote5, c.nvote3,
             isempty(supp) ? NaN : sum(supp) / length(supp), maximum(values(refw); init=0.0),
             isempty(c.refs) ? NaN : minimum(get(mash, r, 1.0) for r in c.refs),
             isnan(tm) ? NaN : clen / tm, tn / max(1, nref_total), c.part, lastp ? 1.0 : 0.0,
             startc, stopok, instops, mod3, acf, gc, reg, dj]

    return feats
end

trna_aa(g) = (m = match(r"^trn(fM|[A-Z])", g); m === nothing ? "" : (m[1] == "fM" ? "M" : String(m[1])))

"""
Possible naming conflict between a candidate and an overlapping label of another gene:
both rRNAs, or both tRNAs with the same amino acid or the same anticodon (trnG-UCC
annotated as trnG-GCC; trnI-CAU / trnM-CAU / trnfM-CAU).
"""
function name_conflict(c, t)
    (c.kind == t.kind && c.kind == "rRNA") && return true
    if c.kind == "tRNA" && t.kind == "tRNA"
        a1, a2 = trna_aa(c.gene), trna_aa(t.gene)
        (!isempty(a1) && a1 == a2) && return true
        x1, x2 = anticodon_dna(c.gene), anticodon_dna(t.gene)
        (x1 !== nothing && x1 == x2) && return true
    end
    return false
end

end # module

# ======================================================================
# PlastomeTrees
# ======================================================================
"""
    PlastomeTrees

Gradient-boosted decision trees for binary classification (logistic loss), histogram
based, Base Julia only. Features are binned at up to `max_bins - 1` quantile cut
points (bin 0 = missing/NaN, which always goes to the left child); trees are grown
level-wise to `max_depth`; a split maximises GL^2/(HL+λ) + GR^2/(HR+λ) - G^2/(H+λ)
with at least `min_h` hessian on each side; leaf value -G/(H+λ) x learning rate. Each
tree uses a row subsample and a column subsample; early stopping on a validation set
(log loss, `patience` trees). Importance = summed split gain per feature.
"""
module PlastomeTrees

export TreeModel, fit_gbdt, predict_gbdt, save_model, load_model, auc

struct Tree
    feat::Vector{Int}      # 0 = leaf
    thr::Vector{Int}       # left child: bin <= thr
    left::Vector{Int}
    right::Vector{Int}
    val::Vector{Float64}
end

struct TreeModel
    features::Vector{String}
    cuts::Vector{Vector{Float64}}
    base::Float64
    trees::Vector{Tree}
    importance::Vector{Float64}
end

# deterministic RNG (splitmix64), no Random dependency
mutable struct SplitMix
    s::UInt64
end
function nextu!(r::SplitMix)
    r.s += 0x9e3779b97f4a7c15
    z = r.s
    z = (z ⊻ (z >> 30)) * 0xbf58476d1ce4e5b9
    z = (z ⊻ (z >> 27)) * 0x94d049bb133111eb
    return z ⊻ (z >> 31)
end
randf!(r::SplitMix) = (nextu!(r) >> 11) * (1.0 / 9007199254740992.0)

"Linear-interpolation quantile of a sorted vector (as numpy's default)."
function qsorted(v::Vector{Float64}, q::Float64)
    n = length(v)
    pos = q * (n - 1)
    lo = floor(Int, pos)
    fr = pos - lo
    return lo + 1 >= n ? v[n] : v[lo + 1] + fr * (v[lo + 2] - v[lo + 1])
end

function make_cuts(X::Matrix{Float64}, max_bins::Int)
    cuts = Vector{Vector{Float64}}()
    for j in 1:size(X, 2)
        v = sort!([x for x in view(X, :, j) if !isnan(x)])
        if isempty(v)
            push!(cuts, Float64[])
            continue
        end
        u = unique(v)
        if length(u) <= max_bins - 1
            push!(cuts, [(u[i] + u[i + 1]) / 2 for i in 1:(length(u) - 1)])
        else
            qs = [qsorted(v, t / (max_bins - 1)) for t in 1:(max_bins - 2)]
            push!(cuts, unique(qs))
        end
    end
    return cuts
end

"Bin matrix: 0 for NaN, otherwise 1 + number of cut points below x."
function binize(X::Matrix{Float64}, cuts::Vector{Vector{Float64}})
    n, p = size(X)
    B = zeros(UInt8, n, p)
    for j in 1:p
        c = cuts[j]
        for i in 1:n
            x = X[i, j]
            B[i, j] = isnan(x) ? 0x00 : UInt8(searchsortedfirst(c, x))
        end
    end
    return B
end

function build_tree(B::Matrix{UInt8}, g::Vector{Float64}, h::Vector{Float64}, rows::Vector{Int},
                    nbins::Vector{Int}, max_depth::Int, lam::Float64, min_h::Float64, cols::Vector{Int},
                    lr::Float64, imp::Vector{Float64})
    feat, thr, left, right, val = Int[], Int[], Int[], Int[], Float64[]
    function new_node()
        push!(feat, 0); push!(thr, 0); push!(left, 0); push!(right, 0); push!(val, 0.0)
        return length(feat)
    end
    root = new_node()
    level = [(root, rows)]
    for depth in 0:max_depth
        nxt = Tuple{Int,Vector{Int}}[]
        for (node, r) in level
            G = sum(g[i] for i in r; init=0.0)
            H = sum(h[i] for i in r; init=0.0)
            val[node] = -G / (H + lam) * lr
            (depth == max_depth || H < 2min_h) && continue
            parent = G^2 / (H + lam)
            bg = zeros(length(cols))
            bt = zeros(Int, length(cols))
            # every variable assigned in the threaded body is declared local: a name also
            # assigned in the enclosing function would be shared by all threads
            Threads.@threads for ci in eachindex(cols)
                local jc = cols[ci]
                local nb = nbins[jc]
                local hg = zeros(nb)
                local hh = zeros(nb)
                for i in r
                    local bb = Int(B[i, jc]) + 1
                    hg[bb] += g[i]
                    hh[bb] += h[i]
                end
                local GL = 0.0
                local HL = 0.0
                local bestg = 0.0
                local bthr = -1
                for t in 0:(nb - 2)
                    GL += hg[t + 1]
                    HL += hh[t + 1]
                    local GR = G - GL
                    local HR = H - HL
                    (HL >= min_h && HR >= min_h) || continue
                    local gain = GL^2 / (HL + lam) + GR^2 / (HR + lam) - parent
                    if gain > bestg
                        bestg, bthr = gain, t
                    end
                end
                bg[ci], bt[ci] = bestg, bthr
            end
            kbest = argmax(bg)
            (bg[kbest] > 0 && bt[kbest] >= 0) || continue
            jbest, tbest = cols[kbest], bt[kbest]
            feat[node], thr[node] = jbest, tbest
            imp[jbest] += bg[kbest]
            rl = [i for i in r if B[i, jbest] <= tbest]
            rr = [i for i in r if B[i, jbest] > tbest]
            l = new_node()
            rn = new_node()
            left[node], right[node] = l, rn
            push!(nxt, (l, rl))
            push!(nxt, (rn, rr))
        end
        level = nxt
    end
    return Tree(feat, thr, left, right, val)
end

function tree_value(t::Tree, B::Matrix{UInt8}, i::Int)
    node = 1
    while t.feat[node] != 0
        node = B[i, t.feat[node]] <= t.thr[node] ? t.left[node] : t.right[node]
    end
    return t.val[node]
end

function add_tree!(F::Vector{Float64}, t::Tree, B::Matrix{UInt8})
    Threads.@threads for i in eachindex(F)
        F[i] += tree_value(t, B, i)
    end
end

sigmoid(x) = 1 / (1 + exp(-x))

function logloss(y::Vector{Float64}, F::Vector{Float64})
    s = 0.0
    for i in eachindex(y)
        p = clamp(sigmoid(F[i]), 1e-12, 1 - 1e-12)
        s -= y[i] * log(p) + (1 - y[i]) * log(1 - p)
    end
    return s / length(y)
end

"""
    fit_gbdt(X, y, features; Xv=nothing, yv=nothing, n_trees=400, lr=0.1, max_depth=6,
             lam=1.0, min_h=1.0, subsample=0.8, colsample=0.8, max_bins=64, patience=30,
             seed=1) -> TreeModel

`X`: samples x features (NaN = missing); `y`: 0/1.
"""
function fit_gbdt(X::Matrix{Float64}, y::Vector{Float64}, features::Vector{String};
                  Xv=nothing, yv=nothing, n_trees::Int=400, lr::Float64=0.1, max_depth::Int=6,
                  lam::Float64=1.0, min_h::Float64=1.0, subsample::Float64=0.8, colsample::Float64=0.8,
                  max_bins::Int=64, patience::Int=30, seed::Int=1, verbose::Bool=true)
    rng = SplitMix(UInt64(seed))
    cuts = make_cuts(X, max_bins)
    B = binize(X, cuts)
    nbins = [length(c) + 2 for c in cuts]
    p0 = clamp(sum(y) / length(y), 1e-6, 1 - 1e-6)
    base = log(p0 / (1 - p0))
    F = fill(base, length(y))
    hasv = Xv !== nothing
    Bv = hasv ? binize(Xv, cuts) : zeros(UInt8, 0, 0)
    Fv = hasv ? fill(base, length(yv)) : Float64[]
    trees = Tree[]
    imp = zeros(size(X, 2))
    best_it, best_loss = 0, Inf
    g = similar(y)
    h = similar(y)
    for it in 1:n_trees
        for i in eachindex(y)
            p = sigmoid(F[i])
            g[i] = p - y[i]
            h[i] = p * (1 - p)
        end
        rows = [i for i in eachindex(y) if randf!(rng) < subsample]
        cols = [j for j in 1:size(X, 2) if randf!(rng) < colsample]
        isempty(cols) && (cols = [Int(nextu!(rng) % size(X, 2)) + 1])
        t = build_tree(B, g, h, rows, nbins, max_depth, lam, min_h, cols, lr, imp)
        push!(trees, t)
        add_tree!(F, t, B)
        if hasv
            add_tree!(Fv, t, Bv)
            ll = logloss(yv, Fv)
            if ll < best_loss - 1e-6
                best_loss, best_it = ll, it
            elseif it - best_it >= patience
                break
            end
            verbose && it % 25 == 0 && @info "tree $it: validation log loss $(round(ll; digits=5)) (best $(round(best_loss; digits=5)) at $best_it)"
        end
    end
    hasv && (trees = trees[1:best_it])
    return TreeModel(features, cuts, base, trees, imp)
end

"Probabilities for the rows of X."
function predict_gbdt(m::TreeModel, X::Matrix{Float64})
    B = binize(X, m.cuts)
    F = fill(m.base, size(X, 1))
    for t in m.trees
        add_tree!(F, t, B)
    end
    return sigmoid.(F)
end

function save_model(path::AbstractString, m::TreeModel)
    open(path, "w") do io
        println(io, "features\t", join(m.features, '\t'))
        println(io, "base\t", m.base)
        println(io, "importance\t", join(m.importance, '\t'))
        for c in m.cuts
            println(io, "cuts\t", join(c, '\t'))
        end
        for t in m.trees
            println(io, "tree\t", length(t.feat))
            for k in eachindex(t.feat)
                println(io, t.feat[k], '\t', t.thr[k], '\t', t.left[k], '\t', t.right[k], '\t', t.val[k])
            end
        end
    end
end

function load_model(path::AbstractString)
    lines = readlines(path)
    features = String.(split(lines[1], '\t')[2:end])
    base = parse(Float64, split(lines[2], '\t')[2])
    imp = parse.(Float64, split(lines[3], '\t')[2:end])
    cuts = Vector{Vector{Float64}}()
    trees = Tree[]
    i = 4
    while i <= length(lines)
        f = split(lines[i], '\t')
        if f[1] == "cuts"
            push!(cuts, [parse(Float64, x) for x in f[2:end] if !isempty(x)])
            i += 1
        elseif f[1] == "tree"
            n = parse(Int, f[2])
            ft, th, le, ri, va = Int[], Int[], Int[], Int[], Float64[]
            for k in 1:n
                r = split(lines[i + k], '\t')
                push!(ft, parse(Int, r[1])); push!(th, parse(Int, r[2])); push!(le, parse(Int, r[3]))
                push!(ri, parse(Int, r[4])); push!(va, parse(Float64, r[5]))
            end
            push!(trees, Tree(ft, th, le, ri, va))
            i += n + 1
        else
            error("bad model line $i")
        end
    end
    return TreeModel(features, cuts, base, trees, imp)
end

"Area under the ROC curve (ties averaged)."
function auc(y::Vector{Float64}, p::Vector{Float64})
    o = sortperm(p)
    ranks = zeros(length(p))
    i = 1
    while i <= length(o)
        j = i
        while j < length(o) && p[o[j + 1]] == p[o[i]]
            j += 1
        end
        for k in i:j
            ranks[o[k]] = (i + j) / 2
        end
        i = j + 1
    end
    npos = sum(y)
    nneg = length(y) - npos
    (npos == 0 || nneg == 0) && return NaN
    return (sum(ranks[k] for k in eachindex(y) if y[k] == 1) - npos * (npos + 1) / 2) / (npos * nneg)
end

end # module

# ======================================================================
# PlastomeProtein
# ======================================================================
"""
    PlastomeProtein

Protein-level search of reference proteins in the six-frame translation of a plastome,
Base Julia only. Used for genes that nucleotide transfer misses in divergent lineages.

- Translation: standard code (plastid code 11 has the same codon assignments).
- Seeds: exact amino-acid 3-mers of the reference protein in an index of the six
  frames; hits binned by diagonal (16 aa); bins with >= `min_seeds` hits are aligned.
- Alignment: local Smith-Waterman with affine gaps (Gotoh), BLOSUM62, gap open 11 /
  extend 1 (first gap residue costs 12), on a window of the frame around the diagonal.
- Significance: E = K m n exp(-lambda S) with the gapped BLOSUM62 11/1 parameters
  lambda = 0.267, K = 0.041; m = protein length, n = total six-frame length.
"""
module PlastomeProtein

export translate_nt, SixFrames, six_frames, ProteinHit, protein_search, protein_search_window, frame_nt

const AA = "ARNDCQEGHILKMFPSTWYV"
# BLOSUM62 (rows/columns in AA order), from Biopython Bio.Align.substitution_matrices
const B62 = Int8[
  4 -1 -2 -2  0 -1 -1  0 -2 -1 -1 -1 -1 -2 -1  1  0 -3 -2  0;
 -1  5  0 -2 -3  1  0 -2  0 -3 -2  2 -1 -3 -2 -1 -1 -3 -2 -3;
 -2  0  6  1 -3  0  0  0  1 -3 -3  0 -2 -3 -2  1  0 -4 -2 -3;
 -2 -2  1  6 -3  0  2 -1 -1 -3 -4 -1 -3 -3 -1  0 -1 -4 -3 -3;
  0 -3 -3 -3  9 -3 -4 -3 -3 -1 -1 -3 -1 -2 -3 -1 -1 -2 -2 -1;
 -1  1  0  0 -3  5  2 -2  0 -3 -2  1  0 -3 -1  0 -1 -2 -1 -2;
 -1  0  0  2 -4  2  5 -2  0 -3 -3  1 -2 -3 -1  0 -1 -3 -2 -2;
  0 -2  0 -1 -3 -2 -2  6 -2 -4 -4 -2 -3 -3 -2  0 -2 -2 -3 -3;
 -2  0  1 -1 -3  0  0 -2  8 -3 -3 -1 -2 -1 -2 -1 -2 -2  2 -3;
 -1 -3 -3 -3 -1 -3 -3 -4 -3  4  2 -3  1  0 -3 -2 -1 -3 -1  3;
 -1 -2 -3 -4 -1 -2 -3 -4 -3  2  4 -2  2  0 -3 -2 -1 -2 -1  1;
 -1  2  0 -1 -3  1  1 -2 -1 -3 -2  5 -1 -3 -1  0 -1 -3 -2 -2;
 -1 -1 -2 -3 -1  0 -2 -3 -2  1  2 -1  5  0 -2 -1 -1 -1 -1  1;
 -2 -3 -3 -3 -2 -3 -3 -3 -1  0  0 -3  0  6 -4 -2 -2  1  3 -1;
 -1 -2 -2 -1 -3 -1 -1 -2 -2 -3 -3 -1 -2 -4  7 -1 -1 -4 -3 -2;
  1 -1  1  0 -1  0  0  0 -1 -2 -2  0 -1 -2 -1  4  1 -3 -2 -2;
  0 -1  0 -1 -1 -1 -1 -2 -2 -1 -1 -1 -1 -2 -1  1  5 -2 -2  0;
 -3 -3 -4 -4 -2 -2 -3 -2 -2 -3 -2 -3 -1  1 -4 -3 -2 11  2 -3;
 -2 -2 -2 -3 -2 -1 -2 -3  2 -1 -1 -2 -1  3 -3 -2 -2  2  7 -1;
  0 -3 -3 -3 -1 -2 -2 -3 -3  3  1 -2  1 -1 -2 -2  0 -3 -1  4]
const AIDX = let t = zeros(Int, 256)
    for (i, c) in enumerate(AA)
        t[Int(c) + 1] = i
    end
    t
end
const STAR = UInt8('*')
const LAMBDA = 0.267
const KAPPA = 0.041

@inline function score(a::UInt8, b::UInt8)
    i, j = AIDX[a + 1], AIDX[b + 1]
    (i > 0 && j > 0) && return Int(B62[i, j])
    return (a == STAR && b == STAR) ? 1 : -4
end

# standard genetic code, bases in TCAG order
const CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"
const NIDX = let t = fill(-1, 256)
    for (i, c) in enumerate("TCAG")
        t[Int(c) + 1] = i - 1
    end
    t
end

"Translation of the codons starting at 1-based `first` (`n` codons) of `s`; X for ambiguous codons."
function translate_nt(s::AbstractVector{UInt8}, first::Int, n::Int)
    out = Vector{UInt8}(undef, max(n, 0))
    for k in 1:n
        p = first + 3(k - 1)
        a, b, c = NIDX[s[p] + 1], NIDX[s[p + 1] + 1], NIDX[s[p + 2] + 1]
        out[k] = (a < 0 || b < 0 || c < 0) ? UInt8('X') : UInt8(CODE[16a + 4b + c + 1])
    end
    return out
end

"Six frames of a circular sequence: frames 1-3 on the + strand (offsets 0-2), 4-6 on the reverse complement."
struct SixFrames
    L::Int
    ext::Int
    prot::Vector{Vector{UInt8}}
    index::Vector{Vector{Int64}}      # 3-mer key + 1 -> packed (frame << 32 | aa position)
end

"Nucleotide strand position (1-based, may exceed L; use mod1) of amino acid `j` of frame `fi`."
frame_nt(fi::Int, j::Int) = mod(fi - 1, 3) + 3(j - 1) + 1

kmer_key(a::UInt8, b::UInt8, c::UInt8) =
    (AIDX[a + 1] > 0 && AIDX[b + 1] > 0 && AIDX[c + 1] > 0) ?
    (AIDX[a + 1] - 1) * 400 + (AIDX[b + 1] - 1) * 20 + AIDX[c + 1] : 0

function six_frames(q::Vector{UInt8}, qrc::Vector{UInt8}; ext::Int=40000)
    L = length(q)
    e = min(L, ext)
    prot = Vector{Vector{UInt8}}()
    for (S, f) in ((q, 0), (q, 1), (q, 2), (qrc, 0), (qrc, 1), (qrc, 2))
        Sx = vcat(S, S[1:e])
        n = (length(Sx) - f) ÷ 3
        push!(prot, translate_nt(Sx, f + 1, n))
    end
    index = [Int64[] for _ in 1:8000]
    for (fi, T) in enumerate(prot), j in 1:(length(T) - 2)
        k = kmer_key(T[j], T[j + 1], T[j + 2])
        k > 0 && push!(index[k], (Int64(fi) << 32) | j)
    end
    return SixFrames(L, e, prot, index)
end

struct ProteinHit
    frame::Int          # 1-3 + strand, 4-6 - strand
    score::Int
    evalue::Float64
    rstart::Int         # aligned protein range
    rend::Int
    tstart::Int         # aligned range in the frame translation (aa positions)
    tend::Int
    m::Int              # protein length
end

const NEG = -1_000_000_000

"Local alignment (Gotoh) of protein P against T[lo:hi] -> (score, rstart, rend, tstart, tend)."
function sw(P::Vector{UInt8}, T::Vector{UInt8}, lo::Int, hi::Int; go::Int=11, ge::Int=1)
    m = length(P)
    gopen = go + ge
    # row i of P is stored at index i + 1; index 1 is the empty-prefix boundary (always 0)
    H = zeros(Int, m + 1); Hi = zeros(Int, m + 1); Hj = zeros(Int, m + 1)
    E = fill(NEG, m + 1); Ei = zeros(Int, m + 1); Ej = zeros(Int, m + 1)
    best, bri, brj, bi, bj = 0, 0, 0, 0, 0
    for j in lo:hi
        hd, hdi, hdj = 0, 0, 0                  # old H of row i-1
        F, Fi, Fj = NEG, 0, 0
        for i in 1:m
            k = i + 1
            # horizontal gap (old H[k] -> this column)
            eo, ee = H[k] - gopen, E[k] - ge
            if eo >= ee
                E[k], Ei[k], Ej[k] = eo, Hi[k], Hj[k]
            else
                E[k] = ee
            end
            # vertical gap (new H of row i-1 -> row i)
            fo, fe = H[k - 1] - gopen, F - ge
            if i > 1 && fo >= fe
                F, Fi, Fj = fo, Hi[k - 1], Hj[k - 1]
            else
                F = fe
            end
            d = hd + score(P[i], T[j])
            di, dj = hd > 0 ? (hdi, hdj) : (i, j)
            hold, holdi, holdj = H[k], Hi[k], Hj[k]
            if d >= E[k] && d >= F && d > 0
                H[k], Hi[k], Hj[k] = d, di, dj
            elseif E[k] >= F && E[k] > 0
                H[k], Hi[k], Hj[k] = E[k], Ei[k], Ej[k]
            elseif F > 0
                H[k], Hi[k], Hj[k] = F, Fi, Fj
            else
                H[k], Hi[k], Hj[k] = 0, 0, 0
            end
            if H[k] > best
                best, bri, brj, bi, bj = H[k], Hi[k], Hj[k], i, j
            end
            hd, hdi, hdj = hold, holdi, holdj
        end
    end
    return best, bri, bi, brj, bj
end

"""
    protein_search(sf, P; min_seeds=3, max_windows=4, flank=40) -> ProteinHit or nothing

Best local alignment of protein `P` (terminal stop removed) in the six frames.
"""
function protein_search(sf::SixFrames, P::Vector{UInt8}; min_seeds::Int=3, max_windows::Int=4, flank::Int=40)
    m = length(P)
    m >= 10 || return nothing
    bins = Dict{Tuple{Int,Int},Int}()
    for i in 1:(m - 2)
        k = kmer_key(P[i], P[i + 1], P[i + 2])
        k > 0 || continue
        for x in sf.index[k]
            fi, j = Int(x >> 32), Int(x & 0xffffffff)
            key = (fi, fld(j - i, 16))
            bins[key] = get(bins, key, 0) + 1
        end
    end
    isempty(bins) && return nothing
    top = sort!([(c, k) for (k, c) in bins if c >= min_seeds]; rev=true)
    n = 2 * sf.L                        # six frames of L/3 amino acids
    best = nothing
    for (_, (fi, b)) in first(top, max_windows)
        T = sf.prot[fi]
        lo = max(1, 16b - flank)
        hi = min(length(T), 16b + 16 + m + flank)
        s, rs, re, ts, te = sw(P, T, lo, hi)
        s > 0 || continue
        ev = KAPPA * m * n * exp(-LAMBDA * s)
        if best === nothing || s > best.score
            best = ProteinHit(fi, s, ev, rs, re, ts, te, m)
        end
    end
    return best
end

"""
    protein_search_window(sf, P, glo, W) -> ProteinHit or nothing

Best local alignment of `P` in the six frames restricted to the genomic window of `W` nt
starting at `glo` (1-based, circular; `W` <= the frame extension). No seeding: the whole
window is aligned in each frame; the E-value uses the window's six-frame length (2W).
"""
function protein_search_window(sf::SixFrames, P::Vector{UInt8}, glo::Int, W::Int)
    m = length(P)
    (m >= 8 && 0 < W <= sf.ext) || return nothing
    L = sf.L
    s0m = L - mod1(glo + W - 1, L) + 1        # 5' end of the window on the reverse strand
    best = nothing
    for fi in 1:6
        T = sf.prot[fi]
        f = mod(fi - 1, 3)
        st = fi <= 3 ? glo : s0m
        lo = max(1, cld(st - f - 1, 3) + 1)              # first codon starting at >= st
        hi = min(length(T), fld(st + W - 1 - f, 3))      # last codon ending at <= st + W - 1
        hi > lo || continue
        sc, rs, re, ts, te = sw(P, T, lo, hi)
        sc > 0 || continue
        ev = KAPPA * m * 2W * exp(-LAMBDA * sc)
        if best === nothing || sc > best.score
            best = ProteinHit(fi, sc, ev, rs, re, ts, te, m)
        end
    end
    return best
end

end # module

# ======================================================================
# annotator
# ======================================================================
#= annotate (annotator) -----------------------------------------------------

annotate - plastome annotation: transfer from references, existence models, gene
models with ORF checks; GFF3 output, or end-to-end evaluation on a dataset split.

  annotate FASTA files:
    julia --threads auto PlastomeAnnotator.jl annotate --fasta genome.fa [--fasta more.fa ...] [--out DIR]
                                     [--index DIR] [--models DIR] [--threshold 0.5]
                                     [--exclude-genus GENUS]
    defaults: DATASET/refindex_all and DATASET/models_final when present (production),
    otherwise refindex_train and models; --exclude-genus drops a genus from the
    references (to test on a genome whose genus is in the index)
  evaluate on held-out records (query genus excluded from the references):
    julia --threads auto PlastomeAnnotator.jl annotate --eval validation [--max-queries N]
                                     [--exclude genus|accession] [--near-ref 0.01]
  --near-ref D: when the nearest reference is closer than Mash D, transferred parts take
  their ends from it instead of the weighted mode of all references (single-exon CDS
  ends are then not moved by the start/stop rules, and passing multi-exon models keep
  its junctions); 0 disables. --exclude accession (eval) keeps the query's congeners in
  the references, as in production.

Pipeline per genome:
  1. IR (PlastomeIR.find_ir_status) for region features.
  2. MinHash shortlist of references, alignment, transfer, candidates (depth >= 0.05),
     baseline CDS start/stop rules (PlastomeAnnot).
  3. Features (PlastomeFeatures) and existence probability: CDS and tRNA models
     (models/existence_<kind>.model); rRNA by rule (depth >= 0.25 and >= 0.5 of the
     deepest candidate of its path), as the rRNA model had too few negatives.
  4. Accepted parts (p >= threshold); overlapping parts of one path, and same-kind parts
     of different genes overlapping >= 50% of the shorter on one strand (paralogous
     tRNAs, name conflicts), keep the higher probability.
  5. Gene models: parts of one gene and strand chained along the transcript (part
     numbers increasing, gap -20..4000 nt); a trans-splicing junction starts a new model.
  6. CDS checks on the spliced model: start (ATG/GTG/ACG) when the model holds part 1;
     stop and length divisible by 3 when it holds part 1 and the last part; in-frame
     internal stops. Editing tolerance (editing_evidence) when one of the 3 best-covering
     references has >= 3 complete annotated CDS with in-frame internal stops or ending
     with CAA/CAG/CGA: internal stops tolerated as in qc_cds.py (<= 5 or <= 0.035 per
     codon), CAA/CAG/CGA count as terminal stops (C-to-U editing).
     Complete multi-exon models: junctions refined (refine_junctions: donor/acceptor
     moves <= 9 nt chosen by gene-specific junction profiles from the shortlisted
     references, junction_pwms; then the 3' end: cut at the first in-frame stop of the
     last exon, or extended to a stop <= 30 nt). Complete single-exon models failing the
     checks: ORF search (refine_orf). Before the checks, missing short exons are
     searched with the references' sequences and trans-spliced pieces are joined
     (build_models). A complete model still failing is reported as pseudo,
     except ycf15 and ycf68, always reported as CDS (flag orf_checks_ignored).
  7. Protein-level rescue (protein_rescue, PlastomeProtein): CDS genes present in >= half
     of the shortlisted references but without a model are searched with the references'
     proteins (single-exon genes) or exon peptides (multi-exon genes) in the window
     between the gene's syntenic neighbours (Smith-Waterman, BLOSUM62); single-exon genes
     also in the whole six frames (E <= 1e-3 either way); models flagged protein_evidence.

Output (annotate mode): DIR/<id>.gff3 and DIR/<id>.models.tsv.
Evaluation mode: DATASET/annotate_<split>/models.tsv, eval_genes.tsv (with the ORF checks
of each labelled CDS: label_pseudo, label_flags; model_segments), queries.tsv (editing setting per query);
per-group summary (gene and part level, label vs model pseudo calls) to stdout.

=#

export Annotator, load_annotator, annotate_record, write_gff, models_rows

using .PlastomeIR: find_ir_status
using .PlastomeRef
using .PlastomeAnnot
using .PlastomeFeatures
using .PlastomeTrees
using .PlastomeProtein
using Printf

acckey(a) = String(first(split(a, '.')))

function read_tsv(f, path)
    open(path) do io
        hdr = split(readline(io), '\t')
        col = Dict(String(h) => i for (i, h) in enumerate(hdr))
        for ln in eachline(io)
            f(col, split(ln, '\t'))
        end
    end
end

function each_fasta(f, path)
    io = open_maybe_gz(path)
    name = ""
    buf = IOBuffer()
    try
        for line in eachline(io)
            if startswith(line, '>')
                isempty(name) || f(name, String(take!(buf)))
                name = String(first(split(line[2:end])))
            else
                print(buf, uppercase(strip(line)))
            end
        end
        isempty(name) || f(name, String(take!(buf)))
    finally
        close(io)
    end
end


struct Part
    c::Candidate
    p::Float64
end

struct GeneModel
    gene::String
    kind::String
    strand::Char
    parts::Vector{Part}          # transcript order
    p::Float64                   # mean part probability
    pseudo::Bool
    flags::Vector{String}
end

"Position of the 5' end along the strand (for sorting in transcript order)."
fivep(c::Candidate, L) = c.strand == '+' ? c.first : L - c.last + 1
threep(c::Candidate, L) = c.strand == '+' ? c.last : L - c.first + 1

function spliced(q, qrc, parts::Vector{Part})
    io = IOBuffer()
    for pt in parts
        print(io, strand_seq(q, qrc, pt.c.first, pt.c.last, pt.c.strand))
    end
    return String(take!(io))
end

const STOP3 = ("TAA", "TAG", "TGA")
const ALWAYS_CDS = Set(["ycf15", "ycf68"])   # never called pseudo (ORF failures kept as flags)
const STOP_CU = ("TAA", "TAG", "TGA", "CAA", "CAG", "CGA")   # C-to-U editing can create a stop

"In-frame stops before the last codon of a spliced CDS."
internal_stops(s) = count(i -> s[i:i+2] in STOP3, 1:3:(3 * (length(s) ÷ 3 - 1) - 2))

stops_ok(n, ncod, editing) = n == 0 || (editing && (n <= 5 || n / max(1, ncod) <= 0.035))

"""
CDS checks of a spliced model -> (flags, pseudo). Start when the model holds part 1;
frame, terminal stop (TAA/TAG/TGA, plus CAA/CAG/CGA with editing) and internal stops
when it holds part 1 and the last part.
"""
function check_cds(s, hasfirst, haslast, editing)
    flags = String[]
    pseudo = false
    if hasfirst && length(s) >= 3 && !(s[1:3] in ("ATG", "GTG", "ACG"))
        push!(flags, "no_start")
    end
    if hasfirst && haslast
        if length(s) % 3 != 0
            push!(flags, "frame")
            pseudo = true
        end
        if length(s) >= 3 && !(s[end-2:end] in (editing ? STOP_CU : STOP3))
            push!(flags, "no_stop")
            pseudo = true
        end
        nst = internal_stops(s)
        if nst > 0
            push!(flags, "internal_stops=$nst")
            stops_ok(nst, length(s) ÷ 3, editing) || (pseudo = true)
        end
    end
    return flags, pseudo
end

# ---------------------------------------------------------------- junction profiles
const BIDX = Dict(UInt8('A') => 1, UInt8('C') => 2, UInt8('G') => 3, UInt8('T') => 4)
const CBASE = Dict(UInt8('A') => UInt8('T'), UInt8('C') => UInt8('G'), UInt8('G') => UInt8('C'), UInt8('T') => UInt8('A'))
const DON = -5:6      # donor window around the last exon base: 6 exon + 6 intron bases
const ACC = -6:5      # acceptor window around the first exon base: 6 intron + 6 exon bases

"Base at offset t (5'->3' along the strand) from genomic position x."
sbase(s::Vector{UInt8}, x::Int, t::Int, strand::Char) =
    strand == '+' ? s[mod1(x + t, length(s))] : get(CBASE, s[mod1(x - t, length(s))], UInt8('N'))

struct JunctionPWM
    donor::Matrix{Float64}      # log-odds vs uniform, length(DON) x 4
    acceptor::Matrix{Float64}
    n::Int
end

"""
Gene- and junction-specific profiles from the shortlisted references: for each CDS exon
(part k) followed by an intron (k+1) and an exon (k+2) of the same copy, the donor window
around the end of exon k and the acceptor window around the start of exon k+2. Key
(gene, k); kept with >= `min_refs` junctions. `rp`: (sequence, parts) per reference.
"""
function junction_pwms(rp; min_refs::Int=3)
    cnt = Dict{Tuple{String,Int},Tuple{Matrix{Float64},Matrix{Float64},Base.RefValue{Int}}}()
    for (s, parts) in rp
        isintron = Set((p.gene, p.copy, p.part) for p in parts if p.kind == "intron")
        cds = Dict((p.gene, p.copy, p.part) => p for p in parts if p.kind == "CDS")
        for ((g, c, k), p) in cds
            ((g, c, k + 1) in isintron && haskey(cds, (g, c, k + 2))) || continue
            nx = cds[(g, c, k + 2)]
            nx.strand == p.strand || continue
            e = p.strand == '+' ? p.last : p.first
            b = nx.strand == '+' ? nx.first : nx.last
            D, A, n = get!(() -> (zeros(length(DON), 4), zeros(length(ACC), 4), Ref(0)), cnt, (g, k))
            for (i, t) in enumerate(DON)
                j = get(BIDX, sbase(s, e, t, p.strand), 0)
                j > 0 && (D[i, j] += 1)
            end
            for (i, t) in enumerate(ACC)
                j = get(BIDX, sbase(s, b, t, nx.strand), 0)
                j > 0 && (A[i, j] += 1)
            end
            n[] += 1
        end
    end
    lo(M) = [log((M[i, j] + 0.5) / (sum(M[i, :]) + 2.0) / 0.25) for i in axes(M, 1), j in 1:4]
    return Dict(k => JunctionPWM(lo(D), lo(A), n[]) for (k, (D, A, n)) in cnt if n[] >= min_refs)
end

"Profile score of a junction: donor = last exon base at strand position e, acceptor = first exon base at b."
function pwm_score(P::JunctionPWM, S::Vector{UInt8}, e::Int, b::Int)
    L = length(S)
    sc = 0.0
    for (i, t) in enumerate(DON)
        j = get(BIDX, S[mod1(e + t, L)], 0)
        j > 0 && (sc += P.donor[i, j])
    end
    for (i, t) in enumerate(ACC)
        j = get(BIDX, S[mod1(b + t, L)], 0)
        j > 0 && (sc += P.acceptor[i, j])
    end
    return sc
end

"""
Editing evidence from the references that cover the query best (top `top` by
coverage): their own complete annotated CDS (parts 1, 3, ... of one copy) with in-frame
internal stops, or ending with CAA/CAG/CGA (stops created by C-to-U editing). Returns
(any reference with >= `min_n` such CDS, max count of CDS with internal stops, max count
of C-to-U terminal stops).
"""
function editing_evidence(idx, rp_by_ri, refw; top::Int=3, min_n::Int=3)
    nint_max, ncu_max = 0, 0
    for ri in first(sort(collect(keys(refw)); by=r -> -refw[r]), top)
        haskey(rp_by_ri, ri) || continue
        s, parts = rp_by_ri[ri]
        rc = revcomp_bytes(s)
        byc = Dict{Tuple{String,Int},Vector{RefPart}}()
        for p in parts
            p.kind == "CDS" && push!(get!(Vector{RefPart}, byc, (p.gene, p.copy)), p)
        end
        nint, ncu = 0, 0
        for (_, ps) in byc
            sort!(ps; by=p -> p.part)
            [p.part for p in ps] == collect(1:2:(2length(ps) - 1)) || continue
            sp = join(strand_seq(s, rc, p.first, p.last, p.strand) for p in ps)
            (length(sp) >= 6 && length(sp) % 3 == 0) || continue
            internal_stops(sp) > 0 && (nint += 1)
            sp[end-2:end] in ("CAA", "CAG", "CGA") && (ncu += 1)
        end
        nint_max, ncu_max = max(nint_max, nint), max(ncu_max, ncu)
    end
    return (nint_max >= min_n || ncu_max >= min_n), nint_max, ncu_max
end

"""
Junction refinement of a complete multi-exon CDS model: exons in strand coordinates; for
each junction in transcript order every donor/acceptor move of up to `w` nt each is
tried; a move is valid when the spliced prefix has no more in-frame stops than allowed
(and, at the last junction, a length divisible by 3 and a terminal stop). Among valid
moves the one with the highest score of the gene's junction profile (junction_pwms) wins,
ties to the smaller move; without a profile, the smallest valid move. When the current
junction is itself valid it is replaced only if the profile score improves by >=
`min_gain` (natural log units). Then the 3' end of the last exon is moved to the first
in-frame stop within `w3` nt if it has none.
"""
function refine_junctions(ch::Vector{Part}, q, qrc, editing, pwms; w::Int=9, w3::Int=30, min_gain::Float64=3.0)
    L = length(q)
    st = ch[1].c.strand
    S = st == '+' ? q : qrc
    ex = [[fivep(pt.c, L), threep(pt.c, L)] for pt in ch]          # strand coordinates
    seg(a, b) = String([S[mod1(x, L)] for x in a:b])
    stopset = editing ? STOP_CU : STOP3
    n = length(ex)
    moves = sort!([(d, a) for d in -w:w for a in -w:w]; by=x -> (abs(x[1]) + abs(x[2]), x[1], x[2]))
    for k in 1:(n - 1)
        P = ch[k + 1].c.part == ch[k].c.part + 2 ? get(pwms, (ch[k].c.gene, ch[k].c.part), nothing) : nothing
        prefix = join(seg(e[1], e[2]) for e in ex[1:(k - 1)])
        best, bsc, sc0 = nothing, -Inf, nothing
        for (d, a) in moves
            e1 = (ex[k][1], ex[k][2] + d)
            e2 = (ex[k + 1][1] + a, ex[k + 1][2])
            (e1[2] >= e1[1] + 2 && e2[2] >= e2[1] + 2 && e2[1] > e1[2]) || continue
            sp = prefix * seg(e1...) * seg(e2...)
            if k + 1 == n
                (length(sp) % 3 == 0 && sp[end-2:end] in stopset) || continue
                stops_ok(internal_stops(sp), length(sp) ÷ 3, editing) || continue
            else
                nst = count(i -> sp[i:i+2] in STOP3, 1:3:(3 * (length(sp) ÷ 3) - 2))
                stops_ok(nst, length(sp) ÷ 3, editing) || continue
            end
            sc = P === nothing ? 0.0 : pwm_score(P, S, e1[2], e2[1])
            (d, a) == (0, 0) && (sc0 = sc)
            if best === nothing || sc > bsc
                best, bsc = (d, a), sc
            end
            P === nothing && break                  # smallest valid move
        end
        best === nothing && continue
        sc0 !== nothing && bsc - sc0 < min_gain && continue
        ex[k][2] += best[1]
        ex[k + 1][1] += best[2]
    end
    # 3' end: without editing, a model running past an in-frame stop in its last exon is
    # cut there when <= 25% of its length is removed
    sp = join(seg(e[1], e[2]) for e in ex)
    if !editing && length(sp) % 3 == 0
        before = length(sp) - (ex[n][2] - ex[n][1] + 1)
        cs = 1:3:(length(sp) - 5)                          # codon starts before the last codon
        ifs = findfirst(y -> sp[y:y+2] in STOP3, cs)
        if ifs !== nothing
            pos = cs[ifs]
            if pos > before && length(sp) - (pos + 2) <= 0.25 * length(sp)
                ex[n][2] = ex[n][1] + (pos + 2 - before) - 1
                sp = join(seg(e[1], e[2]) for e in ex)
            end
        end
    end
    # 3' end: first in-frame stop near the current end when the model has none
    if !(length(sp) % 3 == 0 && length(sp) >= 3 && sp[end-2:end] in stopset)
        before = length(sp) - (ex[n][2] - ex[n][1] + 1)    # spliced bases before the last exon
        for newend in sort(collect((ex[n][2] - w3):(ex[n][2] + w3)); by=x -> abs(x - ex[n][2]))
            tot = before + (newend - ex[n][1] + 1)
            (tot % 3 == 0 && newend - ex[n][1] + 1 >= 3) || continue
            seg(newend - 2, newend) in stopset || continue
            cand = join(seg(e[1], e[2]) for e in ex[1:(n - 1)]) * seg(ex[n][1], newend)
            if stops_ok(internal_stops(cand), length(cand) ÷ 3, editing)
                ex[n][2] = newend
                break
            end
        end
    end
    return to_parts(ch, ex, st, L)
end

function to_parts(ch, ex, st, L)
    out = Part[]
    for (pt, e) in zip(ch, ex)
        f, l = st == '+' ? (mod1(e[1], L), mod1(e[2], L)) : (mod1(L - e[2] + 1, L), mod1(L - e[1] + 1, L))
        push!(out, Part(Candidate(pt.c, f, l), pt.p))
    end
    return out
end

"""
ORF search for a complete single-exon CDS model failing the checks: every start codon
(ATG, GTG, ACG) within `w5` nt (60) of the 5' end, in any frame; for each, the first in-frame
TAA/TAG/TGA; the ORF is kept when its length is within `maxdev` (0.25) of the model length; the
ORF closest to both model ends wins.
"""
function refine_orf(ch::Vector{Part}, q, qrc; w5::Int=60, maxdev::Float64=0.25)
    L = length(q)
    st = ch[1].c.strand
    S = st == '+' ? q : qrc
    a0, b0 = fivep(ch[1].c, L), threep(ch[1].c, L)
    len0 = b0 - a0 + 1
    cod(x) = String([S[mod1(x, L)], S[mod1(x + 1, L)], S[mod1(x + 2, L)]])
    best, bd = nothing, typemax(Int)
    for a in (a0 - w5):(a0 + w5)
        cod(a) in ("ATG", "GTG", "ACG") || continue
        for t in 0:3:floor(Int, (1 + maxdev) * len0)
            if cod(a + t) in STOP3
                e = a + t + 2
                len = e - a + 1
                if t > 0 && abs(len - len0) <= maxdev * len0
                    dd = abs(a - a0) + abs(e - b0)
                    dd < bd && ((best, bd) = ((a, e), dd))
                end
                break
            end
        end
    end
    best === nothing && return nothing
    return to_parts(ch, [[best[1], best[2]]], st, L)
end

"Intron length between consecutive exons a -> b of one copy (strand coordinates), -1 if none."
function gapbetween(a, b)
    (a === nothing || b === nothing || a.strand != b.strand) && return -1
    g = a.strand == '+' ? b.first - a.last - 1 : a.first - b.last - 1
    return 0 < g <= 4000 ? g : -1
end

"""
Short CDS exons (<= `maxlen` nt) of the shortlisted references: (gene, part) ->
[(strand-oriented sequence, intron length to the previous exon, to the next exon)],
-1 when that exon is absent or not cis-spliced.
"""
function short_exon_db(rp; maxlen::Int=80)
    db = Dict{Tuple{String,Int},Vector{Tuple{String,Int,Int}}}()
    for (s, parts) in rp
        cds = Dict((p.gene, p.copy, p.part) => p for p in parts if p.kind == "CDS")
        for ((g, c, k), p) in cds
            n = p.last - p.first + 1
            n <= maxlen || continue
            x = String([sbase(s, p.strand == '+' ? p.first : p.last, t, p.strand) for t in 0:(n - 1)])
            push!(get!(Vector{Tuple{String,Int,Int}}, db, (g, k)),
                  (x, gapbetween(get(cds, (g, c, k - 2), nothing), p), gapbetween(p, get(cds, (g, c, k + 2), nothing))))
        end
    end
    return db
end

"""
Hits of the reference variants of a missing short exon next to an exon of the chain:
downstream of the previous exon (`prev`, strand 3' end) at the reference intron lengths
+- `slack`, or upstream of the next exon (`next`, strand 5' start). Exons < 12 nt must
match exactly, longer ones with <= `maxmm` mismatches. Returns (start, end, mismatches,
|intron - median reference intron|) in strand coordinates, best first.
"""
function exon_hits(S::Vector{UInt8}, variants, prev, next; maxmm::Int=1, slack::Int=50)
    L = length(S)
    hits = Tuple{Int,Int,Int,Int}[]
    for (seq, gp, gn) in unique(variants)
        n = length(seq)
        mm_max = n < 12 ? 0 : maxmm
        if prev !== nothing
            gaps = [v[2] for v in variants if v[2] > 0]
            isempty(gaps) && continue
            med = sort(gaps)[cld(length(gaps), 2)]
            xs = (prev + 1 + minimum(gaps) - slack):(prev + 1 + maximum(gaps) + slack)
        else
            gaps = [v[3] for v in variants if v[3] > 0]
            isempty(gaps) && continue
            med = sort(gaps)[cld(length(gaps), 2)]
            xs = (next - maximum(gaps) - slack - n):(next - minimum(gaps) + slack - n)
        end
        for x in xs
            x > (prev === nothing ? typemin(Int) : prev) || continue
            (next === nothing || x + n - 1 < next) || continue
            mm = 0
            for t in 0:(n - 1)
                S[mod1(x + t, L)] == UInt8(seq[t + 1]) || (mm += 1)
                mm > mm_max && break
            end
            mm <= mm_max || continue
            intr = prev !== nothing ? x - prev - 1 : next - (x + n - 1) - 1
            push!(hits, (x, x + n - 1, mm, abs(intr - med)))
        end
    end
    return sort!(unique(hits); by=h -> (h[3], h[4]))
end

"""
Gene-model assembly from accepted parts:
  1. chains of one gene and strand along the transcript (part numbers increasing, gap
     -20..4000 nt);
  2. missing short exons (<= 80 nt, e.g. the 6-9 nt first exons of petB, petD, rpl16)
     searched next to the neighbouring exon with the references' sequences and intron
     lengths (exon_hits); a hit is kept when the extended chain passes the CDS checks
     (or, for an incomplete chain, starts with a start codon when it gains part 1);
  3. trans-splicing: a chain lacking part 1 whose first part k is common in the
     references (>= 0.5 x part 1) and not preceded by an intron there (intron part
     < 0.5 x part k) is joined to the best
     chain holding parts 1..k-2 of the gene (any strand, e.g. rps12);
  4. CDS checks; complete cis multi-exon models: refine_junctions; complete single-exon
     models failing: refine_orf.
"""
function build_models(acc::Vector{Part}, q, qrc, templates, editing::Bool, pwms, exdb)
    L = length(q)
    out = GeneModel[]
    tn(pth) = get(templates, pth, (0, 0.0))[1]
    complete_first(ch) = minimum(pt.c.part for pt in ch) == 1
    function complete_last(ch, gene, kind)
        m = maximum(pt.c.part for pt in ch)
        tn("$gene/$kind/$(m + 2)") < 0.5 * max(1, tn("$gene/$kind/$m"))
    end
    groups = Dict{Tuple{String,Char},Vector{Part}}()
    for pt in acc
        push!(get!(Vector{Part}, groups, (pt.c.gene, pt.c.strand)), pt)
    end
    bygene = Dict{String,Vector{Vector{Part}}}()
    for ((gene, st), ps) in groups
        sort!(ps; by=pt -> fivep(pt.c, L))
        chains = Vector{Vector{Part}}()
        for pt in ps
            if !isempty(chains)
                prev = chains[end][end]
                gap = fivep(pt.c, L) - threep(prev.c, L) - 1
                if pt.c.part > prev.c.part && -20 <= gap <= 4000
                    push!(chains[end], pt)
                    continue
                end
            end
            push!(chains, [pt])
        end
        # missing short exons
        S = st == '+' ? q : qrc
        for (ci, ch) in enumerate(chains)
            ch[1].c.kind == "CDS" || continue
            recovered = false
            for _ in 1:4
                pn = [pt.c.part for pt in ch]
                want = Int[]
                minimum(pn) > 1 && push!(want, minimum(pn) - 2)
                append!(want, [k for k in minimum(pn):2:maximum(pn) if !(k in pn)])
                complete_last(ch, gene, "CDS") || push!(want, maximum(pn) + 2)
                added = false
                for k in want
                    haskey(exdb, (gene, k)) || continue
                    ip = findfirst(pt -> pt.c.part == k - 2, ch)
                    inx = findfirst(pt -> pt.c.part == k + 2, ch)
                    (ip === nothing && inx === nothing) && continue
                    hits = exon_hits(S, exdb[(gene, k)], ip === nothing ? nothing : threep(ch[ip].c, L),
                                     inx === nothing ? nothing : fivep(ch[inx].c, L))
                    ref = ch[something(ip, inx)]
                    for (x, e, _, _) in first(hits, 20)
                        f, l = st == '+' ? (mod1(x, L), mod1(e, L)) : (mod1(L - e + 1, L), mod1(L - x + 1, L))
                        c = ref.c
                        np = Part(Candidate(c.gene, c.kind, k, c.strand, f, l, c.depth, c.nrefs, NaN, NaN, 0, 0, c.refs), ref.p)
                        ch2 = sort!(vcat(ch, [np]); by=pt -> pt.c.part)
                        sp = spliced(q, qrc, ch2)
                        hf, hl = complete_first(ch2), complete_last(ch2, gene, "CDS")
                        ok = if hf && hl
                            !check_cds(sp, true, true, editing)[2]
                        else
                            !hf || (length(sp) >= 3 && sp[1:3] in ("ATG", "GTG", "ACG"))
                        end
                        if ok
                            ch, added, recovered = ch2, true, true
                            break
                        end
                    end
                    added && break
                end
                added || break
            end
            chains[ci] = recovered ? ch : chains[ci]
            recovered && push!(ch, Part(ch[1].c, -1.0))      # marker, removed below
        end
        for ch in chains
            append!(get!(Vector{Vector{Part}}, bygene, gene), [ch])
        end
    end
    # trans-splicing joins
    joinedby = Dict{String,Vector{Vector{Part}}}()
    for (gene, chs) in bygene
        isrec(ch) = ch[end].p == -1.0
        used = Set{Int}()
        joined = Vector{Vector{Part}}()
        for (i, x) in enumerate(chs)
            x[1].c.kind == "CDS" || continue
            body = isrec(x) ? x[1:end-1] : x
            k = minimum(pt.c.part for pt in body)
            # trans junction in the references: part k common (>= 0.5 x part 1; excludes rare
            # splits, e.g. features crossing the sequence origin) and not preceded by an intron
            (k > 1 && tn("$gene/CDS/$k") >= 0.5 * max(1, tn("$gene/CDS/1")) &&
             tn("$gene/intron/$(k - 1)") < 0.5 * max(1, tn("$gene/CDS/$k"))) || continue
            cands = [j for (j, y) in enumerate(chs) if j != i && y[1].c.kind == "CDS" &&
                     minimum(pt.c.part for pt in (isrec(y) ? y[1:end-1] : y)) == 1 &&
                     maximum(pt.c.part for pt in (isrec(y) ? y[1:end-1] : y)) == k - 2]
            isempty(cands) && continue
            j = cands[argmax([sum(pt.p for pt in (isrec(chs[j]) ? chs[j][1:end-1] : chs[j])) for j in cands])]
            y = isrec(chs[j]) ? chs[j][1:end-1] : chs[j]
            push!(joined, vcat(y, body, [Part(body[1].c, -2.0)]))   # marker: trans-spliced
            push!(used, i, j)
        end
        joinedby[gene] = vcat([chs[i] for i in eachindex(chs) if !(i in used)], joined)
    end
    for (gene, chs) in joinedby, ch0 in chs
        flags = String[]
        ch = ch0
        if ch[end].p == -2.0
            ch = ch[1:end-1]
            push!(flags, "trans_spliced")
        end
        if ch[end].p == -1.0
            ch = ch[1:end-1]
            push!(flags, "exon_recovered")
        end
        trans = "trans_spliced" in flags
        kind = ch[1].c.kind
        # chains made only of rare part paths (< 0.2 x part 1 in the references; e.g. the
        # second piece of features split at the sequence origin) are artifacts
        all(pt -> tn("$gene/$kind/$(pt.c.part)") < 0.2 * tn("$gene/$kind/1"), ch) && continue
        st = ch[1].c.strand
        pseudo = false
        parts_n = [pt.c.part for pt in ch]
        n1 = tn("$gene/$kind/1")
        haslast = complete_last(ch, gene, kind)
        hasfirst = complete_first(ch)
        for k in minimum(parts_n):2:maximum(parts_n)
            (k in parts_n || tn("$gene/$kind/$k") < 0.5 * max(1, n1)) || push!(flags, "missing_part_$k")
        end
        if kind == "CDS"
            f2, pseudo = check_cds(spliced(q, qrc, ch), hasfirst, haslast, editing)
            # junction refinement: failing models, and passing ones unless their ends come
            # from the nearest reference
            if length(ch) >= 2 && hasfirst && haslast && !trans && (pseudo || !all(pt -> pt.c.near, ch))
                # junctions: failing models, and passing models whose profile score improves
                ch2 = refine_junctions(ch, q, qrc, editing, pwms)
                if [(p.c.first, p.c.last) for p in ch2] != [(p.c.first, p.c.last) for p in ch]
                    f3, p3 = check_cds(spliced(q, qrc, ch2), hasfirst, haslast, editing)
                    if !p3
                        ch, f2, pseudo = ch2, vcat(f3, ["junctions_refined"]), false
                    end
                end
            elseif pseudo && length(ch) == 1 && hasfirst && haslast
                ch2 = refine_orf(ch, q, qrc)
                if ch2 !== nothing
                    f3, p3 = check_cds(spliced(q, qrc, ch2), hasfirst, haslast, editing)
                    if !p3
                        ch, f2, pseudo = ch2, vcat(f3, ["orf_refined"]), false
                    end
                end
            end
            append!(flags, f2)
        end
        push!(out, GeneModel(gene, kind, st, ch, sum(pt.p for pt in ch) / length(ch), pseudo, flags))
    end
    return out
end

# ---------------------------------------------------------------- protein-level rescue
const STARTS3 = ("ATG", "GTG", "ACG")

"""
Gene model from a protein hit (PlastomeProtein): start codon (ATG/GTG/ACG) searched in
frame upstream of the aligned start (up to 10 codons beyond where the reference protein
would begin, not across a stop) and in the first 10 aligned codons, nearest to the
expected start; stop = first in-frame stop after the aligned end, within 15 codons of
the expected end (with editing also CAA/CAG/CGA within 2 codons of it).
"""
function hit_model(gene, h, sf, q, qrc, editing)
    L = length(q)
    fi = h.frame
    st = fi <= 3 ? '+' : '-'
    T = sf.prot[fi]
    S = st == '+' ? q : qrc
    star = UInt8('*')
    cod(j) = (x = frame_nt(fi, j); String([S[mod1(x + t, L)] for t in 0:2]))
    exp5 = h.tstart - (h.rstart - 1)
    js = h.tstart
    if h.rstart > 1 || !(cod(h.tstart) in STARTS3)
        cands = Int[]
        j = h.tstart
        while j >= max(1, exp5 - 10)
            (j < h.tstart && T[j] == star) && break
            cod(j) in STARTS3 && push!(cands, j)
            j -= 1
        end
        for j2 in (h.tstart + 1):min(h.tstart + 10, h.tend)
            cod(j2) in STARTS3 && push!(cands, j2)
        end
        isempty(cands) || (js = cands[argmin(abs.(cands .- exp5))])
    end
    expe = h.tend + (h.m - h.rend)
    je = nothing
    for j in (h.tend + 1):min(length(T), expe + 15)
        if T[j] == star || (editing && abs(j - (expe + 1)) <= 2 && cod(j) in ("CAA", "CAG", "CGA"))
            je = j
            break
        end
    end
    je === nothing && (je = min(length(T), expe + 1))
    a, e = frame_nt(fi, js), frame_nt(fi, je) + 2
    f, l = st == '+' ? (mod1(a, L), mod1(e, L)) : (mod1(L - e + 1, L), mod1(L - a + 1, L))
    pt = Part(Candidate(gene, "CDS", 1, st, f, l, 0.0, 0, NaN, NaN, 0, 0, Int[]), 0.0)
    flags, pseudo = check_cds(spliced(q, qrc, [pt]), true, true, editing)
    return GeneModel(gene, "CDS", st, [pt], 0.0, pseudo,
                     vcat(["protein_evidence", @sprintf("evalue=%.1e", h.evalue)], flags))
end

"""
Query window expected to hold `gene`, from its neighbours in the references (highest
coverage first, up to 5): the nearest features upstream and downstream of the gene
within `maxref` nt in the reference whose genes have exactly one model in the query;
the window is the shorter arc between those two models (<= `maxw` nt), padded by `pad`.
-> (start, length) or nothing
"""
function syntenic_window(gene, order, rp_by_ri, gm, L; maxref::Int=15000, maxw::Int=30000, pad::Int=300)
    cnt = Dict{String,Int}()
    span = Dict{String,Tuple{Int,Int}}()
    for m in gm
        m.pseudo && continue
        cnt[m.gene] = get(cnt, m.gene, 0) + 1
        span[m.gene] = (minimum(pt.c.first for pt in m.parts), maximum(pt.c.last for pt in m.parts))
    end
    for ri in first(order, 5)
        s, parts = rp_by_ri[ri]
        Lr = length(s)
        pos = Dict{String,Int}()
        for p in parts
            (p.kind != "intron" && p.copy == 1) || continue
            pos[p.gene] = min(get(pos, p.gene, typemax(Int)), p.first)
        end
        haskey(pos, gene) || continue
        c = pos[gene]
        left, dl, right, dr = "", -typemax(Int), "", typemax(Int)
        for (h, x) in pos
            (h != gene && get(cnt, h, 0) == 1) || continue
            d = mod(x - c + Lr ÷ 2, Lr) - Lr ÷ 2
            abs(d) <= maxref || continue
            if d < 0 && d > dl
                left, dl = h, d
            elseif d > 0 && d < dr
                right, dr = h, d
            end
        end
        (isempty(left) || isempty(right)) && continue
        sl, sr = span[left], span[right]
        w1, w2 = mod(sr[1] - sl[2], L), mod(sl[1] - sr[2], L)
        a, w = w1 <= w2 ? (sl[2], w1) : (sr[2], w2)
        w <= maxw || continue
        return mod1(a - pad, L), w + 2pad
    end
    return nothing
end

"""
Exon peptides of a reference CDS (parts of one copy, transcript order, one strand): for
each exon the amino acids of the codons lying wholly inside it, with the numbers of
leading and trailing bases of split codons; the terminal stop is dropped.
"""
function exon_peptides(s::Vector{UInt8}, ps)
    nt = UInt8[]
    bounds = Int[]
    for p in ps
        n = p.last - p.first + 1
        append!(nt, [sbase(s, p.strand == '+' ? p.first : p.last, t, p.strand) for t in 0:(n - 1)])
        push!(bounds, length(nt))
    end
    aa = translate_nt(nt, 1, length(nt) ÷ 3)
    (!isempty(aa) && aa[end] == UInt8('*')) && pop!(aa)
    out = Tuple{Vector{UInt8},Int,Int}[]
    o = 0
    for b in bounds
        clo, chi = cld(o, 3) + 1, min(fld(b, 3), length(aa))
        push!(out, (chi >= clo ? aa[clo:chi] : UInt8[], 3clo - 2 - (o + 1), b - 3chi))
        o = b
    end
    return out
end

"""
Protein-level rescue of CDS genes that are present in >= half of the shortlisted
references but have no gene model in the query.

Single-exon genes: the proteins of up to `nprot` references (highest coverage first) are
aligned in the syntenic window (syntenic_window) and in the whole six frames; the best
hit with E <= `emax` over its own search space (window or genome) and >= 50% of the
protein aligned gives the model (hit_model). The genome-wide search covers genes moved
out of the window by rearrangements (e.g. psbM in Selaginella).

Multi-exon cis-spliced genes (syntenic window required): the exon peptides of the
best-covering reference holding all exons are aligned in the window (E <= 1e-2, >= 40%
of the peptide); hits on the strand of the best one, in transcript order, aligned over
>= 50% of the protein, give the exons (ends extended to the peptide ends; start codon
within 10 codons of the first exon's start, stop within 15 codons after the last); a
complete model gets refine_junctions.

Models overlapping a model of another gene by more than half of their length are
dropped; models are flagged protein_evidence.
"""
const PROTLOG = Ref{Any}(nothing)          # IO for protein_rescue diagnostics (eval mode)
const PROTLOCK = ReentrantLock()

function plog(fields...)
    PROTLOG[] === nothing && return
    lock(PROTLOCK) do
        println(PROTLOG[], join(fields, '\t'))
    end
end

fmt_hit(h) = h === nothing ? "-\t-\t-" :
    @sprintf("%d\t%.1e\t%.2f", h.score, h.evalue, (h.rend - h.rstart + 1) / h.m)

function protein_rescue(gm, q, qrc, rp_by_ri, refw, single, editing, pwms;
                        nprot::Int=6, mincov::Float64=0.5, emax::Float64=1e-3, id::AbstractString="")
    L = length(q)
    have = Set(m.gene for m in gm if m.kind == "CDS")
    order = sort(collect(keys(rp_by_ri)); by=r -> -get(refw, r, 0.0))
    nref = length(order)
    prots = Dict{String,Vector{Vector{UInt8}}}()
    exons = Dict{String,Vector{Tuple{Vector{UInt8},Int,Int}}}()
    present = Dict{String,Int}()
    for ri in order
        s, parts = rp_by_ri[ri]
        byg = Dict{String,Vector{RefPart}}()
        for p in parts
            (p.kind == "CDS" && p.copy == 1) && push!(get!(Vector{RefPart}, byg, p.gene), p)
        end
        introns = Set((p.gene, p.part) for p in parts if p.kind == "intron" && p.copy == 1)
        for (g, ps) in byg
            present[g] = get(present, g, 0) + 1
            g in have && continue
            sort!(ps; by=p -> p.part)
            if g in single
                p = ps[1]
                p.part == 1 || continue
                v = get!(Vector{Vector{UInt8}}, prots, g)
                n = p.last - p.first + 1
                (length(v) < nprot && n >= 30) || continue
                nt = [sbase(s, p.strand == '+' ? p.first : p.last, t, p.strand) for t in 0:(n - 1)]
                aa = translate_nt(nt, 1, n ÷ 3)
                (!isempty(aa) && aa[end] == UInt8('*')) && pop!(aa)
                push!(v, aa)
            elseif !haskey(exons, g)
                # all exons 1, 3, ..., one strand, cis introns between them
                pn = [p.part for p in ps]
                (pn == collect(1:2:(2length(ps) - 1)) && length(ps) >= 2 &&
                 all(p -> p.strand == ps[1].strand, ps) &&
                 all(k -> (g, k) in introns, 2:2:(2length(ps) - 2))) || continue
                exons[g] = exon_peptides(s, ps)
            end
        end
    end
    todo = sort([g for g in union(keys(prots), keys(exons)) if get(present, g, 0) >= 0.5 * nref])
    isempty(todo) && return GeneModel[]
    sf = six_frames(q, qrc)
    occupied = Tuple{Int,Int,String}[]
    for m in gm
        m.pseudo && continue
        for pt in m.parts
            push!(occupied, (pt.c.first, pt.c.last, m.gene))
        end
    end
    clash(c, gene) = any(o -> o[3] != gene && overlap(c.first, c.last, o[1], o[2], L) > 0.5 * ilen(c.first, c.last, L), occupied)
    out = GeneModel[]
    for gene in todo
        win = syntenic_window(gene, order, rp_by_ri, gm, L)
        if haskey(prots, gene)
            best = nothing
            rawg, raww = nothing, nothing          # best hits regardless of thresholds (log)
            for P in prots[gene]
                hg = protein_search(sf, P)
                hw = win === nothing ? nothing : protein_search_window(sf, P, win[1], win[2])
                hg !== nothing && (rawg === nothing || hg.score > rawg.score) && (rawg = hg)
                hw !== nothing && (raww === nothing || hw.score > raww.score) && (raww = hw)
                # window hit and whole-genome hit are both admissible (the window can be wrong in
                # rearranged plastomes); E <= emax over the respective search space
                for h in (hw, hg)
                    h === nothing && continue
                    (h.evalue <= emax && h.rend - h.rstart + 1 >= mincov * h.m) || continue
                    (best === nothing || h.score > best.score) && (best = h)
                end
            end
            plog(id, gene, "single", length(prots[gene]), maximum(length, prots[gene]),
                 win === nothing ? "none" : "$(win[1])+$(win[2])", fmt_hit(raww), fmt_hit(rawg), best !== nothing)
            best === nothing && continue
            m = hit_model(gene, best, sf, q, qrc, editing)
            any(pt -> clash(pt.c, gene), m.parts) && continue
        elseif haskey(exons, gene) && win !== nothing
            m = exon_model(gene, exons[gene], sf, win, q, qrc, editing, pwms)
            plog(id, gene, "multi", length(exons[gene]), sum(length(p[1]) for p in exons[gene]),
                 "$(win[1])+$(win[2])", "-\t-\t-", "-\t-\t-", m !== nothing)
            (m === nothing || any(pt -> clash(pt.c, gene), m.parts)) && continue
        else
            haskey(exons, gene) && plog(id, gene, "multi", length(exons[gene]), 0, "none",
                                        "-\t-\t-", "-\t-\t-", false)
            continue
        end
        push!(out, m)
        for pt in m.parts
            push!(occupied, (pt.c.first, pt.c.last, gene))
        end
    end
    return out
end

"Multi-exon model from exon-peptide hits in a window (see protein_rescue)."
function exon_model(gene, peps, sf, win, q, qrc, editing, pwms)
    L = length(q)
    n = length(peps)
    hits = Vector{Any}(nothing, n)
    for (k, (P, _, _)) in enumerate(peps)
        length(P) >= 8 || continue
        h = protein_search_window(sf, P, win[1], win[2])
        (h !== nothing && h.evalue <= 1e-2 && h.rend - h.rstart + 1 >= 0.4 * h.m) && (hits[k] = h)
    end
    ks = [k for k in 1:n if hits[k] !== nothing]
    isempty(ks) && return nothing
    kb = ks[argmax([hits[k].score for k in ks])]
    st = hits[kb].frame <= 3 ? '+' : '-'
    keep = Int[]
    lastx = typemin(Int)
    for k in ks
        h = hits[k]
        (h.frame <= 3) == (st == '+') || continue
        x = frame_nt(h.frame, h.tstart)
        x > lastx || continue
        push!(keep, k)
        lastx = frame_nt(h.frame, h.tend)
    end
    sum(hits[k].rend - hits[k].rstart + 1 for k in keep) >= 0.5 * sum(length(p[1]) for p in peps) || return nothing
    S = st == '+' ? q : qrc
    cod(x) = String([S[mod1(x + t, L)] for t in 0:2])
    ex = Vector{Vector{Int}}()
    for k in keep
        h = hits[k]
        P, lead, trail = peps[k]
        a = frame_nt(h.frame, h.tstart) - 3(h.rstart - 1) - lead
        # the last exon's trailing bases are the stop codon, searched below
        e = frame_nt(h.frame, h.tend) + 2 + 3(length(P) - h.rend) + (k == n ? 0 : trail)
        push!(ex, [a, e])
    end
    for i in 2:length(ex)
        ex[i][1] > ex[i - 1][2] || return nothing
    end
    hasfirst, haslast = keep[1] == 1, keep[end] == n
    if hasfirst
        a = ex[1][1]
        cands = [a + 3t for t in -10:10 if cod(a + 3t) in STARTS3]
        isempty(cands) || (ex[1][1] = cands[argmin(abs.(cands .- a))])
    end
    if haslast
        e = ex[end][2]
        for t in 0:15
            c = cod(e + 1 + 3t)
            if c in STOP3 || (editing && t <= 2 && c in ("CAA", "CAG", "CGA"))
                ex[end][2] = e + 3 + 3t
                break
            end
        end
    end
    ch = Part[]
    for (k, x) in zip(keep, ex)
        f, l = st == '+' ? (mod1(x[1], L), mod1(x[2], L)) : (mod1(L - x[2] + 1, L), mod1(L - x[1] + 1, L))
        push!(ch, Part(Candidate(gene, "CDS", 2k - 1, st, f, l, 0.0, 0, NaN, NaN, 0, 0, Int[]), 0.0))
    end
    flags = ["protein_evidence", @sprintf("evalue=%.1e", minimum(hits[k].evalue for k in keep))]
    append!(flags, ["missing_part_$(2k - 1)" for k in 1:n if !(k in keep)])
    f2, pseudo = check_cds(spliced(q, qrc, ch), hasfirst, haslast, editing)
    if hasfirst && haslast && length(ch) >= 2 && length(keep) == n
        ch2 = refine_junctions(ch, q, qrc, editing, pwms)
        f3, p3 = check_cds(spliced(q, qrc, ch2), true, true, editing)
        if !p3
            ch, f2, pseudo = ch2, vcat(f3, ["junctions_refined"]), false
        end
    end
    return GeneModel(gene, "CDS", st, ch, 0.0, pseudo, vcat(flags, f2))
end

"Overlap of two candidates (1-based, circular)."
cov(a::Candidate, b::Candidate, L) = overlap(a.first, a.last, b.first, b.last, L)

function suppress(parts::Vector{Part}, L)
    sort!(parts; by=pt -> -pt.p)
    keep = Part[]
    for pt in parts
        c = pt.c
        clash = any(keep) do k
            d = k.c
            o = cov(c, d, L)
            o == 0 && return false
            same_path = d.gene == c.gene && d.kind == c.kind && d.part == c.part
            paralog = d.gene != c.gene && d.kind == c.kind && d.strand == c.strand &&
                      o >= 0.5 * min(ilen(c.first, c.last, L), ilen(d.first, d.last, L))
            same_path || paralog
        end
        clash || push!(keep, pt)
    end
    return keep
end

"""
Annotates one genome; returns (gene models, IR tuple or nothing, editing flag, accepted parts).
"""
function annotate_genome(id, seq, idx, models, single, threshold; exclude_genus="", exclude_accession="",
                         near_mash::Float64=0.01)
    q = Vector{UInt8}(seq)
    L = length(q)
    st, irr, _ = find_ir_status(seq)
    ir = (st == :found && irr !== nothing) ? (irr.irb[1], irr.irb[2], irr.ira[1], irr.ira[2]) : nothing
    sh = select_references(idx, sketch(seq; k=idx.k, s=idx.s); exclude_genus=exclude_genus,
                           exclude_accession=exclude_accession)
    mash = Dict(ri => d for (ri, d) in sh)
    refs = [(ri, Vector{UInt8}(read_sequence(idx, ri)),
             [(p.gene, p.kind, p.part, p.strand, p.first, p.last)
              for p in get(idx.parts, idx.records[ri].accession, RefPart[])]) for (ri, _) in sh]
    parts_by_path = Dict{String,Int}()
    for (_, _, ps) in refs, pth in unique("$(p[1])/$(p[2])/$(p[3])" for p in ps)
        parts_by_path[pth] = get(parts_by_path, pth, 0) + 1
    end
    # ends from the nearest reference when it is very close (Mash < near_mash)
    near_ref = (!isempty(sh) && sh[1][2] < near_mash) ? sh[1][1] : 0
    cands, refw, qrc = transfer_query(q, refs, idx.templates; min_depth=0.05, near_ref=near_ref)
    # editing tolerance from the references' own annotated CDS (editing_evidence)
    rp_by_ri = Dict(ri => (sq, get(idx.parts, idx.records[ri].accession, RefPart[])) for (ri, sq, _) in refs)
    editing, _, _ = editing_evidence(idx, rp_by_ri, refw)
    # single-exon CDS with the nearest reference's ends: kept when they give an ORF (length
    # divisible by 3, terminal stop, internal stops as allowed); otherwise the weighted-mode
    # ends and the start/stop rules, as for every other candidate
    function near_orf_ok(c)
        s = strand_seq(q, qrc, c.first, c.last, c.strand)
        length(s) >= 6 && length(s) % 3 == 0 && s[end-2:end] in (editing ? STOP_CU : STOP3) &&
            stops_ok(internal_stops(s), length(s) ÷ 3, editing)
    end
    cands = [!c.near ? refine_cds(c, q, qrc, single) :
             (c.kind == "CDS" && c.gene in single && !near_orf_ok(c)) ? refine_cds(mode_ends(c), q, qrc, single) : c
             for c in cands if c.kind != "intron"]
    pwms = junction_pwms(values(rp_by_ri))
    exdb = short_exon_db(values(rp_by_ri))
    nref_total = length(idx.records)
    acc = Part[]
    allc = Part[]                         # every candidate with its probability (diagnostics)
    for kind in ("CDS", "tRNA", "rRNA")
        ck = [c for c in cands if c.kind == kind]
        isempty(ck) && continue
        X = permutedims(reduce(hcat, [candidate_features(c, cands, q, qrc, L, refw, mash, idx.templates, nref_total,
                                                          parts_by_path, ir) for c in ck]))
        if haskey(models, kind)
            p = predict_gbdt(models[kind], X)
        else                                  # rule: depth >= 0.25 and rel_depth >= 0.5
            jd, jr = findfirst(==("depth"), FEATURES), findfirst(==("rel_depth"), FEATURES)
            p = [(X[i, jd] >= 0.25 && X[i, jr] >= 0.5) ? 1.0 : 0.0 for i in axes(X, 1)]
        end
        for (i, c) in enumerate(ck)
            push!(allc, Part(c, p[i]))
            p[i] >= threshold && push!(acc, Part(c, p[i]))
        end
    end
    acc = suppress(acc, L)
    gm = build_models(acc, q, qrc, idx.templates, editing, pwms, exdb)
    append!(gm, protein_rescue(gm, q, qrc, rp_by_ri, refw, single, editing, pwms; id=id))
    # identical models (same gene, kind and segments) reported once
    seen = Set{Any}()
    gm = [m for m in gm if !((m.gene, m.kind, [(pt.c.first, pt.c.last, pt.c.strand) for pt in m.parts]) in seen) &&
          (push!(seen, (m.gene, m.kind, [(pt.c.first, pt.c.last, pt.c.strand) for pt in m.parts])); true)]
    # genes reported as CDS whatever the ORF checks say (user decision): ycf15, ycf68
    gm = [(m.pseudo && m.gene in ALWAYS_CDS) ?
          GeneModel(m.gene, m.kind, m.strand, m.parts, m.p, false, vcat(m.flags, ["orf_checks_ignored"])) : m
          for m in gm]
    return gm, ir, editing, acc, qrc, allc
end

function write_gff(io, id, L, models; header::Bool=true)
    header && println(io, "##gff-version 3")
    println(io, "##sequence-region $id 1 $L")
    # circular sequence: a feature crossing the origin ends at last + L (GFF3 convention)
    println(io, join((id, "PlastomeAnnot", "region", 1, L, ".", "+", ".", "ID=$id;Is_circular=true"), '\t'))
    uend(c) = c.last < c.first ? c.last + L : c.last
    n = 0
    for m in sort(models; by=m -> minimum(pt.c.first for pt in m.parts))
        n += 1
        gid = "$(m.gene)_$n"
        lo = minimum(pt.c.first for pt in m.parts)
        hi = maximum(uend(pt.c) for pt in m.parts)
        attrs = "ID=$gid;Name=$(m.gene);gene=$(m.gene)" * (m.pseudo ? ";pseudo=true" : "") *
                (isempty(m.flags) ? "" : ";note=" * join(m.flags, ","))
        println(io, join((id, "PlastomeAnnot", m.pseudo ? "pseudogene" : "gene", lo, hi,
                          @sprintf("%.3f", m.p), m.strand, ".", attrs), '\t'))
        phase = 0
        for (k, pt) in enumerate(m.parts)
            c = pt.c
            ftype = m.kind == "CDS" ? "CDS" : m.kind
            ph = m.kind == "CDS" ? string(phase) : "."
            println(io, join((id, "PlastomeAnnot", ftype, c.first, uend(c), @sprintf("%.3f", pt.p), c.strand, ph,
                              "Parent=$gid;part=$(c.part);depth=$(@sprintf("%.3f", c.depth));nrefs=$(c.nrefs)"), '\t'))
            if m.kind == "CDS"
                phase = mod(3 - mod(ilen(c.first, c.last, L) - phase, 3), 3)
            end
        end
    end
end

function models_rows(id, models)
    rows = String[]
    for m in models
        segs = join(["$(pt.c.first)..$(pt.c.last):$(pt.c.strand)" for pt in m.parts], ";")
        push!(rows, join((id, m.gene, m.kind, m.strand, segs, join([pt.c.part for pt in m.parts], ","),
                          @sprintf("%.3f", m.p), m.pseudo ? 1 : 0, join(m.flags, ",")), '\t'))
    end
    return rows
end

function parse_args(argv)
    opt = Dict{String,Any}("dataset" => joinpath(homedir(), "Assembler/refseq_plastid/dataset"), "index" => "",
                           "models" => "", "threshold" => "0.5", "out" => "annotations", "eval" => "",
                           "max-queries" => "0", "fasta" => String[], "exclude-genus" => "",
                           "near-ref" => "0.01", "exclude" => "genus")
    i = 1
    while i <= length(argv)
        (startswith(argv[i], "--") && i < length(argv)) || error("bad argument $(argv[i])")
        k = argv[i][3:end]
        k == "fasta" ? push!(opt["fasta"], argv[i + 1]) : (opt[k] = argv[i + 1])
        i += 2
    end
    # defaults: production index and models when built, otherwise the train-split ones
    pick(a, b) = isdir(joinpath(opt["dataset"], a)) ? joinpath(opt["dataset"], a) : joinpath(opt["dataset"], b)
    isempty(opt["index"]) && (opt["index"] = pick("refindex_all", "refindex_train"))
    isempty(opt["models"]) && (opt["models"] = pick("models_final", "models"))
    return opt
end

function load_models(dir)
    ms = Dict{String,TreeModel}()
    for kind in ("CDS", "tRNA")
        f = joinpath(dir, "existence_$kind.model")
        isfile(f) ? (ms[kind] = load_model(f)) : @warn "no model $f: rule used for $kind"
    end
    return ms
end

function single_exon_genes(templates)
    single = Set{String}()
    for (pth, (n, _)) in templates
        g, kind, part = rsplit(pth, '/'; limit=3)
        (kind == "CDS" && part == "1") || continue
        get(templates, "$g/CDS/3", (0, 0.0))[1] < 0.2n && push!(single, String(g))
    end
    return single
end

# ------------------------------------------------------------------ evaluation

const LSEG_RE = r"(\d+)\.\.(\d+):([+-])"

"""
ORF checks of a labelled CDS (segments in transcript order), with the editing setting
used for the query -> (flags, pseudo). Shows whether the record's own CDS passes the
checks applied to the models.
"""
function label_seq(segs, q, qrc)
    io = IOBuffer()
    for m in eachmatch(LSEG_RE, segs)
        print(io, strand_seq(q, qrc, parse(Int, m[1]), parse(Int, m[2]), m[3][1]))
    end
    return String(take!(io))
end

label_orf(segs, q, qrc, editing) = check_cds(label_seq(segs, q, qrc), true, true, editing)

function evaluate(id, group, q, qrc, editing, models, acc, allc, labels, all_labels)
    L = length(q)
    truth = build_ref_parts(labels)
    allp = build_ref_parts(all_labels)
    # gene level: labelled features (unmasked) vs models
    grows = String[]
    gstats = Dict{String,Vector{Int}}()     # kind -> [n, found, exact, within3, pseudo_called]
    pq = Dict{String,Vector{Int}}()         # CDS -> [label ok/model ok, label ok/model pseudo, label pseudo/model ok, both pseudo]
    for (gene, kind, segs) in labels
        kind in ("CDS", "tRNA", "rRNA") || continue
        tp = [p for p in build_ref_parts([(gene, kind, segs)]) if p.kind != "intron"]
        isempty(tp) && continue
        tlen = sum(p.last - p.first + 1 for p in tp)
        best, bo = nothing, 0
        for m in models
            (m.gene == gene && m.kind == kind) || continue
            o = sum(overlap(pt.c.first, pt.c.last, p.first, p.last, L) for pt in m.parts for p in tp if pt.c.strand == p.strand; init=0)
            if o > bo
                best = m
                bo = o
            end
        end
        v = get!(() -> zeros(Int, 6), gstats, kind)
        v[1] += 1
        found = best !== nothing && bo >= 0.5 * tlen
        exact = w3 = same = false
        if found
            v[2] += 1
            pm = sort([(pt.c.first, pt.c.last) for pt in best.parts])
            pt_ = sort([(p.first, p.last) for p in tp])
            if length(pm) == length(pt_)
                exact = all(pm[k] == pt_[k] for k in eachindex(pm))
                w3 = all(abs(cdist(pm[k][1], pt_[k][1], L)) <= 3 && abs(cdist(pm[k][2], pt_[k][2], L)) <= 3 for k in eachindex(pm))
            end
            v[3] += exact
            v[4] += w3
            v[5] += best.pseudo
            # same spliced sequence: junctions placed differently in repeated bases, same product
            same = exact || spliced(q, qrc, best.parts) == label_seq(segs, q, qrc)
            v[6] += same
            if kind == "CDS"
                lp = label_orf(segs, q, qrc, editing)[2]
                w = get!(() -> zeros(Int, 4), pq, "CDS")
                w[1 + 2 * Int(lp) + Int(best.pseudo)] += 1
            end
        end
        lf, lp = kind == "CDS" ? label_orf(segs, q, qrc, editing) : (String[], false)
        # missed labels: best candidate of the gene overlapping the label (any probability)
        mp, md = "", ""
        if !found
            hits = [pt for pt in allc if pt.c.gene == gene && pt.c.kind == kind &&
                    any(p -> p.strand == pt.c.strand && overlap(pt.c.first, pt.c.last, p.first, p.last, L) > 0, tp)]
            if !isempty(hits)
                h = hits[argmax([pt.p for pt in hits])]
                mp, md = @sprintf("%.3f", h.p), @sprintf("%.3f", h.c.depth)
            else
                mp, md = "none", "none"
            end
        end
        push!(grows, join((id, group, gene, kind, segs, Int(found), Int(exact), Int(w3),
                           found ? Int(best.pseudo) : "", found ? join(best.flags, ",") : "",
                           kind == "CDS" ? Int(lp) : "", join(lf, ","),
                           found ? join(["$(pt.c.first)..$(pt.c.last):$(pt.c.strand)" for pt in best.parts], ";") : "",
                           found ? Int(same) : "", mp, md), '\t'))
    end
    # predicted models not overlapping any label of the gene
    nunl = count(m -> !m.pseudo && !any(p -> p.gene == m.gene && p.kind == m.kind &&
                     any(pt -> overlap(pt.c.first, pt.c.last, p.first, p.last, L) > 0, m.parts), allp), models)
    # part level (comparable with annotate_baseline.jl)
    pstats = Dict{String,Vector{Int}}()     # kind -> [n, found, exact5, exact3]
    for t in truth
        t.kind == "intron" && continue
        v = get!(() -> zeros(Int, 4), pstats, t.kind)
        v[1] += 1
        hits = [pt.c for pt in acc if pt.c.gene == t.gene && pt.c.kind == t.kind && pt.c.strand == t.strand &&
                overlap(pt.c.first, pt.c.last, t.first, t.last, L) >= 0.5 * (t.last - t.first + 1)]
        isempty(hits) && continue
        h = hits[argmax([overlap(c.first, c.last, t.first, t.last, L) for c in hits])]
        v[2] += 1
        v[3] += t.strand == '+' ? h.first == t.first : h.last == t.last
        v[4] += t.strand == '+' ? h.last == t.last : h.first == t.first
    end
    return grows, gstats, pstats, nunl, pq
end

function main(argv)
    opt = parse_args(argv)
    idx = load_index(opt["index"])
    models = load_models(opt["models"])
    single = single_exon_genes(idx.templates)
    thr = parse(Float64, opt["threshold"])
    nearm = parse(Float64, opt["near-ref"])
    @info "index $(opt["index"]): $(length(idx.records)) references; models: $(join(keys(models), ", "))"

    if isempty(opt["eval"])
        isempty(opt["fasta"]) && error("give --fasta FILE (or --eval SPLIT)")
        mkpath(opt["out"])
        for f in opt["fasta"]
            each_fasta(f) do id, seq
                ms, ir, editing, _, _, _ = annotate_genome(id, seq, idx, models, single, thr;
                                                           exclude_genus=opt["exclude-genus"], near_mash=nearm)
                # nearest reference (distance ~0: the genome, or a copy of it, is in the index)
                sh = select_references(idx, sketch(seq; k=idx.k, s=idx.s); exclude_genus=opt["exclude-genus"])
                nr = idx.records[sh[1][1]]
                stem = replace(id, r"[^A-Za-z0-9._-]" => "_")
                open(io -> write_gff(io, id, length(seq), ms), joinpath(opt["out"], "$stem.gff3"), "w")
                open(joinpath(opt["out"], "$stem.models.tsv"), "w") do io
                    println(io, "id\tgene\tkind\tstrand\tsegments\tparts\tp\tpseudo\tflags")
                    foreach(r -> println(io, r), models_rows(id, ms))
                end
                @info "$id: $(length(ms)) gene models ($(count(m -> m.pseudo, ms)) pseudo), IR $(ir === nothing ? "not found" : ir), editing-tolerant checks: $editing, nearest reference $(nr.accession) ($(nr.genus), $(nr.group)) Mash $(round(sh[1][2]; digits=4))$(sh[1][2] < nearm ? " (ends taken from it)" : "")"
            end
        end
        return
    end

    # evaluation mode
    ds = opt["dataset"]
    split_ = opt["eval"]
    maxq = parse(Int, opt["max-queries"])
    meta = Dict{String,NTuple{3,String}}()
    read_tsv(joinpath(ds, "records.tsv")) do c, f
        meta[acckey(f[c["accession"]])] = (f[c["group"]], f[c["genus"]], f[c["split"]])
    end
    used = Set{String}()
    read_tsv(joinpath(ds, "records_final.tsv")) do c, f
        f[c["use"]] == "1" && push!(used, acckey(f[c["accession"]]))
    end
    isq(a) = a in used && meta[a][3] == split_
    labels = Dict{String,Vector{Tuple{String,String,String}}}()
    all_labels = Dict{String,Vector{Tuple{String,String,String}}}()
    read_tsv(joinpath(ds, "labels.tsv")) do c, f
        a = acckey(f[c["accession"]])
        (isq(a) && !isempty(f[c["gene"]])) || return
        x = (String(f[c["gene"]]), String(f[c["feature"]]), String(f[c["segments"]]))
        push!(get!(() -> Tuple{String,String,String}[], all_labels, a), x)
        f[c["mask"]] == "0" && push!(get!(() -> Tuple{String,String,String}[], labels, a), x)
    end
    out = joinpath(ds, "annotate_$split_")
    mkpath(out)
    mio = open(joinpath(out, "models.tsv"), "w")
    gio = open(joinpath(out, "eval_genes.tsv"), "w")
    println(mio, "id\tgene\tkind\tstrand\tsegments\tparts\tp\tpseudo\tflags")
    println(gio, "query\tgroup\tgene\tkind\tsegments\tfound\texact\twithin3\tpseudo_called\tflags\tlabel_pseudo\tlabel_flags\tmodel_segments\tsame_seq\tmiss_best_p\tmiss_best_depth")
    qio = open(joinpath(out, "queries.tsv"), "w")
    PROTLOG[] = open(joinpath(out, "protein_log.tsv"), "w")
    println(PROTLOG[], "query\tgene\ttype\tn_prot\tprot_len\twindow\twin_score\twin_evalue\twin_cov\tgenome_score\tgenome_evalue\tgenome_cov\taccepted")
    println(qio, "query\tgroup\tgenus\tediting\tmodels\tpseudo_models")
    G = Dict{Tuple{String,String},Vector{Int}}()
    P = Dict{Tuple{String,String},Vector{Int}}()
    U = Dict{String,Vector{Int}}()               # group -> [queries, unlabelled models, models, editing]
    Q = Dict{Tuple{String,String},Vector{Int}}() # (group, flag) -> pseudo agreement counts
    batch = Tuple{String,String}[]
    nq = 0
    function run_batch!()
        res = Vector{Any}(undef, length(batch))
        Threads.@threads for i in eachindex(batch)
            local id, seq = batch[i]
            local a = acckey(id)
            local g, ge, _ = meta[a]
            local ms, ir, editing, acc, qrc, allc = annotate_genome(id, seq, idx, models, single, thr;
                                                                    exclude_genus=(opt["exclude"] == "genus" ? ge : ""),
                                                                    exclude_accession=id, near_mash=nearm)
            local ev = evaluate(id, g, Vector{UInt8}(seq), qrc, editing, ms, acc, allc,
                                get(labels, a, Tuple{String,String,String}[]),
                                get(all_labels, a, Tuple{String,String,String}[]))
            res[i] = (id, g, ge, ms, editing, ev)
        end
        for (id, g, ge, ms, editing, (grows, gstats, pstats, nunl, pq)) in res
            println(qio, join((id, g, ge, Int(editing), length(ms), count(m -> m.pseudo, ms)), '\t'))
            for (k, v) in pq
                Q[(g, k)] = get(Q, (g, k), zeros(Int, 4)) .+ v
            end
            foreach(r -> println(mio, r), models_rows(id, ms))
            foreach(r -> println(gio, r), grows)
            for (k, v) in gstats
                G[(g, k)] = get(G, (g, k), zeros(Int, 6)) .+ v
            end
            for (k, v) in pstats
                P[(g, k)] = get(P, (g, k), zeros(Int, 4)) .+ v
            end
            U[g] = get(U, g, zeros(Int, 4)) .+ [1, nunl, length(ms), Int(editing)]
        end
        flush(mio)
        flush(gio)
        flush(qio)
        nq_done = sum(v[1] for v in values(U))
        @info "annotated $nq_done queries"
        empty!(batch)
    end
    each_fasta(joinpath(ds, "genomes.fasta.gz")) do id, seq
        (maxq > 0 && nq >= maxq) && return
        isq(acckey(id)) || return
        nq += 1
        push!(batch, (id, seq))
        length(batch) >= 64 && run_batch!()
    end
    isempty(batch) || run_batch!()
    close(mio)
    close(gio)
    close(qio)
    close(PROTLOG[])
    PROTLOG[] = nothing
    groups = sort(collect(keys(U)); by=g -> -U[g][1])
    println("\ngene level (labelled features): found = model of the gene covering >= 50%; exact / within 3 nt = all segment ends; pseudo = found but called pseudo; same_seq = identical spliced sequence")
    println("group\tkind\tlabelled\tfound\texact/found\twithin3/found\tpseudo/found\tsame_seq/found")
    for g in groups, k in ("CDS", "tRNA", "rRNA")
        v = get(G, (g, k), zeros(Int, 6))
        v[1] == 0 && continue
        @printf("%s\t%s\t%d\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n", g, k, v[1], v[2] / v[1], v[3] / max(1, v[2]),
                v[4] / max(1, v[2]), v[5] / max(1, v[2]), v[6] / max(1, v[2]))
    end
    println("\npart level (as annotate_baseline.jl)")
    println("group\tkind\tlabelled\tfound\texact5/found\texact3/found")
    for g in groups, k in ("CDS", "tRNA", "rRNA")
        v = get(P, (g, k), zeros(Int, 4))
        v[1] == 0 && continue
        @printf("%s\t%s\t%d\t%.3f\t%.3f\t%.3f\n", g, k, v[1], v[2] / v[1], v[3] / max(1, v[2]), v[4] / max(1, v[2]))
    end
    println("\ngroup\tqueries\tmodels/query\tmodels without any label of the gene/query\tqueries with editing-tolerant checks")
    for g in groups
        u = U[g]
        @printf("%s\t%d\t%.1f\t%.2f\t%d\n", g, u[1], u[3] / u[1], u[2] / u[1], u[4])
    end
    println("\nfound CDS: ORF checks of the labelled CDS (same editing setting) vs the model's call")
    println("group\tlabel_ok/model_ok\tlabel_ok/model_pseudo\tlabel_pseudo/model_ok\tlabel_pseudo/model_pseudo")
    for g in groups
        v = get(Q, (g, "CDS"), zeros(Int, 4))
        @printf("%s\t%d\t%d\t%d\t%d\n", g, v...)
    end
    println("written $(out)/models.tsv, eval_genes.tsv, queries.tsv")
end

# ------------------------------------------------------------------ library API
"Reference index, existence models and settings, loaded once (load_annotator)."
struct Annotator
    idx::RefIndex
    models::Dict{String,TreeModel}
    single::Set{String}
    threshold::Float64
    near_mash::Float64
end

"""
    load_annotator(dataset=""; index="", models="", threshold=0.5, near_mash=0.01) -> Annotator

`dataset`: directory holding refindex_all/models_final (production) or refindex_train/models;
default ~/Assembler/refseq_plastid/dataset. `index`/`models` override the directories.
"""
function load_annotator(dataset::AbstractString=""; index::AbstractString="", models::AbstractString="",
                        threshold::Float64=0.5, near_mash::Float64=0.01)
    ds = isempty(dataset) ? joinpath(homedir(), "Assembler/refseq_plastid/dataset") : dataset
    pick(a, b) = isdir(joinpath(ds, a)) ? joinpath(ds, a) : joinpath(ds, b)
    ix = isempty(index) ? pick("refindex_all", "refindex_train") : index
    md = isempty(models) ? pick("models_final", "models") : models
    idx = load_index(ix)
    return Annotator(idx, load_models(md), single_exon_genes(idx.templates), threshold, near_mash)
end

"""
    annotate_record(a, id, seq; exclude_genus="") -> (models, ir, editing, nearest)

Gene models of one circular plastome sequence; `nearest` = (reference record, Mash
distance) or nothing.
"""
function annotate_record(a::Annotator, id::AbstractString, seq::AbstractString; exclude_genus::AbstractString="")
    s = uppercase(String(seq))
    ms, ir, editing, _, _, _ = annotate_genome(String(id), s, a.idx, a.models, a.single, a.threshold;
                                               exclude_genus=exclude_genus, near_mash=a.near_mash)
    sh = select_references(a.idx, sketch(s; k=a.idx.k, s=a.idx.s); exclude_genus=exclude_genus)
    return ms, ir, editing, (isempty(sh) ? nothing : (a.idx.records[sh[1][1]], sh[1][2]))
end

# ======================================================================
# BuildIndex (build_refindex.jl)
# ======================================================================
#=
build_refindex.jl - reference index for the plastome annotator (step 1).

  julia --threads auto PlastomeAnnotator.jl build-index [--dataset DIR] [--splits train] [--out DIR]
                                         [--k 21] [--sketch 1000] [--report validation]

References: records with use = 1 in records_final.tsv, split in --splits (comma list;
"all" = every split), structure not anomalous. Features: unmasked rows of labels.tsv
(CDS, tRNA, rRNA), split into parts (exons/introns) by PlastomeRef.build_ref_parts.

Output directory (default DATASET/refindex_<splits>): refs.tsv, refs.fasta,
sketches.bin, parts.tsv, templates.tsv.

--report SPLIT: for every used record of that split (e.g. validation), select references
as the annotator will (query genus excluded) and write neighbours_<SPLIT>.tsv with the
nearest reference, its distance and group, and how many of the shortlist are of the
query's group; a per-group summary goes to stdout.
=#
module BuildIndex
using ..PlastomeRef
using Printf

acckey(a) = String(first(split(a, '.')))

function read_tsv(f, path)
    open(path) do io
        hdr = split(readline(io), '\t')
        col = Dict(String(h) => i for (i, h) in enumerate(hdr))
        for ln in eachline(io)
            f(col, split(ln, '\t'))
        end
    end
end

"Streams a (gzipped) FASTA, calling f(id, seq) for every record."
function foreach_fasta(f, path)
    io = parentmodule(@__MODULE__).open_maybe_gz(path)
    name = ""
    buf = IOBuffer()
    try
        for line in eachline(io)
            if startswith(line, '>')
                isempty(name) || f(name, String(take!(buf)))
                name = String(first(split(line[2:end])))
            else
                print(buf, uppercase(strip(line)))
            end
        end
        isempty(name) || f(name, String(take!(buf)))
    finally
        close(io)
    end
end

function parse_args(argv)
    opt = Dict("dataset" => joinpath(homedir(), "Assembler/refseq_plastid/dataset"), "splits" => "train",
               "out" => "", "k" => "21", "sketch" => "1000", "report" => "")
    i = 1
    while i <= length(argv)
        startswith(argv[i], "--") && i < length(argv) || error("bad argument $(argv[i])")
        opt[argv[i][3:end]] = argv[i + 1]
        i += 2
    end
    isempty(opt["out"]) && (opt["out"] = joinpath(opt["dataset"], "refindex_" * replace(opt["splits"], ',' => '_')))
    return opt
end

function main(argv)
    opt = parse_args(argv)
    ds = opt["dataset"]
    k, s = parse(Int, opt["k"]), parse(Int, opt["sketch"])
    splits = opt["splits"] == "all" ? nothing : Set(split(opt["splits"], ','))

    meta = Dict{String,NTuple{3,String}}()          # key -> (group, genus, split)
    read_tsv(joinpath(ds, "records.tsv")) do c, f
        meta[acckey(f[c["accession"]])] = (f[c["group"]], f[c["genus"]], f[c["split"]])
    end
    use = Dict{String,Tuple{Bool,String}}()          # key -> (use, structure)
    read_tsv(joinpath(ds, "records_final.tsv")) do c, f
        st = haskey(c, "structure") ? f[c["structure"]] : ""
        use[acckey(f[c["accession"]])] = (f[c["use"]] == "1", st)
    end
    isref(a) = (u = get(use, a, (false, "")); u[1] && u[2] != "anomalous" &&
                (splits === nothing || meta[a][3] in splits))

    feats = Dict{String,Vector{Tuple{String,String,String}}}()
    read_tsv(joinpath(ds, "labels.tsv")) do c, f
        a = acckey(f[c["accession"]])
        (isref(a) && f[c["mask"]] == "0" && !isempty(f[c["gene"]])) || return
        push!(get!(() -> Tuple{String,String,String}[], feats, a), (f[c["gene"]], f[c["feature"]], f[c["segments"]]))
    end

    out = opt["out"]
    mkpath(out)
    records = RefRecord[]
    sketches = Vector{Vector{UInt64}}()
    parts = Dict{String,Vector{RefPart}}()
    batch = Tuple{String,String}[]
    fio = open(joinpath(out, "refs.fasta"), "w")
    function flush_batch!()
        sk = Vector{Vector{UInt64}}(undef, length(batch))
        Threads.@threads for i in eachindex(batch)
            sk[i] = sketch(batch[i][2]; k=k, s=s)
        end
        for (i, (id, seq)) in enumerate(batch)
            a = acckey(id)
            g, ge, sp = meta[a]
            print(fio, '>', id, '\n')
            off = position(fio)
            print(fio, seq, '\n')
            push!(records, RefRecord(id, g, ge, sp, use[a][2], length(seq), off))
            push!(sketches, sk[i])
            parts[id] = build_ref_parts(get(feats, a, Tuple{String,String,String}[]))
        end
        empty!(batch)
    end
    foreach_fasta(joinpath(ds, "genomes.fasta.gz")) do id, seq
        isref(acckey(id)) || return
        push!(batch, (id, seq))
        length(batch) >= 256 && flush_batch!()
    end
    flush_batch!()
    close(fio)
    write_index(out, records, sketches, parts; k=k, s=s)
    np = sum(length, values(parts); init=0)
    ntpl = countlines(joinpath(out, "templates.tsv")) - 1
    @info "index $(out): $(length(records)) references, $(np) parts, $(ntpl) part paths"
    bygroup = Dict{String,Int}()
    for r in records
        bygroup[r.group] = get(bygroup, r.group, 0) + 1
    end
    println("references per group: ", join(["$g $n" for (g, n) in sort(collect(bygroup); by=x -> -x[2])], ", "))

    isempty(opt["report"]) || report(opt, ds, meta, use, k, s, out)
end

function report(opt, ds, meta, use, k, s, out)
    idx = load_index(out)
    rsplit = opt["report"]
    rows = String[]
    stats = Dict{String,Vector{Tuple{Float64,Bool,Int}}}()
    lck = ReentrantLock()
    batch = Tuple{String,String}[]
    function run_batch!()
        res = Vector{String}(undef, length(batch))
        Threads.@threads for i in eachindex(batch)
            id, seq = batch[i]
            a = acckey(id)
            g, ge, _ = meta[a]
            sh = select_references(idx, sketch(seq; k=k, s=s); exclude_genus=ge, exclude_accession=id)
            r1 = idx.records[sh[1][1]]
            nsame = count(x -> idx.records[x[1]].group == g, sh)
            res[i] = join((id, g, ge, r1.accession, r1.group, r1.genus, @sprintf("%.4f", sh[1][2]), nsame, length(sh)), '\t')
            lock(lck) do
                push!(get!(() -> Tuple{Float64,Bool,Int}[], stats, g), (sh[1][2], r1.group == g, nsame))
            end
        end
        append!(rows, res)
        empty!(batch)
    end
    foreach_fasta(joinpath(ds, "genomes.fasta.gz")) do id, seq
        a = acckey(id)
        (haskey(use, a) && use[a][1] && meta[a][3] == rsplit) || return
        push!(batch, (id, seq))
        length(batch) >= 256 && run_batch!()
    end
    run_batch!()
    open(joinpath(out, "neighbours_$(rsplit).tsv"), "w") do io
        println(io, "query\tgroup\tgenus\tnearest\tnearest_group\tnearest_genus\tmash\tn_same_group\tn_shortlist")
        foreach(r -> println(io, r), rows)
    end
    println("group\tqueries\tnearest_same_group\tmedian_nearest_mash\tmedian_same_group_in_shortlist")
    for (g, v) in sort(collect(stats); by=x -> -length(x[2]))
        d = sort([x[1] for x in v])
        ns = sort([x[3] for x in v])
        @printf("%s\t%d\t%.3f\t%.4f\t%d\n", g, length(v), count(x -> x[2], v) / length(v), d[(end + 1) ÷ 2], ns[(end + 1) ÷ 2])
    end
end
end # module BuildIndex

# ======================================================================
# BuildTraining (build_training.jl)
# ======================================================================
#=
build_training.jl - step 3a: candidate features and labels for the existence model.

  julia --threads auto PlastomeAnnotator.jl build-training [--dataset DIR] [--index DIR] [--max-angio 3000]
                                         [--min-depth 0.05] [--max-queries N] [--out training]

Production models: --index DATASET/refindex_all --out training_all (references from all
splits; the query's genus is still excluded, so the features stay those of an unseen
genus).

Queries: used records of the train split (all non-angiosperms; angiosperms thinned to
--max-angio, every n-th by accession) and of the validation split (all). References:
the train index; the query's genus (and the query itself) are always excluded, so a
training query is annotated as a genome of an unseen genus.

Candidates come from PlastomeAnnot.transfer_query with --min-depth (lower than the
baseline 0.25: the model decides), then refine_cds as in the baseline.

Label (per candidate; introns are not labelled):
   1  positive   an unmasked label of the same gene, kind and strand overlaps >= 50% of
                 the shorter of the two
   0  negative   pseudo: a feature of the same gene annotated as pseudo (/pseudo in the
                 record) covers >= 50% of the shorter of the two (gene lost or
                 non-functional in the lineage, e.g. infA: Millen et al. 2001);
                 other_gene: overlaps an unmasked label of another gene that is not a
                 possible naming conflict; elsewhere: the query has this gene labelled
                 (unmasked) elsewhere and the candidate is not in an IR copy holding a
                 label of the gene
  -1  excluded   reason in label_reason: same_gene_partial (overlap < 50%, or a masked
                 label of the gene), masked_label (overlaps a masked label of another
                 gene), name_conflict (overlaps a label of a possibly mis-named gene: both
                 rRNA, or tRNAs with the same amino acid or anticodon, e.g. trnG-UCC
                 annotated as trnG-GCC), gene_unlabelled (gene absent from the query's
                 labels: lost, or missing from the annotation), ir_other_copy (gene
                 labelled in an IR copy, candidate in an IR copy)

Features (NaN = not applicable): see FEATURES below.

Output: DATASET/<--out>/candidates.tsv
=#
module BuildTraining
using ..PlastomeRef
using ..PlastomeAnnot
using ..PlastomeFeatures
using Printf

acckey(a) = String(first(split(a, '.')))

function read_tsv(f, path)
    open(path) do io
        hdr = split(readline(io), '\t')
        col = Dict(String(h) => i for (i, h) in enumerate(hdr))
        for ln in eachline(io)
            f(col, split(ln, '\t'))
        end
    end
end

function foreach_fasta(f, path)
    io = parentmodule(@__MODULE__).open_maybe_gz(path)
    name = ""
    buf = IOBuffer()
    try
        for line in eachline(io)
            if startswith(line, '>')
                isempty(name) || f(name, String(take!(buf)))
                name = String(first(split(line[2:end])))
            else
                print(buf, uppercase(strip(line)))
            end
        end
        isempty(name) || f(name, String(take!(buf)))
    finally
        close(io)
    end
end

function features_and_label(c, cands, q, qrc, L, refw, mash, templates, nref_total, parts_by_path,
                            truth, allp, maskedp, pseudop, ir)
    feats = candidate_features(c, cands, q, qrc, L, refw, mash, templates, nref_total, parts_by_path, ir)
    clen = ilen(c.first, c.last, L)
    reg = region_code(ir, c.first, L)
    # label
    sameg = [t for t in truth if t.gene == c.gene && t.kind == c.kind]
    for t in sameg
        t.strand == c.strand || continue
        ov = overlap(c.first, c.last, t.first, t.last, L)
        ov >= 0.5 * min(clen, t.last - t.first + 1) && return feats, 1, "match"
    end
    for p in pseudop
        (p.gene == c.gene && p.kind == c.kind &&
         overlap(c.first, c.last, p.first, p.last, L) >= 0.5 * min(clen, p.last - p.first + 1)) &&
            return feats, 0, "pseudo"
    end
    for p in allp
        (p.gene == c.gene && p.kind == c.kind && overlap(c.first, c.last, p.first, p.last, L) > 0) &&
            return feats, -1, "same_gene_partial"
    end
    any(p -> overlap(c.first, c.last, p.first, p.last, L) > 0, maskedp) && return feats, -1, "masked_label"
    others = [t for t in truth if t.gene != c.gene && overlap(c.first, c.last, t.first, t.last, L) > 0]
    any(t -> name_conflict(c, t), others) && return feats, -1, "name_conflict"
    isempty(others) || return feats, 0, "other_gene"
    isempty(sameg) && return feats, -1, "gene_unlabelled"
    if reg in (1.0, 3.0) && any(t -> region_code(ir, t.first, L) in (1.0, 3.0), sameg)
        return feats, -1, "ir_other_copy"
    end
    return feats, 0, "elsewhere"
end

function parse_args(argv)
    opt = Dict("dataset" => joinpath(homedir(), "Assembler/refseq_plastid/dataset"), "index" => "",
               "max-angio" => "3000", "min-depth" => "0.05", "max-queries" => "0", "out" => "training")
    i = 1
    while i <= length(argv)
        (startswith(argv[i], "--") && i < length(argv)) || error("bad argument $(argv[i])")
        opt[argv[i][3:end]] = argv[i + 1]
        i += 2
    end
    isempty(opt["index"]) && (opt["index"] = joinpath(opt["dataset"], "refindex_train"))
    return opt
end

function main(argv)
    opt = parse_args(argv)
    ds = opt["dataset"]
    idx = load_index(opt["index"])
    nref_total = length(idx.records)
    mind = parse(Float64, opt["min-depth"])
    maxq = parse(Int, opt["max-queries"])
    single = Set{String}()
    for (pth, (n, _)) in idx.templates
        g, kind, part = rsplit(pth, '/'; limit=3)
        (kind == "CDS" && part == "1") || continue
        get(idx.templates, "$g/CDS/3", (0, 0.0))[1] < 0.2n && push!(single, String(g))
    end

    meta = Dict{String,NTuple{3,String}}()
    read_tsv(joinpath(ds, "records.tsv")) do c, f
        meta[acckey(f[c["accession"]])] = (f[c["group"]], f[c["genus"]], f[c["split"]])
    end
    used = Set{String}()
    read_tsv(joinpath(ds, "records_final.tsv")) do c, f
        f[c["use"]] == "1" && push!(used, acckey(f[c["accession"]]))
    end
    # train queries: all non-angiosperms, angiosperms thinned
    maxa = parse(Int, opt["max-angio"])
    angio = sort([a for a in used if meta[a][3] == "train" && meta[a][1] == "angiosperms"])
    step = max(1, cld(length(angio), max(1, maxa)))
    keep_angio = Set(angio[1:step:end])
    isq(a) = a in used && (meta[a][3] == "validation" ||
                           (meta[a][3] == "train" && (meta[a][1] != "angiosperms" || a in keep_angio)))

    irs = Dict{String,NTuple{4,Int}}()
    read_tsv(joinpath(ds, "ir_table.tsv")) do c, f
        f[c["status"]] == "found" || return
        irs[acckey(f[c["accession"]])] = (parse(Int, f[c["irb_start"]]), parse(Int, f[c["irb_end"]]),
                                          parse(Int, f[c["ira_start"]]), parse(Int, f[c["ira_end"]]))
    end
    labels = Dict{String,Vector{Tuple{String,String,String}}}()
    all_labels = Dict{String,Vector{Tuple{String,String,String}}}()
    masked_labels = Dict{String,Vector{Tuple{String,String,String}}}()
    pseudo_labels = Dict{String,Vector{Tuple{String,String,String}}}()
    read_tsv(joinpath(ds, "labels.tsv")) do c, f
        a = acckey(f[c["accession"]])
        (isq(a) && !isempty(f[c["gene"]])) || return
        x = (String(f[c["gene"]]), String(f[c["feature"]]), String(f[c["segments"]]))
        push!(get!(() -> Tuple{String,String,String}[], all_labels, a), x)
        if f[c["mask"]] == "0"
            push!(get!(() -> Tuple{String,String,String}[], labels, a), x)
        else
            push!(get!(() -> Tuple{String,String,String}[], masked_labels, a), x)
            f[c["status"]] == "pseudo" && push!(get!(() -> Tuple{String,String,String}[], pseudo_labels, a), x)
        end
    end

    out = joinpath(ds, opt["out"])
    mkpath(out)
    oio = open(joinpath(out, "candidates.tsv"), "w")
    println(oio, join(vcat(["split", "query", "group", "gene", "kind", "part", "strand", "first", "last",
                            "label", "label_reason"], FEATURES), '\t'))
    batch = Tuple{String,String}[]
    nq = 0
    ndone = 0
    lck = ReentrantLock()
    function run_batch!()
        rows = Vector{Vector{String}}(undef, length(batch))
        Threads.@threads for i in eachindex(batch)
            id, seq = batch[i]
            a = acckey(id)
            g, ge, sp = meta[a]
            q = Vector{UInt8}(seq)
            L = length(q)
            sh = select_references(idx, sketch(seq; k=idx.k, s=idx.s); exclude_genus=ge, exclude_accession=id)
            mash = Dict(ri => d for (ri, d) in sh)
            refs = [(ri, Vector{UInt8}(read_sequence(idx, ri)),
                     [(p.gene, p.kind, p.part, p.strand, p.first, p.last)
                      for p in get(idx.parts, idx.records[ri].accession, RefPart[])]) for (ri, _) in sh]
            parts_by_path = Dict{String,Int}()
            for (_, _, ps) in refs, pth in unique("$(p[1])/$(p[2])/$(p[3])" for p in ps)
                parts_by_path[pth] = get(parts_by_path, pth, 0) + 1
            end
            cands, refw, qrc = transfer_query(q, refs, idx.templates; min_depth=mind)
            cands = [refine_cds(c, q, qrc, single) for c in cands]
            truth = [p for p in build_ref_parts(get(labels, a, Tuple{String,String,String}[])) if p.kind != "intron"]
            allp = build_ref_parts(get(all_labels, a, Tuple{String,String,String}[]))
            maskedp = [p for p in build_ref_parts(get(masked_labels, a, Tuple{String,String,String}[])) if p.kind != "intron"]
            pseudop = [p for p in build_ref_parts(get(pseudo_labels, a, Tuple{String,String,String}[])) if p.kind != "intron"]
            ir = get(irs, a, nothing)
            r = String[]
            for c in cands
                c.kind == "intron" && continue
                feats, lab, why = features_and_label(c, cands, q, qrc, L, refw, mash, idx.templates, nref_total,
                                                     parts_by_path, truth, allp, maskedp, pseudop, ir)
                push!(r, join(vcat([sp, id, g, c.gene, c.kind, string(c.part), string(c.strand), string(c.first),
                                    string(c.last), string(lab), why],
                                   [isnan(x) ? "NaN" : @sprintf("%.5g", x) for x in feats]), '\t'))
            end
            rows[i] = r
        end
        for r in rows
            foreach(x -> println(oio, x), r)
        end
        flush(oio)
        ndone += length(batch)
        @info "processed $(ndone) queries"
        empty!(batch)
    end
    foreach_fasta(joinpath(ds, "genomes.fasta.gz")) do id, seq
        (maxq > 0 && nq >= maxq) && return
        isq(acckey(id)) || return
        nq += 1
        push!(batch, (id, seq))
        length(batch) >= 64 && run_batch!()
    end
    isempty(batch) || run_batch!()
    close(oio)
    println("written $(out)/candidates.tsv ($(ndone) queries)")
end
end # module BuildTraining

# ======================================================================
# TrainExistence (train_existence.jl)
# ======================================================================
#=
train_existence.jl - step 3b: existence models (one per kind: CDS, tRNA, rRNA) from
training/candidates.tsv (build_training.jl), and their evaluation.

  julia --threads auto PlastomeAnnotator.jl train [--dataset DIR] [--threshold 0.5] [--n-trees 400]
                                          [--training training] [--out models] [--final 0|1]

--final 1 (production models): rows of every split are used; 10% of the queries (by
hash) are held out for early stopping and the summary below is computed on them (so it
is slightly optimistic: they also chose the number of trees).

Rows with label 0/1 are used. Train split -> fitting; validation queries are split in two
halves by a hash of the accession: one half for early stopping, the other ("test") for
the evaluation below. The baseline for comparison accepts every candidate with depth >=
0.25 (annotate_baseline.jl).

Evaluation on the test half, per kind and per group: AUC; at --threshold, recall of
positives and rejection of negatives, against the baseline rule. Excluded rows with
reason gene_unlabelled (the gene is absent from the record's labels) are scored too: a
record with several such candidates at p >= 0.9 is listed in
models/unlabelled_high_confidence.tsv (candidates for annotations missing genes).

Output: DATASET/models/existence_<kind>.model, importance_<kind>.tsv, test_predictions.tsv
=#
module TrainExistence
using ..PlastomeTrees
using Printf

const META = ["split", "query", "group", "gene", "kind", "part", "strand", "first", "last", "label", "label_reason"]

function parse_args(argv)
    opt = Dict("dataset" => joinpath(homedir(), "Assembler/refseq_plastid/dataset"), "threshold" => "0.5",
               "n-trees" => "400", "training" => "training", "out" => "models", "final" => "0")
    i = 1
    while i <= length(argv)
        (startswith(argv[i], "--") && i < length(argv)) || error("bad argument $(argv[i])")
        opt[argv[i][3:end]] = argv[i + 1]
        i += 2
    end
    return opt
end

"FNV-1a hash of a string (stable across runs)."
function fnv(s::AbstractString)
    h = 0xcbf29ce484222325
    for b in codeunits(s)
        h = (h ⊻ b) * 0x100000001b3
    end
    return h
end

function main(argv)
    opt = parse_args(argv)
    ds = opt["dataset"]
    thr = parse(Float64, opt["threshold"])
    path = joinpath(ds, opt["training"], "candidates.tsv")
    final = opt["final"] == "1"
    io = open(path)
    hdr = split(readline(io), '\t')
    col = Dict(String(h) => i for (i, h) in enumerate(hdr))
    fnames = String.(hdr[(length(META) + 1):end])
    rows = Dict{String,Vector{Vector{SubString{String}}}}()     # kind -> rows
    for ln in eachline(io)
        f = split(ln, '\t')
        push!(get!(() -> Vector{SubString{String}}[], rows, String(f[col["kind"]])), f)
    end
    close(io)
    out = joinpath(ds, opt["out"])
    mkpath(out)
    pio = open(joinpath(out, "test_predictions.tsv"), "w")
    println(pio, "query\tgroup\tgene\tkind\tpart\tstrand\tfirst\tlast\tlabel\tlabel_reason\tdepth\tp")
    uio = open(joinpath(out, "unlabelled_high_confidence.tsv"), "w")
    println(uio, "query\tgroup\tsplit\tgene\tkind\tpart\tstrand\tfirst\tlast\tdepth\tp")
    jd = col["depth"]
    nmeta = length(META)
    feat(f) = [x == "NaN" ? NaN : parse(Float64, x) for x in f[(nmeta + 1):end]]
    lab(f) = parse(Int, f[col["label"]])
    istest(f) = isodd(fnv(f[col["query"]]) % 2)
    mat(v) = isempty(v) ? zeros(0, length(fnames)) : permutedims(reduce(hcat, [feat(f) for f in v]))
    fr(v, c) = isempty(v) ? NaN : count(c, v) / length(v)
    for kind in ("CDS", "tRNA", "rRNA")
        haskey(rows, kind) || continue
        rs = rows[kind]
        if final
            # production: every split; 10% of the queries (hash) for early stopping, reported below
            tr = [f for f in rs if lab(f) >= 0 && fnv(f[col["query"]]) % 10 != 0]
            es = [f for f in rs if lab(f) >= 0 && fnv(f[col["query"]]) % 10 == 0]
            te = es
        else
            tr = [f for f in rs if f[col["split"]] == "train" && lab(f) >= 0]
            es = [f for f in rs if f[col["split"]] == "validation" && lab(f) >= 0 && !istest(f)]
            te = [f for f in rs if f[col["split"]] == "validation" && lab(f) >= 0 && istest(f)]
        end
        X, y = mat(tr), Float64[lab(f) for f in tr]
        Xe, ye = mat(es), Float64[lab(f) for f in es]
        Xt, yt = mat(te), Float64[lab(f) for f in te]
        @info "$kind: train $(length(y)) rows ($(Int(sum(y))) positive), early-stop $(length(ye)), test $(length(yt))"
        (isempty(y) || sum(y) == 0 || sum(y) == length(y)) && (@warn "$kind: one class only, skipped"; continue)
        m = fit_gbdt(X, y, fnames; Xv=Xe, yv=ye, n_trees=parse(Int, opt["n-trees"]))
        save_model(joinpath(out, "existence_$kind.model"), m)
        open(joinpath(out, "importance_$kind.tsv"), "w") do o
            tot = sum(m.importance)
            for j in sortperm(m.importance; rev=true)
                @printf(o, "%s\t%.4f\n", fnames[j], tot > 0 ? m.importance[j] / tot : 0.0)
            end
        end
        pt = predict_gbdt(m, Xt)
        for (k, f) in enumerate(te)
            println(pio, join((f[col["query"]], f[col["group"]], f[col["gene"]], kind, f[col["part"]], f[col["strand"]],
                               f[col["first"]], f[col["last"]], f[col["label"]], f[col["label_reason"]], f[jd],
                               @sprintf("%.4f", pt[k])), '\t'))
        end
        # summary per group
        @printf("\n%s: %d trees, test AUC %.4f (baseline depth AUC %.4f)\n", kind, length(m.trees), auc(yt, pt),
                auc(yt, [parse(Float64, f[jd]) for f in te]))
        println("group\tpos\tneg\trecall_model\treject_model\trecall_base\treject_base")
        for g in sort(unique(String(f[col["group"]]) for f in te))
            ix = [k for k in eachindex(te) if te[k][col["group"]] == g]
            pos = [k for k in ix if yt[k] == 1]
            neg = [k for k in ix if yt[k] == 0]
            dep = Dict(k => parse(Float64, te[k][jd]) for k in ix)
            @printf("%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\n", g, length(pos), length(neg),
                    fr(pos, k -> pt[k] >= thr), fr(neg, k -> pt[k] < thr),
                    fr(pos, k -> dep[k] >= 0.25), fr(neg, k -> dep[k] < 0.25))
        end
        println("top features: ", join(["$(fnames[j]) $(@sprintf("%.3f", m.importance[j] / sum(m.importance)))"
                                        for j in first(sortperm(m.importance; rev=true), 8)], ", "))
        # candidates of genes absent from the record's labels
        un = [f for f in rs if f[col["label_reason"]] == "gene_unlabelled"]
        if !isempty(un)
            pu = predict_gbdt(m, mat(un))
            for (k, f) in enumerate(un)
                pu[k] >= 0.9 || continue
                println(uio, join((f[col["query"]], f[col["group"]], f[col["split"]], f[col["gene"]], kind, f[col["part"]],
                                   f[col["strand"]], f[col["first"]], f[col["last"]], f[jd], @sprintf("%.4f", pu[k])), '\t'))
            end
            @printf("%s: gene_unlabelled candidates %d, with p >= 0.9: %d\n", kind, length(un), count(>=(0.9), pu))
        end
    end
    close(pio)
    close(uio)
    println("\nwritten $(out)/existence_*.model, importance_*.tsv, test_predictions.tsv, unlabelled_high_confidence.tsv")
end
end # module TrainExistence

# ======================================================================
# IRScan (ir_scan.jl)
# ======================================================================
#=
ir_scan.jl - IR detection for every genome of the land-plant plastid dataset, and
IRscope-style junction plots for chosen records. Needs PlastomeIR.jl in the same
directory. Base Julia; .gz input through open_maybe_gz (CodecZlib, else the gzip command).

  julia --threads auto PlastomeAnnotator.jl ir-scan [--dataset DIR] [--out ir_table.tsv] [--min-len 500]
  julia PlastomeAnnotator.jl ir-scan --plot NC_000932,NC_039750 --svg ir_junctions.svg [--dataset DIR]

Table columns: accession, length, status (found / not_detected / low_identity /
inconsistent / error), irb_start, irb_end, ira_start, ira_end, len_irb, len_ira, lsc,
ssc, gc_lsc, gc_irb, gc_ssc, gc_ira, identity (substitutions only), aligned_frac,
n_unknown, n_snp, n_indel,
indel_bp_irb, indel_bp_ira (coordinates 1-based, circular; IRb follows the LSC), then the
largest direct repeat (both copies on the same strand, >= 1 kb) for every genome
whatever its IR status: dr_len, dr_start1, dr_end1, dr_start2, dr_end2.
Coordinates are written for found and low_identity pairs; only `found` IRs should be
used downstream. Plots use features from labels.tsv (QC failures left out, training
masks ignored) and the organism name from records.tsv.
=#
module IRScan
using ..PlastomeIR
using Printf

const IR_HEADER = ["accession", "length", "status", "irb_start", "irb_end", "ira_start", "ira_end",
                   "len_irb", "len_ira", "lsc", "ssc", "gc_lsc", "gc_irb", "gc_ssc", "gc_ira",
                   "identity", "aligned_frac", "n_unknown", "n_snp", "n_indel", "indel_bp_irb", "indel_bp_ira"]
const DR_HEADER = ["dr_len", "dr_start1", "dr_end1", "dr_start2", "dr_end2"]
const HEADER = vcat(IR_HEADER, DR_HEADER)
const NEMPTY = length(IR_HEADER) - 3
const DR_EMPTY = join(fill("", length(DR_HEADER)), '\t')

function parse_args(argv)
    opt = Dict{String,String}("dataset" => joinpath(homedir(), "Assembler/refseq_plastid/dataset"),
                              "out" => "", "plot" => "", "svg" => "ir_junctions.svg",
                              "min-len" => "500")
    i = 1
    while i <= length(argv)
        k = argv[i]
        startswith(k, "--") || error("unexpected argument $k")
        i + 1 <= length(argv) || error("missing value for $k")
        opt[k[3:end]] = argv[i + 1]
        i += 2
    end
    isempty(opt["out"]) && (opt["out"] = joinpath(opt["dataset"], "ir_table.tsv"))
    return opt
end

acckey(a) = String(first(split(a, '.')))
emptyrow(id, L, status) = string(id, '\t', L, '\t', status, "\t"^NEMPTY)

function dr_part(seq)
    dr = find_direct_repeat(seq)
    return dr === nothing ? DR_EMPTY : join((dr.len, dr.start1, dr.end1, dr.start2, dr.end2), '\t')
end

scan_one(id, seq, min_len) = string(ir_part(id, seq, min_len), '\t', dr_part(seq))

function ir_part(id, seq, min_len)
    st, ir, d = find_ir_status(seq; min_len = min_len)
    ir === nothing && return emptyrow(id, length(seq), st)
    d === nothing && (d = IRDiff[])
    nsnp = count(x -> x.kind == :snp, d)
    nind = count(x -> x.kind == :indel, d)
    bpb = sum((x.len for x in d if x.kind == :indel && x.copy == :IRb); init = 0)
    bpa = sum((x.len for x in d if x.kind == :indel && x.copy == :IRa); init = 0)
    g = ir.gc
    return @sprintf("%s\t%d\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.5f\t%.4f\t%d\t%d\t%d\t%d\t%d",
                    id, ir.L, st, ir.irb[1], ir.irb[2], ir.ira[1], ir.ira[2], ir.len_irb, ir.len_ira,
                    ir.lsc, ir.ssc, g.LSC, g.IRb, g.SSC, g.IRa, ir.identity, ir.aligned_frac,
                    ir.n_unknown, nsnp, nind, bpb, bpa)
end

function scan(opt)
    fasta = joinpath(opt["dataset"], "genomes.fasta.gz")
    min_len = parse(Int, opt["min-len"])
    out = open(opt["out"], "w")
    println(out, join(HEADER, '\t'))
    batch = Pair{String,String}[]
    n = 0
    tally = Dict{String,Int}()
    function flush!()
        res = Vector{String}(undef, length(batch))
        Threads.@threads for t in eachindex(batch)
            res[t] = try
                scan_one(batch[t].first, batch[t].second, min_len)
            catch e
                @warn "IR scan failed" id = batch[t].first exception = e
                string(emptyrow(batch[t].first, length(batch[t].second), :error), '\t', DR_EMPTY)
            end
        end
        for r in res
            println(out, r)
            st = split(r, '\t')[3]
            tally[st] = get(tally, st, 0) + 1
        end
        n += length(batch)
        empty!(batch)
        @info "scanned $n genomes: " * join(("$k $v" for (k, v) in sort!(collect(tally))), ", ")
    end
    foreach_fasta(fasta) do id, seq
        push!(batch, id => seq)
        length(batch) >= 64 * Threads.nthreads() && flush!()
    end
    isempty(batch) || flush!()
    close(out)
    println("written $(opt["out"]): $n genomes; " * join(("$k $v" for (k, v) in sort!(collect(tally))), ", "))
end

function plot(opt)
    want = acckey.(split(opt["plot"], ','))
    feats = read_features_tsv(joinpath(opt["dataset"], "labels.tsv"); accession = want)
    org = Dict{String,String}()
    open(joinpath(opt["dataset"], "records.tsv")) do io
        hdr = split(readline(io), '\t')
        ia, iorg = findfirst(==("accession"), hdr), findfirst(==("organism"), hdr)
        for line in eachline(io)
            f = split(line, '\t')
            org[acckey(f[ia])] = iorg === nothing ? String(f[ia]) : String(f[iorg])
        end
    end
    seqs = Dict{String,String}()
    foreach_fasta(joinpath(opt["dataset"], "genomes.fasta.gz")) do id, seq
        acckey(id) in want && (seqs[acckey(id)] = seq)
    end
    tracks = IRTrack[]
    for a in want
        haskey(seqs, a) || (@warn "not in genomes.fasta.gz" a; continue)
        seq = seqs[a]
        st, ir, d = find_ir_status(seq; min_len = parse(Int, opt["min-len"]))
        name = get(org, a, a)
        if st != :found
            st == :low_identity && @warn "IR pair below identity threshold, drawn as not detected" a ir.identity
            push!(tracks, IRTrack(name, length(seq), nothing, Feature[], nothing))
        else
            push!(tracks, IRTrack(name, ir, get(feats, a, Feature[]), d))
        end
    end
    plot_ir_junctions(tracks, opt["svg"])
    println("written $(opt["svg"])")
end

function main(argv)
    opt = parse_args(argv)
    isempty(opt["plot"]) ? scan(opt) : plot(opt)
end
end # module IRScan

# ======================================================================
# IRCheck (ir_check.jl)
# ======================================================================
#=
ir_check.jl - diagnostics for chosen genomes: length, non-ACGT count, longest N run, the
largest inverted repeat found with relaxed settings (min_len 100, no identity cut), and
the largest direct repeat (k-mers occurring exactly twice on the same strand, runs at a
constant shift with gaps <= 300; exact shift only, so indels split a direct repeat into
pieces).

  julia PlastomeAnnotator.jl ir-check NC_xxxxx,NC_yyyyy [genomes.fasta.gz]
=#
module IRCheck
using ..PlastomeIR

"Largest direct repeat: (length, first copy start, shift to second copy)."
function largest_direct_repeat(seq::AbstractString; k::Int=31, gap_max::Int=300)
    s = Vector{UInt8}(uppercase(seq))
    codes = PlastomeIR.kmer_codes(s, k)
    pos = Dict{UInt64,Vector{Int}}()
    for i in eachindex(codes)
        v = codes[i]
        v == PlastomeIR.INVALID && continue
        push!(get!(pos, v, Int[]), i)
    end
    shifts = Dict{Int,Vector{Int}}()
    for ps in values(pos)
        length(ps) == 2 && push!(get!(shifts, ps[2] - ps[1], Int[]), ps[1])
    end
    best = (0, 0, 0)
    for (d, xs) in shifts
        sort!(xs)
        st = pr = xs[1]
        for q in 2:(length(xs) + 1)
            if q > length(xs) || xs[q] - pr > gap_max
                span = pr - st + k
                span > best[1] && (best = (span, st, d))
                q <= length(xs) && (st = xs[q])
            end
            q <= length(xs) && (pr = xs[q])
        end
    end
    return best
end

function main(argv)
    isempty(argv) && error("usage: ir-check ACC1,ACC2 [genomes.fasta.gz]")
    want = Set(first.(split.(split(argv[1], ','), '.')))
    fasta = length(argv) >= 2 ? argv[2] : joinpath(homedir(), "Assembler/refseq_plastid/dataset/genomes.fasta.gz")
    println(join(["accession", "length", "n_nonACGT", "longest_N_run", "status_relaxed", "len_irb", "len_ira",
                  "lsc", "ssc", "identity", "aligned_frac", "dr_len", "dr_start", "dr_shift"], '\t'))
    foreach_fasta(fasta) do id, seq
        String(first(split(id, '.'))) in want || return
        nn = count(c -> !(c in "ACGTacgt"), seq)
        longest = 0
        for m in eachmatch(r"[Nn]+", seq)
            longest = max(longest, length(m.match))
        end
        st, ir, _ = find_ir_status(seq; min_len = 100, min_identity = 0.0, min_aligned = 0.0)
        cols = ir === nothing ? fill("", 6) :
            string.([ir.len_irb, ir.len_ira, ir.lsc, ir.ssc, round(ir.identity; digits = 4), round(ir.aligned_frac; digits = 3)])
        dr = largest_direct_repeat(seq)
        println(join([id, length(seq), nn, longest, st, cols..., dr...], '\t'))
    end
end
end # module IRCheck

# ======================================================================
# IRDiag (ir_diag.jl)
# ======================================================================
#=
ir_diag.jl - both IR passes of find_ir_status for chosen records, with the rrn16
positions from rna_qc.tsv, to see why a record is inconsistent, low_identity or has
rrn16 copies outside the IR.

  julia PlastomeAnnotator.jl ir-diag NC_085758,NC_069639 [--dataset DIR]
=#
module IRDiag
using ..PlastomeIR
using Printf

function main(argv)
    isempty(argv) && error("usage: julia PlastomeAnnotator.jl ir-diag ACC1,ACC2,... [--dataset DIR]")
    ds = joinpath(homedir(), "Assembler/refseq_plastid/dataset")
    i = findfirst(==("--dataset"), argv)
    i === nothing || (ds = argv[i + 1])
    acckey(a) = String(first(split(a, '.')))
    want = Set(acckey.(split(argv[1], ',')))

    rrn = Dict{String,Vector{String}}()
    open(joinpath(ds, "rna_qc.tsv")) do io
        hdr = split(readline(io), '\t')
        col(n) = findfirst(==(n), hdr)
        ca, cf, cg, cs, cp, cr = col("accession"), col("feature"), col("gene"), col("segments"),
                                 col("rrna_pred"), col("reasons")
        for ln in eachline(io)
            f = split(ln, '\t')
            a = acckey(f[ca])
            (a in want && f[cf] == "rRNA") || continue
            push!(get!(() -> String[], rrn, a), string(f[cg], "(", f[cp], ") ", f[cs],
                  occursin("wrong_strand", f[cr]) ? " wrong_strand" : ""))
        end
    end

    foreach_fasta(joinpath(ds, "genomes.fasta.gz")) do id, seq
        a = acckey(id)
        a in want || return
        println("== ", id, "  length ", length(seq))
        for (lab, ms, cg) in (("pass1 (2000/3000)", 2000, 3000), ("pass2 (10000/10000)", 10000, 10000))
            st, ir, _ = PlastomeIR.ir_pass(seq; max_shift = ms, chain_gap = cg)
            if ir === nothing
                println("  ", lab, ": ", st)
            else
                @printf("  %s: %s IRb %d-%d (%d) IRa %d-%d (%d) LSC %d SSC %d id %.4f aligned %.3f drift %s\n",
                        lab, st, ir.irb..., ir.len_irb, ir.ira..., ir.len_ira, ir.lsc, ir.ssc,
                        ir.identity, ir.aligned_frac, string(ir.drift))
            end
        end
        st, ir, _ = find_ir_status(seq)
        println("  final: ", st, ir === nothing ? "" : "  IRb $(ir.irb) IRa $(ir.ira)")
        dr = find_direct_repeat(seq)
        dr === nothing || println("  direct repeat: ", dr)
        for r in get(rrn, a, String[])
            println("  rRNA ", r)
        end
    end
end
end # module IRDiag

# ======================================================================
# command line
# ======================================================================
const COMMANDS = Dict("annotate" => args -> main(args), "build-index" => args -> BuildIndex.main(args),
                      "build-training" => args -> BuildTraining.main(args),
                      "train" => args -> TrainExistence.main(args), "ir-scan" => args -> IRScan.main(args),
                      "ir-check" => args -> IRCheck.main(args), "ir-diag" => args -> IRDiag.main(args))

function cli(argv)
    if isempty(argv) || argv[1] in ("-h", "--help")
        println("usage: julia --threads auto PlastomeAnnotator.jl <command> [options]")
        println("commands: ", join(sort(collect(keys(COMMANDS))), ", "), " (details at the top of the file)")
        return
    end
    # annotate is the default when the first argument is an option (--fasta / --eval)
    startswith(argv[1], "--") && return main(argv)
    haskey(COMMANDS, argv[1]) || error("unknown command $(argv[1]); one of: $(join(sort(collect(keys(COMMANDS))), ", "))")
    COMMANDS[argv[1]](argv[2:end])
end

end # module PlastomeAnnotator

if abspath(PROGRAM_FILE) == @__FILE__
    PlastomeAnnotator.cli(ARGS)
end
