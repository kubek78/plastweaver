#!/usr/bin/env julia
# AssemblerComplete v16.49.2
# Seed-and-extend organelle assembler for short, substitution-dominated reads.
# Author of v15: user-supplied project. v16: corrected/reworked implementation.
# This is not a port of NOVOPlasty and does not guarantee resolution of long IRs.
# Julia >= 1.10. External dependency: CodecZlib (see Project.toml).

module AssemblerCompleteV16

using Printf
using Statistics
using CodecZlib

# BEGIN EMBEDDED PLASTOME ORDER (identical to OrderPlastome.jl module)
module PlastomeOrder

const BASE_CODE = let t = fill(UInt8(255), 256)
    for (i, b) in enumerate("ACGT")
        t[Int(b)+1] = UInt8(i-1)
        t[Int(lowercase(b))+1] = UInt8(i-1)
    end
    t
end
@inline function complements(a::UInt8, b::UInt8)
    x = BASE_CODE[Int(a)+1]; y = BASE_CODE[Int(b)+1]
    return x < 4 && y < 4 && (x ⊻ 0x03) == y
end
@inline circular_base(seq, p::Int) = seq[mod(p, length(seq))+1]

struct IRPair
    a::Int                         # Zero-based starts in original circular sequence.
    b::Int
    len::Int
    identity::Float64
end

struct OrderedPlastome
    sequence::Vector{UInt8}
    origin0::Int                   # Original coordinate now at FASTA position 1.
    ir_pair::IRPair
    lsc_length::Int
    ssc_length::Int
    method::String
    ssc_inverted::Bool             # SSC reverse-complemented (flip-flop isomer) to put ndhF at its start (next to the first IR)
    ndhf::String                   # SSC_start, SSC_start_after_inversion, not_found or not_checked
end

# ndhF reference proteins (NCBI RefSeq): Marchantia polymorpha YP_009479660.1, Physcomitrium
# patens YP_009477565.1, Anthoceros angustus NP_777459.1, Arabidopsis thaliana NP_051106.1,
# Nicotiana tabacum NP_054554.1. Used only to find which end of the SSC holds ndhF.
const NDHF_PROTEINS = (
    "MELIFQNVWLVPLFPFLASSLLGIGLFFFPNSIKKFRRLSSFISIMFLNIAMFLSFHFFWQQITGSPIHRYLWSWVLYKDFVLEIGYLLD" *
    "PLTSIMLVLVTTVAVMVMIYSDSYMFYDEGYIKFFCYLSLFTASMLGLVLSPNLIQVYIFWELVGMCSYLLIGFWFTRPSAANACQKAFV" *
    "TNRIGDFGLLLGILGFYWITGSFDFQQLSKRFFELLSYNQINLVFATLCALFLFLGPVAKSAQFPLHIWLPDAMEGPTPISALIHAATMV" *
    "AAGIFLVARMFPLFQMLPFVMSIISWTGAITALLGATIALAQKDLKKGLAYSTMSQLGYMMLALGIGSYKAGLFHLITHAYSKALLFLGS" *
    "GSVIHSMEPIVGYHPNKSQNMIFMGGLRQYMPITAITFLFGTLSLCGIPPFACFWSKDEILVNSWLYFPILGSIAFFTAGLTAFYMFRIY" *
    "FLTFEGDFRGHLFDDVKKFSSISIWGSLEFSEEEFKLEKKPTLYPKEANNIMLFPLIILTIPTIFIGFIGILFDDNKMNVDFLSYWLNLS" *
    "INSFNYSNSEKFFEFLFNAIPSVSIAFVGILIAFYLYGPNFSFLKKKIQLKPKIDIVLKSFSNFIYNWSYYRAYIDGFYSFFFIKGLRFL" *
    "IKIVSFIDRWIIDGFINGIGIFSFFGGESLKYIEGGRISSYLFFIIFCMFLFFLYSYII",
    "MEFIYQFVWIVPFFPFIASILIGVNLLFFPKSTKSLREIWAIFSILLLSIAMIFSFNILWQQINDNIIYRYLWSWIFDNKIDFKIGFLID" *
    "PLTSIMLVLITTVGVLVMIYSDSYMSYDQGYVRFFAYLSLFTASMLGLVLSPNLIQIYIFWELVGMCSYLLIGFWFTRPSAANACQKAFV" *
    "TNRIGDFGLLLGILGFYWITGSFEFEILFKRFNDLVINHEVNLYFANFCALLLFLGPIAKSAQFPLHIWLPDAMEGPTPISALIHAATMV" *
    "AAGIFLVARLFPLFQLLPFVMTIISWVGAITAFLGATIALAQTDLKKGLAYSTMSQLGYMMLALGIGSYQAGLFHLITHAYSKALLFLGS" *
    "GSVIHSVESIVGYSPNKCQNMAFMGGLRKYMPITGITFLLGTFSLCGIPPFACFWSKDEIITDSWLYSSTLGSISLVTAGLTAFYMFRIY" *
    "FLTFEGDLRVNLNKVAFTYPVSIWGELNLNKKNSNKNKFVLNKYDIFFETDQNENQTDNLSFFSDIKKLKYPKESDNKMLFPLLVLTLPT" *
    "LFIGFLGAPFPEGQIGSDLLSQWLYPVFKSAEEITSGNWLEFGLNAINSLSVVFSGIFIAFILYGPFSLFPQNLEKNIEFSLEKNLNSFF" *
    "SFIYNWSYFRGYIDVYYNIVFVKGTRLLAQSLSFFDEWIIDGFVNGFGILTFFEGESIKYLEGGRISFYLFGLIFGMIILLFIGFFGAMF",
    "MQLIYQCAWIVPLCPLISSILIGLGLLFFEKSTKSIRRICAISCILSLSVATVICYNLSYEQIIDNSIYQYSWSWISNGDIVLEFGYLID" *
    "PLTCIMLVLVTSVAIIVMIYSDGYMSHDEGYVRFFVYLSLFTASMLGLVLSPNLIQIYIFWELVGMCSYLLIGFWFTRPSAANACQKAFI" *
    "TNRIGDFGLLLGILGFYWITGSFEFENLFKGFKDLLINNEVSPFFAILCASFLFLGPVAKSAQFPLHVWLPDAMEGPTPISALIHAATMV" *
    "AAGIYLVARMFPLFETLPFVMSVISWTGAVTALLGATLAFFQKDLKRGLAYSTMSQLGYMMLALGIGSYRAGLFHLITHAYSKALLFLGS" *
    "GSVIHSMEPIVGYRPDESQNMIFMGGLRKHMPITGTTFLLGTLSLCGIPPFACFWSKDEIIADSWLYSPFIGWIALLTAGLTSFYMFRIY" *
    "FLTFEGEFRANSFKENTPVSSVSLWGEFRFEEFGEKKADSVLQIVEKSSPKEFFRFVQSNRESSETQIGNYLSAYPSKEPIHFLYPKESD" *
    "YTMLIPLLILSIPTLLIGFIGAPLPNGQLGSDLLSHWLNSFGNLSPERLSENWLEFIEDAITSISIASFGIFSAFILYGPASIFPRDLEK" *
    "KIEPQLKGVWGFFINYMYNWSQYRGYIDQYYNKIFVEGTRILAYASSFFDRWIIDGIVNGTGISSFSGGEGMRYGEGGRVSSYLFGLVFG" *
    "MTILLAVILLII",
    "MEHTYQYSWIIPFIPLPVPILLGVGLLLFPTATKNLRRMWTFLSIFLLSIVMIFSIYLSIQQIFLSCIHQNVWSWTINNEFSFEFGYFID" *
    "PLTSIMSILITTVGILVLIYSDNYMSHDQGYLRFFAYMGFFNTSMLGLVTSSNLIQVYFFWELVGMCSYLLIGFWFTRPIAANACQKAFV" *
    "TNRVGDFGLLLGILGLYWITGSFEFQDLFEIFNNLILNNRVNLLFLTLCAFLLFVGPIAKSAQFPLHVWLPDAMEGPTPISALIHAATMV" *
    "AAGIFLVARLLPLFIVIPSIMYIISLIGIITVLLGATLALAQKDIKRGLAYSTMSQLGYMMLALGMGSYRSALFHLITHAYSKALLFLGS" *
    "GSIIHSMEAIVGYSPDKSQNMILMGGLTKHVPITKTAFLIGTLSLCGIPPLACFWSKDEILNDSLLFSPIFAIIACSTAGLTAFYMFRIY" *
    "LLTFEGHLNTYFLNYSGKKSGSFYSLSLWGKEEEKKLNKNFGLVPLLTMNNTKRASFFCNKTYKISNNVRNQIFITVENFGLNTRTFYYP" *
    "HESDNTILFPMLILVLFTLFIGAIGIPFNQEGIDFDILSKFFTPSINLLHKNSQNFVDWYEFLRNATFSVSIAFFGIFIAYCLYKPFYSS" *
    "LLNLTLLNSFQKWNSKRIHWEKLINFVYNWSYNRGYIDSFFKTSLIESIRRLAKQTTFFDKRIIDGITNGVGITSFFVGEVTKYIGGSRI" *
    "SSYLFLYLSYVLIFLMILFFFYFEKF",
    "MEQTYEYAWIIPFIPLPVPMLIGAGLFLFPTATKSFRRMWAFQSVLLLSIVMVFSIYLSIQQINSSSFYQYVWSWIINNDFSLDFGYLID" *
    "PLTSIMSILITTVGIMVLIYSDNYMAHDQGYLRFFAYMSFFSTSMLGLVTSSNLIQIYIFWELVGLCSYLLIGFWFTRPVAANACQKAFV" *
    "TNRVGDFGLLLGILGFYWITGSFEFRDLFEIFNNLIYNNEVDFLFVTLCAVLLFAGAVAKSAQFPLHVWLPDAMEGPTPISALIHAATMV" *
    "AAGIFLVARLLPLFRVIPYIMYLISVIGIITVLLGATLALAQKDIKRGLAYSTMSQLGYMMLALGMGSYRSALFHLITHAYSKALLFLGS" *
    "GSIIHSMETIVGYSPAKSQNMGLMGGLRKHVPISKITFLLGTLSLCGIPPLACFWSKDEILNDSWLYSPIFAIIAWATAGLTAFYMFRIY" *
    "LLTFEGHLNAHFPNYGGKQKTPFYSISLWGKNGVKKNSCLLTMNNNESTYFFAKTKYPIDKNGRKMTRPFMTIAHFEHKAVYSYPYESDN" *
    "TMLFPIFVLGLFTLFVGSIGIPFNQEGGNLDILSKWLAPSINLLHQKSNNSMDWNEFLKDAVLSVSIAYFGIFIASFLYKPIYSSLKNFE" *
    "LINSFVKKGPKRILWDKIINGIYDWSYNRAYIDAFYTRFLVGGIRGLAEFTHFFDRRVIDGMTNGVGVISFIVGEGIKYIGGGRISSYLF" *
    "LYLAYVSIFLLVYYLLFSTL",
)
const NDHF_K = 5
const NDHF_KMERS = let set = Set{String}()
    for prot in NDHF_PROTEINS, i in 1:length(prot)-NDHF_K+1
        push!(set, prot[i:i+NDHF_K-1])
    end
    set
end
const GENETIC_CODE = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"   # codons in TCAG order
const TCAG_CODE = let t = fill(UInt8(255), 256)
    for (i, b) in enumerate("TCAG")
        t[Int(b)+1] = UInt8(i-1); t[Int(lowercase(b))+1] = UInt8(i-1)
    end
    t
end
const RC_BASE = let t = fill(UInt8('N'), 256)
    for (x, y) in zip("ACGTacgt", "TGCAtgca"); t[Int(x)+1] = UInt8(y); end
    t
end
revcomp(s::AbstractVector{UInt8}) = UInt8[RC_BASE[Int(b)+1] for b in Iterators.reverse(s)]

"Translation of `s` from 0-based offset `f` (X for codons with non-ACGT bases)."
function translate_frame(s::AbstractVector{UInt8}, f::Int)
    n = max((length(s) - f) ÷ 3, 0)
    out = Vector{UInt8}(undef, n)
    for i in 1:n
        p = f + 3(i-1)
        a = TCAG_CODE[Int(s[p+1])+1]; b = TCAG_CODE[Int(s[p+2])+1]; c = TCAG_CODE[Int(s[p+3])+1]
        out[i] = (a == 255 || b == 255 || c == 255) ? UInt8('X') : UInt8(GENETIC_CODE[16Int(a)+4Int(b)+Int(c)+1])
    end
    return String(out)
end

"0-based positions in `s` of reference ndhF protein 5-mers in the six-frame translation."
function ndhf_hits(s::AbstractVector{UInt8})
    out = Int[]; n = length(s)
    for (forward, t) in ((true, s), (false, revcomp(s))), f in 0:2
        aa = translate_frame(t, f)
        for i in 1:length(aa)-NDHF_K+1
            aa[i:i+NDHF_K-1] in NDHF_KMERS || continue
            nt = f + 3(i-1)
            push!(out, forward ? nt : n - 1 - nt)
        end
    end
    return sort!(out)
end

"""
End of the SSC holding ndhF in an LSC-IR-SSC-IR sequence: :start, :end or :absent. ndhF is the
densest 3-kb window of reference-protein 5-mer hits over the SSC plus <= 2 kb of each IR, with
>= 100 hits (a complete ndhF gives several hundred; its homologues ndhD and ndhB about 30).
"""
function ndhf_side(seq::AbstractVector{UInt8}, lsc::Int, ir::Int, ssc::Int; flank::Int=2000,
                   window::Int=3000, min_hits::Int=100)
    fl = min(flank, ir)
    lo = lsc + ir - fl; hi = min(lsc + ir + ssc + fl, length(seq))
    h = ndhf_hits(@view seq[lo+1:hi])
    isempty(h) && return :absent
    best = 0; bi = 1; j = 1
    for i in eachindex(h)
        while h[i] - h[j] > window
            j += 1
        end
        if i - j + 1 > best
            best = i - j + 1; bi = i
        end
    end
    best >= min_hits || return :absent
    centre = h[bi - best + cld(best, 2)]
    return centre < fl + ssc ÷ 2 ? :start : :end
end

function pair_gaps(a::Int, b::Int, len::Int, n::Int)
    d = mod(b-a, n)
    return d-len, n-d-len
end

function pair_identity(seq::Vector{UInt8}, a::Int, b::Int, len::Int)
    matches = 0
    for t in 0:len-1
        matches += complements(circular_base(seq,a+t), circular_base(seq,b+len-1-t))
    end
    return matches/len
end

"Extend only exact terminal matches; retain internal substitutions unchanged."
function extend_pair(seq::Vector{UInt8}, a::Int, b::Int, len::Int, min_sc::Int)
    n = length(seq)
    while true
        _, before = pair_gaps(a,b,len,n)
        before >= min_sc+2 || break
        complements(circular_base(seq,a-1), circular_base(seq,b+len)) || break
        a = mod(a-1,n); len += 1
    end
    while true
        after, _ = pair_gaps(a,b,len,n)
        after >= min_sc+2 || break
        complements(circular_base(seq,a+len), circular_base(seq,b-1)) || break
        b = mod(b-1,n); len += 1
    end
    return min(a,b), max(a,b), len
end

"Circular canonical k-mer index, including kmers crossing the FASTA origin."
function circular_seed_index(seq::Vector{UInt8}, k::Int, max_occurrences::Int)
    buckets = Dict{UInt64,Vector{Tuple{Int,Bool}}}()
    n = length(seq)
    mask = (UInt64(1) << (2k))-UInt64(1)
    f::UInt64 = 0; r::UInt64 = 0; valid_run = 0
    for t in 0:n+k-2
        base = BASE_CODE[Int(circular_base(seq,t))+1]
        if base == 255
            f = 0; r = 0; valid_run = 0
        else
            f = ((f << 2) | UInt64(base)) & mask
            r = (r >> 2) | (UInt64(base ⊻ 0x03) << (2(k-1)))
            valid_run += 1
            if valid_run >= k
                p = t-k+1
                key = min(f,r)
                bucket = get!(() -> Tuple{Int,Bool}[], buckets, key)
                # max+1 entries mark a repetitive, unusable seed. Never grow it again.
                length(bucket) <= max_occurrences && push!(bucket,(p,f < r))
            end
        end
    end
    return buckets
end

function anchor_groups(points::Vector{Int}, n::Int, max_gap::Int)
    p = sort!(unique(points))
    length(p) < 2 && return Vector{Vector{Int}}()
    gaps = [p[i+1]-p[i] for i in 1:length(p)-1]
    push!(gaps,p[1]+n-p[end])
    maximum(gaps) > max_gap || return Vector{Vector{Int}}()
    start = mod1(argmax(gaps)+1,length(p))
    first_p = p[start]
    unwrapped = Int[]
    for offset in 0:length(p)-1
        x = p[mod1(start+offset,length(p))]
        push!(unwrapped, x < first_p ? x+n : x)
    end
    groups = Vector{Vector{Int}}()
    current = Int[unwrapped[1]]
    for x in unwrapped[2:end]
        if x-current[end] > max_gap
            push!(groups,current); current = Int[]
        end
        push!(current,x)
    end
    push!(groups,current)
    return groups
end

function circular_intersection(a::Int, la::Int, b::Int, lb::Int, n::Int)
    function pieces(start,len)
        last = start+len
        return last <= n ? [(start,last)] : [(start,n),(0,last-n)]
    end
    return sum(max(0,min(x[2],y[2])-max(x[1],y[1])) for x in pieces(a,la) for y in pieces(b,lb))
end

function same_pair_locus(x::IRPair, y::IRPair, n::Int)
    cutoff = 0.9*min(x.len,y.len)
    direct = circular_intersection(x.a,x.len,y.a,y.len,n) >= cutoff &&
             circular_intersection(x.b,x.len,y.b,y.len,n) >= cutoff
    swapped = circular_intersection(x.a,x.len,y.b,y.len,n) >= cutoff &&
              circular_intersection(x.b,x.len,y.a,y.len,n) >= cutoff
    return direct || swapped
end

"Find a dominant pair of long inverted repeats by circular reflection diagonals."
function detect_ir(seq::Vector{UInt8}; min_ir::Int=1000, min_sc::Int=100,
                   min_identity::Float64=0.995, k::Int=31,
                   max_seed_gap::Int=128, max_occurrences::Int=16,
                   max_diagonals::Int=16)
    1 <= k <= 31 && isodd(k) || error("IR seed k must be odd and in 1:31")
    min_ir >= k && min_sc > 0 || error("Invalid minimum IR or single-copy length")
    0 < min_identity <= 1 || error("IR identity must be in (0,1]")
    n = length(seq)
    n >= 2*min_ir+2*min_sc || return nothing,"sequence_too_short_for_four_regions"
    buckets = circular_seed_index(seq,k,max_occurrences)
    votes = Dict{Int,Int}()
    for bucket in values(buckets)
        2 <= length(bucket) <= max_occurrences || continue
        for i in 1:length(bucket)-1, j in i+1:length(bucket)
            p, strand_p = bucket[i]; q, strand_q = bucket[j]
            strand_p == strand_q && continue
            diagonal = mod(p+q+k-1,n)
            votes[diagonal] = get(votes,diagonal,0)+1
        end
    end
    isempty(votes) && return nothing,"no_inverted_repeat_seeds"
    ranked = sort!(collect(votes); by=x -> (-last(x),first(x)))
    ranked = ranked[1:min(length(ranked),max_diagonals)]
    anchors = Dict(first(x) => Int[] for x in ranked)
    for bucket in values(buckets)
        2 <= length(bucket) <= max_occurrences || continue
        for i in 1:length(bucket)-1, j in i+1:length(bucket)
            p, strand_p = bucket[i]; q, strand_q = bucket[j]
            strand_p == strand_q && continue
            diagonal = mod(p+q+k-1,n)
            haskey(anchors,diagonal) || continue
            push!(anchors[diagonal],p,q)
        end
    end
    candidates = IRPair[]
    seen = Set{Tuple{Int,Int,Int}}()
    # A legal SC interval must always separate anchor groups.
    gap_limit = min(max_seed_gap,min_sc+k-1)
    for (diagonal,points) in sort!(collect(anchors); by=first)
        for group in anchor_groups(points,n,gap_limit)
            length(group) >= 3 || continue
            len = group[end]-group[1]+k
            len >= min_ir && 2*len+2*min_sc <= n || continue
            a = mod(group[1],n)
            b = mod(diagonal-(a+len-1),n)
            minimum(pair_gaps(a,b,len,n)) >= min_sc || continue
            a,b,len = extend_pair(seq,a,b,len,min_sc)
            (a,b,len) in seen && continue
            push!(seen,(a,b,len))
            identity = pair_identity(seq,a,b,len)
            identity+1e-12 >= min_identity || continue
            push!(candidates,IRPair(a,b,len,identity))
        end
    end
    isempty(candidates) && return nothing,"no_long_ir_pair_passing_thresholds"
    sort!(candidates; by=x -> (-x.len,-x.identity,x.a,x.b))
    best = candidates[1]
    for other in candidates[2:end]
        if other.len >= 0.9*best.len && !same_pair_locus(best,other,n)
            return nothing,"multiple_competing_ir_pairs"
        end
    end
    return best,"detected"
end

"Only rotate. Do not duplicate/delete IR, invert SSC or reverse-complement the genome."
function order_plastome(seq::Vector{UInt8}; min_ir::Int=1000, min_sc::Int=100,
                        min_identity::Float64=0.995,
                        manual_starts::Union{Nothing,Tuple{Int,Int}}=nothing,
                        manual_length::Int=0, orient_ssc::Bool=true)
    pair = nothing; method = "detected"
    if manual_starts === nothing
        pair,method = detect_ir(seq;min_ir=min_ir,min_sc=min_sc,min_identity=min_identity)
        pair === nothing && return nothing,method
    else
        n = length(seq)
        all(x -> 1 <= x <= n,manual_starts) || error("Manual IR starts must be 1-based positions within the sequence")
        manual_length >= min_ir || error("Manual IR length is shorter than Min IR length")
        a,b = sort(collect(manual_starts).-1)
        minimum(pair_gaps(a,b,manual_length,n)) >= min_sc || error("Manual IR intervals overlap or leave too short single-copy regions")
        identity = pair_identity(seq,a,b,manual_length)
        identity+1e-12 >= min_identity || error("Manual IR coordinates fail the identity threshold")
        pair = IRPair(a,b,manual_length,identity); method = "manual_coordinates"
    end
    gap1,gap2 = pair_gaps(pair.a,pair.b,pair.len,length(seq))
    gap1 == gap2 && return nothing,"single_copy_regions_have_equal_length"
    if gap1 > gap2
        origin = mod(pair.a+pair.len,length(seq)); lsc = gap1; ssc = gap2
    else
        origin = mod(pair.b+pair.len,length(seq)); lsc = gap2; ssc = gap1
    end
    rotated = vcat(seq[origin+1:end],seq[1:origin])
    length(rotated) == length(seq) || error("Internal error: genome length changed")
    lsc+ssc+2*pair.len == length(seq) || error("Internal error: region lengths do not sum to genome length")
    # This postcondition also verifies IR orientation after the coordinate transform.
    first_ir = @view rotated[lsc+1:lsc+pair.len]
    second_ir = @view rotated[end-pair.len+1:end]
    matches = count(i -> complements(first_ir[i],second_ir[end-i+1]),eachindex(first_ir))
    matches/pair.len+1e-12 >= min_identity || error("Internal error: ordered IRs do not match")
    # SSC orientation: ndhF at the SSC start, next to the first IR (LSC-IRa-SSC-IRb naming: IRa),
    # on the - strand (its 3' end faces the IR/SSC junction). The other orientation is the
    # flip-flop isomer; inverting the SSC between identical IRs leaves the IRs unchanged.
    inverted = false; ndhf = "not_checked"
    if orient_ssc
        side = ndhf_side(rotated, lsc, pair.len, ssc)
        if side == :end
            s0 = lsc + pair.len
            rotated[s0+1:s0+ssc] = revcomp(rotated[s0+1:s0+ssc])
            inverted = true
        end
        ndhf = side == :start ? "SSC_start" : side == :end ? "SSC_start_after_inversion" : "not_found"
    end
    return OrderedPlastome(rotated,origin,pair,lsc,ssc,method,inverted,ndhf),"ordered"
end

const REGIONS_HEADER = "id\tregion\tstart\tend\tlength\toriginal_start\toriginal_end\twraps_original_origin\tir_identity\tmethod\tstrand"

function write_region_rows(io::IO, id::String, item::OrderedPlastome; revcomp::Bool=false)
    n = length(item.sequence); offset = 0
    regions = (("LSC",item.lsc_length),("IR1",item.ir_pair.len),
               ("SSC",item.ssc_length),("IR2",item.ir_pair.len))
    for (name,len) in regions
        old_start = mod(item.origin0+offset,n)+1
        old_end = mod(item.origin0+offset+len-1,n)+1
        wraps = old_start>old_end
        minus = name == "SSC" && item.ssc_inverted     # ordered start maps to the original end
        minus && ((old_start,old_end) = (old_end,old_start))
        if revcomp                                   # item was ordered from the reverse complement
            old_start, old_end = n - old_start + 1, n - old_end + 1; minus = !minus
        end
        println(io,join((id,name,offset+1,offset+len,len,old_start,old_end,
                        wraps,round(item.ir_pair.identity;digits=6),item.method,minus ? "-" : "+"),'\t'))
        offset += len
    end
end

function write_ordered_record(io::IO, id::String, item::OrderedPlastome; closure::AbstractString="",
                              extra::AbstractString="", revcomp::Bool=false)
    n = length(item.sequence)
    tag = (isempty(closure) ? "" : " closure=$closure") * (isempty(extra) ? "" : " $extra")
    start = revcomp ? n - item.origin0 : item.origin0 + 1
    transform = (revcomp ? "reverse_complement+" : "") * (item.ssc_inverted ? "rotation+SSC_inversion" : "rotation")
    println(io,">$id length=$n topology=circular$tag order=LSC-IR-SSC-IR original_start=$start transform=$transform ndhF=$(item.ndhf)")
    for i in 1:80:n
        println(io,String(item.sequence[i:min(i+79,n)]))
    end
end

function fasta_records(path::String)
    records = Pair{String,Vector{UInt8}}[]
    header = ""; seq = UInt8[]; ids = Set{String}()
    for raw in eachline(path)
        line = strip(raw); isempty(line) && continue
        if startswith(line,">")
            if !isempty(header)
                isempty(seq) && error("Empty FASTA record: $header")
                push!(records,header=>seq)
            end
            header = String(strip(line[2:end])); seq = UInt8[]
            isempty(header) && error("Empty FASTA header")
            id = String(first(split(header)))
            id in ids && error("Duplicate FASTA identifier: $id")
            push!(ids,id)
        else
            isempty(header) && error("Sequence before the first FASTA header")
            for b in codeunits(line)
                uppercase(Char(b)) in "ACGTNRYKMSWBDHV" || error("Invalid DNA character")
                push!(seq,b)
            end
        end
    end
    if !isempty(header)
        isempty(seq) && error("Empty FASTA record: $header")
        push!(records,header=>seq)
    end
    isempty(records) && error("No FASTA records")
    return records
end

is_circular_header(header::String) = occursin(r"(?:^|\s)topology=circular(?:\s|$)",header)

function cli(args=ARGS)
    if length(args)<2 || any(x -> x in ("-h","--help"),args)
        println("Usage: julia OrderPlastome.jl INPUT.fasta OUTPUT.fasta [--contig ID]")
        println("Options: --min-ir 1000 --min-sc 100 --identity 0.995")
        println("Manual boundaries: --ir-starts A,B --ir-length L (1-based original coordinates)")
        return
    end
    input,output = abspath(args[1]),abspath(args[2])
    (input == output || (isfile(output) && realpath(input)==realpath(output))) && error("Output must differ from input")
    opts = Dict("--contig"=>"","--min-ir"=>"1000","--min-sc"=>"100","--identity"=>"0.995",
                "--ir-starts"=>"","--ir-length"=>"0")
    seen = Set{String}()
    i = 3
    while i <= length(args)
        key = args[i]
        haskey(opts,key) && i<length(args) || error("Unknown or incomplete option: $key")
        key in seen && error("Duplicate option: $key")
        opts[key] = args[i+1]; push!(seen,key); i += 2
    end
    records = fasta_records(input)
    selected = isempty(opts["--contig"]) ? [r for r in records if is_circular_header(first(r))] :
               [r for r in records if first(split(first(r))) == opts["--contig"]]
    length(selected)==1 || error("Select exactly one circular contig with --contig ID; do not select by length")
    header,seq = only(selected)
    occursin(r"(?:^|\s)topology=linear(?:\s|$)",header) && error("Selected contig is explicitly labelled linear")
    starts = nothing
    if !isempty(opts["--ir-starts"])
        parts = parse.(Int,split(opts["--ir-starts"],','))
        length(parts)==2 || error("--ir-starts requires two coordinates separated by a comma")
        starts = (parts[1],parts[2])
    end
    manual_length = parse(Int,opts["--ir-length"])
    (starts===nothing) == (manual_length==0) || error("Use --ir-starts and --ir-length together")
    min_ir = parse(Int,opts["--min-ir"]); min_sc = parse(Int,opts["--min-sc"])
    identity = parse(Float64,opts["--identity"])
    min_ir >= 31 && min_sc > 0 && 0 < identity <= 1 || error("Invalid ordering thresholds")
    item,status = order_plastome(seq;min_ir=min_ir,min_sc=min_sc,min_identity=identity,
                                 manual_starts=starts,manual_length=manual_length)
    item===nothing && error("Cannot assign LSC/IR/SSC automatically: $status. Supply verified IR coordinates if needed.")
    mkpath(dirname(output))
    id = String(first(split(header)))
    open(output,"w") do io
        write_ordered_record(io,id,item)
    end
    open(output*".regions.tsv","w") do io
        println(io,REGIONS_HEADER); write_region_rows(io,id,item)
    end
    println("Saved: $output")
    println("Regions: LSC=$(item.lsc_length), IR1=$(item.ir_pair.len), SSC=$(item.ssc_length), IR2=$(item.ir_pair.len)")
    println("Original start=$(item.origin0+1); length=$(length(seq)); IR identity=$(round(item.ir_pair.identity;digits=6))")
end

end # module PlastomeOrder

"""
Nearest land-plant reference by FracMinHash: canonical 21-mers hashed with splitmix64, kept
below 2^64/scaled (scaled 500: about 300 hashes per plastome). The reference set is the dataset
directory written by build_landplant_dataset.py (genomes.fasta.gz, records.tsv, features.tsv);
its sketches are built once:  julia --project=. AssemblerComplete_v16.jl --build-reference-sketches <dataset>
"""
module RefSketch

using CodecZlib

const K = 21
const SCALED = 500
const MAGIC = Vector{UInt8}(codeunits("PLSK1\n"))
const SKETCH_FILE = "ref_sketches.bin"
const ENC2 = let t = fill(UInt8(255), 256)
    for (i, b) in enumerate("ACGT")
        t[Int(b)+1] = UInt8(i-1); t[Int(lowercase(b))+1] = UInt8(i-1)
    end
    t
end

@inline function mix64(x::UInt64)
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

"Sorted unique FracMinHash of the canonical k-mers of `seq` (across the origin if circular)."
function sketch(seq::AbstractVector{UInt8}; k::Int=K, scaled::Int=SCALED, circular::Bool=false)
    s = circular && length(seq) >= k ? vcat(seq, seq[1:k-1]) : seq
    thr = typemax(UInt64) ÷ UInt64(scaled)
    mask = (UInt64(1) << (2k)) - UInt64(1)
    f = UInt64(0); r = UInt64(0); run = 0; out = UInt64[]
    @inbounds for b in s
        c = ENC2[Int(b)+1]
        if c == 0xff
            f = UInt64(0); r = UInt64(0); run = 0; continue
        end
        f = ((f << 2) | UInt64(c)) & mask
        r = (r >> 2) | (UInt64(0x03 ⊻ c) << (2(k-1)))
        run += 1
        if run >= k
            h = mix64(min(f, r))
            h < thr && push!(out, h)
        end
    end
    return sort!(unique!(out))
end

"Number of hashes shared by two sorted vectors."
function shared(a::Vector{UInt64}, b::Vector{UInt64})
    i = 1; j = 1; n = 0
    @inbounds while i <= length(a) && j <= length(b)
        x = a[i]; y = b[j]
        if x == y
            n += 1; i += 1; j += 1
        elseif x < y
            i += 1
        else
            j += 1
        end
    end
    return n
end

"Calls f(name, seq) for every FASTA record (gz or plain) until f returns true."
function each_fasta(f::Function, path::AbstractString)
    io = endswith(path, ".gz") ? GzipDecompressorStream(open(path)) : open(path)
    try
        name = ""; buf = UInt8[]
        for line in eachline(io)
            if startswith(line, ">")
                if !isempty(name) && f(name, buf)
                    return
                end
                name = String(first(split(line[2:end]))); buf = UInt8[]
            else
                append!(buf, codeunits(strip(line)))
            end
        end
        isempty(name) || f(name, buf)
    finally
        close(io)
    end
    return
end

const LOCAL_DIR = "local"
const FASTA_EXT = r"\.(fasta|fa|fna)(\.gz)?$"i

"Own plastomes in <dataset>/local/<library>.fasta: (stem, path) pairs."
function local_files(dataset::AbstractString)
    dir = joinpath(dataset, LOCAL_DIR)
    isdir(dir) || return Tuple{String,String}[]
    return [(replace(f, FASTA_EXT => ""), joinpath(dir, f)) for f in sort(readdir(dir)) if occursin(FASTA_EXT, f)]
end

"Name of a local reference: local:<file stem>:<record>."
local_name(stem, record) = "local:" * stem * ":" * record

function build(dataset::AbstractString; k::Int=K, scaled::Int=SCALED)
    fasta = joinpath(dataset, "genomes.fasta.gz")
    isfile(fasta) || error("Missing $fasta (run build_landplant_dataset.py first)")
    out = joinpath(dataset, SKETCH_FILE); n = 0; nlocal = 0
    open(out * ".part", "w") do io
        write(io, MAGIC); write(io, Int32(k)); write(io, Int64(scaled))
        function put(name, seq)
            h = sketch(seq; k=k, scaled=scaled, circular=true)
            a = Vector{UInt8}(codeunits(name))
            write(io, UInt16(length(a))); write(io, a)
            write(io, Int64(length(seq))); write(io, Int32(length(h))); write(io, h)
            n += 1
            n % 2000 == 0 && println("  sketched ", n)
        end
        each_fasta(fasta) do name, seq
            put(name, seq); false
        end
        # Own plastomes (e.g. validated circles of other libraries) as additional references.
        for (stem, path) in local_files(dataset)
            each_fasta(path) do name, seq
                put(local_name(stem, name), seq); nlocal += 1; false
            end
        end
    end
    nlocal > 0 && println("Local references: ", nlocal, " record(s) from ", joinpath(dataset, LOCAL_DIR))
    mv(out * ".part", out; force=true)
    println("Reference sketches: ", n, " genomes, k=", k, ", scaled=", scaled, " -> ", out)
    return out
end

struct SketchDB
    k::Int
    scaled::Int
    names::Vector{String}
    lengths::Vector{Int}
    hashes::Vector{Vector{UInt64}}
end

function load(path::AbstractString)
    open(path) do io
        read(io, length(MAGIC)) == MAGIC || error("Not a reference sketch file: $path")
        k = Int(read(io, Int32)); scaled = Int(read(io, Int64))
        names = String[]; lens = Int[]; hs = Vector{UInt64}[]
        while !eof(io)
            l = Int(read(io, UInt16)); push!(names, String(read(io, l)))
            push!(lens, Int(read(io, Int64))); m = Int(read(io, Int32))
            v = Vector{UInt64}(undef, m); read!(io, v); push!(hs, v)
        end
        return SketchDB(k, scaled, names, lens, hs)
    end
end

"accession => (organism, group) from records.tsv (empty if absent)."
function metadata(dataset::AbstractString)
    path = joinpath(dataset, "records.tsv"); out = Dict{String,Tuple{String,String}}()
    isfile(path) || return out
    for (i, line) in enumerate(eachline(path))
        i == 1 && continue
        f = split(line, '\t')
        length(f) >= 4 && (out[String(f[1])] = (String(f[4]), String(f[3])))
    end
    for (stem, path) in local_files(dataset)
        each_fasta(path) do name, _
            out[local_name(stem, name)] = (stem, "local"); false
        end
    end
    return out
end

"""
Best references for a query (several sequences): containment = shared / query hashes (the share
of the query found in the reference), ani ~ containment^(1/k), reference_covered = shared /
reference hashes (the share of the reference present in the query).
"""
function nearest(db::SketchDB, seqs, circular; top::Int=5, exclude::Function=name -> false)
    q = UInt64[]
    for (s, circ) in zip(seqs, circular)
        append!(q, sketch(s; k=db.k, scaled=db.scaled, circular=circ))
    end
    sort!(unique!(q))
    isempty(q) && return NamedTuple[], 0
    scores = [(i, shared(q, h)) for (i, h) in enumerate(db.hashes) if !exclude(db.names[i])]
    sort!(scores; by=x -> (-x[2], db.names[x[1]]))
    hits = [(accession=db.names[i], shared=n, containment=n/length(q), ani=(n/length(q))^(1/db.k),
             reference_covered=n/max(1, length(db.hashes[i])), reference_length=db.lengths[i])
            for (i, n) in first(scores, top)]
    return hits, length(q)
end

"Sequence of one accession from genomes.fasta.gz, or of a local:<stem>:<record> reference (empty if absent)."
function reference_sequence(dataset::AbstractString, accession::AbstractString)
    out = UInt8[]
    if startswith(accession, "local:")
        for (stem, path) in local_files(dataset)
            each_fasta(path) do name, seq
                local_name(stem, name) == accession || return false
                append!(out, seq); true
            end
            isempty(out) || return out
        end
        return out
    end
    each_fasta(joinpath(dataset, "genomes.fasta.gz")) do name, seq
        name == accession || return false
        append!(out, seq); true
    end
    return out
end

"features.tsv rows of one accession, header first."
function reference_features(dataset::AbstractString, accession::AbstractString)
    path = joinpath(dataset, "features.tsv"); rows = String[]
    isfile(path) || return rows
    tag = accession * "\t"
    for (i, line) in enumerate(eachline(path))
        (i == 1 || startswith(line, tag)) && push!(rows, line)
    end
    return rows
end

end # module RefSketch
# END EMBEDDED PLASTOME ORDER

const VERSION_STRING = "16.49.2"
const DNA = UInt8['A', 'C', 'G', 'T']
const ENCODE = let a = fill(UInt8(255), 256)
    for (i, b) in enumerate(DNA)
        a[Int(b) + 1] = UInt8(i - 1)
        a[Int(lowercase(Char(b))) + 1] = UInt8(i - 1)
    end
    a
end
const COMPLEMENT = let a = fill(UInt8('N'), 256)
    for (x, y) in zip("ACGTN", "TGCAN")
        a[Int(x) + 1] = UInt8(y)
    end
    a
end

# Live reporting is deliberately separate from sequence/branch decisions.
# libuv is bundled with Julia; RSS includes the Julia runtime and native buffers.
function resident_memory_bytes()
    bytes = Ref{Csize_t}(0)
    try
        code = ccall(:uv_resident_set_memory, Cint, (Ref{Csize_t},), bytes)
        return code == 0 ? Int(bytes[]) : nothing
    catch
        return nothing # Unavailable is not zero, and peak RSS is not current RSS.
    end
end

function duration_text(seconds::Real)
    n = floor(Int, max(0, seconds))
    return @sprintf("%02d:%02d:%02d", n ÷ 3600, (n % 3600) ÷ 60, n % 60)
end

mutable struct ProgressMonitor
    lock::ReentrantLock
    io::IO
    terminal::Bool
    interval::Float64
    started::UInt64
    phase_started::UInt64
    last_render::UInt64
    last_sample::UInt64
    phase::String
    done::Int
    total::Int                         # Zero means unknown, not a 0% estimate.
    unit::String
    detail::String
    current_rss::Union{Nothing,Int}
    peak_rss::Int                      # Maximum observed sample, explicitly labelled.
    line_open::Bool
    stopped::Bool
    timer::Union{Nothing,Timer}
    seeds::Dict{String,NTuple{3,Int}}   # longest length, completed paths, spawned paths
    ascii_panel::Bool
    panel_interval::Float64
    last_panel::UInt64
    genome_range::Tuple{Int,Int}
    longest::Int                       # Longest path length reported by any seed.
    paths_done::Int
    paths_found::Int
    closed_paths::Int                  # Paths that passed circular validation (before deduplication).
    circles::Int                       # Final validated circles; -1 until export.
    panel_lines::Int                   # Panel rows currently on screen (terminal only).
    panel_phase::String                # Last phase that received a full panel in redirected output.
    previous_logger::Union{Nothing,Base.CoreLogging.AbstractLogger}
    label::String                      # Sample/project name shown in the panel title.
end

elapsed_seconds(m::ProgressMonitor) = Float64(time_ns()-m.started)/1e9
function progress_line(m::ProgressMonitor, stamp::UInt64)
    elapsed = Float64(stamp-m.started)/1e9
    phase_elapsed = Float64(stamp-m.phase_started)/1e9
    width = m.terminal ? 10 : 16
    if m.total > 0
        filled = clamp(floor(Int,width*m.done/m.total),0,width)
        bar = "["*repeat("=",filled)*repeat(".",width-filled)*"]"
        count = "$(m.done)/$(m.total) $(m.unit)"
    else
        cursor = mod(floor(Int,elapsed),width)
        bar = "["*repeat(".",cursor)*">"*repeat(".",width-cursor-1)*"]"
        count = m.done > 0 ? "$(m.done) $(m.unit)" : "working"
    end
    ram = m.current_rss === nothing ? "n/a" : @sprintf("%.2f GB",m.current_rss/1e9)
    detail = isempty(m.detail) ? "" : " | "*m.detail
    if m.terminal
        return "$(m.phase) $bar | t=$(duration_text(elapsed)) | RAM $ram | $count$detail"
    end
    return "$(m.phase) $bar $count$detail | stage $(duration_text(phase_elapsed)) | total $(duration_text(elapsed)) | RAM $ram"
end

# ---------------------------------------------------------------------------
# ASCII panel. Plain 7-bit ASCII so it survives pipes, log files and any font.
# On a terminal it is redrawn in place (cursor up + erase); in redirected output
# a full panel is written only at phase boundaries and at the end.
# ---------------------------------------------------------------------------
const PANEL_MIN_WIDTH = 60
const PANEL_MAX_WIDTH = 80
const PANEL_LOG_WIDTH = 76
const RING_ROWS = 7
const RING_COLS = 23
# Ring cells (row, column), clockwise from 12 o'clock; 26 cells ~ 3.8% each.
const RING_CELLS = ((1,13),(1,15),(1,17),(2,18),(2,20),(3,22),(4,23),(5,22),(6,20),
                    (6,18),(7,17),(7,15),(7,13),(7,11),(7,9),(7,7),(6,6),(6,4),(5,2),
                    (4,1),(3,2),(2,4),(2,6),(1,7),(1,9),(1,11))

function ring_rows(fraction::Float64, centre::AbstractString)
    grid = [fill(' ', RING_COLS) for _ in 1:RING_ROWS]
    filled = clamp(round(Int, fraction*length(RING_CELLS)), 0, length(RING_CELLS))
    for (i, (r, col)) in enumerate(RING_CELLS)
        grid[r][col] = i <= filled ? '#' : '.'
    end
    text = first(centre, RING_COLS-6)
    start = (RING_COLS-length(text))÷2 + 1
    for (i, ch) in enumerate(text)
        grid[(RING_ROWS+1)÷2][start+i-1] = ch
    end
    return [String(row) for row in grid]
end

fit_text(text::AbstractString, width::Int) =
    width <= 0 ? "" : rpad(first(text, width), width)

function stage_bar_text(m::ProgressMonitor, width::Int, stamp::UInt64)
    width = max(width, 4)
    if m.total > 0
        filled = clamp(floor(Int, width*m.done/m.total), 0, width)
        return "["*repeat("#",filled)*repeat(".",width-filled)*"]"
    end
    # Unknown total: a moving marker, never a fabricated fraction.
    cursor = mod(floor(Int, Float64(stamp-m.started)/1e9), width)
    return "["*repeat(".",cursor)*"<=>"[1:min(3,width-cursor)]*repeat(".",max(0,width-cursor-3))*"]"
end

function panel_status(m::ProgressMonitor)
    phase = m.phase
    phase == "Starting" && return "starting"
    phase == "Indexing" && return "reading FASTQ, indexing minimizers"
    phase == "Filtering" && return "removing high-frequency minimizers"
    phase == "Seeding" && return "selecting initial reads from seed FASTA"
    phase == "Assembly" && return m.closed_paths > 0 ?
        "closure validated; remaining seeds running" : "extending; circular closure pending"
    phase == "Deduplication" && return "removing contained/duplicate paths"
    phase == "Scaffolding" && return isempty(m.detail) ? "scaffolding" : "scaffolding "*m.detail
    phase == "Export" && return "writing FASTA/TSV"
    phase == "LSC/IR/SSC" && return isempty(m.detail) ? "ordering LSC/IR/SSC/IR" : "ordering "*m.detail
    phase == "Annotation" && return isempty(m.detail) ? "annotating genes" : "annotating "*m.detail
    phase == "Finished" && return m.circles > 0 ? "done; $(m.circles) validated circle(s)" :
                                                  "done; no validated circle, see contigs.tsv"
    phase == "Failed" && return "failed; see the error below"
    return m.detail
end

kb_text(n::Int) = n >= 1000 ? @sprintf("%.1f kb", n/1000) : "$(n) bp"

function ascii_panel_lines(m::ProgressMonitor, stamp::UInt64, width::Int=PANEL_LOG_WIDTH)
    width = max(width, PANEL_MIN_WIDTH)
    elapsed = Float64(stamp-m.started)/1e9
    phase_elapsed = Float64(stamp-m.phase_started)/1e9
    low, high = m.genome_range
    closed = m.closed_paths > 0 || m.circles > 0
    target = high > 0 ? high : max(m.longest, 1)
    fraction = closed ? 1.0 : clamp(m.longest/target, 0.0, 1.0)
    centre = closed ? "CLOSED" : m.longest > 0 ? kb_text(m.longest) : "--"
    ring = ring_rows(fraction, centre)

    info_width = width - RING_COLS - 7          # "| " ring "   " info " |"
    count = m.total > 0 ? "$(m.done)/$(m.total) $(m.unit)" :
            (m.done > 0 ? "$(m.done) $(m.unit)" : "")
    bar_width = clamp(info_width-9-3-length(count), 6, 24)
    ram = m.current_rss === nothing ? "n/a" : @sprintf("%.2f GB", m.current_rss/1e9)
    peak = m.peak_rss > 0 ? @sprintf(" (peak %.2f GB)", m.peak_rss/1e9) : ""
    range_text = high > 0 ? "$(low ÷ 1000)-$(high ÷ 1000) kb" : "n/a"
    # Longest raw path (it still contains the closure overlap); after export, the primary result.
    contig = m.longest > 0 ? "$(m.longest) bp | target $(range_text)" : "n/a | target $(range_text)"
    circles = m.circles >= 0 ? "$(m.circles) validated" : "$(m.closed_paths) closed"
    info = (
        "phase    $(m.phase)  ($(duration_text(phase_elapsed)))",
        "stage    $(stage_bar_text(m, bar_width, stamp)) $(count)",
        "time     $(duration_text(elapsed))",
        "RAM      $(ram)$(peak)",
        "longest  $(contig)",
        "paths    $(m.paths_done)/$(m.paths_found) finished | circles $(circles)",
        "status   $(panel_status(m))",
    )
    version = " v$(VERSION_STRING) "
    title = isempty(m.label) ? " CIRCULAR PLASTOME ASSEMBLY " : " CIRCULAR PLASTOME ASSEMBLY | $(m.label) "
    room = width-5-length(version)-1
    length(title) > room && (title = first(title, room-1)*" ")
    dashes = max(1, width-2-2-length(title)-length(version)-1)
    lines = String["+--"*title*repeat("-", dashes)*version*"-+"]
    for r in 1:RING_ROWS
        push!(lines, "| "*ring[r]*"   "*fit_text(info[r], info_width)*" |")
    end
    push!(lines, "+"*repeat("-", width-2)*"+")
    return lines
end

function terminal_columns(io::IO)
    try
        return displaysize(io)[2]
    catch
        return PANEL_MAX_WIDTH+1
    end
end

# Caller owns m.lock. Moves the cursor back over the previous frame and repaints it.
function draw_panel!(m::ProgressMonitor, stamp::UInt64)
    columns = terminal_columns(m.io)
    lines = ascii_panel_lines(m, stamp, clamp(columns-1, PANEL_MIN_WIDTH, PANEL_MAX_WIDTH))
    buf = IOBuffer()
    m.panel_lines > 0 && print(buf, "\r\e[", m.panel_lines, "A")
    for line in lines
        # Never let a row wrap: wrapped rows would break the cursor-up count.
        print(buf, "\e[2K", first(line, max(1, columns-1)), '\n')
    end
    print(buf, "\e[J")
    write(m.io, take!(buf))
    flush(m.io)
    m.panel_lines = length(lines)
    return nothing
end

# Caller owns m.lock. Removes the live frame/line so other output starts on a clean row.
function erase_live_output!(m::ProgressMonitor)
    m.terminal || return nothing
    if m.panel_lines > 0
        print(m.io, "\r\e[", m.panel_lines, "A\e[J")
        m.panel_lines = 0
        flush(m.io)
    elseif m.line_open
        print(m.io, "\r\e[2K")
        m.line_open = false
        flush(m.io)
    end
    return nothing
end

# Caller owns m.lock. Sampling/terminal I/O occur at bounded frequency.
function render_progress!(m::ProgressMonitor; force::Bool=false, final::Bool=false)
    m.stopped && !final && return
    stamp = time_ns()
    if force || Float64(stamp-m.last_sample)/1e9 >= 1.0
        m.current_rss = resident_memory_bytes()
        m.current_rss === nothing || (m.peak_rss=max(m.peak_rss,m.current_rss))
        m.last_sample = stamp
    end
    if m.ascii_panel
        if m.terminal
            !force && !final && Float64(stamp-m.last_panel)/1e9 < m.panel_interval && return
            draw_panel!(m, stamp)
            m.last_panel = stamp
            m.last_render = stamp
            return
        end
        # Redirected output: one full panel per phase and at the end; one-line
        # records in between, at the configured log interval.
        if final || m.phase != m.panel_phase
            for line in ascii_panel_lines(m, stamp, PANEL_LOG_WIDTH)
                println(m.io, line)
            end
            flush(m.io)
            m.panel_phase = m.phase
            m.last_panel = stamp
            m.last_render = stamp
            return
        end
    end
    !force && Float64(stamp-m.last_render)/1e9 < m.interval && return
    line = progress_line(m,stamp)
    if m.terminal
        # Clear one line; crop only on an actual terminal to avoid line wrapping.
        columns = terminal_columns(m.io)
        line = first(line,max(1,columns-1))
        print(m.io,"\r\e[2K",line)
        final && println(m.io)
        m.line_open = !final
    else
        println(m.io,line)              # No ANSI/carriage returns in logs/pipes.
    end
    flush(m.io)
    m.last_render = stamp
end

"Forwards to the previous logger after erasing the live panel, so warnings never interleave with it."
struct ProgressAwareLogger{L<:Base.CoreLogging.AbstractLogger} <: Base.CoreLogging.AbstractLogger
    inner::L
    monitor::ProgressMonitor
end
Base.CoreLogging.min_enabled_level(l::ProgressAwareLogger) = Base.CoreLogging.min_enabled_level(l.inner)
Base.CoreLogging.shouldlog(l::ProgressAwareLogger, args...) = Base.CoreLogging.shouldlog(l.inner, args...)
Base.CoreLogging.catch_exceptions(l::ProgressAwareLogger) = Base.CoreLogging.catch_exceptions(l.inner)
function Base.CoreLogging.handle_message(l::ProgressAwareLogger, args...; kwargs...)
    lock(l.monitor.lock) do
        erase_live_output!(l.monitor)
        Base.CoreLogging.handle_message(l.inner, args...; kwargs...)
    end
end

function start_progress(; io::IO=stderr, terminal::Bool=io isa Base.TTY && get(ENV,"TERM","")!="dumb",
                        interval::Real=1.0, log_interval::Real=30.0,
                        ascii_panel::Bool=false, panel_interval::Real=1.0,
                        genome_range::Tuple{Int,Int}=(0,0), label::AbstractString="")
    stamp=time_ns()
    isfinite(panel_interval) && panel_interval >= 0.1 || error("panel_interval must be finite and >= 0.1 seconds")
    m=ProgressMonitor(ReentrantLock(),io,terminal,Float64(terminal ? interval : log_interval),
        stamp,stamp,UInt64(0),UInt64(0),"Starting",0,0,"","",nothing,0,false,false,nothing,
        Dict{String,NTuple{3,Int}}(),ascii_panel,Float64(panel_interval),UInt64(0),genome_range,
        0,0,0,0,-1,0,"",nothing,String(label))
    if terminal
        # Hide the cursor while the frame is live; stop_progress! (and atexit, for
        # abrupt exits) restores it.
        print(io, "\e[?25l")
        atexit(() -> try print(io, "\e[?25h") catch end)
        m.previous_logger = Base.CoreLogging.global_logger()
        Base.CoreLogging.global_logger(ProgressAwareLogger(m.previous_logger, m))
    end
    lock(m.lock) do
        render_progress!(m;force=true)
    end
    # The timer updates clocks/RSS while waiting for workers or I/O. Explicit
    # checkpoints also work during single-threaded CPU loops between yields.
    m.timer=Timer(min(1.0,m.interval);interval=min(1.0,m.interval)) do _
        lock(m.lock) do
            render_progress!(m)
        end
    end
    return m
end

progress_phase!(::Nothing,args...;kwargs...)=nothing
function progress_phase!(m::ProgressMonitor,phase::String;total::Int=0,unit::String="",detail::String="")
    lock(m.lock) do
        m.stopped && return
        m.phase=phase; m.phase_started=time_ns(); m.done=0; m.total=total
        m.unit=unit; m.detail=detail
        render_progress!(m;force=true)
    end
end
progress_update!(::Nothing,args...;kwargs...)=nothing
function progress_update!(m::ProgressMonitor;done=nothing,total=nothing,detail=nothing,force=false)
    lock(m.lock) do
        m.stopped && return
        done===nothing || (m.done=done)
        total===nothing || (m.total=total)
        detail===nothing || (m.detail=detail)
        render_progress!(m;force=force)
    end
end
progress_seed!(::Nothing,args...;kwargs...)=nothing
function progress_seed!(m::ProgressMonitor,id::String,n::Int,paths::Int,spawned::Int;
                        complete::Bool=false,closed::Bool=false)
    lock(m.lock) do
        m.stopped && return
        previous=get(m.seeds,id,(0,0,0))
        m.seeds[id]=(max(n,previous[1]),paths,spawned)
        m.longest=maximum(x[1] for x in values(m.seeds))
        m.paths_done=sum(x[2] for x in values(m.seeds))
        m.paths_found=sum(x[3] for x in values(m.seeds))
        m.detail="longest=$(m.longest) bp; paths=$(m.paths_done)/$(m.paths_found)"
        closed && (m.closed_paths+=1)
        complete && (m.done+=1)
        render_progress!(m;force=closed)
    end
end
progress_circles!(::Nothing,args...;kwargs...)=nothing
"Final counts at export: validated circles and the primary result (its length replaces the longest raw path)."
function progress_circles!(m::ProgressMonitor,n::Int;primary_length::Int=0,primary_circular::Bool=false)
    lock(m.lock) do
        m.stopped && return
        m.circles=n
        primary_length > 0 && (m.longest=primary_length)
        primary_circular && (m.closed_paths=max(m.closed_paths,1))
    end
end
progress_breakline!(::Nothing)=nothing
function progress_breakline!(m::ProgressMonitor)
    lock(m.lock) do
        erase_live_output!(m)
    end
end
stop_progress!(::Nothing;kwargs...)=nothing
function stop_progress!(m::ProgressMonitor;success::Bool=true)
    m.timer===nothing || close(m.timer)
    lock(m.lock) do
        m.stopped && return
        m.phase=success ? "Finished" : "Failed"
        # Do not fabricate 100% after a failed run.
        success && (m.done=1; m.total=1; m.unit="run"; m.detail="")
        render_progress!(m;force=true,final=true)
        m.stopped=true
        if m.terminal
            print(m.io, "\e[?25h"); flush(m.io)
            m.panel_lines=0             # The final frame stays on screen as the run summary.
            m.previous_logger===nothing || Base.CoreLogging.global_logger(m.previous_logger)
            m.previous_logger=nothing
        end
    end
end

Base.@kwdef mutable struct Config
    project_name::String = "Chloroplast_Assembly"
    assembly_type::String = "chl"
    forward_reads::String = ""
    reverse_reads::String = ""
    seed_input::String = ""
    output_path::String = "results_v16"
    paired::Bool = true
    read_length::Int = 150                 # Informational; actual lengths are used.
    read_pairs_to_process::Int = 0         # First N pairs, not a random sample.
    genome_range::Tuple{Int,Int} = (110000, 130000)
    kmer::Int = 31
    minimizer_window::Int = 10
    max_minimizer_frequency::Float64 = 0.01 # Fraction of actually loaded READS.
    max_postings::Int = 100000             # Hot key removed entirely, never truncated.
    batch_size::Int = 10000
    min_base_quality::Int = 0              # Phred+33; 0 disables masking.
    min_overlap::Int = 50
    min_coverage::Int = 3                  # Unique pairs (unique reads for SE).
    min_identity::Float64 = 0.98
    consensus_fraction::Float64 = 0.80
    consensus_ratio::Float64 = 2.0
    max_extension::Int = 32
    min_seed_hits::Int = 2                 # Distinct minimizers, not occurrences.
    seed_correction::Bool = true           # Local pair-aware consensus of initial read.
    bubble_policy::String = "majority"     # majority: SNV/short-indel bubbles follow the majority
                                           # allele and are reported; branch: fork (v16.7 behaviour).
    bubble_lookahead::Int = 48             # Bases followed on each allele to test reconvergence.
    bubble_max_indel::Int = 20             # Largest allele length difference treated as a bubble.
    bubble_anchor::Int = 16                # Identical bases required after the alleles.
    seeds_per_record::Int = 2
    max_initial_contigs::Int = 8
    max_branch_paths::Int = 4              # Per initial read, including primary path.
    branch_min_fraction::Float64 = 0.20
    max_iterations::Int = 20000
    max_candidates::Int = 100000           # Overflow stops a direction explicitly.
    max_concurrent_assemblies::Int = 4
    insert_size::Int = 350                 # Outer fragment length, not inter-read gap.
    insert_tolerance::Float64 = 0.30
    paired_validation::Bool = true
    min_paired_support::Int = 3            # Required at scaffold/circle junctions in PE.
    scaffold::Bool = true
    scaffold_min_overlap::Int = 100
    scaffold_max_overlap::Int = 5000
    scaffold_identity::Float64 = 0.99
    scaffold_rounds::Int = 8
    circular_min_overlap::Int = 150
    circular_max_overlap::Int = 5000
    coverage_switch_ratio::Float64 = 0.0       # Off by default: contaminants with multicopy rRNA operons have
                                               # plastid-like coverage outside the operon (L03), so coverage
                                               # cannot tell them apart. > 1 enables the experimental switch.
    coverage_trace::Bool = false               # Diagnostic: per-iteration support table (_coverage_trace.tsv).
    coverage_jump_length::Int = 2000           # Bases of sustained > ratio x support before a direction is
                                               # stopped and trimmed back (coverage_jump).
    adaptive_read_pairs::Int = 40_000_000      # No closed circle, path ends stopped for lack of reads and
                                               # the file has more pairs than Read pairs to process: assemble
                                               # again from up to this many pairs (0 = off). RAM grows with it.
    reference_dataset::String = ""             # Land-plant reference set (build_landplant_dataset.py output with
                                               # ref_sketches.bin): nearest reference reported and written.
    reference_orientation::Bool = true         # Orient the ordered output like the nearest reference.
    reference_rescue::Bool = true              # Reference seeds for stretches missing from the paths.
    reference_scaffold::Bool = true            # Reference scaffold when no circle was found (rescue output).
    plastid_panel::String = ""                 # Folder (GenBank/FASTA) or file of reference plastomes, e.g.
                                               # seedref/references/panel: graph nodes with no plastid 21-mers
                                               # (contaminants, e.g. bacteria past a shared rRNA) are removed.
    graph_topology_rescue::Bool = true        # Bounded alternative-path / IR-copy fallback.
    mate_rescue::Bool = true                   # Iterative mate-derived reseeding in graph rescue.
    mate_rescue_rounds::Int = 3
    mate_rescue_seeds_per_round::Int = 6
    reference_rescue_seeds::Int = 6             # Reference seeds per rescue round (one per missing stretch), extra to the above.
    mate_rescue_max_bases::Int = 50000          # Per new path, not a target genome size.
    mate_rescue_branch_paths::Int = 2
    mate_rescue_max_candidates::Int = 4000     # Per graph-side; saturation is reported.
    rescue_seed_fasta::String = ""             # Optional additional seeds; auto-detect PROJECT_seeds_all.fasta.
    island_recovery::Bool = true               # Graph rescue assembles sequence between two obstacles from
                                               # mates anchored at open ends (extra paths for the graph).
    graph_read_bridges::Bool = true          # Read-spanned, graph-unique flanks; graph rescue only.
    inferred_joins::Bool = true                # Graph rescue may join dead ends in the same microsatellite
                                               # motif with Ns when nothing connects them (inferred_join).
    graph_rescue::Bool = true                  # No gap-free circle: search circular walks in the graph of
                                               # all paths, node copy numbers from coverage (graph_path).
    mate_jumps::Int = 2                        # Per path: jumps over an unreadable stretch (e.g. long poly-G)
                                               # assembled from anchored mates; the gap stays as Ns. 0 = off.
    phase_window::Int = 31                     # Extension reads must match the last this-many contig
                                               # bases exactly: keeps nearby variants on one haplotype.
    repeat_anchor::Int = 20                    # A contig end inside a tandem repeat accepts only reads
                                               # crossing the repeat plus this many unique bases.
    unknown_gap_ns::Int = 10                   # Repeat seen from both ends but no length estimate:
                                               # close with this many Ns (NCBI unknown-length gap); 0 = off.
    repeat_gap_closure::Bool = true            # Close a linear genome whose ends stop in the same
                                               # short tandem repeat; length from read pairs, flagged.
    circular_repeat_max_overlap::Int = 30000   # Longer end overlaps, tested only where the
                                               # path's first 31-mer recurs (seam inside an IR).
    circular_identity::Float64 = 0.995
    circular_min_reads::Int = 3
    circular_check_every::Int = 200        # Newly added bases between checks.
    junction_anchor::Int = 20
    min_contig_length::Int = 300
    normalize_plastome::Bool = true       # Additional rotated FASTA for circular chl results.
    annotation_dir::String = ""           # Directory with PlastomeAnnotator.jl:
                                          # genes of the ordered plastomes -> _plastome_LSC_IR_SSC_IR.gff3.
                                          # Empty = no annotation.
    annotation_dataset::String = ""       # Directory with refindex_all and models_final (annotator dataset);
                                          # empty = ~/Assembler/refseq_plastid/dataset.
    ir_min_length::Int = 1000
    ir_min_identity::Float64 = 0.995
    ir_min_sc_length::Int = 100
    debug_mode::Bool = false
    show_progress::Bool = true
    progress_interval::Float64 = 1.0      # Terminal refresh, seconds.
    progress_log_interval::Float64 = 30.0 # Redirected output refresh, seconds.
    ascii_panel::Bool = true              # Plain ASCII panel also works through a pipe.
    ascii_panel_interval::Float64 = 1.0
end

const CONFIG_FIELDS = Dict{String,Symbol}(
    "Graph topology rescue" => :graph_topology_rescue,
    "Mate rescue" => :mate_rescue, "Mate rescue rounds" => :mate_rescue_rounds,
    "Mate rescue seeds per round" => :mate_rescue_seeds_per_round,
    "Mate rescue max bases" => :mate_rescue_max_bases,
    "Mate rescue branch paths" => :mate_rescue_branch_paths,
    "Mate rescue max candidates" => :mate_rescue_max_candidates,
    "Rescue seed FASTA" => :rescue_seed_fasta,
    "Project name" => :project_name, "Assembly type" => :assembly_type,
    "Forward reads" => :forward_reads, "Reverse reads" => :reverse_reads,
    "Seed Input" => :seed_input, "Output path" => :output_path,
    "Read Length" => :read_length, "Read pairs to process" => :read_pairs_to_process,
    "K-mer" => :kmer, "Minimizer window" => :minimizer_window,
    "Max minimizer frequency" => :max_minimizer_frequency,
    "Max postings per minimizer" => :max_postings, "Batch size" => :batch_size,
    "Min base quality" => :min_base_quality, "Min overlap" => :min_overlap,
    "Min coverage" => :min_coverage, "Min identity" => :min_identity,
    "Consensus fraction" => :consensus_fraction, "Consensus ratio" => :consensus_ratio,
    "Max extension" => :max_extension, "Min seed hits" => :min_seed_hits,
    "Seed correction" => :seed_correction,
    "Bubble policy" => :bubble_policy, "Bubble lookahead" => :bubble_lookahead,
    "Bubble max indel" => :bubble_max_indel, "Bubble anchor" => :bubble_anchor,
    "Seeds per record" => :seeds_per_record, "Max initial contigs" => :max_initial_contigs,
    "Max branch paths" => :max_branch_paths, "Branch min fraction" => :branch_min_fraction,
    "Max iterations" => :max_iterations, "Max candidates" => :max_candidates,
    "Max concurrent assemblies" => :max_concurrent_assemblies,
    "Insert size" => :insert_size, "Insert tolerance" => :insert_tolerance,
    "Paired validation" => :paired_validation, "Min paired support" => :min_paired_support,
    "Scaffolding" => :scaffold, "Scaffold min overlap" => :scaffold_min_overlap,
    "Scaffold max overlap" => :scaffold_max_overlap, "Scaffold identity" => :scaffold_identity,
    "Scaffold rounds" => :scaffold_rounds, "Circular min overlap" => :circular_min_overlap,
    "Circular max overlap" => :circular_max_overlap,
    "Circular repeat max overlap" => :circular_repeat_max_overlap,
    "Repeat gap closure" => :repeat_gap_closure, "Unknown gap Ns" => :unknown_gap_ns,
    "Repeat anchor" => :repeat_anchor, "Phase window" => :phase_window, "Mate jumps" => :mate_jumps, "Graph rescue" => :graph_rescue, "Graph read bridges" => :graph_read_bridges, "Inferred joins" => :inferred_joins, "Island recovery" => :island_recovery, "Plastid panel" => :plastid_panel, "Reference dataset" => :reference_dataset, "Reference orientation" => :reference_orientation, "Reference rescue" => :reference_rescue, "Reference scaffold" => :reference_scaffold, "Reference rescue seeds" => :reference_rescue_seeds, "Adaptive read pairs" => :adaptive_read_pairs,
    "Coverage switch ratio" => :coverage_switch_ratio, "Coverage jump length" => :coverage_jump_length,
    "Coverage trace" => :coverage_trace, "Circular identity" => :circular_identity,
    "Circular min reads" => :circular_min_reads,
    "Circular check every" => :circular_check_every, "Junction anchor" => :junction_anchor,
    "Normalize plastome" => :normalize_plastome, "IR min length" => :ir_min_length,
    "Annotation dir" => :annotation_dir, "Annotation dataset" => :annotation_dataset,
    "IR min identity" => :ir_min_identity, "IR min SC length" => :ir_min_sc_length,
    "Min contig length" => :min_contig_length, "Debug mode" => :debug_mode,
    "Show progress" => :show_progress, "Progress interval" => :progress_interval,
    "Progress log interval" => :progress_log_interval, "ASCII genome panel" => :ascii_panel,
    "ASCII panel interval" => :ascii_panel_interval)

function validate_config(c::Config; files::Bool = true)
    c.mate_rescue_rounds >= 0 || error("Mate rescue rounds must be >= 0")
    min(c.mate_rescue_seeds_per_round,c.mate_rescue_branch_paths,c.mate_rescue_max_candidates)>0 || error("Mate rescue budgets must be positive")
    c.mate_rescue_max_bases>=300 || error("Mate rescue max bases must be >= 300")
    occursin(r"^[A-Za-z0-9_.-]+$", c.project_name) || error("Invalid Project name; use letters, digits, _, - or .")
    c.project_name in (".", "..") && error("Invalid Project name")
    c.assembly_type in ("chl", "mito", "organelle") || error("Assembly type: chl, mito or organelle")
    1 <= c.kmer <= 32 || error("K-mer must be in 1:32")
    for field in (:minimizer_window, :max_postings, :batch_size, :min_overlap,
                  :min_coverage, :max_extension, :min_seed_hits, :seeds_per_record,
                  :max_initial_contigs, :max_branch_paths, :max_iterations,
                  :max_candidates, :max_concurrent_assemblies, :insert_size,
                  :scaffold_min_overlap, :scaffold_max_overlap, :scaffold_rounds,
                  :circular_min_overlap, :circular_max_overlap, :circular_min_reads,
                  :circular_check_every, :junction_anchor, :min_contig_length, :read_length)
        getfield(c, field) > 0 || error("$field must be positive")
    end
    c.read_pairs_to_process >= 0 || error("Read pairs to process cannot be negative")
    for field in (:progress_interval,:progress_log_interval,:ascii_panel_interval)
        isfinite(getfield(c,field)) && getfield(c,field)>=0.1 || error("$field must be finite and >= 0.1 seconds")
    end
    0 <= c.min_base_quality <= 93 || error("Min base quality must be in 0:93")
    c.genome_range[1] > 0 && c.genome_range[2] >= c.genome_range[1] || error("Invalid Genome Range")
    for field in (:max_minimizer_frequency, :min_identity, :consensus_fraction,
                  :branch_min_fraction, :scaffold_identity, :circular_identity)
        0 < getfield(c, field) <= 1 || error("$field must be in (0, 1]")
    end
    c.ir_min_length >= 31 && c.ir_min_sc_length > 0 || error("Invalid IR region length thresholds")
    0 < c.ir_min_identity <= 1 || error("IR min identity must be in (0, 1]")
    c.consensus_fraction > 0.5 || error("Consensus fraction must be > 0.5")
    c.bubble_policy in ("majority", "branch") || error("Bubble policy must be majority or branch")
    c.repeat_anchor >= 0 || error("Repeat anchor cannot be negative")
    c.phase_window >= 0 || error("Phase window cannot be negative")
    c.mate_jumps >= 0 || error("Mate jumps cannot be negative")
    c.adaptive_read_pairs >= 0 || error("Adaptive read pairs cannot be negative")
    (c.coverage_switch_ratio == 0 || c.coverage_switch_ratio > 1) || error("Coverage switch ratio must be 0 or > 1")
    c.coverage_jump_length > 0 || error("Coverage jump length must be positive")
    c.bubble_anchor >= 8 && c.bubble_max_indel >= 0 || error("Bubble anchor must be >= 8 and Bubble max indel >= 0")
    c.bubble_lookahead >= c.bubble_anchor + c.bubble_max_indel + 1 ||
        error("Bubble lookahead must be >= Bubble anchor + Bubble max indel + 1")
    c.consensus_ratio > 1 || error("Consensus ratio must be > 1")
    0 <= c.insert_tolerance < 1 || error("Insert tolerance must be in [0, 1)")
    c.min_paired_support > 0 || error("Min paired support must be positive")
    c.scaffold_min_overlap <= c.scaffold_max_overlap || error("Invalid scaffold overlap range")
    c.circular_min_overlap <= c.circular_max_overlap || error("Invalid circular overlap range")
    c.circular_repeat_max_overlap >= c.circular_max_overlap || error("Circular repeat max overlap must be >= Circular max overlap")
    isempty(c.output_path) && error("Output path cannot be empty")
    if files
        for (label, path) in (("Forward reads", c.forward_reads), ("Seed Input", c.seed_input))
            isfile(path) || error("$label: file not found: $path")
        end
        c.paired && !isfile(c.reverse_reads) && error("Reverse reads: file not found: $(c.reverse_reads)")
        c.paired && realpath(c.forward_reads) == realpath(c.reverse_reads) && error("R1 and R2 must be different files")
        isempty(c.annotation_dir) || isfile(joinpath(c.annotation_dir, "PlastomeAnnotator.jl")) ||
            error("Annotation dir: PlastomeAnnotator.jl not found in $(c.annotation_dir)")
    end
    return c
end

# Keys read by the seed/reference stage (SeedRef.jl, BatchAssemble.jl); the
# assembler accepts them in a shared config file and ignores them silently.
const SEED_STAGE_KEYS = Set(["SeedIR reads", "K-mer size", "Min k-mer hits for seed",
    "Input folder", "Marker folder", "Seed markers", "Seed reads per marker",
    "Reference mode", "Reference GenBank", "Genome range margin",
    "BLAST database", "BLAST Entrez query", "BLAST hitlist size", "BLAST timeout",
    "Samples", "Skip completed", "Common seed FASTA", "Marker panel", "Refine seeds with reference"])

"Apply one `key = value` entry to `c`. Returns false for keys the assembler does not own."
function apply_config_entry!(c::Config, key::AbstractString, value::AbstractString)
    if key == "Single/Paired"
        uppercase(value) in ("PE", "SE") || error("Single/Paired must be PE or SE")
        c.paired = uppercase(value) == "PE"
    elseif key == "Genome Range"
        m = match(r"^(\d+)\s*-\s*(\d+)$", value)
        m === nothing && error("Genome Range must be MIN-MAX")
        c.genome_range = (parse(Int, m.captures[1]), parse(Int, m.captures[2]))
    elseif haskey(CONFIG_FIELDS, key)
        field = CONFIG_FIELDS[key]
        T = fieldtype(Config, field)
        parsed = T === String ? String(value) : T === Bool ? parse(Bool, lowercase(value)) : parse(T, value)
        setfield!(c, field, parsed)
    else
        return false
    end
    return true
end

"Resolve relative input/output paths against `base` (the config file's directory)."
function resolve_config_paths!(c::Config, base::AbstractString)
    for field in (:forward_reads, :reverse_reads, :seed_input, :output_path, :plastid_panel, :rescue_seed_fasta, :reference_dataset,
                  :annotation_dir, :annotation_dataset)
        value = getfield(c, field)
        isempty(value) || setfield!(c, field, normpath(joinpath(base, expanduser(value))))
    end
    return c
end

"Write a complete, standalone config (absolute paths) that reproduces `c`."
function write_config(io::IO, c::Config)
    println(io, "# AssemblerComplete ", VERSION_STRING, " - complete configuration used for this run")
    println(io, "Single/Paired = ", c.paired ? "PE" : "SE")
    println(io, "Genome Range = ", c.genome_range[1], "-", c.genome_range[2])
    for key in sort(collect(keys(CONFIG_FIELDS)))
        field = CONFIG_FIELDS[key]
        field === :paired && continue
        value = getfield(c, field)
        value isa AbstractString && isempty(value) && continue
        println(io, key, " = ", value)
    end
end

function parse_config(path::AbstractString)
    c = Config()
    base = dirname(abspath(path))
    seen = Set{String}()
    for (lineno, raw) in enumerate(eachline(path))
        line = strip(first(split(raw, '#'; limit=2)))
        isempty(line) && continue
        occursin('=', line) || error("Config line $lineno: expected key = value")
        key, value = strip.(split(line, '='; limit=2))
        key in seen && error("Duplicate configuration key: $key")
        push!(seen, key)
        key in SEED_STAGE_KEYS && continue
        apply_config_entry!(c, key, value) || error("Unknown configuration key '$key' at line $lineno")
    end
    resolve_config_paths!(c, base)
    validate_config(c)
    c.min_overlap < c.kmer + c.minimizer_window - 1 && @warn "Min overlap is shorter than a complete minimizer window; some overlaps cannot be retrieved."
    c.min_identity < 0.95 && @warn "Min identity < 0.95: permissive for short-read assembly."
    return c
end

const IUPAC_AMBIGUOUS = let a = falses(256)
    for b in codeunits("NRYKMSWBDHV")
        a[Int(b)+1] = true
    end
    a
end

# One pass over the bytes; ASCII case folding replaces uppercase(), which
# allocated a second String for every read.
function normalized_dna(s::AbstractString)
    src = codeunits(s)
    out = Vector{UInt8}(undef, length(src))
    @inbounds for i in eachindex(src)
        b = src[i]
        UInt8('a') <= b <= UInt8('z') && (b -= 0x20)
        if ENCODE[Int(b)+1] == 255
            IUPAC_AMBIGUOUS[Int(b)+1] || error("Invalid DNA byte $(repr(Char(b)))")
            b = UInt8('N')
        end
        out[i] = b
    end
    return out
end

function reverse_complement!(out::Vector{UInt8}, seq::AbstractVector{UInt8})
    n = length(seq)
    resize!(out, n)
    @inbounds for (j, i) in enumerate(eachindex(seq))
        out[n - j + 1] = COMPLEMENT[Int(seq[i]) + 1]
    end
    return out
end
reverse_complement(seq::AbstractVector{UInt8}) = reverse_complement!(Vector{UInt8}(undef, length(seq)), seq)
reverse_complement(s::AbstractString) = String(reverse_complement(normalized_dna(s)))

# Validity is explicit: UInt64(0xffffffffffffffff) is a valid 32-mer (T^32).
function rolling_kmers(seq::AbstractVector{UInt8}, k::Int)
    n = max(0, length(seq) - k + 1)
    keys = zeros(UInt64, n)
    valid = falses(n)
    mask = k == 32 ? typemax(UInt64) : (UInt64(1) << (2k)) - UInt64(1)
    f::UInt64 = 0
    r::UInt64 = 0
    run = 0
    @inbounds for i in eachindex(seq)
        b = ENCODE[Int(seq[i]) + 1]
        if b == 255
            f = 0; r = 0; run = 0
        else
            f = ((f << 2) | UInt64(b)) & mask
            r = (r >> 2) | (UInt64(b ⊻ 0x03) << (2(k-1)))
            run += 1
            if run >= k
                p = i - k + 1
                keys[p] = min(f, r)
                valid[p] = true
            end
        end
    end
    return keys, valid
end

# Fixed SplitMix64 permutation: uniform ordering without session-dependent hash().
@inline function mix64(x::UInt64)
    x ⊻= x >> 30; x *= 0xbf58476d1ce4e5b9
    x ⊻= x >> 27; x *= 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

"Scratch owned by one logical task, never indexed by a migrating task's thread ID."
struct MinimizerWorkspace
    keys::Vector{UInt64}
    ranks::Vector{UInt64}
    positions::Vector{Int}
    selected::Vector{UInt64}
end
MinimizerWorkspace() = MinimizerWorkspace(UInt64[], UInt64[], Int[], UInt64[])

"Returned storage is borrowed until the next call with the same workspace."
function minimizer_keys!(work::MinimizerWorkspace, seq::AbstractVector{UInt8}, k::Int, w::Int)
    1 <= k <= 32 || throw(ArgumentError("k must be in 1:32"))
    w >= 1 || throw(ArgumentError("minimizer window must be positive"))
    selected = work.selected
    empty!(selected)
    n = length(seq)
    n-k+1 < w && return selected
    if length(work.keys) < n
        resize!(work.keys, n); resize!(work.ranks, n); resize!(work.positions, n)
    end
    mask = k == 32 ? typemax(UInt64) : (UInt64(1) << (2k))-UInt64(1)
    f = UInt64(0); r = UInt64(0); run = 0; head = 1; tail = 0
    @inbounds for i in eachindex(seq)
        b = ENCODE[Int(seq[i])+1]
        if b == 255
            f = 0; r = 0; run = 0; head = 1; tail = 0
            continue
        end
        f = ((f << 2) | UInt64(b)) & mask
        r = (r >> 2) | (UInt64(b ⊻ 0x03) << (2(k-1)))
        run += 1
        run < k && continue
        key = min(f,r); rank = mix64(key)
        while head <= tail && work.positions[head] <= i-w
            head += 1
        end
        while head <= tail && work.ranks[tail] >= rank
            tail -= 1
        end
        tail += 1
        work.keys[tail] = key; work.ranks[tail] = rank; work.positions[tail] = i
        if run >= k+w-1
            winner = work.keys[head]
            (isempty(selected) || selected[end] != winner) && push!(selected, winner)
        end
    end
    sort!(selected)
    unique!(selected) # Once per read; duplicate keys in separated windows also collapse.
    return selected
end
minimizer_keys(seq::AbstractVector{UInt8}, k::Int, w::Int) =
    minimizer_keys!(MinimizerWorkspace(), seq, k, w)

"Sorted unique canonical k-mers of the valid windows of `seq`, written into the reused buffer `out` (in-place QuickSort)."
function all_kmer_keys!(out::Vector{UInt64}, seq::AbstractVector{UInt8}, k::Int)
    empty!(out)
    mask = k == 32 ? typemax(UInt64) : (UInt64(1) << (2k)) - UInt64(1)
    f::UInt64 = 0; r::UInt64 = 0; run = 0
    @inbounds for i in eachindex(seq)
        b = ENCODE[Int(seq[i]) + 1]
        if b == 255
            f = 0; r = 0; run = 0
        else
            f = ((f << 2) | UInt64(b)) & mask
            r = (r >> 2) | (UInt64(b ⊻ 0x03) << (2(k-1)))
            run += 1
            run >= k && push!(out, min(f, r))
        end
    end
    sort!(out; alg=QuickSort)
    w = 0
    @inbounds for v in out
        if w == 0 || out[w] != v
            w += 1; out[w] = v
        end
    end
    return resize!(out, w)
end
all_kmer_keys(seq::AbstractVector{UInt8}, k::Int) = all_kmer_keys!(UInt64[], seq, k)

struct CompactReadStore
    data::Vector{UInt64}
    offsets::Vector{UInt32}
    lengths::Vector{UInt16}
    ambiguous::Dict{UInt32,Vector{UInt16}}
end
CompactReadStore() = CompactReadStore(UInt64[], UInt32[], UInt16[], Dict{UInt32,Vector{UInt16}}())

function push_read!(s::CompactReadStore, seq::AbstractVector{UInt8})
    n = length(seq)
    0 < n <= typemax(UInt16) || error("Read length outside supported 1:65535 range")
    length(s.lengths) < typemax(UInt32) || error("Too many reads")
    start = length(s.data) + 1
    words = cld(n, 32)
    start + words - 1 <= typemax(UInt32) || error("Read store offset limit reached")
    resize!(s.data, length(s.data)+words)
    @inbounds for p in start:start+words-1
        s.data[p] = 0
    end
    unknown = nothing
    @inbounds for i in eachindex(seq)
        b = ENCODE[Int(seq[i])+1]
        if b == 255
            unknown === nothing && (unknown = UInt16[])
            push!(unknown, UInt16(i)); b = 0
        end
        s.data[start + ((i-1) >> 5)] |= UInt64(b) << (2((i-1) & 31))
    end
    push!(s.offsets, UInt32(start)); push!(s.lengths, UInt16(n))
    id = UInt32(length(s.lengths))
    unknown === nothing || (s.ambiguous[id] = unknown)
    return id
end

"Decode read `id` into `out` (resized); lets hot loops reuse one buffer."
function get_read!(out::Vector{UInt8}, s::CompactReadStore, id::Integer)
    1 <= id <= length(s.lengths) || throw(BoundsError(s.lengths, id))
    n = Int(s.lengths[id]); start = Int(s.offsets[id])
    resize!(out, n)
    @inbounds for i in 1:n
        b = (s.data[start + ((i-1) >> 5)] >> (2((i-1) & 31))) & 0x03
        out[i] = DNA[Int(b)+1]
    end
    if haskey(s.ambiguous, UInt32(id))
        for p in s.ambiguous[UInt32(id)]
            out[p] = UInt8('N')
        end
    end
    return out
end
get_read(s::CompactReadStore, id::Integer) = get_read!(UInt8[], s, id)

struct PairedReadStore
    r1::CompactReadStore
    r2::CompactReadStore
    paired::Bool
end
PairedReadStore(paired::Bool) = PairedReadStore(CompactReadStore(), CompactReadStore(), paired)
num_pairs(s::PairedReadStore) = length(s.r1.lengths)
@inline pair_id(id::UInt32) = id >> 1
@inline read_id(pair::Integer, first::Bool) = (UInt32(pair) << 1) | UInt32(first)
get_read(s::PairedReadStore, id::UInt32) = get_read(isodd(id) ? s.r1 : s.r2, pair_id(id))
get_read!(out::Vector{UInt8}, s::PairedReadStore, id::UInt32) =
    get_read!(out, isodd(id) ? s.r1 : s.r2, pair_id(id))

function open_fastq(path::AbstractString)
    raw = open(path, "r")
    if endswith(lowercase(path), ".gz")
        try
            return GzipDecompressorStream(raw; bufsize=1 << 20)   # default is 16 KiB
        catch
            close(raw); rethrow()
        end
    end
    return raw
end

"Append reads exactly as repeated push_read! would, packing them on all threads."
function append_reads!(s::CompactReadStore, seqs::Vector{Vector{UInt8}})
    isempty(seqs) && return s
    n0 = length(s.lengths)
    n0 + length(seqs) < typemax(UInt32) || error("Too many reads")
    starts = Vector{Int}(undef, length(seqs))
    pos = length(s.data) + 1
    for (i, seq) in enumerate(seqs)
        0 < length(seq) <= typemax(UInt16) || error("Read length outside supported 1:65535 range")
        starts[i] = pos; pos += cld(length(seq), 32)
    end
    pos - 1 <= typemax(UInt32) || error("Read store offset limit reached")
    resize!(s.data, pos - 1)
    unknown = Vector{Union{Nothing,Vector{UInt16}}}(nothing, length(seqs))
    data = s.data
    chunks = min(Threads.nthreads(), length(seqs))
    Threads.@threads :dynamic for chunk in 1:chunks
        for i in fld((chunk-1)*length(seqs),chunks)+1:fld(chunk*length(seqs),chunks)
            seq = seqs[i]; start = starts[i]
            @inbounds for p in start:start+cld(length(seq), 32)-1     # disjoint word ranges
                data[p] = 0
            end
            u = nothing
            @inbounds for j in eachindex(seq)
                b = ENCODE[Int(seq[j])+1]
                if b == 255
                    u === nothing && (u = UInt16[])
                    push!(u, UInt16(j)); b = 0
                end
                data[start + ((j-1) >> 5)] |= UInt64(b) << (2((j-1) & 31))
            end
            unknown[i] = u
        end
    end
    for (i, seq) in enumerate(seqs)
        push!(s.offsets, UInt32(starts[i])); push!(s.lengths, UInt16(length(seq)))
        unknown[i] === nothing || (s.ambiguous[UInt32(n0+i)] = unknown[i])
    end
    return s
end

@inline function fastq_line(io::IO, path::AbstractString, record::Int, part::Int)
    eof(io) && error("Truncated FASTQ: $path, record $record, line $part")
    return readline(io)
end

"First whitespace-delimited identifier, with a terminal /1 or /2 removed."
function fastq_name(header::String, path::AbstractString, record::Int)
    first = nextind(header, firstindex(header))
    last = lastindex(header)
    while first <= last && isspace(header[first])
        first = nextind(header,first)
    end
    first <= last || error("Empty FASTQ identifier: $path, record $record")
    stop = first
    while stop <= last && !isspace(header[stop])
        stop = nextind(header,stop)
    end
    token = SubString(header,first,prevind(header,stop))
    if endswith(token,"/1") || endswith(token,"/2")
        return String(SubString(token,firstindex(token),prevind(token,lastindex(token),2)))
    end
    return String(token)
end

const RawRecord = NTuple{4,String}

"Four FASTQ lines without interpretation (done by the per-file reader task)."
function read_raw_record(io::IO, path::AbstractString, record::Int)
    eof(io) && return nothing
    return (fastq_line(io,path,record,1), fastq_line(io,path,record,2),
            fastq_line(io,path,record,3), fastq_line(io,path,record,4))
end

function fastq_record(io::IO, path::AbstractString, record::Int, minq::Int)
    raw = read_raw_record(io, path, record)
    raw === nothing && return nothing
    return parse_fastq_record(raw, path, record, minq)
end

"Validation, base normalisation and quality masking of one record (thread-safe)."
function parse_fastq_record(raw::RawRecord, path::AbstractString, record::Int, minq::Int)
    header, text, plus, quality = raw
    startswith(header, "@") && ncodeunits(header) > 1 || error("Invalid FASTQ header: $path, record $record")
    startswith(plus, "+") || error("Missing + line: $path, record $record")
    ncodeunits(text) == ncodeunits(quality) || error("Sequence/quality length mismatch: $path, record $record")
    seq = normalized_dna(text)
    isempty(seq) && error("Empty FASTQ sequence: $path, record $record")
    q = codeunits(quality)
    for i in eachindex(q)
        33 <= q[i] <= 126 || error("Invalid Phred+33 quality: $path, record $record")
        Int(q[i]) - 33 < minq && (seq[i] = UInt8('N'))
    end
    return fastq_name(header,path,record), seq
end

# ---------------------------------------------------------------------------
# Sharded minimizer index.
#
# Build: every minimizer occurrence is appended as one packed UInt128
# (index_hash(key) << 32 | read id) to a buffer owned by (minimizer chunk, shard).
# Minimizer computation and "insertion" therefore run in the same threaded loop,
# with no hashing, no locks and no per-key heap objects.
#
# Freeze (Filtering phase, threaded over shards): concatenate, sort, run-length
# count, drop keys with count > cutoff, and store CSR arrays per shard.
#
# index_hash is bijective (odd-constant add + SplitMix64 finalizer), so equal
# hashes mean equal keys: lookups are exact, never probabilistic. The offset
# decorrelates it from mix64(key), the minimizer rank, whose selected values are
# strongly biased towards small numbers and would unbalance the shards.
#
# Result equivalence with v16.4: the retained key set and per-key read sets are
# the same (count > cutoff removed; cutoff <= Max postings). Ids within a key are
# now strictly ascending (previously R1 before R2 within a pair); every consumer
# sorts or scores ids, so no output depends on that order.
# ---------------------------------------------------------------------------
const INDEX_HASH_OFFSET = 0x9e3779b97f4a7c15
@inline index_hash(key::UInt64) = mix64(key + INDEX_HASH_OFFSET)
@inline pack_posting(h::UInt64, id::UInt32) = (UInt128(h) << 32) | UInt128(id)
@inline packed_hash(x::UInt128) = UInt64(x >> 32)
@inline packed_id(x::UInt128) = UInt32(x & 0xffffffff)

struct PostingShard
    hashes::Vector{UInt64}          # Sorted index_hash of retained keys.
    offsets::Vector{UInt32}         # CSR, zero-based: ids[offsets[j]+1:offsets[j+1]].
    ids::Vector{UInt32}             # Ascending read ids per key.
    directory::Vector{UInt32}       # Bucket t spans hashes[directory[t+1]+1:directory[t+2]].
    dbits::Int
end

struct ShardedPostings
    shards::Vector{PostingShard}
    sbits::Int
    nkeys::Int
end

@inline shard_index(h::UInt64, sbits::Int) = Int(h >> (64 - sbits)) + 1
@inline shard_bucket(h::UInt64, sbits::Int, dbits::Int) = Int((h << sbits) >> (64 - dbits))

"Return (shard, j) with j > 0 if the key is retained, (shard, 0) otherwise."
@inline function posting_lookup(p::ShardedPostings, key::UInt64)
    h = index_hash(key)
    s = @inbounds p.shards[shard_index(h, p.sbits)]
    t = shard_bucket(h, p.sbits, s.dbits)
    @inbounds for j in Int(s.directory[t+1])+1:Int(s.directory[t+2])
        x = s.hashes[j]
        x == h && return s, j
        x > h && break
    end
    return s, 0
end
@inline posting_ids(s::PostingShard, j::Int) =
    @inbounds @view s.ids[Int(s.offsets[j])+1:Int(s.offsets[j+1])]

Base.length(p::ShardedPostings) = p.nkeys
Base.isempty(p::ShardedPostings) = p.nkeys == 0
Base.haskey(p::ShardedPostings, key::UInt64) = posting_lookup(p, key)[2] != 0
function Base.getindex(p::ShardedPostings, key::UInt64)
    s, j = posting_lookup(p, key)
    j == 0 && throw(KeyError(key))
    return posting_ids(s, j)
end
function Base.get(p::ShardedPostings, key::UInt64, default)
    s, j = posting_lookup(p, key)
    return j == 0 ? default : posting_ids(s, j)
end

"Shard bits: about four shards per thread (load balance during the freeze), 16-1024 shards."
index_shard_bits(threads::Int=Threads.nthreads()) = clamp(ceil(Int, log2(4*max(1, threads))), 4, 10)

"Sorted packed occurrences of one shard -> CSR shard. Returns (shard, distinct keys, dropped keys)."
function build_shard(sorted::Vector{UInt128}, cutoff::Int, sbits::Int)
    n = length(sorted)
    nkeys = 0; nids = 0; distinct = 0; dropped = 0
    i = 1
    while i <= n
        h = packed_hash(sorted[i]); j = i
        while j < n && packed_hash(sorted[j+1]) == h
            j += 1
        end
        distinct += 1
        if j-i+1 <= cutoff
            nkeys += 1; nids += j-i+1
        else
            dropped += 1
        end
        i = j+1
    end
    nids < typemax(UInt32) || error("Index shard exceeds 2^32 postings; increase shard count")
    hashes = Vector{UInt64}(undef, nkeys)
    offsets = Vector{UInt32}(undef, nkeys+1)
    ids = Vector{UInt32}(undef, nids)
    offsets[1] = 0
    k = 0; cursor = 0; i = 1
    while i <= n
        h = packed_hash(sorted[i]); j = i
        while j < n && packed_hash(sorted[j+1]) == h
            j += 1
        end
        if j-i+1 <= cutoff
            k += 1
            hashes[k] = h
            for r in i:j
                cursor += 1
                ids[cursor] = packed_id(sorted[r])
            end
            offsets[k+1] = UInt32(cursor)
        end
        i = j+1
    end
    # About four keys per directory bucket: one short scan within a cache line or two.
    dbits = clamp(ceil(Int, log2(max(nkeys, 1)/4 + 1)), 1, 30)
    directory = zeros(UInt32, (1 << dbits) + 1)
    for h in hashes
        directory[shard_bucket(h, sbits, dbits)+2] += 1
    end
    for t in 2:length(directory)
        directory[t] += directory[t-1]
    end
    return PostingShard(hashes, offsets, ids, directory, dbits), distinct, dropped
end

"Threaded over shards; frees each chunk buffer as soon as it has been copied."
function freeze_sharded!(buffers::Vector{Vector{Vector{UInt128}}}, sbits::Int, cutoff::Int; monitor=nothing)
    nshards = 1 << sbits
    shards = Vector{PostingShard}(undef, nshards)
    distinct = zeros(Int, nshards); dropped = zeros(Int, nshards)
    finished = Threads.Atomic{Int}(0)
    Threads.@threads :dynamic for s in 1:nshards
        total = 0
        for chunk in buffers
            total += length(chunk[s])
        end
        occurrences = Vector{UInt128}(undef, total)
        pos = 1
        for chunk in buffers
            part = chunk[s]
            copyto!(occurrences, pos, part, 1, length(part))
            pos += length(part)
            chunk[s] = UInt128[]            # Distinct slot per shard: no data race.
        end
        sort!(occurrences)                  # (hash, id): ids ascending within a key.
        shards[s], distinct[s], dropped[s] = build_shard(occurrences, cutoff, sbits)
        progress_update!(monitor;done=Threads.atomic_add!(finished, 1)+1)
    end
    return ShardedPostings(shards, sbits, sum(x -> length(x.hashes), shards)), sum(distinct), sum(dropped)
end

struct ReadIndex{P}
    store::PairedReadStore
    postings::P                     # ShardedPostings after build_index; any Dict-like map works.
    max_read_length::Int
    cutoff::Int
    dropped_keys::Int
    short_reads::Int
end

"Wall seconds of indexing sub-phases of the last build_index call (written to _run.txt)."
const PHASE_TIMES = Dict{String,Float64}()


"""
Line splitter over large decompressed blocks: one bulk read per 4 MiB instead of one
readline call per line on the decompression stream. Lines match readline(io): the line
break ("\\n" or "\\r\\n") is removed and a last line without a break is returned.
"""
mutable struct ByteLines
    io::IO
    buf::Vector{UInt8}
    chunk::Vector{UInt8}
    start::Int
    stop::Int
    done::Bool
end
ByteLines(io::IO) = ByteLines(io, Vector{UInt8}(undef, 1 << 23), Vector{UInt8}(undef, 1 << 22), 1, 0, false)

"Next line as a range of r.buf (valid until the next call), or nothing at EOF."
function next_line_range!(r::ByteLines)
    while true
        if r.start <= r.stop
            i = GC.@preserve r begin
                p = ccall(:memchr, Ptr{UInt8}, (Ptr{UInt8}, Cint, Csize_t),
                          pointer(r.buf, r.start), 0x0a, r.stop - r.start + 1)
                p == C_NULL ? 0 : Int(p - pointer(r.buf)) + 1
            end
            if i != 0 || r.done
                e = i != 0 ? i - 1 : r.stop
                i != 0 && e >= r.start && r.buf[e] == 0x0d && (e -= 1)   # only "\r\n", as readline
                a = r.start
                r.start = i != 0 ? i + 1 : r.stop + 1
                return (a, e)
            end
        elseif r.done
            return nothing
        end
        # Keep the unfinished line, then append the next decompressed block.
        rest = r.stop - r.start + 1
        rest > 0 && r.start > 1 && copyto!(r.buf, 1, r.buf, r.start, rest)
        r.start = 1; r.stop = max(rest, 0)
        n = readbytes!(r.io, r.chunk, length(r.chunk))
        if n == 0
            r.done = true
        else
            r.stop + n > length(r.buf) && resize!(r.buf, 2*(r.stop + n))
            copyto!(r.buf, r.stop + 1, r.chunk, 1, n)
            r.stop += n
        end
    end
end

function next_line!(r::ByteLines)
    x = next_line_range!(r)
    x === nothing && return nothing
    return String(r.buf[x[1]:x[2]])
end

"""
A batch of FASTQ records kept as bytes: line j (1-4) of record i spans
data[bounds[8(i-1)+2j-1] : bounds[8(i-1)+2j]]. Two allocations per batch instead of
four strings per record.
"""
struct RawBatch
    data::Vector{UInt8}
    bounds::Vector{Int}
end
Base.length(b::RawBatch) = length(b.bounds) ÷ 8
Base.isempty(b::RawBatch) = isempty(b.bounds)

@inline function push_line!(data::Vector{UInt8}, bounds::Vector{Int}, r::ByteLines, x::Tuple{Int,Int})
    s = length(data) + 1
    x[2] >= x[1] && append!(data, @view r.buf[x[1]:x[2]])
    push!(bounds, s); push!(bounds, length(data))
    return nothing
end

"Append one record to a batch; false at EOF (same truncation errors as read_raw_record)."
function read_raw_record!(data::Vector{UInt8}, bounds::Vector{Int}, r::ByteLines, path::AbstractString, record::Int)
    x = next_line_range!(r)
    x === nothing && return false
    push_line!(data, bounds, r, x)
    for part in 2:4
        x = next_line_range!(r)
        x === nothing && error("Truncated FASTQ: $path, record $record, line $part")
        push_line!(data, bounds, r, x)
    end
    return true
end

@inline is_ascii_space(x::UInt8) = x == 0x20 || 0x09 <= x <= 0x0d

"""
parse_fastq_record on record i of a byte batch: same checks, messages, order, base
normalisation and quality masking; returns (sequence, name start, name end) where the name
(token after '@', without a trailing /1 or /2) is a range of b.data, compared without strings.
"""
function parse_fastq_bytes(b::RawBatch, i::Int, path::AbstractString, record::Int, minq::Int)
    o = 8(i-1); d = b.data; bd = b.bounds
    hs, he = bd[o+1], bd[o+2]; ss, se = bd[o+3], bd[o+4]
    ps, pe = bd[o+5], bd[o+6]; qs, qe = bd[o+7], bd[o+8]
    (he - hs + 1 > 1 && d[hs] == UInt8('@')) || error("Invalid FASTQ header: $path, record $record")
    (pe >= ps && d[ps] == UInt8('+')) || error("Missing + line: $path, record $record")
    se - ss == qe - qs || error("Sequence/quality length mismatch: $path, record $record")
    n = se - ss + 1
    seq = Vector{UInt8}(undef, n)
    @inbounds for j in 1:n
        x = d[ss+j-1]
        UInt8('a') <= x <= UInt8('z') && (x -= 0x20)
        if ENCODE[Int(x)+1] == 255
            IUPAC_AMBIGUOUS[Int(x)+1] || error("Invalid DNA byte $(repr(Char(x)))")
            x = UInt8('N')
        end
        seq[j] = x
    end
    n == 0 && error("Empty FASTQ sequence: $path, record $record")
    @inbounds for j in 1:n
        q = d[qs+j-1]
        33 <= q <= 126 || error("Invalid Phred+33 quality: $path, record $record")
        Int(q) - 33 < minq && (seq[j] = UInt8('N'))
    end
    a = hs + 1
    while a <= he && is_ascii_space(d[a])
        a += 1
    end
    a <= he || error("Empty FASTQ identifier: $path, record $record")
    z = a
    while z <= he && !is_ascii_space(d[z])
        z += 1
    end
    z -= 1
    if z - a + 1 >= 2 && d[z-1] == UInt8('/') && (d[z] == UInt8('1') || d[z] == UInt8('2'))
        z -= 2
    end
    return seq, a, z
end

"Four FASTQ lines from a ByteLines reader; same EOF and truncation behaviour as read_raw_record(io)."
function read_raw_record(r::ByteLines, path::AbstractString, record::Int)
    header = next_line!(r)
    header === nothing && return nothing
    text = next_line!(r); text === nothing && error("Truncated FASTQ: $path, record $record, line 2")
    plus = next_line!(r); plus === nothing && error("Truncated FASTQ: $path, record $record, line 3")
    quality = next_line!(r); quality === nothing && error("Truncated FASTQ: $path, record $record, line 4")
    return (header, text, plus, quality)
end

"Reader task: decompress one FASTQ and hand over batches of raw records (bounded buffer)."
struct RawReader
    channel::Channel{RawBatch}
    task::Task
end

"Set by the readers: the input holds more records than Read pairs to process."
const INPUT_MORE = Threads.Atomic{Bool}(false)

function raw_reader(path::AbstractString, limit::Int, batch::Int)
    channel = Channel{RawBatch}(4)
    task = Threads.@spawn begin
        io = open_fastq(path)
        lines = ByteLines(io)
        try
            record = 0
            while limit == 0 || record < limit
                data = UInt8[]; sizehint!(data, batch * 400)
                bounds = Int[]; sizehint!(bounds, 8batch)
                count = 0
                while count < batch && (limit == 0 || record < limit)
                    read_raw_record!(data, bounds, lines, path, record+1) || break
                    record += 1; count += 1
                end
                count > 0 && put!(channel, RawBatch(data, bounds))
                count < batch && break
            end
            limit > 0 && record >= limit && next_line_range!(lines) !== nothing && (INPUT_MORE[] = true)
        finally
            close(io)
        end
    end
    bind(channel, task)                         # channel closes (with the error) when the task ends
    return RawReader(channel, task)
end

"Next batch, or nothing once the file (or the pair limit) is exhausted."
function next_raw_batch(r::RawReader)
    try
        return take!(r.channel)
    catch err
        err isa InvalidStateException && err.state === :closed && (fetch(r.task); return nothing)
        rethrow()
    end
end

finish_reader(r::RawReader) = (fetch(r.task); nothing)                 # surfaces reader errors
stop_reader(r::RawReader) = (isopen(r.channel) && close(r.channel); nothing)

function build_index(c::Config; quiet::Bool=false, monitor=nothing)
    progress_phase!(monitor,"Indexing";total=c.read_pairs_to_process,unit=c.paired ? "pairs limit" : "reads limit")
    store = PairedReadStore(c.paired)
    nthreads = Threads.nthreads()
    sbits = index_shard_bits(nthreads)
    workspaces = [MinimizerWorkspace() for _ in 1:nthreads]
    # buffers[chunk][shard]: written only by the task processing that chunk.
    buffers = [[UInt128[] for _ in 1:(1 << sbits)] for _ in 1:nthreads]
    maxlen = 0; short_reads = 0
    empty!(PHASE_TIMES)
    INPUT_MORE[] = false
    t_wait = 0.0; t_parse = 0.0; t_store = 0.0; t_min = 0.0; gc_min = 0.0
    gc_start = Base.gc_time_ns(); reserved_until = 0; t_reserve = 0.0
    # One reader task per file only decompresses and splits lines; parsing, packing
    # and minimizers run on all threads. Same reads, ids and index as before.
    p1 = raw_reader(c.forward_reads, c.read_pairs_to_process, c.batch_size)
    p2 = c.paired ? raw_reader(c.reverse_reads, c.read_pairs_to_process, c.batch_size) : nothing
    try
        while true
            t0 = time()
            b1 = next_raw_batch(p1)
            b2 = p2 === nothing ? nothing : next_raw_batch(p2)
            t_wait += time()-t0
            if b1 === nothing
                b2 !== nothing && !isempty(b2) && error("R2 has more records than R1")
                break
            end
            first_pair = num_pairs(store) + 1
            npairs = length(b1)
            if p2 !== nothing
                got = b2 === nothing ? 0 : length(b2)
                got < npairs && error("R2 ended before R1 at pair $(first_pair+got)")
                got > npairs && error("R2 has more records than R1")
            end
            first_pair + npairs - 1 <= (typemax(UInt32) >> 1) || error("Pair identifier limit reached")
            t0 = time()
            seqs1 = Vector{Vector{UInt8}}(undef, npairs)
            n1s = Vector{Int}(undef, npairs); n1e = Vector{Int}(undef, npairs)
            m2 = p2 === nothing ? 0 : npairs
            seqs2 = Vector{Vector{UInt8}}(undef, m2)
            n2s = Vector{Int}(undef, m2); n2e = Vector{Int}(undef, m2)
            chunks = min(nthreads, npairs)
            Threads.@threads :dynamic for chunk in 1:chunks
                for i in fld((chunk-1)*npairs,chunks)+1:fld(chunk*npairs,chunks)
                    n = first_pair + i - 1
                    seqs1[i], n1s[i], n1e[i] = parse_fastq_bytes(b1, i, c.forward_reads, n, c.min_base_quality)
                    if p2 !== nothing
                        seqs2[i], n2s[i], n2e[i] = parse_fastq_bytes(b2, i, c.reverse_reads, n, c.min_base_quality)
                    end
                end
            end
            if p2 !== nothing
                for i in 1:npairs
                    @view(b1.data[n1s[i]:n1e[i]]) == @view(b2.data[n2s[i]:n2e[i]]) ||
                        error("R1/R2 names differ at pair $(first_pair+i-1): $(String(b1.data[n1s[i]:n1e[i]])) / $(String(b2.data[n2s[i]:n2e[i]]))")
                end
            end
            t_parse += time()-t0
            t0 = time()
            append_reads!(store.r1, seqs1)
            p2 === nothing || append_reads!(store.r2, seqs2)
            t_store += time()-t0
            batch = Tuple{UInt32,Vector{UInt8}}[]
            sizehint!(batch, npairs*(p2 === nothing ? 1 : 2))
            for i in 1:npairs
                n = first_pair + i - 1
                push!(batch, (read_id(n, true), seqs1[i]))
                p2 === nothing || push!(batch, (read_id(n, false), seqs2[i]))
            end
            for (_, seq) in batch
                maxlen = max(maxlen, length(seq))
                length(seq) < c.kmer + c.minimizer_window - 1 && (short_reads += 1)
            end
            # Each logical chunk owns its workspace and its shard buffers, even if
            # Julia migrates its task; the loop barrier permits reuse next batch.
            t0 = time(); gc_before = Base.gc_time_ns()
            chunks = min(nthreads, length(batch))
            Threads.@threads :dynamic for chunk in 1:chunks
                work = workspaces[chunk]
                out = buffers[chunk]
                first = fld((chunk-1)*length(batch),chunks)+1
                last = fld(chunk*length(batch),chunks)
                for i in first:last
                    id, seq = batch[i]
                    for key in minimizer_keys!(work,seq,c.kmer,c.minimizer_window)
                        h = index_hash(key)
                        push!(out[shard_index(h, sbits)], pack_posting(h, id))
                    end
                end
            end
            t_min += time()-t0
            gc_min += (Base.gc_time_ns() - gc_before)/1e9
            # After the first batch the occurrences per pair are known: every shard buffer is
            # reserved once for the whole read limit instead of growing by repeated doubling.
            # Reserved but untouched pages are not resident (not in RSS), so a file shorter than
            # the limit costs address space only. Staged reservation (16.19.0) re-copied the
            # buffers at each stage and was measured worse (L01: +10 s, peak +3.6 GB).
            if reserved_until == 0 && c.read_pairs_to_process > 0
                reserved_until = 1
                tr = time()
                done_pairs = num_pairs(store)
                occ = sum(length(v) for b in buffers for v in b)
                per = ceil(Int, 1.1 * occ / max(done_pairs, 1) * c.read_pairs_to_process /
                                (length(buffers) * (1 << sbits)))
                for b in buffers, v in b
                    sizehint!(v, per)
                end
                t_reserve += time() - tr
            end
            progress_update!(monitor;done=num_pairs(store),detail="indexing batch")
            !quiet && num_pairs(store) % (10*c.batch_size) < npairs && @printf("Indexed %d pairs\n", num_pairs(store))
        end
        finish_reader(p1)
        p2 === nothing || finish_reader(p2)
    finally
        stop_reader(p1)
        p2 === nothing || stop_reader(p2)
    end
    PHASE_TIMES["input_wait"] = t_wait; PHASE_TIMES["parse"] = t_parse
    PHASE_TIMES["read_store"] = t_store; PHASE_TIMES["minimizers"] = t_min
    PHASE_TIMES["gc_in_minimizers"] = gc_min
    PHASE_TIMES["reserve"] = t_reserve
    PHASE_TIMES["gc_reading_and_minimizers"] = (Base.gc_time_ns() - gc_start)/1e9
    num_pairs(store) > 0 || error("No reads loaded")
    # EOF can precede the configured read limit; only now is the actual total known.
    progress_update!(monitor;done=num_pairs(store),total=num_pairs(store),detail="input complete",force=true)
    progress_phase!(monitor,"Filtering";total=1 << sbits,unit="shards")
    nreads = num_pairs(store) * (c.paired ? 2 : 1)
    cutoff = min(c.max_postings, max(1, floor(Int, nreads*c.max_minimizer_frequency)))
    t0 = time()
    postings, distinct_keys, dropped = freeze_sharded!(buffers, sbits, cutoff; monitor=monitor)
    PHASE_TIMES["index_freeze"] = time()-t0
    empty!(buffers)
    progress_update!(monitor;done=1 << sbits,detail="$(distinct_keys) distinct keys",force=true)
    isempty(postings) && error("All minimizers were filtered; inspect reads or increase Max minimizer frequency")
    !quiet && @printf("Loaded %d pairs; %d minimizers retained; %d hot keys removed; global cutoff %d reads\n",
                      num_pairs(store), length(postings), dropped, cutoff)
    if short_reads > 0 && (!quiet || monitor!==nothing)
        progress_breakline!(monitor)
        @warn "Some reads are too short for the minimizer window" short_reads
    end
    return ReadIndex(store, postings, maxlen, cutoff, dropped, short_reads)
end

"Per-task deduplication scratch: stamps[id] == stamp marks an id already collected."
mutable struct CandidateWorkspace
    stamps::Vector{UInt16}           # 2 bytes per read id; reset on wrap-around.
    stamp::UInt16
    ids::Vector{UInt32}
    keys::Vector{UInt64}             # reused k-mer buffer for all_kmer_keys!
end

# Task-local, so concurrent assemblies never share it and a migrating task keeps its own.
function candidate_workspace(capacity::Int)
    tls = task_local_storage()
    work = get(tls, :assembler_candidate_workspace, nothing)
    if !(work isa CandidateWorkspace)
        work = CandidateWorkspace(zeros(UInt16, capacity), UInt16(0), UInt32[], UInt64[])
        tls[:assembler_candidate_workspace] = work
    elseif length(work.stamps) < capacity
        old = length(work.stamps)
        resize!(work.stamps, capacity)
        fill!(@view(work.stamps[old+1:end]), UInt16(0))
    end
    return work::CandidateWorkspace
end

# Same result as the former Set union (sorted distinct ids; saturated as soon as the
# distinct count exceeds Max candidates), without hashing every posting.
"True for a k-mer with period 1-6 over its whole length ((TA)n, poly(A), (AAT)n ...)."
@inline function low_complexity_kmer(key::UInt64, k::Int)
    for p in 1:min(6, k-1)
        low = k-p == 32 ? typemax(UInt64) : (UInt64(1) << (2(k-p))) - UInt64(1)
        (key >> (2p)) == (key & low) && return true
    end
    return false
end

# Low-complexity k-mers are skipped: they occur at microsatellites all over the (mostly
# nuclear) genome, recruit unrelated reads and can saturate Max candidates.
function candidate_read_ids(seq::AbstractVector{UInt8}, idx::ReadIndex, c::Config)
    work = candidate_workspace(2*num_pairs(idx.store)+1)
    stamp = work.stamp + UInt16(1)
    if stamp == 0                      # Wrapped after 65535 calls: reset marks (one memset).
        fill!(work.stamps, UInt16(0)); stamp = UInt16(1)
    end
    work.stamp = stamp
    stamps = work.stamps
    ids = work.ids
    empty!(ids)
    for key in all_kmer_keys!(work.keys, seq, c.kmer)
        low_complexity_kmer(key, c.kmer) && continue
        bucket = get(idx.postings, key, nothing)
        bucket === nothing && continue
        for id in bucket
            excluded_read(UInt32(id),c) && continue
            if id > length(stamps)
                old = length(stamps)
                resize!(stamps, Int(id)); fill!(@view(stamps[old+1:end]), UInt16(0))
            end
            if stamps[id] != stamp
                stamps[id] = stamp
                push!(ids, id)
                length(ids) > c.max_candidates && return UInt32[], true
            end
        end
    end
    return sort!(copy(ids)), false
end

# N never counts as a matching base, including N versus N.
const GAP_BAND = 3                   # max net shift (bases) of the gapped comparisons

@inline edit_cost(x::UInt8, y::UInt8) = (x != y || x == UInt8('N') || y == UInt8('N')) ? 1 : 0

"""
Edit distance of a and b (global) inside a diagonal band of +-GAP_BAND; lengths may differ by
at most the band. typemax(Int) if outside the band. Used only after an ungapped comparison
failed, so homopolymer or microsatellite length differences (indels) are not mismatches.
"""
function banded_edits(a::AbstractVector{UInt8}, b::AbstractVector{UInt8}; band::Int=GAP_BAND)
    n = length(a); m = length(b)
    abs(n - m) <= band || return typemax(Int)
    W = 2band + 1; INF = typemax(Int) ÷ 2
    prev = fill(INF, W); cur = fill(INF, W)          # column index t: j = i + t - band - 1
    for t in 1:W
        j = t - band - 1
        0 <= j <= m && (prev[t] = j)                  # row i = 0
    end
    for i in 1:n
        fill!(cur, INF)
        for t in 1:W
            j = i + t - band - 1
            (0 <= j <= m) || continue
            if j == 0
                cur[t] = i; continue
            end
            best = prev[t] + edit_cost(a[i], b[j])    # diagonal: (i-1, j-1) has the same t
            t > 1 && (best = min(best, cur[t-1] + 1))  # (i, j-1)
            t < W && (best = min(best, prev[t+1] + 1)) # (i-1, j)
            cur[t] = best
        end
        prev, cur = cur, prev
    end
    t = m - n + band + 1
    return prev[t] >= INF ? typemax(Int) : prev[t]
end

"""
Gapped overlap of a read with the end of the context: the last `ov` context bases align
(global) to a read prefix whose length may differ by +-GAP_BAND. Returns the read length
consumed (0 if more than `allowed` edits, no read base left to extend, or a substitution in
the last `w` aligned context bases: phasing still forbids other SNV alleles, not indels).
"""
function gapped_overlap(ctx::AbstractVector{UInt8}, read::AbstractVector{UInt8}, ov::Int, allowed::Int, w::Int;
                        band::Int=GAP_BAND)
    a = @view ctx[end-ov+1:end]
    mmax = min(length(read) - 1, ov + band)
    mmax >= ov - band || return 0
    W = 2band + 1; INF = typemax(Int32) ÷ Int32(2)
    B = fill(INF, ov + 1, W)                          # B[i+1, t] holds cell (i, j = i + t - band - 1)
    for t in 1:W
        j = t - band - 1
        0 <= j <= mmax && (B[1, t] = Int32(j))
    end
    for i in 1:ov, t in 1:W
        j = i + t - band - 1
        (0 <= j <= mmax) || continue
        if j == 0
            B[i+1, t] = Int32(i); continue
        end
        d = B[i, t] + Int32(edit_cost(a[i], read[j]))
        t > 1 && (d = min(d, B[i+1, t-1] + Int32(1)))
        t < W && (d = min(d, B[i, t+1] + Int32(1)))
        B[i+1, t] = d
    end
    jbest = 0; ebest = INF
    for j in max(1, ov-band):mmax
        e = B[ov+1, j - ov + band + 1]
        e < ebest && (ebest = e; jbest = j)
    end
    (jbest > 0 && ebest <= allowed) || return 0
    # Traceback over the last w context bases: no substitution allowed there.
    i = ov; j = jbest
    while i > ov - w && i > 0 && j > 0
        t = j - i + band + 1
        d = B[i+1, t]
        cost = edit_cost(a[i], read[j])
        if d == B[i, t] + cost
            cost == 1 && return 0
            i -= 1; j -= 1
        elseif t < W && d == B[i, t+1] + 1
            i -= 1
        else
            j -= 1
        end
    end
    return jbest
end

function sequences_similar(a::AbstractVector{UInt8}, b::AbstractVector{UInt8}, identity::Float64)
    length(a) == length(b) || return false
    isempty(a) && return false
    allowed = floor(Int, length(a)*(1-identity) + 1e-9)
    mismatches = 0
    @inbounds for i in eachindex(a, b)
        if a[i] != b[i] || a[i] == UInt8('N') || b[i] == UInt8('N')
            mismatches += 1
            mismatches > allowed && return false
        end
    end
    return true
end

dna_string(seq::AbstractVector{UInt8}) = String(copy(seq))

"Exact forward anchor keys; each linked list enumerates positions in ascending order."
struct ExactAnchorTable
    heads::Dict{UInt64,Int}
    next::Vector{Int}
end

function ExactAnchorTable(text::String,k::Int)
    1 <= k <= 32 || throw(ArgumentError("anchor size must be in 1:32"))
    seq = codeunits(text)
    heads = Dict{UInt64,Int}()
    links = zeros(Int,max(0,length(seq)-k+1))
    key = UInt64(0); run = 0
    # Right-to-left construction prepends progressively smaller coordinates.
    # Enumeration therefore matches findnext, including early placement caps.
    @inbounds for i in length(seq):-1:1
        b = ENCODE[Int(seq[i])+1]
        if b == 255
            key = 0; run = 0
            continue
        end
        key = (key >> 2) | (UInt64(b) << (2(k-1)))
        run += 1
        if run >= k
            links[i] = get(heads,key,0)
            heads[key] = i
        end
    end
    return ExactAnchorTable(heads,links)
end

"One local terminal window; its lazy tables are never shared between assembly tasks."
struct ExactAnchorCache
    text::String
    tables::Dict{Int,ExactAnchorTable}
end
ExactAnchorCache(text::String) = ExactAnchorCache(text,Dict{Int,ExactAnchorTable}())

@inline function exact_anchor_key(anchor::AbstractVector{UInt8})
    key = UInt64(0)
    @inbounds for b in anchor
        encoded = ENCODE[Int(b)+1]
        encoded == 255 && return nothing
        key = (key << 2) | UInt64(encoded)
    end
    return key
end

struct AnchorHits
    next::Vector{Int}
    first::Int
    last::Int
end
@inline function Base.iterate(hits::AnchorHits, pos::Int=hits.first)
    (pos == 0 || pos > hits.last) && return nothing
    return pos, @inbounds(hits.next[pos])
end
Base.IteratorSize(::Type{AnchorHits}) = Base.SizeUnknown()
Base.eltype(::Type{AnchorHits}) = Int

function anchor_hits(cache::ExactAnchorCache,k::Int,key::UInt64,lo::Int,hi::Int)
    table = get!(cache.tables,k) do
        ExactAnchorTable(cache.text,k)
    end
    pos = get(table.heads,key,0)
    while pos != 0 && pos < lo
        pos = table.next[pos]
    end
    return AnchorHits(table.next,pos,hi)
end

function indexed_placements(cache::ExactAnchorCache,read::AbstractVector{UInt8},identity::Float64;
                            lo::Int,hi::Int,anchor_size::Int,cap::Int)
    text = cache.text
    lo = max(1,lo); hi = min(hi,ncodeunits(text)-length(read)+1)
    lo > hi && return Int[]
    a = min(anchor_size,length(read))
    # Preserve the public search API for unusually long anchors.
    if !(1 <= a <= 32)
        return oriented_placements(text,read,identity;lo=lo,hi=hi,anchor_size=anchor_size,cap=cap)
    end
    starts = (1,(length(read)-a)÷2+1,length(read)-a+1)
    seen = falses(hi-lo+1)
    target = codeunits(text)
    found = Int[]
    for anchor_start in starts
        key = exact_anchor_key(@view(read[anchor_start:anchor_start+a-1]))
        key === nothing && continue
        for hit in anchor_hits(cache,a,key,lo+anchor_start-1,hi+anchor_start-1)
            p = hit-anchor_start+1
            seen[p-lo+1] && continue
            seen[p-lo+1] = true
            if sequences_similar(@view(target[p:p+length(read)-1]),read,identity)
                push!(found,p)
                length(found) >= cap && return sort!(found)
            end
        end
    end
    return sort!(found)
end

"Find full-length, ungapped placements using three exact short anchors."
function oriented_placements(text::String, read::AbstractVector{UInt8}, identity::Float64;
                             lo::Int=1, hi::Int=ncodeunits(text)-length(read)+1,
                             anchor_size::Int=15, cap::Int=64,
                             anchors::Union{Nothing,ExactAnchorCache}=nothing)
    anchors === nothing || return indexed_placements(anchors,read,identity;
        lo=lo,hi=hi,anchor_size=anchor_size,cap=cap)
    lo = max(1, lo); hi = min(hi, ncodeunits(text)-length(read)+1)
    lo > hi && return Int[]
    a = min(anchor_size, length(read))
    starts = unique([1, (length(read)-a)÷2+1, length(read)-a+1])
    proposals = Set{Int}()
    target = codeunits(text)
    found = Int[]
    for anchor_start in starts
        anchor = @view read[anchor_start:anchor_start+a-1]
        UInt8('N') in anchor && continue
        needle = dna_string(anchor)
        cursor = lo+anchor_start-1
        while cursor <= ncodeunits(text)-a+1
            hit = findnext(needle, text, cursor)
            hit === nothing && break
            p = first(hit)-anchor_start+1
            p > hi && break
            if p >= lo && !(p in proposals)
                push!(proposals, p)
                if sequences_similar(@view(target[p:p+length(read)-1]), read, identity)
                    push!(found, p)
                    length(found) >= cap && return sort!(found)
                end
            end
            cursor = first(hit)+1
        end
    end
    return sort!(found)
end

struct Placement
    start::Int
    len::Int
    reverse::Bool
end

function read_placements(text::String, read::Vector{UInt8}, c::Config;
                         anchors::Union{Nothing,ExactAnchorCache}=nothing)
    result = Placement[]
    for reverse in (false, true)
        oriented = reverse ? reverse_complement(read) : read
        for p in oriented_placements(text, oriented, c.min_identity;anchors=anchors)
            push!(result, Placement(p, length(read), reverse))
        end
    end
    return result
end

function insert_bounds(c::Config)
    return max(1, floor(Int, c.insert_size*(1-c.insert_tolerance))),
           ceil(Int, c.insert_size*(1+c.insert_tolerance))
end

function mate_supported(context::Vector{UInt8}, text::String, idx::ReadIndex,
                        id::UInt32, reverse::Bool, read_start::Int, read_len::Int,
                        c::Config; anchors::Union{Nothing,ExactAnchorCache}=nothing)
    (!c.paired || !c.paired_validation) && return false
    mate = get_read(idx.store, id ⊻ UInt32(1))
    lower, upper = insert_bounds(c)
    if reverse
        # Reverse-oriented right read: its forward mate must be to the left.
        lo = read_start+read_len-upper
        hi = min(read_start, read_start+read_len-lower)
        oriented = mate
    else
        # Forward-oriented left read: its reverse mate must be to the right.
        lo = max(read_start, read_start+lower-length(mate))
        hi = read_start+upper-length(mate)
        oriented = reverse_complement(mate)
    end
    positions = oriented_placements(text, oriented, c.min_identity; lo=lo, hi=hi, cap=2,anchors=anchors)
    return length(positions) == 1
end

struct ExtensionCandidate
    id::UInt32
    extension::Vector{UInt8}
    overlap::Int
    reverse::Bool
    paired_support::Bool
end

"Candidate overlap lengths from disjoint anchors; exhaustive fallback for loose thresholds."
function proposed_overlaps(text::String, read::Vector{UInt8}, c::Config;
                           anchors::Union{Nothing,ExactAnchorCache}=nothing,
                           scratch::Vector{Int}=Int[])
    empty!(scratch)
    upper = min(ncodeunits(text), length(read))
    upper < c.min_overlap && return scratch
    max_errors = floor(Int, upper*(1-c.min_identity)+1e-9)
    anchor_len = min(15, c.min_overlap÷(max_errors+1))
    # e+1 disjoint anchors ensure at least one exact anchor at <= e mismatches.
    # With extremely permissive thresholds, exhaustive checking is safer.
    anchor_len < 8 && return append!(scratch,upper:-1:c.min_overlap)
    for j in 0:max_errors
        p = 1+j*anchor_len
        anchor = @view read[p:p+anchor_len-1]
        UInt8('N') in anchor && continue
        if anchors !== nothing
            key = exact_anchor_key(anchor)
            key === nothing && continue
            for hit in anchor_hits(anchors,anchor_len,key,
                                   ncodeunits(text)-upper+p,ncodeunits(text)-anchor_len+1)
                overlap = ncodeunits(text)-hit+p
                c.min_overlap <= overlap <= upper && push!(scratch,overlap)
            end
            continue
        end
        needle = dna_string(anchor)
        cursor = ncodeunits(text)-upper+p
        while cursor <= ncodeunits(text)-anchor_len+1
            hit = findnext(needle, text, cursor)
            hit === nothing && break
            overlap = ncodeunits(text)-first(hit)+p
            c.min_overlap <= overlap <= upper && push!(scratch, overlap)
            cursor = first(hit)+1
        end
    end
    sort!(scratch; rev=true, alg=QuickSort)      # short lists: cheaper than the default hybrid
    w = 0
    @inbounds for v in scratch                     # sorted: duplicates are adjacent
        if w == 0 || scratch[w] != v
            w += 1; scratch[w] = v
        end
    end
    return resize!(scratch, w)
end

"Reads flagged as contaminant in the decontamination pass (empty: no filtering)."
const FOREIGN_READS = Ref(BitVector())
@inline is_foreign_read(id::UInt32) = (f = FOREIGN_READS[]; !isempty(f) && id <= length(f) && f[id])
@inline excluded_read(id::UInt32,c::Config) = is_foreign_read(id) || (c.paired && is_foreign_read(id ⊻ UInt32(1)))

function extension_candidates(context::Vector{UInt8}, idx::ReadIndex, c::Config; phase::Int=c.phase_window)
    # A matching mate can recruit an extending read with no exact indexed k-mer
    # in its overlap. In SE, only the terminal read-length window is relevant.
    terminal_len = min(length(context), idx.max_read_length)
    query = c.paired && c.paired_validation ? context : @view(context[end-terminal_len+1:end])
    ids, saturated = candidate_read_ids(query, idx, c)
    saturated && return ExtensionCandidate[], "candidate_limit"
    if c.paired && c.paired_validation
        ids = sort!(unique(vcat(ids, [id ⊻ UInt32(1) for id in ids])))
        length(ids) > c.max_candidates && return ExtensionCandidate[], "candidate_limit"
    end
    candidates = ExtensionCandidate[]
    text = dna_string(context)
    anchors = ExactAnchorCache(text)
    overlaps = Int[]
    # A context ending in a tandem repeat: an overlap lying only in the repeat would match
    # any copy of that microsatellite, so the read must cross it into unique sequence.
    period, run = periodic_suffix(context)
    min_overlap = period > 0 ? max(c.min_overlap, run + c.repeat_anchor) : c.min_overlap
    # Two decode buffers reused across all candidate reads of this call; only the
    # accepted extension suffix is copied into an ExtensionCandidate.
    read = UInt8[]; read_rc = UInt8[]
    for id in ids
        excluded_read(id,c) && continue                 # include flagged mate, not just this read
        get_read!(read, idx.store, id)
        reverse_complement!(read_rc, read)
        for reverse in (false, true)
            oriented = reverse ? read_rc : read
            for overlap in proposed_overlaps(text, oriented, c;anchors=anchors,scratch=overlaps)
                w = min(phase, overlap)
                used = 0                                        # read bases aligned to the context end
                if sequences_similar(@view(context[end-overlap+1:end]), @view(oriented[1:overlap]), c.min_identity)
                    overlap >= min_overlap || break             # largest overlap already too short
                    # Phasing: a read differing within the last Phase window bases carries another
                    # allele of a variant just added; letting it vote for the next bases would build
                    # a combination of alleles that no read contains.
                    w == 0 || @view(context[end-w+1:end]) == @view(oriented[overlap-w+1:overlap]) || break
                    used = overlap
                elseif overlap >= min_overlap
                    # Gapped fallback: the read differs by a short indel (e.g. homopolymer length);
                    # same identity budget, counted as edits; phasing forbids substitutions only.
                    allowed = floor(Int, overlap*(1-c.min_identity) + 1e-9)
                    used = allowed == 0 ? 0 : gapped_overlap(context, oriented, overlap, allowed, w)
                    used == 0 && continue
                else
                    continue
                end
                if used < length(oriented)
                    supported = mate_supported(context, text, idx, id, reverse,
                                               length(context)-overlap+1, length(oriented), c;anchors=anchors)
                    push!(candidates, ExtensionCandidate(id, oriented[used+1:end],
                                                        overlap, reverse, supported))
                end
                # Do not reinterpret a fully contained read as a shorter overlap.
                break
            end
        end
    end
    sort!(candidates; by=x -> (x.id, x.reverse, -x.overlap))
    return candidates, isempty(candidates) ? "no_extension" : "ok"
end

"Resolved bubble inside one extension: alleles in extension orientation."
struct BubbleVariant
    offset::Int                 # 1-based position in the extension where the alleles start
    major::Vector{UInt8}
    minor::Vector{UInt8}
    major_support::Int          # read pairs voting for the major / minor base at divergence
    minor_support::Int
    informative::Int            # all pairs voting at that position (conflicting ones included)
    switch::Bool                # coverage switch: major = coverage-consistent base, minor = majority
end
BubbleVariant(offset, major, minor, ms, mn, informative) = BubbleVariant(offset, major, minor, ms, mn, informative, false)

struct ExtensionChoice
    sequence::Vector{UInt8}
    min_support::Int
    reason::String
    variants::Vector{BubbleVariant}
end
ExtensionChoice(sequence, min_support, reason) = ExtensionChoice(sequence, min_support, reason, BubbleVariant[])

"Pair-level ballots; contradictory placements of one pair abstain."
function extension_votes(candidates::Vector{ExtensionCandidate}, active::Vector{Int}, pos::Int)
    return extension_votes!(VoteWorkspace(),candidates,active,pos)
end

struct VoteWorkspace
    masks::Dict{UInt32,UInt8}
    paired_masks::Dict{UInt32,UInt8}
    counts::Vector{Int}
    paired::Vector{Int}
end
VoteWorkspace() = VoteWorkspace(Dict{UInt32,UInt8}(),Dict{UInt32,UInt8}(),zeros(Int,4),zeros(Int,4))

function extension_votes!(work::VoteWorkspace,candidates::Vector{ExtensionCandidate},active::Vector{Int},pos::Int)
    masks = work.masks; paired_masks = work.paired_masks
    empty!(masks); empty!(paired_masks)
    for j in active
        item = candidates[j]
        pos <= length(item.extension) || continue
        b = ENCODE[Int(item.extension[pos])+1]
        b == 255 && continue
        p = pair_id(item.id)
        bit = UInt8(1) << b
        masks[p] = get(masks, p, UInt8(0)) | bit
        if item.paired_support
            paired_masks[p] = get(paired_masks, p, UInt8(0)) | bit
        end
    end
    counts = work.counts; paired = work.paired
    fill!(counts,0); fill!(paired,0)
    for (p, mask) in masks
        count_ones(mask) == 1 || continue
        b = trailing_zeros(mask)+1
        counts[b] += 1
        get(paired_masks, p, UInt8(0)) == mask && (paired[b] += 1)
    end
    return counts, paired, masks
end

function compatible_active(candidates, active, masks, pos, base_index)
    bit = UInt8(1) << (base_index-1)
    return [j for j in active if pos <= length(candidates[j].extension) &&
            candidates[j].extension[pos] == DNA[base_index] &&
            get(masks, pair_id(candidates[j].id), UInt8(0)) == bit]
end

"Follow the reads carrying `base` at `pos`, taking the plurality base, for up to Bubble lookahead bases."
function bubble_lookahead(candidates, active, masks, pos::Int, base::Int, c::Config, look::VoteWorkspace)
    seq = UInt8[DNA[base]]
    act = compatible_active(candidates, active, masks, pos, base)
    p = pos + 1
    while length(seq) < c.bubble_lookahead
        cnt, _, m = extension_votes!(look, candidates, act, p)
        (length(m) < c.min_coverage || maximum(cnt) < c.min_coverage) && break
        b = 1
        for x in 2:4
            cnt[x] > cnt[b] && (b = x)
        end
        push!(seq, DNA[b])
        act = compatible_active(candidates, act, m, p, b)
        p += 1
    end
    return seq
end

"First shared anchor after two alleles: (sa, sb) with a[sa:...] == b[sb:...], or nothing."
function bubble_reconvergence(a::Vector{UInt8}, b::Vector{UInt8}, anchor::Int, max_indel::Int)
    best = nothing
    for sa in 1:max_indel+1, sb in 1:max_indel+1
        (sa == 1 && sb == 1) && continue                 # alleles differ at their first base
        abs(sa - sb) <= max_indel || continue
        (sa+anchor-1 <= length(a) && sb+anchor-1 <= length(b)) || continue
        @view(a[sa:sa+anchor-1]) == @view(b[sb:sb+anchor-1]) || continue
        key = (sa+sb, abs(sa-sb), sa)
        (best === nothing || key < best[1]) && (best = (key, sa, sb))
    end
    return best === nothing ? nothing : (best[2], best[3])
end

"""
Non-decisive position: follow the majority if every supported alternative is a bubble
(reconverges within Bubble max indel) or too rare to branch on. Returns
(record, major, minor, minor_support) or nothing for a genuine fork.
"""
function resolve_bubble(candidates, active, masks, counts, ranked, pos::Int, total::Int,
                        c::Config, look::VoteWorkspace)
    best = ranked[1]
    counts[best] >= c.min_coverage || return nothing
    alternatives = [b for b in ranked[2:4] if counts[b] >= c.min_coverage]
    if isempty(alternatives)
        # Only conflicting ballots against the plurality base: accept a clear majority.
        return counts[best]/total > 0.5 ? (false, UInt8[], UInt8[], 0) : nothing
    end
    major = bubble_lookahead(candidates, active, masks, pos, best, c, look)
    found = nothing
    for b in alternatives
        minor = bubble_lookahead(candidates, active, masks, pos, b, c, look)
        m = length(major) >= c.bubble_anchor+1 && length(minor) >= c.bubble_anchor+1 ?
            bubble_reconvergence(major, minor, c.bubble_anchor, c.bubble_max_indel) : nothing
        if m !== nothing
            found === nothing && (found = (true, major[1:m[1]-1], minor[1:m[2]-1], counts[b]))
        elseif counts[b]/total >= c.branch_min_fraction
            return nothing                              # divergent and frequent: real branch
        end
    end
    return found === nothing ? (counts[best]/total > 0.5 ? (false, UInt8[], UInt8[], 0) : nothing) : found
end

"Follow coherent read prefixes; never construct a column-wise mosaic."
function consensus_choices(candidates::Vector{ExtensionCandidate}, c::Config; expected::Int=0)
    isempty(candidates) && return ExtensionChoice[]
    pending = [(UInt8[], collect(eachindex(candidates)), typemax(Int), BubbleVariant[])]
    results = ExtensionChoice[]
    # Branches retain only prefixes/active IDs, never these borrowed vote buffers.
    votes = VoteWorkspace()
    look = VoteWorkspace()                      # separate buffers: lookahead must not clobber `votes`
    while !isempty(pending)
        prefix, active, support, variants = popfirst!(pending)
        reason = "max_extension"
        while length(prefix) < c.max_extension
            pos = length(prefix)+1
            counts, paired, masks = extension_votes!(votes,candidates,active,pos)
            total = length(masks) # Includes conflicting pairs in denominator.
            if total < c.min_coverage || maximum(counts) < c.min_coverage
                reason = "low_support"; break
            end
            ranked = sort(collect(1:4); by=b -> (-counts[b], b))
            best, second = ranked[1], ranked[2]
            # Coverage switch: a majority far above this path's support, next to an alternative
            # whose support matches it, is another genome sharing this stretch (contamination).
            if expected > 0 && c.coverage_switch_ratio > 0 && counts[best] > c.coverage_switch_ratio*expected
                alt = 0
                for b in ranked[2:4]
                    if counts[b] >= c.min_coverage && expected/4 <= counts[b] <= 2*expected
                        alt = b; break
                    end
                end
                # Alleles of one genome reconverge (a bubble): that is a mixture or heteroplasmy,
                # handled by the majority rules below, never a reason to switch.
                if alt != 0
                    a = bubble_lookahead(candidates, active, masks, pos, best, c, look)
                    b = bubble_lookahead(candidates, active, masks, pos, alt, c, look)
                    (length(a) > c.bubble_anchor && length(b) > c.bubble_anchor &&
                     bubble_reconvergence(a, b, c.bubble_anchor, c.bubble_max_indel) !== nothing) && (alt = 0)
                end
                if alt != 0
                    push!(variants, BubbleVariant(pos, UInt8[DNA[alt]], UInt8[DNA[best]], counts[alt],
                                                  counts[best], total, true))
                    push!(prefix, DNA[alt])
                    support = min(support, counts[alt])
                    active = compatible_active(candidates, active, masks, pos, alt)
                    continue
                end
            end
            decisive = counts[best]/total >= c.consensus_fraction &&
                       (counts[second] == 0 || counts[best]/counts[second] >= c.consensus_ratio)
            if !decisive && c.paired && c.paired_validation
                pranked = sort(collect(1:4); by=b -> (-paired[b], b))
                pbest, psecond = pranked[1], pranked[2]
                if paired[pbest] >= c.min_paired_support && counts[pbest] >= c.min_coverage &&
                   (paired[psecond] == 0 || paired[pbest]/paired[psecond] >= c.consensus_ratio)
                    best = pbest; decisive = true
                end
            end
            if decisive
                push!(prefix, DNA[best])
                support = min(support, counts[best])
                active = compatible_active(candidates, active, masks, pos, best)
                continue
            end
            if c.bubble_policy == "majority"
                bubble = resolve_bubble(candidates, active, masks, counts, ranked, pos, total, c, look)
                if bubble !== nothing
                    major = ranked[1]
                    bubble[1] && push!(variants, BubbleVariant(pos, bubble[2], bubble[3],
                                                              counts[major], bubble[4], total))
                    push!(prefix, DNA[major])
                    support = min(support, counts[major])
                    active = compatible_active(candidates, active, masks, pos, major)
                    continue
                end
            end
            branches = [b for b in ranked if counts[b] >= c.min_coverage &&
                        counts[b]/total >= c.branch_min_fraction]
            if length(branches) >= 2 && length(results)+length(pending)+length(branches) <= c.max_branch_paths
                for b in branches
                    next_prefix = copy(prefix); push!(next_prefix, DNA[b])
                    push!(pending, (next_prefix, compatible_active(candidates, active, masks, pos, b),
                                    min(support, counts[b]), copy(variants)))
                end
                reason = "forked"
            else
                reason = "ambiguous"
            end
            break
        end
        if reason != "forked" && !isempty(prefix)
            push!(results, ExtensionChoice(prefix, support, reason, variants))
        end
    end
    sort!(results; by=x -> (-x.min_support, dna_string(x.sequence)))
    return results
end

# Two buffers avoid copying a growing contig on every left/right extension.
# left holds the left flank in REVERSED (not complemented) order.
mutable struct ContigBuffer
    left::Vector{UInt8}
    right::Vector{UInt8}
end
ContigBuffer(seed::Vector{UInt8}) = ContigBuffer(UInt8[], copy(seed))
Base.length(b::ContigBuffer) = length(b.left)+length(b.right)
Base.copy(b::ContigBuffer) = ContigBuffer(copy(b.left), copy(b.right))
contig_sequence(b::ContigBuffer) = vcat(reverse(b.left), b.right)

function contig_head(b::ContigBuffer, n::Int)
    n = min(n, length(b))
    left_take = min(n, length(b.left))
    out = reverse(b.left[end-left_take+1:end])
    n > left_take && append!(out, @view b.right[1:n-left_take])
    return out
end

function contig_tail(b::ContigBuffer, n::Int)
    n = min(n, length(b))
    right_take = min(n, length(b.right))
    out = n > right_take ? reverse(b.left[1:n-right_take]) : UInt8[]
    append!(out, @view b.right[end-right_take+1:end])
    return out
end

function append_extension!(b::ContigBuffer, sequence::Vector{UInt8}, direction::Symbol)
    if direction == :right
        append!(b.right, sequence)
    else
        for x in sequence
            push!(b.left, COMPLEMENT[Int(x)+1])
        end
    end
    return b
end

function terminal_context(b::ContigBuffer, idx::ReadIndex, c::Config, direction::Symbol)
    _, high = insert_bounds(c)
    n = max(2*idx.max_read_length, high+idx.max_read_length)
    return direction == :right ? contig_tail(b, n) : reverse_complement(contig_head(b, n))
end

function common_choice_prefix(choices::Vector{ExtensionChoice})
    isempty(choices) && return UInt8[]
    n = minimum(length(x.sequence) for x in choices)
    for i in 1:n
        b = choices[1].sequence[i]
        any(x -> x.sequence[i] != b, choices) && return choices[1].sequence[1:i-1]
    end
    return choices[1].sequence[1:n]
end

function unique_anchor(text::String, anchor::Vector{UInt8})
    UInt8('N') in anchor && return false
    positions = Set{Int}()
    for needle in unique([dna_string(anchor), dna_string(reverse_complement(anchor))])
        cursor = 1
        while cursor <= ncodeunits(text)-length(anchor)+1
            hit = findnext(needle, text, cursor)
            hit === nothing && break
            push!(positions, first(hit))
            length(positions) > 1 && return false
            cursor = first(hit)+1
        end
    end
    return length(positions) == 1
end

struct JunctionEvidence
    reads::Int   # Distinct pairs with at least one uniquely mapped spanning read.
    pairs::Int   # Distinct inward FR fragments crossing the junction in insert range.
    unique::Bool
    reason::String
end

function junction_evidence(seq::Vector{UInt8}, seam::Int, idx::ReadIndex, c::Config;
                           circular::Bool=false)
    a = c.junction_anchor
    n = length(seq)
    if circular
        n >= 4a || return JunctionEvidence(0, 0, false, "short_sequence")
        left = seq[end-a+1:end]; right = seq[1:a]
    else
        a <= seam <= n-a || return JunctionEvidence(0, 0, false, "short_flank")
        left = seq[seam-a+1:seam]; right = seq[seam+1:seam+a]
    end
    genome_text = dna_string(seq)
    if !unique_anchor(genome_text, left) || !unique_anchor(genome_text, right)
        return JunctionEvidence(0, 0, false, "repetitive_junction")
    end
    _, high = insert_bounds(c)
    flank = high+idx.max_read_length
    if circular
        flank = min(flank, n÷2-1)
        segment = vcat(seq[end-flank+1:end], seq[1:flank])
        junction = flank
    else
        low = max(1, seam-flank+1); last = min(n, seam+flank)
        segment = seq[low:last]
        junction = seam-low+1
    end
    ids, saturated = candidate_read_ids(segment, idx, c)
    saturated && return JunctionEvidence(0, 0, false, "candidate_limit")
    text = dna_string(segment)
    anchors = ExactAnchorCache(text)
    spanning = Set{UInt32}(); fragments = Set{UInt32}()
    lower, upper = insert_bounds(c)
    for pid in sort!(unique(pair_id.(ids)))
        maps1 = read_placements(text, get_read(idx.store.r1, pid), c;anchors=anchors)
        maps2 = c.paired ? read_placements(text, get_read(idx.store.r2, pid), c;anchors=anchors) : Placement[]
        for maps in (maps1, maps2)
            if length(maps) == 1
                p = only(maps)
                if p.start <= junction-a+1 && p.start+p.len-1 >= junction+a
                    push!(spanning, pid)
                end
            end
        end
        if length(maps1) == 1 && length(maps2) == 1
            p, q = only(maps1), only(maps2)
            if p.reverse != q.reverse
                forward = p.reverse ? q : p
                reverse = p.reverse ? p : q
                outer = reverse.start+reverse.len-forward.start
                if forward.start <= reverse.start && lower <= outer <= upper &&
                   forward.start <= junction-a+1 && reverse.start+reverse.len-1 >= junction+a
                    push!(fragments, pid)
                end
            end
        end
    end
    return JunctionEvidence(length(spanning), length(fragments), true, "checked")
end

struct Closure
    sequence::Vector{UInt8}
    overlap::Int
    evidence::JunctionEvidence
    kind::String            # "validated" (read evidence at the seam) or "ir_collapse"
end
Closure(sequence, overlap, evidence) = Closure(sequence, overlap, evidence, "validated")

# ---------------------------------------------------------------------------
# Closure attempt log (diagnostic only; never changes a decision).
# ---------------------------------------------------------------------------
const CLOSURE_LOG_LOCK = ReentrantLock()
const CLOSURE_LOG = Vector{Vector{String}}()
const CLOSURE_LOG_HEADER = ["path", "length", "overlap_low", "overlap_high", "start_kmer_at",
    "implied_overlap", "mismatches_at_implied", "first_mismatch_offset", "tested", "outcome"]

"Where the path's first 31-mer recurs, and how well the two ends agree at that offset."
function closure_probe(seq::Vector{UInt8})
    k = 31; n = length(seq)
    n < 2k && return (0, 0, 0, 0)
    at = 0
    for p in n-k+1:-1:2
        @inbounds seq[p] == seq[1] || continue
        @view(seq[p:p+k-1]) == @view(seq[1:k]) && (at = p; break)
    end
    at == 0 && return (0, 0, 0, 0)
    ov = n - at + 1
    mism = 0; first_mm = 0
    for i in 1:ov
        if seq[i] != seq[at+i-1]
            mism += 1; first_mm == 0 && (first_mm = i)
        end
    end
    return (at, ov, mism, first_mm)
end

function log_closure!(key, seq, low, high, tested::Vector{String}, outcome::String)
    key === nothing && return
    at, ov, mism, first_mm = closure_probe(seq)
    row = [String(key), string(length(seq)), string(low), string(high), string(at), string(ov),
           string(mism), string(first_mm), isempty(tested) ? "-" : join(tested, ";"), outcome]
    lock(CLOSURE_LOG_LOCK) do
        push!(CLOSURE_LOG, row)
    end
    return
end

"""
Closure by collapsing a terminal IR copy: end overlap of at least IR min length, identical
up to max(2, overlap/1000) indels or substitutions (repeat length of the two passes),
and the circle orders as LSC-IR-SSC-IR with an IR at least as long as the overlap (the
overlap lies in one IR copy), IR <= 30% of the circle, SSC >= 1 kb and LSC >= SSC.
"""
function ir_collapse_ok(seq::Vector{UInt8}, circle::Vector{UInt8}, overlap::Int, c::Config)
    c.assembly_type == "chl" || return false
    overlap >= c.ir_min_length || return false
    # Ends identical, allowing a few indels (repeat or homopolymer length of the two passes).
    (@view(seq[1:overlap]) == @view(seq[end-overlap+1:end]) ||
     banded_edits(@view(seq[1:overlap]), @view(seq[end-overlap+1:end])) <= max(2, overlap ÷ 1000)) || return false
    item, _ = PlastomeOrder.order_plastome(circle; min_ir=c.ir_min_length, min_sc=c.ir_min_sc_length,
                                           min_identity=c.ir_min_identity)
    item === nothing && return false
    return item.ir_pair.len >= overlap && item.ir_pair.len <= 0.3*length(circle) &&
           item.ssc_length >= 1000 && item.lsc_length >= item.ssc_length
end

"First seam s in from:to whose two flanking Junction-anchor windows are unique in `circle` (0 if none)."
function first_unique_seam(circle::Vector{UInt8}, from::Int, to::Int, a::Int)
    keys, valid = rolling_kmers(circle, a)
    counts = Dict{UInt64,Int}()
    for i in eachindex(keys)
        valid[i] && (counts[keys[i]] = get(counts, keys[i], 0) + 1)
    end
    unique_at(p) = 1 <= p <= length(keys) && valid[p] && counts[keys[p]] == 1
    for s in max(from, a):min(to, length(circle)-a)
        unique_at(s-a+1) && unique_at(s+1) && return s
    end
    return 0
end

"""
Test one closure length. A seam inside a repeat (typically an IR, when a path starts at
an IR boundary) cannot be validated where the ends meet; the identical end overlap then
has to reach unique sequence, and read evidence is taken at the first unique position
inside the overlap, which must still leave Circular min overlap of unique shared bases.
"""
function closure_attempt(seq::Vector{UInt8}, overlap::Int, idx::ReadIndex, c::Config, tested::Vector{String};
                         allow_ir_collapse::Bool=false)
    circle = seq[1:end-overlap]
    evidence = junction_evidence(circle, length(circle), idx, c; circular=true)
    if evidence.reason == "repetitive_junction"
        s = first_unique_seam(circle, c.junction_anchor, overlap - c.circular_min_overlap, c.junction_anchor)
        if s == 0
            # The whole overlap is repeat: typically the path starts and ends in IR copies, at the
            # IR/single-copy boundaries, so no read can cross the seam into unique sequence (IR
            # longer than the insert). Every IR/LSC/SSC boundary is still inside the path and was
            # built from reads; collapsing the terminal IR copy is then the only circle, if the
            # ends are identical and the circle is a plausible quadripartite plastome.
            # Only once the path can grow no further: while it still grows, its end may yet
            # leave the IR and give a seam with read evidence (validated), which comes first.
            if allow_ir_collapse && ir_collapse_ok(seq, circle, overlap, c)
                push!(tested, "$(overlap):repeat_seam_no_unique_overlap:ir_collapse:ok")
                return Closure(circle, overlap, JunctionEvidence(0, 0, false, "ir_collapse"), "ir_collapse")
            end
            push!(tested, "$(overlap):repeat_seam_no_unique_overlap:rejected")
            return nothing
        end
        shifted = junction_evidence(circle, s, idx, c)
        ok = shifted.unique && shifted.reads >= c.circular_min_reads &&
             (!c.paired || shifted.pairs >= c.min_paired_support)
        push!(tested, "$(overlap):repeat_seam_shifted_to_$(s):$(shifted.reason):reads=$(shifted.reads):pairs=$(shifted.pairs):$(ok ? "ok" : "rejected")")
        return ok ? Closure(circle, overlap, shifted) : nothing
    end
    ok = evidence.unique && evidence.reads >= c.circular_min_reads &&
         (!c.paired || evidence.pairs >= c.min_paired_support)
    push!(tested, "$(overlap):$(evidence.reason):reads=$(evidence.reads):pairs=$(evidence.pairs):$(ok ? "ok" : "rejected")")
    return ok ? Closure(circle, overlap, evidence) : nothing
end

function check_circular(seq::Vector{UInt8}, idx::ReadIndex, c::Config; log_key=nothing,
                        allow_ir_collapse::Bool=false)
    n = length(seq)
    low = max(c.circular_min_overlap, n-c.genome_range[2])
    high = min(c.circular_max_overlap, n÷2, n-c.genome_range[1])
    ext_high = min(c.circular_repeat_max_overlap, n÷2, n-c.genome_range[1])
    if low > max(high, ext_high)
        log_closure!(log_key, seq, low, high, String[], "overlap_range_empty")
        return nothing
    end
    closures = Closure[]
    tested = String[]
    candidates = low <= high ? collect(high:-1:low) : Int[]
    # Longer overlaps only where the first 31-mer recurs: cheap, and exact for a seam in an IR.
    k = 31
    if ext_high > max(high, low-1) && n >= 2k
        for ov in ext_high:-1:max(high+1, low)
            p = n - ov + 1
            @inbounds seq[p] == seq[1] || continue
            @view(seq[p:p+k-1]) == @view(seq[1:k]) && push!(candidates, ov)
        end
    end
    for overlap in candidates
        if !sequences_similar(@view(seq[1:overlap]), @view(seq[end-overlap+1:end]), c.circular_identity)
            # Gapped: the two passes through the overlap differ by short indels (e.g. repeat
            # length); tried only where the first 31-mer recurs exactly at the seam start.
            p = n - overlap + 1
            (overlap >= k && @view(seq[p:p+k-1]) == @view(seq[1:k])) || continue
            banded_edits(@view(seq[1:overlap]), @view(seq[end-overlap+1:end])) <=
                floor(Int, overlap*(1-c.circular_identity) + 1e-9) || continue
        end
        closure = closure_attempt(seq, overlap, idx, c, tested; allow_ir_collapse=allow_ir_collapse)
        closure === nothing && continue
        push!(closures, closure)
        # Multiple supported lengths are ambiguous, not a license to choose one.
        if length(closures) > 1
            log_closure!(log_key, seq, low, ext_high, tested, "ambiguous_closure_lengths")
            return nothing
        end
    end
    outcome = !isempty(closures) ? (all(x -> x.kind == "ir_collapse", closures) ? "ir_collapse" : "closed") :
              isempty(tested) ? "no_end_overlap_within_identity" : "read_evidence_rejected"
    log_closure!(log_key, seq, low, ext_high, tested, outcome)
    return isempty(closures) ? nothing : only(closures)
end

function read_seeds(path::AbstractString)
    records = Pair{String,Vector{UInt8}}[]
    header = ""; chunks = String[]; seen = Set{String}()
    for raw in eachline(path)
        line = strip(raw)
        isempty(line) && continue
        if startswith(line, ">")
            if !isempty(header)
                seq = normalized_dna(join(chunks))
                isempty(seq) && error("Empty seed '$header'")
                push!(records, header => seq)
            end
            header = String(strip(line[2:end]))
            isempty(header) && error("Empty FASTA header")
            header in seen && error("Duplicate FASTA header '$header'")
            push!(seen, header); empty!(chunks)
        else
            isempty(header) && error("FASTA sequence before first header")
            push!(chunks, String(line))
        end
    end
    if !isempty(header)
        seq = normalized_dna(join(chunks))
        isempty(seq) && error("Empty seed '$header'")
        push!(records, header => seq)
    end
    isempty(records) && error("No FASTA seeds found")
    return records
end

struct InitialContig
    id::String
    sequence::Vector{UInt8}
    source::String
end

"Correct well-supported seed substitutions without importing reference bases."
function correct_seed_read(seed::Vector{UInt8}, idx::ReadIndex, c::Config)
    c.seed_correction || return seed
    ids, saturated = candidate_read_ids(seed, idx, c)
    saturated && return seed
    text = dna_string(seed)
    ballots = [Dict{UInt32,UInt8}() for _ in seed]
    for id in ids
        read = get_read(idx.store, id)
        hits = Tuple{Int,Int,Vector{UInt8}}[]
        for rc in (false, true)
            oriented = rc ? reverse_complement(read) : read
            a = min(15, length(oriented))
            positions = Set{Int}()
            for p in unique([1, (length(oriented)-a)÷2+1, length(oriented)-a+1])
                anchor = @view oriented[p:p+a-1]
                UInt8('N') in anchor && continue
                cursor = 1
                while cursor <= length(seed)-a+1
                    hit = findnext(dna_string(anchor), text, cursor)
                    hit === nothing && break
                    push!(positions, first(hit)-p+1)
                    cursor = first(hit)+1
                end
            end
            for start in positions
                low = max(1, start); high = min(length(seed), start+length(oriented)-1)
                overlap = high-low+1
                overlap >= c.min_overlap || continue
                if sequences_similar(@view(seed[low:high]), @view(oriented[low-start+1:high-start+1]), c.min_identity)
                    push!(hits, (overlap, start, oriented))
                end
            end
        end
        isempty(hits) && continue
        sort!(hits; by=x -> (-x[1], x[2], dna_string(x[3])))
        length(hits) > 1 && hits[1][1] == hits[2][1] && continue
        _, start, oriented = hits[1]
        for p in max(1, start):min(length(seed), start+length(oriented)-1)
            b = ENCODE[Int(oriented[p-start+1])+1]
            b == 255 && continue
            pid = pair_id(id)
            ballots[p][pid] = get(ballots[p], pid, UInt8(0)) | (UInt8(1) << b)
        end
    end
    corrected = copy(seed)
    for p in eachindex(seed)
        votes = zeros(Int, 4)
        for mask in values(ballots[p])
            count_ones(mask) == 1 && (votes[trailing_zeros(mask)+1] += 1)
        end
        order = sort(collect(1:4); by=b -> (-votes[b], b))
        best, second = order[1], order[2]
        total = length(ballots[p])
        # Strict rule only: per-position majority correction of a read from a mixed library
        # can combine alleles of neighbouring variants from different haplotypes (a mosaic no
        # read contains). The seed is a real read; minority alleles are left to extension.
        if votes[best] >= c.min_coverage && votes[best]/max(1, total) >= c.consensus_fraction &&
           (votes[second] == 0 || votes[best]/votes[second] >= c.consensus_ratio)
            corrected[p] = DNA[best]
        end
    end
    return corrected
end

# Adapter starts (Illumina TruSeq, Nextera, MGI read 1/2). The data are not trimmed; only an
# initial read is cut before its adapter, since degraded (herbarium) DNA gives many
# read-through reads and an adapter tail at the start of a path blocks all extension.
const ADAPTER_PREFIXES = ("AGATCGGAAGAGC", "CTGTCTCTTATACACATCT", "AAGTCGGAGGCCAAGCGG", "AAGTCGGATCGTAGCCATG")
const MIN_INITIAL_LENGTH = 60

"Bases before a full adapter start, or before a >= 8-base adapter start at the read end."
function genomic_length(seq::Vector{UInt8})
    t = String(copy(seq)); cut = length(seq)
    for a in ADAPTER_PREFIXES
        r = findfirst(a, t)
        r === nothing || (cut = min(cut, first(r) - 1))
        for m in length(a)-1:-1:8
            if length(t) >= m && endswith(t, a[1:m])
                cut = min(cut, length(t) - m); break
            end
        end
    end
    return cut
end

const INITIAL_LOG = String[]

function find_initial_reads(c::Config, idx::ReadIndex)
    empty!(INITIAL_LOG)
    initial = InitialContig[]
    seen_sequences = Set{String}()
    for (header, seq) in read_seeds(c.seed_input)
        scores = Dict{UInt32,Int}()
        for key in all_kmer_keys(seq, c.kmer)
            haskey(idx.postings, key) || continue
            for id in idx.postings[key]
                excluded_read(UInt32(id),c) && continue
                scores[id] = get(scores, id, 0)+1
            end
        end
        ranked = sort!(collect(scores); by=x -> (-last(x), first(x)))
        added = 0
        for (id, score) in ranked
            score < c.min_seed_hits && break
            (added >= c.seeds_per_record || length(initial) >= c.max_initial_contigs) && break
            raw = get_read(idx.store, id)
            g = genomic_length(raw)
            g >= MIN_INITIAL_LENGTH || continue                 # mostly adapter: next candidate
            read = correct_seed_read(raw[1:g], idx, c)
            canonical = min(dna_string(read), dna_string(reverse_complement(read)))
            canonical in seen_sequences && continue
            push!(seen_sequences, canonical)
            push!(initial, InitialContig("seed_$(length(initial)+1)", read, header))
            push!(INITIAL_LOG, "seed_$(length(initial))\tread $(id)\tshared minimizers $(score)\tlength $(g) of $(length(raw))\t$(dna_string(read))")
            added += 1
        end
        added == 0 && @warn "No new initial read for FASTA record" header
        length(initial) >= c.max_initial_contigs && break
    end
    isempty(initial) && error("No initial reads: inspect seed, Min seed hits, k-mer size and minimizer filtering")
    return initial
end

"A resolved bubble in output contig coordinates (alleles on the contig's forward strand)."
struct ContigVariant
    position::Int
    major::String
    minor::String
    major_support::Int
    minor_support::Int
    informative::Int
    switch::Bool
end
ContigVariant(position, major, minor, ms, mn, informative) = ContigVariant(position, major, minor, ms, mn, informative, false)

struct AssemblyResult
    id::String
    sequence::Vector{UInt8}
    circular::Bool
    reason::String
    source::String
    closure_overlap::Int
    junction_reads::Int
    junction_pairs::Int
    variants::Vector{ContigVariant}
end
AssemblyResult(id, sequence, circular, reason, source, overlap, reads, pairs) =
    AssemblyResult(id, sequence, circular, reason, source, overlap, reads, pairs, ContigVariant[])

const PathVariant = Tuple{Symbol,Int,BubbleVariant}   # direction, side length before the extension

"Map extension-local bubbles to final contig coordinates (left flank is stored reversed, complemented)."
function contig_variants(buffer::ContigBuffer, path_variants::Vector{PathVariant}, circle_length::Int)
    L = length(buffer.left)
    out = ContigVariant[]
    for (direction, side_length, v) in path_variants
        span = max(length(v.major), 1)
        if direction == :right
            pos = L + side_length + v.offset
            major = v.major; minor = v.minor
        else
            first_index = side_length + v.offset          # index in buffer.left of the allele start
            pos = L - (first_index + span - 1) + 1
            major = reverse_complement(v.major); minor = reverse_complement(v.minor)
        end
        circle_length > 0 && pos > circle_length && (pos -= circle_length)
        push!(out, ContigVariant(pos, String(copy(major)), String(copy(minor)),
                                 v.major_support, v.minor_support, v.informative, v.switch))
    end
    sort!(out; by=x -> (x.position, x.major, x.minor))
    return unique(x -> (x.position, x.major, x.minor), out)
end

function result_from_buffer(id, buffer, reason, source, closure; path_variants=PathVariant[])
    if closure === nothing
        return AssemblyResult(id, contig_sequence(buffer), false, reason, source, 0, 0, 0,
                              contig_variants(buffer, path_variants, 0))
    end
    return AssemblyResult(id, closure.sequence, true,
                          closure.kind == "ir_collapse" ? "closed_ir_collapse" : "validated_closure", source,
                          closure.overlap, closure.evidence.reads, closure.evidence.pairs,
                          contig_variants(buffer, path_variants, length(closure.sequence)))
end

mutable struct PathState
    buffer::ContigBuffer
    iterations::Int
    variants::Vector{PathVariant}
    supports::Dict{Symbol,Vector{Int}}    # per direction: min support of normal iterations
    high_start::Dict{Symbol,Int}          # side length where the current high-support run began
    high_run::Dict{Symbol,Int}            # bases added in that run
    blocked::Set{Symbol}                  # directions stopped by coverage_jump
    jumps::Int                            # mate jumps made on this path
end
PathState(buffer::ContigBuffer, iterations::Int) = PathState(buffer, iterations, PathVariant[])
PathState(buffer::ContigBuffer, iterations::Int, variants::Vector{PathVariant};
          supports=Dict(:left => Int[], :right => Int[])) =
    PathState(buffer, iterations, variants, supports, Dict(:left => 0, :right => 0),
              Dict(:left => 0, :right => 0), Set{Symbol}(), 0)

const COVERAGE_TRACE_LOCK = ReentrantLock()
const COVERAGE_TRACE = Vector{String}()

"Forward code of a k-mer (2 bits per base), or typemax(UInt64) if it contains a non-ACGT base."
function kmer_code(seq::AbstractVector{UInt8})
    v = UInt64(0)
    @inbounds for b in seq
        e = ENCODE[Int(b) + 1]
        e == 255 && return typemax(UInt64)
        v = (v << 2) | UInt64(e)
    end
    return v
end

"Forward k-mer codes of a sequence."
function forward_kmer_set(seq::AbstractVector{UInt8}, k::Int)
    out = Set{UInt64}(); sizehint!(out, max(length(seq) - k + 1, 0))
    for i in 1:length(seq)-k+1
        v = kmer_code(@view seq[i:i+k-1])
        v == typemax(UInt64) || push!(out, v)
    end
    return out
end

"""
Single-path consensus extension of `seq` to the right by up to `maxlen` bases (no forks, the
best-supported choice; phasing off because `seq` may end in a raw read). Local helper for
mate jumps; the main assembly loop is not used.
"""
function local_extend(seq::Vector{UInt8}, idx::ReadIndex, c::Config; maxlen::Int=600)
    out = copy(seq); target = length(seq) + maxlen
    ctxlen = 2*idx.max_read_length
    while length(out) < target
        ctx = out[max(1, end-ctxlen+1):end]
        candidates, status = extension_candidates(ctx, idx, c; phase=0)
        status == "ok" || break
        choices = consensus_choices(candidates, c)
        isempty(choices) && break
        append!(out, choices[1].sequence)
    end
    return out
end

"Mate jumps are allowed only in the fallback pass (see run_assembly)."
const MATE_JUMPS_ON = Ref(false)

const MATE_JUMP_LOG_LOCK = ReentrantLock()
const MATE_JUMP_LOG = String[]

"""
Mate jump over a stretch no read crosses (e.g. a long poly-G the sequencer does not read).
Reads placed wholly and uniquely in the last insert-length of the context, pointing to its
end, have mates beyond it. Their mates (as genome-forward sequence, not placed in the
context) are grouped by shared 31-mers; the best-connected mate must be linked to
>= Min coverage distinct pairs and be new sequence (< 20% of its 31-mers in the path on
either strand). It is corrected like a seed, extended back towards the obstacle and forward
by local consensus, and appended after Unknown gap Ns: the gap length is not known (the insert
distribution only estimates a distance), so the Ns stay in the sequence and are reported as
an internal gap even if the genome later closes with read evidence elsewhere. Used only in a
fallback pass when the normal assembly gave no gap-free circle, so libraries that close
without it are never changed by it.
"""
function mate_jump!(state, context::Vector{UInt8}, built::Int, direction::Symbol, idx::ReadIndex, c::Config, path::String)
    (MATE_JUMPS_ON[] && c.paired && c.mate_jumps > state.jumps && c.unknown_gap_ns > 0) || return false
    k = c.kmer
    _, high = insert_bounds(c)
    L = length(context)
    W = min(L, high)
    W >= 2*idx.max_read_length || return false
    window = context[L-W+1:L]
    text = dna_string(window); anchors = ExactAnchorCache(text)
    ids, saturated = candidate_read_ids(window, idx, c)
    saturated && return false
    mates = Dict{UInt32,Vector{UInt8}}()          # pair -> mate sequence, genome-forward orientation
    for id in ids
        read = get_read(idx.store, id)
        pl = read_placements(text, read, c; anchors=anchors)
        (length(pl) == 1 && !only(pl).reverse) || continue
        q = only(pl)
        q.start + length(read) - 1 + 40 <= W || continue      # read lies before the end region
        mate = get_read(idx.store, id ⊻ UInt32(1))
        isempty(read_placements(text, mate, c; anchors=anchors)) || continue   # mate not in context
        mates[pair_id(id)] = reverse_complement(mate)
    end
    length(mates) >= c.min_coverage || return false
    # Best-connected mate: shares 31-mers with the most other pairs' mates.
    sets = Dict(p => forward_kmer_set(m, k) for (p, m) in mates)
    best = UInt32(0); links = 0
    for (p, sp) in sets
        n = count(q -> q != p && !isdisjoint(sp, sets[q]), keys(sets))
        n > links && (best = p; links = n)
    end
    links + 1 >= c.min_coverage || return false
    far = correct_seed_read(mates[best], idx, c)
    length(unique(far)) >= 3 && longest_run(far) <= 12 || return false
    # New sequence only: not already in the path on either strand.
    pathset = forward_kmer_set(context, k)
    kf = forward_kmer_set(far, k); kr = forward_kmer_set(reverse_complement(far), k)
    isempty(kf) && return false
    known = count(v -> v in pathset, kf) + count(v -> v in pathset, kr)
    known <= 0.2*length(kf) || return false
    # Back towards the obstacle, then forward past the raw read end.
    far = reverse_complement(local_extend(reverse_complement(far), idx, c; maxlen=max(high, 300)))
    far = local_extend(far, idx, c; maxlen=300)
    # Both sides are cut back to their unique flank before a low-complexity run next to the
    # gap (homopolymer, microsatellite): its length and the bases just past it are what no read
    # resolves, so all of that uncertainty is left to the Ns instead of partial, duplicated runs.
    ft = head_run_trim(far)
    ft > 0 && length(far) - ft >= 50 && (far = far[ft+1:end])
    nt = tail_run_trim(context)
    if 0 < nt <= built
        trim_extension!(state.buffer, nt, direction)
        keep = direction == :right ? length(state.buffer.right) : length(state.buffer.left)
        filter!(v -> !(v[1] == direction && v[2] + v[3].offset > keep), state.variants)
    else
        nt = 0
    end
    append_extension!(state.buffer, vcat(fill(UInt8('N'), c.unknown_gap_ns), far), direction)
    state.jumps += 1
    row = join((path, direction, length(state.buffer), length(mates), links + 1, length(far), nt, ft), '\t')
    lock(() -> push!(MATE_JUMP_LOG, row), MATE_JUMP_LOG_LOCK)
    return true
end

"""
Bases to cut from the end of `seq` so that it ends before a low-complexity run (homopolymer
or tandem repeat, period 1-6, as periodic_suffix) that ends within the last `slack` bases;
0 if there is none.
"""
function tail_run_trim(seq::AbstractVector{UInt8}; slack::Int=30)
    L = length(seq)
    for kk in 0:slack
        L - kk > 0 || break
        p, r = periodic_suffix(@view seq[1:L-kk])
        p > 0 && return kk + r
    end
    return 0
end

"Same at the start of `seq` (a reversed sequence keeps its periodicity)."
head_run_trim(seq::AbstractVector{UInt8}; slack::Int=40) = tail_run_trim(reverse(seq); slack=slack)

"""
Runs of N (>= 5) in a sequence: gaps of unknown length (mate jumps, N closures).
Returns vector of (start, length).
"""
function n_gaps(seq::AbstractVector{UInt8})
    out = Tuple{Int,Int}[]; i = 1; n = length(seq)
    while i <= n
        if seq[i] == UInt8('N')
            j = i
            while j < n && seq[j+1] == UInt8('N'); j += 1; end
            j - i + 1 >= 5 && push!(out, (i, j - i + 1))
            i = j + 1
        else
            i += 1
        end
    end
    return out
end

"""
Long inverted self-matches of a sequence: segments A and B (>= min_len) with B the reverse
complement of A, from unique 31-mers (anti-diagonal chains). A plastome path shows one such
pair, its IR; a path that turned back at an IR exit shows a much longer one (IR plus the
single-copy region walked backwards). Diagnostic only.
"""
function self_inversions(seq::Vector{UInt8}, k::Int; min_len::Int=1000)
    n = length(seq)
    fwd = Dict{UInt64,Int}(); dup = Set{UInt64}()
    for i in 1:n-k+1
        v = kmer_code(@view seq[i:i+k-1]); v == typemax(UInt64) && continue
        haskey(fwd, v) ? push!(dup, v) : (fwd[v] = i)
    end
    rcs = reverse_complement(seq)
    blocks = NTuple{4,Int}[]          # (i_first, i_last, j_first, j_last), i decreasing as j increases
    cur = nothing
    for j in 1:n-k+1
        v = kmer_code(@view rcs[n-j-k+2:n-j+1])
        i = (v == typemax(UInt64) || v in dup) ? 0 : get(fwd, v, 0)
        if i > 0 && i < j && cur !== nothing && j == cur[4] + 1 && i == cur[2] - 1
            cur = (cur[1], i, cur[3], j)
        else
            cur !== nothing && cur[4] - cur[3] + k >= min_len && push!(blocks, cur)
            cur = (i > 0 && i < j) ? (i, i, j, j) : nothing
        end
    end
    cur !== nothing && cur[4] - cur[3] + k >= min_len && push!(blocks, cur)
    # A = seq[i_last : i_first+k-1], B = seq[j_first : j_last+k-1]
    return [(b[2], b[1] + k - 1, b[3], b[4] + k - 1) for b in blocks]
end

"""
Extension stopped with the end inside a tandem repeat: try resolve_tandem_end and, if it
succeeds, trim the overshoot, append the rest of the repeat plus the exit flank, and log it.
"""
function tandem_step!(state, context::Vector{UInt8}, built::Int, direction::Symbol,
                      idx::ReadIndex, c::Config, path::String)
    r = resolve_tandem_end(context, idx, c; phase=clamp(built,0,c.phase_window), genome=contig_sequence(state.buffer))
    r === nothing && return false
    r[1] <= built || return false                 # never trim into the seed read
    tandem_apply!(state, context, r, direction, path)
    return true
end

"Apply a resolve_tandem_end result: trim the overshoot, append repeat rest + exit flank, log."
function tandem_apply!(state, context::Vector{UInt8}, r, direction::Symbol, path::String)
    trim, append, best, second, rlen = r
    if trim > 0
        trim_extension!(state.buffer, trim, direction)
        keep = direction == :right ? length(state.buffer.right) : length(state.buffer.left)
        filter!(v -> !(v[1] == direction && v[2] + v[3].offset > keep), state.variants)
    end
    append_extension!(state.buffer, append, direction)
    row = join((path, direction, length(state.buffer), periodic_suffix(context)[1],
                periodic_suffix(context)[2], rlen, trim, length(append), best, second), '\t')
    lock(() -> push!(TANDEM_LOG, row), TANDEM_LOG_LOCK)
    return nothing
end

"""
Choices diverge at the end of a homopolymer or short tandem run: the sequence up to the
divergence (last 30 context bases + common prefix of the choices) ends in a periodic run of
period 1-6 covering >= 6 bases and >= 2 motif copies.
"""
function repeat_length_fork(context::Vector{UInt8}, choices::Vector{ExtensionChoice})
    common = common_choice_prefix(choices)
    window = vcat(context[max(1, end-29):end], common)
    period, run = periodic_suffix(window)
    return period > 0 && run >= max(6, 2period)
end

"Pairs voting at the first extension position: the path's support on the scale votes are compared on."
function first_position_support(candidates::Vector{ExtensionCandidate}, work::VoteWorkspace)
    _, _, masks = extension_votes!(work, candidates, collect(eachindex(candidates)), 1)
    return length(masks)
end

"Median first-position support of the last 50 normal iterations in this direction (0 while < 10 known)."
function expected_support(state::PathState, direction::Symbol)
    h = state.supports[direction]
    length(h) < 10 && return 0
    return round(Int, median(@view h[max(1, end-49):end]))
end

function assemble_seed(seed::InitialContig, idx::ReadIndex, c::Config; monitor=nothing, length_limit::Int=c.genome_range[2]+c.circular_repeat_max_overlap)
    pending = [PathState(ContigBuffer(seed.sequence), 0)]
    results = AssemblyResult[]
    spawned = 1
    progress_seed!(monitor,seed.id,length(seed.sequence),0,spawned)
    support_votes = VoteWorkspace()
    while !isempty(pending)
        state = popfirst!(pending)
        last_checked = 0
        closure = nothing
        # IR-collapse circle seen while the path still grew: kept with the buffer and variants of
        # that moment, used only if no validated closure follows (the end may leave the IR by the
        # other exit, after which the ends never match again).
        collapse = nothing
        reason = "max_iterations"
        while state.iterations < c.max_iterations
            state.iterations += 1
            progress = false
            stops = String[]
            for direction in (:right, :left)
                if direction in state.blocked
                    push!(stops, "$(direction):coverage_jump"); continue
                end
                context = terminal_context(state.buffer, idx, c, direction)
                # Phasing applies only to bases built by extension: the seed is one raw read, and a
                # sequencing error near its end would otherwise match no other read exactly.
                built = direction == :right ? length(state.buffer.right) - length(seed.sequence) :
                                              length(state.buffer.left)
                candidates, status = extension_candidates(context, idx, c;
                                                          phase=clamp(built, 0, c.phase_window))
                if status != "ok"
                    if tandem_step!(state, context, built, direction, idx, c, "$(seed.id)_path_$(length(results)+1)")
                        progress = true; continue
                    end
                    if status == "no_extension" &&
                       mate_jump!(state, context, built, direction, idx, c, "$(seed.id)_path_$(length(results)+1)")
                        progress = true; continue
                    end
                    push!(stops, "$(direction):$status"); continue
                end
                track = c.coverage_switch_ratio > 0 || c.coverage_trace
                expected = track ? expected_support(state, direction) : 0
                here = track ? first_position_support(candidates, support_votes) : 0
                choices = consensus_choices(candidates, c; expected=expected)
                if isempty(choices)
                    if tandem_step!(state, context, built, direction, idx, c, "$(seed.id)_path_$(length(results)+1)")
                        progress = true; continue
                    end
                    if mate_jump!(state, context, built, direction, idx, c, "$(seed.id)_path_$(length(results)+1)")
                        progress = true; continue
                    end
                    push!(stops, "$(direction):low_support_or_ambiguity"); continue
                end
                # A fork that opens at the end of a homopolymer or short tandem run is a length
                # variant of that run (or sequencing stutter), not a second genome structure: it is
                # resolved from reads spanning the whole run (distinct pairs, entry flank -> exit
                # flank); if they do not decide, the best-supported choice is taken without
                # spawning a branch, so the branch budget stays for real forks (IR exits etc.).
                if length(choices) > 1 && repeat_length_fork(context, choices)
                    nchoices = length(choices)
                    supports = join((x.min_support for x in choices), ",")
                    common = common_choice_prefix(choices)
                    ctx2 = vcat(context, common)
                    r = resolve_tandem_end(ctx2, idx, c; phase=clamp(built+length(common),0,c.phase_window), genome=contig_sequence(state.buffer))
                    label = "$(seed.id)_path_$(length(results)+1)"
                    if r !== nothing && r[1] <= built + length(common)
                        append_extension!(state.buffer, common, direction)
                        tandem_apply!(state, ctx2, r, direction, label)
                        action = "repeat_resolved"
                        choices = ExtensionChoice[]
                    else
                        action = "repeat_majority"
                        choices = choices[1:1]
                    end
                    row = join((label, direction, length(state.buffer), nchoices, action, spawned, supports, ""), '\t')
                    lock(() -> push!(FORK_LOG, row), FORK_LOG_LOCK)
                    if isempty(choices)
                        progress = true; continue
                    end
                end
                if length(choices) > 1
                    budget_ok = spawned+length(choices)-1 <= c.max_branch_paths
                    row = join(("$(seed.id)_path_$(length(results)+1)", direction, length(state.buffer),
                                length(choices), budget_ok ? "branched" : "branch_limit", spawned,
                                join((x.min_support for x in choices), ","),
                                join((String(copy(x.sequence[1:min(40, length(x.sequence))])) for x in choices), ",")), '\t')
                    lock(() -> push!(FORK_LOG, row), FORK_LOG_LOCK)
                    if budget_ok
                        for choice in choices[2:end]
                            alternate = copy(state.buffer)
                            side = direction == :right ? length(alternate.right) : length(alternate.left)
                            append_extension!(alternate, choice.sequence, direction)
                            alt_variants = copy(state.variants)
                            for v in choice.variants
                                push!(alt_variants, (direction, side, v))
                            end
                            branch = PathState(alternate, state.iterations, alt_variants;
                                               supports=deepcopy(state.supports))
                            # A branch inherits the coverage-control state, including stopped directions.
                            branch.high_start = copy(state.high_start); branch.high_run = copy(state.high_run)
                            branch.blocked = copy(state.blocked)
                            branch.jumps = state.jumps
                            push!(pending, branch)
                            spawned += 1
                        end
                    else
                        common = common_choice_prefix(choices)
                        if isempty(common)
                            push!(stops, "$(direction):branch_limit"); continue
                        end
                        kept = [v for v in choices[1].variants if v.offset + max(length(v.major), 1) - 1 <= length(common)]
                        choices = [ExtensionChoice(common, minimum(x.min_support for x in choices), "branch_limit", kept)]
                    end
                end
                chosen = choices[1]
                side = direction == :right ? length(state.buffer.right) : length(state.buffer.left)
                append_extension!(state.buffer, chosen.sequence, direction)
                for v in chosen.variants
                    push!(state.variants, (direction, side, v))
                end
                progress = true
                if c.coverage_trace
                    switches = count(v -> v.switch, chosen.variants)
                    row = join(("$(seed.id)_path_$(length(results)+1)", direction, length(state.buffer), here, expected,
                                chosen.min_support, state.high_run[direction], switches, length(candidates)), '\t')
                    lock(() -> push!(COVERAGE_TRACE, row), COVERAGE_TRACE_LOCK)
                end
                # Sustained support far above the path's own: the direction has left the
                # plastome (e.g. into a contaminant genome). Trim back to where it began, stop.
                if expected > 0 && c.coverage_switch_ratio > 0 && here > c.coverage_switch_ratio*expected
                    state.high_run[direction] == 0 && (state.high_start[direction] = side)
                    state.high_run[direction] += length(chosen.sequence)
                    if state.high_run[direction] >= c.coverage_jump_length
                        keep = state.high_start[direction]
                        resize!(direction == :right ? state.buffer.right : state.buffer.left, keep)
                        filter!(v -> !(v[1] == direction && v[2] + v[3].offset > keep), state.variants)
                        push!(state.blocked, direction)
                        push!(stops, "$(direction):coverage_jump")
                    end
                else
                    state.high_run[direction] = 0
                    track && push!(state.supports[direction], here)
                end
                length(state.buffer) >= length_limit && break
            end
            n = length(state.buffer)
            progress_seed!(monitor,seed.id,n,length(results),spawned)
            if n >= c.genome_range[1]+c.circular_min_overlap &&
               (n-last_checked >= c.circular_check_every || !progress || n >= length_limit)
                closure = check_circular(contig_sequence(state.buffer), idx, c;
                                         log_key="$(seed.id)_path_$(length(results)+1)",
                                         allow_ir_collapse=true)
                last_checked = n
                if closure !== nothing && closure.kind == "ir_collapse" && progress && n < length_limit
                    collapse === nothing && (collapse = (closure, copy(state.buffer), copy(state.variants)))
                    closure = nothing
                end
                if closure !== nothing
                    reason = "validated_closure"; break
                end
            end
            if n >= length_limit
                reason = "length_limit_unresolved"; break
            elseif !progress
                reason = join(stops, ";"); break
            end
        end
        # Also check the final state after the iteration limit.
        if closure === nothing && length(state.buffer) >= c.genome_range[1]+c.circular_min_overlap
            closure = check_circular(contig_sequence(state.buffer), idx, c;
                                     log_key="$(seed.id)_path_$(length(results)+1)")
        end
        id = "$(seed.id)_path_$(length(results)+1)"
        if closure === nothing && collapse !== nothing
            closure, buffer_then, variants_then = collapse
            push!(results, result_from_buffer(id, buffer_then, reason, seed.source, closure;
                                              path_variants=variants_then))
        else
            push!(results, result_from_buffer(id, state.buffer, reason, seed.source, closure;
                                              path_variants=state.variants))
        end
        progress_seed!(monitor,seed.id,length(state.buffer),length(results),spawned;
                       closed=closure!==nothing)
    end
    progress_seed!(monitor,seed.id,maximum(length(x.sequence) for x in results),length(results),spawned;complete=true)
    return results
end

"Rank of closure evidence (lower = stronger); linear results last."
closure_rank(x::AssemblyResult) = something(findfirst(==(closure_status(x)),
    ("validated", "ir_collapse", "read_bridge", "graph_path", "estimated_repeat", "gap_unknown", "inferred_join", "topology_inferred", "unresolved", "ref_scaffold")), 11)

function deduplicate(results::Vector{AssemblyResult})
    # Among identical sequences the one with the strongest closure evidence is kept: two paths
    # may give the same circle, one validated by reads at the seam and one by IR collapse.
    ordered = sort(results; by=x -> (!x.circular, closure_rank(x), -length(x.sequence), x.id))
    kept = AssemblyResult[]
    for item in ordered
        text = dna_string(item.sequence)
        rc = dna_string(reverse_complement(item.sequence))
        duplicate = false
        for old in kept
            length(item.sequence) <= length(old.sequence) || continue
            oldtext = dna_string(old.sequence)
            # A circle may contain a fragment crossing its arbitrary FASTA origin.
            haystack = old.circular ? oldtext*oldtext : oldtext
            if occursin(text, haystack) || occursin(rc, haystack)
                duplicate = true; break
            end
        end
        duplicate || push!(kept, item)
    end
    return kept
end

struct ContigLink
    first::Int
    second::Int
    reverse_first::Bool
    reverse_second::Bool
    overlap::Int
    evidence::JunctionEvidence
end

function find_contig_links(contigs::Vector{AssemblyResult}, idx::ReadIndex, c::Config)
    links = ContigLink[]
    for i in eachindex(contigs), j in i+1:length(contigs)
        (contigs[i].circular || contigs[j].circular) && continue
        for ri in (false, true), rj in (false, true)
            first = ri ? reverse_complement(contigs[i].sequence) : contigs[i].sequence
            second = rj ? reverse_complement(contigs[j].sequence) : contigs[j].sequence
            high = min(length(first)-1, length(second)-1, c.scaffold_max_overlap)
            for overlap in high:-1:c.scaffold_min_overlap
                sequences_similar(@view(first[end-overlap+1:end]), @view(second[1:overlap]), c.scaffold_identity) || continue
                merged = vcat(first, @view(second[overlap+1:end]))
                length(merged) <= c.genome_range[2]+c.circular_max_overlap || break
                evidence = junction_evidence(merged, length(first), idx, c)
                if evidence.unique && evidence.reads >= c.min_coverage &&
                   (!c.paired || evidence.pairs >= c.min_paired_support)
                    push!(links, ContigLink(i, j, ri, rj, overlap, evidence))
                end
                break
            end
        end
    end
    return links
end

link_endpoints(link::ContigLink) = ((link.first, link.reverse_first ? :left : :right),
                                  (link.second, link.reverse_second ? :right : :left))

function next_scaffold_id!(occupied::Set{String}, counter::Base.RefValue{Int}, round::Int)
    while true
        id = "scaffold_r$(round)_$(counter[])"
        counter[] += 1
        if !(id in occupied)
            push!(occupied, id)
            return id
        end
    end
end

function scaffold_contigs(contigs::Vector{AssemblyResult}, idx::ReadIndex, c::Config; quiet::Bool=false, monitor=nothing)
    c.scaffold || return contigs
    occupied = Set(x.id for x in contigs)
    for round in 1:c.scaffold_rounds
        progress_update!(monitor;detail="round $round/$(c.scaffold_rounds); $(length(contigs)) contigs")
        links = find_contig_links(contigs, idx, c)
        isempty(links) && break
        degree = Dict{Tuple{Int,Symbol},Int}()
        for link in links, endpoint in link_endpoints(link)
            degree[endpoint] = get(degree, endpoint, 0)+1
        end
        # Only mutually unambiguous ends; other links remain unresolved.
        filter!(link -> all(e -> degree[e] == 1, link_endpoints(link)), links)
        isempty(links) && break
        sort!(links; by=x -> (-x.evidence.pairs, -x.evidence.reads, -x.overlap,
                             x.first, x.second, x.reverse_first, x.reverse_second))
        used = Set{Int}(); next = AssemblyResult[]; counter = Ref(1)
        for link in links
            (link.first in used || link.second in used) && continue
            a, b = contigs[link.first], contigs[link.second]
            first = link.reverse_first ? reverse_complement(a.sequence) : a.sequence
            second = link.reverse_second ? reverse_complement(b.sequence) : b.sequence
            seq = vcat(first, @view(second[link.overlap+1:end]))
            id = next_scaffold_id!(occupied, counter, round)
            closure = check_circular(seq, idx, c; log_key=id)
            result = result_from_buffer(id, ContigBuffer(seq), "supported_overlap_merge",
                                        a.source*";"*b.source, closure)
            push!(next, result)
            push!(used, link.first); push!(used, link.second)
        end
        isempty(used) && break
        for i in eachindex(contigs)
            i in used || push!(next, contigs[i])
        end
        contigs = deduplicate(next)
        !quiet && @printf("Scaffolding round %d: %d contigs\n", round, length(contigs))
    end
    return contigs
end

# ---------------------------------------------------------------------------
# Closure across a short tandem repeat (e.g. (TA)n) that no read spans.
# Both path ends stop in the same repeat; read pairs with one mate in unique
# sequence before the right end and the other after the left end span the gap.
# Gap = median library insert - (mate distances to the two ends); the repeat is
# continued in phase for that length. The result is circular but explicitly
# flagged (stop_reason closed_estimated_repeat, _estimated_gaps.tsv).
# ---------------------------------------------------------------------------
struct GapClosure
    status::String            # closed, or the reason the attempt was abandoned
    contig::String
    junction::Int             # first filled base in the circle (0 when bases were trimmed)
    motif::String
    change::Int               # bases added (> 0) or removed (< 0) at the junction
    gap_median::Float64
    gap_min::Int
    gap_max::Int
    spanning_pairs::Int
    median_insert::Float64
    right_overhang::String    # read consensus beyond the right end (>= 1 read), for inspection
    left_overhang::String     # same beyond the left start, forward strand
end
const GAP_CLOSURES_LOCK = ReentrantLock()
const GAP_CLOSURES = GapClosure[]

"Longest p-periodic suffix (p = 1..6, >= 8 bases and >= 3 units): (period, length) or (0, 0)."
function periodic_suffix(seq::AbstractVector{UInt8})
    n = length(seq); best = (0, 0)
    for p in 1:6
        n > p || continue
        r = p
        while n - r >= 1 && seq[n-r] == seq[n-r+p]
            r += 1
        end
        (r >= 8 && r >= 3p) || continue
        (best[1] == 0 || r > best[2]) && (best = (p, r))
    end
    return best
end
periodic_prefix(seq::AbstractVector{UInt8}) = periodic_suffix(reverse(seq))

"Unique FR placements of both mates: pair id => (forward placement, reverse placement)."
function pair_placements(text::String, ids, idx::ReadIndex, c::Config)
    anchors = ExactAnchorCache(text)
    out = Dict{UInt32,Tuple{Vector{Placement},Vector{Placement}}}()
    for pid in sort!(unique(pair_id.(ids)))
        excluded_read(read_id(pid,true),c) && continue
        out[pid] = (read_placements(text, get_read(idx.store.r1, pid), c; anchors=anchors),
                    read_placements(text, get_read(idx.store.r2, pid), c; anchors=anchors))
    end
    return out
end

"Median outer fragment length from FR pairs inside a central segment of the contig."
function contig_insert_median(seq::Vector{UInt8}, idx::ReadIndex, c::Config)
    span = min(length(seq), 6000)
    start = (length(seq) - span) ÷ 2 + 1
    segment = seq[start:start+span-1]
    ids, saturated = candidate_read_ids(segment, idx, c)
    saturated && return NaN
    outers = Int[]
    for (_, (p1, p2)) in pair_placements(dna_string(segment), ids, idx, c)
        (length(p1) == 1 && length(p2) == 1) || continue
        p, q = only(p1), only(p2)
        p.reverse != q.reverse || continue
        f = p.reverse ? q : p; r = p.reverse ? p : q
        f.start <= r.start && push!(outers, r.start + r.len - f.start)
    end
    return length(outers) < 30 ? NaN : median(outers)
end

"""
Per spanning pair: bases of the forward mate inside the contig before its right end plus
bases of the reverse mate inside the contig after its left start. Mates may lie wholly in
the flanks (ungapped placement) or overhang into the gap (terminal overlap >= Min overlap);
with inserts near 2x read length only overhanging mates can span a gap at all.
"""
function end_spanning_sums(seq::Vector{UInt8}, idx::ReadIndex, c::Config)
    c.paired || return Int[]
    W = min(length(seq)÷3,insert_bounds(c)[2]+idx.max_read_length)
    right_d = Dict{UInt32,Int}(); left_d = Dict{UInt32,Int}()
    genome_text = dna_string(seq); genome_anchors = ExactAnchorCache(genome_text)
    for (side,window) in ((:right,seq[end-W+1:end]),(:left,seq[1:W]))
        ids,saturated = candidate_read_ids(window,idx,c)
        saturated && continue
        for id in ids
            # A unique placement in a terminal window alone is not locus-specific.
            pl = read_placements(genome_text,get_read(idx.store,id),c;anchors=genome_anchors)
            length(pl)==1 || continue
            q=only(pl)
            if side==:right && !q.reverse && q.start>=length(seq)-W+1
                right_d[id]=length(seq)-q.start+1
            elseif side==:left && q.reverse && q.start+q.len-1<=W
                left_d[id]=q.start+q.len-1
            end
        end
    end
    for (d,context) in ((right_d,seq[end-W+1:end]),(left_d,reverse_complement(seq[1:W])))
        candidates,status=extension_candidates(context,idx,c)
        status=="ok" || continue
        for x in candidates
            (x.reverse || haskey(d,x.id)) && continue
            # Partial placement must include a unique, complex anchor before the repeat.
            read=get_read(idx.store,x.id)
            a=min(25,x.overlap,length(read))
            a>=20 || continue
            anchor=read[1:a]
            length(unique(anchor))>=3 && longest_run(anchor)<=6 || continue
            @view(context[end-x.overlap+1:end-x.overlap+a])==anchor || continue
            unique_anchor(genome_text,anchor) || continue
            d[x.id]=x.overlap
        end
    end
    sums=Dict{UInt32,Set{Int}}()
    for id in sort!(collect(keys(right_d)))
        mate=id ⊻ UInt32(1)
        haskey(left_d,mate) || continue
        push!(get!(Set{Int},sums,pair_id(id)),right_d[id]+left_d[mate])
    end
    return [only(sums[p]) for p in sort!(collect(keys(sums))) if length(sums[p])==1]
end

"Column consensus of read overhangs beyond a context end (>= min_depth reads, >= 2/3 agreement)."
function overhang_consensus(candidates::Vector{ExtensionCandidate}, min_depth::Int)
    out=UInt8[]; active=collect(eachindex(candidates)); work=VoteWorkspace()
    while !isempty(active)
        pos=length(out)+1
        cnt,_,masks=extension_votes!(work,candidates,active,pos)
        best=argmax(cnt)
        cnt[best]>=min_depth && cnt[best]>=2length(masks)/3 || break
        push!(out,DNA[best])
        active=compatible_active(candidates,active,masks,pos,best)
    end
    return out
end

"Last `k` path bases + read overhang; its periodic suffix (if any) is the repeat the end ran into."
struct EndModel
    x::Vector{UInt8}          # path tail (k bases) followed by the overhang consensus
    k::Int
    period::Int               # 0: no repeat seen from this end
    run_start::Int            # first index of the periodic suffix in x (period > 0)
end

function end_model(tail::Vector{UInt8}, candidates::Vector{ExtensionCandidate}, min_depth::Int)
    x = vcat(tail, overhang_consensus(candidates, min_depth))
    p, r = periodic_suffix(x)
    return EndModel(x, length(tail), p, length(x) - r + 1)
end

function end_candidates(context::Vector{UInt8}, idx::ReadIndex, c::Config)
    candidates, status = extension_candidates(context, idx, c)
    return status == "ok" ? candidates : ExtensionCandidate[]
end

"""
First tandem run (period 1-6, >= 8 bases, >= 3 units) starting at or after `from` in x,
as (period, start) or (0, 0). Unlike periodic_suffix, the run need not reach the end of
x: single-read overhangs often leave the repeat or carry errors at their tails.
"""
function first_tandem_run(x::Vector{UInt8}, from::Int)
    for q in max(1, from):length(x), p in 1:6
        q + p <= length(x) || continue
        r = p
        while q + r <= length(x) && x[q+r] == x[q+r-p]
            r += 1
        end
        (r >= 8 && r >= 3p) && return (p, q)
    end
    return (0, 0)
end

"Rotation-invariant key of the repeat an end runs into (forward strand), \"\" if none."
function end_repeat_key(m::EndModel, forward::Bool)
    p, q = first_tandem_run(m.x, m.k - 5)
    p == 0 && return ""
    u = m.x[q:q+p-1]
    forward || (u = reverse_complement(u))
    return minimum(String(copy(circshift(u, i))) for i in 0:p-1)
end

"Repeat unit of an end in forward orientation, as a rotation-invariant key (\"\" if none)."
function unit_key(m::EndModel, forward::Bool)
    m.period == 0 && return ""
    u = m.x[end-m.period+1:end]
    forward || (u = reverse_complement(u))
    return minimum(String(copy(circshift(u, i))) for i in 0:m.period-1)
end

observed(m::EndModel) = m.x[m.k+1:end]                    # read-supported bases beyond the end
"Observed bases between the end and the repeat (all observed bases if no repeat was seen)."
flank(m::EndModel) = m.period > 0 ? m.x[m.k+1:max(m.k, m.run_start-1)] : observed(m)

"L repeat bases continuing the periodic run in phase, starting right after the flank."
function repeat_bases(m::EndModel, L::Int)
    a = max(m.run_start, m.k + 1)
    return UInt8[m.x[m.run_start + mod(a - m.run_start + i - 1, m.period)] for i in 1:L]
end

"Gap of length G = right flank + repeat + left flank (forward strand), or nothing."
function gap_candidate(right::EndModel, left::EndModel, G::Int)
    rf = flank(right); lf = reverse_complement(flank(left))
    L = G - length(rf) - length(lf)
    L >= 0 || return nothing
    R = right.period > 0 ? repeat_bases(right, L) :
        left.period > 0 ? reverse_complement(repeat_bases(left, L)) : nothing
    R === nothing && return nothing
    return vcat(rf, R, lf)
end

"The candidate must reproduce every base the reads showed beyond both ends."
function gap_consistent(g::Vector{UInt8}, right::EndModel, left::EndModel)
    ro = observed(right); m = min(length(ro), length(g))
    g[1:m] == ro[1:m] || return false
    lo = observed(left); rg = reverse_complement(g); m = min(length(lo), length(g))
    return rg[1:m] == lo[1:m]
end

"""
Read-observed bridge: the overhang beyond one end runs into the other end of the path.
Returns (gap bases, reads carrying the whole bridge) or nothing. `target` is the far end
in the overhang's orientation; the bridge is accepted only if every overlapping base agrees.
"""
function read_bridge(obs::Vector{UInt8}, target::Vector{UInt8}, cands::Vector{ExtensionCandidate}; anchor::Int=20)
    (length(obs) >= anchor && length(target) >= anchor) || return nothing
    for j in 1:length(obs)-anchor+1
        @view(obs[j:j+anchor-1]) == @view(target[1:anchor]) || continue
        m = min(length(obs) - j + 1, length(target))
        @view(obs[j:j+m-1]) == @view(target[1:m]) || continue
        full = obs[1:j+anchor-1]
        support = count(x -> length(x.extension) >= length(full) &&
                             @view(x.extension[1:length(full)]) == full, cands)
        support >= 1 || continue
        return obs[1:j-1], support
    end
    return nothing
end

"""
Close a linear genome whose two ends stop at (or in) the same short tandem repeat.
The repeat and the bases next to it come from read overhangs at both ends; the gap
length from read pairs spanning it. The gap sequence must read identically from both
ends, which fixes the repeat phase. Every attempt, closed or not, is reported.
"""
longest_run(a::AbstractVector{UInt8}) = (best = 0; run = 0; prev = 0x00;
    for x in a; run = x == prev ? run + 1 : 1; prev = x; best = max(best, run); end; best)

"Start index of `needle` in `hay` at or after `from` (0 if absent)."
function find_sub(hay::AbstractVector{UInt8}, needle::AbstractVector{UInt8}, from::Int=1)
    m = length(needle)
    @inbounds for i in max(from, 1):length(hay)-m+1
        ok = true
        for j in 1:m
            hay[i+j-1] == needle[j] || (ok = false; break)
        end
        ok && return i
    end
    return 0
end

function count_sub(hay::AbstractVector{UInt8}, needle::AbstractVector{UInt8})
    n = 0; i = find_sub(hay, needle)
    while i != 0
        n += 1; i = find_sub(hay, needle, i + 1)
    end
    return n
end

"""
Bridge the two ends of a path through a low-complexity junction (e.g. adjacent homopolymers)
that stopped extension on both sides, so the ends end in the same tract or overlap inside it.
Anchors are the complex (>= 3 bases, no run > 6) 20-mers nearest to each end that occur once
in the path on both strands; reads containing both anchors in order (either strand) show the
junction directly. Support is counted in distinct read pairs (both mates of one fragment
are one observation). The most frequent junction is used if >= 2 pairs agree and no other
junction is as frequent; the low-complexity end bases outside the anchors are replaced.
Returns (circle, pairs, first kept base, last kept base) or nothing.
"""
function anchored_bridge(seq::Vector{UInt8}, idx::ReadIndex, c::Config; alen::Int=20, reach::Int=120)
    n = length(seq)
    n > 4*(reach + alen) || return nothing
    rc = reverse_complement(seq)
    good(a) = length(unique(a)) >= 3 && longest_run(a) <= 6 && !any(==(UInt8('N')), a) &&
              count_sub(seq, a) + count_sub(rc, a) == 1
    pr = 0
    for L in 0:reach
        good(@view seq[n-L-alen+1:n-L]) && (pr = n - L; break)
    end
    pl = 0
    for L in 0:reach
        good(@view seq[1+L:alen+L]) && (pl = 1 + L; break)
    end
    (pr > 0 && pl > 0) || return nothing
    ra = seq[pr-alen+1:pr]; la = seq[pl:pl+alen-1]
    ids = UInt32[]
    for region in (seq[max(1, pr-99):pr], seq[pl:min(n, pl+99)])
        got, saturated = candidate_read_ids(region, idx, c)
        saturated || append!(ids, got)
    end
    unique!(sort!(ids))
    middles = Dict{Vector{UInt8},Set{UInt32}}()
    read = UInt8[]
    for id in ids
        get_read!(read, idx.store, id)
        for oriented in (read, reverse_complement(read))
            i = find_sub(oriented, ra); i == 0 && continue
            j = find_sub(oriented, la, i + alen); j == 0 && continue
            push!(get!(Set{UInt32}, middles, oriented[i+alen:j-1]), id >> 1)
        end
    end
    isempty(middles) && return nothing
    ranked = sort([(m, length(p)) for (m, p) in middles]; by=x -> -x[2])
    junction, reads = ranked[1]
    (reads >= 2 && (length(ranked) == 1 || ranked[2][2] < reads)) || return nothing
    return vcat(seq[pl:pr], junction), reads, pl, pr
end

const FORK_LOG_LOCK = ReentrantLock()
const FORK_LOG = String[]
const TANDEM_LOG_LOCK = ReentrantLock()
const TANDEM_LOG = String[]

"""
A path end stopped inside a tandem repeat (period 1-6): resolve the whole repeat from reads
that contain the unique entry flank (20-mer before the repeat), then the repeat, then >= 20
exit bases; support = distinct read pairs. Where the path goes next is decided by the exit
flank: the best exit needs >= Min coverage pairs and at least twice the pairs of any other
exit. The repeat length is then the most frequent one among the pairs with that exit (short
reads vary in homopolymer and microsatellite length, so the length is a mode, not a vote
that must be decisive). Returns (bases to trim from the end, bases to append, pairs with the
exit, pairs with the next exit, repeat length) or nothing.
"""
function resolve_tandem_end(context::Vector{UInt8}, idx::ReadIndex, c::Config; alen::Int=20,
                            phase::Int=c.phase_window, genome::Vector{UInt8}=context)
    period, run = periodic_suffix(context)
    period > 0 || return nothing
    L = length(context); before = L-run
    before >= alen+30 || return nothing
    entry = context[before-alen+1:before]
    complex(a) = length(unique(a)) >= 3 && longest_run(a) <= 6 && !(UInt8('N') in a)
    complex(entry) && periodic_suffix(entry)[1] == 0 || return nothing
    unique_anchor(dna_string(genome),entry) || return nothing
    motif = context[before+1:before+period]
    ids,saturated = candidate_read_ids(context[max(1,before-99):before],idx,c)
    saturated && return nothing
    # A pair with incompatible placements/alleles abstains entirely.
    ballots = Dict{UInt32,Set{Tuple{Int,String}}}()
    read = UInt8[]
    for id in ids
        get_read!(read,idx.store,id)
        for oriented in (read,reverse_complement(read))
            count_sub(oriented,entry) == 1 || continue
            i = find_sub(oriented,entry); p = i+alen
            available = min(before,p-1)
            required = max(alen,max(0,phase-run))
            available >= required || continue
            sequences_similar(@view(context[before-available+1:before]),
                              @view(oriented[p-available:p-1]),c.min_identity) || continue
            @view(context[before-required+1:before]) == @view(oriented[p-required:p-1]) || continue
            r = 0
            while p+r <= length(oriented) && oriented[p+r] == motif[mod(r,period)+1]
                r += 1
            end
            r >= period || continue
            p+r+alen-1 <= length(oriented) || continue
            e = oriented[p+r:p+r+alen-1]
            complex(e) || continue
            push!(get!(Set{Tuple{Int,String}},ballots,pair_id(id)),(r,dna_string(e)))
        end
    end
    variants = Dict{Tuple{Int,String},Int}()
    for votes in values(ballots)
        length(votes) == 1 || continue
        v = only(votes); variants[v] = get(variants,v,0)+1
    end
    isempty(variants) && return nothing
    exits = Dict{String,Int}()
    for ((_,e),n) in variants
        exits[e] = get(exits,e,0)+n
    end
    ranked = sort(collect(exits);by=x -> (-last(x),first(x)))
    exit,best = first(ranked); second = length(ranked)>1 ? last(ranked[2]) : 0
    best >= c.min_coverage && best >= 2second || return nothing
    lengths = sort([(r,n) for ((r,e),n) in variants if e==exit];by=x -> (-x[2],x[1]))
    r,support = first(lengths)
    support >= c.min_coverage || return nothing
    length(lengths)>1 && lengths[2][2]==support && return nothing
    return max(0,run-r),vcat(motif_bases(motif,run,r),Vector{UInt8}(codeunits(exit))),best,second,r
end

"Repeat bases from position run+1 to r of a repeat with this motif (empty if r <= run)."
motif_bases(motif::Vector{UInt8}, run::Int, r::Int) = UInt8[motif[k % length(motif) + 1] for k in run:r-1]

"Remove the last n extension bases of one side (never seed bases: caller checks)."
function trim_extension!(b::ContigBuffer, n::Int, direction::Symbol)
    side = direction == :right ? b.right : b.left
    resize!(side, length(side) - n)
    return b
end

"""
Safety net for circles closed without read evidence at the seam (read bridge, repeat-length
estimate): for plastomes (Assembly type chl) the circle must order as LSC-IR-SSC-IR with
IR <= 30% of the circle, SSC >= 1 kb and LSC >= SSC. A path that skipped or duplicated a
region can otherwise be bridged into a circle of plausible length but wrong structure.
"""
function plausible_plastome(circle::Vector{UInt8}, c::Config)
    c.assembly_type == "chl" || return true
    item, _ = PlastomeOrder.order_plastome(circle; min_ir=c.ir_min_length, min_sc=c.ir_min_sc_length,
                                           min_identity=c.ir_min_identity)
    item === nothing && return false
    return item.ir_pair.len <= 0.3*length(circle) && item.ssc_length >= 1000 && item.lsc_length >= item.ssc_length
end

function repeat_gap_closure(res::AssemblyResult, idx::ReadIndex, c::Config)
    ro = ""; lo = ""
    log(status, args...) = (lock(() -> push!(GAP_CLOSURES, GapClosure(status, res.id, args..., ro, lo)), GAP_CLOSURES_LOCK); nothing)
    seq = res.sequence; n = length(seq)
    ctx = max(2*idx.max_read_length, insert_bounds(c)[2] + idx.max_read_length)
    n > 2ctx || return log("contig_too_short", 0, "", 0, NaN, 0, 0, 0, NaN)
    k = 40
    rc_cand = end_candidates(seq[n-ctx+1:n], idx, c)
    lc_cand = end_candidates(reverse_complement(seq[1:ctx]), idx, c)
    right = end_model(seq[n-k+1:n], rc_cand, 2)                 # flanks: >= 2 reads per base
    left = end_model(reverse_complement(seq[1:k]), lc_cand, 2)
    right1 = end_model(seq[n-k+1:n], rc_cand, 1)                # repeat identity: >= 1 read
    left1 = end_model(reverse_complement(seq[1:k]), lc_cand, 1)
    ro = String(copy(first(observed(right1), 60)))
    lo = String(reverse_complement(first(observed(left1), 60)))
    # 1. (16.22) The former bridge anchored on the first bases of the path, which inside a
    #    microsatellite occur many times, and accepted a single read. Only the anchored bridge
    #    below remains: unique complex anchors on both sides and >= 2 distinct read pairs.
    # Ends stopped in the same low-complexity tract (homopolymer, microsatellite): pairs whose
    #    reads contain the nearest unique anchor on both sides show the junction directly.
    ab = anchored_bridge(seq, idx, c)
    if ab !== nothing
        circle, reads, pl, pr = ab
        if c.genome_range[1] <= length(circle) <= c.genome_range[2] && plausible_plastome(circle, c)
            log("closed_read_bridge(pairs=$(reads),anchored,trimmed=$(pl-1)+$(n-pr))", n+1, "",
                length(circle) - (pr - pl + 1), NaN, 0, 0, 0, NaN)
            shifted = [ContigVariant(v.position - (pl-1), v.major, v.minor, v.major_support, v.minor_support,
                                     v.informative, v.switch) for v in res.variants if pl <= v.position <= pr]
            return AssemblyResult(res.id, circle, true, "closed_read_bridge", res.source,
                                  0, reads, 0, shifted)
        end
    end
    rkey = end_repeat_key(right1, true); lkey = end_repeat_key(left1, false)
    motif = isempty(rkey) ? lkey : rkey
    isempty(motif) && return log("no_tandem_repeat_at_ends", 0, "", 0, NaN, 0, 0, 0, NaN)
    # Fallback when the repeat length cannot be estimated: only if both ends ran into the
    # same repeat (the locus is shared), close with observed bases + Unknown gap Ns.
    # Same locus: both ends run into the same repeat unit (one read suffices here).
    same_repeat = !isempty(rkey) && rkey == lkey
    function unknown_gap(reason, pairs, gmed, gmin, gmax, insert)
        (c.unknown_gap_ns > 0 && same_repeat && pairs >= c.min_paired_support) ||
            return log(reason, 0, motif, 0, gmed, gmin, gmax, pairs, insert)
        gapfill = vcat(observed(right), fill!(Vector{UInt8}(undef, c.unknown_gap_ns), UInt8('N')),
                       reverse_complement(observed(left)))
        ncircle = vcat(seq, gapfill)
        c.genome_range[1] <= length(ncircle) <= c.genome_range[2] ||
            return log("circle_length_outside_genome_range", n+1, motif, length(gapfill), gmed, gmin, gmax, pairs, insert)
        # Without a length estimate the circle must at least be a complete plastome:
        # one IR pair with both single-copy regions (nothing else can be missing).
        if c.assembly_type == "chl"
            item, _ = PlastomeOrder.order_plastome(ncircle; min_ir=c.ir_min_length, min_sc=c.ir_min_sc_length,
                                                   min_identity=c.ir_min_identity)
            item === nothing &&
                return log("unknown_gap_circle_not_quadripartite", n+1, motif, length(gapfill), gmed, gmin, gmax, pairs, insert)
            # Ordering alone accepts chimeric paths whose duplicated halves look like a huge IR.
            (item.ir_pair.len <= 0.3*length(ncircle) && item.ssc_length >= 1000 && item.lsc_length >= item.ssc_length) ||
                return log("unknown_gap_structure_implausible(IR=$(item.ir_pair.len),SSC=$(item.ssc_length))",
                           n+1, motif, length(gapfill), gmed, gmin, gmax, pairs, insert)
        end
        log("closed_unknown_gap_length($(reason))", n+1, motif, length(gapfill), gmed, gmin, gmax, pairs, insert)
        return AssemblyResult(res.id, ncircle, true, "closed_gap_unknown_length", res.source,
                              0, 0, pairs, res.variants)
    end
    sums = end_spanning_sums(seq, idx, c)
    length(sums) >= c.min_paired_support ||
        return unknown_gap("spanning_pairs_below_min", length(sums), NaN, 0, 0, NaN)
    insert = contig_insert_median(seq, idx, c)
    isnan(insert) && return unknown_gap("insert_size_unknown", length(sums), NaN, 0, 0, NaN)
    gaps = [insert - x for x in sums]
    gmed = median(gaps); target = round(Int, gmed)
    gmin = round(Int, minimum(gaps)); gmax = round(Int, maximum(gaps))
    if target < 0
        # Negative gap: the ends overlap (both ran through the same repeat). The pairs spanning
        # the ends tie the two loci together; the overlap length is searched near the estimate
        # and must be an exact end match. Its copy number is estimated, not observed.
        win = max(12, 2*max(right.period, left.period))
        O = 0
        for o in sort(max(1, -target-win):(-target+win); by=o -> (abs(o+target), o))
            o <= n ÷ 4 || continue
            @view(seq[n-o+1:n]) == @view(seq[1:o]) || continue
            O = o; break
        end
        O > 0 || return unknown_gap("negative_gap_overlap_inconsistent", length(sums), gmed, gmin, gmax, insert)
        circle = seq[1:n-O]
        c.genome_range[1] <= length(circle) <= c.genome_range[2] ||
            return log("circle_length_outside_genome_range", n+1, motif, -O, gmed, gmin, gmax, length(sums), insert)
        plausible_plastome(circle, c) ||
            return log("closed_circle_structure_implausible", n+1, motif, -O, gmed, gmin, gmax, length(sums), insert)
        log("closed_overlap_collapsed", n+1, motif, -O, gmed, gmin, gmax, length(sums), insert)
        return AssemblyResult(res.id, circle, true, "closed_estimated_repeat", res.source,
                              O, 0, length(sums), [v for v in res.variants if v.position <= n - O])
    end
    # Search lengths near the estimate; the phase of the repeat and both flanks are fixed
    # by requiring the gap to reproduce all bases observed from both ends.
    window = max(12, 2*max(right.period, left.period))
    fill = nothing
    for G in sort(max(0, target-window):target+window; by=g -> (abs(g-target), g))
        g = gap_candidate(right, left, G)
        g === nothing && continue
        gap_consistent(g, right, left) || continue
        fill = g; break
    end
    fill === nothing && return unknown_gap("repeat_phase_or_flanks_inconsistent", length(sums), gmed, gmin, gmax, insert)
    circle = vcat(seq, fill)
    c.genome_range[1] <= length(circle) <= c.genome_range[2] ||
        return log("circle_length_outside_genome_range", n+1, motif, length(fill), gmed, gmin, gmax, length(sums), insert)
    # A side from which no read went beyond the end contributes no flank: its few
    # boundary bases (if any) are inferred from the repeat alone - flagged separately.
    status = isempty(observed(right)) || isempty(observed(left)) ? "closed_one_flank_unobserved" : "closed"
    log(status, n+1, motif, length(fill), gmed, gmin, gmax, length(sums), insert)
    plausible_plastome(circle, c) ||
        return log("closed_circle_structure_implausible", n+1, motif, length(fill), gmed, gmin, gmax, length(sums), insert)
    return AssemblyResult(res.id, circle, true, "closed_estimated_repeat", res.source,
                          0, 0, length(sums), res.variants)
end

"""
Try repeat_gap_closure on the linear paths when none is circular: shortest plausible genome
first (longer paths sharing its ends usually repeat the IR); at most one is closed.
"""
function rescue_closure!(results::Vector{AssemblyResult}, idx::ReadIndex, c::Config)
    (c.repeat_gap_closure && !any(x -> x.circular, results)) || return results
    for i in sortperm(results; by=x -> (length(x.sequence), x.id))
        x = results[i]
        length(x.sequence) >= c.genome_range[1] - 5000 || continue
        closed = repeat_gap_closure(x, idx, c)
        closed === nothing && continue              # reason is in _estimated_gaps.tsv
        results[i] = closed
        break
    end
    return results
end

function assemble(c::Config, idx::ReadIndex; quiet::Bool=false, monitor=nothing)
    lock(CLOSURE_LOG_LOCK) do
        empty!(CLOSURE_LOG)
    end
    lock(() -> empty!(GAP_CLOSURES), GAP_CLOSURES_LOCK)
    lock(() -> empty!(COVERAGE_TRACE), COVERAGE_TRACE_LOCK)
    lock(() -> empty!(TANDEM_LOG), TANDEM_LOG_LOCK)
    lock(() -> empty!(FORK_LOG), FORK_LOG_LOCK)
    lock(() -> empty!(MATE_JUMP_LOG), MATE_JUMP_LOG_LOCK)
    progress_phase!(monitor,"Seeding")
    seeds = find_initial_reads(c, idx)
    progress_phase!(monitor,"Assembly";total=length(seeds),unit="seeds")
    all_results = AssemblyResult[]
    # A bounded worker pool consumes the next seed as soon as a worker finishes.
    # Result slots preserve input order, independent of completion/thread order.
    parallel = min(c.max_concurrent_assemblies,Threads.nthreads(),length(seeds))
    jobs = Channel{Int}(length(seeds))
    for i in eachindex(seeds)
        put!(jobs,i)
    end
    close(jobs)
    per_seed = Vector{Vector{AssemblyResult}}(undef,length(seeds))
    @sync for _ in 1:parallel
        Threads.@spawn begin
            for i in jobs
                per_seed[i] = assemble_seed(seeds[i],idx,c;monitor=monitor)
            end
        end
    end
    for (i,results) in enumerate(per_seed)
        append!(all_results,results)
        !quiet && @printf("Completed %s: %d paths, longest %d bp\n",seeds[i].id,
                          length(results),maximum(length(x.sequence) for x in results))
    end
    progress_update!(monitor;done=length(seeds),force=true)
    # Closure rescue before deduplication: deduplication drops a shorter linear path contained
    # in a longer one (e.g. one that repeats the IR), and the shorter one may be the closable one.
    rescued_results = copy(all_results)            # raw results keep their own stop reasons
    rescue_closure!(rescued_results, idx, c)
    progress_phase!(monitor,"Deduplication";detail="$(length(all_results)) paths")
    unique_results = deduplicate(rescued_results)
    progress_phase!(monitor,"Scaffolding";detail=c.scaffold ? "checking overlaps" : "disabled")
    final = deduplicate(scaffold_contigs(unique_results, idx, c;quiet=quiet,monitor=monitor))
    rescue_closure!(final, idx, c)                  # scaffolds may have become closable
    sort!(final; by=x -> (-length(x.sequence), x.id))
    return final, all_results
end

"""
How a circle was closed. Only `validated` has read and pair evidence across the junction
(Circular min reads, Min paired support); the others are closed by inference and must not
be reported as validated circles.
"""
function closure_status(x::AssemblyResult)
    x.circular || return "linear"
    x.reason == "closed_ir_collapse" && return "ir_collapse"            # terminal IR copy collapsed
    x.reason == "closed_topology_inferred" && return "topology_inferred"
    x.reason == "closed_graph_path" && return "graph_path"              # circular walk of the path graph
    x.reason == "closed_inferred_join" && return "inferred_join"        # joined by motif, not by data
    x.reason == "closed_ref_scaffold" && return "ref_scaffold"          # paths in reference order, Ns between
    x.reason == "closed_read_bridge" && return "read_bridge"            # gap observed in reads from >= 2 distinct pairs
    x.reason == "closed_estimated_repeat" && return "estimated_repeat"  # repeat length from inserts or few spanning reads
    x.reason == "closed_gap_unknown_length" && return "gap_unknown"     # N gap of unknown length
    if x.reason == "validated_closure"                                   # read and pair evidence at the junction
        # A validated junction does not validate an internal N gap (e.g. from a mate jump).
        return isempty(n_gaps(x.sequence)) ? "validated" : "gap_unknown"
    end
    warn_unknown_closure(x.reason)
    return "unresolved"                                                  # unknown reason: never validated
end

const UNKNOWN_CLOSURES = Set{String}()
const UNKNOWN_CLOSURES_LOCK = ReentrantLock()
"Warns once per unknown closure reason (a circle with it gets status unresolved)."
function warn_unknown_closure(reason::AbstractString)
    first_time = lock(UNKNOWN_CLOSURES_LOCK) do
        reason in UNKNOWN_CLOSURES ? false : (push!(UNKNOWN_CLOSURES, String(reason)); true)
    end
    first_time && @warn "Unknown closure reason $(repr(reason)): circle reported as unresolved"
    return nothing
end

"""
Output order (rank 1 = primary result): validated (gap-free) circles in Genome Range, then circles in
range closed by IR collapse, read bridge, repeat-length estimate or N gap, then circles
outside the range, then linear contigs; longer first within a tier. Length alone never decides.
"""
function result_tier(x::AssemblyResult, c::Config)
    in_range = c.genome_range[1] <= length(x.sequence) <= c.genome_range[2]
    st = closure_status(x)
    x.circular || return 9
    in_range || return 8
    st == "topology_inferred" && return 6
    st == "inferred_join" && return 6                # adjacency inferred, not observed
    st == "unresolved" && return 6                   # unknown closure reason
    st == "ref_scaffold" && return 7                 # adjacency from the reference only (rescue)
    isempty(n_gaps(x.sequence)) || return 5         # internal gap of unknown length
    return st == "validated" ? 1 : st == "ir_collapse" ? 2 : st in ("read_bridge", "graph_path") ? 3 :
           st == "estimated_repeat" ? 4 : 5
end

const READ_STATS = Ref{Any}(nothing)

"""
Approximate plastid share of the library: reads sharing at least one indexed minimizer with
the primary result (both strands, across the circle origin), their share of all reads, and
mean coverage = their bases / genome length. Reads of plastid-like nuclear or mitochondrial
copies (NUPT/NUMT) are included; hot minimizers removed by the frequency cutoff are not.
"""
function plastid_read_stats(idx::ReadIndex, x::AssemblyResult, c::Config)
    seq = x.circular ? vcat(x.sequence, x.sequence[1:min(c.kmer-1, length(x.sequence))]) : x.sequence
    store = idx.store
    marked = falses(2*num_pairs(store) + 2)
    for key in all_kmer_keys(seq, c.kmer)
        low_complexity_kmer(key, c.kmer) && continue
        bucket = get(idx.postings, key, nothing)
        bucket === nothing && continue
        for id in bucket
            marked[id] = true
        end
    end
    n = 0; bases = 0
    for id in findall(marked)
        pair = Int(id) >> 1
        (pair < 1 || pair > num_pairs(store)) && continue
        n += 1
        bases += isodd(id) ? store.r1.lengths[pair] : store.r2.lengths[pair]
    end
    total = num_pairs(store) * (store.paired ? 2 : 1)
    return (id=x.id, reads=n, total=total, percent=100*n/max(total, 1), coverage=bases/length(x.sequence))
end

function write_fasta(path::String, results::Vector{AssemblyResult}; roles::Bool=false, extra::AbstractString="")
    open(path, "w") do io
        for (rank, item) in enumerate(results)
            topology = item.circular ? "circular" : "linear"
            closure = item.circular ? " closure=$(closure_status(item))" : ""
            role = roles ? (rank == 1 ? " role=primary" : " role=alternative") : ""
            gaps = length(n_gaps(item.sequence))
            gaptag = gaps > 0 ? " internal_gaps=$gaps" : ""
            println(io, ">$(item.id) length=$(length(item.sequence)) topology=$topology$closure$gaptag rank=$rank$role$extra")
            for i in 1:80:length(item.sequence)
                println(io, dna_string(@view item.sequence[i:min(i+79, end)]))
            end
        end
    end
end

clean_tsv(s::AbstractString) = replace(s, '\t' => ' ', '\n' => ' ', '\r' => ' ')

"Forward and reverse-complement exact k-mer matches of `a` against `b` (both read 5'->3')."
function strand_votes(a::AbstractVector{UInt8}, b::AbstractVector{UInt8}; k::Int=21)
    function kmers(s)
        out = UInt64[]; mask = (UInt64(1) << (2k)) - UInt64(1); f = UInt64(0); run = 0
        for x in s
            code = PlastomeOrder.BASE_CODE[Int(x)+1]
            if code > 3
                run = 0; f = UInt64(0); continue
            end
            f = ((f << 2) | UInt64(code)) & mask; run += 1
            run >= k && push!(out, f)
        end
        out
    end
    ref = Set(kmers(b))
    fwd = count(in(ref), kmers(a)); rev = count(in(ref), kmers(PlastomeOrder.revcomp(a)))
    return fwd, rev
end
decisive(f, r) = max(f, r) >= 100 && max(f, r) >= 3 * min(f, r)

"""
Orientation from the nearest reference (both put in LSC-IR-SSC-IR): the molecule strand whose
LSC runs the same way as the reference LSC (so the FASTA starts at the same LSC end); for an
SSC without ndhF, the SSC orientation of the reference. Needs >= 100 matching 21-mers and a
3x majority; otherwise the assembly orientation is kept. Returns (item, label): label assembly,
reference, or reference_revcomp (the whole molecule was reverse-complemented first).
"""
function orient_by_reference(seq::Vector{UInt8}, item, ref, c::Config)
    ref === nothing && return item, "assembly"
    rl = ref.lsc_length
    f, r = strand_votes(@view(item.sequence[1:item.lsc_length]), @view(ref.sequence[1:rl]))
    decisive(f, r) || return item, "assembly"
    label = "reference"
    if r > f
        flipped, _ = PlastomeOrder.order_plastome(PlastomeOrder.revcomp(seq); min_ir=c.ir_min_length,
                         min_sc=c.ir_min_sc_length, min_identity=c.ir_min_identity)
        flipped === nothing && return item, "assembly"
        item = flipped; label = "reference_revcomp"    # coordinates of item now refer to revcomp(seq)
    end
    if item.ndhf == "not_found"
        s0 = item.lsc_length + item.ir_pair.len; r0 = ref.lsc_length + ref.ir_pair.len
        f2, r2 = strand_votes(@view(item.sequence[s0+1:s0+item.ssc_length]), @view(ref.sequence[r0+1:r0+ref.ssc_length]))
        if decisive(f2, r2) && r2 > f2
            q = copy(item.sequence); q[s0+1:s0+item.ssc_length] = PlastomeOrder.revcomp(q[s0+1:s0+item.ssc_length])
            item = PlastomeOrder.OrderedPlastome(q, item.origin0, item.ir_pair, item.lsc_length, item.ssc_length,
                                                 item.method, !item.ssc_inverted, "not_found_SSC_by_reference")
        end
    end
    return item, label
end

# ------------------------------------------------------------------ molecule orientation by rpoB
# Convention: rpoB (and the rpoC genes next to it) on the + strand of the ordered FASTA. The
# rpoB proteins come from the CDS features of the plastid panel (GenBank files spanning land
# plants); rpoB is found in our LSC by protein 5-mers in the six-frame translation, as ndhF is.
const RPOB_KMERS = Ref{Any}(nothing)            # (panel, Set{String}) of the last panel read

"Protein sequences of the CDS features named `gene` in a GenBank file (complement handled)."
function genbank_cds_proteins(path::AbstractString, gene::AbstractString)
    seq = UInt8[]; feats = Tuple{String,String}[]          # (location, gene)
    state = :header; key = ""; loc = IOBuffer(); g = ""; inloc = false
    flush!() = (key == "CDS" && push!(feats, (String(take!(loc)), g)); take!(loc); key = ""; g = "")
    for line in eachline(path)
        if state == :header
            startswith(line, "FEATURES") && (state = :features)
        elseif state == :features
            if startswith(line, "ORIGIN")
                flush!(); state = :origin; continue
            end
            length(line) >= 22 || continue
            fkey = strip(line[6:min(20, end)]); body = strip(line[22:end])
            if !isempty(fkey)
                flush!(); key = fkey; print(loc, body); inloc = true
            elseif startswith(body, "/")
                inloc = false
                m = match(r"^/gene=\"?([^\"]*)\"?", body)
                m === nothing || isempty(g) && (g = String(m.captures[1]))
            elseif inloc
                print(loc, body)
            end
        else
            startswith(line, "//") && break
            for b in codeunits(line)
                UInt8('a') <= b <= UInt8('z') && (b -= 0x20)
                UInt8('A') <= b <= UInt8('Z') && push!(seq, b)
            end
        end
    end
    out = String[]
    for (l, name) in feats
        lowercase(name) == lowercase(gene) || continue
        occursin(':', l) && continue                            # remote location
        ranges = sort([(parse(Int, m.captures[1]), parse(Int, m.captures[2]))
                       for m in eachmatch(r"<?(\d+)\.\.>?(\d+)", l)])
        (isempty(ranges) || !all(r -> 1 <= r[1] <= r[2] <= length(seq), ranges)) && continue
        nt = reduce(vcat, [seq[a:b] for (a, b) in ranges])
        occursin("complement", l) && (nt = PlastomeOrder.revcomp(nt))
        push!(out, PlastomeOrder.translate_frame(nt, 0))
    end
    return out
end

"rpoB protein 5-mers from the GenBank files of the plastid panel (read once per panel)."
function rpob_kmers(panel::AbstractString)
    (RPOB_KMERS[] !== nothing && RPOB_KMERS[][1] == panel) && return RPOB_KMERS[][2]
    set = Set{String}(); k = PlastomeOrder.NDHF_K
    files = isdir(panel) ? [joinpath(panel, f) for f in readdir(panel) if occursin(r"\.(gb|gbk|gbff)$"i, f)] :
            (isfile(panel) && occursin(r"\.(gb|gbk|gbff)$"i, panel)) ? [panel] : String[]
    for f in files, prot in genbank_cds_proteins(f, "rpoB"), i in 1:length(prot)-k+1
        w = prot[i:i+k-1]
        occursin('*', w) || occursin('X', w) || push!(set, w)
    end
    RPOB_KMERS[] = (panel, set)
    return set
end

"""
rpoB strand in the LSC of an ordered plastome: :plus, :minus or :absent. Protein 5-mer hits are
counted per strand in the densest `window`-bp stretch (whole-sequence counts are dominated by
chance hits in translated non-coding sequence); the strand with >= 100 hits there and a 3x
majority over the other strand's densest stretch wins.
"""
function rpob_strand(seq::AbstractVector{UInt8}, lsc::Int, kmers::Set{String}; window::Int=4000)
    isempty(kmers) && return :absent
    k = PlastomeOrder.NDHF_K
    s = Vector{UInt8}(seq[1:lsc]); n = length(s)
    function positions(t, forward)
        out = Int[]
        for f in 0:2
            aa = PlastomeOrder.translate_frame(t, f)
            for i in 1:length(aa)-k+1
                aa[i:i+k-1] in kmers || continue
                nt = f + 3(i-1)
                push!(out, forward ? nt : n - 1 - nt)
            end
        end
        return sort!(out)
    end
    function densest(h)
        best = 0; j = 1
        for i in eachindex(h)
            while h[i] - h[j] > window
                j += 1
            end
            best = max(best, i - j + 1)
        end
        return best
    end
    fw = densest(positions(s, true)); rv = densest(positions(PlastomeOrder.revcomp(s), false))
    decisive(fw, rv) || return :absent
    return fw > rv ? :plus : :minus
end

"SSC without ndhF: the SSC orientation of the reference (as in orient_by_reference)."
function ssc_by_reference(item, ref)
    (ref === nothing || item.ndhf != "not_found") && return item
    s0 = item.lsc_length + item.ir_pair.len; r0 = ref.lsc_length + ref.ir_pair.len
    f2, r2 = strand_votes(@view(item.sequence[s0+1:s0+item.ssc_length]), @view(ref.sequence[r0+1:r0+ref.ssc_length]))
    if decisive(f2, r2) && r2 > f2
        q = copy(item.sequence); q[s0+1:s0+item.ssc_length] = PlastomeOrder.revcomp(q[s0+1:s0+item.ssc_length])
        item = PlastomeOrder.OrderedPlastome(q, item.origin0, item.ir_pair, item.lsc_length, item.ssc_length,
                                             item.method, !item.ssc_inverted, "not_found_SSC_by_reference")
    end
    return item
end

"""
Molecule orientation: rpoB on the + strand (label rpoB, or rpoB_revcomp when the whole molecule
was reverse-complemented and ordered again); without a decisive rpoB call, the nearest dataset
reference (orient_by_reference). The SSC keeps ndhF at its start, next to the first IR, in every case
(order_plastome).
"""
function orient_molecule(seq::Vector{UInt8}, item, ref, c::Config, kmers::Set{String})
    st = rpob_strand(item.sequence, item.lsc_length, kmers)
    st == :absent && return orient_by_reference(seq, item, ref, c)
    label = "rpoB"
    if st == :minus
        flipped, _ = PlastomeOrder.order_plastome(PlastomeOrder.revcomp(seq); min_ir=c.ir_min_length,
                         min_sc=c.ir_min_sc_length, min_identity=c.ir_min_identity)
        flipped === nothing && return orient_by_reference(seq, item, ref, c)
        item = flipped; label = "rpoB_revcomp"         # coordinates of item now refer to revcomp(seq)
    end
    return ssc_by_reference(item, ref), label
end

"Write a separate rotated FASTA; closure evidence keeps its original coordinates."
function write_plastome_outputs(c::Config, accepted, prefix::String; monitor=nothing, reference=nothing)
    # Reference orientation from the nearest dataset record, never from a local reference (our own
    # earlier circles keep the orientation of their assembly, and two libraries that are each
    # other's nearest local reference would orient each other). No containment threshold: the
    # strand vote itself needs >= 100 shared 21-mers and a 3x majority (orient_by_reference),
    # so a distant reference orients only when the vote is decisive.
    ref_item = nothing; ref_tag = ""
    if c.reference_orientation && reference !== nothing && !isempty(reference.orient_sequence)
        ref_item, _ = PlastomeOrder.order_plastome(reference.orient_sequence; min_ir=c.ir_min_length,
                          min_sc=c.ir_min_sc_length, min_identity=c.ir_min_identity)
        ref_tag = ref_item === nothing ? "" : reference.orient_accession
    end
    written = Dict{Vector{UInt8},String}()     # ordered sequence -> first record written with it
    rpk = isempty(c.plastid_panel) ? Set{String}() : rpob_kmers(c.plastid_panel)
    progress_phase!(monitor,"LSC/IR/SSC";total=length(accepted),unit="contigs")
    output = prefix*"_plastome_LSC_IR_SSC_IR.fasta"
    ordered_count = 0
    # Always truncate these outputs, including when disabled, to avoid stale genomes.
    open(output,"w") do fasta
        open(output*".regions.tsv","w") do regions
            open(prefix*"_plastome_ordering.tsv","w") do status_io
                println(regions,PlastomeOrder.REGIONS_HEADER)
                println(status_io,"id\tstatus\toriginal_start\tLSC_length\tIR_length\tSSC_length\tir_identity\tclosure_status\tndhF\torientation\treference")
                for (number,x) in enumerate(accepted)
                    progress_update!(monitor;detail=x.id)
                    status = !c.normalize_plastome ? "disabled" : c.assembly_type != "chl" ?
                             "not_chloroplast" : !x.circular ? "not_circular" : "pending"
                    item = nothing
                    if status == "pending"
                        item,status = PlastomeOrder.order_plastome(x.sequence;
                            min_ir=c.ir_min_length,min_sc=c.ir_min_sc_length,
                            min_identity=c.ir_min_identity)
                    end
                    if item === nothing
                        println(status_io,join((x.id,status,".",".",".",".",".",closure_status(x),".",".","."),'\t'))
                        if c.normalize_plastome && c.assembly_type == "chl" && x.circular
                            progress_breakline!(monitor)
                            @warn "Plastome ordering unresolved; inspect ordering TSV or provide verified IR coordinates with OrderPlastome.jl" contig=x.id status=status
                        end
                        progress_update!(monitor;done=number)
                        continue
                    end
                    item, orientation = orient_molecule(x.sequence, item, ref_item, c, rpk)
                    rc = orientation in ("reference_revcomp", "rpoB_revcomp")
                    # the same molecule ordered twice (e.g. two SSC isomers of one assembly): written once
                    if haskey(written, item.sequence)
                        println(status_io,join((x.id,"duplicate_of:"*written[item.sequence],".",".",".",".",".",
                                closure_status(x),".",".","."),'\t'))
                        progress_update!(monitor;done=number)
                        continue
                    end
                    written[copy(item.sequence)] = x.id
                    extra = "orientation=$orientation" * (startswith(orientation, "reference") ? " ref=$ref_tag" : "")
                    PlastomeOrder.write_ordered_record(fasta,x.id,item;closure=closure_status(x),extra=extra,revcomp=rc)
                    PlastomeOrder.write_region_rows(regions,x.id,item;revcomp=rc)
                    println(status_io,join((x.id,status,rc ? length(item.sequence)-item.origin0 : item.origin0+1,item.lsc_length,
                            item.ir_pair.len,item.ssc_length,item.ir_pair.identity,closure_status(x),item.ndhf,
                            orientation,startswith(orientation, "reference") ? ref_tag : "."),'\t'))
                    ordered_count += 1
                    progress_update!(monitor;done=number)
                end
            end
        end
    end
    progress_update!(monitor;done=length(accepted),force=true)
    ordered_count > 0 && monitor===nothing && println("Ordered plastome FASTA: ",output)
    return ordered_count
end

# ------------------------------------------------------------------ plastome annotation
const ANNOTATOR_MODULE = Ref{Any}(nothing)
const ANNOTATOR_CACHE = Ref{Any}(nothing)      # ((dir, dataset), Annotator)

"Loads PlastomeAnnotator (PlastomeAnnotator.jl in `dir`) once per process."
function annotator_module(dir::AbstractString)
    ANNOTATOR_MODULE[] === nothing || return ANNOTATOR_MODULE[]
    path = joinpath(dir, "PlastomeAnnotator.jl")
    isfile(path) || error("Annotation dir: PlastomeAnnotator.jl not found in $dir")
    Base.include(@__MODULE__, path)
    # the module was defined in a newer world: look it up there
    ANNOTATOR_MODULE[] = Base.invokelatest(getfield, @__MODULE__, :PlastomeAnnotator)
    return ANNOTATOR_MODULE[]
end

"""
Gene annotation of the ordered plastomes (prefix_plastome_LSC_IR_SSC_IR.fasta) with
PlastomeAnnotator (Annotation dir: directory of PlastomeAnnotator.jl; Annotation dataset: the
directory with refindex_all and models_final). Writes prefix_plastome_LSC_IR_SSC_IR.gff3 in
the coordinates of the ordered FASTA and prefix_plastome_genes.tsv (column editing: 1 when the
annotator applied RNA-editing tolerance to the record). Off when Annotation dir
is empty; a failure is reported and does not stop the run. Returns (records, summary).
"""
function annotate_plastomes(c::Config, prefix::AbstractString; monitor=nothing)
    isempty(c.annotation_dir) && return 0, ""
    fasta = prefix*"_plastome_LSC_IR_SSC_IR.fasta"
    gff = prefix*"_plastome_LSC_IR_SSC_IR.gff3"; tsv = prefix*"_plastome_genes.tsv"
    rm(gff; force=true); rm(tsv; force=true)               # no stale annotation
    (isfile(fasta) && filesize(fasta) > 0) || return 0, "no ordered plastome"
    try
        A = annotator_module(c.annotation_dir)
        key = (c.annotation_dir, c.annotation_dataset)
        if ANNOTATOR_CACHE[] === nothing || ANNOTATOR_CACHE[][1] != key
            ANNOTATOR_CACHE[] = (key, Base.invokelatest(() -> A.load_annotator(c.annotation_dataset)))
        end
        ann = ANNOTATOR_CACHE[][2]
        recs = PlastomeOrder.fasta_records(fasta)
        progress_phase!(monitor, "Annotation"; total=length(recs), unit="plastomes")
        lines = String[]
        open(gff, "w") do g
            open(tsv, "w") do t
                println(g, "##gff-version 3")
                println(t, "id\tgene\tkind\tstrand\tsegments\tparts\tp\tpseudo\tflags\tediting")
                for (i, rec) in enumerate(recs)
                    id = String(first(split(first(rec))))
                    seq = String(copy(last(rec)))
                    progress_update!(monitor; detail=id)
                    # functions of the annotator are called in the latest world (defined by include)
                    ms, ir, editing, nearest = Base.invokelatest(() -> A.annotate_record(ann, id, seq))
                    Base.invokelatest(() -> A.write_gff(g, id, length(seq), ms; header=false))
                    foreach(r -> println(t, r, '\t', editing ? 1 : 0), Base.invokelatest(() -> A.models_rows(id, ms)))
                    near = nearest === nothing ? "-" :
                           @sprintf("%s (%s) Mash %.4f", nearest[1].accession, nearest[1].genus, nearest[2])
                    push!(lines, "$id: $(length(ms)) gene models ($(count(m -> m.pseudo, ms)) pseudo), " *
                                 "editing tolerance $(editing ? "on" : "off"), nearest annotation reference $near")
                    progress_update!(monitor; done=i)
                end
            end
        end
        return length(recs), join(lines, "; ")
    catch err
        progress_breakline!(monitor)
        @warn "Plastome annotation failed; assembly outputs are unaffected" exception=(err, catch_backtrace())
        rm(gff; force=true); rm(tsv; force=true)
        return 0, "failed: " * sprint(showerror, err)
    end
end

"Maximal tandem runs (period 1-6, >= max(12, 3 periods), no N) as (start, stop, period); smallest period first."
function tandem_runs(seq::AbstractVector{UInt8})
    n = length(seq); taken = falses(n); out = Tuple{Int,Int,Int}[]
    for p in 1:6
        i = 1
        while i + p <= n
            if seq[i] == seq[i+p] && seq[i] != UInt8('N')
                j = i
                while j + p < n && seq[j+1] == seq[j+1+p] && seq[j+1] != UInt8('N')
                    j += 1
                end
                start, stop = i, j + p; len = stop - start + 1
                if len >= max(12, 3p) && count(@view taken[start:stop]) <= len ÷ 2
                    push!(out, (start, stop, p)); taken[start:stop] .= true
                end
                i = j + 1
            else
                i += 1
            end
        end
    end
    return sort!(out)
end

"Bases between the first exact `L` and the next exact `R` in `s` (at most `maxlen`), or -1."
function flanked_length(s::AbstractVector{UInt8}, L::AbstractVector{UInt8}, R::AbstractVector{UInt8}, maxlen::Int)
    nl = length(L); nr = length(R); n = length(s)
    for i in 1:n-nl-nr+1
        @view(s[i:i+nl-1]) == L || continue
        a = i + nl
        for j in a:min(n-nr+1, a+maxlen)
            @view(s[j:j+nr-1]) == R && return j - a
        end
        return -1
    end
    return -1
end

"""
Repeat-length variants of tandem runs shorter than the reads. Reads (either strand) holding
both 12-bp flanks exactly are tallied by run length. Reported, with the assembled sequence
unchanged:
  repeat_length (>= 20 such reads): another length with >= 2 reads and >= 20% of the most frequent length
    (>= 35% when it differs by at most one repeat unit: stutter of homopolymers and
    dinucleotides reaches 20%);
  repeat_length_unsupported (>= 10 such reads): no read has the assembled length while
    another has >= 3 reads.
In an IR both copies get the same counts: reads of the two copies cannot be told apart.
"""
function repeat_length_variants(x::AssemblyResult, idx::ReadIndex, c::Config; anchor::Int=12, window::Int=60,
                                min_reads::Int=10, min_reads_variant::Int=20)
    out = Tuple{ContigVariant,String}[]
    seq = x.sequence; n = length(seq)
    n >= 4window || return out
    maxrun = idx.max_read_length - 2anchor - 4
    for (start, stop, p) in tandem_runs(seq)
        len = stop - start + 1
        len <= maxrun || continue
        (start - window >= 1 && stop + window <= n) || continue
        L = seq[start-anchor:start-1]; R = seq[stop+1:stop+anchor]
        (UInt8('N') in L || UInt8('N') in R) && continue
        ids = Set{UInt32}()
        for w in (@view(seq[start-window:start-1]), @view(seq[stop+1:stop+window]))
            for key in all_kmer_keys(w, c.kmer)
                bucket = get(idx.postings, key, nothing)
                bucket === nothing && continue
                for id in bucket
                    excluded_read(UInt32(id), c) || push!(ids, UInt32(id))
                end
            end
        end
        tally = Dict{Int,Int}()
        for id in ids
            r = get_read(idx.store, id)
            for s in (r, reverse_complement(r))
                m = flanked_length(s, L, R, maxrun + 20)
                m < 0 && continue
                tally[m] = get(tally, m, 0) + 1
                break
            end
        end
        total = sum(values(tally); init=0)
        total >= min_reads || continue
        ranked = sort(collect(tally); by=t -> (-last(t), first(t)))
        top = last(first(ranked))
        kind = "repeat_length"
        if get(tally, len, 0) == 0 && top >= 3
            alt = 1; kind = "repeat_length_unsupported"
        else
            total >= min_reads_variant || continue
            alt = findfirst(t -> first(t) != len &&
                                 last(t) >= max(2, (abs(first(t) - len) <= p ? 0.35 : 0.2) * top), ranked)
            alt === nothing && continue
        end
        al, an = ranked[alt]
        motif = seq[start:start+p-1]
        minor = String(UInt8[motif[k % p + 1] for k in 0:al-1])
        push!(out, (ContigVariant(start, String(seq[start:stop]), minor, get(tally, len, 0), an, total), kind))
    end
    return out
end

"The library's own plastome is never its reference (local:<Project name>:...)."
not_self(c::Config) = name -> startswith(name, "local:" * c.project_name * ":")

"""
Nearest reference of the result: the primary circle, else all accepted contigs >= 5 kb. Writes
_nearest_reference.tsv (top 5), _nearest_reference.fasta and _nearest_reference_features.tsv;
returns (summary for _run.txt, "" when no Reference dataset is set; accession; sequence; containment;
and the orientation reference: the nearest record of the dataset itself, local references excluded).
"""
function write_nearest_reference(c::Config, accepted, prefix::String)
    none(msg) = (summary=msg, accession="", sequence=UInt8[], containment=0.0,
                 orient_accession="", orient_sequence=UInt8[], orient_containment=0.0)
    isempty(c.reference_dataset) && return none("")
    sk = joinpath(c.reference_dataset, RefSketch.SKETCH_FILE)
    isfile(sk) || return none("not available (missing $(sk))")
    query = !isempty(accepted) && accepted[1].circular ? accepted[1:1] : [x for x in accepted if length(x.sequence) >= 5000]
    isempty(query) && return none("not available (no contig >= 5 kb)")
    db = RefSketch.load(sk); meta = RefSketch.metadata(c.reference_dataset)
    hits, nq = RefSketch.nearest(db, [x.sequence for x in query], [x.circular for x in query]; exclude=not_self(c))
    isempty(hits) && return none("not available (no k-mers)")
    open(prefix*"_nearest_reference.tsv", "w") do io
        println(io, "rank\taccession\torganism\tgroup\tshared_hashes\tquery_hashes\tcontainment\tani_estimate\treference_covered\treference_length")
        for (i, h) in enumerate(hits)
            org, grp = get(meta, h.accession, ("", ""))
            @printf(io, "%d\t%s\t%s\t%s\t%d\t%d\t%.4f\t%.4f\t%.4f\t%d\n", i, h.accession, org, grp, h.shared, nq,
                    h.containment, h.ani, h.reference_covered, h.reference_length)
        end
    end
    best = hits[1]; org, grp = get(meta, best.accession, ("", ""))
    seq = RefSketch.reference_sequence(c.reference_dataset, best.accession)
    open(prefix*"_nearest_reference.fasta", "w") do io
        println(io, ">", best.accession, " ", org)
        for i in 1:80:length(seq)
            println(io, String(seq[i:min(i+79, end)]))
        end
    end
    open(prefix*"_nearest_reference_features.tsv", "w") do io
        foreach(l -> println(io, l), RefSketch.reference_features(c.reference_dataset, best.accession))
    end
    # orientation reference: nearest dataset record, local references excluded
    ohits, _ = RefSketch.nearest(db, [x.sequence for x in query], [x.circular for x in query]; top=1,
                                 exclude=name -> startswith(name, "local:"))
    ob = isempty(ohits) ? nothing : ohits[1]
    oseq = ob === nothing ? UInt8[] : ob.accession == best.accession ? seq :
           RefSketch.reference_sequence(c.reference_dataset, ob.accession)
    summary = @sprintf("%s %s (%s): containment %.3f, ANI ~%.3f, reference covered %.3f (query: %s)", best.accession, org, grp,
                       best.containment, best.ani, best.reference_covered,
                       length(query) == 1 && query[1].circular ? "primary circle" : "$(length(query)) contigs >= 5 kb")
    ob === nothing || (summary *= @sprintf("; orientation reference %s (containment %.3f)", ob.accession, ob.containment))
    return (summary=summary, accession=best.accession, sequence=seq, containment=best.containment,
            orient_accession=ob === nothing ? "" : ob.accession, orient_sequence=oseq,
            orient_containment=ob === nothing ? 0.0 : ob.containment)
end

# ---------------------------------------------------------------------------
# Reference scaffold (last resort: no circle at all)
# ---------------------------------------------------------------------------
const REF_SCAFFOLD_REPORT = String[]
const REF_SCAFFOLD_REJECT = Ref("")       # why the last reference_scaffold call returned nothing
const REF_SCAFFOLD_HEADER = "junction\tfrom\tto\treference_from\treference_to\tpredicted_bp\tjoin\tbases"

"(1-based start, forward 2-bit k-mer) of every valid window of `s`."
function forward_kmers(s::AbstractVector{UInt8}, k::Int)
    out = Tuple{Int,UInt64}[]; mask = (UInt64(1) << (2k)) - UInt64(1); f = UInt64(0); run = 0
    @inbounds for i in eachindex(s)
        v = ENCODE[Int(s[i])+1]
        if v == 255
            run = 0; f = UInt64(0); continue
        end
        f = ((f << 2) | UInt64(v)) & mask; run += 1
        run >= k && push!(out, (i - k + 1, f))
    end
    return out
end

"""
Placements of a sequence on a circular reference: clusters of exact k-mer anchors on one
diagonal (reference position - query position, modulo n, in `bin`-bp bins; indels move the
diagonal little). One placement per cluster with >= `min_anchors` anchors: the anchored span of
the query and its reference start. A path in the IR gets one placement in each IR copy.
"""
function ref_placements(q::Vector{UInt8}, refidx::Dict{UInt64,Vector{Int32}}, n::Int, k::Int;
                        bin::Int=200, min_anchors::Int=100)
    out = NamedTuple{(:rstart,:seq,:anchors),Tuple{Int,Vector{UInt8},Int}}[]
    for t in (q, reverse_complement(q))
        anchors = Tuple{Int,Int}[]                                   # (query pos, diagonal)
        for (cp, key) in forward_kmers(t, k)
            for rp in get(refidx, key, Int32[])
                push!(anchors, (cp, mod(Int(rp) - cp, n)))
            end
        end
        isempty(anchors) && continue
        nb = cld(n, bin); counts = zeros(Int, nb)
        for (_, d) in anchors; counts[d ÷ bin + 1] += 1; end
        seen = falses(nb)
        for b0 in sortperm(counts; rev=true)
            counts[b0] < min_anchors && break
            seen[b0] && continue
            cluster = Set([b0]); seen[b0] = true
            for dir in (-1, 1)                                       # neighbouring bins (circular)
                b = b0
                while true
                    b = mod1(b + dir, nb)
                    (seen[b] || counts[b] == 0) && break
                    push!(cluster, b); seen[b] = true
                end
            end
            members = [(cp, d) for (cp, d) in anchors if (d ÷ bin + 1) in cluster]
            length(members) >= min_anchors || continue
            d0 = members[1][2]
            unwrapped = sort([d0 + mod(d - d0 + n ÷ 2, n) - n ÷ 2 for (_, d) in members])
            dm = unwrapped[cld(length(unwrapped), 2)]
            cmin = minimum(first, members); cmax = maximum(first, members) + k - 1
            push!(out, (rstart=mod1(cmin + dm, n), seq=t[cmin:cmax], anchors=length(members)))
        end
    end
    return out
end

"Pieces of `q` between N runs (n_gaps: runs of >= 5 N as (start, length))."
function split_at_ns(q::Vector{UInt8})
    out = Vector{UInt8}[]; a = 1
    for (x, len) in n_gaps(q)
        x > a && push!(out, q[a:x-1]); a = x + len
    end
    a <= length(q) && push!(out, q[a:end])
    return out
end

"Exact suffix(q) = prefix(b) overlap of >= 30 bp nearest to `want` (within +-50), or 0."
function scaffold_overlap(q::AbstractVector{UInt8}, b::AbstractVector{UInt8}, want::Int)
    best = 0
    for o in max(30, want - 50):min(length(b) - 1, length(q) - 1, want + 50)
        @view(q[end-o+1:end]) == @view(b[1:o]) || continue
        (best == 0 || abs(o - want) < abs(best - want)) && (best = o)
    end
    return best
end

"""
Reference scaffold (rescue when no circle was found). Paths are split at N runs; pieces >= 1 kb
are taken longest first unless they add < 300 new 21-mers, and placed on the ordered nearest
reference. Placements are tiled greedily (longest first, each must add >= 1 kb of uncovered
reference); in reference order a long overlap with the previous piece is trimmed to 100 bp.
Neighbours are merged where their sequences overlap exactly by >= 30 bp near the predicted
overlap (sequence from the reads); otherwise an N gap of the predicted length (>= 10) is put
between them, and a predicted overlap not found in the sequences is trimmed from the second
piece. Returns the circle (closure ref_scaffold) or nothing when the placements cover
< 80% of the reference, the circle has > `max_n` of N, or it falls outside Genome Range.
"""
function reference_scaffold(paths::Vector{Vector{UInt8}}, ref::Vector{UInt8}, c::Config; k::Int=21, max_n::Float64=0.02)
    empty!(REF_SCAFFOLD_REPORT); REF_SCAFFOLD_REJECT[] = ""
    n = length(ref)
    n >= 10_000 || (REF_SCAFFOLD_REJECT[] = "reference shorter than 10 kb"; return nothing)
    refidx = Dict{UInt64,Vector{Int32}}()
    for (p, key) in forward_kmers(ref, k)
        push!(get!(Vector{Int32}, refidx, key), Int32(p))
    end
    pieces = Vector{UInt8}[]
    for q in paths, piece in split_at_ns(q)
        length(piece) >= 1000 && push!(pieces, piece)
    end
    taken = Set{UInt64}(); chosen = Vector{UInt8}[]
    for q in sort(pieces; by=q -> -length(q))
        keys = rescue_keys(q)
        count(!in(taken), keys) >= 300 || continue
        push!(chosen, q); union!(taken, keys)
    end
    pl = NamedTuple{(:rstart,:seq,:anchors),Tuple{Int,Vector{UInt8},Int}}[]
    for q in chosen; append!(pl, ref_placements(q, refidx, n, k)); end
    isempty(pl) && (REF_SCAFFOLD_REJECT[] = "no path placed on the reference ($(length(chosen)) pieces >= 1 kb tried)"; return nothing)
    # Greedy tiling: longest placements first; one is taken only if it adds >= 1 kb of reference
    # not covered yet (near-duplicate paths add nothing and are skipped).
    sort!(pl; by=x -> (-length(x.seq), -x.anchors))
    covered = falses(n); keep = similar(pl, 0)
    for x in pl
        newbp = count(i -> !covered[mod1(x.rstart + i, n)], 0:length(x.seq)-1)
        newbp >= 1000 || continue
        push!(keep, x)
        for i in 0:length(x.seq)-1; covered[mod1(x.rstart + i, n)] = true; end
    end
    if count(covered) < 0.8n
        REF_SCAFFOLD_REJECT[] = @sprintf("placements cover %.1f%% of the reference (< 80%%; %d placements from %d pieces)",
                                         100count(covered)/n, length(keep), length(chosen))
        return nothing
    end
    junction_gap(a, b) = (d = mod(b.rstart - (a.rstart + length(a.seq)), n); d > n ÷ 2 ? d - n : d)
    # In reference order, a long overlap with the previous piece is trimmed from the next one down
    # to `margin` bp (the junction is then checked over that stretch); a piece ending inside the
    # previous one is dropped.
    margin = 100
    sort!(keep; by=x -> x.rstart)
    tiled = [keep[1]]
    for x in keep[2:end]
        d = junction_gap(tiled[end], x)
        if d < -margin
            cut = -d - margin
            cut >= length(x.seq) - 100 && continue
            x = (rstart=mod1(x.rstart + cut, n), seq=x.seq[cut+1:end], anchors=x.anchors)
        end
        push!(tiled, x)
    end
    if length(tiled) > 1                                   # closing pair: trim the end of the last piece
        d = junction_gap(tiled[end], tiled[1])
        if d < -margin
            cut = -d - margin
            if cut < length(tiled[end].seq) - 100
                y = tiled[end]; tiled[end] = (rstart=y.rstart, seq=y.seq[1:end-cut], anchors=y.anchors)
            else
                pop!(tiled)
            end
        end
    end
    keep = tiled
    function join_to!(q, a, b, idx_a, idx_b, closing)
        d = junction_gap(a, b); bseq = b.seq; how = ""; used = 0
        if d < 0
            want = -d; best = scaffold_overlap(q, bseq, want)
            if best > 0
                how = "overlap"; used = best
                closing || append!(q, @view bseq[best+1:end])
                push!(REF_SCAFFOLD_REPORT, join((length(REF_SCAFFOLD_REPORT)+1, idx_a, idx_b, mod1(a.rstart + length(a.seq) - 1, n), b.rstart, d, how, best), '\t'))
                return best
            end
            how = "N_gap_unverified_overlap"; used = min(want, length(bseq) - 1); gap = c.unknown_gap_ns > 0 ? c.unknown_gap_ns : 10
        else
            how = "N_gap"; gap = max(d, 10)
        end
        append!(q, fill(UInt8('N'), gap))
        closing || append!(q, @view bseq[used+1:end])
        push!(REF_SCAFFOLD_REPORT, join((length(REF_SCAFFOLD_REPORT)+1, idx_a, idx_b, mod1(a.rstart + length(a.seq) - 1, n), b.rstart, d, how, gap), '\t'))
        return closing ? -used : 0
    end
    q = copy(keep[1].seq)
    for i in 2:length(keep)
        join_to!(q, keep[i-1], keep[i], i-1, i, false)
    end
    # Closing junction: last placement back to the first.
    r = join_to!(q, keep[end], keep[1], length(keep), 1, true)
    if r > 0
        q = q[r+1:end]                       # verified overlap with the start: drop the duplicate
    elseif r < 0
        q = q[1-r:end]                       # unverified predicted overlap: trimmed from the first path
    end
    if !(c.genome_range[1] <= length(q) <= c.genome_range[2])
        REF_SCAFFOLD_REJECT[] = "length $(length(q)) outside Genome Range $(c.genome_range[1])-$(c.genome_range[2])"
        return nothing
    end
    nbp = count(==(UInt8('N')), q)
    if nbp > max_n * length(q)
        REF_SCAFFOLD_REJECT[] = @sprintf("%d N bp = %.2f%% of %d bp (> %.0f%%)", nbp, 100nbp/length(q), length(q), 100max_n)
        return nothing
    end
    return q
end

"""
Rescue step after graph rescue: if no circle exists and a Reference dataset is set, the nearest
reference of the linear results (not the library itself) is ordered LSC-IR-SSC-IR and the paths
are scaffolded on it. Returns (result or nothing, note).
"""
function reference_scaffold_rescue(results::Vector{AssemblyResult}, c::Config)
    (isempty(c.reference_dataset) || !c.reference_scaffold) && return nothing, ""
    sk = joinpath(c.reference_dataset, RefSketch.SKETCH_FILE); isfile(sk) || return nothing, ""
    paths = vcat([x.sequence for x in results if !x.circular], REF_RESCUE_PATHS)
    query = [q for q in paths if length(q) >= 5000]
    isempty(query) && return nothing, "Reference scaffold: no path >= 5 kb"
    hits, _ = RefSketch.nearest(RefSketch.load(sk), query, fill(false, length(query)); top=1, exclude=not_self(c))
    isempty(hits) && return nothing, "Reference scaffold: no reference"
    best = hits[1]
    best.containment >= 0.2 || return nothing, @sprintf("Reference scaffold: nearest %s too distant (containment %.3f)", best.accession, best.containment)
    ref = RefSketch.reference_sequence(c.reference_dataset, best.accession)
    item, _ = PlastomeOrder.order_plastome(ref; min_ir=c.ir_min_length, min_sc=c.ir_min_sc_length, min_identity=c.ir_min_identity)
    item === nothing || (ref = item.sequence)
    q = reference_scaffold(paths, ref, c)
    q === nothing && return nothing, "Reference scaffold: $(best.accession): not accepted: $(REF_SCAFFOLD_REJECT[])"
    joins = REF_SCAFFOLD_REPORT
    nov = count(r -> split(r, '\t')[7] == "overlap", joins)
    nbp = count(==(UInt8('N')), q)
    note = @sprintf("Reference scaffold: %s (containment %.3f); %d junction(s): %d sequence overlap(s), %d N gap(s), %d N bp; length %d (see _ref_scaffold.tsv)",
                    best.accession, best.containment, length(joins), nov, length(joins) - nov, nbp, length(q))
    return AssemblyResult("ref_scaffold_1", q, true, "closed_ref_scaffold", "reference_scaffold:" * best.accession, 0, 0, 0), note
end

"Copy of `x` whose id starts with the Project name (unchanged if it already does)."
function with_library_name(x::AssemblyResult, c::Config)
    tag = c.project_name * "_"
    startswith(x.id, tag) && return x
    return AssemblyResult(tag * x.id, x.sequence, x.circular, x.reason, x.source, x.closure_overlap,
                          x.junction_reads, x.junction_pairs, x.variants)
end

"0 if the circle, ordered LSC-IR-SSC-IR without SSC inversion, has ndhF at the SSC start (next to the first IR); else 1."
function ndhf_rank(x::AssemblyResult, c::Config)
    (x.circular && c.normalize_plastome && c.assembly_type == "chl" &&
     c.genome_range[1] <= length(x.sequence) <= c.genome_range[2]) || return 1
    item, _ = PlastomeOrder.order_plastome(x.sequence; min_ir=c.ir_min_length, min_sc=c.ir_min_sc_length,
                                           min_identity=c.ir_min_identity, orient_ssc=false)
    item === nothing && return 1
    return PlastomeOrder.ndhf_side(item.sequence, item.lsc_length, item.ir_pair.len, item.ssc_length) == :start ? 0 : 1
end

function write_results(c::Config, idx::ReadIndex, results, raw, seconds::Float64; monitor=nothing)
    progress_phase!(monitor,"Export")
    READ_STATS[] = nothing
    mkpath(c.output_path)
    prefix = joinpath(c.output_path, c.project_name)
    # Every output record is named <Project name>_<id> (Project name = library name), so
    # sequences from different libraries stay distinguishable after files are merged.
    results = AssemblyResult[with_library_name(x, c) for x in results]
    raw = AssemblyResult[with_library_name(x, c) for x in raw]
    # Among equal circles (SSC flip-flop isomers) the primary is the one that, ordered
    # LSC-IR-SSC-IR, already has ndhF at the SSC start (next to the first IR).
    pref = IdDict{Any,Int}(x => ndhf_rank(x, c) for x in results)
    results = sort(results; by=x -> (result_tier(x, c), -length(x.sequence), pref[x], x.id))
    accepted = [x for x in results if length(x.sequence) >= c.min_contig_length]
    rejected = [x for x in results if length(x.sequence) < c.min_contig_length]
    validated_count = count(x -> closure_status(x) == "validated", accepted)
    closure_counts = Dict(k => count(x -> closure_status(x) == k, accepted)
                          for k in ("validated", "ir_collapse", "read_bridge", "graph_path", "estimated_repeat", "gap_unknown", "inferred_join", "topology_inferred", "unresolved", "ref_scaffold"))
    progress_circles!(monitor,validated_count;primary_length=isempty(accepted) ? 0 : length(accepted[1].sequence),
                      primary_circular=!isempty(accepted) && accepted[1].circular)
    progress_update!(monitor;detail="validated circles=$(validated_count)",force=true)
    # Always write both, so a shorter contig is not silently lost or a stale file left behind.
    # Run-level reference use that can change the sequence: seeds (paths grown from reference
    # seeds were retained) and scaffold (a reference scaffold is among the results). The ordered
    # output adds orientation=reference... when the reference set its orientation.
    ruse = String[]
    isempty(REF_RESCUE_PATHS) || push!(ruse, "seeds")
    any(x -> closure_status(x) == "ref_scaffold", accepted) && push!(ruse, "scaffold")
    reference_use = isempty(c.reference_dataset) ? "" : " reference_use=" * (isempty(ruse) ? "none" : join(ruse, "+"))
    write_fasta(prefix*"_assembly.fasta", accepted; roles=true, extra=reference_use)
    write_fasta(prefix*"_short_contigs.fasta", rejected)
    open(prefix*"_contigs.tsv", "w") do io
        println(io, "id\tlength\ttopology\toutput\tstop_reason\tclosure_overlap\tclosure_spanning_read_pairs\tclosure_FR_pairs\tseed\tclosure_status")
        for x in results
            println(io, join((x.id, length(x.sequence), x.circular ? "circular" : "linear",
                             length(x.sequence) >= c.min_contig_length ? "assembly" : "short_contigs",
                             clean_tsv(x.reason), x.closure_overlap, x.junction_reads,
                             x.junction_pairs, clean_tsv(x.source), closure_status(x)), '\t'))
        end
    end
    open(prefix*"_closure_log.tsv", "w") do io
        println(io, join(CLOSURE_LOG_HEADER, '\t'))
        rows = lock(() -> copy(CLOSURE_LOG), CLOSURE_LOG_LOCK)
        for row in sort(rows; by=r -> (r[1], parse(Int, r[2])))
            println(io, join(row, '\t'))
        end
    end
    open(prefix*"_mate_jumps.tsv", "w") do io
        println(io, "path\tdirection\tcontig_length_after\tanchored_pairs\tpairs_linked_to_far_mate\tfar_side_bases\tnear_trimmed\tfar_trimmed")
        for row in lock(() -> copy(MATE_JUMP_LOG), MATE_JUMP_LOG_LOCK)
            println(io, row)
        end
    end
    if !isempty(raw)
        g = GRAPH_RESCUE_GRAPH[] === nothing ? assembly_graph([x.sequence for x in raw], idx, c; read_bridges=false) : GRAPH_RESCUE_GRAPH[]
        write_graph_gfa(prefix*"_graph.gfa", g)
        write_graph_nodes(prefix*"_graph_nodes.tsv", g)
        open(prefix*"_graph_read_bridges.tsv", "w") do io
            println(io, GRAPH_READ_BRIDGE_HEADER)
            rows=isempty(GRAPH_READ_BRIDGE_REPORT[]) ? g.read_bridge_report : GRAPH_READ_BRIDGE_REPORT[]
            foreach(r -> println(io,r), rows)
        end
        open(prefix*"_graph_links.tsv", "w") do io
            println(io, "node\tside\tnode_length\tterminal_run_excluded\tbest_partner:pairs\tnext_partner_pairs\tresult")
            foreach(r -> println(io, r), g.report)
            if !isempty(GRAPH_RESCUE_LINKS[])
                println(io, "# graph rescue (all passes' paths, with inferred joins if tried):")
                foreach(r -> println(io, r), GRAPH_RESCUE_LINKS[])
            end
            if !isempty(ISLAND_SEQS)
                open(prefix*"_islands.fasta", "w") do fio
                    for (name, q) in ISLAND_SEQS
                        println(fio, ">", name, " length=", length(q)); println(fio, dna_string(q))
                    end
                end
            end
            if !isempty(ISLAND_REPORT)
                println(io, "# island recovery: end\tside\twindow\tterminal_run\tanchored_reads\tmates_beyond\tnew_mates\tlinked_mates\tisland_bp\tresult")
                foreach(r -> println(io, "# ", r), ISLAND_REPORT)
            end
        end
    end
    open(prefix*"_graph_topology.tsv","w") do io
        println(io,GRAPH_TOPOLOGY_HEADER)
        foreach(row->println(io,row),GRAPH_TOPOLOGY_REPORT)
    end
    # Paths grown from reference seeds (also when a later pass replaced the results).
    open(prefix*"_reference_seed_paths.fasta", "w") do io
        for (i, q) in enumerate(REF_RESCUE_PATHS)
            println(io, ">", c.project_name, "_reference_seed_path_", i, " length=", length(q))
            for j in 1:80:length(q)
                println(io, dna_string(@view q[j:min(j+79, end)]))
            end
        end
    end
    open(prefix*"_ref_scaffold.tsv", "w") do io
        println(io, REF_SCAFFOLD_HEADER)
        foreach(r -> println(io, r), REF_SCAFFOLD_REPORT)
    end
    open(prefix*"_mate_rescue.tsv", "w") do io
        println(io,MATE_RESCUE_HEADER)
        foreach(row->println(io,row),MATE_RESCUE_REPORT)
    end
    write_fasta(prefix*"_mate_rescue.fasta",MATE_RESCUE_RESULTS)
    open(prefix*"_self_inversions.tsv", "w") do io
        println(io, "id\tlength\ttopology\tA_start\tA_end\tB_start\tB_end\tlength_AB\tN_gaps")
        for x in accepted
            gaps = join(("$(a)+$(l)" for (a, l) in n_gaps(x.sequence)), ",")
            for (a1, a2, b1, b2) in self_inversions(x.sequence, c.kmer)
                println(io, join((x.id, length(x.sequence), x.circular ? "circular" : "linear",
                                  a1, a2, b1, b2, a2 - a1 + 1, gaps), '\t'))
            end
        end
    end
    open(prefix*"_forks.tsv", "w") do io
        println(io, "path\tdirection\tcontig_length\tchoices\taction\tbranches_already_spawned\tmin_supports\tfirst_40_bases")
        for row in lock(() -> copy(FORK_LOG), FORK_LOG_LOCK)
            println(io, row)
        end
    end
    open(prefix*"_tandem_repeats.tsv", "w") do io
        println(io, "path\tdirection\tcontig_length_after\tperiod\trepeat_bases_in_path\trepeat_bases_in_reads\ttrimmed\tappended\tsupporting_pairs\tnext_best_pairs")
        for row in lock(() -> copy(TANDEM_LOG), TANDEM_LOG_LOCK)
            println(io, row)
        end
    end
    if c.coverage_trace
        open(prefix*"_coverage_trace.tsv", "w") do io
            println(io, "path\tdirection\tcontig_length\tfirst_position_pairs\texpected\tmin_support\thigh_run\tswitches\tcandidates")
            for row in lock(() -> copy(COVERAGE_TRACE), COVERAGE_TRACE_LOCK)
                println(io, row)
            end
        end
    end
    open(prefix*"_estimated_gaps.tsv", "w") do io
        println(io, "status\tcontig\tjunction_position\tmotif\tbases_added\tgap_estimate_median\tgap_min\tgap_max\tspanning_pairs\tmedian_insert\tright_overhang\tleft_overhang")
        for g in lock(() -> copy(GAP_CLOSURES), GAP_CLOSURES_LOCK)
            @printf(io, "%s\t%s\t%d\t%s\t%d\t%.1f\t%d\t%d\t%d\t%.1f\t%s\t%s\n", g.status, g.contig, g.junction, g.motif, g.change,
                    g.gap_median, g.gap_min, g.gap_max, g.spanning_pairs, g.median_insert,
                    isempty(g.right_overhang) ? "-" : g.right_overhang, isempty(g.left_overhang) ? "-" : g.left_overhang)
        end
    end
    # Repeat-length variants: primary result only (positions in its _assembly.fasta record).
    repeat_variants = [i == 1 ? repeat_length_variants(x, idx, c) : Tuple{ContigVariant,String}[]
                       for (i, x) in enumerate(accepted)]
    open(prefix*"_variants.tsv", "w") do io
        println(io, "contig\tposition\tmajor\tminor\ttype\tmajor_pairs\tminor_pairs\tinformative_pairs\tminor_fraction")
        for x in accepted, v in x.variants
            kind = v.switch ? "coverage_switch" :
                   length(v.major) == length(v.minor) ? (length(v.major) == 1 ? "SNV" : "MNV") : "indel"
            @printf(io, "%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.3f\n", x.id, v.position,
                    isempty(v.major) ? "-" : v.major, isempty(v.minor) ? "-" : v.minor, kind,
                    v.major_support, v.minor_support, v.informative,
                    v.minor_support/max(1, v.major_support+v.minor_support))
        end
        # repeat_length rows count reads, not pairs (columns kept for one table).
        for (x, rv) in zip(accepted, repeat_variants), (v, kind) in rv
            @printf(io, "%s\t%d\t%s\t%s\t%s\t%d\t%d\t%d\t%.3f\n", x.id, v.position, v.major, v.minor,
                    kind, v.major_support, v.minor_support, v.informative,
                    v.minor_support/max(1, v.major_support+v.minor_support))
        end
    end
    ordering_t0 = time()
    nearest_reference = write_nearest_reference(c, accepted, prefix)
    ordered_count = write_plastome_outputs(c,accepted,prefix;monitor=monitor,reference=nearest_reference)
    ordering_seconds = time()-ordering_t0
    annotation_t0 = time()
    annotated_count, annotation_note = annotate_plastomes(c, prefix; monitor=monitor)
    annotation_seconds = time()-annotation_t0
    open(prefix*"_run.txt", "w") do io
        println(io, "AssemblerComplete ", VERSION_STRING)
        println(io, "Julia: ", VERSION, "; threads: ", Threads.nthreads())
        println(io, "Loaded pairs: ", num_pairs(idx.store))
        println(io, "Reads shorter than minimizer window: ", idx.short_reads)
        println(io, "Minimizers: ", length(idx.postings), "; filtered: ", idx.dropped_keys,
                    "; global read-frequency cutoff: ", idx.cutoff)
        println(io, "Raw paths: ", length(raw), "; output contigs: ", length(accepted),
                    "; short contigs: ", length(rejected))
        println(io, "Validated circles: ", validated_count)
        isempty(nearest_reference.summary) || println(io, "Nearest reference: ", nearest_reference.summary)
        if !isempty(accepted)
            p = accepted[1]
            println(io, "Primary result: ", p.id, " length=", length(p.sequence), " topology=",
                    p.circular ? "circular" : "linear", " closure=", closure_status(p))
            st = plastid_read_stats(idx, p, c)
            READ_STATS[] = st
            @printf(io, "Plastid reads (approx., sharing a minimizer with %s): %d of %d (%.2f%%)\n",
                    st.id, st.reads, st.total, st.percent)
            @printf(io, "Mean plastid coverage (approx.): %.1fx\n", st.coverage)
        end
        println(io, "Circles by closure: ", join(("$(k)=$(closure_counts[k])" for k in
                    ("validated", "ir_collapse", "read_bridge", "graph_path", "estimated_repeat", "gap_unknown", "inferred_join", "topology_inferred", "unresolved", "ref_scaffold")), " "))
        for line in INITIAL_LOG
            println(io, "Initial read: ", line)
        end
        println(io, "Bubble policy: ", c.bubble_policy, "; resolved bubbles reported: ",
                sum(x -> length(x.variants), accepted; init=0), "; repeat-length variants: ",
                sum(length, repeat_variants; init=0), " (", prefix, "_variants.tsv)")
        println(io, "Plastomes ordered LSC-IR-SSC-IR: ", ordered_count)
        @printf(io, "Ordering seconds: %.2f\n", ordering_seconds)
        if !isempty(c.annotation_dir)
            println(io, "Plastomes annotated: ", annotated_count, isempty(annotation_note) ? "" : " ($annotation_note)")
            @printf(io, "Annotation seconds: %.2f\n", annotation_seconds)
        end
        @printf(io, "Elapsed seconds: %.2f\n", seconds)
        println(io, "Reads: first N input pairs; no random subsampling or automatic insert estimation.")
        println(io, "Alignment: ungapped; long repeats and structural alternatives can remain unresolved.")
        for field in fieldnames(Config)
            println(io, field, " = ", repr(getfield(c, field)))
        end
    end
    if c.debug_mode
        write_fasta(prefix*"_all_paths.fasta", raw)
    end
    if isempty(accepted)
        progress_breakline!(monitor)
        @warn "No contig reached Min contig length; inspect short_contigs FASTA and contigs TSV"
    end
    return prefix*"_assembly.fasta"
end

# ------------------------------------------------------------------------------------------
# Assembly graph of the paths (16.30): compacted de Bruijn graph (k = 63) of all path
# sequences, node multiplicity from read coverage, and circular walks that use every node as
# many times as its multiplicity (as GetOrganelle does on SPAdes graphs, here without SPAdes).
# A path that walked a single-copy region twice, or two paths that each cover part of the
# genome, collapse onto the same nodes; the IR becomes one node of multiplicity 2.
# ------------------------------------------------------------------------------------------
const GRAPH_K = 63

"2-bit forward code of a k-mer (UInt128, k <= 63)."
@inline function code128(seq::AbstractVector{UInt8}, i::Int, k::Int)
    v = UInt128(0)
    @inbounds for j in i:i+k-1
        e = ENCODE[Int(seq[j]) + 1]
        e == 255 && return typemax(UInt128)
        v = (v << 2) | UInt128(e)
    end
    return v
end

function rc128(x::UInt128, k::Int)
    r = UInt128(0)
    for _ in 1:k
        r = (r << 2) | (UInt128(3) - (x & UInt128(3)))
        x >>= 2
    end
    return r
end

@inline canon128(x::UInt128, k::Int) = min(x, rc128(x, k))

function decode128(x::UInt128, k::Int)
    out = Vector{UInt8}(undef, k)
    for j in k:-1:1
        out[j] = DNA[Int(x & UInt128(3)) + 1]; x >>= 2
    end
    return out
end

struct GraphReadBridge
    trim_from::Int
    trim_to::Int
    sequence::Vector{UInt8}       # complete inter-flank sequence; Ns if unresolved
    supporting_pairs::Int
    unresolved::Bool
    estimated::Bool               # estimated_repeat: few or discordant spanning reads (length a best guess)
end
GraphReadBridge(trim_from, trim_to, sequence, pairs, unresolved) =
    GraphReadBridge(trim_from, trim_to, sequence, pairs, unresolved, false)

const PANEL_K = 21
const PANEL_CACHE = Dict{String,Set{UInt64}}()

"Canonical 21-mers of all plastomes in a GenBank/FASTA file or a folder of them (cached)."
function load_plastid_panel(path::AbstractString)
    haskey(PANEL_CACHE, path) && return PANEL_CACHE[path]
    files = isdir(path) ? [joinpath(path, f) for f in readdir(path) if occursin(r"\.(gb|gbk|gbff|fa|fasta|fna)$"i, f)] : [path]
    keys = Set{UInt64}()
    for f in files
        txt = read(f, String)
        seqs = if occursin(r"^LOCUS"m, txt)
            [uppercase(replace(split(part, "//")[1], r"[^A-Za-z]" => "")) for part in split(txt, "ORIGIN")[2:end]]
        else
            [uppercase(replace(join(split(r, '\n')[2:end]), r"[^A-Za-z]" => "")) for r in split(txt, '>')[2:end]]
        end
        for q in seqs
            isempty(q) || union!(keys, all_kmer_keys(Vector{UInt8}(codeunits(q)), PANEL_K))
        end
    end
    PANEL_CACHE[path] = keys
    return keys
end

"Share of a sequence's canonical 21-mers present in the plastid panel."
function plastid_fraction(seq::AbstractVector{UInt8}, panel::Set{UInt64})
    ks = all_kmer_keys(seq, PANEL_K)
    return isempty(ks) ? 0.0 : count(in(panel), ks) / length(ks)
end

mutable struct PathGraph
    k::Int
    kset::Set{UInt128}                                   # canonical k-mers
    paths::Vector{Vector{UInt128}}                        # per node: oriented k-mers
    seqs::Vector{Vector{UInt8}}                           # per node: sequence (forward)
    kinfo::Dict{UInt128,Tuple{Int,Int}}                   # canonical -> (node, position)
    gaps::Dict{Tuple{Int,Bool},Set{Tuple{Int,Bool,Int}}}  # (node, plus) -> (node, plus, N count)
    cov::Vector{Float64}
    base::Float64
    mult::Vector{Int}
    report::Vector{String}                                # dead-end linking, one row per end
    inferred::Set{NTuple{4,Any}}                          # (node, plus, node, plus) of inferred joins
    lo::Vector{Int}                                       # allowed uses of each node in a walk
    hi::Vector{Int}
    triples::Set{NTuple{6,Any}}                           # node triples seen contiguously in the paths
    read_bridges::Dict{NTuple{4,Any},GraphReadBridge}
    read_bridge_report::Vector{String}
    plastid::Vector{Float64}                              # share of a node's 21-mers in the plastid panel (NaN: no panel)
end

function graph_succ(g::PathGraph, x::UInt128)
    mask = (UInt128(1) << (2g.k)) - UInt128(1)
    out = UInt128[]
    for b in 0:3
        y = ((x << 2) | UInt128(b)) & mask
        canon128(y, g.k) in g.kset && push!(out, y)
    end
    return out
end

function graph_pred(g::PathGraph, x::UInt128)
    out = UInt128[]
    for b in 0:3
        y = (x >> 2) | (UInt128(b) << (2(g.k-1)))
        canon128(y, g.k) in g.kset && push!(out, y)
    end
    return out
end

"Node, position and orientation (true = as stored) of an oriented k-mer of the graph."
function graph_locate(g::PathGraph, y::UInt128)
    n, pos = g.kinfo[canon128(y, g.k)]
    return n, pos, y == g.paths[n][pos]
end

graph_tail(g::PathGraph, n::Int, plus::Bool) = plus ? g.paths[n][end] : rc128(g.paths[n][1], g.k)

"Oriented nodes following (n, plus) through a (k-1)-overlap."
function graph_nexts(g::PathGraph, n::Int, plus::Bool)
    out = Tuple{Int,Bool}[]
    for y in graph_succ(g, graph_tail(g, n, plus))
        v, pos, o = graph_locate(g, y)
        ((o && pos == 1) || (!o && pos == length(g.paths[v]))) && push!(out, (v, o))
    end
    return out
end

graph_prevs(g::PathGraph, n::Int, plus::Bool) = [(v, !o) for (v, o) in graph_nexts(g, n, !plus)]

"Normal (gap = -1) and N-gap successors."
graph_all_nexts(g::PathGraph, n::Int, plus::Bool) =
    vcat([(v, o, -1) for (v, o) in graph_nexts(g, n, plus)], sort(collect(get(g.gaps, (n, plus), Set{Tuple{Int,Bool,Int}}()))))

graph_seq(g::PathGraph, n::Int, plus::Bool) = plus ? g.seqs[n] : reverse_complement(g.seqs[n])

"""
Build the graph from path sequences. `covfun(seq)` gives a node's read coverage (any unit);
multiplicity = coverage / single-copy coverage (length-weighted median over nodes >= 500 bp),
nodes under 0.25x dropped, others at least 1; of parallel nodes (same neighbours on both
sides: SNV/indel alleles, mixture haplotypes) only the best covered is kept. N runs inside
a path become gap edges between the nodes on both sides.
"""
function PathGraph(seqs::Vector{Vector{UInt8}}, covfun; k::Int=GRAPH_K, panel=nothing)
    kset = Set{UInt128}()
    for s in seqs, i in 1:length(s)-k+1
        v = code128(s, i, k)
        v == typemax(UInt128) || push!(kset, canon128(v, k))
    end
    g = PathGraph(k, kset, Vector{UInt128}[], Vector{UInt8}[], Dict{UInt128,Tuple{Int,Int}}(),
                  Dict{Tuple{Int,Bool},Set{Tuple{Int,Bool,Int}}}(), Float64[], 0.0, Int[], String[],
                  Set{NTuple{4,Any}}(), Int[], Int[], Set{NTuple{6,Any}}(), Dict{NTuple{4,Any},GraphReadBridge}(), String[],
                  Float64[])
    visited = Set{UInt128}()
    for c0 in sort!(collect(kset))
        c0 in visited && continue
        path = [c0]; seen = Set{UInt128}([c0])
        cur = c0
        while true
            nx = graph_succ(g, cur); length(nx) == 1 || break
            y = nx[1]; length(graph_pred(g, y)) == 1 || break
            cy = canon128(y, k); (cy in seen || cy in visited) && break
            push!(path, y); push!(seen, cy); cur = y
        end
        cur = c0
        while true
            pv = graph_pred(g, cur); length(pv) == 1 || break
            y = pv[1]; length(graph_succ(g, y)) == 1 || break
            cy = canon128(y, k); (cy in seen || cy in visited) && break
            pushfirst!(path, y); push!(seen, cy); cur = y
        end
        n = length(g.paths) + 1
        for (pos, x) in enumerate(path)
            cx = canon128(x, k); push!(visited, cx); g.kinfo[cx] = (n, pos)
        end
        seq = decode128(path[1], k)
        for x in @view path[2:end]
            push!(seq, DNA[Int(x & UInt128(3)) + 1])
        end
        push!(g.paths, path); push!(g.seqs, seq)
    end
    # gap edges from N runs
    for s in seqs
        parts = UnitRange{Int}[]; i = 1
        while i <= length(s)
            if s[i] == UInt8('N')
                i += 1
            else
                j = i
                while j < length(s) && s[j+1] != UInt8('N'); j += 1; end
                push!(parts, i:j); i = j + 1
            end
        end
        for t in 1:length(parts)-1
            a, b = parts[t], parts[t+1]
            (length(a) >= k && length(b) >= k) || continue
            ya = code128(s, last(a)-k+1, k); yb = code128(s, first(b), k)
            (ya == typemax(UInt128) || yb == typemax(UInt128)) && continue
            u, pu, ou = graph_locate(g, ya); v, pv, ov = graph_locate(g, yb)
            ((ou && pu == length(g.paths[u])) || (!ou && pu == 1)) || continue
            ((ov && pv == 1) || (!ov && pv == length(g.paths[v]))) || continue
            gap = first(b) - last(a) - 1
            push!(get!(Set{Tuple{Int,Bool,Int}}, g.gaps, (u, ou)), (v, ov, gap))
            push!(get!(Set{Tuple{Int,Bool,Int}}, g.gaps, (v, !ov)), (u, !ou, gap))
        end
    end
    # coverage and multiplicity
    g.cov = [Float64(covfun(s)) for s in g.seqs]
    L = length.(g.seqs)
    # Plastid label (GetOrganelle-like, by k-mers instead of BLAST): nodes >= 500 bp with < 0.5%
    # of their 21-mers in the plastid panel are not plastid (e.g. a bacterial genome entered
    # through a shared rRNA) and are excluded from coverage baseline and walks.
    g.plastid = panel === nothing ? fill(NaN, length(L)) : [plastid_fraction(q, panel) for q in g.seqs]
    foreign = [!isnan(g.plastid[i]) && L[i] >= 500 && g.plastid[i] < 0.005 for i in eachindex(L)]
    pool = [(g.cov[i], L[i]) for i in eachindex(L) if L[i] >= 500 && !foreign[i]]
    isempty(pool) && (pool = [(g.cov[i], L[i]) for i in eachindex(L)])
    sort!(pool); total = sum(last, pool; init=0); acc = 0; g.base = isempty(pool) ? 0.0 : pool[end][1]
    for (cv, l) in pool
        acc += l
        if acc >= total / 2
            g.base = cv; break
        end
    end
    g.mult = [g.base <= 0 || foreign[i] ? 0 : max(1, round(Int, g.cov[i] / g.base)) for i in eachindex(g.cov)]
    low = [g.base <= 0 || cv < 0.25g.base for cv in g.cov]
    # parallel short nodes (< 1 kb: SNV/indel alleles, mixture haplotypes): keep the best
    # covered. Long nodes with the same neighbours are not alleles (e.g. both single-copy
    # stretches next to a collapsed hairpin) and all stay.
    groups = Dict{Any,Vector{Int}}()
    for n in eachindex(g.seqs)
        g.mult[n] == 0 && continue
        ka = (sort(graph_prevs(g, n, true)), sort(graph_nexts(g, n, true)))
        kb = (sort(graph_prevs(g, n, false)), sort(graph_nexts(g, n, false)))
        key = min(ka, kb)
        (isempty(key[1]) || isempty(key[2])) && continue
        push!(get!(Vector{Int}, groups, key), n)
    end
    for ns0 in values(groups)
        ns = [i for i in ns0 if length(g.seqs[i]) < 1000]
        length(ns) < 2 && continue
        best = ns[argmax([g.cov[i] for i in ns])]
        for i in ns
            i == best || (g.mult[i] = 0)
        end
    end
    # Low-coverage nodes (< 0.25x) are dropped only as dead ends; a low-coverage node that
    # connects two parts of the graph (AT-rich, poorly sequenced) is genome and stays.
    changed = true
    while changed
        changed = false
        for n in eachindex(g.seqs)
            (g.mult[n] > 0 && low[n]) || continue
            a = any(x -> g.mult[x[1]] > 0, graph_nexts(g, n, true))
            b = any(x -> g.mult[x[1]] > 0, graph_nexts(g, n, false))
            if !(a && b)
                g.mult[n] = 0; changed = true
            end
        end
    end
    # Allowed uses in a walk. Coverage of short nodes is noisy (a collapsed hairpin is used
    # twice): 1 to max(2, mult+1). Long nodes: the integers within +-0.25 of their relative
    # coverage (IR at 1.34x may be used once or twice; structure and paths decide).
    # Low-coverage nodes that stay (connectors) are optional (0 uses allowed): needed when they
    # are the only way through (an AT-rich stretch), skipped when they are error side loops.
    for n in eachindex(g.seqs)
        m = g.mult[n]; rel = g.base > 0 ? g.cov[n] / g.base : 0.0
        if m == 0
            push!(g.lo, 0); push!(g.hi, 0)
        elseif low[n]
            push!(g.lo, 0); push!(g.hi, length(g.seqs[n]) < 500 ? 2 : 1)
        elseif length(g.seqs[n]) < 500
            push!(g.lo, 1); push!(g.hi, max(2, m + 1))
        else
            lo = max(1, floor(Int, rel + 0.25)); push!(g.lo, lo); push!(g.hi, max(lo, ceil(Int, rel - 0.25)))
        end
    end
    # Node triples seen contiguously in the input paths: walks are ranked by them.
    for s in seqs
        visits = Tuple{Int,Bool}[]
        for i in 1:length(s)-k+1
            y = code128(s, i, k); y == typemax(UInt128) && (empty!(visits); continue)
            n, _, o = graph_locate(g, y)
            (isempty(visits) || visits[end] != (n, o)) && push!(visits, (n, o))
            if length(visits) >= 3
                a, b, cc = visits[end-2], visits[end-1], visits[end]
                push!(g.triples, (a..., b..., cc...))
                push!(g.triples, (cc[1], !cc[2], b[1], !b[2], a[1], !a[2]))
            end
        end
    end
    return g
end

"""
Circular walks using every kept node between `lo` and `hi` times, starting at the longest
single-copy node; depth-first with a step budget; the walks best supported by node triples
seen in the paths are kept (up to `maxsol`). Returns the circle sequences (gap edges filled
with that many Ns, min 1).
"""
function graph_circles(g::PathGraph; maxsol::Int=4, budget::Int=200_000)
    return [walk_sequence(g, w) for w in graph_walks(g; maxsol=maxsol, budget=budget)]
end

"Circular walks of graph_circles, as lists of (node, plus, gap of the edge leading to it)."
function graph_walks(g::PathGraph; maxsol::Int=4, budget::Int=200_000, accept_walk=w->true, stats=nothing, prefer_coverage::Bool=false)
    keep = [i for i in eachindex(g.mult) if g.mult[i] > 0]
    isempty(keep) && return Vector{Tuple{Int,Bool,Int}}[]
    start = keep[argmax([(g.lo[i] == 1 && g.hi[i] == 1, g.lo[i] > 0, length(g.seqs[i])) for i in keep])]
    used = zeros(Int, length(g.seqs)); used[start] = 1
    sols = Vector{Vector{Tuple{Int,Bool,Int}}}()
    steps = Ref(0); closures = Ref(0); rejected = Ref(0)
    walk = Tuple{Int,Bool,Int}[(start, true, -1)]
    function dfs(n::Int, plus::Bool)
        (length(sols) >= 4maxsol || steps[] > budget) && return
        steps[] += 1
        nexts = graph_all_nexts(g, n, plus)
        # Explore observed bridges before optional repeat walks. No change to
        # traversal order when the new evidence is absent.
        isempty(g.read_bridges) || sort!(nexts;by=x ->
            (haskey(g.read_bridges,(n,plus,x[1],x[2])) ? 0 : 1,x[1],x[2],x[3]))
        prefer_coverage && sort!(nexts;by=x -> (haskey(g.read_bridges,(n,plus,x[1],x[2])) ? 0 : 1,-g.cov[x[1]],x[1],x[2],x[3]))
        if all(i -> used[i] >= g.lo[i], keep)
            for (v, o, gap) in nexts
                if v == start && o
                    candidate=vcat(copy(walk),[(start,true,gap)]);closures[]+=1
                    if accept_walk(candidate);push!(sols,candidate);else;rejected[]+=1;end
                    break
                end
            end
        end
        for (v, o, gap) in nexts
            (v != start && used[v] < g.hi[v]) || continue
            used[v] += 1; push!(walk, (v, o, gap))
            dfs(v, o)
            pop!(walk); used[v] -= 1
        end
    end
    dfs(start, true)
    stats===nothing || merge!(stats,Dict("steps"=>steps[],"budget"=>budget,"budget_hit"=>Int(steps[]>budget),"closed_candidates"=>closures[],"rejected"=>rejected[],"accepted"=>length(sols)))
    # Ranking of closed walks (larger is better), lexicographic:
    #   1. fewer inferred joins (adjacency from a shared motif only, no data);
    #   2. fewer estimated repeat bridges (length from few or discordant spanning pairs);
    #   3. more observed read bridges (independent adjacencies, not repeated traversals);
    #   4. stronger weakest observed bridge (minimum spanning pairs; no bridge = no weak link);
    #   5. more node triples seen contiguously in the paths.
    # Each is counted once per edge, so a weak link is not hidden behind strong ones.
    function score(w)
        v = [(x[1], x[2]) for x in w[1:end-1]]
        cyc = vcat(v, v[1:min(2, length(v))])
        triples=count(i -> (cyc[i]..., cyc[i+1]..., cyc[i+2]...) in g.triples, 1:length(cyc)-2)
        observed=Set{NTuple{4,Any}}(); estimated=Set{NTuple{4,Any}}(); inferred=Set{NTuple{4,Any}}()
        weakest=typemax(Int)
        for i in 1:length(w)-1
            u,o,_=w[i]; v,p,_=w[i+1]; key=(u,o,v,p); edge=min(key,(v,!p,u,!o))
            b=get(g.read_bridges,key,nothing)
            if b!==nothing
                if b.estimated
                    push!(estimated,edge)
                else
                    push!(observed,edge); weakest=min(weakest,b.supporting_pairs)
                end
            end
            key in g.inferred && push!(inferred,edge)
        end
        return (-length(inferred),-length(estimated),length(observed),weakest,triples)
    end
    if !isempty(sols)
        best = maximum(score.(sols))
        sols = [w for w in sols if score(w) == best]
    end
    return first(sols, maxsol)
end

"Sequence of a circular walk (gap edges filled with that many Ns, min 1)."
function walk_sequence(g::PathGraph, w::Vector{Tuple{Int,Bool,Int}})
    # Trims are node-local, not additive guessed overlaps. Check both ends of
    # every visit before concatenation (including the circular first node).
    N=length(w)-1
    for i in 1:N
        u=w[i]; prev=w[i==1 ? N : i-1]; nxt=w[i+1]
        incoming=get(g.read_bridges,(prev[1],prev[2],u[1],u[2]),nothing)
        outgoing=get(g.read_bridges,(u[1],u[2],nxt[1],nxt[2]),nothing)
        left=incoming===nothing ? 0 : incoming.trim_to
        right=outgoing===nothing ? 0 : outgoing.trim_from
        left+right < length(g.seqs[u[1]]) || return UInt8[]
    end
    q=copy(graph_seq(g,w[1][1],w[1][2]))
    for i in 2:length(w)
        u=w[i-1]; v,o,gap=w[i]; closing=i==length(w)
        bridge=get(g.read_bridges,(u[1],u[2],v,o),nothing)
        if bridge!==nothing
            length(q)>=bridge.trim_from || return UInt8[]
            resize!(q,length(q)-bridge.trim_from)
            append!(q,bridge.sequence)
            if closing
                q=q[bridge.trim_to+1:end]
            else
                append!(q,@view graph_seq(g,v,o)[bridge.trim_to+1:end])
            end
        elseif gap<0
            if closing
                length(q)>=g.k-1 || return UInt8[]
                resize!(q,length(q)-(g.k-1))
            else
                append!(q,@view graph_seq(g,v,o)[g.k:end])
            end
        else
            append!(q,fill(UInt8('N'),max(gap,1)))
            closing || append!(q,graph_seq(g,v,o))
        end
    end
    return q
end

"Does a circular walk use an inferred join?"
walk_inferred(g::PathGraph, w) = any(i -> (w[i][1], w[i][2], w[i+1][1], w[i+1][2]) in g.inferred, 1:length(w)-1)

"Does a circular walk use an estimated repeat bridge?"
function walk_estimated(g::PathGraph, w)
    for i in 1:length(w)-1
        b = get(g.read_bridges, (w[i][1], w[i][2], w[i+1][1], w[i+1][2]), nothing)
        b !== nothing && b.estimated && return true
    end
    return false
end

"Motif of the low-complexity run ending within the last `slack` bases (rotation and strand canonical), or empty."
function terminal_motif(q::AbstractVector{UInt8}; slack::Int=30)
    L = length(q)
    for kk in 0:slack
        L - kk > 0 || break
        p, r = periodic_suffix(@view q[1:L-kk])
        if p > 0
            m = q[L-kk-p+1:L-kk]
            rots = [vcat(m[i:end], m[1:i-1]) for i in 1:p]
            mr = reverse_complement(m)
            append!(rots, [vcat(mr[i:end], mr[1:i-1]) for i in 1:p])
            return minimum(dna_string.(rots))
        end
    end
    return ""
end

"""
Length of the longest suffix of `q` that is a run of `motif` (any rotation, either strand).
An end that stops right at the start of a microsatellite shows only a few units of it.
"""
function motif_suffix(q::AbstractVector{UInt8}, motif::String)
    p = length(motif); best = 0
    m = Vector{UInt8}(codeunits(motif))
    for mm in (m, reverse_complement(m)), i in 1:p
        rot = vcat(mm[i:end], mm[1:i-1])
        j = 0
        while j < length(q) && q[end-j] == rot[p - (j % p)]
            j += 1
        end
        best = max(best, j)
    end
    return best
end

"""
Inferred joins: dead ends that both stop in a run of the same microsatellite or homopolymer
motif (e.g. (TA)n) are joined by Unknown gap Ns when no read or pair connects them, because
the repeat is longer than the inserts. Weaker than every observed closure: the adjacency is
inferred from the shared motif and the graph, not seen in the data (closure inferred_join).
"""
function infer_joins!(g::PathGraph, c::Config)
    ends = open_ends(g)
    # An end with a detected run (>= 8 bases) sets the motif; any other end that ends in that
    # run (detected) or in >= 2 units of it (stopped at its start) may be joined to it.
    motifs = Dict(e => terminal_motif(graph_seq(g, e[1], e[2])) for e in ends)
    pairs = Set{Tuple{Tuple{Int,Bool},Tuple{Int,Bool}}}()
    for e in ends
        key = motifs[e]; isempty(key) && continue
        for f in ends
            f == e && continue
            ok = motifs[f] == key || motif_suffix(graph_seq(g, f[1], f[2]), key) >= max(4, 2length(key))
            ok && push!(pairs, e < f ? (e, f) : (f, e))
        end
    end
    for (e, f) in pairs
        e[1] == f[1] && length(g.seqs[e[1]]) < 10_000 && continue
        push!(get!(Set{Tuple{Int,Bool,Int}}, g.gaps, (e[1], e[2])), (f[1], !f[2], c.unknown_gap_ns))
        push!(get!(Set{Tuple{Int,Bool,Int}}, g.gaps, (f[1], f[2])), (e[1], !e[2], c.unknown_gap_ns))
        push!(g.inferred, (e[1], e[2], f[1], !f[2])); push!(g.inferred, (f[1], f[2], e[1], !e[2]))
        push!(g.report, join(("n$(e[1])", e[2] ? "end" : "start", length(g.seqs[e[1]]), 0,
                              "n$(f[1])$(f[2] ? "+" : "-"):0", 0, "inferred_join"), '\t'))
    end
    return g
end

"Read coverage of a sequence: median number of indexed reads per minimizer key it contains."
function node_coverage(seq::AbstractVector{UInt8}, idx::ReadIndex, c::Config)
    counts = Int[]
    filtered = !isempty(FOREIGN_READS[])
    for key in all_kmer_keys(seq, c.kmer)
        bucket = get(idx.postings, key, nothing)
        bucket === nothing && continue
        push!(counts, filtered ? count(id -> !excluded_read(UInt32(id),c), bucket) : length(bucket))
    end
    return isempty(counts) ? 0.0 : Float64(median(counts))
end

"""
Open ends: sides of required nodes (lo > 0) from which no other required node follows, by
overlap or gap. Optional neighbours (low-coverage side loops) and the node's own reverse
(a hairpin at its end) do not count: the path does not really continue there.
"""
function open_ends(g::PathGraph)
    out = Tuple{Int,Bool}[]
    for n in eachindex(g.seqs)
        g.lo[n] > 0 || continue
        for s in (true, false)
            any(x -> x[1] != n && g.lo[x[1]] > 0, graph_all_nexts(g, n, s)) || push!(out, (n, s))
        end
    end
    return out
end

const ISLAND_REPORT = String[]
const ISLAND_SEQS = Tuple{String,Vector{UInt8}}[]

"Graph-rescue state written at export (graph, its links and bridges, islands, report line)."
graph_rescue_state() = (graph=GRAPH_RESCUE_GRAPH[], report=GRAPH_RESCUE_REPORT[], links=copy(GRAPH_RESCUE_LINKS[]),
                        bridges=copy(GRAPH_READ_BRIDGE_REPORT[]), islands=copy(ISLAND_REPORT),
                        island_seqs=copy(ISLAND_SEQS), mate_last=copy(MATE_RESCUE_LAST_RESULTS[]))

"Restores a graph_rescue_state(), so the exported GFA and evidence match the restored sequences."
function restore_graph_rescue_state!(st)
    GRAPH_RESCUE_GRAPH[] = st.graph; GRAPH_RESCUE_REPORT[] = st.report; GRAPH_RESCUE_LINKS[] = copy(st.links)
    GRAPH_READ_BRIDGE_REPORT[] = copy(st.bridges)
    empty!(ISLAND_REPORT); append!(ISLAND_REPORT, st.islands)
    empty!(ISLAND_SEQS); append!(ISLAND_SEQS, st.island_seqs)
    MATE_RESCUE_LAST_RESULTS[] = copy(st.mate_last)
    return nothing
end

"Logs written by assemble() (each call empties them first)."
logs_state() = (closure=copy(CLOSURE_LOG), gaps=copy(GAP_CLOSURES), coverage=copy(COVERAGE_TRACE),
                tandem=copy(TANDEM_LOG), forks=copy(FORK_LOG), jumps=copy(MATE_JUMP_LOG))

"Restores a logs_state()."
function restore_logs!(st)
    lock(() -> (empty!(CLOSURE_LOG); append!(CLOSURE_LOG, st.closure)), CLOSURE_LOG_LOCK)
    lock(() -> (empty!(GAP_CLOSURES); append!(GAP_CLOSURES, st.gaps)), GAP_CLOSURES_LOCK)
    lock(() -> (empty!(COVERAGE_TRACE); append!(COVERAGE_TRACE, st.coverage)), COVERAGE_TRACE_LOCK)
    lock(() -> (empty!(TANDEM_LOG); append!(TANDEM_LOG, st.tandem)), TANDEM_LOG_LOCK)
    lock(() -> (empty!(FORK_LOG); append!(FORK_LOG, st.forks)), FORK_LOG_LOCK)
    lock(() -> (empty!(MATE_JUMP_LOG); append!(MATE_JUMP_LOG, st.jumps)), MATE_JUMP_LOG_LOCK)
    return nothing
end

"Everything the export reads besides the results: graph-rescue state and assemble() logs."
export_state() = (graph=graph_rescue_state(), logs=logs_state())

"Restores an export_state(): reports, GFA and logs then match the results restored with it."
function restore_export_state!(st)
    restore_graph_rescue_state!(st.graph); restore_logs!(st.logs)
    return nothing
end

"""
Islands: sequence lying between two obstacles that no read crosses (e.g. (TA)27 and a poly-A
1 kb apart) is never reached by extension from outside and has no seed inside. From each open
end, reads placed wholly and uniquely in the insert-length window before the end's own
low-complexity run, pointing to the end, whose mates are not in the graph, give mates from
beyond the obstacle. With >= Min coverage such pairs whose mate is new sequence (>= 30
complex 31-mers absent from the graph on either strand), the mate reaching furthest into new
sequence is corrected like a seed and extended by local consensus: back towards the obstacle
(up to 2 insert lengths) and forward into the island
(up to 5 kb). The seed is the mate with most new complex 31-mers; the island is kept as an extra
path if forward read consensus grew it by >= 200 bp.
"""
function recover_islands(g::PathGraph, idx::ReadIndex, c::Config)
    c.paired || return Vector{UInt8}[]
    k = c.kmer
    _, high = insert_bounds(c)
    known = Set{UInt64}()
    for q in g.seqs
        union!(known, forward_kmer_set(q, k)); union!(known, forward_kmer_set(reverse_complement(q), k))
    end
    islands = Vector{UInt8}[]
    empty!(ISLAND_REPORT)
    for (n, s) in open_ends(g)
        q = graph_seq(g, n, s)
        r = tail_run_trim(q); r < length(q) - 100 || (r = 0)
        window = q[max(1, end-r-high+1):end-r]
        row(result, a=0, b=0, cnov=0, l=0, len=0) = push!(ISLAND_REPORT, join(("n$n", s ? "end" : "start",
            length(window), r, a, b, cnov, l, len, result), '\t'))
        length(window) >= 2*idx.max_read_length || (row("window_too_short"); continue)
        text = dna_string(window); anchors = ExactAnchorCache(text)
        ids, saturated = candidate_read_ids(window, idx, c)
        saturated && (row("candidates_saturated"); continue)
        mates = Dict{UInt32,Vector{UInt8}}(); newk = Dict{UInt32,Int}()
        anchored = 0; beyond = 0
        for id in ids
            read = get_read(idx.store, id)
            pl = read_placements(text, read, c; anchors=anchors)
            (length(pl) == 1 && !only(pl).reverse) || continue
            anchored += 1
            mate = get_read(idx.store, id ⊻ UInt32(1))
            isempty(read_placements(text, mate, c; anchors=anchors)) || continue
            beyond += 1
            m = reverse_complement(mate)
            # New sequence = >= 30 complex 31-mers absent from the graph; repeat k-mers do not
            # count (a mate straddling the (TA)n and the island is new despite its TA part).
            km = [v for v in forward_kmer_set(m, k) if !low_complexity_kmer(v, k)]
            nnew = count(v -> !(v in known), km)
            nnew >= 30 || continue
            mates[pair_id(id)] = m; newk[pair_id(id)] = nnew
        end
        length(mates) >= c.min_coverage || (row("too_few_new_mates", anchored, beyond, length(mates)); continue)
        groups=mate_seed_groups(mates,known,c;k=k)
        isempty(groups) && (row("no_coherent_mate_cluster",anchored,beyond,length(mates));continue)
        group=first(groups);best=group.best;links=length(group.ids)-1
        far = correct_seed_read(mates[best], idx, c)
        # Back towards the obstacle only far enough to reach known sequence (overlap links);
        # forward into the island up to 5 kb, which must grow by >= 200 bp of read consensus.
        back = reverse_complement(local_extend(reverse_complement(far), idx, c; maxlen=2high))
        island = local_extend(back, idx, c; maxlen=5000)
        length(island) - length(back) >= 200 ||
            (row("island_not_extendable", anchored, beyond, length(mates), links + 1, length(island)); continue)
        row("island", anchored, beyond, length(mates), links + 1, length(island))
        push!(ISLAND_SEQS, ("island_n$(n)_$(s ? "end" : "start")_seed", far))
        push!(ISLAND_SEQS, ("island_n$(n)_$(s ? "end" : "start")", island))
        push!(islands, island)
        union!(known, forward_kmer_set(island, k)); union!(known, forward_kmer_set(reverse_complement(island), k))
    end
    return islands
end

"""
Scaffold links between dead ends of the graph (kept nodes with no kept neighbour on one side),
as SPAdes/GetOrganelle scaffolds do: a read placed wholly and uniquely in the last insert-length
of one dead end, pointing to it, whose mate is placed wholly and uniquely in the last
insert-length of another dead end, pointing to that one, shows that the two ends face each
other across a stretch no read crosses (e.g. a long (TA)n). Ends are linked when they are each
other's best partner with >= Min paired support distinct pairs and at least twice the pairs of
any other partner; the link is a gap edge of Unknown gap Ns (its length is not observed).
"""
function link_dead_ends!(g::PathGraph, idx::ReadIndex, c::Config)
    (c.paired && c.unknown_gap_ns > 0) || return g
    _, high = insert_bounds(c)
    ends = open_ends(g)
    length(ends) >= 2 || return g
    # The insert-length window before the end, without the end's own low-complexity run: a
    # node that stops in a long (TA)n would otherwise offer only repeat reads, none unique.
    term = Dict{Tuple{Int,Bool},String}(); runs = Dict{Tuple{Int,Bool},Int}()
    for e in ends
        q = graph_seq(g, e[1], e[2])
        r = tail_run_trim(q)
        r < length(q) - 100 || (r = 0)
        runs[e] = r
        term[e] = dna_string(q[max(1, end-r-high+1):end-r])
    end
    anchors = Dict(e => ExactAnchorCache(t) for (e, t) in term)
    links = Dict{Tuple{Tuple{Int,Bool},Tuple{Int,Bool}},Set{UInt32}}()
    for e in ends
        ids, saturated = candidate_read_ids(Vector{UInt8}(codeunits(term[e])), idx, c)
        saturated && continue
        for id in ids
            read = get_read(idx.store, id)
            pl = read_placements(term[e], read, c; anchors=anchors[e])
            (length(pl) == 1 && !only(pl).reverse) || continue
            mate = get_read(idx.store, id ⊻ UInt32(1))
            for f in ends
                f == e && continue
                f[1] == e[1] && length(g.seqs[e[1]]) < 10_000 && continue
                pm = read_placements(term[f], mate, c; anchors=anchors[f])
                (length(pm) == 1 && !only(pm).reverse) || continue
                push!(get!(Set{UInt32}, links, (e, f)), pair_id(id))
            end
        end
    end
    support(e, f) = length(union(get(links, (e, f), Set{UInt32}()), get(links, (f, e), Set{UInt32}())))
    function partner(e)
        ranked = sort([(f, support(e, f)) for f in ends if f != e]; by=x -> -x[2])
        isempty(ranked) && return nothing
        best, n1 = ranked[1]; n2 = length(ranked) > 1 ? ranked[2][2] : 0
        (n1 >= c.min_paired_support && n1 >= 2n2) ? best : nothing
    end
    done = Set{Tuple{Int,Bool}}()
    for e in ends
        ranked = sort([(f, support(e, f)) for f in ends if f != e]; by=x -> -x[2])
        best = isempty(ranked) ? "-" : "n$(ranked[1][1][1])$(ranked[1][1][2] ? "+" : "-"):$(ranked[1][2])"
        second = length(ranked) > 1 ? ranked[2][2] : 0
        e in done && continue
        f = partner(e)
        linked = f !== nothing && partner(f) == e
        push!(g.report, join(("n$(e[1])", e[2] ? "end" : "start", length(g.seqs[e[1]]), runs[e], best, second,
                              linked ? "linked" : "-"), '\t'))
        linked || continue
        push!(get!(Set{Tuple{Int,Bool,Int}}, g.gaps, (e[1], e[2])), (f[1], !f[2], c.unknown_gap_ns))
        push!(get!(Set{Tuple{Int,Bool,Int}}, g.gaps, (f[1], f[2])), (e[1], !e[2], c.unknown_gap_ns))
        push!(done, e); push!(done, f)
    end
    # Short isolated nodes (no kept neighbour on either side, < 1 kb), typically reads lying
    # wholly inside a long microsatellite, cannot be part of one circle with the rest.
    for n in eachindex(g.seqs)
        (g.mult[n] > 0 && length(g.seqs[n]) < 1000) || continue
        any(x -> g.mult[x[1]] > 0, graph_all_nexts(g, n, true)) && continue
        any(x -> g.mult[x[1]] > 0, graph_all_nexts(g, n, false)) && continue
        g.mult[n] = 0
        push!(g.report, join(("n$n", "isolated", length(g.seqs[n]), 0, "-", 0, "dropped"), '\t'))
    end
    return g
end

# Graph rescue only: discover whole-read bridges across terminal tandem repeats.
# Recruitment deliberately bypasses the minimizer index: it uses small exact
# graph-unique flank anchors on the reads already loaded by the assembler.
const GRAPH_READ_BRIDGE_HEADER = "end_a\tend_b\tspanning_pairs\tdistinct_layouts\tvariants_bp_pairs\ttrim_a\ttrim_b\tlength_resolved\tstatus\tloaded_pair_ids"
const GRAPH_READ_BRIDGE_REPORT = Ref(String[])
struct MRAnchor
    endpoint::Int
    pos::Int
    forward::Bool
end
struct MRPlacement
    endpoint::Int
    reverse::Bool
    start::Int
end
function mr_kmers(f,s,k)
    mask=(UInt64(1)<<(2k))-1; x=UInt64(0); y=UInt64(0); run=0
    @inbounds for i in eachindex(s)
        v=ENCODE[Int(s[i])+1]
        if v==255
            x=0; y=0; run=0; continue
        end
        x=((x<<2)|UInt64(v))&mask
        y=(y>>2)|(UInt64(3-v)<<(2(k-1))); run+=1
        run>=k && f(min(x,y),x<=y,i-k+1)
    end
end
function mr_ends(g)
    out=NamedTuple{(:key,:q,:boundary,:trim,:motif),Tuple{Tuple{Int,Bool},Vector{UInt8},Int,Int,Vector{UInt8}}}[]
    for key in open_ends(g)
        q=graph_seq(g,key...); p,r=periodic_suffix(q)
        # A trim across a non-repeat trailing sequence needs a different model;
        # do not silently discard such sequence in this narrowly scoped rescue.
        (p>0 && length(q)-r>=21) || continue
        push!(out,(key=key,q=q,boundary=length(q)-r,trim=r,motif=q[end-r+1:end-r+p]))
    end
    out
end
function mr_anchors(g,ends;k=21,flank=600)
    wanted=Dict{UInt64,Vector{MRAnchor}}()
    for (i,e) in enumerate(ends)
        lo=max(1,e.boundary-flank+1)
        mr_kmers(@view(e.q[lo:e.boundary]),k) do key,forward,p
            pos=lo+p-1
            (e.boundary-pos<60 || mod(p-1,5)==0) || return
            low_complexity_kmer(key,k) && return
            push!(get!(wanted,key,MRAnchor[]),MRAnchor(i,pos,forward))
        end
    end
    counts=Dict(key=>0 for key in keys(wanted))
    for q in g.seqs
        mr_kmers(q,k) do key,_,pos
            haskey(counts,key) && (counts[key]+=1)
        end
    end
    filter!(kv->counts[first(kv)]==1,wanted)
    wanted
end
function mr_placements(s,anchors,ends,identity;k=21)
    proposals=nothing
    mr_kmers(s,k) do key,forward,p
        refs=get(anchors,key,nothing); refs===nothing && return
        proposals===nothing && (proposals=Set{Tuple{Int,Bool,Int}}())
        for a in refs
            rev=a.forward!=forward
            rp=rev ? length(s)-p-k+2 : p
            push!(proposals,(a.endpoint,rev,a.pos-rp+1))
        end
    end
    proposals===nothing && return nothing
    out=MRPlacement[]; sr=reverse_complement(s)
    for (ei,rev,st) in proposals
        e=ends[ei]; r=rev ? sr : s
        a=max(1,2-st); b=min(length(r),e.boundary-st+1)
        b-a+1>=k || continue
        sequences_similar(@view(e.q[st+a-1:st+b-1]),@view(r[a:b]),identity) || continue
        push!(out,MRPlacement(ei,rev,st))
    end
    out
end
function mr_read_bridge(s,hs,ends)
    hs===nothing && return nothing
    length(hs)==2 || return nothing
    a,b=hs
    a.reverse && ((a,b)=(b,a))
    (a.reverse || !b.reverse || a.endpoint==b.endpoint) && return nothing
    ea,eb=ends[a.endpoint],ends[b.endpoint]
    ea.key[1]==eb.key[1] && return nothing # self-circle is left to existing closure checks
    l=ea.boundary-a.start+1; r=length(s)-(eb.boundary-b.start+1)+1
    1<=l<r<=length(s) || return nothing
    middle=s[l+1:r-1]
    length(middle)>=length(ea.motif) || return nothing
    all(i->middle[i]==ea.motif[mod1(i,length(ea.motif))],eachindex(middle)) || return nothing
    # Normalize direction once so the two mates cannot cast separate votes.
    if a.endpoint<b.endpoint
        return (a.endpoint,b.endpoint,dna_string(middle),l,length(s)-r+1)
    end
    return (b.endpoint,a.endpoint,dna_string(reverse_complement(middle)),length(s)-r+1,l)
end
function mr_pair_bridge(s1,s2,anchors,ends,c)
    votes=Tuple{Int,Int,String,Int,Int}[]
    for s in (s1,s2)
        b=mr_read_bridge(s,mr_placements(s,anchors,ends,c.min_identity),ends)
        b===nothing || push!(votes,b)
    end
    isempty(votes) && return nothing
    signatures=unique((v[1],v[2],v[3]) for v in votes)
    length(signatures)==1 || return nothing # discordant mates abstain
    first(votes)
end
function mr_install!(g,ends,key,variants,layouts,ids,c)
    a,b=key; ea,eb=ends[a],ends[b]; n=sum(values(variants))
    ranked=sort(collect(variants);by=x->(-last(x),first(x)))
    best,top=first(ranked); second=length(ranked)>1 ? last(ranked[2]) : 0
    # A mode alone is not enough: require >=3 fragments, majority and 2x runner-up.
    resolved=top>=max(3,c.min_paired_support) && top>second && top>=2second && 2top>n
    status=resolved ? "read_bridge" : "read_bridge_gap"
    if !resolved && c.unknown_gap_ns<=0; status="unresolved_length_disabled"; end
    lengths=Dict{Int,Int}()
    for (s,nv) in variants; lengths[length(s)]=get(lengths,length(s),0)+nv; end
    summary=join(("$l:$nv" for (l,nv) in sort(collect(lengths);by=first)),",")
    name(e)="n$(e.key[1]):$(e.key[2] ? "end" : "start")"
    push!(g.read_bridge_report,join((name(ea),name(eb),n,length(layouts),summary,ea.trim,eb.trim,resolved,status,join(sort(ids),",")),'\t'))
    status=="unresolved_length_disabled" && return
    q=resolved ? Vector{UInt8}(codeunits(best)) : fill(UInt8('N'),c.unknown_gap_ns)
    u,ou=ea.key; v,ov=eb.key
    g.read_bridges[(u,ou,v,!ov)]=GraphReadBridge(ea.trim,eb.trim,q,n,!resolved)
    g.read_bridges[(v,ov,u,!ou)]=GraphReadBridge(eb.trim,ea.trim,reverse_complement(q),n,!resolved)
    # -2 is a typed bridge edge; walk_sequence MUST use read_bridges metadata.
    push!(get!(Set{Tuple{Int,Bool,Int}},g.gaps,(u,ou)),(v,!ov,-2))
    push!(get!(Set{Tuple{Int,Bool,Int}},g.gaps,(v,ov)),(u,!ou,-2))
    push!(g.report,join(("n$u",ou ? "end" : "start",length(ea.q),ea.trim,"n$v$(ov ? "+" : "-"):$n",0,status),'\t'))
end
"Same tandem motif up to rotation and strand."
function same_motif(a::AbstractVector{UInt8}, b::AbstractVector{UInt8})
    (length(a) == length(b) && !isempty(a)) || return false
    p = length(a); rb = reverse_complement(b)
    for s in 0:p-1
        rot = vcat(a[s+1:end], a[1:s])
        (rot == b || rot == rb) && return true
    end
    return false
end

"""
Estimated repeat bridge between two repeat ends with the same motif whose spanning reads are
too few or too discordant for read_bridge (>= 3 fragments, 2x runner-up), e.g. a microsatellite
shorter than the reads with coverage dropout. The two ends must be each other's only partner
(2x the support of any other end). Length: the most frequent spanning variant; ties go to the
variant nearest the median length, then the shorter. Without spanning reads the ends stay open
(inferred_join may still join them with Ns). Closure estimated_repeat, never validated.
"""
function mr_estimated_bridge!(g,ends,key,variants,layouts,ids,support,c)
    a,b=key; ea,eb=ends[a],ends[b]
    same_motif(ea.motif,eb.motif) || return false
    n=sum(values(variants)); n>=1 || return false
    function other(x,y)
        m=0
        for z in eachindex(ends)
            (z==x || z==y) && continue
            m=max(m,support(x,z))
        end
        return m
    end
    (n>=2*other(a,b) && n>=2*other(b,a)) || return false
    lens=Int[]
    for (s,nv) in variants; append!(lens,fill(length(s),nv)); end
    sort!(lens); med=lens[cld(length(lens),2)]
    best=first(first(sort(collect(variants);by=x->(-last(x),abs(length(first(x))-med),length(first(x)),first(x)))))
    lengths=Dict{Int,Int}()
    for (s,nv) in variants; lengths[length(s)]=get(lengths,length(s),0)+nv; end
    summary=join(("$l:$nv" for (l,nv) in sort(collect(lengths);by=first)),",")
    name(e)="n$(e.key[1]):$(e.key[2] ? "end" : "start")"
    push!(g.read_bridge_report,join((name(ea),name(eb),n,length(layouts),summary,ea.trim,eb.trim,false,"estimated_repeat",join(sort(ids),",")),'\t'))
    q=Vector{UInt8}(codeunits(best))
    u,ou=ea.key; v,ov=eb.key
    g.read_bridges[(u,ou,v,!ov)]=GraphReadBridge(ea.trim,eb.trim,q,n,false,true)
    g.read_bridges[(v,ov,u,!ou)]=GraphReadBridge(eb.trim,ea.trim,reverse_complement(q),n,false,true)
    push!(get!(Set{Tuple{Int,Bool,Int}},g.gaps,(u,ou)),(v,!ov,-2))
    push!(get!(Set{Tuple{Int,Bool,Int}},g.gaps,(v,ov)),(u,!ou,-2))
    push!(g.report,join(("n$u",ou ? "end" : "start",length(ea.q),ea.trim,"n$v$(ov ? "+" : "-"):$n",0,"estimated_repeat"),'\t'))
    return true
end
function graph_read_bridges!(g::PathGraph,idx::ReadIndex,c::Config)
    ends=mr_ends(g); length(ends)>=2 || return g
    anchors=mr_anchors(g,ends); isempty(anchors) && return g
    np=num_pairs(idx.store); chunks=min(8,Threads.nthreads(),max(1,np))
    # Store only positive whole-read bridges, not all recruited reads.
    found=[Tuple{Int,Int,String,Int,Int,UInt32}[] for _ in 1:chunks]
    Threads.@threads :dynamic for chunk in 1:chunks
        for pid in fld((chunk-1)*np,chunks)+1:fld(chunk*np,chunks)
            excluded_read(read_id(pid,true),c) && continue
            s1=get_read(idx.store,read_id(pid,true))
            s2=c.paired ? get_read(idx.store,read_id(pid,false)) : UInt8[]
            b=mr_pair_bridge(s1,s2,anchors,ends,c)
            b===nothing || push!(found[chunk],(b...,UInt32(pid)))
        end
    end
    variants=Dict{Tuple{Int,Int},Dict{String,Int}}()
    layouts=Dict{Tuple{Int,Int},Set{Tuple{Int,Int}}}()
    ids=Dict{Tuple{Int,Int},Vector{UInt32}}()
    for rows in found, (a,b,s,l,r,id) in rows
        key=(a,b); vv=get!(variants,key,Dict{String,Int}()); vv[s]=get(vv,s,0)+1
        push!(get!(layouts,key,Set{Tuple{Int,Int}}()),(l,r))
        push!(get!(ids,key,UInt32[]),id)
    end
    support(a,b)=sum(values(get(variants,minmax(a,b),Dict{String,Int}()));init=0)
    function partner(a)
        ranked=sort([(b,support(a,b)) for b in eachindex(ends) if b!=a];by=x->(-x[2],x[1]))
        b,n=first(ranked); next=length(ranked)>1 ? ranked[2][2] : 0
        key=minmax(a,b)
        n>=max(3,c.min_paired_support) && n>=2next && length(get(layouts,key,Set{Tuple{Int,Int}}()))>=2 ? b : 0
    end
    for key in sort(collect(keys(variants)))
        a,b=key
        if partner(a)==b && partner(b)==a
            mr_install!(g,ends,key,variants[key],layouts[key],ids[key],c)
        elseif mr_estimated_bridge!(g,ends,key,variants[key],layouts[key],ids[key],support,c)
            # installed as estimated_repeat
        else
            n=sum(values(variants[key])); ea,eb=ends[a],ends[b]
            push!(g.read_bridge_report,join(("n$(ea.key[1]):$(ea.key[2] ? "end" : "start")","n$(eb.key[1]):$(eb.key[2] ? "end" : "start")",n,length(layouts[key]),".",ea.trim,eb.trim,false,"insufficient_or_competing_support",join(sort(ids[key]),",")),'\t'))
        end
    end
    return g
end

"""
Decontamination: graph nodes >= 500 bp with < 0.5% of their 21-mers in the plastid panel are
foreign (e.g. a bacterium entered through a shared rRNA). Reads carrying >= 2 indexed 31-mers
of those nodes that are absent from all other nodes (so shared rRNA stretches do not count)
are flagged; extension and node coverage then ignore them. Returns (flagged reads, foreign bp).
"""
function flag_foreign_reads!(g::PathGraph, idx::ReadIndex, c::Config)
    (isempty(g.plastid) || all(isnan, g.plastid)) && return 0, 0
    foreign = [n for n in eachindex(g.seqs) if !isnan(g.plastid[n]) && length(g.seqs[n]) >= 500 && g.plastid[n] < 0.005]
    total = sum((length(g.seqs[n]) for n in foreign); init=0)
    total >= 5000 || return 0, total
    others = Set{UInt64}()
    for n in eachindex(g.seqs)
        n in foreign || union!(others, all_kmer_keys(g.seqs[n], c.kmer))
    end
    hits = zeros(UInt8, 2*num_pairs(idx.store) + 2)
    previous = FOREIGN_READS[]
    for n in foreign, key in all_kmer_keys(g.seqs[n], c.kmer)
        key in others && continue
        bucket = get(idx.postings, key, nothing)
        bucket === nothing && continue
        for id in bucket
            hits[id] < 0xff && (hits[id] += 0x01)
        end
    end
    flags = BitVector(hits .>= 0x02)
    isempty(previous) || (flags .|= previous)            # rounds accumulate
    FOREIGN_READS[] = flags
    return count(flags), total
end

"Graph of the given paths with read coverage from the index, dead ends linked by read pairs."
function assembly_graph(paths::Vector{Vector{UInt8}}, idx::ReadIndex, c::Config; read_bridges::Bool=c.graph_read_bridges)
    panel = isempty(c.plastid_panel) ? nothing : load_plastid_panel(c.plastid_panel)
    g = PathGraph(paths, s -> node_coverage(s, idx, c); panel=panel)
    read_bridges && graph_read_bridges!(g, idx, c)
    return link_dead_ends!(g, idx, c)
end

"GFA 1 of the graph (Bandage): S with LN and DP (coverage), L for k-1 overlaps and N gaps (0M, NG tag)."
function write_graph_gfa(path::String, g::PathGraph)
    open(path, "w") do io
        println(io, "H\tVN:Z:1.0")
        for n in eachindex(g.seqs)
            println(io, "S\tn$n\t", dna_string(g.seqs[n]), "\tLN:i:", length(g.seqs[n]),
                    "\tDP:f:", round(g.cov[n]; digits=2))
        end
        seen = Set{NTuple{4,Any}}()
        sgn(p) = p ? "+" : "-"
        for n in eachindex(g.seqs), plus in (true, false)
            for (v, o, gap) in graph_all_nexts(g, n, plus)
                key = min((n, plus, v, o), (v, !o, n, !plus))
                key in seen && continue
                push!(seen, key)
                bridge=get(g.read_bridges,(n,plus,v,o),nothing)
                if bridge!==nothing
                    println(io,"L\tn$n\t",sgn(plus),"\tn$v\t",sgn(o),"\t0M\tRB:Z:",dna_string(bridge.sequence),"\tTF:i:",bridge.trim_from,"\tTT:i:",bridge.trim_to,"\tRP:i:",bridge.supporting_pairs,"\tRU:i:",Int(bridge.unresolved))
                    continue
                end
                println(io, "L\tn$n\t", sgn(plus), "\tn$v\t", sgn(o), "\t",
                        gap < 0 ? "$(g.k-1)M" : "0M\tNG:i:$gap")
            end
        end
    end
end

"""
Node table. links_start / links_end: kept neighbours (multiplicity > 0) before the node's
start and after its end, overlaps and N gaps. A kept node with 0 on one side is a dead end
(the graph is broken there and no circle can pass); its terminal 40 bases are listed.
"""
function write_graph_nodes(path::String, g::PathGraph)
    open(path, "w") do io
        println(io, "node\tlength\tcoverage\trelative_coverage\tmultiplicity\tlinks_start\tlinks_end\tstart_40\tend_40\tplastid_kmers\tmin_visits\tmax_visits")
        for n in eachindex(g.seqs)
            rel = g.base > 0 ? g.cov[n] / g.base : NaN
            kept(list) = count(x -> g.mult[x[1]] > 0, list)
            ls = kept(graph_all_nexts(g, n, false)); le = kept(graph_all_nexts(g, n, true))
            s = g.seqs[n]
            println(io, join(("n$n", length(s), round(g.cov[n]; digits=2), round(rel; digits=3), g.mult[n], ls, le,
                              dna_string(s[1:min(40, end)]), dna_string(s[max(1, end-39):end]),
                              isempty(g.plastid) || isnan(g.plastid[n]) ? "-" : round(g.plastid[n]; digits=4),g.lo[n],g.hi[n]), '\t'))
        end
    end
end

# Mate-derived reseeding. A recruited seed discovers sequence; it is NOT a gap edge.
const MATE_RESCUE_HEADER = "call\tround\tsite\tkind\tcandidates\tcluster_pairs\tseed_id\tlength_bp\tnew_complex_kmers\tstatus"
const MATE_RESCUE_REPORT = String[]
const MATE_RESCUE_RESULTS = AssemblyResult[]
const MATE_RESCUE_LAST_RESULTS = Ref(AssemblyResult[])
const MATE_RESCUE_CALL = Ref(0)
function rescue_row(round,site,kind,candidates,support,id,len,novel,status)
    push!(MATE_RESCUE_REPORT,join((MATE_RESCUE_CALL[],round,site,kind,candidates,support,id,len,novel,status),'\t'))
end
function rescue_keys(s,k=21)
    out=Set{UInt64}()
    mr_kmers(s,k) do key,_,p
        low_complexity_kmer(key,k) || push!(out,key)
    end
    out
end
rescue_novel(s,known;k=21)=count(x->!(x in known),rescue_keys(s,k))

"Oriented overlap must support shared new kmers; a common repeat alone cannot cluster mates."
function rescue_overlap(a,b,c;k=21)
    positions=Dict{UInt64,Vector{Int}}()
    for (key,p) in forward_kmer_positions(a,k)
        low_complexity_kmer(key,k) || push!(get!(positions,key,Int[]),p)
    end
    offsets=Dict{Int,Int}()
    for (key,p) in forward_kmer_positions(b,k), q in get(positions,key,Int[])
        d=q-p; offsets[d]=get(offsets,d,0)+1
    end
    for (d,hits) in offsets
        hits>=2 || continue
        lo=max(1,1+d); hi=min(length(a),length(b)+d)
        hi-lo+1>=max(k+1,min(c.min_overlap,40)) || continue
        sequences_similar(@view(a[lo:hi]),@view(b[lo-d:hi-d]),c.min_identity) && return true
    end
    false
end
# A compact position iterator, separate from minimizer sampling.
function forward_kmer_positions(s,k)
    out=Tuple{UInt64,Int}[]; mask=(UInt64(1)<<(2k))-1; key=UInt64(0); run=0
    for i in eachindex(s)
        v=ENCODE[Int(s[i])+1]
        if v==255; key=0;run=0;continue;end
        key=((key<<2)|UInt64(v))&mask;run+=1
        run>=k && push!(out,(key,i-k+1))
    end
    out
end
function mate_seed_groups(mates,known,c;k=21)
    ids=sort!(collect(keys(mates))); n=length(ids); parent=collect(1:n)
    root(i)=begin
        while parent[i]!=i;parent[i]=parent[parent[i]];i=parent[i];end
        i
    end
    novel=[setdiff(rescue_keys(mates[id],k),known) for id in ids]
    postings=Dict{UInt64,Vector{Int}}(); degree=zeros(Int,n)
    for i in 1:n
        shared=Dict{Int,Int}()
        for key in novel[i], j in get(postings,key,Int[])
            shared[j]=get(shared,j,0)+1
        end
        for j in sort!(collect(keys(shared)))
            shared[j]>=2 || continue
            root(i)==root(j) && continue
            rescue_overlap(mates[ids[i]],mates[ids[j]],c;k=k) || continue
            parent[root(i)]=root(j);degree[i]+=1;degree[j]+=1
        end
        for key in novel[i];push!(get!(postings,key,Int[]),i);end
    end
    components=Dict{Int,Vector{Int}}()
    for i in 1:n;push!(get!(components,root(i),Int[]),i);end
    groups=NamedTuple[]
    for members in values(components)
        length(members)>=max(3,c.min_coverage) || continue
        counts=Dict{UInt64,Int}()
        for i in members,key in novel[i];counts[key]=get(counts,key,0)+1;end
        supported(i)=count(key->counts[key]>=2,novel[i])
        best=first(sort(members;by=i->(-supported(i),-degree[i],ids[i])))
        supported(best)>=12 || continue
        push!(groups,(best=ids[best],ids=ids[members],novel=supported(best)))
    end
    sort!(groups;by=g->(-length(g.ids),-g.novel,g.best))
end

"Include occupied branch/IR exits, not only graph dead ends. IR anchoring identifies a node, not its copy."
function rescue_sites(g,c)
    opens=Set(open_ends(g)); sites=NamedTuple[]
    for n in eachindex(g.seqs), side in (true,false)
        # link_dead_ends! excludes isolated nodes <1 kb from circular walks,
        # but leaves their positive lower bound: they can still recruit mates.
        (g.mult[n]>0 || (g.lo[n]>0 && length(g.seqs[n])<1000)) || continue
        nexts=graph_all_nexts(g,n,side)
        kind=(n,side) in opens ? "dead_end" :
             length(nexts)>1 ? "branch" :
             g.mult[n]>1 && length(g.seqs[n])>=c.ir_min_length ? "repeat_exit" : ""
        isempty(kind) && continue
        q=Vector{UInt8}(graph_seq(g,n,side)); trim=tail_run_trim(q)
        trim<length(q)-60 || (trim=0)
        boundary=length(q)-trim;boundary>=21 || continue
        push!(sites,(key=(n,side),q=q,boundary=boundary,trim=trim,motif=UInt8[],kind=kind))
    end
    sites
end
function recruit_rescue_mates(g,idx,c,known;monitor=nothing)
    sites=rescue_sites(g,c); pools=[Dict{UInt32,Vector{UInt8}}() for _ in sites]
    saturated=falses(length(sites))
    (c.paired && !isempty(sites)) || return sites,pools,saturated
    _,high=insert_bounds(c)
    anchors=mr_anchors(g,sites;k=21,flank=high+idx.max_read_length)
    isempty(anchors) && return sites,pools,saturated
    np=num_pairs(idx.store); chunks=min(8,Threads.nthreads()); cap=c.mate_rescue_max_candidates
    localp=[[Dict{UInt32,Vector{UInt8}}() for _ in sites] for _ in 1:chunks]
    localfull=[falses(length(sites)) for _ in 1:chunks]
    Threads.@threads for chunk in 1:chunks
        readbuf=UInt8[]
        for pid in fld((chunk-1)*np,chunks)+1:fld(chunk*np,chunks)
            excluded_read(read_id(pid,true),c) && continue
            votes=Tuple{Int,Vector{UInt8}}[]
            for first in (true,false)
                id=read_id(pid,first);read=get_read!(readbuf,idx.store,id)
                hits=mr_placements(read,anchors,sites,c.min_identity)
                hits===nothing && continue
                filter!(h->!h.reverse && 0<length(sites[h.endpoint].q)-h.start+1<=high,hits)
                length(hits)==1 || continue
                hit=only(hits)
                raw=get_read(idx.store,id⊻UInt32(1)); cut=genomic_length(raw)
                cut>=MIN_INITIAL_LENGTH || continue
                mate=reverse_complement(raw[1:cut])
                rescue_novel(mate,known)>=12 || continue
                push!(votes,(hit.endpoint,mate))
            end
            # The two mates cannot propose unrelated destinations with one pair vote.
            length(votes)==1 || continue
            ei,mate=only(votes); pool=localp[chunk][ei]
            if length(pool)<cap;pool[UInt32(pid)]=mate;else;localfull[chunk][ei]=true;end
        end
    end
    for ei in eachindex(sites)
        combined=Dict{UInt32,Vector{UInt8}}()
        for chunk in 1:chunks
            merge!(combined,localp[chunk][ei]);saturated[ei]|=localfull[chunk][ei]
        end
        ids=sort!(collect(keys(combined)));saturated[ei]|=length(ids)>cap
        for pid in first(ids,cap);pools[ei][pid]=combined[pid];end
    end
    sites,pools,saturated
end
function rescue_seed_file(c)
    !isempty(c.rescue_seed_fasta) && return c.rescue_seed_fasta
    isempty(c.seed_input) && return ""
    path=joinpath(dirname(c.seed_input),c.project_name*"_seeds_all.fasta")
    isfile(path) ? path : ""
end
const RefTile = NamedTuple{(:stretch,:header,:seq),Tuple{Int,String,Vector{UInt8}}}
const REF_RESCUE_TILES = Ref(RefTile[])
const REF_RESCUE_NOTE = Ref("")
const REF_RESCUE_ACTIVE = Ref(false)      # only the first graph rescue of a run uses reference seeds
const REF_RESCUE_PATHS = Vector{UInt8}[]   # paths grown from reference seeds, kept for the reference scaffold

"""
Reference rescue seeds, set before graph rescue when there is no gap-free circle. The nearest
reference of the paths (k-mer containment >= 0.2) is compared with the paths by 21-mers; a
reference stretch longer than `min_gap` bp without any shared 21-mer is missing from the paths
(divergence alone leaves shared 21-mers every few dozen bases at ANI >= 0.95). Tiles of `tile` bp
every `step` bp across each stretch become reference seeds (reference_round_seeds). Like SeedRef's
seeds they only select a library read, which is then extended from the reads: no reference base
enters the result.
"""
function reference_rescue_seeds!(paths::Vector{Vector{UInt8}}, c::Config; min_gap::Int=300, tile::Int=300,
                                 step::Int=1500, max_tiles::Int=200)
    REF_RESCUE_TILES[] = RefTile[]; REF_RESCUE_NOTE[] = ""
    (isempty(c.reference_dataset) || !c.reference_rescue) && return
    sk = joinpath(c.reference_dataset, RefSketch.SKETCH_FILE)
    isfile(sk) || (REF_RESCUE_NOTE[] = "Reference rescue: missing $(sk)"; return)
    query = [q for q in paths if length(q) >= 5000]
    isempty(query) && (REF_RESCUE_NOTE[] = "Reference rescue: no path >= 5 kb"; return)
    hits, _ = RefSketch.nearest(RefSketch.load(sk), query, fill(false, length(query)); top=1, exclude=not_self(c))
    isempty(hits) && return
    best = hits[1]
    if best.containment < 0.2
        REF_RESCUE_NOTE[] = @sprintf("Reference rescue: nearest %s too distant (containment %.3f), not used",
                                     best.accession, best.containment)
        return
    end
    ref = RefSketch.reference_sequence(c.reference_dataset, best.accession); n = length(ref)
    n >= 2tile || return
    have = Set{UInt64}()
    for q in paths; union!(have, rescue_keys(q)); end
    hitpos = Int[]
    mr_kmers(ref, 21) do key, _, pos
        key in have && push!(hitpos, pos)
    end
    sort!(hitpos)
    if isempty(hitpos)
        REF_RESCUE_NOTE[] = "Reference rescue: $(best.accession) shares no 21-mer with the paths"; return
    end
    circ(a, b) = a <= b ? ref[a:b] : vcat(ref[a:end], ref[1:b])      # 1-based, circular
    gaps = Tuple{Int,Int}[]
    for i in eachindex(hitpos)
        a = hitpos[i]; b = i < length(hitpos) ? hitpos[i+1] : hitpos[1] + n
        b - a > min_gap + 21 && push!(gaps, (a + 21, b - 1))           # positions not covered by any shared 21-mer
    end
    tiles = RefTile[]
    for (gi, (a, b)) in enumerate(gaps)
        len = b - a + 1; starts = collect(a:step:max(a, b - tile + 1))
        last(starts) < b - tile + 1 && push!(starts, b - tile + 1)
        for st in starts
            length(tiles) < max_tiles || break
            s0 = mod1(st, n); s1 = mod1(st + min(tile, len) - 1, n)
            push!(tiles, (stretch=gi, header="reference_gap $(best.accession):$(s0)-$(s1)", seq=circ(s0, s1)))
        end
    end
    REF_RESCUE_TILES[] = tiles
    missing_bp = sum((b - a + 1 for (a, b) in gaps); init=0)
    REF_RESCUE_NOTE[] = @sprintf("Reference rescue: %s (containment %.3f); %d stretch(es) of %d bp missing from the paths; %d seed tile(s)",
                                 best.accession, best.containment, length(gaps), missing_bp, length(tiles))
    return
end

"""
Reference seeds of one rescue round, outside the mate seed budget (a rescue for libraries that do
not assemble): per missing stretch at most one tile, the first not tried yet whose sequence is still
absent from the paths (no 21-mer in `known`; stretches filled by earlier rounds are skipped), turned
into the library read with most shared minimizers; at most `c.reference_rescue_seeds` per round.
"""
function reference_round_seeds(idx, c, known, tried::Set{String})
    out = NamedTuple[]; used = Set{Int}()
    REF_RESCUE_ACTIVE[] || return out
    for t in REF_RESCUE_TILES[]
        length(out) >= c.reference_rescue_seeds && break
        (t.stretch in used || t.header in tried) && continue
        push!(tried, t.header)
        any(in(known), rescue_keys(t.seq)) && continue          # already in the paths
        scores = Dict{UInt32,Int}()
        for key in all_kmer_keys(t.seq, c.kmer), id in get(idx.postings, key, UInt32[])
            excluded_read(UInt32(id), c) && continue
            scores[id] = get(scores, id, 0) + 1
        end
        for (id, score) in sort!(collect(scores); by=x -> (-last(x), first(x)))
            score >= c.min_seed_hits || break
            raw = get_read(idx.store, id); cut = genomic_length(raw); cut >= MIN_INITIAL_LENGTH || continue
            q = raw[1:cut]
            rescue_novel(q, known) >= 12 || continue
            push!(used, t.stretch)
            push!(out, (q=q, site=t.header, kind="reference_seed", candidates=length(scores), support=0))
            break
        end
    end
    return out
end

function extra_rescue_seeds(idx,c,known,attempted)
    records=Pair{String,Vector{UInt8}}[]
    path=rescue_seed_file(c)
    if !isempty(path)
        isfile(path) || error("Rescue seed FASTA missing: $path")
        append!(records,read_seeds(path))
    end
    isempty(records) && return NamedTuple[]
    out=NamedTuple[]; seen=copy(attempted)
    for (header,seq) in records
        scores=Dict{UInt32,Int}()
        for key in all_kmer_keys(seq,c.kmer), id in get(idx.postings,key,UInt32[])
            excluded_read(UInt32(id),c) && continue
            scores[id]=get(scores,id,0)+1
        end
        # One actual library read per record gives every marker a fair first chance.
        for (id,score) in sort!(collect(scores);by=x->(-last(x),first(x)))
            score>=c.min_seed_hits || break
            raw=get_read(idx.store,id);cut=genomic_length(raw);cut>=MIN_INITIAL_LENGTH || continue
            q=raw[1:cut];key=min(dna_string(q),dna_string(reverse_complement(q)))
            key in seen && continue
            rescue_novel(q,known)>=12 || continue
            push!(seen,key)
            push!(out,(q=q,site=replace(header,r"[\t\r\n]"=>" "),kind="extra_seed",candidates=length(scores),support=0))
            break
        end
    end
    out
end
# Topology fallback changes traversal constraints, never invents a sequence edge.
const GRAPH_TOPOLOGY_REPORT=String[]
const GRAPH_TOPOLOGY_HEADER="attempt\ttype\tstart\tend\tdetail"
const GRAPH_TOPOLOGY_ATTEMPT=Ref(0)
const OrientedNode=Tuple{Int,Bool}
function topology_row(kind,a,b,detail)
    push!(GRAPH_TOPOLOGY_REPORT,join((GRAPH_TOPOLOGY_ATTEMPT[],kind,a,b,detail),'\t'))
end
function topology_arm_sequence(g,path)
    # Only a boundary k-mer from each flank; large flanks must not hide divergence.
    q=copy(graph_seq(g,path[1]...)[end-g.k+1:end])
    for v in path[2:end-1];append!(q,@view graph_seq(g,v...)[g.k:end]);end
    append!(q,@view graph_seq(g,path[end]...)[g.k:g.k])
    q
end
function topology_bubbles(g;max_bases=2000,max_steps=512,max_nodes=24)
    groups=NamedTuple[];seen=Set{String}()
    name(p)=join(("$n$(o ? "+" : "-")" for (n,o) in p),",")
    for n in eachindex(g.seqs), o in (false,true)
        g.mult[n]>0 || continue
        start=(n,o)
        nx=filter(x->g.mult[x[1]]>0 && x[1]!=n,graph_nexts(g,n,o))
        length(nx)>=2 || continue
        routes=Dict{OrientedNode,Vector{Vector{OrientedNode}}}()
        steps=Ref(0);saturated=Ref(false)
        function visit(path,bases)
            steps[]+=1
            steps[]>max_steps && (saturated[]=true;return)
            v=last(path)
            push!(get!(routes,v,Vector{OrientedNode}[]),copy(path))
            (length(path)>=max_nodes || bases+length(g.seqs[v[1]])-g.k+1>max_bases) && return
            # Large sequences are boundary nodes, never made optional by this rule.
            length(g.seqs[v[1]])>max_bases && return
            for w in graph_nexts(g,v...)
                g.mult[w[1]]>0 || continue
                any(x->x[1]==w[1],path) && continue
                push!(path,w);visit(path,bases+length(g.seqs[v[1]])-g.k+1);pop!(path)
                saturated[] && return
            end
        end
        for v in nx
            visit(OrientedNode[start,v],0);saturated[] && break
        end
        saturated[] && (topology_row("bubble_search_limit","n$n","-","steps=$max_steps; no relaxation from this start");continue)
        for target in sort!(collect(keys(routes)))
            paths=routes[target]
            for i in eachindex(paths),j in 1:i-1
                a,b=paths[i],paths[j]
                (length(a)>=3 && length(b)>=3 && a[2]!=b[2]) || continue
                ia=Set(first.(a[2:end-1]));ib=Set(first.(b[2:end-1]))
                isdisjoint(ia,ib) || continue
                # Both physical flanks must provide complex sequence.
                min(length(rescue_keys(graph_seq(g,start...)[end-g.k+1:end])),
                    length(rescue_keys(graph_seq(g,target...)[1:g.k])))>=8 || continue
                qa=topology_arm_sequence(g,a);qb=topology_arm_sequence(g,b)
                # A local structural alternative can include indels and divergent
                # haplotypes. This is a bounded hypothesis, never a validated allele.
                allowed=min(150,max(3,floor(Int,0.15max(length(qa),length(qb)))))
                abs(length(qa)-length(qb))<=allowed || continue
                edits=banded_edits(qa,qb;band=allowed)
                edits<=allowed || continue
                ra=[(v,!p) for (v,p) in reverse(a)];rb=[(v,!p) for (v,p) in reverse(b)]
                key=min(join(sort([name(a),name(b)]),"|"),join(sort([name(ra),name(rb)]),"|"))
                key in seen && continue;push!(seen,key)
                nodes=union(ia,ib)
                push!(groups,(a=copy(a),b=copy(b),nodes=nodes,edits=edits))
                topology_row("alternative_paths",name([start]),name([target]),
                    "arms=$(name(a))|$(name(b)); bp=$(length(qa)),$(length(qb)); edits=$edits")
            end
        end
    end
    groups
end
function topology_long_exits(g,start,c;max_bases=3000,max_nodes=32)
    exits=Set{OrientedNode}();best=Dict{OrientedNode,Int}();queue=[(start,0,0)];head=1
    while head<=length(queue)
        v,dist,depth=queue[head];head+=1
        depth>=max_nodes && continue
        for w in graph_nexts(g,v...)
            (g.mult[w[1]]>0 && w[1]!=start[1]) || continue
            d=dist+(v==start ? 0 : length(g.seqs[v[1]])-g.k+1)
            d<=max_bases || continue
            if length(g.seqs[w[1]])>=c.ir_min_length && g.lo[w[1]]>0
                push!(exits,w);continue
            end
            d<get(best,w,typemax(Int)) || continue
            best[w]=d;push!(queue,(w,d,depth+1))
        end
    end
    exits
end
function topology_repeat_candidates(g,c)
    out=Int[]
    for n in eachindex(g.seqs)
        (g.mult[n]>0 && g.hi[n]==1 && length(g.seqs[n])>=c.ir_min_length) || continue
        a=topology_long_exits(g,(n,true),c);b=topology_long_exits(g,(n,false),c)
        (length(a)>=2 && length(b)>=2 && isdisjoint(Set(first.(a)),Set(first.(b)))) || continue
        push!(out,n)
        topology_row("repeat_candidate","n$n","-","length=$(length(g.seqs[n])); exits=$(join(sort!(collect(a)),','))|$(join(sort!(collect(b)),',')); coverage=$(g.cov[n]); base=$(g.base)")
    end
    out
end
function topology_path_in_walk(path,walknodes)
    N=length(walknodes);length(path)>N+1 && return false
    for p in (path,[(n,!o) for (n,o) in reverse(path)])
        for i in 1:N
            all(j->walknodes[mod1(i+j-1,N)]==p[j],eachindex(p)) && return true
        end
    end
    false
end
function topology_constraints_ok(original,adjusted,groups,lifted,w)
    visits=w[1:end-1];counts=zeros(Int,length(original.seqs));nodes=[(n,o) for (n,o,_) in visits]
    for (n,_,_) in visits;counts[n]+=1;end
    # Every waived obligatory node needs an actually traversed alternative arm.
    covered=Set{Int}()
    for group in groups
        (topology_path_in_walk(group.a,nodes)||topology_path_in_walk(group.b,nodes)) && union!(covered,group.nodes)
    end
    for n in eachindex(counts)
        counts[n]<original.lo[n] && !(n in covered) && return false
    end
    for n in lifted
        counts[n]<=1 && continue
        counts[n]==2 || return false
        count(x->x[1]==n && x[2],visits)==1 || return false
    end
    true
end
function topology_tandem_limits!(h,source_paths)
    loops=Set(n for n in eachindex(h.seqs) if length(h.paths[n])==1 && isempty(rescue_keys(h.seqs[n])))
    isempty(loops) && return
    maxima=Dict{Int,Int}()
    for q in source_paths
        previous=0;run=0
        for i in 1:length(q)-h.k+1
            key=code128(q,i,h.k)
            node=key==typemax(UInt128) ? 0 : first(get(h.kinfo,canon128(key,h.k),(0,0)))
            if node in loops
                run=node==previous ? run+1 : 1
                maxima[node]=max(get(maxima,node,0),run)
            else
                run=0
            end
            previous=node
        end
    end
    for (n,observed) in maxima
        old=h.hi[n];h.hi[n]=max(h.lo[n],min(old,observed))
        h.hi[n]<old && topology_row("observed_tandem_limit","n$n","-","old_max=$old; observed_max=$observed; new_max=$(h.hi[n])")
    end
end
function topology_graph_rescue(g,c;source_paths=Vector{UInt8}[])
    (c.graph_topology_rescue && c.assembly_type=="chl") || return AssemblyResult[],g
    GRAPH_TOPOLOGY_ATTEMPT[]+=1
    groups=topology_bubbles(g);lifted=topology_repeat_candidates(g,c)
    isempty(groups) && isempty(lifted) && return AssemblyResult[],g
    # Bound the hypothesis space; never lift an arbitrary number of long sequences.
    length(lifted)>8 && (topology_row("repeat_limit","-","-","$(length(lifted)) candidates; no topology rescue");return AssemblyResult[],g)
    h=deepcopy(g)
    for group in groups,n in group.nodes
        # An alternative arm does not waive a mandatory multicopy requirement.
        g.lo[n]==1 && (h.lo[n]=0)
    end
    for n in lifted;h.hi[n]=2;end
    topology_tandem_limits!(h,source_paths)
    stats=Dict{String,Int}();cache=Dict{String,Vector{UInt8}}()
    failures=Dict("constraints"=>0,"length"=>0,"structure"=>0,"IR"=>0)
    function accept(w)
        topology_constraints_ok(g,h,groups,lifted,w) || (failures["constraints"]+=1;return false)
        q=walk_sequence(h,w)
        c.genome_range[1]<=length(q)<=c.genome_range[2] || (failures["length"]+=1;return false)
        plausible_plastome(q,c) || (failures["structure"]+=1;return false)
        twice=[n for n in lifted if count(x->x[1]==n,w[1:end-1])==2]
        if !isempty(twice)
            item,_=PlastomeOrder.order_plastome(q;min_ir=c.ir_min_length,min_sc=c.ir_min_sc_length,min_identity=c.ir_min_identity)
            item===nothing && (failures["IR"]+=1;return false)
            ir=item.sequence[item.lsc_length+1:item.lsc_length+item.ir_pair.len]
            text=dna_string(ir);rev=dna_string(reverse_complement(ir))
            all(n->occursin(dna_string(h.seqs[n]),text)||occursin(dna_string(h.seqs[n]),rev),twice) ||
                (failures["IR"]+=1;return false)
        end
        cache[string(w)]=q
        true
    end
    walks=graph_walks(h;accept_walk=accept,stats=stats,prefer_coverage=true)
    topology_row("search","-","-",join(["$k=$(stats[k])" for k in sort!(collect(keys(stats)))],"; ")*"; "*
        join(["rejected_$k=$(failures[k])" for k in sort!(collect(keys(failures)))],"; "))
    out=AssemblyResult[];seen=Set{String}()
    for w in walks
        q=cache[string(w)];key=dna_string(q);key in seen && continue;push!(seen,key)
        reason=walk_inferred(h,w) ? "closed_inferred_join" : "closed_topology_inferred"
        push!(out,AssemblyResult("graph_topology_$(length(out)+1)",q,true,reason,"assembly_graph_topology",0,0,0))
        topology_row("candidate","-","-","id=$(last(out).id); length=$(length(q)); closure=$(closure_status(last(out))); gaps=$(length(n_gaps(q))); walk=$(join(["n"*string(n)*(o ? "+" : "-") for (n,o,_) in w],','))")
    end
    isempty(out) ? (out,g) : (out,h)
end

"""
accept_walk for plastome circles: length in Genome Range, then a plausible LSC-IR-SSC-IR
structure (cached per sequence). Used inside the walk search, so that its solution limit counts
only walks that pass these checks.
"""
function plastome_walk_filter(g::PathGraph, c::Config)
    cache = Dict{UInt64,Bool}()
    return function (w)
        q = walk_sequence(g, w)
        c.genome_range[1] <= length(q) <= c.genome_range[2] || return false
        return get!(() -> plausible_plastome(q, c), cache, hash(q))
    end
end

function rescue_graph_circles(g,c)
    out=AssemblyResult[]; walks=graph_walks(g;accept_walk=plastome_walk_filter(g,c))
    for (i,w) in enumerate(walks)
        q=walk_sequence(g,w)
        c.genome_range[1]<=length(q)<=c.genome_range[2] || continue
        plausible_plastome(q,c) || continue
        reason=walk_inferred(g,w) ? "closed_inferred_join" : walk_estimated(g,w) ? "closed_estimated_repeat" : "closed_graph_path"
        push!(out,AssemblyResult("graph_$i",q,true,reason,"assembly_graph",0,0,0))
    end
    out
end
function iterative_mate_graph_rescue(paths,idx,c;monitor=nothing)
    isempty(paths) && return AssemblyResult[]
    MATE_RESCUE_CALL[]+=1
    g=assembly_graph(paths,idx,c);baseline=g
    current=copy(paths);partial=AssemblyResult[];direct=AssemblyResult[]
    known=Set{UInt64}();for q in current;union!(known,rescue_keys(q));end
    attempted=Set{String}();out=rescue_graph_circles(g,c);rounds=0;ref_tried=Set{String}()
    rc=deepcopy(c);rc.max_branch_paths=c.mate_rescue_branch_paths
    # A graph walk containing gaps is not a reason to stop looking for missing sequence.
    if !has_gap_free_circle(out) && c.island_recovery
        for round in 1:c.mate_rescue_rounds
            rounds=round
            progress_phase!(monitor,"Mate recruitment";detail="round $round / $(c.mate_rescue_rounds)")
            sites,pools,saturated=recruit_rescue_mates(g,idx,c,known;monitor=monitor)
            queues=Vector{NamedTuple}[]
            for ei in eachindex(sites)
                e=sites[ei];site="n$(e.key[1]):$(e.key[2] ? "end" : "start")"
                groups=mate_seed_groups(pools[ei],known,c)
                rescue_row(round,site,e.kind,length(pools[ei]),0,"-",0,0,
                    saturated[ei] ? "candidate_limit" : isempty(groups) ? "no_coherent_cluster" : "clusters_found")
                push!(queues,[(q=pools[ei][gr.best],site=site,kind=e.kind,candidates=length(pools[ei]),support=length(gr.ids)) for gr in groups])
            end
            # Round-robin across origins; a high-coverage site cannot consume every seed slot.
            push!(queues,extra_rescue_seeds(idx,c,known,attempted))
            tasks=NamedTuple[];depth=1
            while any(q->length(q)>=depth,queues)
                for q in queues;length(q)>=depth && push!(tasks,q[depth]);end
                depth+=1
            end
            # Reference seeds first and on top of the mate budget.
            reftasks=reference_round_seeds(idx,c,known,ref_tried)
            tasks=vcat(reftasks,tasks); cap=c.mate_rescue_seeds_per_round+length(reftasks)
            progress_phase!(monitor,"Mate reseeding";total=min(length(tasks),cap),unit="seeds")
            added=0;used=0
            for task in tasks
                used<cap || break
                q=task.q;key=min(dna_string(q),dna_string(reverse_complement(q)))
                key in attempted && continue
                rescue_novel(q,known)>=12 || continue
                push!(attempted,key);used+=1
                id="mate_rescue_$(MATE_RESCUE_CALL[])_r$(round)_s$(used)"
                seed=InitialContig(id,correct_seed_read(q,idx,rc),task.site)
                assembled=assemble_seed(seed,idx,rc;monitor=monitor,length_limit=c.mate_rescue_max_bases)
                kept=0
                for result in assembled
                    novel=rescue_novel(result.sequence,known)
                    accept=length(result.sequence)>=length(q)+100 && novel>=30
                    rescue_row(round,task.site,task.kind,task.candidates,task.support,result.id,length(result.sequence),novel,
                        accept ? "retained" : "no_new_extension")
                    accept || continue
                    push!(partial,result);push!(current,result.sequence);union!(known,rescue_keys(result.sequence))
                    # Kept even if a decontamination pass replaces the results later.
                    task.kind=="reference_seed" && push!(REF_RESCUE_PATHS,result.sequence)
                    result.circular && push!(direct,result)
                    added+=1;kept+=1
                end
                isempty(assembled) && rescue_row(round,task.site,task.kind,task.candidates,task.support,id,0,0,"no_path")
            end
            if added==0
                # Unproductive seeds must not starve later sites/FASTA records merely
                # because this round exhausted its bounded number of attempts.
                deferred=any(tasks) do task
                    key=min(dna_string(task.q),dna_string(reverse_complement(task.q)))
                    !(key in attempted) && rescue_novel(task.q,known)>=12
                end
                # Reference tiles not tried yet also keep the rounds going.
                deferred=deferred || (REF_RESCUE_ACTIVE[] && any(t -> !(t.header in ref_tried), REF_RESCUE_TILES[]))
                deferred || break
                rescue_row(round,"-","budget",length(tasks),0,"-",0,0,"seeds_deferred")
                continue
            end
            # Keep the accumulated paths even when this rebuilt graph has no circle.
            g=assembly_graph(current,idx,c);out=rescue_graph_circles(g,c)
            has_gap_free_circle(vcat(out,direct)) && break
        end
    end
    if isempty(out) && c.inferred_joins
        infer_joins!(g,c);out=rescue_graph_circles(g,c)
    end
    if isempty(out)
        progress_phase!(monitor,"Graph topology";detail="checking alternatives and IR copies")
        out,g=topology_graph_rescue(g,c;source_paths=current)
    end
    # A newly discovered partial island must not erase an older circle.
    if isempty(out) && g!==baseline
        fallback=rescue_graph_circles(baseline,c)
        if isempty(fallback) && c.inferred_joins;infer_joins!(baseline,c);fallback=rescue_graph_circles(baseline,c);end
        if isempty(fallback);fallback,baseline=topology_graph_rescue(baseline,c;source_paths=paths);end
        if !isempty(fallback);g=baseline;out=fallback;end
    end
    append!(MATE_RESCUE_RESULTS,partial);MATE_RESCUE_LAST_RESULTS[]=partial
    GRAPH_RESCUE_GRAPH[]=g;GRAPH_RESCUE_LINKS[]=copy(g.report)
    GRAPH_READ_BRIDGE_REPORT[]=copy(g.read_bridge_report)
    GRAPH_RESCUE_REPORT[]="mate rescue rounds $rounds, retained paths $(length(partial)), retained bp $(sum(length(x.sequence) for x in partial;init=0)), " *
        "nodes in use $(count(>(0),g.mult)), read bridges $(length(g.read_bridges)÷2), inferred joins $(length(g.inferred)÷2); see _mate_rescue.tsv and _graph_topology.tsv"
    vcat(direct,out)
end

"Graph fallback with optional iterative reseeding; circle range and structure checks still apply."
function graph_rescue(paths::Vector{Vector{UInt8}}, idx::ReadIndex, c::Config;monitor=nothing)
    MATE_RESCUE_LAST_RESULTS[]=AssemblyResult[]
    c.mate_rescue && return iterative_mate_graph_rescue(paths,idx,c;monitor=monitor)
    isempty(paths) && return AssemblyResult[]
    g = assembly_graph(paths, idx, c)
    out = AssemblyResult[]
    nlinks = count(r -> endswith(r, "linked"), g.report)
    walks = graph_walks(g; accept_walk=plastome_walk_filter(g, c)); w0 = length(walks)
    island_lengths = Int[]; island_use = ""
    # Islands are an extra attempt, never a replacement: if the graph with them gives no circuit
    # (even with inferred joins), the graph without them is tried as before, so a recovered
    # island cannot take away a circle that existed without it.
    if isempty(walks) && c.island_recovery
        islands = recover_islands(g, idx, c)
        if !isempty(islands)
            island_lengths = length.(islands)
            gi = assembly_graph(vcat(paths, islands), idx, c)
            wi = graph_walks(gi; accept_walk=plastome_walk_filter(gi, c))
            if isempty(wi) && c.inferred_joins
                infer_joins!(gi, c); wi = graph_walks(gi; accept_walk=plastome_walk_filter(gi, c))
            end
            if !isempty(wi)
                g, walks = gi, wi; island_use = " (used)"
                nlinks = count(r -> endswith(r, "linked"), g.report)
            else
                island_use = " (discarded: no circuit with them)"
            end
        end
    end
    if isempty(walks) && c.inferred_joins
        infer_joins!(g, c)
        walks = graph_walks(g; accept_walk=plastome_walk_filter(g, c))
    end
    njoins = count(r -> endswith(r, "inferred_join"), g.report)
    out_range = 0; implausible = 0
    for (i, w) in enumerate(walks)
        q = walk_sequence(g, w)
        c.genome_range[1] <= length(q) <= c.genome_range[2] || (out_range += 1; continue)
        plausible_plastome(q, c) || (implausible += 1; continue)
        reason = walk_inferred(g, w) ? "closed_inferred_join" : walk_estimated(g, w) ? "closed_estimated_repeat" : "closed_graph_path"
        push!(out, AssemblyResult("graph_$(i)", q, true, reason, "assembly_graph", 0, 0, 0))
    end
    isempty(out) && ((out,g)=topology_graph_rescue(g,c;source_paths=paths))
    kept = count(>(0), g.mult)
    deadends = count(r -> endswith(r, "\t-") || endswith(r, "\tlinked"), g.report)
    GRAPH_RESCUE_REPORT[] = "nodes in use $(kept), dead ends $(deadends), islands $(isempty(island_lengths) ? 0 : join(island_lengths, "+")) bp$(island_use), " *
        "pair links $(nlinks), read bridges $(length(g.read_bridges)÷2), inferred joins $(njoins), " *
        "walks without joins $(w0), walks total $(length(walks)), rejected: length $(out_range), structure $(implausible)"
    GRAPH_RESCUE_GRAPH[] = g
    GRAPH_RESCUE_LINKS[] = copy(g.report)
    GRAPH_READ_BRIDGE_REPORT[] = copy(g.read_bridge_report)
    return out
end

const GRAPH_RESCUE_GRAPH = Ref{Union{Nothing,PathGraph}}(nothing)
const GRAPH_RESCUE_REPORT = Ref("")
const GRAPH_RESCUE_LINKS = Ref(String[])

"""
Circle strong enough to stop further attempts: no N and closure validated, ir_collapse,
read_bridge or graph_path. estimated_repeat (a repeat length from as little as one pair),
inferred_join, topology_inferred, ref_scaffold and unresolved circles stay candidates for the
output ranking but do not stop the search for a better-supported circle.
"""
has_gap_free_circle(results) = any(x -> x.circular && closure_status(x) in CONCLUSIVE_CLOSURES && isempty(n_gaps(x.sequence)), results)
const CONCLUSIVE_CLOSURES = ("validated", "ir_collapse", "read_bridge", "graph_path")

"""
Second pass with more reads: only if the first used a read limit and the file holds more,
no conclusive circle was closed (has_gap_free_circle, the same rule as the other stages), and
some linear path end stopped for lack of reads (low support or no extension, not branching) or
a circle was closed only weakly (N gap, estimated repeat length).
"""
function adaptive_needed(c::Config, results::Vector{AssemblyResult}; raw::Vector{AssemblyResult}=results)
    c.adaptive_read_pairs > c.read_pairs_to_process > 0 || return false
    INPUT_MORE[] || return false
    has_gap_free_circle(results) && return false
    return any(x -> !x.circular && (occursin("low_support", x.reason) || occursin("no_extension", x.reason)), raw) ||
           any(x -> x.circular, results)
end

"Run one assembly for a validated Config; returns the main FASTA path and timings."
function run_assembly(c::Config; label::AbstractString=c.project_name)
    GRAPH_RESCUE_GRAPH[]=nothing; GRAPH_READ_BRIDGE_REPORT[]=String[]; empty!(ISLAND_REPORT); empty!(ISLAND_SEQS)
    REF_RESCUE_TILES[]=RefTile[]; REF_RESCUE_NOTE[]=""; empty!(REF_SCAFFOLD_REPORT); empty!(REF_RESCUE_PATHS)
    empty!(MATE_RESCUE_REPORT);empty!(MATE_RESCUE_RESULTS);MATE_RESCUE_CALL[]=0
    empty!(GRAPH_TOPOLOGY_REPORT);GRAPH_TOPOLOGY_ATTEMPT[]=0
    MATE_RESCUE_LAST_RESULTS[]=AssemblyResult[]
    FOREIGN_READS[]=BitVector()
    mkpath(c.output_path)
    @printf("AssemblerComplete %s | Julia %s | %d threads\n", VERSION_STRING, string(VERSION), Threads.nthreads())
    @printf("k=%d window=%d overlap=%d coverage=%d identity=%.3f\n",
            c.kmer, c.minimizer_window, c.min_overlap, c.min_coverage, c.min_identity)
    t0 = time_ns()
    monitor = c.show_progress ? start_progress(interval=c.progress_interval,log_interval=c.progress_log_interval,
        ascii_panel=c.ascii_panel,panel_interval=c.ascii_panel_interval,genome_range=c.genome_range,
        label=label) : nothing
    output = ""; assembly_seconds = 0.0; success = false; adaptive_note = ""; mate_note = ""; graph_note = ""; decon_note = ""
    try
        idx = build_index(c;quiet=c.show_progress,monitor=monitor)
        assembly_t0 = time_ns()
        results, raw = assemble(c, idx;quiet=c.show_progress,monitor=monitor)
        assembly_seconds = Float64(time_ns()-assembly_t0)/1e9
        if adaptive_needed(c, results; raw=raw)
            first_pairs = num_pairs(idx.store)
            first_phase_times = copy(PHASE_TIMES)
            idx = nothing; results = nothing; raw = nothing; GC.gc()
            limit = c.read_pairs_to_process
            c.read_pairs_to_process = c.adaptive_read_pairs
            try
                idx = build_index(c;quiet=c.show_progress,monitor=monitor)
                assembly_t0 = time_ns()
                results, raw = assemble(c, idx;quiet=c.show_progress,monitor=monitor)
                assembly_seconds += Float64(time_ns()-assembly_t0)/1e9
                for (key,value) in first_phase_times
                    PHASE_TIMES[key] = get(PHASE_TIMES,key,0.0) + value
                end
            finally
                c.read_pairs_to_process = limit
            end
            adaptive_note = "Adaptive read pairs: $(first_pairs) pairs gave no closed circle and path ends " *
                            "stopped for lack of reads; assembled again from $(num_pairs(idx.store)) pairs " *
                            "(limit $(c.adaptive_read_pairs))"
            println(adaptive_note)
        end
        graph_paths = [x.sequence for x in raw]              # every pass contributes its paths
        if c.mate_jumps > 0 && c.paired && !has_gap_free_circle(results)
            # Logs of the normal pass are restored if its result is kept.
            saved_logs = logs_state()
            MATE_JUMPS_ON[] = true
            jump_t0 = time_ns()
            results2, raw2 = try
                assemble(c, idx;quiet=c.show_progress,monitor=monitor)
            finally
                MATE_JUMPS_ON[] = false
            end
            assembly_seconds += Float64(time_ns()-jump_t0)/1e9
            append!(graph_paths, [x.sequence for x in raw2])
            if any(x -> x.circular, results2) && !any(x -> x.circular, results)
                results, raw = results2, raw2
                mate_note = "Mate-jump pass: no gap-free circle from normal extension; jumps over stretches " *
                            "no read crosses gave a circle with internal N gap(s) (see _mate_jumps.tsv)"
            else
                mate_note = "Mate-jump pass: tried, did not give a circle; normal result kept"
                restore_logs!(saved_logs)
                lock(() -> empty!(MATE_JUMP_LOG), MATE_JUMP_LOG_LOCK)
            end
            println(mate_note)
        end
        if c.graph_rescue && !has_gap_free_circle(results)
            reference_rescue_seeds!(graph_paths, c)
            graph_t0 = time_ns()
            GRAPH_RESCUE_REPORT[] = ""; GRAPH_RESCUE_LINKS[] = String[]; GRAPH_READ_BRIDGE_REPORT[] = String[]
            # Reference seeds (if any) in this first graph rescue only, not in decontamination rounds.
            REF_RESCUE_ACTIVE[] = !isempty(REF_RESCUE_TILES[])
            circles = try
                graph_rescue(graph_paths, idx, c;monitor=monitor)
            finally
                REF_RESCUE_ACTIVE[] = false
            end
            append!(results,filter(x->!x.circular,MATE_RESCUE_LAST_RESULTS[]))
            assembly_seconds += Float64(time_ns()-graph_t0)/1e9
            if !isempty(circles) && (!any(x -> x.circular, results) ||
                                     any(x -> isempty(n_gaps(x.sequence)), circles))
                results = vcat(circles, results)
                graph_note = "Graph rescue: $(length(circles)) circular walk(s) of the path graph with node " *
                             "copy numbers from coverage (closure $(closure_status(circles[1]))); see _graph.gfa"
            else
                graph_note = "Graph rescue: no plausible circular walk of the path graph"
            end
            graph_note *= " [" * GRAPH_RESCUE_REPORT[] * "]"
            isempty(REF_RESCUE_NOTE[]) || (graph_note *= "\n" * REF_RESCUE_NOTE[])
            println(graph_note)
        end
        if !has_gap_free_circle(results) && !isempty(c.plastid_panel) && GRAPH_RESCUE_GRAPH[] !== nothing
            nflag, fbp = flag_foreign_reads!(GRAPH_RESCUE_GRAPH[], idx, c)
            if nflag > 0
                saved_state = export_state()             # restored if the previous result is kept
                decon_t0 = time_ns()
                rounds = String[]
                # A bacterium with many rRNA operons offers a new flank once the first is
                # filtered: up to 4 rounds, each flagging the foreign nodes of the last graph.
                # Output tier of the best circle in a result set (9 = no circle; lower is better).
                best_tier = rs -> minimum((result_tier(x, c) for x in rs if x.circular); init=9)
                results5, raw5 = try
                    r5 = AssemblyResult[]; w5 = AssemblyResult[]
                    kept5 = nothing                      # round with the best weak circle so far
                    for rnd in 1:4
                        r5, w5 = assemble(c, idx;quiet=c.show_progress,monitor=monitor)
                        if !has_gap_free_circle(r5) && c.graph_rescue
                            GRAPH_RESCUE_REPORT[] = ""
                            circ5 = graph_rescue([x.sequence for x in w5], idx, c;monitor=monitor)
                            r5 = vcat(circ5,r5,filter(x->!x.circular,MATE_RESCUE_LAST_RESULTS[]))
                        end
                        longest = isempty(r5) ? 0 : maximum(length(x.sequence) for x in r5)
                        pf = isempty(r5) ? NaN : plastid_fraction(sort(r5; by=x -> -length(x.sequence))[1].sequence,
                                                                  load_plastid_panel(c.plastid_panel))
                        push!(rounds, "round $(rnd): $(count(FOREIGN_READS[])) reads filtered, longest $(longest) bp " *
                                      "(plastid 21-mers $(round(pf; digits=3))), circle $(any(x -> x.circular, r5)) " *
                                      "[$(GRAPH_RESCUE_REPORT[])]")
                        # Only a conclusive circle ends the rounds; a weak one is kept in case
                        # later rounds find nothing better.
                        has_gap_free_circle(r5) && break
                        if any(x -> x.circular, r5) && (kept5 === nothing || best_tier(r5) < best_tier(kept5.r5))
                            kept5 = (r5=r5, w5=w5, round=rnd, state=export_state())
                        end
                        GRAPH_RESCUE_GRAPH[] === nothing && break
                        before = count(FOREIGN_READS[])
                        flag_foreign_reads!(GRAPH_RESCUE_GRAPH[], idx, c)
                        count(FOREIGN_READS[]) > before || break       # nothing new to filter
                    end
                    if !has_gap_free_circle(r5) && kept5 !== nothing && best_tier(kept5.r5) < best_tier(r5)
                        r5, w5 = kept5.r5, kept5.w5
                        restore_export_state!(kept5.state)              # GFA, evidence and logs of that round
                        push!(rounds, "round $(kept5.round) restored: its circle is better supported than the later rounds' results")
                    end
                    r5, w5
                finally
                    FOREIGN_READS[] = BitVector()
                end
                assembly_seconds += Float64(time_ns()-decon_t0)/1e9
                # Same ordering as the output (result_tier): a better-supported circle wins,
                # also over a weaker circle of the previous pass.
                better = any(x -> x.circular, results5) && best_tier(results5) < best_tier(results)
                # No circle either way: the decontaminated paths are still the better output, since
                # the previous ones run from the plastome into the contaminant (chimeras).
                cleaner = !better && !any(x -> x.circular, results) && !isempty(results5)
                if cleaner
                    results, raw = results5, raw5
                    decon_note = "Decontamination pass: $(nflag) reads flagged as contaminant from $(fbp) bp of " *
                                 "non-plastid graph nodes; no circle, decontaminated linear paths reported instead " *
                                 "of paths running into the contaminant"
                elseif better
                    results, raw = results5, raw5
                    decon_note = "Decontamination pass: $(nflag) reads flagged as contaminant from $(fbp) bp of " *
                                 "non-plastid graph nodes; reassembly without them gave a circle " *
                                 "(closure $(closure_status(first(filter(x -> x.circular, results5)))))"
                else
                    restore_export_state!(saved_state)           # graph, evidence and logs of the kept result
                    decon_note = "Decontamination pass: $(nflag) reads flagged from $(fbp) bp of non-plastid nodes; " *
                                 "no better circle, previous result kept"
                end
                decon_note *= "\n" * join(("  " * r for r in rounds), "\n")
                println(decon_note)
            end
        end
        scaffold_note = ""
        if !any(x -> x.circular, results)
            sc, scaffold_note = reference_scaffold_rescue(results, c)
            sc === nothing || pushfirst!(results, sc)
            isempty(scaffold_note) || println(scaffold_note)
            isempty(scaffold_note) || (graph_note = isempty(graph_note) ? scaffold_note : graph_note * "\n" * scaffold_note)
        end
        export_t0 = time_ns()
        output = write_results(c, idx, results, raw, Float64(time_ns()-t0)/1e9;monitor=monitor)
        PHASE_TIMES["export"] = Float64(time_ns()-export_t0)/1e9
        success = true
    finally
        stop_progress!(monitor;success=success)
    end
    total_seconds = Float64(time_ns()-t0)/1e9
    rss = resident_memory_bytes()
    sampled_peak = monitor===nothing ? nothing : monitor.peak_rss
    @printf("Assembly time: %s (%.2f s) | Total time: %s (%.2f s)\n",
        duration_text(assembly_seconds),assembly_seconds,duration_text(total_seconds),total_seconds)
    rss===nothing ? println("Current RAM (RSS): n/a") : @printf("Current RAM (RSS): %.3f GB\n",rss/1e9)
    st = READ_STATS[]
    st === nothing || @printf("Plastid reads (approx.): %d of %d (%.2f%%) | mean coverage %.1fx\n",
                              st.reads, st.total, st.percent, st.coverage)
    if sampled_peak!==nothing && sampled_peak>0
        @printf("Peak sampled RAM (RSS): %.3f GB\n",sampled_peak/1e9)
    end
    open(joinpath(c.output_path,c.project_name*"_run.txt"),"a") do io
        @printf(io,"Assembly seconds (seeding, extension, deduplication, scaffolding): %.2f\n",assembly_seconds)
        @printf(io,"Total seconds (including indexing and export): %.2f\n",total_seconds)
        for key in ("input_wait", "parse", "read_store", "minimizers", "reserve", "gc_in_minimizers",
                    "gc_reading_and_minimizers", "index_freeze", "export")
            haskey(PHASE_TIMES, key) && @printf(io,"Phase seconds %s: %.2f\n",key,PHASE_TIMES[key])
        end
        println(io,"Current RSS GB: ",rss===nothing ? "n/a" : @sprintf("%.6f",rss/1e9))
        println(io,"Peak sampled RSS GB: ",sampled_peak===nothing || sampled_peak==0 ? "n/a" : @sprintf("%.6f",sampled_peak/1e9))
        println(io,"Memory unit: 1 GB = 1000000000 bytes; RSS for this Julia process, all threads.")
        isempty(adaptive_note) || println(io, adaptive_note)
        isempty(mate_note) || println(io, mate_note)
        isempty(graph_note) || println(io, graph_note)
        isempty(decon_note) || println(io, decon_note)
    end
    println("Saved: ", output)
    peak = sampled_peak===nothing ? 0 : sampled_peak
    return (output=output, assembly_seconds=assembly_seconds, total_seconds=total_seconds, peak_rss=peak)
end

function main(args=ARGS)
    if length(args) == 1 && args[1] == "--self-test"
        testfile = joinpath(@__DIR__, "tests_v16.jl")
        isfile(testfile) || error("tests_v16.jl must be beside this script; download the complete package")
        Base.include(Main, testfile)
        return
    end
    if length(args) == 2 && args[1] == "--build-reference-sketches"
        RefSketch.build(args[2])
        return
    end
    if length(args) in (4, 5) && args[1] == "--orient"
        # Re-order an existing _assembly.fasta (no reassembly): --orient <dataset> <fasta> <prefix> [plastid panel]
        # writes <prefix>_plastome_LSC_IR_SSC_IR.fasta (+ regions, ordering) and <prefix>_nearest_reference.*
        dataset, fasta, prefix = args[2], args[3], args[4]
        c = Config(); c.reference_dataset = dataset
        length(args) == 5 && (c.plastid_panel = abspath(args[5]))
        recs = AssemblyResult[]; name = ""; buf = UInt8[]; hdr = ""
        function push_rec!()
            isempty(name) && return
            m = match(r"closure=(\S+)", hdr); cl = m === nothing ? "validated" : m.captures[1]
            reason = cl == "validated" ? "validated_closure" : cl == "gap_unknown" ? "closed_gap_unknown_length" : "closed_" * cl
            circ = occursin("topology=circular", hdr)
            push!(recs, AssemblyResult(name, buf, circ, circ ? reason : "linear", "external", 0, 0, 0))
        end
        for line in eachline(fasta)
            if startswith(line, ">")
                push_rec!(); hdr = line; name = String(first(split(line[2:end]))); buf = UInt8[]
            else
                append!(buf, codeunits(strip(line)))
            end
        end
        push_rec!()
        mkpath(dirname(abspath(prefix)))
        ref = write_nearest_reference(c, recs, prefix)
        println("Nearest reference: ", ref.summary)
        n = write_plastome_outputs(c, recs, prefix; reference=ref)
        println("Ordered: ", n, " -> ", prefix, "_plastome_LSC_IR_SSC_IR.fasta")
        return
    end
    if length(args) in (3, 4) && args[1] == "--annotate"
        # Annotate an existing ordered plastome: --annotate <annotation dir> <prefix> [annotation dataset]
        # reads <prefix>_plastome_LSC_IR_SSC_IR.fasta, writes .gff3 and <prefix>_plastome_genes.tsv
        c = Config(); c.annotation_dir = abspath(args[2])
        length(args) == 4 && (c.annotation_dataset = abspath(args[4]))
        n, note = annotate_plastomes(c, args[3])
        println("Annotated: ", n, isempty(note) ? "" : " ($note)")
        return
    end
    if length(args) == 3 && args[1] == "--nearest-reference"
        # Check an existing assembly FASTA: first circular record, else all records >= 5 kb.
        dataset, fasta = args[2], args[3]
        recs = Tuple{String,Vector{UInt8},Bool}[]
        name = ""; buf = UInt8[]; circ = false
        for line in eachline(fasta)
            if startswith(line, ">")
                isempty(name) || push!(recs, (name, buf, circ))
                name = line[2:end]; buf = UInt8[]; circ = occursin("topology=circular", line)
            else
                append!(buf, codeunits(strip(line)))
            end
        end
        isempty(name) || push!(recs, (name, buf, circ))
        i = findfirst(r -> r[3], recs)
        query = i === nothing ? [r for r in recs if length(r[2]) >= 5000] : [recs[i]]
        db = RefSketch.load(joinpath(dataset, RefSketch.SKETCH_FILE)); meta = RefSketch.metadata(dataset)
        hits, nq = RefSketch.nearest(db, [r[2] for r in query], [r[3] for r in query])
        println("query: ", length(query), " record(s), ", nq, " hashes")
        for h in hits
            org, grp = get(meta, h.accession, ("", ""))
            @printf("%-14s %-40s %-12s containment %.3f  ANI~%.3f  reference covered %.3f\n", h.accession,
                    first(org, 40), grp, h.containment, h.ani, h.reference_covered)
        end
        return
    end
    if length(args) != 1 || args[1] in ("-h", "--help")
        println("Usage: julia --project=. --threads auto AssemblerComplete_v16.jl config_v16.txt")
        println("Reference sketches: julia --project=. AssemblerComplete_v16.jl --build-reference-sketches <dataset dir>")
        println("Annotation only: julia --project=. --threads auto AssemblerComplete_v16.jl --annotate <annotation dir> <prefix> [annotation dataset]")
        println("Tests: julia --project=. --threads 4 AssemblerComplete_v16.jl --self-test")
        return
    end
    run_assembly(parse_config(args[1]))
    return
end

end # module AssemblerCompleteV16

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    AssemblerCompleteV16.main()
end
