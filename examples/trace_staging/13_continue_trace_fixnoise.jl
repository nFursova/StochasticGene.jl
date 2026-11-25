# Stage B: continue the stage-A fits with ALL noise parameters frozen; fit only kinetic rates.
#
# Reads every info_<key>.toml / rates_<key>.txt in SG_FROM (written by
# 12_stage_trace_separate_noise.jl, classic hierarchical path) and, for each model, from ITS OWN
# fit:
#   * freezes the per-trace noise (offsets 1, 3) in every trace block,
#   * freezes the per-sample noise (offsets 2, 4),
#   * keeps fitting the kinetic rates with the same priors and the same ties
#     (all Rshift equal, all splice rates equal, Eject separate).
# Frozen values come from row SG_RATETYPE (median by default) of rates_<key>.txt and are written
# into all four rows of a new starting-rate file, so no later ratetype choice can change them.
# Kinetic rates start from the same row.
#
# Run:
#   export SG_FROM=trace-separate-enhancer-noise-<tag>
#   export SG_INTO=$SG_FROM-fixnoise-24h          (default)
#   julia --project=$SG_PROJECT scripts/13_continue_trace_fixnoise.jl
#
# Environment (defaults): SG_ROOT, SG_PROJECT=/data/$USER/sg201, SG_RATETYPE=median,
#   SG_MAXTIME=86400, SG_NCHAINS=16, SG_FIXNOISE_PROPCV=0.05, SG_ONLY (comma list of key
#   substrings), SG_SWARMDIR=swarm-<SG_INTO>

using StochasticGene

const ROOT = abspath(get(ENV, "SG_ROOT", "/data/Larsonlab/Nadya/StochasticGene"))
const PROJECT = get(ENV, "SG_PROJECT", "/data/$(get(ENV, "USER", "unknown"))/sg201")
const FROM = strip(get(ENV, "SG_FROM", ""))
isempty(FROM) && error("SG_FROM (stage-A results folder) must be set")
const INTO = strip(get(ENV, "SG_INTO", FROM * "-fixnoise-24h"))
FROM == INTO && error("SG_FROM and SG_INTO must differ")
const RATETYPE = strip(get(ENV, "SG_RATETYPE", "median"))
RATETYPE in ("ml", "mean", "median", "last") || error("SG_RATETYPE must be ml|mean|median|last")
const MAXTIME = parse(Float64, get(ENV, "SG_MAXTIME", "86400"))
const NCHAINS = parse(Int, get(ENV, "SG_NCHAINS", "16"))
const PROPCV = parse(Float64, get(ENV, "SG_FIXNOISE_PROPCV", "0.05"))
const ONLY = [strip(x) for x in split(get(ENV, "SG_ONLY", ""), ',') if !isempty(strip(x))]
const SWARMDIR = joinpath(ROOT, get(ENV, "SG_SWARMDIR", "swarm-" * INTO))

const SRCDIR = joinpath(ROOT, "results", FROM)
const DESTDIR = joinpath(ROOT, "results", INTO)
const INITDIR = joinpath(ROOT, "results", INTO * "-initrates")
isdir(SRCDIR) || error("source results folder not found: $SRCDIR")
if isdir(DESTDIR)
    stale = filter(f -> startswith(f, "rates_") && endswith(f, ".txt"), readdir(DESTDIR))
    isempty(stale) || error("$DESTDIR already has rates files; use a fresh SG_INTO")
end

function resolve_data(path)
    for p in (path, joinpath(ROOT, path), joinpath(ROOT, "data", path))
        isdir(p) && return p
    end
    error("trace folder not found: $path")
end

specs = Dict{Symbol,Any}[]
report = NamedTuple[]
skipped = String[]
for f in sort(readdir(SRCDIR))
    (startswith(f, "info_") && endswith(f, ".toml")) || continue
    toml = joinpath(SRCDIR, f)
    isfile(StochasticGene.info_jld2_path(toml)) || (push!(skipped, "$f (no .jld2)"); continue)
    spec = StochasticGene.read_run_spec(toml)
    key = String(spec[:key])
    isempty(ONLY) || any(o -> occursin(o, key), ONLY) || continue
    ratefile = joinpath(SRCDIR, "rates_" * key * ".txt")
    isfile(ratefile) || (push!(skipped, "$key (no rates file: run not finished?)"); continue)

    h = spec[:hierarchical]
    haskey(spec, :shared_parameters) && spec[:shared_parameters] !== nothing &&
        error("$key: this is a shared-parameter (joint) run; use 15_continue_trace_joint_fixnoise.jl")
    (h isa Tuple && length(h) == 3) || error("$key: expected classic hierarchical=(nhyper, individual, fixed)")
    nhyper = Int(h[1]); individual = Int.(collect(h[2])); fixedind = h[3]
    n = num_rates(spec[:transitions], spec[:R], spec[:S], spec[:insertstep])
    nall = n + length(spec[:noisepriors])
    noise = collect(n+1:nall)

    lines = readlines(ratefile)
    hdr = split(lines[1], ',')
    rem(length(hdr), nall) == 0 || error("$key: $(length(hdr)) columns is not a multiple of block size $nall")
    nblocks = length(hdr) ÷ nall
    ntr = nblocks - nhyper
    nload = length(StochasticGene.read_tracefiles(resolve_data(String(spec[:datapath])), String(spec[:datacond]), 1, -1, 3))
    nload == ntr || error("$key: rates file has $ntr trace blocks but $nload traces load from $(spec[:datapath])")

    # starting-rate file: noise columns of every block set to the chosen row in all four rows
    ref = parse.(Float64, split(lines[1 + StochasticGene.get_row(RATETYPE)], ','))
    frozen = [(b - 1) * nall + j for b in 1:nblocks for j in noise]
    mkpath(INITDIR)
    out = joinpath(INITDIR, "rates_" * key * ".txt")
    open(out, "w") do io
        println(io, lines[1])
        for row in 1:4
            v = parse.(Float64, split(lines[1+row], ','))
            v[frozen] = ref[frozen]
            all(isfinite, v) || error("$key: NaN/Inf in row $row")
            println(io, join(v, ","))
        end
    end

    kinetic = sort(setdiff(Int.(spec[:fittedparam]), noise))
    spec[:resultfolder] = INTO
    spec[:maxtime] = MAXTIME
    spec[:nchains] = NCHAINS
    spec[:ratetype] = RATETYPE
    spec[:init_ratefile] = abspath(out)
    spec[:fittedparam] = kinetic
    spec[:hierarchical] = (nhyper, Int[], fixedind)
    spec[:propcv] = PROPCV
    for k in (:propcv_rate, :propcv_noise, :proposal_cv_rates, :proposal_cv_noise)
        spec[k] = 0.0
    end
    for k in (:propcv_levels, :proposal_cv_levels)
        spec[k] = nothing
    end
    push!(specs, spec)
    push!(report, (; key, ntr, kinetic, ties=spec[:fixedeffects], nfrozen=length(frozen), out))
end
isempty(specs) && error("nothing to stage from $SRCDIR (check SG_FROM / SG_ONLY)")

println("from      : $SRCDIR")
println("into      : $DESTDIR")
println("frozen    : all noise parameters (per trace and per sample), row = $RATETYPE")
for r in report
    println("  ", rpad(r.key, 36), " traces=", lpad(r.ntr, 4), " fitted kinetic=", r.kinetic,
            " ties=", r.ties, " frozen values=", r.nfrozen)
end
isempty(skipped) || (println("skipped:"); foreach(s -> println("  - ", s), skipped))

out = StochasticGene.stage_write_run_specs(
    specs; resultfolder=INTO, root=ROOT, filedir=SWARMDIR,
    juliafile="fitscript", swarmfile="fit", commandfile="fit.commands",
    nchains=NCHAINS, nthreads=1, project=PROJECT,
)

problems = String[]
for (spec, r) in zip(specs, report)
    got = StochasticGene.read_run_spec(joinpath(DESTDIR, "info_" * r.key * ".toml"))
    for k in (:fittedparam, :hierarchical, :init_ratefile, :maxtime, :priormean, :priorcv, :fixedeffects)
        (haskey(got, k) && got[k] == spec[k]) || push!(problems, "$(r.key): $k not as staged")
    end
    isempty(got[:hierarchical][2]) || push!(problems, "$(r.key): per-trace parameters still fitted")
    any(i -> i > num_rates(got[:transitions], got[:R], got[:S], got[:insertstep]), got[:fittedparam]) &&
        push!(problems, "$(r.key): a noise parameter is still fitted")
end
isempty(problems) || error("staged specs are wrong, do not submit:\n  " * join(problems, "\n  "))
println("spec check: all ", length(specs), " specs ok")

for f in ("fit.swarm", "fit.commands")
    p = joinpath(SWARMDIR, f); isfile(p) || continue
    lines = [startswith(l, "julia ") ? "JULIA_WORKER_TIMEOUT=600 JULIA_PKG_PRECOMPILE_AUTO=0 " * l : l for l in readlines(p)]
    write(p, join(lines, "\n") * "\n")
end
mkpath(joinpath(SWARMDIR, "logs"))
wall = ceil(Int, MAXTIME / 3600) + 3
println("submit    : cd $SWARMDIR && swarm -f fit.swarm -g 32 -t $(NCHAINS + 2) --time $(lpad(wall, 2, '0')):00:00 --module julialang --logdir logs")
