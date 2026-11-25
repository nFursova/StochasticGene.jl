# Stage B for the JOINT fits of 14_stage_trace_joint_gene.jl (shared-parameter path):
# freeze ALL noise parameters from each model's own fit and refit only the shared kinetic rates.
#
# (13_continue_trace_fixnoise.jl is for the separate fits of the 12 scripts, which use the
#  classic hierarchical path; it refuses these shared-parameter runs.)
#
# For every finished run in SG_FROM:
#   * noise means per trace and noise SDs per sample (every column whose label contains
#     "noiseparam") are copied from row SG_RATETYPE (median by default) into all four rows of a
#     new starting-rate file, so no later ratetype choice can change them;
#   * fittedparam keeps only the kinetic rates (same priors, same Rshift ties), so the noise
#     stays at those values;
#   * the kinetic rates start from the same row.
#
# Run:
#   export SG_FROM=3and5Prime-joint-gene-<stamp>
#   export SG_INTO=$SG_FROM-fixnoise-24h            (default)
#   julia --project=$SG_PROJECT scripts/15_continue_trace_joint_fixnoise.jl
#
# Environment (defaults): SG_ROOT, SG_PROJECT=/data/$USER/sg201, SG_RATETYPE=median,
#   SG_MAXTIME=86400, SG_NCHAINS=16, SG_ONLY (comma list of key substrings), SG_SWARMDIR=swarm-<SG_INTO>

using StochasticGene

const ROOT = abspath(get(ENV, "SG_ROOT", "/data/Larsonlab/Nadya/StochasticGene"))
const PROJECT = get(ENV, "SG_PROJECT", "/data/$(get(ENV, "USER", "unknown"))/sg201")
const FROM = strip(get(ENV, "SG_FROM", ""))
isempty(FROM) && error("SG_FROM (results folder of 14_stage_trace_joint_gene.jl) must be set")
const INTO = strip(get(ENV, "SG_INTO", FROM * "-fixnoise-24h"))
FROM == INTO && error("SG_FROM and SG_INTO must differ")
const RATETYPE = strip(get(ENV, "SG_RATETYPE", "median"))
RATETYPE in ("ml", "mean", "median", "last") || error("SG_RATETYPE must be ml|mean|median|last")
const MAXTIME = parse(Float64, get(ENV, "SG_MAXTIME", "86400"))
const NCHAINS = parse(Int, get(ENV, "SG_NCHAINS", "16"))
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
    sp = get(spec, :shared_parameters, nothing)
    sp === nothing && error("$key: not a shared-parameter run; use 13_continue_trace_fixnoise.jl for the separate fits")
    ratefile = joinpath(SRCDIR, "rates_" * key * ".txt")
    isfile(ratefile) || (push!(skipped, "$key (no rates file: run not finished?)"); continue)

    lines = readlines(ratefile)
    length(lines) >= 5 || error("$key: $ratefile should have a header and 4 rows (ml, mean, median, last)")
    hdr = String.(split(lines[1], ','))
    noisecols = findall(h -> occursin("noiseparam", h), hdr)
    isempty(noisecols) && error("$key: no noiseparam columns in $ratefile")
    ntrace = count(h -> startswith(h, "noiseparam_1_") && occursin("_trace", h), hdr)

    ref = parse.(Float64, split(lines[1 + StochasticGene.get_row(RATETYPE)], ','))
    length(ref) == length(hdr) || error("$key: row $RATETYPE has $(length(ref)) values for $(length(hdr)) columns")
    mkpath(INITDIR)
    out = joinpath(INITDIR, "rates_" * key * ".txt")
    open(out, "w") do io
        println(io, lines[1])
        for row in 1:4
            v = parse.(Float64, split(lines[1+row], ','))
            v[noisecols] = ref[noisecols]
            all(isfinite, v) || error("$key: NaN/Inf in row $row")
            println(io, join(v, ","))
        end
    end

    n = num_rates(spec[:transitions], spec[:R], spec[:S], spec[:insertstep])
    kinetic = sort([i for i in Int.(spec[:fittedparam]) if i <= n])
    spec[:resultfolder] = INTO
    spec[:maxtime] = MAXTIME
    spec[:nchains] = NCHAINS
    spec[:ratetype] = RATETYPE
    spec[:init_ratefile] = abspath(out)
    spec[:fittedparam] = kinetic
    push!(specs, spec)
    push!(report, (; key, ntrace, kinetic, ties=spec[:fixedeffects], nfrozen=length(noisecols), ncols=length(hdr)))
end
isempty(specs) && error("nothing to stage from $SRCDIR (check SG_FROM / SG_ONLY)")

println("from      : $SRCDIR")
println("into      : $DESTDIR")
println("frozen    : all noise parameters (means per trace, SDs per sample), row = $RATETYPE")
for r in report
    println("  ", rpad(r.key, 34), " traces=", lpad(r.ntrace, 4), " fitted kinetic=", r.kinetic,
            " ties=", r.ties, " frozen=", r.nfrozen, " of ", r.ncols, " columns")
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
    for k in (:fittedparam, :init_ratefile, :maxtime, :priormean, :priorcv, :fixedeffects, :ratetype)
        (haskey(got, k) && got[k] == spec[k]) || push!(problems, "$(r.key): $k not as staged")
    end
    haskey(got, :shared_parameters) || push!(problems, "$(r.key): shared_parameters lost")
    n = num_rates(got[:transitions], got[:R], got[:S], got[:insertstep])
    any(>(n), got[:fittedparam]) && push!(problems, "$(r.key): a noise parameter is still fitted")
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
println("submit    : cd $SWARMDIR && swarm -f fit.swarm -g 48 -t $(NCHAINS + 2) --time $(lpad(wall, 2, '0')):00:00 --module julialang --logdir logs")
