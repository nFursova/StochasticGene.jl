# Stage A: separate (unshared) trace fits for 3Prime and 5Prime, per-trace noise means.
#
# Each sample is its own fit: nothing is shared between 3Prime and 5Prime. Within one sample:
#
#   hierarchy   classic hierarchical=(2, individual, ()) : population block + per-trace blocks
#   noise       offsets 1, 3 (background mean, reporter mean)  -> fitted PER TRACE
#               offsets 2, 4 (background SD, reporter SD)      -> fitted once per sample
#   kinetics    one value per sample (shared by all traces of that sample)
#   ties        all Rshift rates equal (Rshift_1 fitted, the others tied to it);
#               Eject separate; all splice rates equal
#   decay       not fitted
#
# Priors (rates per minute):
#   G transitions  0.01, CV 1.0            Initiate 0.2, CV 0.25
#   Rshift         1/SG_RSHIFT_MINUTES (default 1/2.5 min), CV 0.25
#   Eject          SG_EJECT_RATE (default 1.0 /min), CV 0.25
#   splice         0.05, CV 0.2            noise: per-sample priors from stage 05, CV 0.1
#
# Models (name = G R S insertstep topology): every base in SG_BASES x every topology in
# SG_TOPOLOGIES, i.e. by default 9 bases x 3 = 27 models per sample:
#   bases       3301 3302 3303 3401 3402 3403 3404 3331 3441
#   topologies  none = linear 1<->2<->3
#               cyc  = linear plus 3<->1 (fully reversible)
#               kp   = 1<->2, 2->3, 3->1 (kinetic proofreading)
# SG_MODELS (full names, e.g. 3301kp,3441cyc) overrides the bases x topologies list.
# Splicing is staged only with insertstep 1: with insertstep > 1, fit() resets S (reset_S) and
# the splice rates never enter the likelihood. The script checks every model and refuses
# dead splice rates.
#
# Run length: SG_MAXTIME (default 43200 s = 12 h) is the stop condition; SG_SAMPLESTEPS
# (default 10^12) is never reached.
#
# Then run 13_continue_trace_fixnoise.jl on the result folder to freeze the noise and fit only
# the kinetic rates.
#
# Environment (defaults):
#   SG_ROOT, SG_PROJECT=/data/$USER/sg201, SG_CONDITION=enhancer, SG_SAMPLES_LIST=3Prime,5Prime
#   SG_DATA_3PRIME=traces/3Prime, SG_DATA_5PRIME=traces/5Prime
#   SG_BASES=3301,3302,3303,3401,3402,3403,3404,3331,3441  SG_TOPOLOGIES=none,cyc,kp
#   SG_MODELS (overrides bases x topologies), SG_TAG=<date>, SG_RESULTFOLDER=trace-separate-<condition>-noise-<tag>,
#   SG_SWARMDIR=swarm-<resultfolder>
#   SG_MAXTIME=43200 SG_SAMPLESTEPS=10^12 SG_SAMPLE_STRIDE=1000 SG_NCHAINS=16
#   SG_PROPCV_RATE=0.05 SG_PROPCV_NOISE=0.001 SG_WRITESAMPLES=true
#   SG_RSHIFT_MINUTES=2.5 SG_PRIORCV_R=0.25 SG_EJECT_RATE=1.0 SG_PRIORCV_EJECT=0.25
#   SG_PRIOR_G=0.01 SG_PRIORCV_G=1.0 SG_PRIOR_INITIATE=0.2 SG_PRIORCV_INITIATE=0.25
#   SG_SPLICE_PRIOR=0.05 SG_PRIORCV_SPLICE=0.2 SG_NOISE_INDIVIDUAL=1,3 SG_NOISE_SHARED=2,4

using StochasticGene
using Dates

envstr(v, d) = get(ENV, v, d)
envint(v, d) = parse(Int, get(ENV, v, string(d)))
envfloat(v, d) = parse(Float64, get(ENV, v, string(d)))
envbool(v, d) = lowercase(get(ENV, v, string(d))) in ("1", "true", "yes", "on")
csvlist(v, d) = [String(strip(s)) for s in split(get(ENV, v, d), ',') if !isempty(strip(s))]

const USER = get(ENV, "USER", "unknown")
const ROOT = abspath(envstr("SG_ROOT", "/data/Larsonlab/Nadya/StochasticGene"))
const PROJECT = envstr("SG_PROJECT", "/data/$USER/sg201")
const CONDITION = envstr("SG_CONDITION", "enhancer")
CONDITION in ("gene", "enhancer") || error("SG_CONDITION must be gene or enhancer")
const SAMPLES = csvlist("SG_SAMPLES_LIST", "3Prime,5Prime")
const TAG = envstr("SG_TAG", string(Dates.today()))
const RESULTFOLDER = envstr("SG_RESULTFOLDER", "trace-separate-$CONDITION-noise-$TAG")
const SWARMDIR = joinpath(ROOT, envstr("SG_SWARMDIR", "swarm-" * RESULTFOLDER))

const MAXTIME = envfloat("SG_MAXTIME", 43200.0)
const SAMPLESTEPS = envint("SG_SAMPLESTEPS", 10^12)
const STRIDE = envint("SG_SAMPLE_STRIDE", 1000)
const NCHAINS = envint("SG_NCHAINS", 16)
const PROPCV_RATE = envfloat("SG_PROPCV_RATE", 0.05)
const PROPCV_NOISE = envfloat("SG_PROPCV_NOISE", 0.001)
const WRITESAMPLES = envbool("SG_WRITESAMPLES", true)
const INTERVAL = 1.6666666666666667

const PRIOR_G = envfloat("SG_PRIOR_G", 0.01)
const PRIORCV_G = envfloat("SG_PRIORCV_G", 1.0)
const PRIOR_INITIATE = envfloat("SG_PRIOR_INITIATE", 0.2)
const PRIORCV_INITIATE = envfloat("SG_PRIORCV_INITIATE", 0.25)
const RSHIFT_MINUTES = envfloat("SG_RSHIFT_MINUTES", 2.5)
const PRIORCV_R = envfloat("SG_PRIORCV_R", 0.25)
const EJECT_RATE = envfloat("SG_EJECT_RATE", 1.0)
const PRIORCV_EJECT = envfloat("SG_PRIORCV_EJECT", 0.25)
const SPLICE_PRIOR = envfloat("SG_SPLICE_PRIOR", 0.05)
const PRIORCV_SPLICE = envfloat("SG_PRIORCV_SPLICE", 0.2)
const DECAY_PLACEHOLDER = 0.03165055618995184

# Noise offsets past the last rate: 1 background mean, 2 background SD, 3 reporter mean, 4 reporter SD
const NOISE_INDIVIDUAL = parse.(Int, csvlist("SG_NOISE_INDIVIDUAL", "1,3"))
const NOISE_SHARED = parse.(Int, csvlist("SG_NOISE_SHARED", "2,4"))
isempty(intersect(NOISE_INDIVIDUAL, NOISE_SHARED)) || error("a noise offset cannot be both per-trace and shared")
all(in(1:4), vcat(NOISE_INDIVIDUAL, NOISE_SHARED)) || error("noise offsets must be in 1:4")

# Per-sample noise priors (stage 05 values)
const NOISE_PRIORS = Dict(
    ("3Prime", "gene") => [-0.018007, 0.214861, 1.256466, 0.405046],
    ("5Prime", "gene") => [-0.004966, 0.118037, 0.607264, 0.223634],
    ("3Prime", "enhancer") => [-0.015733, 0.088230, 0.532594, 0.198461],
    ("5Prime", "enhancer") => [-0.011295, 0.053569, 0.322406, 0.124341],
)
datapath(sample) = envstr("SG_DATA_" * uppercase(sample), "traces/$sample")

const TRANSITIONS = Dict(
    (2, "none") => ([1, 2], [2, 1]),
    (3, "none") => ([1, 2], [2, 1], [2, 3], [3, 2]),
    (3, "cyc") => ([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1]),
    (3, "kp") => ([1, 2], [2, 1], [2, 3], [3, 1]),
)

function model_spec(name)
    m = match(r"^(\d)(\d)(\d)(\d)([a-z]*)$", name)
    m === nothing && error("model name must be GRSI<topology>, got $name")
    G, R, S, I = parse.(Int, m.captures[1:4])
    topo = isempty(m.captures[5]) ? "none" : m.captures[5]
    haskey(TRANSITIONS, (G, topo)) || error("$name: no transitions for G=$G topology=$topo")
    1 <= I <= R || error("$name: insertstep must be in 1:R")
    return (name="$(G)$(R)$(S)$(I)$(topo)", G=G, R=R, S=S, insertstep=I, transitions=TRANSITIONS[(G, topo)])
end

const BASES = csvlist("SG_BASES", "3301,3302,3303,3401,3402,3403,3404,3331,3441")
const TOPOLOGIES = csvlist("SG_TOPOLOGIES", "none,cyc,kp")
const MODELS = [model_spec(n) for n in (haskey(ENV, "SG_MODELS") ? csvlist("SG_MODELS", "") :
                                         [b * t for b in BASES for t in TOPOLOGIES])]
isempty(MODELS) && error("no models selected")
length(unique(m.name for m in MODELS)) == length(MODELS) || error("duplicate models")

nrates(m) = num_rates(m.transitions, m.R, m.S, m.insertstep)
nsplice(m) = m.S == 0 ? 0 : min(m.S, m.R - m.insertstep + 1)

# Rate positions (G | Initiate | Rshift_1..Rshift_{R-1} | Eject | splice | Decay | noise 1-4)
function layout(m)
    nt = length(m.transitions); n = nrates(m)
    (G=collect(1:nt), initiate=nt + 1, rshift=collect(nt+2:nt+m.R), eject=nt + 1 + m.R,
     splice=collect(nt+2+m.R:nt+1+m.R+nsplice(m)), decay=n, noise=collect(n+1:n+4), n=n, nall=n + 4)
end

# Refuse splice rates that would never enter the likelihood (the reset_S / builder mismatch).
function check_splice_live(m)
    m.S == 0 && return
    S_eff = StochasticGene.reset_S(m.S, m.R, m.insertstep)
    idx = StochasticGene.set_indices(length(m.transitions), m.R, S_eff, m.insertstep)
    elements, _ = StochasticGene.set_elements_TGRS(m.transitions, m.G, m.R, S_eff, m.insertstep, idx, "")
    used = Set(e.index for e in elements)
    dead = [i for i in layout(m).splice if !(i in used)]
    isempty(dead) || error("""
        $(m.name): splice rate(s) $dead never enter the likelihood (fit() resets S to $S_eff and
        the matrix builder only splices at steps insertstep..S). Use insertstep 1 for splicing.""")
end

function priors(m, sample)
    L = layout(m)
    mean = zeros(L.nall); cv = zeros(L.nall)
    mean[L.G] .= PRIOR_G; cv[L.G] .= PRIORCV_G
    mean[L.initiate] = PRIOR_INITIATE; cv[L.initiate] = PRIORCV_INITIATE
    mean[L.rshift] .= 1 / RSHIFT_MINUTES; cv[L.rshift] .= PRIORCV_R
    mean[L.eject] = EJECT_RATE; cv[L.eject] = PRIORCV_EJECT
    mean[L.splice] .= SPLICE_PRIOR; cv[L.splice] .= PRIORCV_SPLICE
    mean[L.decay] = DECAY_PLACEHOLDER; cv[L.decay] = 0.1
    mean[L.noise] = NOISE_PRIORS[(sample, CONDITION)]; cv[L.noise] .= 0.1
    # Second hyperparameter block (population spread of the per-trace parameters), the
    # StochasticGene default values: 1.0 for rates, 0.1 for noise. Only the entries of the
    # per-trace parameters are fitted; their prior CV is 1.0 (weak).
    hyper = vcat(fill(1.0, L.n), fill(0.1, 4))
    return vcat(mean, hyper), vcat(cv, fill(1.0, L.nall))
end

function fitted_and_ties(m)
    L = layout(m)
    ties = Vector{Int}[]
    length(L.rshift) > 1 && push!(ties, L.rshift)          # all Rshift rates equal
    length(L.splice) > 1 && push!(ties, L.splice)          # all splice rates equal
    kinetic = vcat(L.G, L.initiate, isempty(L.rshift) ? Int[] : [L.rshift[1]], L.eject,
                   isempty(L.splice) ? Int[] : [L.splice[1]])
    individual = L.n .+ NOISE_INDIVIDUAL
    shared_noise = L.n .+ NOISE_SHARED
    fitted = sort(vcat(kinetic, shared_noise, individual))
    return fitted, Tuple(ties), individual
end

function resolve_data(path)
    for p in (path, joinpath(ROOT, path), joinpath(ROOT, "data", path))
        isdir(p) && return p
    end
    error("trace folder not found: $path (tried as-is, under $ROOT and $ROOT/data)")
end

const NTRACES = Dict(s => length(StochasticGene.read_tracefiles(resolve_data(datapath(s)), CONDITION, 1, -1, 3))
                     for s in SAMPLES)
for (s, n) in NTRACES
    n > 0 || error("no '*$CONDITION*' trace files found for $s in $(datapath(s))")
end

const TRACE_SPECS = [(unit=1, interval=INTERVAL, start=1.0, t_end=-1.0, zeromedian=true)]

function run_row(m, sample)
    pm, pcv = priors(m, sample)
    fitted, ties, individual = fitted_and_ties(m)
    return (
        root=ROOT,
        key="$(sample)-tr-$(m.name)-$(CONDITION)",
        datatype="trace",
        datapath=datapath(sample),
        datacond=CONDITION,
        datacol=3,
        zeromedian=true,
        trace_specs=TRACE_SPECS,
        gene="MYC", cell="HBEC", nalleles=1, ejectnumber=1,
        probfn=prob_Gaussian,
        G=m.G, R=m.R, S=m.S, insertstep=m.insertstep, transitions=m.transitions,
        fittedparam=fitted,
        fixedeffects=ties,
        hierarchical=(2, individual, tuple()),
        priormean=pm, priorcv=pcv,
        noisepriors=NOISE_PRIORS[(sample, CONDITION)],
        decayrate=-1.0,
        elongationtime=m.R * RSHIFT_MINUTES,
        method=(Tsit5(), true),
        resultfolder=RESULTFOLDER,
        ratetype="median",
        nchains=NCHAINS, samplesteps=SAMPLESTEPS, sample_stride=STRIDE, warmupsteps=0,
        maxtime=MAXTIME, temp=1.0, burst=false, writesamples=WRITESAMPLES,
        propcv=PROPCV_RATE, propcv_rate=PROPCV_RATE, propcv_noise=PROPCV_NOISE,
    )
end

foreach(check_splice_live, MODELS)
const RUNS = [run_row(m, s) for s in SAMPLES for m in MODELS]

println("condition         : ", CONDITION, " (each sample fitted separately, nothing shared)")
println("traces            : ", join(["$s=$(NTRACES[s]) ($(datapath(s)))" for s in SAMPLES], ", "))
println("noise per trace   : offsets ", NOISE_INDIVIDUAL, "   per sample: offsets ", NOISE_SHARED)
println("priors            : G $PRIOR_G (CV $PRIORCV_G), Initiate $PRIOR_INITIATE (CV $PRIORCV_INITIATE), ",
        "Rshift 1/$RSHIFT_MINUTES min (CV $PRIORCV_R, all tied), Eject $EJECT_RATE (CV $PRIORCV_EJECT), ",
        "splice $SPLICE_PRIOR (CV $PRIORCV_SPLICE, all tied)")
println("run length        : maxtime $MAXTIME s, samplesteps $SAMPLESTEPS (not limiting), stride $STRIDE")
for m in MODELS
    fitted, ties, individual = fitted_and_ties(m)
    L = layout(m)
    println("  ", rpad(m.name, 9), " nrates=", L.n, " fitted=", fitted, " ties=", ties, " per-trace=", individual)
end

isdefined(StochasticGene, :stage_write_run_specs) || error("needs StochasticGene.stage_write_run_specs")
out = StochasticGene.stage_write_run_specs(
    Any[r for r in RUNS];
    resultfolder=RESULTFOLDER, root=ROOT, filedir=SWARMDIR,
    juliafile="fitscript", swarmfile="fit", commandfile="fit.commands",
    nchains=NCHAINS, nthreads=1, project=PROJECT,
)

specdir = joinpath(ROOT, "results", RESULTFOLDER)
for r in RUNS
    spec = StochasticGene.read_run_spec(joinpath(specdir, "info_$(r.key).toml"))
    for f in (:fittedparam, :fixedeffects, :hierarchical, :priormean, :priorcv, :G, :R, :S, :insertstep,
              :transitions, :noisepriors, :datapath, :maxtime, :samplesteps)
        haskey(spec, f) || error("$(r.key): $f absent after staging")
        spec[f] == getproperty(r, f) || error("$(r.key): $f changed during staging")
    end
end
println("spec check        : all ", length(RUNS), " specs round-tripped ok")

for f in ("fit.swarm", "fit.commands")
    p = joinpath(SWARMDIR, f); isfile(p) || continue
    lines = [startswith(l, "julia ") ? "JULIA_WORKER_TIMEOUT=600 JULIA_PKG_PRECOMPILE_AUTO=0 " * l : l for l in readlines(p)]
    write(p, join(lines, "\n") * "\n")
end
mkpath(joinpath(SWARMDIR, "logs"))
write(joinpath(SWARMDIR, "expected_commands.txt"), string(length(RUNS)), "\n")
wall = ceil(Int, MAXTIME / 3600) + 3
println("keys              : ", join(out.keys, ", "))
println("results           : ", specdir)
println("submit            : cd $SWARMDIR && swarm -f fit.swarm -g 32 -t $(NCHAINS + 2) --time $(lpad(wall, 2, '0')):00:00 --module julialang --logdir logs")
println("next              : SG_FROM=$RESULTFOLDER julia --project=\$SG_PROJECT scripts/13_continue_trace_fixnoise.jl")
