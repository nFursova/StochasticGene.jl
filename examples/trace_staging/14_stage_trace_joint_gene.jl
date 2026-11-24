# Joint 3Prime + 5Prime MYC (gene) trace fits: ALL kinetic rates shared (:top).
#
# Shared-parameter path (as in 05_swarm_joint_condition_expanded_topology.jl):
#   hierarchy   top -> group(:sample = 3Prime, 5Prime) -> individual(:trace_id)
#   kinetics    G, Initiate, Rshift/Eject, Decay       :top         one value for both samples
#   noise       offsets 1, 3 (background / reporter mean) :individual  per trace
#               offsets 2, 4 (background / reporter SD)   :group       per sample
#   ties        all Rshift rates equal; Eject separate
#   decay       not fitted
#
# Priors and models as in 12_stage_trace_separate_noise_gene.jl:
#   G 0.01 (CV 1), Initiate 0.2 (CV 0.25), Rshift 1/SG_RSHIFT_MINUTES (CV 0.25),
#   Eject SG_EJECT_RATE (CV 0.25), noise = per-sample stage-05 MYC priors (CV 0.1)
#   bases 3301 3302 3303 3401 3402 3403 3404 x topologies none cyc kp = 21 runs
#
# Run:
#   export SG_STAMP=$(date +%Y-%m-%d)
#   julia --project="$SG_PROJECT" scripts/14_stage_trace_joint_gene.jl
#
# Environment (defaults): SG_ROOT, SG_PROJECT=/data/$USER/sg201,
#   SG_RESULTFOLDER=3and5Prime-joint-gene-$SG_STAMP, SG_SWARMDIR=swarm-<resultfolder>,
#   SG_DATA_3PRIME=traces/3Prime, SG_DATA_5PRIME=traces/5Prime,
#   SG_BASES, SG_TOPOLOGIES=none,cyc,kp, SG_MODELS (full names; overrides),
#   SG_MAXTIME=43200 SG_SAMPLESTEPS=10^12 SG_SAMPLE_STRIDE=1000 SG_NCHAINS=16 SG_WRITESAMPLES=true
#   SG_PROPCV_TOP=0.05 SG_PROPCV_NOISE=0.001, SG_RSHIFT_MINUTES=2.5 SG_PRIORCV_R=0.25
#   SG_EJECT_RATE=1.0 SG_PRIORCV_EJECT=0.25 SG_PRIOR_G=0.01 SG_PRIORCV_G=1.0
#   SG_PRIOR_INITIATE=0.2 SG_PRIORCV_INITIATE=0.25

using StochasticGene
using Dates

envstr(v, d) = get(ENV, v, d)
envint(v, d) = parse(Int, get(ENV, v, string(d)))
envfloat(v, d) = parse(Float64, get(ENV, v, string(d)))
envbool(v, d) = lowercase(get(ENV, v, string(d))) in ("1", "true", "yes", "on")
csvlist(v, d) = [String(strip(s)) for s in split(get(ENV, v, d), ',') if !isempty(strip(s))]

const CONDITION = "gene"
const ROOT = abspath(envstr("SG_ROOT", "/data/Larsonlab/Nadya/StochasticGene"))
const PROJECT = envstr("SG_PROJECT", "/data/$(get(ENV, "USER", "unknown"))/sg201")
const STAMP = envstr("SG_STAMP", Dates.format(Dates.today(), "yyyy-mm-dd"))
const RESULTFOLDER = envstr("SG_RESULTFOLDER", "3and5Prime-joint-gene-$STAMP")
const SWARMDIR = joinpath(ROOT, envstr("SG_SWARMDIR", "swarm-" * RESULTFOLDER))

const MAXTIME = envfloat("SG_MAXTIME", 43200.0)
const SAMPLESTEPS = envint("SG_SAMPLESTEPS", 10^12)
const STRIDE = envint("SG_SAMPLE_STRIDE", 1000)
const NCHAINS = envint("SG_NCHAINS", 16)
const WRITESAMPLES = envbool("SG_WRITESAMPLES", true)
const PROPCV_TOP = envfloat("SG_PROPCV_TOP", 0.05)
const PROPCV_NOISE = envfloat("SG_PROPCV_NOISE", 0.001)
const INTERVAL = 1.6666666666666667

const PRIOR_G = envfloat("SG_PRIOR_G", 0.01)
const PRIORCV_G = envfloat("SG_PRIORCV_G", 1.0)
const PRIOR_INITIATE = envfloat("SG_PRIOR_INITIATE", 0.2)
const PRIORCV_INITIATE = envfloat("SG_PRIORCV_INITIATE", 0.25)
const RSHIFT_MINUTES = envfloat("SG_RSHIFT_MINUTES", 2.5)
const PRIORCV_R = envfloat("SG_PRIORCV_R", 0.25)
const EJECT_RATE = envfloat("SG_EJECT_RATE", 1.0)
const PRIORCV_EJECT = envfloat("SG_PRIORCV_EJECT", 0.25)
const DECAY_PLACEHOLDER = 0.03165055618995184

const SAMPLES = ["3Prime", "5Prime"]
datapath(sample) = envstr("SG_DATA_" * uppercase(sample), "traces/$sample")
# Per-sample MYC noise priors (stage 05). The shared path has one base prior block, seeded from
# 3Prime; each sample starts from its own values through initial_rates.
const NOISE_PRIORS = Dict("3Prime" => [-0.018007, 0.214861, 1.256466, 0.405046],
                          "5Prime" => [-0.004966, 0.118037, 0.607264, 0.223634])
const PRIOR_SAMPLE = "3Prime"

const TRANSITIONS = Dict(
    "none" => ([1, 2], [2, 1], [2, 3], [3, 2]),
    "cyc" => ([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1]),
    "kp" => ([1, 2], [2, 1], [2, 3], [3, 1]),
)
function model_spec(name)
    m = match(r"^(\d)(\d)(\d)(\d)([a-z]+)$", name)
    m === nothing && error("model name must be GRSI<topology>, got $name")
    G, R, S, I = parse.(Int, m.captures[1:4])
    (G == 3 && S == 0) || error("$name: this script stages 3-state models without splicing")
    haskey(TRANSITIONS, m.captures[5]) || error("$name: unknown topology $(m.captures[5])")
    1 <= I <= R || error("$name: insertstep must be in 1:R")
    (name=name, G=G, R=R, S=S, insertstep=I, transitions=TRANSITIONS[m.captures[5]])
end
const BASES = csvlist("SG_BASES", "3301,3302,3303,3401,3402,3403,3404")
const TOPOLOGIES = csvlist("SG_TOPOLOGIES", "none,cyc,kp")
const MODELS = [model_spec(n) for n in (haskey(ENV, "SG_MODELS") ? csvlist("SG_MODELS", "") :
                                         [b * t for b in BASES for t in TOPOLOGIES])]
isempty(MODELS) && error("no models selected")

nrates(m) = num_rates(m.transitions, m.R, m.S, m.insertstep)

function rate_families(m)
    nt = length(m.transitions); n = nrates(m)
    Dict{Any,Vector{Int}}(
        :G => collect(1:nt),
        :initiation => [nt + 1],
        :R => collect(nt+2:nt+1+m.R),              # Rshift_1..Rshift_{R-1}, Eject
        :decay => [n],
        :noise_mean => [n + 1, n + 3],
        :noise_sd => [n + 2, n + 4],
    )
end
parameter_scope(m) = Dict{Any,Symbol}(:G => :top, :initiation => :top, :R => :top, :decay => :top,
                                      :noise_mean => :individual, :noise_sd => :group)

const TRACE_SPECS = [(unit=1, interval=INTERVAL, start=1.0, t_end=-1.0, zeromedian=true)]

function priors(m)
    nt = length(m.transitions); n = nrates(m)
    mean = vcat(fill(PRIOR_G, nt), PRIOR_INITIATE, fill(1 / RSHIFT_MINUTES, m.R - 1), EJECT_RATE,
                DECAY_PLACEHOLDER, NOISE_PRIORS[PRIOR_SAMPLE])
    cv = vcat(fill(PRIORCV_G, nt), PRIORCV_INITIATE, fill(PRIORCV_R, m.R - 1), PRIORCV_EJECT, 0.1, fill(0.1, 4))
    return mean, cv
end
function initial_rates(m)
    n = nrates(m); base = priors(m)[1]
    d = Dict{Any,Vector{Float64}}()
    for s in SAMPLES
        v = copy(base); v[n+1:n+4] = NOISE_PRIORS[s]
        d[Symbol(s)] = v; d[Symbol("$(s)_$(CONDITION)")] = v
    end
    d
end
function fitted_and_ties(m)
    nt = length(m.transitions); n = nrates(m)
    rshift = collect(nt+2:nt+m.R)                  # Rshift_1..Rshift_{R-1}
    kinetic = vcat(1:nt, nt + 1, rshift[1], nt + 1 + m.R)
    ties = length(rshift) > 1 ? (rshift,) : tuple()
    return sort(vcat(kinetic, n+1:n+4)), ties
end

hierarchy(m) = HierarchySpec(
    HierarchyNode(:top, children=[
        HierarchyNode(:group, key=:sample, levels=Symbol.(SAMPLES),
            children=[HierarchyNode(:individual, key=:trace_id)])]);
    parameter_scope=parameter_scope(m))

datasets() = [DatasetSpec(Symbol("$(s)_$(CONDITION)"), "trace", datapath(s), CONDITION;
                          metadata=(sample=Symbol(s), dataset=Symbol("$(s)_$(CONDITION)")),
                          trace_specs=TRACE_SPECS) for s in SAMPLES]

shared_parameters(m) = (kind=:shared, levels=hierarchy(m), datasets=datasets(),
                        transition_families=[:G, :initiation, :R, :decay],
                        emission_families=[:noise_mean, :noise_sd],
                        rate_families=rate_families(m), initial_rates=initial_rates(m))

const PROPCV_LEVELS = Dict{Any,Float64}(:top => PROPCV_TOP, :group => PROPCV_NOISE, :individual => PROPCV_NOISE)

function resolve_data(path)
    for p in (path, joinpath(ROOT, path), joinpath(ROOT, "data", path))
        isdir(p) && return p
    end
    error("trace folder not found: $path")
end
const NTRACES = Dict(s => length(StochasticGene.read_tracefiles(resolve_data(datapath(s)), CONDITION, 1, -1, 3)) for s in SAMPLES)
all(>(0), values(NTRACES)) || error("no '*$CONDITION*' traces found: $NTRACES")

function run_row(m)
    pm, pcv = priors(m)
    fitted, ties = fitted_and_ties(m)
    (
        root=ROOT, key="3Prime-5Prime-trs-$(m.name)-$(CONDITION)",
        datatype="trace", datapath="", datacond=CONDITION, datacol=3, zeromedian=true,
        gene="MYC", cell="HBEC", nalleles=1, ejectnumber=1, probfn=prob_Gaussian,
        G=m.G, R=m.R, S=m.S, insertstep=m.insertstep, transitions=m.transitions,
        fittedparam=fitted, fixedeffects=ties, hierarchical=tuple(),
        shared_parameters=shared_parameters(m),
        priormean=pm, priorcv=pcv, noisepriors=NOISE_PRIORS[PRIOR_SAMPLE],
        decayrate=-1.0, elongationtime=m.R * RSHIFT_MINUTES, trace_specs=TRACE_SPECS,
        method=Tsit5(), resultfolder=RESULTFOLDER, ratetype="median",
        nchains=NCHAINS, samplesteps=SAMPLESTEPS, sample_stride=STRIDE, warmupsteps=0,
        maxtime=MAXTIME, temp=1.0, burst=false, writesamples=WRITESAMPLES,
        propcv=0.0, propcv_rate=PROPCV_TOP, propcv_noise=PROPCV_NOISE,
        propcv_levels=PROPCV_LEVELS, proposal_cv_levels=PROPCV_LEVELS,
    )
end

const RUNS = [run_row(m) for m in MODELS]
println("condition   : gene, 3Prime + 5Prime JOINT, all kinetic rates shared (:top)")
println("traces      : ", join(["$s=$(NTRACES[s])" for s in SAMPLES], ", "))
println("noise       : means per trace (:individual), SDs per sample (:group)")
println("run length  : maxtime $MAXTIME s, samplesteps $SAMPLESTEPS (not limiting), stride $STRIDE")
for (m, r) in zip(MODELS, RUNS)
    println("  ", rpad(m.name, 9), " nrates=", nrates(m), " fitted=", r.fittedparam, " ties=", r.fixedeffects)
end

out = StochasticGene.stage_write_run_specs(Any[r for r in RUNS];
    resultfolder=RESULTFOLDER, root=ROOT, filedir=SWARMDIR, juliafile="fitscript",
    swarmfile="fit", commandfile="fit.commands", nchains=NCHAINS, nthreads=1, project=PROJECT)

specdir = joinpath(ROOT, "results", RESULTFOLDER)
for (m, r) in zip(MODELS, RUNS)
    spec = StochasticGene.read_run_spec(joinpath(specdir, "info_$(r.key).toml"))
    for f in (:fittedparam, :fixedeffects, :priormean, :priorcv, :G, :R, :S, :insertstep, :transitions, :maxtime)
        haskey(spec, f) && spec[f] == getproperty(r, f) || error("$(r.key): $f changed during staging")
    end
    sp = spec[:shared_parameters]
    sp.rate_families == rate_families(m) || error("$(r.key): rate_families changed")
    sp.levels.parameter_scope == parameter_scope(m) || error("$(r.key): parameter_scope changed")
end
println("spec check  : all ", length(RUNS), " specs round-tripped ok")

for f in ("fit.swarm", "fit.commands")
    p = joinpath(SWARMDIR, f); isfile(p) || continue
    lines = [startswith(l, "julia ") ? "JULIA_WORKER_TIMEOUT=600 JULIA_PKG_PRECOMPILE_AUTO=0 " * l : l for l in readlines(p)]
    write(p, join(lines, "\n") * "\n")
end
mkpath(joinpath(SWARMDIR, "logs"))
write(joinpath(SWARMDIR, "expected_commands.txt"), string(length(RUNS)), "\n")
wall = ceil(Int, MAXTIME / 3600) + 3
println("results     : ", specdir)
println("submit      : cd $SWARMDIR && swarm -f fit.swarm -g 48 -t $(NCHAINS + 2) --time $(lpad(wall, 2, '0')):00:00 --module julialang --logdir logs")
