# Shared code for the 12-hour ON/OFF dwell-time staging scripts (08-11).
#
# Included by:
#   08_stage_dwelltime_types1.jl          MYC + CCAT1, model types 1
#   09_stage_dwelltime_types2_ccat1.jl    CCAT1,       model types 2 (splicing)
#   10_stage_crispr_types1.jl             CRISPR conditions, MYC + CCAT1, model types 1
#   11_stage_crispr_types2_ccat1.jl       CRISPR conditions, CCAT1,       model types 2
#
# Model names are <G><R><S><I><topology>, the numbers being exactly the G, R, S, insertstep
# passed to StochasticGene:
#   none  linear                      2-state: 1<->2          3-state: 1<->2<->3
#   cyc   3-state, fully reversible   1<->2, 2<->3, 1<->3
#   kp    3-state, kinetic proofreading 1<->2, 2->3, 3->1
#   gior  double-regime ("Giorgetti") model, staged as a COUPLED model:
#           unit 1 (hidden, G=2, R=0): regime low <-> high (k_f, k_b)
#           unit 2 (observed):          promoter OFF <-> ON, ON initiates, R steps as in the name
#           coupling (1,2,2,1) :activate  k_on(high) = k_on(low) * (1 + gamma)
#         The G digit of a gior name is the promoter's G (=2). One gior model is staged for each
#         R/S/I combination of the 3-state list.
#
# Splicing (model types 2) is staged only for insertstep = 1:
#   S = 1 -> one splice site at R1 ("simple");  S = R -> splicing at every step R1..R ("recursive")
# Splicing with insertstep > 1 is not staged: fit() resets S to at most R - insertstep + 1
# (reset_S), while the matrix builder only splices at steps j <= S (and at the last step only if
# S == R). With insertstep = R no S value reaches the last step, so the splice rate would never
# enter the likelihood (checked: the loglikelihood of 2333none is identical for splice rates
# 0.001, 1 and 10).
#
# Dwell types: dttype ON/OFF with onstates = [] => ON = at least one unspliced reporter-bearing
# R step occupied (RNA visible), OFF = none. For gior the ON/OFF times are those of unit 2.
#
# Run length: maxtime is the stop condition. MH grows its sample storage as it goes (no
# preallocation by samplesteps), so samplesteps defaults to 10^12 and is never reached.

using StochasticGene
using Dates

# ----------------------------------------------------------------------------------------
# environment helpers
# ----------------------------------------------------------------------------------------

envstr(var, default) = get(ENV, var, default)
envint(var, default) = parse(Int, get(ENV, var, string(default)))
envfloat(var, default) = parse(Float64, get(ENV, var, string(default)))
envbool(var, default) = lowercase(get(ENV, var, string(default))) in ("1", "true", "yes", "on")
csvlist(var, default) = [String(strip(s)) for s in split(get(ENV, var, default), ',') if !isempty(strip(s))]

const USER = get(ENV, "USER", "unknown")
const ROOT = abspath(envstr("SG_ROOT", "/data/Larsonlab/Nadya/StochasticGene"))
const PROJECT = envstr("SG_PROJECT", "/data/$USER/sg201")

# Run length: stop on time, not on steps.
const MAXTIME = envfloat("SG_MAXTIME", 43200.0)            # 12 h of sampling
const SAMPLESTEPS = envint("SG_SAMPLESTEPS", 10^12)        # effectively unlimited
const STRIDE = envint("SG_SAMPLE_STRIDE", 1000)
const NCHAINS = envint("SG_NCHAINS", 16)
const PROPCV = envfloat("SG_PROPCV", 0.02)
const WRITESAMPLES = envbool("SG_WRITESAMPLES", true)
# Extra wall time on top of maxtime for startup, compilation, merging and writing.
const WALL_PAD_HOURS = envint("SG_WALL_PAD_HOURS", 2)
const MEM_GB = envint("SG_MEM_GB", 32)
const TAG = envstr("SG_TAG", "$(Dates.today())-12h")
const DWELL_ROOT = envstr("SG_DWELL_ROOT", "dwelltimes")

const DTTYPE = csvlist("SG_DTTYPE", "ON,OFF")
(!isempty(DTTYPE) && all(in(("ON", "OFF")), DTTYPE)) ||
    error("SG_DTTYPE must be a nonempty subset of ON,OFF, got $(join(DTTYPE, ','))")

# Priors (same values as the earlier dwell-time staging unless overridden)
const ELONGATION = Dict(
    "gene" => envfloat("SG_ELONGATION_GENE", 3.2),          # MYC, total minutes
    "enhancer" => envfloat("SG_ELONGATION_ENHANCER", 22.0),  # CCAT1, total minutes
)
const PRIOR_G = envfloat("SG_PRIOR_G", 0.01)
const PRIORCV_G = envfloat("SG_PRIORCV_G", 1.0)
const PRIOR_INITIATE = envfloat("SG_PRIOR_INITIATE", 0.2)
const PRIORCV_INITIATE = envfloat("SG_PRIORCV_INITIATE", 0.25)
const PRIORCV_R = envfloat("SG_PRIORCV_R", 0.25)
const SPLICE_PRIOR = envfloat("SG_SPLICE_PRIOR", 0.05)
const PRIORCV_SPLICE = envfloat("SG_PRIORCV_SPLICE", 0.2)
const GIOR_GAMMA = envfloat("SG_GIOR_GAMMA", 0.1)          # coupling gamma: k_on(high) = k_on(low)*(1+gamma)
const GIOR_GAMMA_CV = envfloat("SG_GIOR_GAMMA_CV", 10.0)
const DECAY_PLACEHOLDER = 0.03165055618995184
const HIDDEN_EJECT_PLACEHOLDER = 1.0                       # unused rate of the R=0 regime unit

# ----------------------------------------------------------------------------------------
# model lists
# ----------------------------------------------------------------------------------------

const TRANSITIONS = Dict(
    (2, "none") => ([1, 2], [2, 1]),
    (3, "none") => ([1, 2], [2, 1], [2, 3], [3, 2]),
    (3, "cyc") => ([1, 2], [2, 1], [2, 3], [3, 2], [1, 3], [3, 1]),
    (3, "kp") => ([1, 2], [2, 1], [2, 3], [3, 1]),
)
const GIOR_REGIME_TRANSITIONS = ([1, 2], [2, 1])   # k_f (low->high), k_b (high->low)
const GIOR_PROMOTER_TRANSITIONS = ([1, 2], [2, 1]) # k_on (OFF->ON), k_off (ON->OFF)

# (R, insertstep) combinations
const RI_TYPES1 = [(3, 1), (3, 2), (3, 3), (4, 1), (4, 3), (4, 4), (5, 1), (5, 4), (5, 5)]
const RI_TYPES2 = [(3, 1), (4, 1), (5, 1)]   # insertstep 1 only (see the splicing note above)

splice_options(R, I) = I == 1 ? [1, R] : error("splicing is staged only for insertstep 1 (got R=$R, insertstep=$I)")

struct ModelSpec
    name::String
    topology::String     # none | cyc | kp | gior
    G::Int               # promoter G (gior: 2)
    R::Int
    S::Int
    insertstep::Int
end

modelname(G, R, S, I, topo) = "$(G)$(R)$(S)$(I)$(topo)"

function models_for(RSI; topologies)
    out = ModelSpec[]
    for (R, S, I) in RSI
        "none" in topologies && push!(out, ModelSpec(modelname(2, R, S, I, "none"), "none", 2, R, S, I))
        for topo in ("none", "cyc", "kp")
            topo in topologies && push!(out, ModelSpec(modelname(3, R, S, I, topo), topo, 3, R, S, I))
        end
        "gior" in topologies && push!(out, ModelSpec(modelname(2, R, S, I, "gior"), "gior", 2, R, S, I))
    end
    return out
end

function model_list(types::Int)
    topologies = csvlist("SG_TOPOLOGIES", "none,cyc,kp,gior")
    all(in(("none", "cyc", "kp", "gior")), topologies) || error("SG_TOPOLOGIES entries must be none,cyc,kp,gior")
    RSI = types == 1 ? [(R, 0, I) for (R, I) in RI_TYPES1] :
          [(R, S, I) for (R, I) in RI_TYPES2 for S in splice_options(R, I)]
    models = models_for(RSI; topologies=topologies)
    wanted = csvlist("SG_MODELS", "")
    if !isempty(wanted)
        names = Set(m.name for m in models)
        bad = setdiff(wanted, names)
        isempty(bad) || error("SG_MODELS has names not in model types $types: $(join(bad, ", "))")
        models = [m for m in models if m.name in wanted]
    end
    isempty(models) && error("no models selected")
    return models
end

# ----------------------------------------------------------------------------------------
# rates, priors, fitted parameters
# ----------------------------------------------------------------------------------------

nsplice(R, S, I) = S == 0 ? 0 : min(S, R - I + 1)

# promoter block: G transitions | Initiate | R steps (Rshift_1..Rshift_{R-1}, Eject) | splice | Decay
function promoter_priors(transitions, R, S, I, elongation)
    mean = vcat(fill(PRIOR_G, length(transitions)), PRIOR_INITIATE,
                fill(R / elongation, R), fill(SPLICE_PRIOR, nsplice(R, S, I)), DECAY_PLACEHOLDER)
    cv = vcat(fill(PRIORCV_G, length(transitions)), PRIORCV_INITIATE,
              fill(PRIORCV_R, R - 1), max(PRIORCV_R, 0.25), fill(PRIORCV_SPLICE, nsplice(R, S, I)), 0.1)
    return mean, cv
end

function single_unit_fields(m::ModelSpec, datacond)
    tr = TRANSITIONS[(m.G, m.topology)]
    n = num_rates(tr, m.R, m.S, m.insertstep)
    pm, pcv = promoter_priors(tr, m.R, m.S, m.insertstep, ELONGATION[datacond])
    length(pm) == n || error("$(m.name): $(length(pm)) priors for $n rates")
    return (
        G=m.G, R=m.R, S=m.S, insertstep=m.insertstep,
        transitions=tr,
        onstates=[Int[] for _ in DTTYPE],
        fittedparam=collect(1:n-1),           # Decay never enters the dwell-time T matrix
        priormean=pm, priorcv=pcv,
        elongationtime=ELONGATION[datacond],
    )
end

function gior_fields(m::ModelSpec, datacond)
    tr1, tr2 = GIOR_REGIME_TRANSITIONS, GIOR_PROMOTER_TRANSITIONS
    n1 = num_rates(tr1, 0, 0, 1)                        # k_f, k_b, Eject, Decay (last two unused)
    n2 = num_rates(tr2, m.R, m.S, m.insertstep)
    pm2, pcv2 = promoter_priors(tr2, m.R, m.S, m.insertstep, ELONGATION[datacond])
    pm = vcat(PRIOR_G, PRIOR_G, HIDDEN_EJECT_PLACEHOLDER, DECAY_PLACEHOLDER, pm2, GIOR_GAMMA)
    pcv = vcat(PRIORCV_G, PRIORCV_G, 0.1, 0.1, pcv2, GIOR_GAMMA_CV)
    length(pm) == n1 + n2 + 1 || error("$(m.name): prior length mismatch")
    fitted = vcat(1, 2, collect(n1+1:n1+n2-1), n1 + n2 + 1)
    return (
        G=(2, m.G), R=(0, m.R), S=(0, m.S), insertstep=(1, m.insertstep),
        transitions=(tr1, tr2),
        coupling=((1, 2), [(1, 2, 2, 1)], [:activate]),
        onstates=[[[2] for _ in DTTYPE], [Int[] for _ in DTTYPE]],
        dwell_specs=[(unit=2, onstates=[Int[] for _ in DTTYPE], dttype=copy(DTTYPE))],
        noisepriors=[Float64[], Float64[]],     # coupled models need one (empty) entry per unit
        fittedparam=fitted,
        priormean=pm, priorcv=pcv,
        elongationtime=ELONGATION[datacond],
    )
end

model_fields(m::ModelSpec, datacond) = m.topology == "gior" ? gior_fields(m, datacond) : single_unit_fields(m, datacond)

function run_row(m::ModelSpec; key, datapath, datacond, resultfolder)
    common = (
        root=ROOT,
        key=key,
        datatype="dwelltime",
        datapath=datapath,
        datacond=datacond,
        dttype=copy(DTTYPE),
        gene="MYC",
        cell="HBEC",
        nalleles=1,
        ejectnumber=1,
        fixedeffects=tuple(),
        hierarchical=tuple(),
        noisepriors=Float64[],                  # gior overrides with one entry per unit
        decayrate=-1.0,
        method=Tsit5(),
        resultfolder=resultfolder,
        ratetype="ml",
        nchains=NCHAINS,
        samplesteps=SAMPLESTEPS,
        sample_stride=STRIDE,
        warmupsteps=0,
        maxtime=MAXTIME,
        temp=1.0,
        burst=false,
        writesamples=WRITESAMPLES,
        propcv=PROPCV,
    )
    return merge(common, model_fields(m, datacond))
end

# ----------------------------------------------------------------------------------------
# data checks
# ----------------------------------------------------------------------------------------

# Mirror how fit() resolves datapath entries: as-is, then under ROOT, then under ROOT/data.
function resolve_data(path)
    ispath(path) && return path
    p = joinpath(ROOT, path)
    ispath(p) && return p
    p = joinpath(ROOT, "data", path)
    ispath(p) && return p
    return nothing
end

# Bins must be uniform and run d, 2d, 3d, ... (dwelltimeCDF prepends 2*bins[1]-bins[2] = 0).
function check_histogram(file)
    c = StochasticGene.readfile(file)
    size(c, 2) >= 2 || error("$file must have two columns (bin, count)")
    bins = Float64.(c[:, 1])
    vals = Float64.(c[:, 2])
    length(bins) >= 2 || error("$file needs at least two bins")
    d = diff(bins)
    all(abs.(d .- d[1]) .<= 1e-6 * max(d[1], 1.0)) || error("$file bins are not uniform")
    abs(2 * bins[1] - bins[2]) <= 1e-6 * max(bins[2], 1.0) ||
        error("$file first bin is $(bins[1]) but spacing is $(d[1]); bins must run d, 2d, 3d, ...")
    all(vals .>= 0) || error("$file has negative values")
    total = sum(vals)
    total > 0 || error("$file sums to zero")
    isapprox(total, 1.0; atol=1e-6) &&
        @warn "$file sums to 1.0 - it looks normalized; write raw event counts so sample size weights the fit."
    return (nbins=length(bins), width=d[1], tmax=bins[end], total=total)
end

# ----------------------------------------------------------------------------------------
# staging
# ----------------------------------------------------------------------------------------

describe(m::ModelSpec) = m.topology == "gior" ?
    "hidden regime (G2 R0) x promoter G2 R$(m.R) S$(m.S) I$(m.insertstep), k_on x (1+gamma) in high regime" :
    "G$(m.G) R$(m.R) S$(m.S) I$(m.insertstep) $(m.topology)"

function print_models(models)
    println("models            : ", length(models))
    for m in models
        println("  ", rpad(m.name, 12), describe(m))
    end
end

# Write info_<key>.jld2 + fit scripts + fit.swarm/fit.commands, check the round trip,
# make the swarm robust to worker start-up, and print the submit line.
function stage!(runs; swarmdir, default_resultfolder)
    isempty(runs) && error("nothing to stage")
    keys_ = [r.key for r in runs]
    length(unique(keys_)) == length(keys_) || error("run keys are not unique")
    for r in runs
        length(r.priormean) == length(r.priorcv) || error("$(r.key): priorcv length must match priormean")
        all(>(0), r.priormean) || error("$(r.key): priormean must be positive")
        length(r.dttype) == length(r.datapath) || error("$(r.key): one data file per dwell type")
        all(isempty, r.noisepriors) || error("$(r.key): dwell-time fits carry no noise parameters")
    end

    isdefined(StochasticGene, :stage_write_run_specs) ||
        error("this script requires StochasticGene.stage_write_run_specs")
    out = StochasticGene.stage_write_run_specs(
        Any[r for r in runs];
        resultfolder=default_resultfolder,
        root=ROOT,
        filedir=swarmdir,
        juliafile="fitscript",
        swarmfile="fit",
        commandfile="fit.commands",
        nchains=NCHAINS,
        nthreads=1,
        project=PROJECT,
    )

    # Round-trip check of what fit(; key=...) will read.
    for r in runs
        specfile = joinpath(ROOT, "results", r.resultfolder, "info_$(r.key).toml")
        spec = StochasticGene.read_run_spec(specfile)
        for field in (:fittedparam, :priormean, :priorcv, :G, :R, :S, :insertstep, :transitions,
                      :dttype, :onstates, :datapath, :datatype, :datacond, :maxtime, :samplesteps,
                      :sample_stride)
            haskey(spec, field) || error("$(r.key): $field absent after staging")
            spec[field] == getproperty(r, field) || error("$(r.key): $field changed during staging")
        end
        if haskey(r, :coupling)
            spec[:coupling] == r.coupling || error("$(r.key): coupling changed during staging")
            haskey(spec, :dwell_specs) && !isempty(spec[:dwell_specs]) ||
                error("$(r.key): dwell_specs absent after staging")
        end
    end
    println("spec check        : all ", length(runs), " specs round-tripped ok")

    # gior (coupled) runs: load the coupled dwell-time patch on every worker, so they also run
    # on an unpatched StochasticGene (e.g. 2.0.7). The patch skips itself if the fix is present.
    gior_keys = Set(r.key for r in runs if haskey(r, :coupling))
    if !isempty(gior_keys) && envbool("SG_GIOR_PATCH", true)
        patch_src = joinpath(@__DIR__, "stochasticgene_coupled_dwell_patch.jl")
        isfile(patch_src) || error("patch file not found: $patch_src")
        patch_dst = joinpath(abspath(swarmdir), basename(patch_src))
        cp(patch_src, patch_dst; force=true)
        npatched = 0
        for (k, script) in zip(out.keys, out.scripts)
            k in gior_keys || continue
            lines = readlines(script)
            i = findfirst(l -> occursin("using StochasticGene", l), lines)
            i === nothing && error("$script: no `using StochasticGene` line to patch after")
            insert!(lines, i + 1, "@everywhere include($(repr(patch_dst)))")
            write(script, join(lines, "\n") * "\n")
            npatched += 1
        end
        npatched == length(gior_keys) || error("patched $npatched of $(length(gior_keys)) gior fit scripts")
        println("gior patch        : ", npatched, " fit scripts load ", patch_dst)
    end

    # Same worker-startup guard as the earlier batches.
    if envbool("SG_WORKER_GUARD", true)
        for f in ("fit.swarm", "fit.commands")
            p = joinpath(swarmdir, f)
            isfile(p) || continue
            lines = readlines(p)
            lines = [startswith(l, "julia ") ? "JULIA_WORKER_TIMEOUT=600 JULIA_PKG_PRECOMPILE_AUTO=0 " * l : l for l in lines]
            write(p, join(lines, "\n") * "\n")
        end
    end
    mkpath(joinpath(swarmdir, "logs"))
    write(joinpath(swarmdir, "expected_commands.txt"), string(length(runs)), "\n")

    wall = ceil(Int, MAXTIME / 3600) + WALL_PAD_HOURS
    folders = sort(unique(r.resultfolder for r in runs))
    println("runs              : ", length(runs))
    println("result folders    : ", join(folders, ", "))
    println("swarm             : ", joinpath(swarmdir, "fit.swarm"))
    println("maxtime           : ", MAXTIME, " s (", round(MAXTIME / 3600, digits=1), " h); samplesteps ",
            SAMPLESTEPS, " (not limiting); stride ", STRIDE)
    println("submit            : cd $swarmdir && swarm -f fit.swarm -g $MEM_GB -t $(NCHAINS + 2) ",
            "--time $(lpad(wall, 2, '0')):00:00 --module julialang --logdir logs")
    return out
end

# ----------------------------------------------------------------------------------------
# data units -> runs
# ----------------------------------------------------------------------------------------

# A data unit is one set of ON/OFF histograms. `units` holds NamedTuples with
#   label, datacond, datapath (unresolved, as fit() will see it), resultfolder, keyfn(model)
# Units whose files are missing are skipped with a note.
function build_runs(units, models)
    runs = Any[]
    found = 0
    println("dwell types       : ", join(DTTYPE, ", "))
    println("elongation prior  : gene=", ELONGATION["gene"], " min, enhancer=", ELONGATION["enhancer"],
            " min; R-step priorcv=", PRIORCV_R)
    println("data units:")
    for u in units
        got = map(resolve_data, u.datapath)
        if any(isnothing, got)
            println("  skip (no data)  ", u.label)
            continue
        end
        found += 1
        for (dt, f) in zip(DTTYPE, got)
            s = check_histogram(f)
            println("  ", rpad("$(u.label) $(dt)", 40), " nbins=", lpad(s.nbins, 4),
                    " d=", round(s.width, digits=4), " tmax=", lpad(round(s.tmax, digits=1), 7),
                    " events=", lpad(round(Int, s.total), 6))
        end
        for m in models
            push!(runs, run_row(m; key=u.keyfn(m), datapath=u.datapath, datacond=u.datacond,
                                resultfolder=u.resultfolder))
        end
    end
    found == 0 && error("no dwell-time files found under $ROOT or $ROOT/data ($DWELL_ROOT)")
    println("units with data   : ", found, "; models per unit ", length(models))
    return runs
end
