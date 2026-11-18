# stage_coupled_hidden.jl
#
# Table-driven staging of coupled models with a HIDDEN unit that regulates both genes. Same seeding as
# stage_coupled_reciprocal.jl / stage_coupled_from_separate_final.jl: unit 1 = enhancer, unit 2 = gene, both
# seeded from two SEPARATE hierarchical fits, all four noise parameters frozen PER TRACE, full rinit built here.
# Unit 3 = the hidden unit: no reporter, no observations, no noise parameters (StochasticGene's hidden-latent
# layout, unit_model = (1, 2, 3)).
#
# =======================================================================================
# MODEL SPEC SYNTAX
#
#   [<name> = ] <kind>:<targets> [@ <mode>]
#
#   kind     H2    two hidden states 1 <-> 2; state 2 regulates the target rates of BOTH genes
#            H3s   three hidden states 1 <-> 2 <-> 3; state 1 regulates the enhancer, state 3 the gene;
#                  symmetric: h23 = h21, h32 = h12, and each target rate has the SAME gamma in both genes
#            H3a   the same, asymmetric: four free hidden rates and separate gammas for the two genes
#   targets  comma separated transition indices (or names) regulated in BOTH genes:
#                1 = Rate_12  2 = Rate_21  3 = Rate_23  4 = Rate_32  5 = Initiate
#   mode     free (default) | activate | inhibit, applied to every gamma
#
#   H3s-R12R23 = H3s:1,3@free
#
# GAMMA STRUCTURE. One gamma per (gene, target rate), free; H3s ties each rate's two gammas.
# HIDDEN RATES. Start at SG_HIDDEN_RATE0 (0.1 per minute) and are fitted (SG_FIT_HIDDEN=1, the default) with
# prior cv SG_PRIORCV_HIDDEN (10). Rate vector of the hidden unit: its switching rates (h12, h21 for H2;
# h12, h21, h23, h32 for H3), then one slot StochasticGene keeps for a unit without reporter steps (held fixed),
# then decay.
#
# =======================================================================================
# RUN (environment as for stage_coupled_reciprocal.jl, plus SG_MODELFILE or SG_MODELS)
#
#   julia --project=$SG_PROJECT stage_coupled_hidden.jl

using StochasticGene
using DelimitedFiles
using Statistics

println("StochasticGene: version=", pkgversion(StochasticGene),
        "  project=", Base.active_project(), "  Julia=", VERSION)

let missing = [f for f in (:read_run_spec, :stage_write_run_specs, :info_jld2_path)
               if !isdefined(StochasticGene, f)]
    isempty(missing) || error("""
        This StochasticGene build ($(pkgversion(StochasticGene))) lacks $(join(missing, ", ")).
        Launch with an explicit project: julia --project=/data/fursovan2/sg111 ...""")
end

# ------------------------------------------------------------------------------- CONFIG
const ROOT     = abspath(get(ENV, "SG_ROOT", "/data/Larsonlab/Nadya/StochasticGene"))
const RUN_ROOT = abspath(get(ENV, "SG_RUN_ROOT", ROOT))
const PROJECT  = get(ENV, "SG_PROJECT", "")
isempty(strip(PROJECT)) && error("SG_PROJECT must be set")

const SAMPLE     = strip(get(ENV, "SG_SAMPLE", "3Prime"))
const SRC_FOLDER = strip(get(ENV, "SG_SRC_FOLDER", ""))
const SRC_ROW    = strip(get(ENV, "SG_SRC_ROW", "median"))
const DATA       = strip(get(ENV, "SG_DATA", ""))
isempty(DATA) && error("SG_DATA must point at the folder holding the PAIRED gene/enhancer traces")
const INTO       = strip(get(ENV, "SG_INTO", "coupled-hidden-" * SAMPLE))
const SWARMROOT  = joinpath(ROOT, get(ENV, "SG_SWARMROOT", "swarm_" * INTO))

const G, R, S, INSERTSTEP = 3, 3, 0, 1
const TRANSITIONS = ([1, 2], [2, 1], [2, 3], [3, 2])
const NRATES  = num_rates(TRANSITIONS, R, S, INSERTSTEP)     # 9
const NNOISE  = 4
const NUNIT   = NRATES + NNOISE                              # 13: rates + noise of one observed unit
const NHYPER  = 2
# Hidden unit transitions: H2 1 <-> 2; H3 the chain 1 <-> 2 <-> 3 (h12, h21, h23, h32)
const HIDDEN_TRANSITIONS = Dict(:H2 => ([1, 2], [2, 1]), :H3 => ([1, 2], [2, 1], [2, 3], [3, 2]))
const HIDDEN_SOURCE = Dict(:H2 => (2, 2), :H3 => (1, 3))    # hidden state regulating (enhancer, gene)

const INTERVAL   = parse(Float64, get(ENV, "SG_INTERVAL", "1.6666666666666667"))
const ELONG      = (parse(Float64, get(ENV, "SG_ELONG_ENH",  "20.0")),
                    parse(Float64, get(ENV, "SG_ELONG_GENE", "5.0")))
const NCHAINS    = parse(Int, get(ENV, "SG_NCHAINS", "16"))
const MAXTIME    = parse(Float64, get(ENV, "SG_MAXTIME", "86400"))
const SAMPLES    = parse(Int, get(ENV, "SG_SAMPLES", "10000000"))
const STRIDE     = parse(Int, get(ENV, "SG_STRIDE", "100"))
const PROPCV     = parse(Float64, get(ENV, "SG_PROPCV", "0.1"))

const COUPLE_MODE = strip(get(ENV, "SG_COUPLE_MODE", "target"))
COUPLE_MODE in ("target", "rates") ||
    error("SG_COUPLE_MODE must be \"target\" or \"rates\", got $(repr(COUPLE_MODE))")

const PRIORCV_TARGET   = parse(Float64, get(ENV, "SG_PRIORCV_TARGET",   "1.0"))
const PRIORCV_G        = parse(Float64, get(ENV, "SG_PRIORCV_G",        "1.0"))
const PRIORCV_INITIATE = parse(Float64, get(ENV, "SG_PRIORCV_INITIATE", "0.25"))
const PRIORCV_COUPLING = parse(Float64, get(ENV, "SG_PRIORCV_COUPLING", "10.0"))
const PRIORCV_HIDDEN   = parse(Float64, get(ENV, "SG_PRIORCV_HIDDEN",   "10.0"))
const HIDDEN_RATE0     = parse(Float64, get(ENV, "SG_HIDDEN_RATE0",     "0.1"))
const FIT_HIDDEN       = get(ENV, "SG_FIT_HIDDEN", "1") == "1"

const GAMMA0 = Dict(:free     => parse(Float64, get(ENV, "SG_GAMMA0",      "0.0")),
                    :activate => parse(Float64, get(ENV, "SG_GAMMA0_ACT",  "0.1")),
                    :inhibit  => parse(Float64, get(ENV, "SG_GAMMA0_INH", "-0.5")))
GAMMA0[:activate] > 0 || error("SG_GAMMA0_ACT must be > 0 (log transform); got $(GAMMA0[:activate])")
-1 < GAMMA0[:inhibit] < 0 || error("SG_GAMMA0_INH must be in (-1, 0); got $(GAMMA0[:inhibit])")
GAMMA0[:free] > -1 || error("SG_GAMMA0 must be > -1; got $(GAMMA0[:free])")
# --------------------------------------------------------------------------- END CONFIG

const RATE_LABELS = ["Rate_12", "Rate_21", "Rate_23", "Rate_32", "Initiate"]
const TARGET_ALIAS = Dict("R12" => 1, "R21" => 2, "R23" => 3, "R32" => 4, "INIT" => 5, "INITIATE" => 5)
const MODE_SUFFIX = Dict(:free => "free", :activate => "act", :inhibit => "inh")
global_index(u, t) = (u - 1) * NUNIT + t                     # observed units 1 and 2

"""Parse `[name = ] <kind>:<targets> [@ mode]`."""
function parse_spec(entry)
    name, specstr = occursin('=', entry) ? String.(strip.(split(entry, '=', limit = 2))) : ("", String(strip(entry)))
    body, md = occursin('@', specstr) ? String.(strip.(rsplit(specstr, '@', limit = 2))) : (specstr, "free")
    mode = Symbol(lowercase(md))
    mode in (:free, :activate, :inhibit) || error("$(repr(entry)): mode must be free, activate or inhibit")
    parts = String.(strip.(split(body, ':')))
    length(parts) == 2 || error("$(repr(entry)): spec must be <kind>:<targets>")
    kind = uppercase(parts[1])
    kind in ("H2", "H3S", "H3A") || error("$(repr(entry)): kind must be H2, H3s or H3a")
    hidden = kind == "H2" ? :H2 : :H3
    symmetric = kind == "H3S"
    targets = map(split(parts[2], ',')) do t
        t = uppercase(strip(t))
        haskey(TARGET_ALIAS, t) ? TARGET_ALIAS[t] :
            (v = tryparse(Int, t); v === nothing || !(1 <= v <= 5) ?
             error("$(repr(entry)): bad target $(repr(t)); use 1..5 or R12/R21/R23/R32/Init") : v)
    end
    allunique(targets) || error("$(repr(entry)): repeated target")

    nh = length(HIDDEN_TRANSITIONS[hidden])                  # hidden switching rates
    hblock = nh + 2                                          # + 1 fixed slot + decay = num_rates(R = 0)
    offset = 2 * NUNIT + hblock                              # gammas start after the hidden unit
    sE, sG = HIDDEN_SOURCE[hidden]
    conns = NTuple{4,Int}[]
    for t in targets
        push!(conns, (3, sE, 1, t))                          # hidden -> enhancer rate t
        push!(conns, (3, sG, 2, t))                          # hidden -> gene rate t
    end
    gamma_index(k) = offset + k
    # tied groups (first index free): H3s ties the two genes' gammas of each rate, and h23 = h21, h32 = h12
    gamma_groups = symmetric ? [[gamma_index(2j - 1), gamma_index(2j)] for j in eachindex(targets)] :
                               [[gamma_index(k)] for k in eachindex(conns)]
    h(i) = 2 * NUNIT + i                                     # index of hidden rate i (h12 = 1, h21 = 2, h23 = 3, h32 = 4)
    hidden_groups = hidden == :H2 ? [[h(1)], [h(2)]] :
                    symmetric     ? [[h(1), h(4)], [h(2), h(3)]] : [[h(1)], [h(2)], [h(3)], [h(4)]]
    key = (isempty(name) ? string(kind == "H3S" ? "H3s" : kind == "H3A" ? "H3a" : "H2", "-",
                                  join(RATE_LABELS[t] == "Initiate" ? "Init" : replace(RATE_LABELS[t], "Rate_" => "R")
                                       for t in targets), "-", MODE_SUFFIX[mode]) : name) * "-h"
    occursin('_', key) && error("$key: keys must not contain '_' (it separates filename fields)")
    return (; key, spec = specstr, kind, hidden, symmetric, targets, conns, nh, hblock, gamma_groups, hidden_groups, mode)
end

const MODEL_SPECS = let
    txt = strip(get(ENV, "SG_MODELS", ""))
    file = strip(get(ENV, "SG_MODELFILE", ""))
    raw = !isempty(file) ? read(file, String) : replace(txt, ';' => '\n')
    [l for l in strip.(split(raw, '\n')) if !isempty(l) && !startswith(l, "#")]
end
isempty(MODEL_SPECS) && error("no models: set SG_MODELFILE (one spec per line) or SG_MODELS (separated by ;)")
const MODELS = [parse_spec(s) for s in MODEL_SPECS]
allunique([m.key for m in MODELS]) || error("model keys are not unique: $([m.key for m in MODELS])")

# --- read the two source fits (identical to stage_coupled_reciprocal.jl) ----------------
function src_dir(env)
    isempty(SRC_FOLDER) && error("""
        SG_SRC_FOLDER is not set, so there is no folder to resolve $env against.
        Either set SG_SRC_FOLDER or give SG_SRC_RATES_GENE and SG_SRC_RATES_ENH as absolute paths.""")
    dir = joinpath(ROOT, "results", SRC_FOLDER)
    isdir(dir) || error("source results folder not found: $dir")
    return dir
end

function src_ratefile(transcript)
    env = transcript == "gene" ? "SG_SRC_RATES_GENE" : "SG_SRC_RATES_ENH"
    explicit = strip(get(ENV, env, ""))
    if !isempty(explicit)
        p = isabspath(explicit) ? explicit : joinpath(src_dir(env), explicit)
        isfile(p) || error("$env points at a missing file: $p")
        return p
    end
    dir = src_dir(env)
    suffix = "-" * SAMPLE * ".txt"
    c = filter(f -> startswith(f, "rates_") && endswith(f, suffix) && occursin("-" * transcript * "-", f), readdir(dir))
    length(c) == 1 || error("expected exactly one rates_*-$transcript-*$suffix in $dir, found $(length(c)). Set $env.")
    return joinpath(dir, c[1])
end

function read_unit(transcript)
    path = src_ratefile(transcript)
    raw  = readdlm(path, ',')
    row  = 1 + StochasticGene.get_row(SRC_ROW)
    size(raw, 1) >= row || error("$path has $(size(raw,1)) lines; need labels + 4 rows")
    hdr  = String.(strip.(string.(raw[1, :])))
    vals = Float64.(raw[row, :])
    labels = ["Rate_12", "Rate_21", "Rate_23", "Rate_32", "Initiate", "Rshift_1", "Rshift_2", "Eject", "Decay"]
    rates = map(labels) do nm
        i = findfirst(==(nm), hdr)
        i === nothing && error("column $nm not in $path")
        vals[i]
    end
    found = Dict{Tuple{Int,Int},Float64}()
    rx = Regex("^noiseparam_(\\d+)_.*_trace(\\d+)\$")
    for (i, h) in enumerate(hdr)
        m = match(rx, h)
        m === nothing && continue
        found[(parse(Int, m[1]), parse(Int, m[2]))] = vals[i]
    end
    isempty(found) && error("no per-trace noiseparam_*_trace<N> columns in $path")
    ntr = maximum(k[2] for k in keys(found))
    noise = Matrix{Float64}(undef, NNOISE, ntr)
    for p in 1:NNOISE, t in 1:ntr
        haskey(found, (p, t)) || error("$path is missing noiseparam_$(p)_..._trace$t")
        noise[p, t] = found[(p, t)]
    end
    println("  $transcript  <- ", path, "\n                 ($ntr traces, row \"$SRC_ROW\")")
    return (; rates, noise, ntraces = ntr, path)
end

println("\nsource fits:")
const ENH  = read_unit("enhancer")
const GENE = read_unit("gene")
ENH.ntraces == GENE.ntraces || error("enhancer has $(ENH.ntraces) traces, gene has $(GENE.ntraces); cannot pair")
const NTR = ENH.ntraces
isapprox(ENH.rates[NRATES], GENE.rates[NRATES]; rtol = 1e-8) ||
    @warn("Decay differs between units ($(ENH.rates[NRATES]) vs $(GENE.rates[NRATES]))")
const DECAY = parse(Float64, get(ENV, "SG_DECAY", string(ENH.rates[NRATES])))

const DATADIR = joinpath(ROOT, "data", DATA)
isdir(DATADIR) || error("paired trace folder does not exist: $DATADIR")
const PAIRS = let
    files = readdir(DATADIR)
    by = [sort(filter(f -> occursin(l, f), files)) for l in ("enhancer", "gene")]
    any(isempty, by) && error("no files matching \"enhancer\" and/or \"gene\" in $DATADIR")
    (; n = minimum(length.(by)), enh = by[1], gene = by[2])
end
println("\npaired traces in $DATADIR : $(PAIRS.n)")
if PAIRS.n != NTR && get(ENV, "SG_ALLOW_TRACE_COUNT_MISMATCH", "") != "1"
    error("""The source fits have $NTR traces but this folder yields $(PAIRS.n) pairs. Per-trace noise is assigned
        BY POSITION. Point SG_DATA at the folder the separate fits used, or set SG_ALLOW_TRACE_COUNT_MISMATCH=1.""")
end

# --- rate vector -----------------------------------------------------------------------
hyper2_unit() = [1.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25, 1.0, 1.0, 0.25, 0.25, 0.25, 0.25]
hidden_block(m) = vcat(fill(HIDDEN_RATE0, m.nh), [1.0, DECAY])          # switching rates, 1 fixed slot, decay
nall(m) = 2 * NUNIT + m.hblock + length(m.conns)

function block(m, enh_noise::AbstractVector, gene_noise::AbstractVector, g0::Vector{Float64})
    r = Float64[]
    append!(r, ENH.rates);  r[NRATES] = DECAY;  append!(r, enh_noise)
    append!(r, GENE.rates); r[NUNIT + NRATES] = DECAY; append!(r, gene_noise)
    append!(r, hidden_block(m))
    append!(r, g0)
    return r
end

function rinit_for(m, g0)
    enh_med  = [median(ENH.noise[p, :])  for p in 1:NNOISE]
    gene_med = [median(GENE.noise[p, :]) for p in 1:NNOISE]
    r = block(m, enh_med, gene_med, g0)
    append!(r, vcat(hyper2_unit(), hyper2_unit(), fill(1.0, m.hblock), fill(1.0, length(g0))))
    for t in 1:NTR
        append!(r, block(m, ENH.noise[:, t], GENE.noise[:, t], g0))
    end
    length(r) == (NHYPER + NTR) * nall(m) || error("internal: rinit is $(length(r)), expected $((NHYPER + NTR) * nall(m))")
    return r
end

priormean_for(m, g0) = rinit_for(m, g0)[1:(2 * nall(m))]

free_firsts(groups) = [first(g) for g in groups]
target_indices(m) = sort([global_index(u, t) for u in (1, 2) for t in m.targets])

function fittedparam_for(m)
    hidden = FIT_HIDDEN ? free_firsts(m.hidden_groups) : Int[]
    base = COUPLE_MODE == "target" ? target_indices(m) : vcat(collect(1:5), collect((NUNIT + 1):(NUNIT + 5)))
    return sort(vcat(base, hidden, free_firsts(m.gamma_groups)))
end

fixedeffects_for(m) = tuple(filter(g -> length(g) > 1, vcat(FIT_HIDDEN ? m.hidden_groups : Vector{Int}[], m.gamma_groups))...)

function priorcv_for(m)
    unit = [1.0, 1.0, 1.0, 1.0, 0.2, 0.25, 0.25, 0.25, 0.1, 0.2, 0.2, 0.1, 0.1]
    cv = vcat(unit, unit, fill(PRIORCV_HIDDEN, m.nh), [1.0, 0.1], fill(PRIORCV_COUPLING, length(m.conns)))
    if COUPLE_MODE == "target"
        for i in target_indices(m)
            cv[i] = PRIORCV_TARGET
        end
    else
        for u in (0, NUNIT)
            cv[u+1:u+4] .= PRIORCV_G
            cv[u+5]      = PRIORCV_INITIATE
        end
    end
    return vcat(cv, fill(2.0, nall(m)))
end

ts() = [(unit = 1, interval = INTERVAL, start = 1.0, t_end = -1.0, zeromedian = true),
        (unit = 2, interval = INTERVAL, start = 1.0, t_end = -1.0, zeromedian = true)]

function spec_for(m)
    nc = length(m.conns)
    g0 = fill(GAMMA0[m.mode], nc)
    Gh = m.hidden == :H2 ? 2 : 3
    return Dict{Symbol,Any}(
        :key => m.key, :label => m.key,
        :datatype => "tracejoint", :datacond => ["enhancer", "gene"],
        :datapath => DATA, :datacol => 3, :zeromedian => [true, true],
        :gene => get(ENV, "SG_GENE", "MYC"), :cell => get(ENV, "SG_CELL", "HBEC"),
        :nalleles => 1, :ejectnumber => 1,
        :probfn => (prob_Gaussian, prob_Gaussian),
        :G => (G, G, Gh), :R => (R, R, 0), :S => (S, S, 0), :insertstep => (INSERTSTEP, INSERTSTEP, 0),
        :transitions => (TRANSITIONS, TRANSITIONS, HIDDEN_TRANSITIONS[m.hidden]),
        # unit 3 is the hidden unit: unit_model (1, 2, 3); direction lives in the connection tuples
        :coupling => ((1, 2, 3), m.conns, fill(m.mode, nc)),
        :noisepriors => [ENH.noise[:, 1], GENE.noise[:, 1], Float64[]],     # the hidden unit has no noise parameters
        :decayrate => DECAY,
        :elongationtime => [ELONG[1], ELONG[2]],
        :hierarchical => (NHYPER, Int[], tuple()),
        :fittedparam => fittedparam_for(m),
        :fixedeffects => fixedeffects_for(m),
        :rinit => rinit_for(m, g0),
        :priormean => priormean_for(m, g0),
        :priorcv => priorcv_for(m),
        :propcv => PROPCV,
        :propcv_rate => 0.0, :propcv_noise => 0.0,
        :proposal_cv_rates => 0.0, :proposal_cv_noise => 0.0,
        :propcv_levels => nothing, :proposal_cv_levels => nothing,
        :method => (Tsit5(), true),
        :trace_specs => ts(), :dwell_specs => [],
        :resultfolder => INTO, :root => RUN_ROOT, :ratetype => "median",
        :nchains => NCHAINS, :samplesteps => SAMPLES, :sample_stride => STRIDE,
        :warmupsteps => parse(Int, get(ENV, "SG_WARMUP", "0")),
        :maxtime => MAXTIME, :temp => 1.0, :burst => false, :writesamples => false,
        :gradient => nothing, :grid => nothing,
    )
end

# --- report ------------------------------------------------------------------------------
const STORED = cld(SAMPLES, STRIDE)
peak_gib(nf) = 2 * NCHAINS * nf * STORED * 8 / 2^30 + 0.5 * (NCHAINS + 1)
const MEMORY_GB_ENV = strip(get(ENV, "SG_MEMORY_GB", ""))
gb_for(nf) = isempty(MEMORY_GB_ENV) ? max(16, ceil(Int, peak_gib(nf) * 1.15 / 8) * 8) : parse(Int, MEMORY_GB_ENV)
const SLURM_TIME = let s = strip(get(ENV, "SG_SLURM_TIME", ""))
    isempty(s) ? string(lpad(ceil(Int, MAXTIME / 3600) + 2, 2, '0'), ":00:00") : s
end

println("\ncouple mode : $COUPLE_MODE     models: $(length(MODELS))     fit hidden rates: $FIT_HIDDEN")
println("decay       : $DECAY     traces: $NTR     blocks: $(NHYPER + NTR)     hidden rates start at $HIDDEN_RATE0")
for m in MODELS
    println(rpad(m.key, 30), rpad(m.spec, 22), " connections ", length(m.conns), "  fitted ", fittedparam_for(m))
    println("    hidden ", m.hidden, m.symmetric ? " symmetric" : "", ": state ", HIDDEN_SOURCE[m.hidden][1], " -> enhancer, state ",
            HIDDEN_SOURCE[m.hidden][2], " -> gene; rates ", join(RATE_LABELS[m.targets], ", "))
    tied = fixedeffects_for(m)
    isempty(tied) || println("    tied: ", join(string.(tied), "  "))
end
println("walltime: $SLURM_TIME")

# --- stage --------------------------------------------------------------------------------
let dest = joinpath(ROOT, "results", INTO)
    if isdir(dest)
        stale = filter(f -> startswith(f, "rates_") && endswith(f, ".txt"), readdir(dest))
        isempty(stale) || error("$dest already holds $(length(stale)) rates_*.txt; use a fresh SG_INTO")
    end
end

println("\nstaging:")
swarmlines = Tuple{String,Int}[]
for m in MODELS
    spec = spec_for(m)
    dir  = joinpath(SWARMROOT, m.key)
    StochasticGene.stage_write_run_specs(
        [spec]; resultfolder = INTO, root = RUN_ROOT, filedir = dir,
        juliafile = "fitscript", swarmfile = "fit", commandfile = "fit.commands",
        nchains = NCHAINS, nthreads = 1, project = PROJECT)
    push!(swarmlines, (dir, gb_for(length(spec[:fittedparam]))))

    got = StochasticGene.read_run_spec(joinpath(ROOT, "results", INTO, "info_" * m.key * ".toml"))
    bad = String[]
    got[:coupling][1] == (1, 2, 3) || push!(bad, "coupling[1] is $(got[:coupling][1]), must be (1, 2, 3)")
    length(got[:coupling][2]) == length(m.conns) || push!(bad, "$(length(got[:coupling][2])) connections, expected $(length(m.conns))")
    all(c -> c[1] == 3, got[:coupling][2]) || push!(bad, "a connection does not start at the hidden unit")
    for k in (:fittedparam, :fixedeffects, :hierarchical, :coupling, :rinit, :priormean, :priorcv, :propcv, :maxtime,
              :datapath, :datacond, :G, :R, :transitions, :method, :trace_specs, :decayrate, :sample_stride, :samplesteps)
        haskey(got, k) || (push!(bad, "$k absent"); continue)
        got[k] == spec[k] || push!(bad, "$k mismatch")
    end
    length(got[:rinit]) == (NHYPER + NTR) * nall(m) || push!(bad, "rinit length $(length(got[:rinit]))")
    println("  ", rpad(m.key, 30), isempty(bad) ? "ok" : "FAILED  " * join(bad, "; "))
    isempty(bad) || error("staged spec for $(m.key) does not match; not submitting")
end

open(joinpath(SWARMROOT, "submit_all.sh"), "w") do io
    println(io, "#!/usr/bin/env bash\nset -euo pipefail")
    for (d, gb) in swarmlines
        rd = ROOT == RUN_ROOT ? d : joinpath(RUN_ROOT, relpath(d, ROOT))
        println(io, "( cd $rd && mkdir -p logs && swarm -f fit.swarm -g $gb -t $(NCHAINS + 2) " *
                    "--time $SLURM_TIME --module julialang --logdir logs )")
    end
end
chmod(joinpath(SWARMROOT, "submit_all.sh"), 0o755)

println("""

results : $(joinpath(ROOT, "results", INTO))
swarms  : $SWARMROOT/<key>/fit.swarm   ($(length(MODELS)) dirs)

Submit with:

    bash $SWARMROOT/submit_all.sh
""")
