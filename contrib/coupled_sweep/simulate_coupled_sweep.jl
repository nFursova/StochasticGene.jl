#!/usr/bin/env julia

"""
Simulate two observed 3301 transcription units (gene1, gene2) over a Cartesian
parameter sweep with StochasticGene's own coupled simulator.

Extends `simulate_shared_hidden_sweep.jl`. The simulator, rate layout,
time-average rate convention, seeds, analysis and output files are the same.
Additions:

- `[hidden] topology = "none" | "H2" | "H3"`. Without a topology key a
  `[hidden]` table means H2, so configurations written for the original script
  run unchanged.
- `[[coupling]]` blocks take `source = "hidden" | "gene1" | "gene2"` and
  `source_state`. Gene sources give direct gene-to-gene and reciprocal coupling.
- Several blocks may regulate the same rate. As in StochasticGene, their
  effects add: rate = base * (1 + sum of active gammas). A point whose combined
  factor is not positive in some promoter state is an error before simulation.
- `rate_reference = "time_average"` is solved on the joint promoter chain.
  For hidden-only coupling this equals the original closed form.
- Burn-in is a multiple of the model's real relaxation time.
- Six FCS multi-tau curves (binary/reporter x EE/GG/EG) computed the same way
  as experimental.R, optional experimental trace lengths, optional scores.

Run:

    julia -t 16 --project=/path/to/StochasticGene simulate_coupled_sweep.jl CONFIG.toml
    julia ... simulate_coupled_sweep.jl CONFIG.toml --check        # validate only
    julia ... simulate_coupled_sweep.jl CONFIG.toml FIRST LAST     # scenario subset

Cross-correlation convention: gene1(t) versus gene2(t + lag). A positive lag
therefore means that gene2 follows gene1.
"""

using StochasticGene
using LinearAlgebra
using Random
using Statistics
using TOML
using Printf

const TRANS_3301 = ((1, 2), (2, 1), (2, 3), (3, 2))
const TRANS_H2 = ((1, 2), (2, 1))
const SIMULATOR_TOTALSTEPS = 100_000_000
const RATE_NAMES = ("R12", "R21", "R23", "R32", "Initiate",
                    "Rshift1", "Rshift2", "Eject", "Decay")
const RATE_INDEX = Dict(name => i for (i, name) in enumerate(RATE_NAMES))
const RATE_REFERENCE_MODES = ("off", "on", "time_average")
const GENES = ("gene1", "gene2")
# StochasticGene's generic R=0 layout reserves two additional slots for RNA
# production and decay. The hidden regulator is not an RNA-producing unit, so
# both are fixed to zero internally and their event clocks remain disabled.
const HIDDEN_NON_RNA_RATES = (0.0, 0.0)
const MIN_MULTITAU_FRAMES = 17

hidden_transitions(nh) = nh == 2 ? TRANS_H2 : TRANS_3301

fmt(x::Real) = isfinite(x) ? @sprintf("%.12g", Float64(x)) : "NA"
csvcell(x) = begin
    s = string(x)
    occursin(r"[,\"\n]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end

asfloat(x, label) = begin
    x isa Real || error("$label must be numeric; got $(repr(x))")
    v = Float64(x)
    isfinite(v) || error("$label must be finite; got $v")
    v
end

asvalues(x, label) = x isa AbstractVector ? [asfloat(v, label) for v in x] :
                     [asfloat(x, label)]

finite_mean(v) = begin
    values = filter(isfinite, Float64.(v))
    isempty(values) ? NaN : mean(values)
end

finite_sd(v) = begin
    values = filter(isfinite, Float64.(v))
    isempty(values) ? NaN : length(values) > 1 ? std(values) : 0.0
end

finite_se(v) = begin
    values = filter(isfinite, Float64.(v))
    length(values) > 1 ? std(values) / sqrt(length(values)) : NaN
end

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

function sweepvalues(key, value)
    if key == "coupling_scheme"
        values = value isa AbstractVector ? String.(value) : [String(value)]
        all(v -> !isempty(strip(v)), values) ||
            error("every sweep.coupling_scheme value must be a nonempty string")
        length(unique(values)) == length(values) ||
            error("sweep.coupling_scheme contains duplicate names")
        return values
    end
    values = asvalues(value, "sweep.$key")
    all(>(0), values) || error("every sweep.$key value must be positive")
    values
end

function require_table(config, name)
    haskey(config, name) || error("configuration is missing [$name]")
    config[name] isa AbstractDict || error("[$name] must be a TOML table")
    config[name]
end

function getpositive(table, key, default=nothing)
    value = haskey(table, key) ? asfloat(table[key], key) :
            default === nothing ? error("missing required value $key") : Float64(default)
    value > 0 || error("$key must be positive; got $value")
    value
end

function base_gene_rates(config, unit)
    table = require_table(config, unit)
    rates = Float64[]
    for name in RATE_NAMES
        haskey(table, name) || error("[$unit] is missing $name")
        push!(rates, getpositive(table, name))
    end
    rates
end

function validate_sweep_key(key)
    key in ("coupling_scheme", "coupling_fold", "hidden_R12", "hidden_R21",
            "hidden_R23", "hidden_R32",
            "common_gene_scale", "gene1_scale", "gene2_scale") && return
    occursin(r"^gene[12]_(R12|R21|R23|R32|Initiate|Rshift1|Rshift2|Eject|Decay)_scale$", key) && return
    endswith(key, "_fold") && return
    error("unsupported [sweep] parameter '$key'")
end

function sweep_points(config)
    table = get(config, "sweep", Dict{String,Any}())
    table isa AbstractDict || error("[sweep] must be a TOML table")
    preferred = ["coupling_scheme", "coupling_fold", "hidden_R12", "hidden_R21",
                 "hidden_R23", "hidden_R32",
                 "common_gene_scale", "gene1_scale", "gene2_scale"]
    keys_all = String.(collect(keys(table)))
    foreach(validate_sweep_key, keys_all)
    keys_ordered = vcat(filter(k -> k in keys_all, preferred),
                        sort(setdiff(keys_all, preferred)))
    points = [Dict{String,Any}()]
    for key in keys_ordered
        values = sweepvalues(key, table[key])
        points = [merge(copy(point), Dict{String,Any}(key => value))
                  for point in points for value in values]
    end
    return points, keys_ordered
end

function hidden_settings(config)
    table = get(config, "hidden", Dict{String,Any}())
    table isa AbstractDict || error("[hidden] must be a TOML table")
    # A [hidden] table without a topology key is the original two-state unit.
    topology = String(get(table, "topology", isempty(table) ? "none" : "H2"))
    topology in ("none", "H2", "H3") ||
        error("hidden.topology must be none, H2 or H3; got $(repr(topology))")
    nh = topology == "none" ? 0 : topology == "H2" ? 2 : 3
    symmetric = Bool(get(table, "symmetric", false))
    symmetric && nh != 3 && error("hidden.symmetric applies only to topology = \"H3\"")
    (; topology, nh, symmetric, table)
end

function check_hidden_sweep_keys(hidden, sweep_keys)
    allowed = hidden.nh == 0 ? String[] :
              hidden.nh == 2 || hidden.symmetric ? ["hidden_R12", "hidden_R21"] :
              ["hidden_R12", "hidden_R21", "hidden_R23", "hidden_R32"]
    for key in sweep_keys
        startswith(key, "hidden_R") && !(key in allowed) &&
            error("sweep key $key does not apply to hidden topology $(hidden.topology)" *
                  (hidden.symmetric ? " with symmetric = true (R23 = R21, R32 = R12)" : ""))
    end
end

# Unit order is [hidden, gene1, gene2], or [gene1, gene2] without a hidden unit.
function unit_layout(nh)
    if nh == 0
        return (; G=(3, 3), R=(3, 3), S=(0, 0), insertstep=(1, 1),
                transitions=(TRANS_3301, TRANS_3301), unit_model=(1, 2),
                observed_units=[1, 2], noise_counts=[4, 4],
                unit=Dict("gene1" => 1, "gene2" => 2))
    end
    (; G=(nh, 3, 3), R=(0, 3, 3), S=(0, 0, 0), insertstep=(1, 1, 1),
     transitions=(hidden_transitions(nh), TRANS_3301, TRANS_3301), unit_model=(1, 2, 3),
     observed_units=[2, 3], noise_counts=[0, 4, 4],
     unit=Dict("hidden" => 1, "gene1" => 2, "gene2" => 3))
end

function parse_coupling(config, sweep_keys, hidden)
    blocks = get(config, "coupling", Any[])
    blocks isa AbstractVector || error("[[coupling]] entries must be TOML array tables")
    isempty(blocks) && error("at least one [[coupling]] block is required; fold 1 is uncoupled")
    defaults = get(config, "coupling_defaults", Dict{String,Any}())
    defaults isa AbstractDict || error("[coupling_defaults] must be a TOML table")
    default_reference = lowercase(String(get(defaults, "rate_reference", "off")))
    default_reference in RATE_REFERENCE_MODES || error(
        "coupling_defaults.rate_reference must be one of $(join(RATE_REFERENCE_MODES, ", "))")
    edges = NamedTuple[]
    for (block_number, block) in enumerate(blocks)
        target = lowercase(String(get(block, "target", "")))
        target in GENES || error("coupling block $block_number target must be gene1 or gene2")
        source = lowercase(String(get(block, "source", "hidden")))
        source in ("hidden", GENES...) ||
            error("coupling block $block_number source must be hidden, gene1 or gene2")
        source == target && error("coupling block $block_number: a gene cannot regulate itself")
        source == "hidden" && hidden.nh == 0 && error(
            "coupling block $block_number uses source = \"hidden\" but hidden.topology is none")
        nstates = source == "hidden" ? hidden.nh : 3
        source_state = if haskey(block, "source_state")
            Int(block["source_state"])
        elseif source == "hidden" && hidden.nh == 2
            2  # the active state of the original two-state hidden unit
        else
            error("coupling block $block_number needs source_state (1 to $nstates)")
        end
        1 <= source_state <= nstates || error(
            "coupling block $block_number source_state must be between 1 and $nstates")
        mode = lowercase(String(get(block, "mode", "")))
        mode in ("activate", "inhibit") ||
            error("coupling block $block_number mode must be activate or inhibit")
        scheme = haskey(block, "scheme") ? strip(String(block["scheme"])) : nothing
        scheme === nothing || !isempty(scheme) ||
            error("coupling block $block_number scheme must be a nonempty string")
        names = get(block, "rates", Any[])
        names isa AbstractVector && !isempty(names) ||
            error("coupling block $block_number needs a nonempty rates array")
        parameter = String(get(block, "fold_parameter", "coupling_fold"))
        fixed_fold = haskey(block, "fold") ? getpositive(block, "fold") : nothing
        rate_reference = lowercase(String(get(block, "rate_reference", default_reference)))
        rate_reference in RATE_REFERENCE_MODES || error(
            "coupling block $block_number rate_reference must be one of " *
            join(RATE_REFERENCE_MODES, ", "))
        if fixed_fold === nothing && !(parameter in sweep_keys) && parameter != "coupling_fold"
            error("coupling block $block_number references absent sweep parameter '$parameter'")
        end
        for raw_name in names
            name = String(raw_name)
            haskey(RATE_INDEX, name) || error(
                "coupling block $block_number has unknown rate '$name'; use $(join(RATE_NAMES, ", "))")
            push!(edges, (; target, source, source_state, rate=name, reaction=RATE_INDEX[name],
                          mode, parameter, fixed_fold, scheme, rate_reference))
        end
    end
    edges
end

function read_csv_rows(path)
    isfile(path) || error("file not found: $path")
    lines = filter(!isempty, strip.(readlines(path)))
    isempty(lines) && error("empty file: $path")
    clean(s) = String(strip(replace(s, "\"" => "")))
    header = clean.(split(lines[1], ','))
    [Dict(zip(header, clean.(split(line, ',')))) for line in lines[2:end]]
end

function read_trace_lengths(path)
    rows = read_csv_rows(path)
    isempty(rows) && error("$path lists no traces")
    haskey(rows[1], "nframes") || error("$path needs an nframes column")
    lengths = [round(Int, parse(Float64, row["nframes"])) for row in rows]
    all(>=(MIN_MULTITAU_FRAMES), lengths) ||
        error("every trace in $path needs at least $MIN_MULTITAU_FRAMES frames")
    lengths
end

# experimental.csv from experimental.R: observable, pair, lag_frames, mean, ...
function read_experimental(path)
    targets = Dict{Tuple{String,String,Int},Float64}()
    for row in read_csv_rows(path)
        value = tryparse(Float64, row["mean"])
        (value === nothing || !isfinite(value)) && continue
        targets[(row["observable"], row["pair"], round(Int, parse(Float64, row["lag_frames"])))] = value
    end
    isempty(targets) && error("$path contains no finite target values")
    targets
end

function simulation_settings(config)
    table = require_table(config, "simulation")
    lengths = nothing
    if haskey(table, "trace_lengths_csv")
        lengths = read_trace_lengths(abspath(String(table["trace_lengths_csv"])))
        haskey(table, "ntraces") && Int(table["ntraces"]) != length(lengths) && error(
            "simulation.ntraces = $(table["ntraces"]) but trace_lengths_csv lists $(length(lengths)) traces")
    end
    nsim = Int(get(table, "nsim", 100))
    ntraces = lengths === nothing ? Int(get(table, "ntraces", 50)) : length(lengths)
    nframes = lengths === nothing ? Int(get(table, "nframes", 216)) : maximum(lengths)
    interval = getpositive(table, "interval_min", 1.0)
    maxlag_min = getpositive(table, "max_lag_min", 60.0)
    max_delay_min = if haskey(table, "burst_max_delay_min")
        asfloat(table["burst_max_delay_min"], "burst_max_delay_min")
    elseif haskey(table, "burst_match_window_min")
        @warn "burst_match_window_min is deprecated; use burst_max_delay_min"
        asfloat(table["burst_match_window_min"], "burst_match_window_min")
    else
        interval
    end
    max_delay_min >= 0 || error("burst_max_delay_min must be nonnegative")
    seed = Int(get(table, "seed", 20260903))
    equilibration = getpositive(table, "equilibration_timescales", 10.0)
    warmup = haskey(table, "warmup_min") ? asfloat(table["warmup_min"], "warmup_min") : nothing
    warmup === nothing || warmup >= 0 || error("warmup_min must be nonnegative")
    nsim > 0 && ntraces > 0 || error("nsim and ntraces must be positive")
    nframes >= MIN_MULTITAU_FRAMES || error("nframes must be at least $MIN_MULTITAU_FRAMES for multi-tau")
    maxlag = min(floor(Int, maxlag_min / interval), nframes - 3)
    max_delay = round(Int, max_delay_min / interval)
    noise = Float64.(get(table, "noise_params", [50.0, 15.0, 200.0, 20.0]))
    length(noise) == 4 || error("simulation.noise_params must contain four values")
    all(isfinite, noise) || error("simulation.noise_params contains a nonfinite value")
    epsilon = asfloat(get(table, "epsilon", 1e-7), "epsilon")
    epsilon >= 0 || error("epsilon must be nonnegative")
    output_dir = abspath(String(get(table, "output_dir", "coupled_sweep_output")))
    save_baseline_binary = Bool(get(table, "save_baseline_binary", true))
    return (; nsim, ntraces, nframes, lengths, interval, maxlag, max_delay, seed, equilibration,
            warmup, noise, epsilon, output_dir, save_baseline_binary)
end

function scoring_settings(config)
    haskey(config, "scoring") || return nothing
    table = require_table(config, "scoring")
    haskey(table, "experimental_csv") || error("[scoring] needs experimental_csv")
    path = abspath(String(table["experimental_csv"]))
    (; path, targets=read_experimental(path), max_lag_min=getpositive(table, "max_lag_min", 90.0))
end

pointvalue(point, key, default) = get(point, key, Float64(default))

function effective_gene_rates(base, unit, point)
    common = pointvalue(point, "common_gene_scale", 1.0)
    unit_scale = pointvalue(point, "$(unit)_scale", 1.0)
    [base[i] * common * unit_scale *
     pointvalue(point, "$(unit)_$(RATE_NAMES[i])_scale", 1.0) for i in eachindex(base)]
end

hidden_rate(point, table, key) =
    haskey(point, "hidden_$key") ? Float64(point["hidden_$key"]) : getpositive(table, key)

function hidden_rates_for_point(hidden, point)
    hidden.nh == 0 && return Float64[]
    r12 = hidden_rate(point, hidden.table, "R12")
    r21 = hidden_rate(point, hidden.table, "R21")
    hidden.nh == 2 && return [r12, r21]
    hidden.symmetric && return [r12, r21, r21, r12]
    [r12, r21, hidden_rate(point, hidden.table, "R23"), hidden_rate(point, hidden.table, "R32")]
end

# ---------------------------------------------------------------------------
# Joint promoter chain: hidden state h, gene1 G state a, gene2 G state b.
# Promoter transitions never depend on reporter occupancy, so this small chain
# (at most 27 states) gives exact stationary source occupancies, time averages
# of coupled rates, and the slowest promoter relaxation time.
# ---------------------------------------------------------------------------

npromoter(nh) = 9 * max(nh, 1)
pindex(h, a, b) = a + 3 * (b - 1) + 9 * (h - 1)
pstate(p) = (div(p - 1, 9) + 1, mod(p - 1, 3) + 1, mod(div(p - 1, 3), 3) + 1)

function coupling_factor(p, gene, reaction, edges, gamma)
    h, a, b = pstate(p)
    factor = 1.0
    for (edge, g) in zip(edges, gamma)
        edge.target == gene && edge.reaction == reaction || continue
        state = edge.source == "hidden" ? h : edge.source == "gene1" ? a : b
        state == edge.source_state && (factor += g)
    end
    factor
end

function promoter_generator(nh, hidden_rates, basal, edges, gamma)
    n = npromoter(nh)
    Q = zeros(n, n)
    for p in 1:n
        h, a, b = pstate(p)
        for (k, (from, to)) in enumerate(TRANS_3301)
            if a == from
                rate = basal["gene1"][k] * coupling_factor(p, "gene1", k, edges, gamma)
                Q[pindex(h, to, b), p] += rate
                Q[p, p] -= rate
            end
            if b == from
                rate = basal["gene2"][k] * coupling_factor(p, "gene2", k, edges, gamma)
                Q[pindex(h, a, to), p] += rate
                Q[p, p] -= rate
            end
        end
        nh == 0 && continue
        for (k, (from, to)) in enumerate(hidden_transitions(nh))
            h == from || continue
            Q[pindex(to, a, b), p] += hidden_rates[k]
            Q[p, p] -= hidden_rates[k]
        end
    end
    Q
end

function stationary(Q)
    n = size(Q, 1)
    A = copy(Q)
    A[end, :] .= 1.0
    b = zeros(n)
    b[end] = 1.0
    prob = A \ b
    all(isfinite, prob) && minimum(prob) > -1e-10 ||
        error("could not solve for the stationary promoter distribution")
    prob = max.(prob, 0.0)
    prob ./ sum(prob)
end

time_average(prob, gene, k, basal, edges, gamma) =
    basal[gene][k] * sum(prob[p] * coupling_factor(p, gene, k, edges, gamma) for p in eachindex(prob))

"""
Return base (all sources off) rates for each gene.

- off: base = supplied rate.
- on: base = supplied / multiplier (one regulating block only).
- time_average: base chosen so the stationary time average of the coupled rate
  equals the supplied rate. With hidden sources only, one pass gives the
  original closed form supplied / ((1 - p) + p * multiplier). With gene sources
  the source occupancy depends on the base rates, so the pass is repeated.
"""
function calibrate(nh, hidden_rates, input, edges, gamma, multipliers)
    basal = Dict(g => copy(input[g]) for g in GENES)
    reference = Dict{Tuple{String,Int},String}()
    for edge in edges
        key = (edge.target, edge.reaction)
        get(reference, key, edge.rate_reference) == edge.rate_reference || error(
            "$(edge.target):$(edge.rate) is regulated by blocks with different rate_reference values")
        reference[key] = edge.rate_reference
    end
    for (key, mode) in reference
        mode == "on" || continue
        idx = findall(e -> (e.target, e.reaction) == key, edges)
        length(idx) == 1 || error("rate_reference = \"on\" is ambiguous for " *
            "$(key[1]):$(RATE_NAMES[key[2]]), which $(length(idx)) blocks regulate; use off or time_average")
        basal[key[1]][key[2]] = input[key[1]][key[2]] / multipliers[only(idx)]
    end
    averaged = [key for (key, mode) in reference if mode == "time_average"]
    iterations = 0
    if !isempty(averaged)
        converged = false
        for iteration in 1:500
            iterations = iteration
            prob = stationary(promoter_generator(nh, hidden_rates, basal, edges, gamma))
            change = 0.0
            damping = iteration <= 100 ? 1.0 : 0.5
            for (g, k) in averaged
                ratio = input[g][k] / time_average(prob, g, k, basal, edges, gamma)
                change = max(change, abs(log(ratio)))
                basal[g][k] *= ratio^damping
            end
            change < 1e-12 && (converged = true; break)
        end
        converged || error("time-average calibration did not converge in 500 iterations")
    end
    prob = stationary(promoter_generator(nh, hidden_rates, basal, edges, gamma))
    average = Dict(g => [time_average(prob, g, k, basal, edges, gamma) for k in eachindex(RATE_NAMES)]
                   for g in GENES)
    for (g, k) in averaged
        abs(average[g][k] / input[g][k] - 1) <= 1e-8 ||
            error("time-average calibration missed $g:$(RATE_NAMES[k])")
    end
    (; basal, average, probabilities=prob, iterations)
end

function model_for_point(point, hidden, layout, edge_specs, base)
    hidden_rates = hidden_rates_for_point(hidden, point)
    input = Dict(g => effective_gene_rates(base[g], g, point) for g in GENES)
    selected_scheme = get(point, "coupling_scheme", nothing)
    edges = selected_scheme === nothing ? edge_specs :
            filter(edge -> edge.scheme === nothing || edge.scheme == selected_scheme, edge_specs)
    isempty(edges) && error("coupling scheme $(repr(selected_scheme)) selects no rates")
    folds = [edge.fixed_fold === nothing ? Float64(pointvalue(point, edge.parameter, 1.0)) :
             edge.fixed_fold for edge in edges]
    all(>(0), folds) || error("coupling folds must be positive")
    multipliers = [edge.mode == "activate" ? fold : inv(fold) for (edge, fold) in zip(edges, folds)]
    gamma = multipliers .- 1
    np = npromoter(hidden.nh)
    for (g, k) in unique([(edge.target, edge.reaction) for edge in edges]), p in 1:np
        factor = coupling_factor(p, g, k, edges, gamma)
        factor > 0 || error("combined coupling factor on $g:$(RATE_NAMES[k]) is $(fmt(factor)) " *
            "in promoter state (hidden, gene1, gene2) = $(pstate(p)); every coupled rate must stay positive")
    end
    calibration = calibrate(hidden.nh, hidden_rates, input, edges, gamma, multipliers)
    basal = calibration.basal
    decay_rates = sort(abs.(real.(eigvals(promoter_generator(hidden.nh, hidden_rates, basal, edges, gamma)))))
    decay_rates[2] > 0 || error("the promoter chain has no finite relaxation time")
    # Time for a polymerase to initiate and leave the 3 reporter positions at the
    # slowest coupled rates: an upper bound on how long reporters remember the start.
    relaxation_reporter = maximum(
        sum(1 / minimum(basal[g][k] * coupling_factor(p, g, k, edges, gamma) for p in 1:np) for k in 5:8)
        for g in GENES)
    (; hidden_rates, input, basal, average=calibration.average,
       promoter_probabilities=calibration.probabilities,
       calibration_iterations=calibration.iterations, edges, folds, multipliers, gamma,
       connections=[(layout.unit[e.source], e.source_state, layout.unit[e.target], e.reaction) for e in edges],
       edge_names=["$(e.source)$(e.source_state)->$(e.target):$(e.rate):$(e.mode)" for e in edges],
       relaxation_promoter=1 / decay_rates[2], relaxation_reporter)
end

# The trace starts in promoter state 1 with no reporters; the deviation from
# stationarity decays as exp(-t / slowest relaxation time).
warmup_for(model, settings) = settings.warmup !== nothing ? settings.warmup :
    settings.equilibration * max(model.relaxation_promoter, model.relaxation_reporter)

# ---------------------------------------------------------------------------
# Simulation (StochasticGene)
# ---------------------------------------------------------------------------

function simulate_one(model, layout, settings)
    rates = Float64[]
    if !isempty(model.hidden_rates)
        append!(rates, model.hidden_rates)
        append!(rates, HIDDEN_NON_RNA_RATES)
    end
    for g in GENES
        append!(rates, model.basal[g])
        append!(rates, settings.noise)
    end
    append!(rates, model.gamma)
    totaltime = (settings.nframes + 2) * settings.interval
    StochasticGene.simulate_trace_vector(
        rates, layout.transitions, layout.G, layout.R, layout.S, layout.insertstep,
        (layout.unit_model, model.connections), settings.interval, totaltime,
        settings.ntraces; onstates=[Int[] for _ in layout.G],
        probfn=StochasticGene.prob_Gaussian, noiseparams=layout.noise_counts,
        observed_units=layout.observed_units, col=3, totalsteps=SIMULATOR_TOTALSTEPS,
        warmuptime=warmup_for(model, settings), equilibration_timescales=1.0)
end

# Column 3 of StochasticGene traces is the reporter count. Truncate trace i to
# its experimental length when trace_lengths_csv is given.
function reporter_traces(raw, settings)
    traces = Matrix{Int}[]
    for (i, trace) in enumerate(raw)
        size(trace, 2) == 2 || error("trace $i has $(size(trace, 2)) columns; expected 2")
        n = settings.lengths === nothing ? settings.nframes : settings.lengths[i]
        size(trace, 1) >= n || error("trace $i has $(size(trace, 1)) frames; expected $n")
        push!(traces, round.(Int, trace[1:n, :]))
    end
    traces
end

# ---------------------------------------------------------------------------
# Analysis (unchanged from the original unless noted)
# ---------------------------------------------------------------------------

function pearson_at_lag(x, y, lag)
    n = length(x)
    if lag >= 0
        m, xi, yi = n - lag, 1, 1 + lag
    else
        m, xi, yi = n + lag, 1 - lag, 1
    end
    m >= 3 || return NaN
    sx = 0.0; sy = 0.0
    @inbounds for j in 0:(m - 1)
        sx += x[xi + j]; sy += y[yi + j]
    end
    mx = sx / m; my = sy / m
    sxx = 0.0; syy = 0.0; sxy = 0.0
    @inbounds for j in 0:(m - 1)
        a = x[xi + j] - mx; b = y[yi + j] - my
        sxx += a * a; syy += b * b; sxy += a * b
    end
    sxx > 0 && syy > 0 ? sxy / sqrt(sxx * syy) : NaN
end

function connected_at_lag(x, y, lag)
    n = length(x)
    if lag >= 0
        m, xi, yi = n - lag, 1, 1 + lag
    else
        m, xi, yi = n + lag, 1 - lag, 1
    end
    m > 0 || return NaN
    sx = 0.0; sy = 0.0; sxy = 0.0
    @inbounds for j in 0:(m - 1)
        a = Float64(x[xi + j]); b = Float64(y[yi + j])
        sx += a; sy += b; sxy += a * b
    end
    mx, my = sx / m, sy / m
    mx > 0 && my > 0 ? sxy / m / (mx * my) - 1 : NaN
end

# FCS multitau implementation ported from fcs_multitau_xcor_final.R via
# run_coupled_xcor_sim_final_v2.jl (verbatim from simulate_shared_hidden_sweep.jl).
idl_round(x) = sign(x) * floor(abs(x) + 0.5)

struct MTCurve
    lag::Vector{Float64}
    Gn::Vector{Float64}
    Md::Vector{Float64}
    Me::Vector{Float64}
    n::Vector{Float64}
end

function bg_level(x::Vector{Float64}, cutoff::Float64)
    values = filter(t -> isfinite(t) && t >= 0, x)
    isempty(values) && return 0.0
    mx = maximum(values)
    mx <= 0 && return 0.0
    hi = ceil(mx) + 1.0
    nb = min(Int(ceil(mx)) + 1, 1_000_000)
    width = hi / nb
    counts = zeros(Int, nb)
    for t in values
        k = clamp(Int(ceil(t / width)), 1, nb)
        counts[k] += 1
    end
    run = 0
    index = 0
    for k in eachindex(counts)
        run += counts[k]
        run / length(values) <= cutoff ? (index = k) : break
    end
    index == 0 ? 0.0 : index * width
end

function multi_tau_corr(left::Vector{Float64}, right::Vector{Float64};
                        M::Int=4, bg_correct::Bool=true,
                        bg_cutoff::Float64=0.001, lag_sign::Int=1)
    n = length(left)
    length(right) == n || error("left and right traces differ in length")
    selected = n - 1
    selected < 4 * M && return nothing
    N = Int(idl_round(log2((selected + 1) / M)))
    N < 2 && return nothing
    nchan = M * N
    epsilon = bg_correct ? 1e-7 : 0.0
    L = left[1:selected] .+ epsilon
    Rv = right[1:selected] .+ epsilon
    background_factor = 1.0
    if bg_correct
        lb, rb = bg_level(L, bg_cutoff), bg_level(Rv, bg_cutoff)
        mL, mR = mean(L), mean(Rv)
        if mL - lb > 0 && mR - rb > 0
            background_factor = (mL / (mL - lb)) * (mR / (mR - rb))
        end
    end
    Go = zeros(nchan); Md = zeros(nchan); Me = zeros(nchan)
    np = zeros(nchan); lags = zeros(nchan)
    k = 1
    lag = 1
    for _ in 1:(2 * M)
        j_end = selected - lag
        if j_end >= 1
            @inbounds for j in 1:j_end
                a, b = L[j], Rv[j + lag]
                Go[k] += a * b; Md[k] += a; Me[k] += b
            end
            np[k] = j_end
        end
        lags[k] = lag_sign * lag
        k += 1; lag += 1
    end
    Lb, Rb = copy(L), copy(Rv)
    points = selected
    for octave in 1:(N - 2)
        points -= points % 2
        half = div(points, 2)
        Lb2 = Vector{Float64}(undef, half)
        Rb2 = Vector{Float64}(undef, half)
        @inbounds for j in 1:half
            Lb2[j] = Lb[2j - 1] + Lb[2j]
            Rb2[j] = Rb[2j - 1] + Rb[2j]
        end
        Lb, Rb = Lb2, Rb2
        points = half
        for i in 1:M
            lag = M + i
            j_end = points - lag
            if j_end >= 1
                @inbounds for j in 1:j_end
                    a, b = Lb[j], Rb[j + lag]
                    Go[k] += a * b; Md[k] += a; Me[k] += b
                end
                np[k] = j_end
            end
            lags[k] = lag_sign * lag * 2.0^octave
            k += 1
        end
    end
    used = k - 1
    Gn = fill(NaN, used)
    @inbounds for j in 1:used
        if np[j] > 0 && isfinite(Go[j]) && Md[j] > 0 && Me[j] > 0
            Gn[j] = ((np[j] * Go[j]) / (Md[j] * Me[j]) - 1) * background_factor
        end
    end
    MTCurve(lags[1:used], Gn, Md[1:used], Me[1:used], np[1:used])
end

function multi_tau_xcorr(left::Vector{Float64}, right::Vector{Float64}; kwargs...)
    forward = multi_tau_corr(left, right; lag_sign=1, kwargs...)
    backward = multi_tau_corr(right, left; lag_sign=-1, kwargs...)
    (forward === nothing || backward === nothing) && return nothing
    MTCurve(vcat(forward.lag, backward.lag), vcat(forward.Gn, backward.Gn),
            vcat(forward.Md, backward.Md), vcat(forward.Me, backward.Me),
            vcat(forward.n, backward.n))
end

function coulon_multitau_mean(curves::Vector{MTCurve})
    isempty(curves) && return NamedTuple[]
    common = copy(first(curves).lag)
    for curve in Iterators.drop(curves, 1)
        common = intersect(common, curve.lag)
    end
    sort!(common)
    isempty(common) && return NamedTuple[]
    ncurves = length(curves)
    weights = fill(NaN, length(common), ncurves)
    Md1 = zeros(ncurves); Me1 = zeros(ncurves); n1 = zeros(ncurves)
    for i in eachindex(curves)
        curve = curves[i]
        Md1[i], Me1[i], n1[i] = curve.Md[1], curve.Me[1], curve.n[1]
        positions = Dict(lag => j for (j, lag) in enumerate(curve.lag))
        n1[i] > 0 || continue
        for j in eachindex(common)
            index = get(positions, common[j], 0)
            index == 0 && continue
            weights[j, i] = curve.Gn[index] * Md1[i] * Me1[i] / n1[i]^2
        end
    end
    sn = sum(n1)
    denominator = sn > 0 ? (sum(Md1) / sn) * (sum(Me1) / sn) : 0.0
    Gn = fill(NaN, length(common))
    if denominator != 0
        for j in eachindex(common)
            values = filter(isfinite, weights[j, :])
            isempty(values) || (Gn[j] = mean(values) / denominator)
        end
    end
    [(; lag_frames=common[j], Gn=Gn[j], n_curves=ncurves) for j in eachindex(common)]
end

# Original output: binary gene1 x gene2 curve, background correction on, |lag| <= max_lag_min.
function fcs_multitau(trace_on, settings)
    curves = MTCurve[]
    for (gene1, gene2) in trace_on
        curve = multi_tau_xcorr(Float64.(gene1), Float64.(gene2))
        curve === nothing || push!(curves, curve)
    end
    rows = coulon_multitau_mean(curves)
    [(; row..., lag_min=row.lag_frames * settings.interval)
     for row in rows if abs(row.lag_frames) <= settings.maxlag]
end

# New: the six curves exactly as experimental.R computes its targets
# (epsilon added to the trace, background correction off, all lags).
function multitau_all(traces, settings)
    rows = NamedTuple[]
    for observable in ("binary", "reporter"), (pair, i, j) in (("EE", 1, 1), ("GG", 2, 2), ("EG", 1, 2))
        curves = MTCurve[]
        for trace in traces
            x = Float64.(trace[:, i]); y = Float64.(trace[:, j])
            if observable == "binary"
                x = Float64.(x .> 0); y = Float64.(y .> 0)
            end
            x .+= settings.epsilon; y .+= settings.epsilon
            curve = multi_tau_xcorr(x, y; bg_correct=false)
            curve === nothing || push!(curves, curve)
        end
        for row in coulon_multitau_mean(curves)
            push!(rows, (; observable, pair, lag_frames=round(Int, row.lag_frames),
                         lag_min=row.lag_frames * settings.interval, Gn=row.Gn,
                         n_curves=row.n_curves))
        end
    end
    rows
end

burst_starts(on) = [i for i in 2:length(on) if on[i] && !on[i - 1]]

"""
For every gene-1 anchor start, select the nearest gene-2 reader start. Reader
events may be selected by more than one anchor, matching the event-aligned R
analysis. Ties select the earlier reader event. Only signed delays within
`max_delay` frames are classified; all other anchors are unpaired.
"""
function classify_nearest_burst_starts(starts1, starts2, max_delay)
    before = same = after = within1 = within2 = unpaired1 = 0
    selected2 = Int[]
    for anchor in starts1
        if isempty(starts2)
            unpaired1 += 1
            continue
        end
        best = firstindex(starts2)
        delta = starts2[best] - anchor
        for j in eachindex(starts2)
            candidate = starts2[j] - anchor
            if abs(candidate) < abs(delta)
                best, delta = j, candidate
            end
        end
        if abs(delta) > max_delay
            unpaired1 += 1
        else
            push!(selected2, best)
            abs(delta) <= 1 && (within1 += 1)
            abs(delta) <= 2 && (within2 += 1)
            delta > 0 ? (before += 1) : delta < 0 ? (after += 1) : (same += 1)
        end
    end
    (; before, same, after, within1, within2, matched=before + same + after, unpaired1,
     unpaired2=length(starts2) - length(unique(selected2)))
end

percent(count, denominator) = denominator > 0 ? 100 * count / denominator : NaN

function analyze_simulation(traces, settings)
    ntraces = length(traces)
    on1_total = on2_total = joint_total = 0
    reporters1_total = reporters2_total = 0
    nframes_total = 0
    bursts1 = bursts2 = before = same = after = within1 = within2 = matched = 0
    unpaired1 = unpaired2 = 0
    xcor = NamedTuple[]
    trace_on = Tuple{Vector{Bool},Vector{Bool}}[]
    for trace in traces
        on1 = Vector{Bool}(trace[:, 1] .> 0)
        on2 = Vector{Bool}(trace[:, 2] .> 0)
        push!(trace_on, (on1, on2))
        on1_total += count(on1); on2_total += count(on2)
        joint_total += count(on1 .& on2); nframes_total += length(on1)
        reporters1_total += sum(trace[:, 1]); reporters2_total += sum(trace[:, 2])
        starts1, starts2 = burst_starts(on1), burst_starts(on2)
        pairing = classify_nearest_burst_starts(starts1, starts2, settings.max_delay)
        bursts1 += length(starts1); bursts2 += length(starts2)
        before += pairing.before; same += pairing.same; after += pairing.after
        within1 += pairing.within1; within2 += pairing.within2
        matched += pairing.matched; unpaired1 += pairing.unpaired1
        unpaired2 += pairing.unpaired2
    end
    for lag in -settings.maxlag:settings.maxlag
        pearson = [pearson_at_lag(x, y, lag) for (x, y) in trace_on]
        connected = [connected_at_lag(x, y, lag) for (x, y) in trace_on]
        push!(xcor, (; lag_frames=lag, lag_min=lag * settings.interval,
                     pearson_mean=finite_mean(pearson), pearson_sd=finite_sd(pearson),
                     connected_mean=finite_mean(connected),
                     connected_sd=finite_sd(connected),
                     n_pearson=count(isfinite, pearson),
                     n_connected=count(isfinite, connected)))
    end
    metrics = (; n_traces=ntraces, n_frames=nframes_total,
        on_fraction_gene1=on1_total / nframes_total,
        on_fraction_gene2=on2_total / nframes_total,
        joint_on_fraction=joint_total / nframes_total,
        gene1_bursts=bursts1, gene2_bursts=bursts2, matched_bursts=matched,
        gene1_before_gene2_pct=percent(before, bursts1),
        gene1_same_as_gene2_pct=percent(same, bursts1),
        gene1_after_gene2_pct=percent(after, bursts1),
        gene1_start_within_1_frame_gene2_pct=percent(within1, bursts1),
        gene1_start_within_2_frames_gene2_pct=percent(within2, bursts1),
        gene1_unpaired_pct=percent(unpaired1, bursts1),
        before_count=before, same_count=same, after_count=after,
        within_1_frame_count=within1, within_2_frames_count=within2,
        unpaired_gene1_count=unpaired1, unpaired_gene2_count=unpaired2,
        mean_reporters_gene1=reporters1_total / nframes_total,
        mean_reporters_gene2=reporters2_total / nframes_total)
    return (; metrics, xcor, multitau=fcs_multitau(trace_on, settings),
            multitau_all=multitau_all(traces, settings))
end

function multitau_summary(results)
    order = Tuple{String,String,Int}[]
    values = Dict{Tuple{String,String,Int},Vector{Float64}}()
    lag_min = Dict{Tuple{String,String,Int},Float64}()
    for rows in results, row in rows
        key = (row.observable, row.pair, row.lag_frames)
        if !haskey(values, key)
            push!(order, key); values[key] = Float64[]; lag_min[key] = row.lag_min
        end
        push!(values[key], row.Gn)
    end
    [(; observable=k[1], pair=k[2], lag_frames=k[3], lag_min=lag_min[k],
        mean=finite_mean(values[k]), sd=finite_sd(values[k]), se=finite_se(values[k]),
        n=count(isfinite, values[k])) for k in order]
end

function score_curve(curve, targets, observable, pair, max_lag)
    target(r) = get(targets, (observable, pair, r.lag_frames), NaN)
    usable = [r for r in curve if isfinite(r.mean) && isfinite(target(r))]
    region(inside) = (finite_mean([r.mean for r in curve if inside(r)]),
                      finite_mean([target(r) for r in usable if inside(r)]))
    center, center_target = region(r -> abs(r.lag_frames) == 1)
    left, left_target = region(r -> -60 - 1e-9 <= r.lag_min <= -3 + 1e-9)
    right, right_target = region(r -> 3 - 1e-9 <= r.lag_min <= 60 + 1e-9)
    both = [r for r in usable if abs(r.lag_min) <= max_lag + 1e-9]
    x = [r.mean for r in both]; y = [target(r) for r in both]
    nrmse = length(x) >= 2 && sum(abs2, y) > 0 ? sqrt(sum(abs2, x .- y) / sum(abs2, y)) : NaN
    ccc = NaN
    if length(x) >= 2
        mx, my = mean(x), mean(y)
        vx, vy = mean((x .- mx) .^ 2), mean((y .- my) .^ 2)
        denominator = vx + vy + (mx - my)^2
        ccc = denominator > 0 ? 2 * mean((x .- mx) .* (y .- my)) / denominator : NaN
    end
    (; n_bins=length(x), nrmse, ccc, center, center_target, left, left_target, right, right_target)
end

# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------

function write_binary_traces(path, traces, settings)
    open(path, "w") do io
        println(io, "trace,frame,time_min,gene1_on,gene2_on")
        for (trace_number, trace) in enumerate(traces), frame in axes(trace, 1)
            println(io, trace_number, ",", frame, ",", fmt(frame * settings.interval), ",",
                    Int(trace[frame, 1] > 0), ",", Int(trace[frame, 2] > 0))
        end
    end
end

function write_sweep_points(path, points, keys_ordered)
    open(path, "w") do io
        println(io, join(vcat("scenario_id", keys_ordered), ','))
        for (i, point) in enumerate(points)
            values = [point[key] isa Real ? fmt(point[key]) : csvcell(point[key])
                      for key in keys_ordered]
            println(io, join(vcat(string(i), values), ','))
        end
    end
end

function write_design(path, points, models, keys_ordered, hidden, settings)
    header = vcat("scenario_id", keys_ordered,
        ["hidden_topology", "hidden_R12", "hidden_R21", "hidden_R23", "hidden_R32"],
        ["hidden_state$(h)_occupancy" for h in 1:hidden.nh],
        ["coupling_edges", "coupling_folds", "coupling_gammas", "coupling_rate_references",
         "calibration_iterations", "promoter_relaxation_min", "reporter_relaxation_min",
         "warmup_min"],
        ["$(g)_$(name)_$(kind)" for kind in ("input_reference", "base", "time_average")
         for g in GENES for name in RATE_NAMES])
    open(path, "w") do io
        println(io, join(header, ','))
        for (i, (point, model)) in enumerate(zip(points, models))
            prob = model.promoter_probabilities
            values = String[string(i)]
            append!(values, [point[key] isa Real ? fmt(point[key]) : csvcell(point[key])
                             for key in keys_ordered])
            push!(values, hidden.topology)
            append!(values, fmt.(vcat(model.hidden_rates, fill(NaN, 4 - length(model.hidden_rates)))))
            append!(values, [fmt(sum(prob[p] for p in eachindex(prob) if pstate(p)[1] == h))
                             for h in 1:hidden.nh])
            append!(values, csvcell.([join(model.edge_names, ';'), join(fmt.(model.folds), ';'),
                                      join(fmt.(model.gamma), ';'),
                                      join([edge.rate_reference for edge in model.edges], ';')]))
            append!(values, [string(model.calibration_iterations), fmt(model.relaxation_promoter),
                             fmt(model.relaxation_reporter), fmt(warmup_for(model, settings))])
            for rates in (model.input, model.basal, model.average), g in GENES
                append!(values, fmt.(rates[g]))
            end
            println(io, join(values, ','))
        end
    end
end

const METRIC_FIELDS = (:n_traces, :n_frames, :on_fraction_gene1, :on_fraction_gene2,
    :joint_on_fraction, :gene1_bursts, :gene2_bursts, :matched_bursts,
    :gene1_before_gene2_pct, :gene1_same_as_gene2_pct, :gene1_after_gene2_pct,
    :gene1_start_within_1_frame_gene2_pct,
    :gene1_start_within_2_frames_gene2_pct,
    :gene1_unpaired_pct, :before_count, :same_count, :after_count,
    :within_1_frame_count, :within_2_frames_count,
    :unpaired_gene1_count, :unpaired_gene2_count,
    :mean_reporters_gene1, :mean_reporters_gene2)

const REPORTED_METRICS = (:on_fraction_gene1, :on_fraction_gene2, :joint_on_fraction,
    :gene1_before_gene2_pct, :gene1_same_as_gene2_pct,
    :gene1_after_gene2_pct, :gene1_start_within_1_frame_gene2_pct,
    :gene1_start_within_2_frames_gene2_pct, :gene1_unpaired_pct,
    :mean_reporters_gene1, :mean_reporters_gene2)

function write_metrics_rows(io, scenario, results)
    for sim in eachindex(results)
        values = [fmt(Float64(getproperty(results[sim], field))) for field in METRIC_FIELDS]
        println(io, scenario, ",", sim, ",", join(values, ','))
    end
end

function write_metric_summary_row(io, scenario, results)
    values = String[]
    for field in REPORTED_METRICS
        sample = [Float64(getproperty(result, field)) for result in results]
        append!(values, [fmt(finite_mean(sample)), fmt(finite_sd(sample))])
    end
    println(io, scenario, ",", join(values, ','))
end

function write_xcor_rows(io, scenario, results)
    for sim in eachindex(results), row in results[sim]
        println(io, scenario, ",", sim, ",", row.lag_frames, ",",
                fmt(row.lag_min), ",", fmt(row.pearson_mean), ",",
                fmt(row.pearson_sd), ",", fmt(row.connected_mean), ",",
                fmt(row.connected_sd), ",", row.n_pearson, ",", row.n_connected)
    end
end

function write_xcor_summary_rows(io, scenario, results)
    for lag_index in eachindex(first(results))
        rows = [result[lag_index] for result in results]
        p = [row.pearson_mean for row in rows]
        c = [row.connected_mean for row in rows]
        row = first(rows)
        println(io, scenario, ",", row.lag_frames, ",", fmt(row.lag_min), ",",
                fmt(finite_mean(p)), ",", fmt(finite_sd(p)), ",",
                fmt(finite_mean(c)), ",", fmt(finite_sd(c)))
    end
end

function write_multitau_rows(io, scenario, results)
    for sim in eachindex(results), row in results[sim]
        println(io, scenario, ",", sim, ",", fmt(row.lag_frames), ",",
                fmt(row.lag_min), ",", fmt(row.Gn), ",", row.n_curves)
    end
end

function write_multitau_summary_rows(io, scenario, results)
    isempty(first(results)) && return
    for lag_index in eachindex(first(results))
        rows = [result[lag_index] for result in results]
        all(row -> row.lag_frames == first(rows).lag_frames, rows) ||
            error("FCS multitau lag grids differ among simulations in scenario $scenario")
        values = [row.Gn for row in rows]
        row = first(rows)
        println(io, scenario, ",", fmt(row.lag_frames), ",", fmt(row.lag_min), ",",
                fmt(finite_mean(values)), ",", fmt(finite_sd(values)), ",",
                count(isfinite, values))
    end
end

function write_multitau_all_rows(io, scenario, results)
    for sim in eachindex(results), row in results[sim]
        println(io, join((scenario, sim, row.observable, row.pair, row.lag_frames,
                          fmt(row.lag_min), fmt(row.Gn), row.n_curves), ','))
    end
end

function write_multitau_all_summary_rows(io, scenario, summary)
    for row in summary
        println(io, join((scenario, row.observable, row.pair, row.lag_frames, fmt(row.lag_min),
                          fmt(row.mean), fmt(row.sd), fmt(row.se), row.n), ','))
    end
end

function write_score_rows(io, scenario, summary, scoring)
    for observable in ("binary", "reporter"), pair in ("EE", "GG", "EG")
        curve = [row for row in summary if row.observable == observable && row.pair == pair]
        s = score_curve(curve, scoring.targets, observable, pair, scoring.max_lag_min)
        println(io, join((scenario, observable, pair, s.n_bins, fmt(s.nrmse), fmt(s.ccc),
                          fmt(s.center), fmt(s.center_target), fmt(s.left), fmt(s.left_target),
                          fmt(s.right), fmt(s.right_target)), ','))
    end
end

const OUTPUT_HEADERS = [
    :metrics => ("metrics", "scenario_id,sim," * join(string.(METRIC_FIELDS), ',')),
    :xcor => ("xcor", "scenario_id,sim,lag_frames,lag_min,pearson_mean,pearson_sd," *
                      "connected_mean,connected_sd,n_pearson,n_connected"),
    :metrics_summary => ("metrics_summary", "scenario_id," *
        join(reduce(vcat, (["$(f)_mean", "$(f)_sd"] for f in REPORTED_METRICS)), ',')),
    :xcor_summary => ("xcor_summary", "scenario_id,lag_frames,lag_min,pearson_mean,pearson_sd," *
                                      "connected_mean,connected_sd"),
    :multitau => ("xcor_multitau", "scenario_id,sim,lag_frames,lag_min,Gn,n_curves"),
    :multitau_summary => ("xcor_multitau_summary",
                          "scenario_id,lag_frames,lag_min,Gn_mean,Gn_sd,n_simulations"),
    :multitau_all => ("multitau_all", "scenario_id,sim,observable,pair,lag_frames,lag_min,Gn,n_curves"),
    :multitau_all_summary => ("multitau_all_summary",
        "scenario_id,observable,pair,lag_frames,lag_min,Gn_mean,Gn_sd,Gn_se,n_simulations"),
    :scores => ("scores", "scenario_id,observable,pair,n_bins,nrmse,ccc,center,center_target," *
                          "left,left_target,right,right_target"),
]

function open_output_files(output_dir, suffix, scoring)
    ios = Dict{Symbol,IO}()
    try
        for (key, (name, header)) in OUTPUT_HEADERS
            key == :scores && scoring === nothing && continue
            ios[key] = open(joinpath(output_dir, name * suffix * ".csv"), "w")
            println(ios[key], header)
        end
    catch
        close_output_files(ios)
        rethrow()
    end
    ios
end

function close_output_files(ios)
    for io in values(ios)
        isopen(io) && close(io)
    end
end

# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

function run_config(config; scenario_range=nothing, check_only=false, label="")
    version = pkgversion(StochasticGene)
    version !== nothing && version >= v"2.0.4" || error(
        "StochasticGene >= 2.0.4 is required for hidden observed-unit filtering and " *
        "multi-reaction coupling; loaded $(something(version, "unknown"))")
    settings = simulation_settings(config)
    hidden = hidden_settings(config)
    layout = unit_layout(hidden.nh)
    base = Dict(g => base_gene_rates(config, g) for g in GENES)
    points, sweep_keys = sweep_points(config)
    check_hidden_sweep_keys(hidden, sweep_keys)
    max_scenarios = Int(get(require_table(config, "simulation"), "max_scenarios", 100_000))
    length(points) <= max_scenarios || error(
        "the Cartesian sweep has $(length(points)) scenarios, above max_scenarios=$max_scenarios")
    edges = parse_coupling(config, sweep_keys, hidden)
    scoring = scoring_settings(config)
    range = scenario_range === nothing ? (1:length(points)) : scenario_range
    !isempty(range) && first(range) >= 1 && last(range) <= length(points) ||
        error("scenario range $range is outside 1:$(length(points))")

    # Build and validate every scenario before simulating any of them.
    models = Vector{Any}(undef, length(points))
    for s in eachindex(points)
        models[s] = try
            model_for_point(points[s], hidden, layout, edges, base)
        catch err
            error("scenario $s ($(points[s])): " * sprint(showerror, err))
        end
    end

    mkpath(settings.output_dir)
    suffix = scenario_range === nothing ? "" : "_scenarios_$(first(range))-$(last(range))"
    shared = scenario_range === nothing || 1 in range
    if shared
        write_sweep_points(joinpath(settings.output_dir, "sweep_points.csv"), points, sweep_keys)
        write_design(joinpath(settings.output_dir, "sweep_design.csv"),
                     points, models, sweep_keys, hidden, settings)
    end

    warmups = [warmup_for(models[s], settings) for s in range]
    println(isempty(label) ? "StochasticGene $version" : "$label — StochasticGene $version")
    println("Julia threads: ", Threads.nthreads())
    println("hidden topology: ", hidden.topology, "; ", length(range), " of ", length(points),
            " scenarios x ", settings.nsim, " simulations x ", settings.ntraces, " traces")
    println("burn-in per trace (min): ", fmt(minimum(warmups)), " to ", fmt(maximum(warmups)),
            "; recorded time (min): ", fmt((settings.nframes + 2) * settings.interval))
    println("output: ", settings.output_dir)
    if check_only
        println("--check: configuration valid; no simulations run")
        return settings.output_dir
    end

    if settings.save_baseline_binary && shared
        Random.seed!(settings.seed)
        # Fold 1 everywhere gives gamma 0: uncoupled, with the first point's hidden rates.
        baseline_edges = [merge(edge, (fixed_fold=1.0,)) for edge in edges]
        baseline = model_for_point(points[1], hidden, layout, baseline_edges, base)
        baseline_traces = reporter_traces(simulate_one(baseline, layout, settings), settings)
        write_binary_traces(joinpath(settings.output_dir, "baseline_binary_traces$(suffix).csv"),
                            baseline_traces, settings)
    end

    ios = open_output_files(settings.output_dir, suffix, scoring)
    completed = Threads.Atomic{Int}(0)
    total = length(range) * settings.nsim
    # Thread over (scenario, simulation) pairs so that all threads are used even
    # when nsim is small; results are written in scenario order after each block.
    block_size = max(1, cld(4 * Threads.nthreads(), settings.nsim))
    try
        for block in Iterators.partition(range, block_size)
            jobs = [(scenario, sim) for scenario in block for sim in 1:settings.nsim]
            results = Vector{Any}(undef, length(jobs))
            Threads.@threads :dynamic for j in eachindex(jobs)
                local scenario, sim
                scenario, sim = jobs[j]
                Random.seed!(settings.seed + 1_000_003 * scenario + sim)
                local traces = reporter_traces(simulate_one(models[scenario], layout, settings), settings)
                results[j] = analyze_simulation(traces, settings)
                local done = Threads.atomic_add!(completed, 1) + 1
                (done == total || done % max(1, total ÷ 20) == 0) &&
                    println("completed ", done, "/", total)
            end
            for (i, scenario) in enumerate(block)
                r = results[(i - 1) * settings.nsim .+ (1:settings.nsim)]
                metrics = [x.metrics for x in r]
                xcors = [x.xcor for x in r]
                multitaus = [x.multitau for x in r]
                alls = [x.multitau_all for x in r]
                write_metrics_rows(ios[:metrics], scenario, metrics)
                write_xcor_rows(ios[:xcor], scenario, xcors)
                write_metric_summary_row(ios[:metrics_summary], scenario, metrics)
                write_xcor_summary_rows(ios[:xcor_summary], scenario, xcors)
                write_multitau_rows(ios[:multitau], scenario, multitaus)
                write_multitau_summary_rows(ios[:multitau_summary], scenario, multitaus)
                write_multitau_all_rows(ios[:multitau_all], scenario, alls)
                summary = multitau_summary(alls)
                write_multitau_all_summary_rows(ios[:multitau_all_summary], scenario, summary)
                scoring === nothing || write_score_rows(ios[:scores], scenario, summary, scoring)
            end
            foreach(flush, values(ios))
        end
    finally
        close_output_files(ios)
    end
    println("done")
    settings.output_dir
end

const USAGE = "usage: julia simulate_coupled_sweep.jl CONFIG.toml [FIRST LAST] [--check]"

function main(args)
    check_only = "--check" in args
    rest = filter(!=("--check"), args)
    length(rest) in (1, 3) || error(USAGE)
    config_path = abspath(rest[1])
    isfile(config_path) || error("configuration file does not exist: $config_path")
    scenario_range = length(rest) == 3 ? (parse(Int, rest[2]):parse(Int, rest[3])) : nothing
    run_config(TOML.parsefile(config_path); scenario_range, check_only)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
