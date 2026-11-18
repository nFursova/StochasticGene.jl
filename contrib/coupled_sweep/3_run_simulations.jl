#!/usr/bin/env julia
#=
Step 3. Simulate a list of simulations and compare them with the experimental curves.

    julia 3_run_simulations.jl config.toml --swarm        write <run_dir>/simulations.swarm (simulates nothing)
    julia -t 12 3_run_simulations.jl config.toml 1 12     simulate rows 1 to 12 (this is one swarm line)
    julia -t 12 3_run_simulations.jl config.toml          simulate every row on this machine
    julia 3_run_simulations.jl config.toml --export-traces 57,1155,58429
                                                          write the true promoter states and reporter counts of these
                                                          simulations (one dataset each), in StochasticGene's
                                                          predictedtraces format, for 6_trace_relationships.py

Two kinds of run, recognised from the config file:
    regular    config.toml with [output] run_dir: reads <run_dir>/simulations.csv from 1_generate_models.py
    expanded   config_expanded.toml with [profiles]: reads expanded_models.csv and expanded_simulations.csv from
               1_generate_models_expanded.py, in <run_root>/<profile>. Add --profile NAME (default: full), e.g.
               julia -t 12 3_run_simulations.jl config_expanded.toml --profile limited 1 12

Each run writes two files into <run_dir>/results/:
    curves_<first>-<last>.csv    regular runs: the six correlation curves of every simulation, mean and SE over
                                 replicates, and every replicate (columns rep1, rep2, ...)
    curves_<first>-<last>.f32    expanded runs, to keep millions of simulations small: mean and SE only, as 32-bit
                                 numbers [simulation, row, (mean, se)]; the rows are listed in curve_layout.csv
    metrics_<first>-<last>.csv   one row per simulation: ON fractions, mean reporter counts, NRMSE and CCC,
                                 the base rates (rates while no source is active), the ON fractions and
                                 reporter counts of every replicate (columns on_E_rep1, ...), and clipped (1 if a
                                 rate multiplier would have gone below MIN_FACTOR and was held there)

Theory curves ([simulation] curves = "theory"): the curves are the expected value of the correlation estimator,
computed from the exact joint chain (see theory_curves), instead of the mean over simulated replicates; their SE
is 0 and the curves files have no replicate columns. The replicates are still simulated, for the ON fractions,
reporter means and trace statistics (10 replicates are plenty for those).

Insert steps ([simulation] insert_steps, e.g. ["11", "13", "33"]): the same simulations are also analysed as if
the reporter were inserted at a later position: "13" counts CCAT1 transcripts from position 1 and MYC transcripts
only from position 3 (the first digit is CCAT1, the second MYC). Every variant is written to its own folder,
results/insert11, results/insert13, ... (fast simulator only). Without insert_steps: results/, all positions counted.

Names: E = CCAT1 (enhancer), G = MYC (gene), H = hidden unit. Both genes are the 3301 model:
promoter states 1 <-> 2 <-> 3, initiation from state 3, three reporter positions.
--export-traces takes sim_ids separated by commas, or a CSV file with a column sim_id or start_sim_id (for example
fit_candidates.csv from 7_compare_folds.py). It always uses the fast simulator, whose states are exact, and writes
<run_dir>/exported_traces/sim_<id>/CCAT1_predictedtraces.csv and MYC_predictedtraces.csv (columns Gstate_<i>,
Rstate_<i>, Reporters_<i>; trace lengths as in the run), plus hidden_states.csv for hidden-unit models.

Sections: 1 inputs, 2 model, 3 fast simulator, 4 StochasticGene simulator, 5 analysis, 6 run and swarm.
=#
using LinearAlgebra, SparseArrays, Random, Statistics, TOML, Printf

const RATE_NAMES = ["R12", "R21", "R23", "R32", "Initiate", "Rshift1", "Rshift2", "Eject", "Decay"]
const GENE_TRANSITIONS = ((1, 2), (2, 1), (2, 3), (3, 2))    # promoter transitions, in the order of R12 R21 R23 R32
const OBSERVABLES = ("binary", "reporter")                   # binary = ON/OFF state, reporter = reporter count
const PAIRS = (("EE", 1, 1), ("GG", 2, 2), ("EG", 1, 2))     # autocorrelation of E, of G, and E x G cross-correlation
const MIN_FACTOR = 1e-6       # a coupled rate is never multiplied by less than this (keeps every rate positive)

# ---------------------------------------------------------------------------
# 1. Inputs
# ---------------------------------------------------------------------------

"Rows of a comma-separated file as dictionaries (column name => text)."
function read_csv(path)
    lines = filter(!isempty, readlines(path))
    clean(text) = String(strip(replace(text, "\"" => "")))
    header = clean.(split(lines[1], ','))
    [Dict(zip(header, clean.(split(line, ',')))) for line in lines[2:end]]
end

"The nine fitted rates of one gene from a StochasticGene rates_*.txt file (data rows: ml, mean, median, last)."
function read_rates(path, statistic)
    lines = filter(!isempty, readlines(path))
    tryparse(Float64, first(split(lines[1], ','))) === nothing && popfirst!(lines)    # drop the header line
    row = Dict("ml" => 1, "mean" => 2, "median" => 3, "last" => 4)[statistic]
    parse.(Float64, split(lines[row], ',')[1:9])
end

function read_config(path, profile)
    config = TOML.parsefile(path)
    folder = dirname(abspath(path))
    expanded = haskey(config, "profiles")                  # config_expanded.toml of 1_generate_models_expanded.py
    if expanded
        haskey(config["profiles"], profile) || error("no profile $profile in $path; available: $(join(keys(config["profiles"]), ", "))")
        run_dir = joinpath(folder, config["output"]["run_root"], profile)
        isfile(joinpath(run_dir, "COMPLETE.json")) || error("$run_dir has no COMPLETE.json: run 1_generate_models_expanded.py first")
        replicates = config["profiles"][profile]["replicates"]
    else
        run_dir = joinpath(folder, config["output"]["run_dir"])
        replicates = config["simulation"]["replicates"]
    end
    input(name) = joinpath(folder, config["inputs"][name])
    sim = config["simulation"]
    statistic = config["inputs"]["rate_statistic"]
    fitted = permutedims(hcat(read_rates(input("rates_E"), statistic), read_rates(input("rates_G"), statistic)))
    # Number of frames of every trace: all nframes, or the experimental lengths.
    lengths = isempty(config["inputs"]["trace_lengths"]) ? fill(sim["nframes"], sim["ntraces"]) :
              [parse(Int, row["nframes"]) for row in read_csv(input("trace_lengths"))]
    minimum(lengths) >= 17 || error("every trace needs at least 17 frames for the multi-tau correlation")
    get(sim, "curves", "simulated") in ("simulated", "theory") || error("[simulation] curves must be \"simulated\" or \"theory\"")
    # Experimental curves: (observable, pair, lag in frames) => value.
    targets = Dict{Tuple{String,String,Int},Float64}()
    for row in read_csv(input("experimental_curves"))
        value = tryparse(Float64, row["mean"])
        value === nothing || !isfinite(value) ||
            (targets[(row["observable"], row["pair"], parse(Int, row["lag_frames"]))] = value)
    end
    (; run_dir, expanded, profile, fitted,                                   # 2 x 9 fitted rates: row 1 = E, row 2 = G
       lengths, ntraces=length(lengths), nframes=maximum(lengths), targets,
       simulator=sim["simulator"], replicates, interval=sim["interval_min"],
       rate_matching=get(sim, "rate_matching", "time_average"),
       curves=get(sim, "curves", "simulated"),                          # "simulated" or "theory"
       warmup_timescales=get(sim, "warmup_timescales", 3.0), epsilon=sim["epsilon"],
       insert_steps=String.(get(sim, "insert_steps", [""])),           # "" = every position counted, results/
       max_lag=sim["score_max_lag_min"], seed=sim["seed"], swarm=config["swarm"])
end

"""
Rows first..last of the simulation list, as dictionaries with the columns of a regular simulations.csv:
sim_id, model_id, hidden, EG_state, EG_rates, GE_state, GE_rates, HE_state, HE_rates, HG_state, HG_rates,
fold_EG, fold_GE, fold_HE, fold_HG, h12, h21, h23, h32 (and seed, for expanded runs); with one fold per
regulated rate, columns such as fold_EG_R12 replace the connection's fold.
A state is one state ("3"), several states that each regulate ("1+2"), or "R" (Rsum).
"""
function read_rows(cfg, first, last)
    cfg.expanded || return read_csv(joinpath(cfg.run_dir, "simulations.csv"))[first:last]
    models = Dict(row["architecture_id"] => row for row in read_csv(joinpath(cfg.run_dir, "expanded_models.csv")))
    rows = Dict{String,String}[]
    open(joinpath(cfg.run_dir, "expanded_simulations.csv")) do io
        clean(text) = String(strip(replace(text, "\"" => "")))
        header = clean.(split(readline(io), ','))
        for (i, line) in enumerate(eachline(io))             # line i is sim_id i
            i < first && continue
            i > last && break
            point = Dict(zip(header, clean.(split(line, ','))))
            model = models[point["architecture_id"]]
            state(edge) = model["$(edge)_source_kind"] == "reporter_count" ? "R" : model["$(edge)_states"]
            group = model["group"]
            push!(rows, merge(Dict(name => value for (name, value) in point            # sim_id, folds, hidden rates
                                   if name == "sim_id" || startswith(name, "fold_") || startswith(name, "h")),
                              Dict("model_id" => model["architecture_id"], "seed" => point["point_seed"],
                                   "hidden" => group == "H2" ? "H2" : startswith(group, "H3") ? "H3" : "none"),
                              Dict("$(edge)_state" => state(edge) for edge in ("EG", "GE", "HE", "HG")),
                              Dict("$(edge)_rates" => model["$(edge)_rates"] for edge in ("EG", "GE", "HE", "HG"))))
        end
    end
    rows
end

"""
Write the true promoter states and reporter counts of the simulations `sim_ids` (one dataset each, fast simulator)
in StochasticGene's predictedtraces format, so that 6_trace_relationships.py can analyse them like the data.
"""
function export_traces(cfg, sim_ids)
    for id in sim_ids
        row = read_rows(cfg, id, id)[1]
        model = build_model(row, cfg.fitted, cfg.rate_matching)
        states = simulate_fast_states(model, cfg, Xoshiro(cfg.seed + 7919 * id))
        folder = mkpath(joinpath(cfg.run_dir, "exported_traces", @sprintf("sim_%06d", id)))
        decoded = [sstate(model, s) for s in states]               # (joint promoter state, bits_E, bits_G)
        within(frame, trace) = frame <= cfg.lengths[trace]
        for (gene, name) in ((1, "CCAT1"), (2, "MYC"))
            open(joinpath(folder, "$(name)_predictedtraces.csv"), "w") do io
                println(io, join(["$(column)_$trace" for trace in 1:cfg.ntraces for column in ("Gstate", "Rstate", "Reporters")], ","))
                for frame in 1:cfg.nframes
                    cells = String[]
                    for trace in 1:cfg.ntraces
                        if !within(frame, trace)
                            append!(cells, ("", "", ""))
                            continue
                        end
                        p, bits_E, bits_G = decoded[frame, trace]
                        bits = gene == 1 ? bits_E : bits_G
                        append!(cells, (string(pstate(p)[gene]), "\"[$(bits & 1), $((bits >> 1) & 1), $((bits >> 2) & 1)]\"",
                                        string(count_ones(bits))))
                    end
                    println(io, join(cells, ","))
                end
            end
        end
        if model.nh > 1
            open(joinpath(folder, "hidden_states.csv"), "w") do io
                println(io, join(["Hstate_$trace" for trace in 1:cfg.ntraces], ","))
                for frame in 1:cfg.nframes
                    println(io, join([within(frame, trace) ? string(pstate(decoded[frame, trace][1])[3]) : "" for trace in 1:cfg.ntraces], ","))
                end
            end
        end
        counts = zeros(Int, size(states)..., 2)
        for index in CartesianIndices(states)
            counts[index, 1] = count_ones(decoded[index][2]); counts[index, 2] = count_ones(decoded[index][3])
        end
        promoter, first = promoter_and_first(model, states)
        open(joinpath(folder, "statistics.csv"), "w") do io      # the trace statistics of this dataset
            println(io, "statistic,value")
            for (name, value) in zip(STAT_NAMES, trace_statistics(promoter, first, counts, cfg))
                println(io, name, ",", fmt(value))
            end
        end
        println("sim $id ($(row["model_id"])): wrote $folder")
    end
end

"The sim_ids after --export-traces: a comma-separated list, or a CSV file with a column sim_id or start_sim_id."
function export_ids(text)
    isfile(text) || return parse.(Int, split(text, ","))
    rows = read_csv(text)
    column = haskey(rows[1], "sim_id") ? "sim_id" : "start_sim_id"
    unique(parse(Int, row[column]) for row in rows if !isempty(row[column]))
end

"Number of simulations in the list (for the swarm file)."
count_rows(cfg) = countlines(joinpath(cfg.run_dir, cfg.expanded ? "expanded_simulations.csv" : "simulations.csv")) - 1

# ---------------------------------------------------------------------------
# 2. Model: connections, states and base rates
# ---------------------------------------------------------------------------

"One connection: while unit `source` is in one of `states`, rate `rate` of gene `target` is multiplied by the fold."
struct Edge
    source::Int    # 1 = E, 2 = G, 3 = H
    states::Vector{Int}  # promoter states 1 to 3 (hidden states for H) that regulate, any of them (e.g. [1, 2]);
                   # 4 = "R": StochasticGene's Rsum, every occupied reporter position adds one coupling effect;
                   # [3, 4] = "3+R" with one shared fold (lists made before "3+R" had a reporter fold;
                   # newer lists split "3+R" into a [3] and a [4] connection with their own folds)
    target::Int    # 1 = E, 2 = G
    rate::Int      # index into RATE_NAMES
    gamma::Float64 # fold - 1 of this connection (fold_EG, fold_GE, fold_HE or fold_HG)
end

struct Model
    id::String
    nh::Int                         # hidden states: 1 = no hidden unit, 2 = H2, 3 = H3
    hidden_rates::Vector{Float64}   # h12 h21 (H2) or h12 h21 h23 h32 (H3), per minute
    edges::Vector{Edge}
    base::Matrix{Float64}           # 2 x 9 rates while no source is active: row 1 = E, row 2 = G
    clipped::Bool                   # does the multiplier reach MIN_FACTOR in some state?
end

# A state of the whole system is (p, bits_E, bits_G):
#   p        joint promoter state of (E promoter state, G promoter state, hidden state)
#   bits_E   occupied reporter positions of E: bit 1, 2, 4 = position 1, 2, 3 (so 0 to 7); bits_G likewise for G
pindex(e, g, h) = e + 3 * (g - 1) + 9 * (h - 1)
pstate(p) = (mod(p - 1, 3) + 1, mod(div(p - 1, 3), 3) + 1, div(p - 1, 9) + 1)
npromoter(model) = 9 * model.nh
sindex(model, p, bits_E, bits_G) = p + npromoter(model) * (bits_E + 8 * bits_G)
function sstate(model, s)
    code = div(s - 1, npromoter(model))
    mod(s - 1, npromoter(model)) + 1, code % 8, div(code, 8)
end

"How many coupling effects a connection contributes in this state: 1 if its source is in one of the regulating
promoter (or hidden) states, else 0; for R, the number of occupied reporter positions of the source gene (0 to 3)."
function activity(edge, p, bits_E, bits_G)
    effects = Int(pstate(p)[edge.source] in edge.states)          # in a regulating promoter (or hidden) state
    4 in edge.states && (effects += count_ones(edge.source == 1 ? bits_E : bits_G))   # R: + each occupied position
    effects
end

"""
Multiplier of rate k of gene u in state (p, bits_E, bits_G): 1 + gamma for every active coupling effect onto it,
but never below MIN_FACTOR. Clipping only happens with inhibition from source state R, where every occupied
reporter position adds fold - 1: at fold 0.5 three occupied positions would give 1 - 1.5 < 0.
"""
factor(model, p, bits_E, bits_G, u, k) = max(MIN_FACTOR, unclipped_factor(model.edges, p, bits_E, bits_G, u, k))

function unclipped_factor(edges, p, bits_E, bits_G, u, k)
    f = 1.0
    for edge in edges
        edge.target == u && edge.rate == k && (f += edge.gamma * activity(edge, p, bits_E, bits_G))
    end
    f
end

"Every transition of the system for given base rates: for each state, the states it can jump to and at what rates."
function system_transitions(model, base)
    n = 64 * npromoter(model)
    destinations = [Int[] for _ in 1:n]
    rates = [Float64[] for _ in 1:n]
    for bits_G in 0:7, bits_E in 0:7, p in 1:npromoter(model)
        s = sindex(model, p, bits_E, bits_G)
        e, g, h = pstate(p)
        function add(target, value)
            push!(destinations[s], target)
            push!(rates[s], value)
        end
        rate(u, k) = base[u, k] * factor(model, p, bits_E, bits_G, u, k)
        for (k, (from, to)) in enumerate(GENE_TRANSITIONS)                                 # promoter transitions
            e == from && add(sindex(model, pindex(to, g, h), bits_E, bits_G), rate(1, k))
            g == from && add(sindex(model, pindex(e, to, h), bits_E, bits_G), rate(2, k))
        end
        for (k, (from, to)) in enumerate(GENE_TRANSITIONS[1:length(model.hidden_rates)])   # hidden unit
            h == from && add(sindex(model, pindex(e, g, to), bits_E, bits_G), model.hidden_rates[k])
        end
        for u in 1:2                                                                       # reporter steps of gene u
            bits = u == 1 ? bits_E : bits_G
            with(b) = u == 1 ? sindex(model, p, b, bits_G) : sindex(model, p, bits_E, b)
            (u == 1 ? e : g) == 3 && bits & 1 == 0 && add(with(bits | 1), rate(u, 5))      # initiate into empty position 1
            bits & 1 != 0 && bits & 2 == 0 && add(with(xor(bits, 3)), rate(u, 6))          # move 1 -> 2
            bits & 2 != 0 && bits & 4 == 0 && add(with(xor(bits, 6)), rate(u, 7))          # move 2 -> 3
            bits & 4 != 0 && add(with(xor(bits, 4)), rate(u, 8))                           # eject from position 3
        end
    end
    destinations, rates
end

"Long-run probability of every state: balance equations, with the last one replaced by 'probabilities sum to 1'."
function stationary(destinations, rates)
    n = length(rates)
    I, J, V = Int[], Int[], Float64[]
    for s in 1:n
        for (target, rate) in zip(destinations[s], rates[s])
            target < n && (push!(I, target); push!(J, s); push!(V, rate))
        end
        s < n && (push!(I, s); push!(J, s); push!(V, -sum(rates[s])))
        push!(I, n); push!(J, s); push!(V, 1.0)
    end
    prob = max.(sparse(I, J, V, n, n) \ [zeros(n - 1); 1.0], 0.0)
    prob ./ sum(prob)
end

"Can rate k of gene u act in this state? A promoter transition needs its starting promoter state;
initiation needs promoter state 3 and an empty first reporter position."
function possible(p, bits_E, bits_G, u, k)
    state = pstate(p)[u]
    k <= 4 && return state == GENE_TRANSITIONS[k][1]
    state == 3 && (u == 1 ? bits_E : bits_G) & 1 == 0
end

"""
Base rates that keep every regulated rate's average equal to its fitted value: base = fitted / (1 + gamma * a),
where a is the average number of active coupling effects. `matching` chooses which average:
- "time_average": over all time, as `rate_reference = "time_average"` did in simulate_shared_hidden_sweep.jl.
  Coupling can then still shift a gene's ON fraction.
- "state_average": only over the time when the gene is in the state where that rate can act. This keeps each
  gene's promoter-state occupancies equal to those of the uncoupled fit.
The average depends on the coupling itself, so the base rates are solved with Newton's method.
"""
function base_rates(model, fitted, matching)
    matching in ("time_average", "state_average") || error("rate_matching must be time_average or state_average")
    regulated = unique([(edge.target, edge.rate) for edge in model.edges])
    # log(average / fitted) of every regulated rate, for log base rates x
    function mismatch(x)
        base = copy(fitted)
        for ((u, k), value) in zip(regulated, x)
            base[u, k] = exp(value)
        end
        prob = stationary(system_transitions(model, base)...)
        counts(s, u, k) = matching == "time_average" || possible(sstate(model, s)..., u, k)     # is state s averaged over?
        average(u, k) = base[u, k] * sum(prob[s] * factor(model, sstate(model, s)..., u, k) for s in eachindex(prob) if counts(s, u, k)) /
                        sum(prob[s] for s in eachindex(prob) if counts(s, u, k))
        [log(average(u, k) / fitted[u, k]) for (u, k) in regulated], base
    end
    x = [log(fitted[u, k]) for (u, k) in regulated]
    miss, base = mismatch(x)
    for _ in 1:100
        (isempty(miss) || maximum(abs, miss) < 1e-10) && return base
        jacobian = hcat([(mismatch(x .+ 1e-6 .* (eachindex(x) .== j))[1] .- miss) ./ 1e-6 for j in eachindex(x)]...)
        step = jacobian \ miss
        step .*= min(1.0, 4.0 / maximum(abs, step))
        for _ in 1:30                                   # shorten the step until the mismatch decreases
            sum(abs2, mismatch(x .- step)[1]) < sum(abs2, miss) && break
            step ./= 2
        end
        x .-= step
        miss, base = mismatch(x)
    end
    error("could not match the fitted average rates of $(model.id)")
end

"Build the model of one row of the simulation list (see read_rows)."
function build_model(row, fitted, matching)
    nh = row["hidden"] == "none" ? 1 : row["hidden"] == "H2" ? 2 : 3
    edges = Edge[]
    # The fold of one regulated rate: its own column (fold_EG_R12) if set, else the connection's (fold_EG),
    # else the single "fold" column of lists made before the folds were split per connection.
    function fold(edge, rate)
        value = get(row, "fold_$(edge)_$(rate)", "")
        isempty(value) ? (haskey(row, "fold_$edge") ? row["fold_$edge"] : row["fold"]) : value
    end
    # "3+R": the reporter fold is the same column with "_reporter" appended; lists without it share one fold.
    function reporter_fold(edge, rate)
        value = get(row, "fold_$(edge)_$(rate)_reporter", "")
        isempty(value) ? get(row, "fold_$(edge)_reporter", "") : value
    end
    function connect(edge, source, states, target)
        for name in split(row["$(edge)_rates"], '+'; keepempty=false)
            rate = findfirst(==(name), RATE_NAMES)
            reporter = states == [3, 4] ? reporter_fold(edge, name) : ""
            if isempty(reporter)
                push!(edges, Edge(source, states, target, rate, parse(Float64, fold(edge, name)) - 1))
            else                                        # 3+R with independent folds: two connections onto the rate
                push!(edges, Edge(source, [3], target, rate, parse(Float64, fold(edge, name)) - 1))
                push!(edges, Edge(source, [4], target, rate, parse(Float64, reporter) - 1))
            end
        end
    end
    states(text) = [part == "R" ? 4 : parse(Int, part) for part in split(text, '+')]     # "3", "1+2", "R" or "3+R"
    # Hidden units regulate from fixed states unless the list says otherwise:
    # H2: state 2 regulates both genes; H3: state 1 regulates E, state 3 regulates G.
    hidden_states(edge, default) = isempty(get(row, "$(edge)_state", "")) ? [default] : states(row["$(edge)_state"])
    isempty(row["EG_rates"]) || connect("EG", 1, states(row["EG_state"]), 2)               # E regulates G
    isempty(row["GE_rates"]) || connect("GE", 2, states(row["GE_state"]), 1)               # G regulates E
    isempty(row["HE_rates"]) || connect("HE", 3, hidden_states("HE", nh == 2 ? 2 : 1), 1)
    isempty(row["HG_rates"]) || connect("HG", 3, hidden_states("HG", nh == 2 ? 2 : 3), 2)
    hidden_rates = Float64[parse(Float64, row[name]) for name in ("h12", "h21", "h23", "h32")[1:2*(nh-1)]]
    clipped = any(unclipped_factor(edges, p, bits_E, bits_G, edge.target, edge.rate) < MIN_FACTOR
                  for edge in edges, p in 1:9*nh, bits_E in 0:7, bits_G in 0:7)
    model = Model(row["model_id"], nh, hidden_rates, edges, copy(fitted), clipped)
    model.base .= base_rates(model, fitted, matching)
    model
end

# ---------------------------------------------------------------------------
# 3. Fast simulator
# ---------------------------------------------------------------------------

"""
Uses the list of all states and transitions of the system (9 * nh * 8 * 8 states), starts each trace from
the exact long-run distribution, so no burn-in is needed, and runs a standard Gillespie simulation.
Returns the system state (see sindex) at every frame of every trace, [frame, trace].
"""
function simulate_fast_states(model, cfg, rng)
    destinations, rates = system_transitions(model, model.base)
    start = cumsum(stationary(destinations, rates))         # used to draw the starting state of each trace
    cumulative = cumsum.(rates)                             # used to draw which event happens next
    states = zeros(Int, cfg.nframes, cfg.ntraces)
    for trace in 1:cfg.ntraces
        s = min(searchsortedfirst(start, rand(rng)), length(start))
        next_event = randexp(rng) / cumulative[s][end]
        for frame in 1:cfg.nframes
            while next_event <= (frame - 1) * cfg.interval
                choice = searchsortedfirst(cumulative[s], rand(rng) * cumulative[s][end])
                s = destinations[s][min(choice, end)]
                next_event += randexp(rng) / cumulative[s][end]
            end
            states[frame, trace] = s
        end
    end
    states
end

# ---------------------------------------------------------------------------
# 4. StochasticGene simulator (the original way)
# ---------------------------------------------------------------------------

const NOISE = [0.0, 0.1, 1.0, 0.1]    # required by the call below but never used: traces are reporter counts (col = 3)

"""
StochasticGene.simulate_trace_vector, called as in simulate_shared_hidden_sweep.jl.
Every trace starts in state 1 and runs for a burn-in of warmup_timescales / (slowest rate) before
recording; the slowest rate is the smallest of the fitted gene rates and the hidden rates.
Returns reporter counts [frame, trace, gene].
"""
function simulate_stochasticgene(model, cfg)
    hidden = model.nh > 1
    unit(u) = hidden ? (u == 3 ? 1 : u + 1) : u            # StochasticGene unit order: [H, E, G] or [E, G]
    # StochasticGene starts every unit in state 1 and draws its first reaction times without the coupling,
    # so no unit may start in a state that regulates. Three adjustments, for the StochasticGene call only:
    #   H3: the hidden states are renumbered. Its state 1 is our middle state 2, its state 2 is our state 1
    #       (regulates E), state 3 stays state 3 (regulates G). H2 needs nothing (state 2 regulates).
    #   A connection whose regulating states include state 1 (after renumbering) is written the equivalent
    #       way round: the rate is fold times higher, and every other state divides it by fold. This is exact
    #       only when no other connection regulates that rate (3+R's two connections never include state 1).
    #   R is StochasticGene's Rsum: one connection from each reporter position (source states 4, 5, 6).
    renumber(state) = model.nh == 3 ? (2, 1, 3)[state] : state
    hidden_transitions = model.nh == 3 ? ((2, 1), (1, 2), (1, 3), (3, 1)) : GENE_TRANSITIONS[1:2]    # order of h12 h21 h23 h32
    base = copy(model.base)
    connections = NTuple{4,Int}[]       # (source unit, source state, target unit, target rate)
    gammas = Float64[]                  # one coupling strength per connection
    for edge in model.edges
        states = edge.source == 3 ? renumber.(edge.states) :
                 vcat(filter(<(4), edge.states), 4 in edge.states ? [4, 5, 6] : Int[])     # R: one connection per position
        gamma = edge.gamma
        if 1 in states                                  # the unit starts in a regulating state: write it the other way round
            count(other -> (other.target, other.rate) == (edge.target, edge.rate), model.edges) == 1 ||
                error("$(model.id): $(RATE_NAMES[edge.rate]) is regulated by more than one connection")
            base[edge.target, edge.rate] *= 1 + gamma
            states = setdiff(1:(edge.source == 3 ? model.nh : 3), states)
            gamma = 1 / (1 + gamma) - 1
        end
        for state in states
            push!(connections, (unit(edge.source), state, unit(edge.target), edge.rate))
            push!(gammas, gamma)
        end
    end
    rates = Float64[]
    hidden && append!(rates, model.hidden_rates, (0.0, 0.0))    # two unused RNA slots of the hidden unit
    for u in 1:2
        append!(rates, base[u, :], NOISE)
    end
    append!(rates, gammas)
    transitions = hidden ? (hidden_transitions, GENE_TRANSITIONS, GENE_TRANSITIONS) : (GENE_TRANSITIONS, GENE_TRANSITIONS)
    G = hidden ? (model.nh, 3, 3) : (3, 3)
    R = hidden ? (0, 3, 3) : (3, 3)
    slowest = minimum(vcat(model.hidden_rates, vec(cfg.fitted[:, 1:8])))
    traces = StochasticGene.simulate_trace_vector(
        rates, transitions, G, R, map(_ -> 0, G), map(_ -> 1, G),
        (Tuple(1:length(G)), connections), cfg.interval, (cfg.nframes + 2) * cfg.interval, cfg.ntraces;
        onstates=[Int[] for _ in G], probfn=StochasticGene.prob_Gaussian,
        noiseparams=hidden ? [0, 4, 4] : [4, 4], observed_units=hidden ? [2, 3] : [1, 2], col=3,
        totalsteps=100_000_000, warmuptime=cfg.warmup_timescales / slowest, equilibration_timescales=1.0)
    counts = zeros(Int, cfg.nframes, cfg.ntraces, 2)
    for (trace, values) in enumerate(traces)
        counts[:, trace, :] = round.(Int, values[1:cfg.nframes, :])
    end
    counts
end

# ---------------------------------------------------------------------------
# 5. Analysis
# ---------------------------------------------------------------------------

# FCS multi-tau correlation and Coulon average, copied from simulate_shared_hidden_sweep.jl
# (ported there from fcs_multitau_xcor_final.R). Only the background-correction branch is removed:
# the traces here are reporter counts, and experimental.R also runs with bg_correct = FALSE.
idl_round(x) = sign(x) * floor(abs(x) + 0.5)

struct MTCurve
    lag::Vector{Float64}
    Gn::Vector{Float64}
    Md::Vector{Float64}
    Me::Vector{Float64}
    n::Vector{Float64}
end

function multi_tau_corr(left::Vector{Float64}, right::Vector{Float64}; M::Int=4, lag_sign::Int=1)
    n = length(left)
    selected = n - 1
    N = Int(idl_round(log2((selected + 1) / M)))
    nchan = M * N
    L = left[1:selected]
    Rv = right[1:selected]
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
            Gn[j] = (np[j] * Go[j]) / (Md[j] * Me[j]) - 1
        end
    end
    MTCurve(lags[1:used], Gn, Md[1:used], Me[1:used], np[1:used])
end

function coulon_multitau_mean(curves::Vector{MTCurve})
    common = copy(first(curves).lag)
    for curve in Iterators.drop(curves, 1)
        common = intersect(common, curve.lag)
    end
    sort!(common)
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
    [(; lag_frames=common[j], Gn=Gn[j]) for j in eachindex(common)]
end

"The six curves of one dataset, in a fixed order: (observable, pair, lag in frames, value)."
function correlation_curves(counts, cfg)
    rows = NamedTuple[]
    for observable in OBSERVABLES, (pair, i, j) in PAIRS
        curves = MTCurve[]
        for trace in axes(counts, 2)
            frames = 1:cfg.lengths[trace]
            x = Float64.(counts[frames, trace, i])
            y = Float64.(counts[frames, trace, j])
            if observable == "binary"
                x = Float64.(x .> 0)
                y = Float64.(y .> 0)
            end
            x .+= cfg.epsilon
            y .+= cfg.epsilon
            forward = multi_tau_corr(x, y; lag_sign=1)       # positive lag: the second gene follows the first
            backward = multi_tau_corr(y, x; lag_sign=-1)
            push!(curves, MTCurve(vcat(forward.lag, backward.lag), vcat(forward.Gn, backward.Gn),
                                  vcat(forward.Md, backward.Md), vcat(forward.Me, backward.Me),
                                  vcat(forward.n, backward.n)))
        end
        for row in coulon_multitau_mean(curves)
            push!(rows, (; observable, pair, lag=round(Int, row.lag_frames), value=row.Gn))
        end
    end
    rows
end

"""
Theoretical curves ([simulation] curves = "theory"): the expected value of the estimator above (multi_tau_corr and
coulon_multitau_mean) for this model, from the exact moments m(d) = E[x(t) y(t + d)] of the joint chain (d frames
apart; x and y are the observable of the two genes, plus epsilon, as in correlation_curves), instead of the mean
over simulated replicates.

For one trace the estimator at a lag is close to a covariance centred on the trace's own means,
    (mean of the products x y at that lag)  -  (mean of x) (mean of y),
and the curve is its average over traces divided by (long-run mean x) (long-run mean y). The expectation of this is a
weighted sum of m(d) (see expected_multitau); the second term is what centring each trace on its own mean subtracts.
The approximation: the ratios of sums inside the estimator are replaced by the ratios of their expected values
(exact for long traces and many traces). The lags are those of the simulated curves.
Returns the curves of every insert-step variant, each ordered as correlation_curves.
"""
function theory_curves(model, cfg)
    destinations, rates = system_transitions(model, model.base)
    n = length(rates)
    Q = zeros(n, n)
    for s in 1:n, (target, rate) in zip(destinations[s], rates[s])
        Q[s, target] += rate
        Q[s, s] -= rate
    end
    Pt = permutedims(exp(Q .* cfg.interval))            # transposed one-frame transition matrix
    prob = stationary(destinations, rates)
    states = [sstate(model, s) for s in 1:n]
    dmax = cfg.nframes - 1
    map(cfg.insert_steps) do variant
        first_E, first_G = variant == "" ? (1, 1) : parse.(Int, collect(variant))
        counts = (Float64[count_ones(b >> (first_E - 1)) for (_, b, _) in states],
                  Float64[count_ones(b >> (first_G - 1)) for (_, _, b) in states])
        rows = NamedTuple[]
        for observable in OBSERVABLES
            f = [(observable == "binary" ? Float64.(c .> 0) : c) .+ cfg.epsilon for c in counts]
            means = [dot(prob, g) for g in f]
            m = [zeros(dmax + 1, 2) for _ in 1:2]           # m[i][d + 1, j] = E[f_i(t) f_j(t + d)]
            for i in 1:2
                w = prob .* f[i]
                for d in 0:dmax
                    m[i][d + 1, 1] = dot(w, f[1]); m[i][d + 1, 2] = dot(w, f[2])
                    w = Pt * w
                end
            end
            for (pair, i, j) in PAIRS
                moment(d) = d >= 0 ? m[i][d + 1, j] : m[j][-d + 1, i]
                for (lag, value) in expected_multitau(moment, cfg.lengths; window_from=i != j ? 0 : cfg.max_lag / cfg.interval)
                    push!(rows, (; observable, pair, lag, value=value / (means[i] * means[j])))
                end
            end
        end
        rows
    end
end

"""(lag in frames, expected trace-centred covariance) of the multi-tau estimator, averaged over traces of the given
lengths, for moments moment(d) = E[x(t) y(t + d)]; lags as multi_tau_corr makes them (lags present in every trace,
averaged over the traces where the lag has products), sorted. At a lag of k bins of B frames with p products:
    (average of moment over the bin's frame offsets)  -  centring
centring = E[(mean of x) (mean of y)], over the n - 2 frames of the estimator's trace means, or, at lags beyond
window_from frames, over the two windows of the p products. Checked against 100-replicate simulated means of 137
traces of 216 frames: window centring matches the cross-correlation at every lag (window_from = 0), trace centring
the autocorrelations out to 60 min, both within about 1 SE; the autocorrelations use window centring beyond the
scoring range (score_max_lag_min), where it is closer, and are less exact from 60 min on (the estimator is a ratio
of sums, which these expectations only approximate)."""
function expected_multitau(moment, lengths; window_from::Real, M::Int=4)
    totals, counts, lag_lists = Dict{Int,Float64}(), Dict{Int,Int}(), Vector{Int}[]
    # E[x(bin i) y(bin i + k)] / B^2 for bins of B frames: the average of moment over the frame offsets
    binned(k, B) = sum((B - abs(o)) * moment(k * B + o) for o in (1 - B):(B - 1)) / B^2
    for (n, copies) in [(n, count(==(n), lengths)) for n in unique(lengths)]
        selected = n - 1
        n1 = n - 2
        trace_centring = sum((n1 - abs(d - 1)) * moment(d) for d in (2 - n1):n1) / n1^2
        lags = Int[]
        function add(lag, B, k, p)                      # lag in frames, bin size, lag in bins, products
            push!(lags, lag)
            p >= 1 || return
            centring = abs(lag) > window_from ? sum((p - abs(d)) * binned(k + d, B) for d in (1 - p):(p - 1)) / p^2 : trace_centring
            totals[lag] = get(totals, lag, 0.0) + copies * (binned(k, B) - centring)
            counts[lag] = get(counts, lag, 0) + copies
        end
        for sign in (1, -1)
            for lag in 1:(2 * M)                                       # linear part: lags 1 .. 2M frames
                add(sign * lag, 1, sign * lag, selected - lag)
            end
            points = selected
            for octave in 1:(Int(idl_round(log2((selected + 1) / M))) - 2)   # binned part: bins of 2^octave frames
                points = div(points - points % 2, 2)
                for i in 1:M
                    add(sign * (M + i) * 2^octave, 2^octave, sign * (M + i), points - (M + i))
                end
            end
        end
        push!(lag_lists, lags)
    end
    common = sort(intersect(lag_lists...))
    [(lag, haskey(counts, lag) ? totals[lag] / counts[lag] : NaN) for lag in common]
end

"ON fraction and mean reporter count of each gene, over all recorded frames."
function on_and_reporters(counts, cfg)
    total(f, gene) = sum(sum(f, view(counts, 1:cfg.lengths[trace], trace, gene)) for trace in axes(counts, 2)) /
                     sum(cfg.lengths)
    [total(>(0), 1), total(>(0), 2), total(identity, 1), total(identity, 2)]    # on_E, on_G, reporters_E, reporters_G
end

"NRMSE (0 = identical) and Lin's concordance correlation coefficient (1 = identical) of simulated x against experimental y."
function nrmse_ccc(x, y)
    length(x) < 2 && return NaN, NaN
    mx, my = mean(x), mean(y)
    nrmse = sqrt(sum(abs2, x .- y) / sum(abs2, y))
    ccc = 2 * mean((x .- mx) .* (y .- my)) / (mean(abs2, x .- mx) + mean(abs2, y .- my) + (mx - my)^2)
    nrmse, ccc
end

# ---------------------------------------------------------------------------
# 6. Run and swarm
# ---------------------------------------------------------------------------

"One replicate dataset of one simulation: its curves and its ON fractions / reporter means."
function simulate_replicate(row, model, replicate, cfg)
    seed = haskey(row, "seed") ? parse(UInt64, row["seed"]) + UInt64(replicate) :      # expanded runs: the generator's point_seed
           cfg.seed + 1000 * parse(Int, row["sim_id"]) + replicate
    if cfg.simulator != "fast"
        cfg.insert_steps == [""] || error("insert_steps needs simulator = \"fast\"")
        Random.seed!(seed)
        counts = simulate_stochasticgene(model, cfg)
        return [(; curves=cfg.curves == "theory" ? nothing : correlation_curves(counts, cfg),
                   metrics=on_and_reporters(counts, cfg), stats=trace_statistics(nothing, nothing, counts, cfg))]
    end
    states = simulate_fast_states(model, cfg, Xoshiro(seed))
    promoter, first = promoter_and_first(model, states)
    map(cfg.insert_steps) do variant
        first_E, first_G = variant == "" ? (1, 1) : parse.(Int, collect(variant))
        counts = zeros(Int, cfg.nframes, cfg.ntraces, 2)
        for index in CartesianIndices(states)            # count the occupied positions from the insert position on
            _, bits_E, bits_G = sstate(model, states[index])
            counts[index, 1] = count_ones(bits_E >> (first_E - 1))
            counts[index, 2] = count_ones(bits_G >> (first_G - 1))
        end
        (; curves=cfg.curves == "theory" ? nothing : correlation_curves(counts, cfg),
           metrics=on_and_reporters(counts, cfg), stats=trace_statistics(promoter, first, counts, cfg))
    end
end

# ---------------------------------------------------------------------------
# Trace statistics, computed exactly as paper_figures.R computes them from the decoded experimental traces
# ---------------------------------------------------------------------------
const STAT_SOURCES = ("1", "2", "3", "R", "3R")          # source conditions ("3R" = 3+R)
const STAT_RATES = ("R12", "R21", "R23", "R32", "Initiate")
const STAT_MIN_FRAMES = 5      # per track: frames needed in and out of the condition
const STAT_MIN_TRACKS = 5      # tracks needed for a median
const STAT_NAMES = vcat(["mfold_$(d)_$(c)_$(r)" for d in ("EG", "GE") for c in STAT_SOURCES for r in STAT_RATES],
                        ["cob_onset_G_given_E", "cob_onset_E_given_G", "cob_coON"])

"""
The median over tracks of each trace statistic of one dataset (STAT_NAMES):
  mfold_<EG|GE>_<condition>_<rate>: the fold measured like in the data: the rate of a transition of the target gene
      (G for EG) is base x (1 + (fold - 1) x n), n = number of coupling effects of the source condition at that frame
      (states 1, 2, 3: 1 while in it; R: the number of reporters; 3R: 1 in state 3 plus the number of reporters);
      fold = 1 + ((k1 + 0.5) T0 / (k0 + 0.5) - T1) / N1 (k = events, T = frames at risk, 1 = n > 0, 0 = n = 0,
      N1 = sum of n), at least 0.001; NA with fewer than STAT_MIN_FRAMES frames on a side or no event.
      Needs promoter states (NaN with the StochasticGene simulator).
  cob_onset_G_given_E: G burst-onset rate (reporters 0 -> above 0) while E has reporters / while it has none
      ((k1 + 0.5) / T1) / ((k0 + 0.5) / T0); cob_onset_E_given_G likewise;
  cob_coON: P(both ON) / (P(E ON) P(G ON)) = (n_both + 0.5) L / ((n_E + 0.5)(n_G + 0.5)); NA if a gene is never ON.
promoter[frame, trace, gene] (or nothing), first[frame, trace, gene] = first reporter position occupied,
counts[frame, trace, gene] = reporter counts.
"""
function trace_statistics(promoter, first, counts, cfg)
    values = [Float64[] for _ in STAT_NAMES]
    for trace in axes(counts, 2)
        L = cfg.lengths[trace]
        if promoter !== nothing
            for (d, (src, tgt)) in enumerate(((1, 2), (2, 1)))
                T1 = zeros(5, 5); k1 = zeros(5, 5); N1 = zeros(5, 5); T = zeros(5); K = zeros(5)
                for t in 1:(L - 1)
                    ps, rs = promoter[t, trace, src], counts[t, trace, src]
                    n = (Int(ps == 1), Int(ps == 2), Int(ps == 3), rs, Int(ps == 3) + rs)
                    pt, pn = promoter[t, trace, tgt], promoter[t + 1, trace, tgt]
                    risk = (pt == 1, pt == 2, pt == 2, pt == 3, pt == 3 && !first[t, trace, tgt])
                    event = (pn == 2, pn == 1, pn == 3, pn == 2, first[t + 1, trace, tgt])
                    for r in 1:5
                        risk[r] || continue
                        e = Int(event[r])
                        T[r] += 1; K[r] += e
                        for c in 1:5
                            n[c] > 0 || continue
                            T1[c, r] += 1; k1[c, r] += e; N1[c, r] += n[c]
                        end
                    end
                end
                for c in 1:5, r in 1:5
                    T0, k0 = T[r] - T1[c, r], K[r] - k1[c, r]
                    (T1[c, r] < STAT_MIN_FRAMES || T0 < STAT_MIN_FRAMES || k1[c, r] + k0 == 0) && continue
                    fold = 1 + ((k1[c, r] + 0.5) * T0 / (k0 + 0.5) - T1[c, r]) / N1[c, r]
                    push!(values[(d - 1) * 25 + (c - 1) * 5 + r], max(fold, 1e-3))
                end
            end
        end
        for (k, (cond, tgt)) in enumerate(((1, 2), (2, 1)))          # onsets of tgt while cond has reporters or not
            T1 = k1 = T0 = k0 = 0
            for t in 1:(L - 1)
                counts[t, trace, tgt] == 0 || continue
                e = Int(counts[t + 1, trace, tgt] > 0)
                if counts[t, trace, cond] > 0
                    T1 += 1; k1 += e
                else
                    T0 += 1; k0 += e
                end
            end
            (T1 < STAT_MIN_FRAMES || T0 < STAT_MIN_FRAMES || k1 + k0 == 0) && continue
            push!(values[50 + k], ((k1 + 0.5) / T1) / ((k0 + 0.5) / T0))
        end
        onE = count(>(0), view(counts, 1:L, trace, 1)); onG = count(>(0), view(counts, 1:L, trace, 2))
        both = count(t -> counts[t, trace, 1] > 0 && counts[t, trace, 2] > 0, 1:L)
        (onE > 0 && onG > 0) && push!(values[53], (both + 0.5) * L / ((onE + 0.5) * (onG + 0.5)))
    end
    [length(v) >= STAT_MIN_TRACKS ? median(v) : NaN for v in values]
end

"Promoter states and first-position occupancy [frame, trace, gene] of fast-simulator states."
function promoter_and_first(model, states)
    promoter = zeros(Int, size(states)..., 2); first = falses(size(states)..., 2)
    for index in CartesianIndices(states)
        p, bits_E, bits_G = sstate(model, states[index])
        e, g, _ = pstate(p)
        promoter[index, 1] = e; promoter[index, 2] = g
        first[index, 1] = isodd(bits_E); first[index, 2] = isodd(bits_G)
    end
    promoter, first
end

fmt(x) = isfinite(x) ? @sprintf("%.8g", x) : "NA"
mean_se(values) = (mean(values), length(values) > 1 ? std(values) / sqrt(length(values)) : NaN)

"Simulate rows first_id to last_id of simulations.csv and write their curves and metrics."
function run_simulations(rows, cfg, first_id, last_id)
    models = Vector{Model}(undef, length(rows))
    Threads.@threads :dynamic for i in eachindex(rows)
        models[i] = build_model(rows[i], cfg.fitted, cfg.rate_matching)
    end
    nclipped = count(model -> model.clipped, models)
    nclipped > 0 && println("$nclipped of $(length(models)) simulations have a rate multiplier clipped at $MIN_FACTOR")
    jobs = [(i, replicate) for i in eachindex(rows) for replicate in 1:cfg.replicates]
    results = Vector{Any}(undef, length(jobs))
    Threads.@threads :dynamic for j in eachindex(jobs)
        results[j] = simulate_replicate(rows[jobs[j][1]], models[jobs[j][1]], jobs[j][2], cfg)
    end

    theory = nothing
    if cfg.curves == "theory"
        theory = Vector{Any}(undef, length(rows))
        Threads.@threads :dynamic for i in eachindex(rows)
            theory[i] = theory_curves(models[i], cfg)
        end
    end
    for (v, variant) in enumerate(cfg.insert_steps)
        folder = mkpath(joinpath(cfg.run_dir, "results", isempty(variant) ? "" : "insert$variant"))
        write_results(rows, models, [result[v] for result in results], cfg, folder, first_id, last_id,
                      theory === nothing ? nothing : [t[v] for t in theory])
    end
end

"Write the curves and metrics of one insert-step variant into `folder`."
function write_results(rows, models, results, cfg, folder, first_id, last_id, theory=nothing)
    stem = @sprintf("%06d-%06d", first_id, last_id)
    curves_path = joinpath(folder, cfg.expanded ? "curves_$stem.f32" : "curves_$stem.csv")
    metrics_path = joinpath(folder, "metrics_$stem.csv")
    score_names = ["$(name)_$(observable)_$(pair)" for observable in OBSERVABLES for (pair, _, _) in PAIRS
                   for name in ("nrmse", "ccc")]
    base_names = ["base_$(gene)_$(RATE_NAMES[k])" for gene in ("E", "G") for k in 1:5]      # R12 R21 R23 R32 Initiate
    metric_names = ("on_E", "on_G", "reporters_E", "reporters_G")
    replicate_names = ["$(name)_rep$(r)" for name in metric_names for r in 1:cfg.replicates]
    if cfg.expanded                                     # the rows of every .f32 curves file, written by every job
        layout_path = joinpath(folder, "curve_layout.csv")
        open(layout_path * ".tmp.$(getpid())", "w") do io
            println(io, "observable,pair,lag_frames,lag_min")
            for bin in (theory === nothing ? results[1].curves : theory[1])
                println(io, join((bin.observable, bin.pair, bin.lag, fmt(bin.lag * cfg.interval)), ','))
            end
        end
        mv(layout_path * ".tmp.$(getpid())", layout_path; force=true)
    end
    # Write to temporary files first, so that a killed job never leaves a half-written result.
    open(curves_path * ".tmp", "w") do curves_io
        open(metrics_path * ".tmp", "w") do metrics_io
            reps = theory === nothing ? ["rep$r" for r in 1:cfg.replicates] : String[]      # theory: no replicate columns
            cfg.expanded || println(curves_io, join(vcat(["sim_id", "observable", "pair", "lag_frames", "lag_min", "mean", "se"], reps), ','))
            println(metrics_io, "sim_id,model_id,on_E,on_E_se,on_G,on_G_se,reporters_E,reporters_E_se," *
                                "reporters_G,reporters_G_se," * join(vcat(score_names, base_names, replicate_names, "clipped",
                                                                          STAT_NAMES), ','))
            for (i, row) in enumerate(rows)
                replicates = results[(i - 1) * cfg.replicates .+ (1:cfg.replicates)]
                if theory === nothing
                    bins = replicates[1].curves                               # same lag bins in every replicate
                    values = hcat([[bin.value for bin in r.curves] for r in replicates]...)
                    means = [mean_se(values[b, :]) for b in eachindex(bins)]
                else                                                          # the expected curve: SE 0
                    bins = theory[i]
                    values = zeros(length(bins), 0)
                    means = [(bin.value, 0.0) for bin in bins]
                end
                if cfg.expanded                                                # mean and SE as 32-bit numbers
                    write(curves_io, Float32[x for b in eachindex(bins) for x in means[b]])
                else
                    for (b, bin) in enumerate(bins)                            # mean, SE, then every replicate
                        println(curves_io, join((row["sim_id"], bin.observable, bin.pair, bin.lag, fmt(bin.lag * cfg.interval),
                                                 fmt.(means[b])..., fmt.(values[b, :])...), ','))
                    end
                end
                cells = String[row["sim_id"], row["model_id"]]
                for k in 1:4                                                  # on_E, on_G, reporters_E, reporters_G
                    append!(cells, fmt.(mean_se([r.metrics[k] for r in replicates])))
                end
                for observable in OBSERVABLES, (pair, _, _) in PAIRS          # score the mean curve
                    use = [b for b in eachindex(bins) if bins[b].observable == observable && bins[b].pair == pair &&
                           abs(bins[b].lag * cfg.interval) <= cfg.max_lag && isfinite(means[b][1]) &&
                           haskey(cfg.targets, (observable, pair, bins[b].lag))]
                    append!(cells, fmt.(nrmse_ccc([means[b][1] for b in use],
                                                  [cfg.targets[(observable, pair, bins[b].lag)] for b in use])))
                end
                append!(cells, fmt.(models[i].base[1, 1:5]), fmt.(models[i].base[2, 1:5]))
                for k in 1:4                                                  # the same four metrics, every replicate
                    append!(cells, fmt.([r.metrics[k] for r in replicates]))
                end
                push!(cells, models[i].clipped ? "1" : "0")
                for k in eachindex(STAT_NAMES)                                # trace statistics: mean over replicates
                    finite = filter(isfinite, [r.stats[k] for r in replicates])
                    push!(cells, fmt(isempty(finite) ? NaN : mean(finite)))
                end
                println(metrics_io, join(cells, ','))
            end
        end
    end
    mv(curves_path * ".tmp", curves_path; force=true)
    mv(metrics_path * ".tmp", metrics_path; force=true)
    println("simulations $first_id to $last_id done: ", metrics_path)
end

"Write the swarm file: one line per job of [swarm] simulations_per_job simulations."
function write_swarm(config_path, cfg, nrows)
    s = cfg.swarm
    project = isempty(s["julia_project"]) ? "" : " --project=$(s["julia_project"])"
    profile = cfg.expanded ? " --profile $(cfg.profile)" : ""
    command = "$(s["julia"]) --threads=$(s["threads"])$project $(abspath(@__FILE__)) $(abspath(config_path))$profile"
    per_job = s["simulations_per_job"]
    lines = ["$command $start $(min(start + per_job - 1, nrows))" for start in 1:per_job:nrows]
    path = joinpath(cfg.run_dir, "simulations.swarm")
    write(path, join(lines, "\n") * "\n")
    println("$(length(lines)) swarm lines of up to $per_job simulations x $(cfg.replicates) replicates. Submit with:")
    println("swarm -f $path -t $(s["threads"]) -g $(s["memory_gb"]) --time $(s["walltime"]) ",
            "--module $(s["modules"]) --logdir $(joinpath(cfg.run_dir, "swarm_logs"))")
end

# ---------------------------------------------------------------------------

const USAGE = "usage: julia 3_run_simulations.jl CONFIG [--profile NAME] [--swarm | FIRST LAST | --export-traces IDS]"
const PROFILE_AT = findfirst(==("--profile"), ARGS)
const EXPORT_AT = findfirst(==("--export-traces"), ARGS)
const VALUE_POSITIONS = vcat([[at, at + 1] for at in (PROFILE_AT, EXPORT_AT) if at !== nothing]...)
const POSITIONAL = [a for (i, a) in enumerate(ARGS) if !(i in VALUE_POSITIONS)]
length(POSITIONAL) in (1, 2, 3) || error(USAGE)
const CONFIG = read_config(POSITIONAL[1], PROFILE_AT === nothing ? "full" : ARGS[PROFILE_AT + 1])
if CONFIG.simulator == "stochasticgene" && !("--swarm" in ARGS) && EXPORT_AT === nothing
    using StochasticGene
    # StochasticGene stops when a coupled rate's multiplier (1 + sum of the active gammas) is not positive.
    # Replace that one internal function, for this run only, so that the multiplier is clipped at MIN_FACTOR,
    # exactly as in the fast simulator. The installed package is not changed.
    isdefined(StochasticGene, :_coupling_factor_for_target) && isdefined(StochasticGene, :_source_occupied) || error(
        "this StochasticGene version does not have the internal function that this script replaces to clip multipliers")
    @eval StochasticGene function _coupling_factor_for_target(state, target::Int, coupling_strength,
            connections, G, R, S; target_transition::Union{Nothing,Int}=nothing,
            override_unit::Int=0, override_state=nothing)
        factor = 1.0
        for (k, (source, source_state, connection_target, connection_transition)) in enumerate(connections)
            connection_target == target || continue
            isnothing(target_transition) || connection_transition == target_transition || continue
            source_vector = source == override_unit ? override_state : state[source, 1]
            occupied = _source_occupied(findall(!iszero, vec(source_vector)), [source_state],
                                        G[source], R[source], S[source] > 0 ? 3 : 2)
            occupied && (factor += coupling_strength[k])
        end
        return max(factor, $MIN_FACTOR)
    end
end
if EXPORT_AT !== nothing
    export_traces(CONFIG, export_ids(ARGS[EXPORT_AT + 1]))
elseif "--swarm" in ARGS
    write_swarm(POSITIONAL[1], CONFIG, count_rows(CONFIG))
else
    first_row, last_row = length(POSITIONAL) == 3 ? parse.(Int, POSITIONAL[2:3]) : (1, count_rows(CONFIG))
    run_simulations(read_rows(CONFIG, first_row, last_row), CONFIG, first_row, last_row)
end
