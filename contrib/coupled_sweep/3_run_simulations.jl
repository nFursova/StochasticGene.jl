#!/usr/bin/env julia
#=
Step 3. Simulate the rows of simulations.csv and compare them with the experimental curves.

    julia 3_run_simulations.jl config.toml --swarm        write <run_dir>/simulations.swarm (simulates nothing)
    julia -t 12 3_run_simulations.jl config.toml 1 12     simulate rows 1 to 12 (this is one swarm line)
    julia -t 12 3_run_simulations.jl config.toml          simulate every row on this machine

Each run writes two files into <run_dir>/results/:
    curves_<first>-<last>.csv    the six correlation curves of every simulation: mean and SE over replicates,
                                 and every replicate (columns rep1, rep2, ...)
    metrics_<first>-<last>.csv   one row per simulation: ON fractions, mean reporter counts, NRMSE and CCC,
                                 the base rates (rates while no source is active), the ON fractions and
                                 reporter counts of every replicate (columns on_E_rep1, ...), and clipped (1 if a
                                 rate multiplier would have gone below MIN_FACTOR and was held there)

Names: E = CCAT1 (enhancer), G = MYC (gene), H = hidden unit. Both genes are the 3301 model:
promoter states 1 <-> 2 <-> 3, initiation from state 3, three reporter positions.
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

function read_config(path)
    config = TOML.parsefile(path)
    folder = dirname(abspath(path))
    input(name) = joinpath(folder, config["inputs"][name])
    sim = config["simulation"]
    statistic = config["inputs"]["rate_statistic"]
    fitted = permutedims(hcat(read_rates(input("rates_E"), statistic), read_rates(input("rates_G"), statistic)))
    # Number of frames of every trace: all nframes, or the experimental lengths.
    lengths = isempty(config["inputs"]["trace_lengths"]) ? fill(sim["nframes"], sim["ntraces"]) :
              [parse(Int, row["nframes"]) for row in read_csv(input("trace_lengths"))]
    minimum(lengths) >= 17 || error("every trace needs at least 17 frames for the multi-tau correlation")
    # Experimental curves: (observable, pair, lag in frames) => value.
    targets = Dict{Tuple{String,String,Int},Float64}()
    for row in read_csv(input("experimental_curves"))
        value = tryparse(Float64, row["mean"])
        value === nothing || !isfinite(value) ||
            (targets[(row["observable"], row["pair"], parse(Int, row["lag_frames"]))] = value)
    end
    (; run_dir=joinpath(folder, config["output"]["run_dir"]),
       fitted,                                   # 2 x 9 fitted rates: row 1 = E, row 2 = G
       lengths, ntraces=length(lengths), nframes=maximum(lengths), targets,
       simulator=sim["simulator"], replicates=sim["replicates"], interval=sim["interval_min"],
       rate_matching=get(sim, "rate_matching", "time_average"),
       warmup_timescales=sim["warmup_timescales"], epsilon=sim["epsilon"],
       max_lag=sim["score_max_lag_min"], seed=sim["seed"], swarm=config["swarm"])
end

# ---------------------------------------------------------------------------
# 2. Model: connections, states and base rates
# ---------------------------------------------------------------------------

"One connection: while unit `source` is in `state`, rate `rate` of gene `target` is multiplied by the fold."
struct Edge
    source::Int    # 1 = E, 2 = G, 3 = H
    state::Int     # promoter state 1 to 3 (hidden state for H), or 4 = "R": StochasticGene's Rsum,
                   # every occupied reporter position of the source gene adds one coupling effect
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

"How many coupling effects a connection contributes in this state: 1 if its source is in the regulating promoter
(or hidden) state, else 0; for R, the number of occupied reporter positions of the source gene (0 to 3)."
function activity(edge, p, bits_E, bits_G)
    edge.state == 4 && return count_ones(edge.source == 1 ? bits_E : bits_G)
    Int(pstate(p)[edge.source] == edge.state)
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

"Build the model of one simulations.csv row."
function build_model(row, fitted, matching)
    nh = row["hidden"] == "none" ? 1 : row["hidden"] == "H2" ? 2 : 3
    edges = Edge[]
    function connect(column, source, state, target, fold)        # fold: the simulations.csv column with this connection's fold
        for name in split(row[column], '+'; keepempty=false)
            push!(edges, Edge(source, state, target, findfirst(==(name), RATE_NAMES), parse(Float64, row[fold]) - 1))
        end
    end
    gene_state(text) = text == "R" ? 4 : parse(Int, text)                               # EG_state / GE_state: 1, 2, 3 or R (Rsum)
    isempty(row["EG_rates"]) || connect("EG_rates", 1, gene_state(row["EG_state"]), 2, "fold_EG")   # E regulates G
    isempty(row["GE_rates"]) || connect("GE_rates", 2, gene_state(row["GE_state"]), 1, "fold_GE")   # G regulates E
    connect("HE_rates", 3, nh == 2 ? 2 : 1, 1, "fold_HE")        # H2: hidden state 2 regulates both genes
    connect("HG_rates", 3, nh == 2 ? 2 : 3, 2, "fold_HG")        # H3: hidden state 1 regulates E, state 3 regulates G
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
Returns reporter counts [frame, trace, gene].
"""
function simulate_fast(model, cfg, rng)
    destinations, rates = system_transitions(model, model.base)
    start = cumsum(stationary(destinations, rates))         # used to draw the starting state of each trace
    cumulative = cumsum.(rates)                             # used to draw which event happens next
    counts = zeros(Int, cfg.nframes, cfg.ntraces, 2)
    for trace in 1:cfg.ntraces
        s = min(searchsortedfirst(start, rand(rng)), length(start))
        next_event = randexp(rng) / cumulative[s][end]
        for frame in 1:cfg.nframes
            while next_event <= (frame - 1) * cfg.interval
                choice = searchsortedfirst(cumulative[s], rand(rng) * cumulative[s][end])
                s = destinations[s][min(choice, end)]
                next_event += randexp(rng) / cumulative[s][end]
            end
            _, bits_E, bits_G = sstate(model, s)
            counts[frame, trace, 1] = count_ones(bits_E)
            counts[frame, trace, 2] = count_ones(bits_G)
        end
    end
    counts
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
    #   A gene regulating from promoter state 1: "state 1 multiplies the rate by fold" is written the
    #       equivalent way round: the rate is fold times higher, and states 2 and 3 divide it by fold.
    #   R is StochasticGene's Rsum: one connection from each reporter position (source states 4, 5, 6).
    renumber(state) = model.nh == 3 ? (2, 1, 3)[state] : state
    hidden_transitions = model.nh == 3 ? ((2, 1), (1, 2), (1, 3), (3, 1)) : GENE_TRANSITIONS[1:2]    # order of h12 h21 h23 h32
    base = copy(model.base)
    connections = NTuple{4,Int}[]       # (source unit, source state, target unit, target rate)
    gammas = Float64[]                  # one coupling strength per connection
    for edge in model.edges
        if edge.source != 3 && edge.state == 1
            base[edge.target, edge.rate] *= 1 + edge.gamma
            for state in (2, 3)
                push!(connections, (unit(edge.source), state, unit(edge.target), edge.rate))
                push!(gammas, 1 / (1 + edge.gamma) - 1)
            end
        else
            states = edge.source == 3 ? (renumber(edge.state),) : edge.state == 4 ? (4, 5, 6) : (edge.state,)
            for state in states
                push!(connections, (unit(edge.source), state, unit(edge.target), edge.rate))
                push!(gammas, edge.gamma)
            end
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
    seed = cfg.seed + 1000 * parse(Int, row["sim_id"]) + replicate
    counts = if cfg.simulator == "fast"
        simulate_fast(model, cfg, Xoshiro(seed))
    else
        Random.seed!(seed)
        simulate_stochasticgene(model, cfg)
    end
    (; curves=correlation_curves(counts, cfg), metrics=on_and_reporters(counts, cfg))
end

fmt(x) = isfinite(x) ? @sprintf("%.8g", x) : "NA"
mean_se(values) = (mean(values), std(values) / sqrt(length(values)))

"Simulate rows first_id to last_id of simulations.csv and write their curves and metrics."
function run_simulations(rows, cfg, first_id, last_id)
    models = [build_model(row, cfg.fitted, cfg.rate_matching) for row in rows]
    nclipped = count(model -> model.clipped, models)
    nclipped > 0 && println("$nclipped of $(length(models)) simulations have a rate multiplier clipped at $MIN_FACTOR")
    jobs = [(i, replicate) for i in eachindex(rows) for replicate in 1:cfg.replicates]
    results = Vector{Any}(undef, length(jobs))
    Threads.@threads :dynamic for j in eachindex(jobs)
        results[j] = simulate_replicate(rows[jobs[j][1]], models[jobs[j][1]], jobs[j][2], cfg)
    end

    folder = mkpath(joinpath(cfg.run_dir, "results"))
    stem = @sprintf("%06d-%06d", first_id, last_id)
    curves_path = joinpath(folder, "curves_$stem.csv")
    metrics_path = joinpath(folder, "metrics_$stem.csv")
    score_names = ["$(name)_$(observable)_$(pair)" for observable in OBSERVABLES for (pair, _, _) in PAIRS
                   for name in ("nrmse", "ccc")]
    base_names = ["base_$(gene)_$(RATE_NAMES[k])" for gene in ("E", "G") for k in 1:5]      # R12 R21 R23 R32 Initiate
    metric_names = ("on_E", "on_G", "reporters_E", "reporters_G")
    replicate_names = ["$(name)_rep$(r)" for name in metric_names for r in 1:cfg.replicates]
    # Write to temporary files first, so that a killed job never leaves a half-written result.
    open(curves_path * ".tmp", "w") do curves_io
        open(metrics_path * ".tmp", "w") do metrics_io
            println(curves_io, "sim_id,observable,pair,lag_frames,lag_min,mean,se," * join(["rep$r" for r in 1:cfg.replicates], ','))
            println(metrics_io, "sim_id,model_id,on_E,on_E_se,on_G,on_G_se,reporters_E,reporters_E_se," *
                                "reporters_G,reporters_G_se," * join(vcat(score_names, base_names, replicate_names, "clipped"), ','))
            for (i, row) in enumerate(rows)
                replicates = results[(i - 1) * cfg.replicates .+ (1:cfg.replicates)]
                bins = replicates[1].curves                                   # same lag bins in every replicate
                values = hcat([[bin.value for bin in r.curves] for r in replicates]...)
                means = [mean_se(values[b, :]) for b in eachindex(bins)]
                for (b, bin) in enumerate(bins)                                # mean, SE, then every replicate
                    println(curves_io, join((row["sim_id"], bin.observable, bin.pair, bin.lag, fmt(bin.lag * cfg.interval),
                                             fmt.(means[b])..., fmt.(values[b, :])...), ','))
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
    command = "$(s["julia"]) --threads=$(s["threads"])$project $(abspath(@__FILE__)) $(abspath(config_path))"
    per_job = s["simulations_per_job"]
    lines = ["$command $start $(min(start + per_job - 1, nrows))" for start in 1:per_job:nrows]
    path = joinpath(cfg.run_dir, "simulations.swarm")
    write(path, join(lines, "\n") * "\n")
    println("$(length(lines)) swarm lines of up to $per_job simulations x $(cfg.replicates) replicates. Submit with:")
    println("swarm -f $path -t $(s["threads"]) -g $(s["memory_gb"]) --time $(s["walltime"]) ",
            "--module $(s["modules"]) --logdir $(joinpath(cfg.run_dir, "swarm_logs"))")
end

# ---------------------------------------------------------------------------

length(ARGS) in (1, 2, 3) || error("usage: julia 3_run_simulations.jl config.toml [--swarm | FIRST LAST]")
const CONFIG = read_config(ARGS[1])
const ROWS = read_csv(joinpath(CONFIG.run_dir, "simulations.csv"))
if CONFIG.simulator == "stochasticgene" && !("--swarm" in ARGS)
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
if "--swarm" in ARGS
    write_swarm(ARGS[1], CONFIG, length(ROWS))
else
    first_row, last_row = length(ARGS) == 3 ? parse.(Int, ARGS[2:3]) : (1, length(ROWS))
    run_simulations(ROWS[first_row:last_row], CONFIG, first_row, last_row)
end
