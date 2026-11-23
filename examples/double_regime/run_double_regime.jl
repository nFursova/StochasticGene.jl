# Load a hand-written run-spec TOML and call `fit`.
#
# `fit(key=...)` only reloads the companion info_<key>.jld2 that a previous run wrote.
# A hand-written TOML is therefore parsed here and its [run] table is passed as
# keyword arguments.
#
# Usage:  julia --project run_double_regime.jl [path/to/spec.toml]

using StochasticGene
using TOML

spec_file = isempty(ARGS) ? joinpath(@__DIR__, "double_regime_two_state.toml") : ARGS[1]
raw = TOML.parsefile(spec_file)["run"]

intvecs(v) = [Vector{Int}(x) for x in v]

spec = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in raw)

# TOML arrays come back as Vector{Any}; convert to the types fit expects.
spec[:transitions] = Tuple(intvecs(raw["transitions"]))
spec[:fixedeffects] = Tuple(intvecs(get(raw, "fixedeffects", Any[])))
spec[:fittedparam] = Vector{Int}(get(raw, "fittedparam", Int[]))
spec[:priormean] = Float64.(raw["priormean"])
spec[:onstates] = intvecs(raw["onstates"])
spec[:datapath] = String.(raw["datapath"])
spec[:inference_method] = Symbol(get(raw, "inference_method", "mh"))
if haskey(raw, "dwell_specs")
    spec[:dwell_specs] = [(unit=Int(d["unit"]),
                           onstates=intvecs(d["onstates"]),
                           dttype=String.(d["dttype"])) for d in raw["dwell_specs"]]
end

fit(; spec...)
