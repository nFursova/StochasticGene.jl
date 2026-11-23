# Stage a key-based run from a hand-written run-spec TOML.
#
# Converts the TOML [run] table to Julia types, then calls `stage_write_run_specs`, which writes
#   <root>/results/<resultfolder>/info_<key>.jld2   (read by fit(; key=...))
#   <root>/results/<resultfolder>/info_<key>.toml   (marker file)
#   <filedir>/fitscript_<key>.jl                    (calls fit(; key=..., resultfolder, root))
#   <filedir>/fit.swarm                             (julia -p <nchains> ... fitscript_<key>.jl)
#
# Usage:
#   julia --project stage_double_regime.jl [spec.toml] [filedir]
#   cd <filedir> && bash fit.swarm        # or: swarm -f fit.swarm  (Biowulf)

using StochasticGene
using TOML

spec_file = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "double_regime_coupled.toml")
filedir = length(ARGS) >= 2 ? ARGS[2] : "swarm-doubleregime"

raw = TOML.parsefile(spec_file)["run"]

ints(v) = Vector{Int}(v)
intvecs(v) = [ints(x) for x in v]
tup(v) = Tuple(Int.(v))

spec = Dict{Symbol,Any}(Symbol(k) => v for (k, v) in raw)

# ---- model structure: per-model tuples for coupled runs ----
coupled = raw["G"] isa AbstractVector
if coupled
    spec[:G] = tup(raw["G"])
    spec[:R] = tup(raw["R"])
    spec[:S] = tup(raw["S"])
    spec[:insertstep] = tup(raw["insertstep"])
    spec[:transitions] = Tuple(Tuple(intvecs(t)) for t in raw["transitions"])
    spec[:onstates] = [intvecs(o) for o in raw["onstates"]]
    haskey(raw, "elongationtime") && (spec[:elongationtime] = Tuple(Float64.(raw["elongationtime"])))
else
    spec[:transitions] = Tuple(intvecs(raw["transitions"]))
    spec[:onstates] = intvecs(raw["onstates"])
end

# ---- coupling = (unit_model, [(beta, s, alpha, t), ...], sign_modes) ----
if haskey(raw, "coupling_unit_model")
    unit_model = tup(raw["coupling_unit_model"])
    connections = [Tuple(Int.(c)) for c in raw["coupling_connections"]]
    modes = Symbol.(get(raw, "coupling_modes", fill("free", length(connections))))
    spec[:coupling] = (unit_model, connections, modes)
    foreach(k -> delete!(spec, k), (:coupling_unit_model, :coupling_connections, :coupling_modes))
end

# ---- parameters ----
spec[:fittedparam] = ints(raw["fittedparam"])
spec[:fixedeffects] = Tuple(intvecs(get(raw, "fixedeffects", Any[])))
spec[:priormean] = Float64.(raw["priormean"])
spec[:priorcv] = raw["priorcv"] isa AbstractVector ? Float64.(raw["priorcv"]) : Float64(raw["priorcv"])

# ---- data ----
spec[:datapath] = String.(raw["datapath"])
if haskey(raw, "dwell_specs")
    spec[:dwell_specs] = [(unit=Int(d["unit"]),
                           onstates=intvecs(d["onstates"]),
                           dttype=String.(d["dttype"])) for d in raw["dwell_specs"]]
end

# ---- inference ----
spec[:inference_method] = Symbol(get(raw, "inference_method", "mh"))

# Sanity check: one prior per rate (+ coupling constants)
nr = StochasticGene.num_rates(spec[:transitions], spec[:R], spec[:S], spec[:insertstep])
ncoup = haskey(spec, :coupling) ? length(spec[:coupling][2]) : 0
ntotal = sum(nr) + ncoup
length(spec[:priormean]) == ntotal ||
    error("priormean has $(length(spec[:priormean])) entries, model needs $ntotal (rates per model $nr + $ncoup coupling)")
spec[:priorcv] isa AbstractVector && length(spec[:priorcv]) != ntotal &&
    error("priorcv has $(length(spec[:priorcv])) entries, model needs $ntotal")
println("rates per model = $nr, coupling constants = $ncoup, total = $ntotal")

key = String(spec[:key])
out = stage_write_run_specs(Dict(key => spec);
    resultfolder=String(spec[:resultfolder]),
    root=String(spec[:root]),
    filedir=filedir,
    swarmfile="fit",
    nchains=Int(spec[:nchains]),
    nthreads=1,
    script_override_keys=(),   # everything comes from info_<key>.jld2
)

println("staged key:  ", key)
println("run spec:    ", out.specs[1].jld2)
println("fit script:  ", out.scripts[1])
println("swarm file:  ", out.commandfiles[1])
