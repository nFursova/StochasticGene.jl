# Runtime patch: coupled dwell-time fits (the "gior" models) on an unpatched StochasticGene.
#
# Tested on StochasticGene 2.0.7. Load it AFTER `using StochasticGene`, on every worker:
#
#     @everywhere using StochasticGene
#     @everywhere include("/path/to/stochasticgene_coupled_dwell_patch.jl")
#
# The staging scripts add these lines to the gior fit scripts automatically. Nothing is
# installed or changed on disk; the methods below replace StochasticGene's in the running
# Julia session only. If StochasticGene already contains the fix (it defines
# `_coupled_dwell_layout`), the patch does nothing.
#
# What it fixes (both only affect coupled models fitted to dwell-time data):
#  1. num_all_parameters counted each unit's ON/OFF state lists as noise parameters, so the
#     coupling strength was looked up at the wrong position of the rate vector (BoundsError).
#  2. fit() reads ON/OFF histograms as one flat list, while the coupled predictor indexes them
#     per unit and also tried to compute dwell times for the hidden unit. The flat data are now
#     mapped onto the observed unit(s), hidden units are skipped, and predictions come back in
#     the order of the data files.

if isdefined(StochasticGene, :_coupled_dwell_layout)
    @info "StochasticGene already includes the coupled dwell-time fix; patch not applied"
else
    @eval StochasticGene begin

        function num_all_parameters(transitions, R::Int, S, insertstep, reporter, coupling=tuple(), grid=nothing)
            if typeof(reporter) <: HMMReporter
                n = reporter.n
            elseif typeof(reporter) <: Vector
                # numbers = noise priors; a vector of state sets (dwell sojourn lists) has none
                n = all(x -> x isa Real, reporter) ? length(reporter) : 0
            elseif typeof(reporter) <: Int
                n = reporter
            else
                n = 0
            end
            c = isempty(coupling) ? 0 : length(coupling[2])
            g = isnothing(grid) ? 0 : 1
            num_rates(transitions, R, S, insertstep) + n + c + g
        end

        function _coupled_dwell_layout(data::DwellTimeData, reporter)
            sojourn = reporter[1]
            nunits = length(sojourn)
            observed = isempty(data.units) ? collect(1:nunits) : collect(Int, data.units)
            nested = !isempty(data.bins) && first(data.bins) isa AbstractVector && !isempty(first(data.bins)) &&
                     first(first(data.bins)) isa AbstractVector
            nested && return data.bins, data.DTtypes, collect(1:nunits), true
            all(u -> 1 <= u <= nunits, observed) ||
                throw(ArgumentError("dwell-time data units $(observed) are outside 1:$nunits"))
            counts = [length(sojourn[u]) for u in observed]
            length(data.bins) == sum(counts) ||
                throw(ArgumentError("dwell-time data have $(length(data.bins)) histograms but observed units $(observed) need $(sum(counts))"))
            bins = [Vector{Vector{Float64}}() for _ in 1:nunits]
            dttype = [String[] for _ in 1:nunits]
            k = 0
            for (u, n) in zip(observed, counts)
                for i in 1:n
                    k += 1
                    push!(bins[u], Float64.(data.bins[k]))
                    push!(dttype[u], String(length(data.DTtypes) == sum(counts) ? data.DTtypes[k] : data.DTtypes[i]))
                end
            end
            return bins, dttype, observed, false
        end

        function predictedarray(r, data::DwellTimeData, model::AbstractGRSMmodel{T}; steady_state_solver::Symbol=:default) where {T<:NamedTuple{(:coupling,)}}
            r, coupling_strength = prepare_rates_coupled(r, model.nrates, model.trait.coupling.couplingindices)
            components = get_components(model, data)
            reporter = get_reporter(model, data)
            components isa TDCoupledFullComponents ||
                return predictedarray(r, coupling_strength, components, data.bins, reporter, data.DTtypes; steady_state_solver=steady_state_solver, method=model.method)
            bins, dttype, observed, nested = _coupled_dwell_layout(data, reporter)
            hists = predictedarray(r, coupling_strength, components, bins, reporter, dttype; steady_state_solver=steady_state_solver, method=model.method)
            nested && return hists
            return reduce(vcat, [hists[u] for u in observed])
        end

        function predictedarray(r, coupling_strength, components::TDCoupledFullComponents, bins, reporter, dttype; steady_state_solver::Symbol=:default, method=Tsit5())
            sojourn_full, _ = reporter
            coupling_rates = [coupling_strength[k] * r[components.targets[k][1]][components.targets[k][2]]
                              for k in eachindex(components.targets)]
            T = make_mat_TC(components, r, coupling_rates)
            pss = steady_state_vector(T; solver=steady_state_solver)
            hists = Vector{Vector}[]
            for α in eachindex(sojourn_full)
                h = Vector[]
                if α > length(bins) || isempty(bins[α])   # hidden unit: nothing to predict
                    push!(hists, h)
                    continue
                end
                for i in eachindex(sojourn_full[α])
                    soj = sojourn_full[α][i]
                    TCD = make_mat_TCD(components, α, i, r, coupling_rates)
                    init = init_S(soj, T, pss, dttype[α][i])
                    init = init[soj]
                    s = sum(init)
                    (iszero(s) || !isfinite(s)) && continue
                    init = init / s
                    push!(h, dwelltimePDF(bins[α][i], TCD[soj, soj], 1:length(soj), init, method))
                end
                push!(hists, h)
            end
            return hists
        end

    end
    @info "Applied coupled dwell-time patch to StochasticGene $(pkgversion(StochasticGene))"
end
