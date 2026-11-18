#=
Compare the ON/OFF dwell-time histograms with the fitted models (step 1 of 2; dwelltime_report.py draws the report).

    julia --project=$SG_PROJECT dwelltime_report.jl RESULTFOLDER [RESULTFOLDER ...]

RESULTFOLDER: a results folder of dwell-time fits, under $SG_ROOT/results (default root
/data/Larsonlab/Nadya/StochasticGene) or an absolute path. Every fit with an info_<key>.jld2 and a rates_<key>.txt
is rebuilt from its saved run spec with its maximum-likelihood rates, and evaluated with StochasticGene's own
dwell-time likelihood. Written into <RESULTFOLDER>/dwelltime_report/:

    fits.csv           one row per fit: key, data (the histogram files), sample, model, parameters, log-likelihood,
                       AIC, number of events, Kolmogorov-Smirnov distance of each dwell type, largest Rhat, smallest ESS
    distributions.csv  one row per fit, dwell type and bin: time, data probability, model probability

Log-likelihood = sum over the bins of both dwell types of (events in the bin) x log(model probability of the bin);
the model probabilities of each type sum to 1 over the bins of the histogram, as in the fit.
AIC = 2 x (fitted parameters) - 2 x log-likelihood. AIC compares only fits of the same histograms (same source,
sample and condition); dwelltime_report.py compares within those groups.
=#
using StochasticGene, DelimitedFiles

const ROOT = abspath(get(ENV, "SG_ROOT", "/data/Larsonlab/Nadya/StochasticGene"))

"""The dwell types of a fit: the spec's dttype if it is there, else the end of each data file name
(<sample>_<condition>_ON.csv -> ON), as the staging scripts name them. fit() drops dttype from the run spec it
saves, and dwell_specs stays empty, so a finished fit does not record its dwell types."""
function dwell_types(spec)
    haskey(spec, :dttype) && return String.(collect(spec[:dttype]))
    [String(match(r"_([A-Za-z]+)\.[A-Za-z]+$", basename(string(path)))[1]) for path in spec[:datapath]]
end

"""Data and model of one fit, rebuilt from its run spec with the maximum-likelihood rates (rates file row "ml"),
as fit(; key) builds them."""
function data_and_model(resultdir, key, spec)
    s = merge(StochasticGene.fit_default_spec(), spec)
    rinit = StochasticGene.readrates(joinpath(resultdir, "rates_$key.txt"), StochasticGene.get_row("ml"))
    data, model, _ = StochasticGene.make_structures(rinit, s[:datatype], s[:datapath], s[:gene], s[:cell], s[:datacond],
        resultdir, key, s[:fittedparam], s[:fixedeffects], s[:transitions], s[:G], s[:R], s[:S], s[:insertstep],
        s[:coupling], s[:grid], ROOT, s[:maxtime], s[:elongationtime], s[:priormean], s[:priorcv], s[:nalleles],
        s[:onstates], s[:decayrate], s[:splicetype], s[:probfn], s[:noisepriors], s[:hierarchical], "ml", s[:propcv],
        1, 0, s[:temp], s[:temprna], s[:method], s[:zeromedian], s[:datacol], s[:ejectnumber], s[:yieldfactor],
        s[:trace_specs], s[:dwell_specs]; legacy_dttype=dwell_types(spec))
    data, model
end

"Largest Rhat and smallest ESS over the fitted parameters, from measures_<key>.txt (rows 5 and 6), or NaN."
function convergence(resultdir, key)
    file = joinpath(resultdir, "measures_$key.txt")
    isfile(file) || return NaN, NaN
    rows = readlines(file)
    row(i) = length(rows) >= i ? filter(isfinite, something.(tryparse.(Float64, split(rows[i], ',')), NaN)) : Float64[]
    rhat, ess = row(5), row(6)
    (isempty(rhat) ? NaN : maximum(rhat)), (isempty(ess) ? NaN : minimum(ess))
end

csv(x::AbstractString) = occursin(r"[,\"]", x) ? "\"" * replace(x, "\"" => "\"\"") * "\"" : x
csv(x::Real) = isfinite(x) ? string(round(x; sigdigits=10)) : "NA"

function report(folder)
    resultdir = isabspath(folder) ? folder : joinpath(ROOT, "results", folder)
    isdir(resultdir) || error("results folder not found: $resultdir")
    keys = sort([m[1] for m in (match(r"^info_(.+)\.jld2$", f) for f in readdir(resultdir))
                 if m !== nothing && isfile(joinpath(resultdir, "rates_$(m[1]).txt"))])
    isempty(keys) && error("no finished fits (info_<key>.jld2 with rates_<key>.txt) in $resultdir")
    out = mkpath(joinpath(resultdir, "dwelltime_report"))
    open(joinpath(out, "fits.csv"), "w") do fits_io
        open(joinpath(out, "distributions.csv"), "w") do dist_io
            println(fits_io, "folder,key,data,sample,model,G,R,S,insertstep,parameters,loglik,AIC,events_ON,events_OFF,KS_ON,KS_OFF,max_rhat,min_ess")
            println(dist_io, "key,dttype,time,data,model")
            for key in keys
                spec = StochasticGene.read_run_spec(joinpath(resultdir, "info_$key.jld2"))
                get(spec, :datatype, "") == "dwelltime" || (println("skip $key: not a dwell-time fit"); continue)
                data, model = data_and_model(resultdir, String(key), spec)
                param = StochasticGene.get_param(model)
                ll = StochasticGene.loglikelihood(param, data, model)[1]
                predicted = StochasticGene.predictedfn(param, data, model)
                events, ks = Dict{String,Float64}(), Dict{String,Float64}()
                start = 1
                for (i, dttype) in enumerate(data.DTtypes)
                    counts = Float64.(data.DwellTimes[i])
                    p = predicted[start:start + length(counts) - 1]
                    start += length(counts)
                    q = counts ./ sum(counts)
                    events[dttype] = sum(counts)
                    ks[dttype] = maximum(abs.(cumsum(q) .- cumsum(p)))
                    for b in eachindex(counts)
                        println(dist_io, join((csv(key), dttype, csv(data.bins[i][b]), csv(q[b]), csv(p[b])), ','))
                    end
                end
                parts = split(key, '-')
                rhat, ess = convergence(resultdir, key)
                println(fits_io, join((csv(folder), csv(key), csv(join(string.(spec[:datapath]), ';')), csv(parts[1]),
                                      csv(join(parts[3:end-1], '-')), spec[:G], spec[:R], spec[:S], spec[:insertstep],
                                      length(param), csv(ll), csv(2 * length(param) - 2 * ll),
                                      csv(get(events, "ON", NaN)), csv(get(events, "OFF", NaN)),
                                      csv(get(ks, "ON", NaN)), csv(get(ks, "OFF", NaN)), csv(rhat), csv(ess)), ','))
                println("$key: log-likelihood $(round(ll; digits=1)), AIC $(round(2 * length(param) - 2 * ll; digits=1))")
            end
        end
    end
    println("wrote $(joinpath(out, "fits.csv")) and distributions.csv")
end

isempty(ARGS) && error("usage: julia --project=\$SG_PROJECT dwelltime_report.jl RESULTFOLDER [RESULTFOLDER ...]")
foreach(report, ARGS)
