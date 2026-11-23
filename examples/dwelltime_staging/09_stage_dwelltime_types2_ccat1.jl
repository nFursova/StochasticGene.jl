# Stage 12-hour ON/OFF dwell-time fits for CCAT1 (enhancer), model types 2 (splicing).
#
# Model types 2, one splice step or splicing at every step:
#   R/I = 301, 401, 501 (insertstep 1): S = 1 (splice at R1)  and  S = R (splice at R1..R)
#   R/I = 303, 404, 505 (insertstep R): S = R (one splice at the last step; the only option
#                                       StochasticGene allows when insertstep = R)
#   -> 9 layouts x (2none, 3none, 3cyc, 3kp, 2gior) = 45 models per data unit, e.g.
#      2311none, 2331none, 2333none, 3441kp, 2551gior
#
# Layout:
#   data    <SG_DWELL_ROOT>/<source>/<sample>_enhancer_<DT>.csv
#   folder  results/dwelltime-<source>-enhancer-splicing-<SG_TAG>
#   key     <sample>-dt-<model>-enhancer
#   swarm   swarm-dwelltime-<source>-enhancer-splicing-<SG_TAG>/fit.swarm  (one per folder)
#
# Environment: as 08_stage_dwelltime_types1.jl, plus SG_SPLICE_PRIOR=0.05 SG_PRIORCV_SPLICE=0.2.

include(joinpath(@__DIR__, "dwelltime_staging_common.jl"))

const SOURCES = csvlist("SG_SOURCES", "model,segm")
const DATACONDS = csvlist("SG_DATACONDS", "enhancer")
const SAMPLES = csvlist("SG_SAMPLES_LIST", "3Prime,5Prime")

const MODELS = model_list(2)
print_models(MODELS)

units = [(label="$source/$(sample)_$datacond",
          source=source,
          datacond=datacond,
          datapath=[joinpath(DWELL_ROOT, source, "$(sample)_$(datacond)_$(dt).csv") for dt in DTTYPE],
          resultfolder="dwelltime-$source-$datacond-splicing-$TAG",
          keyfn=m -> "$(sample)-dt-$(m.name)-$(datacond)")
         for source in SOURCES for datacond in DATACONDS for sample in SAMPLES]

# Keys do not contain the source, so each source x datacond gets its own result folder AND
# its own swarm folder (as in the earlier 06 batches).
for source in SOURCES, datacond in DATACONDS
    group = [u for u in units if u.source == source && u.datacond == datacond]
    folder = "dwelltime-$source-$datacond-splicing-$TAG"
    println("\n================  ", folder, "  ================")
    runs = try
        build_runs(group, MODELS)
    catch err
        err isa ErrorException && occursin("no dwell-time files", err.msg) || rethrow()
        println("  no data for this group, skipped"); continue
    end
    stage!(runs; swarmdir=joinpath(ROOT, "swarm-" * folder), default_resultfolder=folder)
end
