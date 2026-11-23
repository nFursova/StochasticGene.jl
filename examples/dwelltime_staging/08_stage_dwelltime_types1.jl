# Stage 12-hour ON/OFF dwell-time fits for MYC (gene) and CCAT1 (enhancer), model types 1.
#
# Model types 1 (no splicing), for R/I in 301 302 303 401 403 404 501 504 505:
#   2-state linear (2RSInone); 3-state linear, cyclic, kinetic proofreading (3RSInone/cyc/kp);
#   double-regime coupled model with a hidden 2-state regime (2RSIgior).
#   = 9 x 5 = 45 models per data unit.
#
# Layout (same as 06_swarm_dwelltime_expanded_topology.jl):
#   data    <SG_DWELL_ROOT>/<source>/<sample>_<datacond>_<DT>.csv   e.g. dwelltimes/model/3Prime_gene_ON.csv
#   folder  results/dwelltime-<source>-<datacond>-types1-<SG_TAG>
#   key     <sample>-dt-<model>-<datacond>                           e.g. 3Prime-dt-3301cyc-gene
#   swarm   swarm-dwelltime-<source>-<datacond>-types1-<SG_TAG>/fit.swarm  (one per folder)
#
# Environment (defaults):
#   SG_SOURCES=model  SG_DATACONDS=gene,enhancer  SG_SAMPLES_LIST=3Prime,5Prime
#   SG_TOPOLOGIES=none,cyc,kp,gior  SG_MODELS=(all)  SG_TAG=<date>-12h
#   SG_MAXTIME=43200  SG_SAMPLESTEPS=10^12  SG_SAMPLE_STRIDE=1000  SG_NCHAINS=16
#   SG_ELONGATION_GENE=3.2  SG_ELONGATION_ENHANCER=22  SG_PRIORCV_R=0.25
#   SG_GIOR_GAMMA=0.1  SG_GIOR_GAMMA_CV=10   (see dwelltime_staging_common.jl for the rest)

include(joinpath(@__DIR__, "dwelltime_staging_common.jl"))

const SOURCES = csvlist("SG_SOURCES", "model")
const DATACONDS = csvlist("SG_DATACONDS", "gene,enhancer")
const SAMPLES = csvlist("SG_SAMPLES_LIST", "3Prime,5Prime")

const MODELS = model_list(1)
print_models(MODELS)

units = [(label="$source/$(sample)_$datacond",
          source=source,
          datacond=datacond,
          datapath=[joinpath(DWELL_ROOT, source, "$(sample)_$(datacond)_$(dt).csv") for dt in DTTYPE],
          resultfolder="dwelltime-$source-$datacond-types1-$TAG",
          keyfn=m -> "$(sample)-dt-$(m.name)-$(datacond)")
         for source in SOURCES for datacond in DATACONDS for sample in SAMPLES]

# Keys do not contain the source, so each source x datacond gets its own result folder AND
# its own swarm folder (as in the earlier 06 batches).
for source in SOURCES, datacond in DATACONDS
    group = [u for u in units if u.source == source && u.datacond == datacond]
    folder = "dwelltime-$source-$datacond-types1-$TAG"
    println("\n================  ", folder, "  ================")
    runs = try
        build_runs(group, MODELS)
    catch err
        err isa ErrorException && occursin("no dwell-time files", err.msg) || rethrow()
        println("  no data for this group, skipped"); continue
    end
    stage!(runs; swarmdir=joinpath(ROOT, "swarm-" * folder), default_resultfolder=folder)
end
