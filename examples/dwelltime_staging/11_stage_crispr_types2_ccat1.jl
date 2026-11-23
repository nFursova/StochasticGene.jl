# Stage 12-hour ON/OFF dwell-time fits across CRISPR conditions, CCAT1 only, model types 2
# (splicing; same 30 models as 09_stage_dwelltime_types2_ccat1.jl).
#
# Layout:
#   data    <SG_DWELL_ROOT>/<source>/<condition>/<sample>_enhancer_<DT>.csv
#   folder  results/CRISPR-dwelltime-splicing-<SG_TAG>
#   key     <source>-<condition>-<sample>-enhancer-<model>
# Units without files (e.g. KRAB, which has no enhancer data) are skipped.

include(joinpath(@__DIR__, "dwelltime_staging_common.jl"))

const SOURCES = csvlist("SG_SOURCES", "model")
const CONDITIONS = csvlist("SG_CRISPR_CONDITIONS", "UNT,VPR,shKRAB,KRAB")
const SAMPLES = csvlist("SG_SAMPLES_LIST", "5Prime")
const DATACONDS = csvlist("SG_DATACONDS", "enhancer")
const RESULTFOLDER = envstr("SG_RESULTFOLDER", "CRISPR-dwelltime-splicing-$TAG")
const SWARMDIR = joinpath(ROOT, envstr("SG_SWARMDIR", "swarm-CRISPR-dwelltime-splicing-$TAG"))

const MODELS = model_list(2)
print_models(MODELS)

units = [(label="$source/$cond/$(sample)_$datacond",
          datacond=datacond,
          datapath=[joinpath(DWELL_ROOT, source, cond, "$(sample)_$(datacond)_$(dt).csv") for dt in DTTYPE],
          resultfolder=RESULTFOLDER,
          keyfn=m -> "$source-$cond-$sample-$datacond-$(m.name)")
         for source in SOURCES for cond in CONDITIONS for sample in SAMPLES for datacond in DATACONDS]

stage!(build_runs(units, MODELS); swarmdir=SWARMDIR, default_resultfolder=RESULTFOLDER)
