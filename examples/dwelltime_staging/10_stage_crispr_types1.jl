# Stage 12-hour ON/OFF dwell-time fits across CRISPR conditions, MYC + CCAT1, model types 1.
#
# Same models as 08_stage_dwelltime_types1.jl (45 per data unit).
#
# Layout (same as 07_swarm_dwelltime_conditions.jl):
#   data    <SG_DWELL_ROOT>/<source>/<condition>/<sample>_<datacond>_<DT>.csv
#   folder  results/CRISPR-dwelltime-types1-<SG_TAG>          (one folder; keys separate units)
#   key     <source>-<condition>-<sample>-<datacond>-<model>   e.g. model-UNT-5Prime-gene-3301kp
# Units without files (e.g. KRAB enhancer) are skipped.
#
# Environment (defaults):
#   SG_SOURCES=model  SG_CRISPR_CONDITIONS=UNT,VPR,shKRAB,KRAB  SG_SAMPLES_LIST=5Prime
#   SG_DATACONDS=gene,enhancer  plus everything listed in 08_stage_dwelltime_types1.jl

include(joinpath(@__DIR__, "dwelltime_staging_common.jl"))

const SOURCES = csvlist("SG_SOURCES", "model")
const CONDITIONS = csvlist("SG_CRISPR_CONDITIONS", "UNT,VPR,shKRAB,KRAB")
const SAMPLES = csvlist("SG_SAMPLES_LIST", "5Prime")
const DATACONDS = csvlist("SG_DATACONDS", "gene,enhancer")
const RESULTFOLDER = envstr("SG_RESULTFOLDER", "CRISPR-dwelltime-types1-$TAG")
const SWARMDIR = joinpath(ROOT, envstr("SG_SWARMDIR", "swarm-CRISPR-dwelltime-types1-$TAG"))

const MODELS = model_list(1)
print_models(MODELS)

units = [(label="$source/$cond/$(sample)_$datacond",
          datacond=datacond,
          datapath=[joinpath(DWELL_ROOT, source, cond, "$(sample)_$(datacond)_$(dt).csv") for dt in DTTYPE],
          resultfolder=RESULTFOLDER,
          keyfn=m -> "$source-$cond-$sample-$datacond-$(m.name)")
         for source in SOURCES for cond in CONDITIONS for sample in SAMPLES for datacond in DATACONDS]

stage!(build_runs(units, MODELS); swarmdir=SWARMDIR, default_resultfolder=RESULTFOLDER)
