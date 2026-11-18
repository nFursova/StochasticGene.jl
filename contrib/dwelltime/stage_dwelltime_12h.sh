#!/usr/bin/env bash
# Restage every ON/OFF dwell-time fit from scratch, 12 hours each, into new folders:
#   expanded topologies (scripts/06_swarm_dwelltime_expanded_topology.jl):
#       {model, segm} histograms x {gene, enhancer}, samples 3Prime and 5Prime, the 16 default topologies
#   enhancer splicing models (scripts/06_swarm_dwelltime_enhancer_splicing.jl):
#       {model, segm} histograms x {3Prime, 5Prime}
# Every batch sets all the variables it uses, so nothing carries over from the batch before it.
#
#   SG_PROJECT=/data/fursovan2/sg111 bash stage_dwelltime_12h.sh                  stage only
#   SG_PROJECT=/data/fursovan2/sg111 SG_SUBMIT=true bash stage_dwelltime_12h.sh   stage and submit
set -euo pipefail
: "${SG_ROOT:=/data/Larsonlab/Nadya/StochasticGene}"
: "${SG_PROJECT:?set SG_PROJECT, e.g. /data/fursovan2/sg111}"
: "${SG_TAG:=$(date +%F)-12h}"          # part of every new folder name
: "${SG_SUBMIT:=false}"
export SG_ROOT SG_PROJECT

# The same for every batch
export SG_MAXTIME=43200                 # 12 hours: the time limit ends each run ...
export SG_SAMPLESTEPS=1000000000        # ... because this step limit is never reached (samples are stored as they
                                        # come, so a large limit costs no memory)
export SG_SAMPLE_STRIDE=1000            # twice the 6-hour runs' 500: about as many stored samples in twice the time
export SG_NCHAINS=16
export SG_PROPCV=0.02
export SG_WRITESAMPLES=true
export SG_DTTYPE=ON,OFF
WALLTIME=14:00:00                       # SG_MAXTIME plus start-up and writing the results
MEMORY=32

cd "$SG_ROOT"
julia --project="$SG_PROJECT" -e 'using Pkg; Pkg.precompile()'   # once, before many workers start

# Write the swarm's worker settings and submit it (if SG_SUBMIT=true).
submit() {
    local swarmdir="$SG_ROOT/$1"
    sed -i 's|^julia |JULIA_WORKER_TIMEOUT=600 JULIA_PKG_PRECOMPILE_AUTO=0 julia |' "$swarmdir/fit.swarm"
    mkdir -p "$swarmdir/logs"
    if [[ "$SG_SUBMIT" == true ]]; then
        (cd "$swarmdir" && swarm -f fit.swarm -g "$MEMORY" -t $((SG_NCHAINS + 2)) --time "$WALLTIME" \
            --module julialang --logdir logs)
    fi
}

# Expanded topologies: one results folder per histogram source and condition.
stage_expanded() {
    local source=$1 condition=$2 elongation=$3
    unset SG_MODELS SG_SAMPLE                               # the 16 default topologies, both samples
    export SG_SAMPLES_LIST=3Prime,5Prime
    export SG_DWELL_DIR="dwelltimes/$source"
    export SG_CONDITION="$condition"
    export SG_ELONGATION="$elongation"
    export SG_PRIORCV_R=0.25
    export SG_RESULTFOLDER="dwelltime-$source-$condition-$SG_TAG"
    export SG_SWARMDIR="swarm-$SG_RESULTFOLDER"
    echo; echo "== results/$SG_RESULTFOLDER"
    julia --project="$SG_PROJECT" scripts/06_swarm_dwelltime_expanded_topology.jl
    submit "$SG_SWARMDIR"
}

# Enhancer splicing models: one results folder per histogram source (both samples in it), one swarm per sample.
stage_splicing() {
    local source=$1 sample=$2 models=$3
    unset SG_SAMPLES_LIST SG_CONDITION
    export SG_SAMPLE="$sample"
    export SG_MODELS="$models"
    export SG_DWELL_DIR="dwelltimes/$source"
    export SG_ELONGATION=20.0
    export SG_SPLICE_PRIOR=0.05
    export SG_PRIORCV_R=0.25
    export SG_PRIORCV_SPLICE=0.2
    export SG_RESULTFOLDER="dwelltime-$source-enhancer-splicing-$SG_TAG"
    export SG_SWARMDIR="swarm-dwelltime-$source-enhancer-splicing-$sample-$SG_TAG"
    echo; echo "== results/$SG_RESULTFOLDER ($sample)"
    julia --project="$SG_PROJECT" scripts/06_swarm_dwelltime_enhancer_splicing.jl
    submit "$SG_SWARMDIR"
}

for source in model segm; do
    stage_expanded "$source" gene 3.2
    stage_expanded "$source" enhancer 22
    stage_splicing "$source" 3Prime 3301,3302,3303,3331,3322
    stage_splicing "$source" 5Prime 3301,3331,3322,3311
done

echo
[[ "$SG_SUBMIT" == true ]] && echo "Submitted." || echo "Staged, not submitted: rerun with SG_SUBMIT=true to submit."
