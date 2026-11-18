#!/usr/bin/env bash
# Stage the 25 sweep candidates as coupled fits: 3Prime, target mode, 24 hours.
#   9 direct models (models_direct.txt) with your scripts/stage_coupled_reciprocal.jl
#   16 hidden-unit models (models_hidden.txt) with stage_coupled_hidden.jl (this folder)
# Everything as in your earlier coupled runs except the architectures and propcv = 0.1.
# Run from anywhere:  SG_PROJECT=/data/fursovan2/sg111 bash stage_candidates_24h.sh
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export SG_ROOT=/data/Larsonlab/Nadya/StochasticGene
: "${SG_PROJECT:?set SG_PROJECT to your StochasticGene project, e.g. /data/fursovan2/sg111}"
export SG_PROJECT

export SG_SRC_GENE_FOLDER=3and5Prime-stage05-gene-2026-08-25-fixnoise-R-median-continue-24h
export SG_SRC_ENH_FOLDER=3and5Prime-stage05-enhancer-2026-08-25-fixnoise-R-median-continue-24h
export SG_SRC_ROW=median
export SG_MAXTIME=86400          # 24 hours (the continuation adds 72)
export SG_MEMORY_GB=60
export SG_NCHAINS=16
export SG_SAMPLES=10000000
export SG_STRIDE=100
export SG_PROPCV=0.1
export SG_WARMUP=0
export SG_PRIORCV_TARGET=1.0
export SG_PRIORCV_G=1.0
export SG_PRIORCV_INITIATE=0.25
export SG_PRIORCV_COUPLING=10.0
export SG_GAMMA0=0.0             # every gamma free, starting at 0 (no coupling)
export SG_ELONG_ENH=20.0
export SG_ELONG_GENE=5.0
export SG_INTERVAL=1.6666666666666667
# hidden-unit models only
export SG_HIDDEN_RATE0=0.1       # hidden switching rates start at 0.1 per minute ...
export SG_FIT_HIDDEN=1           # ... and are fitted
export SG_PRIORCV_HIDDEN=10.0

export SG_SAMPLE=3Prime
export SG_COUPLE_MODE=target
export SG_DATA=traces/3Prime
export SG_SRC_FOLDER="$SG_SRC_GENE_FOLDER"
export SG_SRC_RATES_GENE="$SG_ROOT/results/$SG_SRC_GENE_FOLDER/rates_3Prime-5Prime-h2-3301none-gene-3Prime.txt"
export SG_SRC_RATES_ENH="$SG_ROOT/results/$SG_SRC_ENH_FOLDER/rates_3Prime-5Prime-h2-3301none-enhancer-3Prime.txt"
[ -d "$SG_ROOT/data/$SG_DATA" ] || { echo "FAIL: missing $SG_ROOT/data/$SG_DATA" >&2; exit 1; }

cd "$SG_ROOT"

export SG_MODELFILE="$HERE/models_direct.txt"
export SG_INTO=3Prime-coupled-target-candidates-direct-24h
export SG_SWARMROOT="swarm_$SG_INTO"
julia --project="$SG_PROJECT" scripts/stage_coupled_reciprocal.jl

export SG_MODELFILE="$HERE/models_hidden.txt"
export SG_INTO=3Prime-coupled-target-candidates-hidden-24h
export SG_SWARMROOT="swarm_$SG_INTO"
julia --project="$SG_PROJECT" "$HERE/stage_coupled_hidden.jl"

echo
echo "Submit both sets:"
echo "  bash $SG_ROOT/swarm_3Prime-coupled-target-candidates-direct-24h/submit_all.sh"
echo "  bash $SG_ROOT/swarm_3Prime-coupled-target-candidates-hidden-24h/submit_all.sh"
