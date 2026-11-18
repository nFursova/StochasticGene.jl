#!/usr/bin/env bash
# After the 24-hour runs end: restage every model from where it stopped and run 72 hours more,
# with your scripts/30_continue_coupled.jl, as for your earlier "-continue-72h" runs.
#   SG_PROJECT=/data/fursovan2/sg111 bash continue_candidates_72h.sh
set -euo pipefail
export SG_ROOT=/data/Larsonlab/Nadya/StochasticGene
: "${SG_PROJECT:?set SG_PROJECT}"
export SG_PROJECT

FOLDERS=(
  3Prime-coupled-target-candidates-direct-24h
  3Prime-coupled-target-candidates-hidden-24h
)
export SG_FROM="$(IFS=,; echo "${FOLDERS[*]}")"
export SG_INTO_SUFFIX='-continue-72h'
export SG_MAXTIME=259200         # 72 hours
export SG_RATETYPE=median
export SG_MEMORY_GB=90           # 30_continue_coupled.jl sizes the swarm memory itself; this replaces it below
unset SG_SWARMROOT SG_ONLY SG_PROPCV SG_WARMUPSTEPS SG_NCHAINS

cd "$SG_ROOT"
julia --project="$SG_PROJECT" scripts/30_continue_coupled.jl

for from in "${FOLDERS[@]}"; do
  submit="$SG_ROOT/swarm_${from}-continue-72h/submit_all.sh"
  sed -i -E "s/-g [0-9]+/-g $SG_MEMORY_GB/" "$submit"
  bash "$submit"
done
