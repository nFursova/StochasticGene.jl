# Stage A for MYC (gene): same as 12_stage_trace_separate_noise.jl, without splicing.
#
# Separate 3Prime and 5Prime fits (nothing shared), noise means per trace, noise SDs per sample,
# all Rshift rates tied, Eject separate, same priors as the enhancer script.
# Models: 3301 3302 3303 3401 3402 3403 3404  x  none cyc kp  = 21 per sample (42 runs).
#
# Defaults set here (any of them can still be overridden by exporting it):
#   SG_CONDITION=gene
#   SG_BASES=3301,3302,3303,3401,3402,3403,3404
#   SG_RESULTFOLDER=3and5Prime-stage05-gene-$SG_STAMP
#   SG_SWARMDIR=swarm-3and5Prime-stage05-gene-$SG_STAMP
# SG_STAMP defaults to today's date (YYYY-MM-DD).
#
# Run:
#   export SG_STAMP=$(date +%Y-%m-%d)
#   julia --project="$SG_PROJECT" scripts/12_stage_trace_separate_noise_gene.jl
# then, after the 12 h runs:
#   SG_FROM=3and5Prime-stage05-gene-$SG_STAMP julia --project="$SG_PROJECT" scripts/13_continue_trace_fixnoise.jl

using Dates

const STAMP = get(ENV, "SG_STAMP", Dates.format(Dates.today(), "yyyy-mm-dd"))
get!(ENV, "SG_CONDITION", "gene")
get!(ENV, "SG_BASES", "3301,3302,3303,3401,3402,3403,3404")
get!(ENV, "SG_RESULTFOLDER", "3and5Prime-stage05-gene-$STAMP")
get!(ENV, "SG_SWARMDIR", "swarm-3and5Prime-stage05-gene-$STAMP")

ENV["SG_CONDITION"] == "gene" || error("this wrapper stages gene fits; use 12_stage_trace_separate_noise.jl for enhancer")
any(b -> length(b) == 4 && b[3] != '0', split(ENV["SG_BASES"], ',')) &&
    error("SG_BASES contains a splicing model; the gene staging is without splicing")

include(joinpath(@__DIR__, "12_stage_trace_separate_noise.jl"))
