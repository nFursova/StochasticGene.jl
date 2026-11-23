# 12-hour dwell-time staging (08-11)

Four staging scripts for ON/OFF dwell-time fits, sharing `dwelltime_staging_common.jl`.
Copy the folder next to your other staging scripts (e.g. `scripts/`).

| script | data | models | runs (default) |
|---|---|---|---|
| `08_stage_dwelltime_types1.jl` | MYC (gene) + CCAT1 (enhancer), 3Prime + 5Prime | types 1 | 2 folders x 90 = 180 |
| `09_stage_dwelltime_types2_ccat1.jl` | CCAT1, 3Prime + 5Prime | types 2 | 1 folder x 60 = 60 |
| `10_stage_crispr_types1.jl` | CRISPR UNT/VPR/shKRAB/KRAB, MYC + CCAT1, 5Prime | types 1 | 7 units x 45 = 315 |
| `11_stage_crispr_types2_ccat1.jl` | CRISPR, CCAT1 (KRAB has none), 5Prime | types 2 | 3 units x 30 = 90 |

All scripts use only the model-based ON/OFF times (`dwelltimes/model/...`), not the raw
segmentation. `SG_SOURCES=model,segm` brings segm back if ever needed.

## Models (types 1: 45, types 2: 30 per data unit)

Names are `<G><R><S><I><topology>`, with G, R, S, insertstep exactly as passed to StochasticGene.

**Types 1**: R/I in 301 302 303 401 403 404 501 504 505, each as
`2RSInone`, `3RSInone`, `3RSIcyc` (fully reversible), `3RSIkp` (kinetic proofreading) and
`2RSIgior`.

**Types 2** (splicing), one splice step or all steps:

| base | one step | all steps |
|---|---|---|
| x301 | S=1 (`x311`) | S=3 (`x331`) |
| x401 | S=1 (`x411`) | S=4 (`x441`) |
| x501 | S=1 (`x511`) | S=5 (`x551`) |

x303, x404 and x505 are **not** staged with splicing. `fit()` resets `S` to at most
`R - insertstep + 1` (`reset_S`), while the matrix builder splices only at steps `j <= S` (and at
the last step only if `S == R`). With insertstep = R no setting reaches the last step, so the
splice rate would be a dead parameter (the loglikelihood of `2333none` is identical for splice
rates 0.001, 1 and 10). Fixing this needs a change in StochasticGene.

**`gior`** (double-regime model) is staged as a coupled model: a hidden 2-state regime unit
(G=2, R=0; k_f, k_b) and the observed promoter unit (G=2 with the R/S/I of the name). In the
high regime OFF->ON is multiplied by (1 + gamma). ON/OFF times are those of the promoter unit.
Rate vector: `k_f, k_b, Eject1, Decay1 | promoter rates | gamma`; Eject1/Decay1 and the promoter
Decay are not fitted. Priors: G rates 0.01 (CV 1), gamma 0.1 (CV 10).

**gior runs on an old StochasticGene (patch).** Released StochasticGene (checked 2.0.4 and
2.0.7) cannot fit coupled models to dwell-time data: `fit()` crashes with
`BoundsError ... index [[16]]`. `stochasticgene_coupled_dwell_patch.jl` fixes this at run time
without reinstalling anything. Staging copies it into each swarm folder and adds
`@everywhere include(".../stochasticgene_coupled_dwell_patch.jl")` to the gior fit scripts only,
so every worker loads it. Each worker log then shows
`Applied coupled dwell-time patch to StochasticGene <version>`. On a StochasticGene that already
has the fix, the patch does nothing. Tested on 2.0.7 with `julia -p 2`. Set `SG_GIOR_PATCH=false`
to leave the fit scripts untouched. Non-gior runs never load the patch.

## Run length

`SG_MAXTIME=43200` (12 h of sampling) is the stop condition. `SG_SAMPLESTEPS=10^12` is never
reached: MH grows its sample storage as it goes. The previous gene runs stopped at their
2,000,000-step limit (`mcmc_steps = 3.2e7 = 16 x 2e6`), not at maxtime.
Samples are stored every `SG_SAMPLE_STRIDE=1000` steps. The printed submit line pads the wall
time by `SG_WALL_PAD_HOURS=2`.

## Usage

```bash
export SG_ROOT=/data/Larsonlab/Nadya/StochasticGene
export SG_PROJECT=/data/$USER/sg201          # your existing environment (gior runs load the patch)
export SG_TAG=$(date +%F)-12h
julia --project="$SG_PROJECT" -e 'using Pkg; Pkg.precompile()'

julia --project="$SG_PROJECT" scripts/08_stage_dwelltime_types1.jl
julia --project="$SG_PROJECT" scripts/09_stage_dwelltime_types2_ccat1.jl
julia --project="$SG_PROJECT" scripts/10_stage_crispr_types1.jl
julia --project="$SG_PROJECT" scripts/11_stage_crispr_types2_ccat1.jl
```

Each script prints one `submit` line per swarm folder, e.g.

```bash
cd $SG_ROOT/swarm-dwelltime-model-gene-types1-$SG_TAG && swarm -f fit.swarm -g 32 -t 18 --time 14:00:00 --module julialang --logdir logs
```

Useful filters: `SG_MODELS=3301kp,2301gior`, `SG_TOPOLOGIES=none,kp`,
`SG_DATACONDS=gene`, `SG_SAMPLES_LIST=5Prime`, `SG_CRISPR_CONDITIONS=UNT,VPR`.
Prior settings: `SG_ELONGATION_GENE=3.2`, `SG_ELONGATION_ENHANCER=22`, `SG_PRIORCV_R=0.25`,
`SG_PRIORCV_G=1`, `SG_SPLICE_PRIOR=0.05`, `SG_PRIORCV_SPLICE=0.2`, `SG_GIOR_GAMMA=0.1`,
`SG_GIOR_GAMMA_CV=10`.
