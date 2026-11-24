# Separate 3Prime / 5Prime trace fits with per-trace noise (12, 13)

Two-stage workflow. Copy both scripts next to your other staging scripts (`scripts/`).

## Stage A: `12_stage_trace_separate_noise.jl`, fit the noise (12 h)

- Each sample (3Prime, 5Prime) is its own fit: nothing is shared between them.
- Classic hierarchical path, `hierarchical = (2, per-trace params, ())`:
  - noise offsets 1 and 3 (background mean, reporter mean): fitted **per trace**
  - noise offsets 2 and 4 (background SD, reporter SD): one value **per sample**
  - kinetic rates: one value per sample
- Ties: all Rshift rates equal (Rshift_1 fitted); Eject separate; all splice rates equal.
- Priors: G 0.01 (CV 1), Initiate 0.2 (CV 0.25), Rshift 1/2.5 min (CV 0.25), Eject 1.0 /min
  (CV 0.25), splice 0.05 (CV 0.2), noise = the per-sample stage-05 priors (CV 0.1).
- Models: bases `3301 3302 3303 3401 3402 3403 3404 3331 3441` x topologies `none cyc kp`
  = 27 per sample, 54 runs. `cyc` = linear + 3<->1, `kp` = 1<->2, 2->3, 3->1.
- Splicing only with insertstep 1: with insertstep > 1, `fit()` resets S and the splice rates
  never enter the likelihood. The script refuses such models with an explanation.
- `SG_MAXTIME=43200` is the stop condition; `SG_SAMPLESTEPS=10^12` is never reached.

```bash
export SG_ROOT=/data/Larsonlab/Nadya/StochasticGene
export SG_PROJECT=/data/$USER/sg201
export SG_TAG=$(date +%F)
julia --project="$SG_PROJECT" scripts/12_stage_trace_separate_noise.jl
# prints: cd <swarm dir> && swarm -f fit.swarm -g 32 -t 18 --time 15:00:00 --module julialang --logdir logs
```

## Stage B: `13_continue_trace_fixnoise.jl`, freeze the noise, fit the kinetics

For every finished stage-A run, from **its own** fit: all noise parameters (per trace and per
sample) are frozen at the `SG_RATETYPE` row (median by default), written into all four rows of a
new starting-rate file; the kinetic rates keep their priors and ties and start from the same row.

```bash
export SG_FROM=trace-separate-enhancer-noise-$SG_TAG
julia --project="$SG_PROJECT" scripts/13_continue_trace_fixnoise.jl     # SG_INTO defaults to $SG_FROM-fixnoise-24h
```

## Options

`SG_SAMPLES_LIST=3Prime,5Prime`, `SG_CONDITION=enhancer` (or gene), `SG_DATA_3PRIME`,
`SG_DATA_5PRIME`, `SG_BASES`, `SG_TOPOLOGIES=none,cyc,kp`, `SG_MODELS` (full names, overrides),
`SG_RSHIFT_MINUTES=2.5`, `SG_PRIORCV_R=0.25`, `SG_EJECT_RATE=1.0`, `SG_PRIORCV_EJECT=0.25`,
`SG_NOISE_INDIVIDUAL=1,3`, `SG_NOISE_SHARED=2,4`, `SG_PROPCV_RATE=0.05`, `SG_PROPCV_NOISE=0.001`,
stage B: `SG_INTO`, `SG_RATETYPE=median`, `SG_MAXTIME=86400`, `SG_ONLY`, `SG_FIXNOISE_PROPCV=0.05`.

Tested on synthetic traces (24 + 20 traces): stage A recovers per-trace reporter means with the
SDs shared and ties held; stage B leaves every noise value identical to stage A and refits the
kinetic rates with the ties still held (including `cyc`, whose extra transitions shift the
indices).
