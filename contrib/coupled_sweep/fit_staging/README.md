# Coupled fits of the sweep candidates

The 25 architectures chosen from the sweep, staged as coupled StochasticGene fits exactly like the earlier coupled
runs (two separate hierarchical source fits, per-trace noise frozen, the same priors), with only the architectures
new and `propcv = 0.1`. Every gamma is free and starts at 0.

| File | What it is |
|---|---|
| `models_direct.txt` | 9 CCAT1 -> MYC, MYC -> CCAT1 and reciprocal models, for your `scripts/stage_coupled_reciprocal.jl` |
| `models_hidden.txt` | 16 hidden-unit models (H2, H3 symmetric, H3 asymmetric), for `stage_coupled_hidden.jl` |
| `stage_coupled_hidden.jl` | Stager for hidden-unit models: like `stage_coupled_reciprocal.jl` plus an unobserved unit 3 |
| `stage_candidates_24h.sh` | Stages both sets, 3Prime, target mode, 24 hours |
| `continue_candidates_72h.sh` | Restages both sets from where they stopped and runs 72 hours more (`scripts/30_continue_coupled.jl`) |

Steps on Biowulf (copy this folder next to your scripts, or call the scripts by their path):

```bash
SG_PROJECT=/data/fursovan2/sg111 bash stage_candidates_24h.sh
```

Then submit the two `submit_all.sh` files it prints. When the 24-hour runs have ended:

```bash
SG_PROJECT=/data/fursovan2/sg111 bash continue_candidates_72h.sh
```

Hidden-unit models (StochasticGene's hidden-latent layout, `unit_model = (1, 2, 3)`): unit 3 has no reporter and
no noise parameters. H2: hidden states 1 <-> 2, state 2 regulates the given rates of both genes. H3: hidden states
1 <-> 2 <-> 3, state 1 regulates CCAT1 and state 3 MYC; H3s ties h23 = h21, h32 = h12 and the two genes' gammas of
each rate; H3a leaves them separate. Hidden rates start at `SG_HIDDEN_RATE0`, which is also their prior mean, and
are fitted (`SG_FIT_HIDDEN=1`) with prior cv `SG_PRIORCV_HIDDEN` (stager defaults 0.1 per minute and 10;
stage_candidates_24h.sh uses 1 per minute and 1.0, the cv of the G rates in the individual fits). In target mode the fitted parameters are the regulated
rates of both genes, the hidden rates and the gammas.
