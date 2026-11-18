# ON/OFF dwell-time fits

Scripts around the dwell-time staging in `$SG_ROOT/scripts` (`06_swarm_dwelltime_expanded_topology.jl`,
`06_swarm_dwelltime_enhancer_splicing.jl`).

## Restage for 12 hours

```bash
SG_PROJECT=/data/fursovan2/sg111 SG_SUBMIT=true bash stage_dwelltime_12h.sh
```

Copy the script to `$SG_ROOT/scripts` or run it from here (it changes into `$SG_ROOT`). New folders, named with
`SG_TAG` (default `<date>-12h`):

| results folder | fits |
|---|---|
| `dwelltime-model-gene-<tag>`, `dwelltime-segm-gene-<tag>` | 16 topologies x 3Prime, 5Prime (elongation prior 3.2 min) |
| `dwelltime-model-enhancer-<tag>`, `dwelltime-segm-enhancer-<tag>` | 16 topologies x 3Prime, 5Prime (22 min) |
| `dwelltime-model-enhancer-splicing-<tag>`, `dwelltime-segm-enhancer-splicing-<tag>` | 3Prime: 3301, 3302, 3303, 3331, 3322; 5Prime: 3301, 3331, 3322, 3311 (20 min) |

Each run stops at `SG_MAXTIME` = 12 hours: the step limit (10^9) is never reached, and samples are stored as they
come, so it costs no memory. Each batch sets every variable it uses, so nothing carries over between batches.

## Report: data against the fits, models ranked

```bash
julia --project=$SG_PROJECT dwelltime_report.jl dwelltime-model-enhancer-<tag> dwelltime-model-enhancer-splicing-<tag>
python3 dwelltime_report.py dwelltime-model-enhancer-<tag> dwelltime-model-enhancer-splicing-<tag> --out model_enhancer.pdf
```

The Julia step rebuilds every finished fit from its saved run spec with its maximum-likelihood rates and writes
`<folder>/dwelltime_report/fits.csv` and `distributions.csv`. The Python step draws, for every set of histograms
(source, sample, condition), the models ranked by AIC and the ON and OFF distributions of data and models.
Folders that fit the same histograms are compared together.

- Log-likelihood: sum over the bins of both dwell types of (events in the bin) x log(model probability of the bin),
  with the model probabilities of each type summing to 1 over the histogram's bins (StochasticGene's likelihood).
- AIC = 2 x (fitted parameters) - 2 x log-likelihood, as StochasticGene's `aic`. Only fits of the same histograms
  can be compared: delta AIC below 2 = as good as the best, above 10 = clearly worse.
- KS: the largest gap between the data and model cumulative distributions, per dwell type.
- Rhat above 1.1: the chains have not converged; that fit's AIC is not reliable yet.
