# Coupled CCAT1–MYC model sweep

E = CCAT1 (enhancer), G = MYC (gene), H = hidden unit.

| File | What it is |
|---|---|
| `config.toml` | The one file you edit: input files, sweep, simulation and swarm settings |
| `0_experimental.R` | Computes the experimental curves and ON fractions from your predicted traces |
| `1_generate_models.py` | Lists every model and every parameter set |
| `2_model_report.py` | Draws the models and parameters as an HTML page, to check before simulating |
| `3_run_simulations.jl` | Runs the simulations (either simulator) and writes the swarm file |
| `4_plot_results.py` | Writes the results report (HTML), every plot as PDF and PNG, and the ranking |
| `5_nearest_simulations.py` | Optional: the simulations of one architecture nearest a given fold or StochasticGene coupling γ (fold = 1 + γ), e.g. the value a direct fit found |
| `6_trace_relationships.py` | Optional: how CCAT1 and MYC relate in decoded traces (burst-aligned averages, model-free fold of every architecture, between- vs within-cell correlation), no model fit |
| `paper_figures.R` | The paper figures: shape map, folds in the data against the simulations, heatmaps of the two selections, ranking and candidates |
| `7_compare_folds.py` | Optional: the folds the shape-passing simulations need against the folds in the decoded traces, and `fit_candidates.csv`: the architectures to fit, with starting values |

Needs Python 3.11+ with numpy, pandas and matplotlib, Julia, and R for step 0. StochasticGene is only
needed for `simulator = "stochasticgene"`.

## Short path for the paper figures

Three commands; everything else in this folder is optional.

1. Generate the model panel (`config_v5_semi.toml`: one CCAT1 or MYC state, 1, 2, 3, R or 3+R, regulating one set of
   rates in the other gene; the same in both directions for reciprocal models; hidden units regulating the same rates
   in both genes; every rate its own fold, 0.01 to 100; hidden units activate: 1 to 100 for R12, R23 and Initiate,
   0.01 to 1 for R21 and R32; 100 parameter sets per architecture):
   `python3 1_generate_models.py config_v5_semi.toml`
2. Simulate (step 5 and 6 below write and submit the swarm). With `insert_steps = ["11", "13", "33"]` every simulation
   is analysed three ways and written to `results/insert11`, `results/insert13`, `results/insert33`.
3. Figures, once per insert variant:

```bash
Rscript paper_figures.R runs/semi_extensive_insert_steps 11 data/experimental.csv CCAT1_predictedtraces.csv MYC_predictedtraces.csv figures_insert11
```

Earlier runs without insert steps: give `none` as the insert, e.g. `Rscript paper_figures.R runs/grid none ...`;
expanded runs: give `<run_root>/<profile>` as the run folder.

`paper_figures.R` (R with ggplot2; settings at its top) writes, in order: `1_shape_map` (every simulation's
cross-correlation centre, within +/-5 min, against its flanks, 20 to 60 min, coloured by NRMSE; shaded = Mexican-hat
shape = selection 1), `2a_folds_data_vs_permuted` (for every CCAT1-MYC connection, the fold in each experimental
track, estimated with the model's form, rate x (1 + (fold - 1) x n) with n = number of reporters for R and
1 + reporters in state 3 for 3+R, blue, against tracks paired with another cell, 1,000 permutations, grey),
`2b_folds_data_vs_simulations` (the same blue violins with the folds of the hat-shaped simulations as points; row
labels give the number of tracks and how many simulations use the connection and are hat-shaped), `3_heatmap_selection1`, `4_heatmap_selection2` and `4b_folds_data_vs_selection2` (as 2b, points = selection 2 only; selection 1 with every
fold inside the middle 50% of the experimental tracks' folds; hidden-unit models have no observable fold and stay),
and `5_ranking` with `5_candidates_to_fit.csv` (best parameter set per architecture of selection 2, ranked by NRMSE;
candidates = top 25%), and `6_selection2_architectures_insert_steps` (every architecture of selection 2, in ranking order, 24 per
page: experimental cross-correlation and the simulated one, mean +/- SE, of their best parameter set with insert
steps 11, 13 and 33 in three colours). Every selection is also saved as CSV.

## Two kinds of run

- **Regular**: `config.toml` and `1_generate_models.py`, as in the steps below.
- **Expanded**: `config_expanded.toml` and `1_generate_models_expanded.py`: subsets of source states (for example
  1+2: regulates in either state), subsets of 1 to 4 regulated rates, and separate folds for every connection,
  sampled over ranges. Generate with
  `python3 1_generate_models_expanded.py config_expanded.toml --profile limited` (or `full`; add `--count-only` to
  see the size first). Then run steps 3 to 7 below with `config_expanded.toml --profile limited` (or `full`) in
  place of `config.toml`; `--profile` goes after the config file, for example
  `julia -t 12 3_run_simulations.jl config_expanded.toml --profile full --swarm`.
  Expanded runs save each simulation's curves as mean and SE only, in compact binary files
  (`results/curves_*.f32`, rows listed in `results/curve_layout.csv`); ON fractions of every replicate are in the
  metrics files as usual.

## Steps

**0. Compute the experimental data** (once). Arguments: your correlation R file, the CCAT1 and MYC
predicted-trace files, the output folder, bootstrap replicates, epsilon, frame interval in minutes.

```bash
Rscript 0_experimental.R fcs_multitau_xcor_final.R CCAT1_predictedtraces.csv MYC_predictedtraces.csv data 200 1e-7 1.6666666666666667
```

**1. Edit `config.toml`.**

- `[inputs]`: the two fitted-rate files, and `experimental.csv` and `experimental_metrics.csv` from step 0.
- `[output] run_dir`: the folder for this run.
- `[sweep] mode`: `"grid"` or `"extensive"`.
- `[simulation] simulator`: `"fast"` or `"stochasticgene"`.
- `[report]`: which curves the report draws, e.g. `observables = ["reporter"]` and `curves = ["EG"]` to skip the
  ON/OFF-state and autocorrelation plots (leave the section out to draw everything).
- `[swarm] julia_project`: your StochasticGene project folder (only for `"stochasticgene"`).

**2. Generate the models.**

```bash
python3 1_generate_models.py config.toml
```

**3. Check what will be simulated.** Open `<run_dir>/models_report.html`.

```bash
python3 2_model_report.py config.toml
```

**4. Test one job** (simulations 1 to 12) before submitting everything.

```bash
julia -t 12 --project=$SG_PROJECT 3_run_simulations.jl config.toml 1 12
```

**5. Write the swarm file.** This prints the `swarm` command to submit it.

```bash
julia 3_run_simulations.jl config.toml --swarm
```

**6. Submit the swarm** with the printed command. It looks like this:

```bash
swarm -f runs/grid/simulations.swarm -t 12 -g 8 --time 08:00:00 --module julia --logdir runs/grid/swarm_logs
```

**7. Report and rank.** Open `<run_dir>/results_report.html`.

```bash
python3 4_plot_results.py config.toml
```

For the extensive sweep, set `mode = "extensive"` and a new `run_dir` in `config.toml`, then repeat from step 2.

**8. (Optional) Compare with a fitted coupling.** The simulations of one architecture nearest a coupling value,
for example γ = −0.062 from a StochasticGene fit, drawn against the experiment and the uncoupled control:

```bash
python3 5_nearest_simulations.py config.toml E_to_G__state2__R23 --gamma -0.062
```

Add `--param fold_EG_R23` when the architecture has more than one fold column, `--n 10` for more simulations,
`--observable binary` for the ON/OFF state. Writes `nearest_*.csv`, `.pdf` and `.png` into the run folder.

**9. (Optional) Relationships in the decoded traces, without a coupled fit.** From the two `predictedtraces_*.csv`
files of the uncoupled fits:

```bash
python3 6_trace_relationships.py CCAT1_predictedtraces.csv MYC_predictedtraces.csv relationships
```

Writes `relationships/relationships_report.html` and its tables:
- activity over the recording: drifts break every stationary model
- burst-aligned averages of each gene around the other's burst onsets, in every window width (no cutoff)
- the rate of every MYC (and CCAT1) transition while the other gene is in state 1, 2, 3 or has reporters, divided by
  the rate otherwise: a model-free fold for every architecture, from a Poisson regression with a separate level for
  every cell and a smooth time course shared by all cells (so differences between cells and drifts do not count),
  plus a check table of the same fold by other methods
- per-cell activity of the two genes across cells
- the multi-tau cross-correlation of `fcs_multitau_xcor_final.R` both around each trace's own means (Coulon, as in
  step 0) and around the global means; their difference is the part that comes from differences between cells

The results report (step 7) also shows, for every single-connection architecture, its cross-correlation at every
fold tested, coloured from inhibition (blue) to activation (red), in a grid of source states by regulated rates.

**10. (Optional) Which architectures to fit.** After step 7 (it writes `shape_passing_reporter.csv`) and step 9:

```bash
python3 7_compare_folds.py config.toml relationships
```

For every connection component (for example CCAT1 state 2 -> MYC R23) it plots the fold of every simulation inside the
shape selection against the fold in the decoded traces, and writes `fit_candidates.csv`: every architecture with a
simulation inside the selection, whether one of them is consistent with the data, and the parameters of its best
simulation as starting values (StochasticGene coupling gamma = fold - 1).

**11. (Optional) Check simulations the way the data are checked.** Folds measured from frames are pulled towards 1
even from true states (CCAT1 changes state within a frame), and hidden-unit models have no fold to compare. So
export the best simulation of each candidate and analyse it exactly like the data:

```bash
julia 3_run_simulations.jl config.toml --export-traces runs/grid/fit_candidates.csv
```

```bash
for d in runs/grid/exported_traces/sim_*; do python3 6_trace_relationships.py $d/CCAT1_predictedtraces.csv $d/MYC_predictedtraces.csv $d/relationships --reference relationships; done
```

Each `relationships_report.html` then shows the simulation (colours) and the data (black): burst-aligned averages,
rate ratios, and the correlations within and between cells. `--export-traces` also takes sim_ids, e.g. `12,57`.

## What ends up in the run folder

| File | Written by | Contents |
|---|---|---|
| `models.csv` | step 2 | One row per model: what regulates what, and the parameter values tested |
| `simulations.csv` | step 2 | One row per simulation: a model plus one parameter set. The simulator's work list |
| `models_report.html` | step 3 | Schematic, regulated rates and parameter values of every model group |
| `simulations.swarm` | step 5 | One line per job: 12 simulations, each with all its replicates, on 12 threads |
| `results/curves_*.csv` | steps 4, 6 | The six correlation curves of every simulation: mean, SE, and every replicate (`rep1`, `rep2`, ...) |
| `results/metrics_*.csv` | steps 4, 6 | One row per simulation: ON fractions and mean reporters (mean, SE and every replicate), NRMSE and CCC of the six curves, base rates |
| `results_report.html` | step 7 | Interactive explorer (pick a model, vary one parameter), summary tables, and every plot. Needs `plots/` and `explorer/` next to it |
| (in the report) | step 7 | A selector at the top switches every fit-based plot and table between NRMSE, CCC and the mean NRMSE of the three curves; a second one shows all architectures or the top 10, 100 or 200. One row per architecture with its best parameter set, rows grouped by model group and sub-type and sorted by fit, with a sub-type bar and a fit bar. A shape map: every simulation's cross-correlation centre against its flanks, per group, and heatmaps of every architecture that passes the shape thresholds (labelled, 40 per page) |
| `plots/` | step 7 | Every plot as a PDF and a 600 dpi PNG |
| `explorer/` | step 7 | The data the report's explorer loads when you pick a model: one small file per model (mean and SE of its curves) |
| `ranked_scores.csv` | step 7 | Every simulation with parameters, metrics and score, best first |

Each swarm job writes one `curves_` and one `metrics_` file, named by its first and last simulation.

## The models

| Group | Coupling | Parameters swept |
|---|---|---|
| `uncoupled` | none (control) | none |
| `E_to_G` | CCAT1 regulates MYC from one source state: promoter state 1, 2, 3, or R | fold_EG |
| `G_to_E` | MYC regulates CCAT1 from one source state: promoter state 1, 2, 3, or R | fold_GE |
| `reciprocal` | both directions at once | fold_EG, fold_GE |
| `H2` | hidden state 2 regulates both genes | fold_HE, fold_HG, h12, h21 |
| `H3_symmetric` | hidden state 1 regulates CCAT1, state 3 regulates MYC, state 2 nothing; h32 = h12, h23 = h21 | fold_HE, fold_HG, h12, h21 |
| `H3_asymmetric` | same; four independent rates | fold_HE, fold_HG, h12, h21, h23, h32 |

- Regulated rates: each connection multiplies one of nine rate sets by its fold: R12, R21, R23, R32,
  Initiate, R12+R23, R21+R32, R12+R23+Initiate, R21+R32+Initiate. Reciprocal models take every
  pairing of the two directions' rate sets (81).
- Fold: a regulated rate is multiplied by 1 + (fold − 1) for every active coupling effect. A fold above 1
  activates, below 1 inhibits, 1 does nothing. Each kind of connection has its own fold column in
  `simulations.csv`: `fold_EG` (E regulates G), `fold_GE`, `fold_HE`, `fold_HG`.
- Gene-to-gene models (`E_to_G`, `G_to_E`, `reciprocal`) are run with every fold and its inverse.
  Hidden models activate ([sweep] hidden_effect = "activating", the default): they speed R12, R23 and Initiate
  (fold above 1) and slow R21 and R32 (fold below 1). "both" allows either direction; "increase" puts every hidden
  fold above 1, also for R21 and R32 (the behaviour of runs made before hidden_effect existed).
- Reciprocal models, `[sweep] reciprocal_states`: `"same"` uses the same source state and the same fold in
  both directions; `"all_pairs"` uses every pair of source states and every pair of folds, so one direction
  can activate while the other inhibits.
- Source state R is StochasticGene's `Rsum`: every occupied reporter position of the source gene adds one
  coupling effect, so the rate is multiplied by 1 + (fold − 1) × (number of occupied positions).
- Clipping: a multiplier that would go below 0.000001 (inhibition from R: fold 0.5 with three positions gives
  1 − 1.5) is held at 0.000001 (`MIN_FACTOR` at the top of `3_run_simulations.jl`), in the rate matching and in
  both simulators. For StochasticGene the script replaces the package's internal check when it starts; the
  installed package is not changed. Such simulations are marked `clipped = 1` in the metrics files.
- Controls: every group also has one model with all folds 1 (gamma = 0).
- Rate matching (`[simulation] rate_matching`): base rates are rescaled so that each regulated rate still
  averages to its fitted value. `"time_average"` averages over all time, as in the original sweep; coupling
  can then still shift a gene's ON fraction. `"state_average"` averages over the time when the gene is in
  the state where that rate can act, which keeps ON fractions close to the uncoupled fit.
- `grid`: every combination of `coupling_fold` (and its inverses, for gene-to-gene models) and `hidden_rate`.
  One fold per connection, multiplying all its regulated rates; hidden models use one fold for both connections.
- `extensive`: `n_extensive` random combinations per model. Every regulated rate has its own fold (columns
  `fold_EG_R12`, `fold_EG_R23`, ...), sampled over the range of its connection in `[sweep.extensive]`; hidden
  rates likewise. All are independent, also in `"same"` reciprocal models (there `"same"` only means the same
  source state). Without `[sweep.extensive]` the ranges come from the grid lists. The expanded generator also
  gives every regulated rate its own fold.

## The two simulators

- `fast`: lists all states of the system, starts each trace from the exact long-run distribution
  (no burn-in) and runs a Gillespie simulation.
- `stochasticgene`: `StochasticGene.simulate_trace_vector`, as in `simulate_shared_hidden_sweep.jl`.
  Each trace runs for a burn-in of `warmup_timescales` (3) / slowest rate before recording.

## Scoring

`score()` at the top of `4_plot_results.py` ranks the simulations (lower is better). It is the
NRMSE of the reporter cross-correlation by default; edit it to use any saved metric.

The report also computes the Mexican-hat shape of each cross-correlation: its mean within ±10 min
(`center_*`), its mean in the flanks from 20 to 60 min on each side (`left_*`, `right_*`), the depth
`hat_*` = centre minus the mean of the flanks, and `hat_shape_*` = true when the centre is at least
`center_min` and the mean of the two flanks is at most `flank_max` (`[report.shape]` in the config; defaults 0.01 and −0.01). The
windows are `HAT_CENTER_MIN` and `HAT_FLANKS_MIN` at the top of the script.

**Shape selection.** The report has two selections, switched by the "Shape selection" menu at the top:
`shape` (the thresholds above) and `shape_nrmse` (also NRMSE below `nrmse_max`, default 0.9; `nrmse_of` picks the
NRMSE: `"cross"` = of the cross-correlation (default), `"combined"` = the combined score, `"each"` = of each of the
three curves; `nrmse_max = ""` drops this selection). For each, the report draws one row per architecture (source
states and regulated rates) with its best passing parameter set, grouped by model group and sub-type and sorted by
the selected fit metric, and (with `labelled_pages = true`) the same with every model labelled, 40 per page.
Files in the run folder, with `<selection>` = `shape` or `shape_nrmse`:
- `<selection>_passing_<observable>.csv`: every passing parameter set
- `<selection>_selected_<observable>_by_<metric>.csv`: the best passing set of each architecture, by `nrmse`, `ccc`
  or `nrmse3`
All are columns of `ranked_scores.csv` and can be used in `score()`.

`nrmse3_reporter` is the fit of all three curves: the mean of the NRMSEs of the cross-correlation, the CCAT1
autocorrelation and the MYC autocorrelation (reporter count, equal weights; `[report] combine = "sum"` gives their
sum, which ranks the same). It is the third option of the report's fit-metric menu and ranks the lists `top_10_architectures.csv` and `top_100_architectures.csv` in the run
folder: the best architectures, each with its best parameter set (all its parameters and fits in one row).

## Notes

- The report shows the mean ± SE over the replicates. Every replicate is kept in `results/curves_*.csv`
  (columns `rep1`, `rep2`, ...) and `results/metrics_*.csv` for your own plots.
- If a simulation is in more than one results file (a test job and the swarm, say), the newest file is used.
