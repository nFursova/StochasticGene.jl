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

Needs Python 3.11+ with numpy, pandas and matplotlib, Julia, and R for step 0. StochasticGene is only
needed for `simulator = "stochasticgene"`.

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

## What ends up in the run folder

| File | Written by | Contents |
|---|---|---|
| `models.csv` | step 2 | One row per model: what regulates what, and the parameter values tested |
| `simulations.csv` | step 2 | One row per simulation: a model plus one parameter set. The simulator's work list |
| `models_report.html` | step 3 | Schematic, regulated rates and parameter values of every model group |
| `simulations.swarm` | step 5 | One line per job: 12 simulations, each with all its replicates, on 12 threads |
| `results/curves_*.csv` | steps 4, 6 | The six correlation curves of every simulation: mean, SE, and every replicate (`rep1`, `rep2`, ...) |
| `results/metrics_*.csv` | steps 4, 6 | One row per simulation: ON fractions and mean reporters (mean, SE and every replicate), NRMSE and CCC of the six curves, base rates |
| `results_report.html` | step 7 | Interactive explorer (pick a model, vary one parameter), summary tables, and every plot |
| `plots/` | step 7 | Every plot as a PDF and a 600 dpi PNG |
| `ranked_scores.csv` | step 7 | Every simulation with parameters, metrics and score, best first |

Each swarm job writes one `curves_` and one `metrics_` file, named by its first and last simulation.

## The models

| Group | Coupling | Parameters swept |
|---|---|---|
| `uncoupled` | none (control) | none |
| `E_to_G` | CCAT1 in promoter state 3 regulates MYC | fold |
| `G_to_E` | MYC in promoter state 3 regulates CCAT1 | fold |
| `reciprocal` | both, same fold | fold |
| `H2` | hidden state 2 regulates both genes | fold, h12, h21 |
| `H3_symmetric` | hidden state 1 regulates CCAT1, state 3 regulates MYC, state 2 nothing; h32 = h12, h23 = h21 | fold, h12, h21 |
| `H3_asymmetric` | same; four independent rates | fold, h12, h21, h23, h32 |

- Regulated rates: each model multiplies one of nine rate sets by the fold: R12, R21, R23, R32,
  Initiate, R12+R23, R21+R32, R12+R23+Initiate, R21+R32+Initiate. Reciprocal models take every
  pairing of the two directions (81 models).
- Controls: every group also has one model with fold 1 (gamma = 0).
- Rate matching: base rates are rescaled so that each regulated rate's long-run average equals its
  fitted value, as `rate_reference = "time_average"` did in the original sweep.
- `grid`: every combination of the listed fold and hidden-rate values.
- `extensive`: 100 random combinations per model. Each parameter takes 100 log-spaced values over
  the same range, each used once.

## The two simulators

- `fast`: lists all states of the system, starts each trace from the exact long-run distribution
  (no burn-in) and runs a Gillespie simulation.
- `stochasticgene`: `StochasticGene.simulate_trace_vector`, as in `simulate_shared_hidden_sweep.jl`.
  Each trace runs for a burn-in of `warmup_timescales` (3) / slowest rate before recording.

## Scoring

`score()` at the top of `4_plot_results.py` ranks the simulations (lower is better). It is the
NRMSE of the reporter cross-correlation by default; edit it to use any saved metric.
