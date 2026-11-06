# Coupled CCAT1–MYC sweeps with StochasticGene's simulator

`simulate_coupled_sweep.jl` extends `simulate_shared_hidden_sweep.jl`. Traces are
still simulated by `StochasticGene.simulate_trace_vector`. The rate layout, the
`time_average` rate convention, the seeds, the Cartesian `[sweep]` and the
original output files are unchanged.

## Run

```bash
julia -t 16 --project=$SG_PROJECT contrib/coupled_sweep/simulate_coupled_sweep.jl CONFIG.toml --check
```

`--check` builds and validates every scenario, writes `sweep_points.csv` and
`sweep_design.csv`, and prints the burn-in range. It simulates nothing. Then:

```bash
julia -t 16 --project=$SG_PROJECT contrib/coupled_sweep/simulate_coupled_sweep.jl CONFIG.toml
```

To split a large sweep across jobs, give a scenario range. Its output files get the
suffix `_scenarios_FIRST-LAST`:

```bash
julia -t 16 --project=$SG_PROJECT contrib/coupled_sweep/simulate_coupled_sweep.jl CONFIG.toml 1 96
```

## Direct and reciprocal coupling configs

`direct_CCAT1_MYC_sweep.toml` is generated; to change it, edit `direct_base.toml`
and regenerate:

```bash
python3 contrib/coupled_sweep/make_direct_sweep.py contrib/coupled_sweep/direct_base.toml contrib/coupled_sweep/direct_CCAT1_MYC_sweep.toml
```

Default rate sets: R12, R21, R23, R32, Initiate, R12+R23, R21+R32,
R12+R23+Initiate and R21+R32+Initiate. The source gene regulates from promoter
state 3. The generated schemes are:

- 1 uncoupled control
- 9 CCAT1 → MYC
- 9 MYC → CCAT1
- 81 reciprocal (every pairing of the two directions)

With 4 folds that makes 400 scenarios, each with 3 replicate datasets
(`nsim = 3`). Useful options:

- `--source-states 1 2 3`
- `--rate-sets R12 Initiate`
- `--mode inhibit`
- `--no-one-way`

## Configuration additions

| Setting | Meaning |
|---|---|
| `[hidden] topology` | `"none"`, `"H2"` (default when `[hidden]` exists) or `"H3"` (chain 1 ⇄ 2 ⇄ 3, rates `R12 R21 R23 R32`; `symmetric = true` sets R23 = R21, R32 = R12) |
| `[[coupling]] source` | `"hidden"` (default), `"gene1"` or `"gene2"` |
| `[[coupling]] source_state` | Hidden state, or the source gene's promoter (G) state. Defaults to 2 for H2 |
| `[sweep] hidden_R23`, `hidden_R32` | Sweep the H3 rates |
| `[simulation] equilibration_timescales` | Burn-in = this × the slowest relaxation time (default 10) |
| `[simulation] trace_lengths_csv` | `experimental_lengths.csv` from `experimental.R`: simulate the experimental trace lengths (replaces `ntraces`, `nframes`) |
| `[simulation] epsilon` | Added to traces for the six multi-tau curves, as in `experimental.R` (default 1e-7) |
| `[scoring] experimental_csv`, `max_lag_min` | Score each scenario against `experimental.csv` (matched on `lag_frames`) |

Coupling rules are StochasticGene's: blocks regulating the same rate add,
rate = base × (1 + Σ active gammas), with gamma = fold − 1 (activate) or 1/fold − 1
(inhibit). A scenario whose combined factor is ≤ 0 in any promoter state stops the
run before anything is simulated.

`rate_reference = "time_average"` chooses base rates so that each coupled rate's
stationary time average equals the supplied rate, solved on the joint promoter
chain (at most 27 states). With hidden sources this is your original closed form
OFF = supplied / ((1 − p) + pF). With gene sources it is iterated, because the
source's occupancy depends on the coupling. `"on"` is only allowed when one block
regulates the rate.

## Burn-in

Each trace starts in promoter state 1 with no reporters, and StochasticGene runs it
for `warmup` minutes before recording. The original rule was
10 × max(hidden relaxation, 1 / slowest average rate). For CCAT1 that is about 9,500
minutes, because R12 = 0.00105/min. The real slowest relaxation time of the CCAT1
promoter is 401 min, and the script now uses it: computed exactly from the promoter
chain, and from the reporter traversal time. The default of 10 timescales gives
about 4,000 min with the same safety margin. 7–8 timescales leaves an
initial-state bias of about 0.1%.

## Outputs (per configuration or scenario range)

The original files are unchanged:

- `metrics.csv`
- `metrics_summary.csv` (adds mean reporter counts)
- `xcor.csv`
- `xcor_summary.csv`
- `xcor_multitau.csv` and `xcor_multitau_summary.csv` (binary gene1 × gene2)
- `sweep_points.csv`
- `baseline_binary_traces.csv`

New files:

- `sweep_design.csv`: hidden occupancies, coupling edges, input/base/time-average
  rates, relaxation times and burn-in per scenario.
- `multitau_all.csv`, `multitau_all_summary.csv`: binary and reporter-count curves,
  EE/GG/EG, all lags, computed exactly as `experimental.R` computes its targets.
- `scores.csv` (with `[scoring]`): NRMSE, Lin's CCC and the center/left/right means
  of each mean curve against the experimental target.
