# Review of the delivered coupled-model screen

This review covers the delivered screen:

- `workflow.py`, `prepare.py`
- `simulate_model_sweep.jl`, `sweep_multitau_helpers.jl`, `run_batch.jl`
- `experimental.R`, `render_report.py`, `merge_runs.py`, `report_viewer.js`

They were compared with the original sweep: `simulate_shared_hidden_sweep.jl`, its
TOML and `hidden_unit_sweep.pl`. Files not seen:

- `simulate_shared_hidden_sweep_looping.jl`
- `fcs_multitau_xcor_final.R`
- `plot_shared_hidden_sweep_v2.R`
- the delivered `sweep.toml`

## Delivered screen compared with the original sweep

| | Original | Delivered |
|---|---|---|
| Simulator | `StochasticGene.simulate_trace_vector`, with burn-in | Its own joint-chain Gillespie simulator; StochasticGene is not used |
| Coupling | Hidden unit → both genes | Adds CCAT1 → MYC, MYC → CCAT1, reciprocal and a 3-state hidden unit (H3) |
| Fold | One `coupling_fold` for both genes, 1–100,000, activate or inhibit | Separate fold per connection, 0.01–100 |
| Hidden rates | 0.01–10 /min | 0.001–1000 /min |
| Same rate regulated twice | Error | Effects add, clipped at 1e-12. This happens in every hidden + direct design, because those designs reuse one rate set for all their connections |
| Sweep | Cartesian, 384 scenarios | ≥ 239,401 stratified points over 2,395 architectures (1,764 of them reciprocal) |
| `time_average` | Closed form OFF = supplied / ((1−p) + pF) | Joint promoter-chain solve, identical for hidden sources |
| Outputs | Binary multi-tau, Pearson, connected correlation, burst-timing %, joint ON | Six multi-tau curves, ON fractions. Burst timing, Pearson/connected and joint ON were dropped |

### Verified correct

- The multi-tau code is copied verbatim. Turning background correction off and
  adding ε matches the original on 0/1 traces: the background step is a no-op
  unless a trace is ON in more than 99.9% of frames.
- Lag sign: positive lag means MYC follows CCAT1.
- The reporter logic of the delivered simulator matches G3 R3 S0 insertstep 1.
- The hidden-state occupancy and dwell-time formulas in the viewer are correct.

### Problems found

1. Simulated traces all have 216 frames, but the experimental traces have variable
   lengths (from 17 frames). The per-trace mean subtraction biases correlations by
   roughly 1/(trace length).
2. `coulon_multitau_mean` keeps only lags shared by every trace. If the R
   `coulon_average` does the same, a 17-frame trace limits the experimental target
   to ±8 frames (±13 min). Check this in `fcs_multitau_xcor_final.R`.
3. `merge_runs.py` always reports `missing=[]`, so a chunk that was never included
   goes unnoticed.
4. One bad point (calibration failure, truncated file) makes `render_report.py`
   refuse to render the whole report.
5. The viewer matches lags in minutes rounded to 1e-7, while the Julia side matches
   on `lag_frames`. If the archived figure file stores lags with fewer digits, its
   scores become unavailable.
6. `read_targets` crashes on `NA` in `experimental.csv`.
7. A killed job leaves stale `.lock` and `.partial.*` folders behind, and the
   point fails on every rerun until they are deleted by hand.
8. `merge_runs.py` re-renders reports inside the source runs; hard links tie the
   merged copy to the originals.
9. Every architecture gets the same number of points, so high-dimensional families
   are sampled sparsely and reciprocal designs take 74% of the compute.
10. The search box re-scores every model on each keystroke.

## StochasticGene simulator compared with the delivered simulator

| | StochasticGene | Delivered |
|---|---|---|
| Method | Next-reaction method: one clock per reaction; coupling rescales the target's clock when its source changes state | Lists all joint states (≤ 1,728) and their transitions, then runs a Gillespie direct-method loop |
| Start | State 1, then burn-in | Exact stationary distribution, no burn-in |
| Coupling | Additive; error if the factor is ≤ 0 | Additive; silently clipped |
| Scope | Any G, R, S, insertstep, alleles, noise | Only G3 R3 S0 insertstep 1 with two genes |
| Validation | Package code, same rules as the fits | Never cross-checked against StochasticGene |
| Speed | Slower; the burn-in dominates | Much faster per point |

Both are exact stochastic simulations of the same model. StochasticGene is the
better reference, because it is general, tested and consistent with the fits.
The delivered simulator is faster for this one model, but should be validated
before use. Validation: compare ON fractions and mean reporter counts with the
exact stationary values, and compare multi-tau curves, on about 5 representative
points.

StochasticGene already supports every architecture the screen added, through
`connections` from any unit to any unit. A custom simulator was not needed.

## What this folder does instead

`simulate_coupled_sweep.jl` extends the original sweep on StochasticGene's
simulator; see README.md.

Burn-in now uses the exact relaxation time:

- CCAT1 promoter: 401 min, against 950 min from the slowest-rate rule.
- Reciprocal coupling at fold 1,000: about 980 min.

Python checks of the rate-matching maths (not the Julia code itself):

- With hidden sources it reproduces the original closed form to 12 digits, up to
  fold 100,000.
- Reciprocal coupling converges to a relative error of 1e-13.

Reciprocal designs are generated from 9 rate sets with 3 replicate datasets
(`make_direct_sweep.py`).

Not yet run in Julia. Before a full sweep:

1. Run `--check`.
2. Run one or two scenarios of the original H2 configuration and compare them
   with earlier results.
