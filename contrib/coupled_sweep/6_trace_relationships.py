#!/usr/bin/env python3
"""Optional. How CCAT1 and MYC relate in decoded traces, without fitting a coupled model.

    python3 6_trace_relationships.py CCAT1_predictedtraces.csv MYC_predictedtraces.csv OUTPUT_FOLDER
            [--interval 1.6666666666666667] [--bootstrap 1000] [--max-lag-min 120] [--epsilon 1e-7]

Inputs: the two predictedtraces_*.csv files of StochasticGene (columns Gstate_<i>, Rstate_<i>, Reporters_<i> per
trace; trace i of one file is the same cell as trace i of the other). Any traces in this format work, for example
simulated ones.

Writes into OUTPUT_FOLDER: relationships_report.html, plots/ (PDF and PNG of every plot) and these tables:

1. Burst-aligned averages (burst_aligned.csv, burst_aligned_windows.csv)
   Reference events: burst onsets of one gene (decoded reporter count goes from 0 to above 0). Response: the other
   gene's (and the same gene's) burst onsets, time in promoter state 3, and reporter count, at every lag before and
   after. "Change" = response at that lag / its expected value - 1. The expected value is the cell's own activity
   level times the population's time course at that frame (smoothed over --trend-window frames), so neither cells
   that are simply more active nor a drift over the recording counts as an effect (--trend-window 0: no time
   course, the cell's average only). "Level" = how much more active the response gene is, on average, in cells
   with more reference bursts (between-cell effect). No window cutoff: the windows table and plot give the average
   change within every window from 1 frame to the maximum lag, before, after and around the onset.

2. Rate ratios for every architecture (rate_ratios.csv)
   For every target rate (R12, R21, R23, R32 from the decoded promoter states; Initiate = first reporter position
   becoming occupied while in state 3) and every source condition (other gene in state 1, 2, 3, or with reporters
   = R): the rate while the condition holds divided by the rate while it does not. This is a model-free estimate
   of that architecture's fold. Within-cell (Mantel-Haenszel over strata = each cell's trace cut into
   --time-blocks parts: compares times inside the same cell and the same part of the recording) and pooled (all
   cells together; includes between-cell differences and drifts). 95% intervals from bootstrapping cells.
   Source sets of two states are the reciprocal of the third (state 1+2 = not 3).

3. Between cells versus within cells (cell_statistics.csv, cell_correlations.csv, correlations.csv)
   Per-cell activity of each gene and their rank correlation across cells. The multi-tau correlation of
   fcs_multitau_xcor_final.R (Coulon average: each trace around its own mean = within cells) and the same
   correlation around the global means (pooled over all cells = within + between cells). Their difference is the
   between-cell part.

0. Stationarity (activity_over_time.csv): each gene's activity against time in the recording. A drift breaks the
   assumption of every stationary model (the sweep and StochasticGene fits) and shapes the correlations.
"""
import argparse, html, re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GENES = {"E": "CCAT1", "G": "MYC"}
INK, MUTED, BLUE, RED = "#0b0b0b", "#898781", "#2a78d6", "#e34948"
RATES = [("R12", 1, 2), ("R21", 2, 1), ("R23", 2, 3), ("R32", 3, 2), ("Initiate", None, None)]
SOURCES = ["1", "2", "3", "R"]

# ---------------------------------------------------------------------------
# Traces
# ---------------------------------------------------------------------------

def read_traces(path):
    """{trace id: dict(state=promoter state, first=first reporter position occupied, reporters=count)}."""
    frame = pd.read_csv(path, na_values=["", "NA", "NaN", "missing"], low_memory=False)
    ids = sorted(int(m.group(1)) for c in frame.columns if (m := re.fullmatch(r"Reporters_(\d+)", c)))
    traces = {}
    for i in ids:
        reporters = frame[f"Reporters_{i}"].to_numpy(dtype=float)
        n = int(np.isfinite(reporters).sum())
        positions = [[int(v) for v in re.findall(r"\d+", str(text))] for text in frame[f"Rstate_{i}"].iloc[:n]]
        traces[i] = dict(state=frame[f"Gstate_{i}"].to_numpy(dtype=float)[:n].astype(int),
                         first=np.array([p[0] if p else 0 for p in positions]), reporters=reporters[:n].astype(int))
    return traces


def onsets(reporters):
    """1 at frames where the reporter count goes from 0 to above 0."""
    out = np.zeros(len(reporters))
    out[1:] = (reporters[:-1] == 0) & (reporters[1:] > 0)
    return out


def responses(trace):
    return {"burst onset": onsets(trace["reporters"]), "in state 3": (trace["state"] == 3).astype(float),
            "reporter count": trace["reporters"].astype(float)}


def resample(n_cells, reps, rng):
    """Bootstrap weights: how often each cell is drawn, one row per replicate."""
    return np.stack([np.bincount(rng.integers(0, n_cells, n_cells), minlength=n_cells) for _ in range(reps)]).astype(float)


def interval(values):
    return np.nanpercentile(values, [2.5, 97.5], axis=0)

# ---------------------------------------------------------------------------
# 1. Burst-aligned averages
# ---------------------------------------------------------------------------

def time_course(cells, gene, kind, window):
    """The population's response against frame, relative to its mean, smoothed over `window` frames (1 if window = 0)."""
    longest = max(len(cell[gene]["state"]) for cell in cells)
    total, count = np.zeros(longest), np.zeros(longest)
    for cell in cells:
        r = responses(cell[gene])[kind]
        total[:len(r)] += r
        count[:len(r)] += 1
    if window <= 1:
        return np.ones(longest)
    kernel = np.ones(window)
    smooth = np.convolve(total, kernel, "same") / np.maximum(np.convolve(count, kernel, "same"), 1)
    return smooth / (total.sum() / count.sum()) if total.sum() > 0 else np.ones(longest)


def burst_aligned(cells, reference, response_gene, K, weights, window):
    """Per response kind: lag, change (within-cell, with 95% interval) and level (between-cell).
    Expected response of cell c at frame f = a_c * course(f), with a_c = the cell's total / its sum of course."""
    lags = np.arange(-K, K + 1)
    out = {}
    for kind in responses(cells[0][response_gene]):
        course = time_course(cells, response_gene, kind, window)
        S = np.zeros((len(cells), len(lags))); N = np.zeros_like(S); B = np.zeros_like(S)
        total, expected, events_c, level_c = (np.zeros(len(cells)) for _ in range(4))
        for c, cell in enumerate(cells):
            r = responses(cell[response_gene])[kind]
            events = np.flatnonzero(onsets(cell[reference]["reporters"]))
            a = r.sum() / course[:len(r)].sum()
            total[c], expected[c], events_c[c], level_c[c] = r.sum(), course[:len(r)].sum(), len(events), len(events) * a
            for j, lag in enumerate(lags):
                at = events + lag
                at = at[(at >= 0) & (at < len(r))]
                S[c, j], N[c, j], B[c, j] = r[at].sum(), len(at), a * course[at].sum()
        if reference == response_gene and kind == "burst onset":
            S[:, K] = np.nan                                  # an onset coincides with itself
        def stats(w):
            return (w @ S) / (w @ B) - 1, (w @ level_c) / (w @ events_c) / ((w @ total) / (w @ expected)) - 1
        change, level = stats(np.ones(len(cells)))
        boot_change, boot_level = stats(weights)
        windows = []
        for w_ in range(1, K + 1):
            for side, keep in (("before", (lags >= -w_) & (lags < 0)), ("after", (lags > 0) & (lags <= w_)),
                               ("around", np.abs(lags) <= w_)):
                s, b = np.nansum(S[:, keep], axis=1), np.nansum(np.where(np.isnan(S[:, keep]), 0, B[:, keep]), axis=1)
                value, boot = s.sum() / b.sum() - 1, (weights @ s) / (weights @ b) - 1
                windows.append(dict(window_frames=w_, side=side, change=value, low=interval(boot)[0], high=interval(boot)[1]))
        out[kind] = dict(lags=lags, change=change, low=interval(boot_change)[0], high=interval(boot_change)[1],
                         level=level, level_low=interval(boot_level)[0], level_high=interval(boot_level)[1],
                         events=int(np.nansum(N[:, K])), windows=pd.DataFrame(windows))
    return out

# ---------------------------------------------------------------------------
# 2. Rate ratios
# ---------------------------------------------------------------------------

def rate_counts(cells, source, target, blocks):
    """Per stratum (cell and part of its trace), condition and rate: events and frames at risk, while the source
    condition holds (1) or not (0)."""
    rows = []
    for c, cell in enumerate(cells):
      s_all, t_all = cell[source], cell[target]
      n_all = len(t_all["state"]) - 1                     # transitions from frame f to f + 1
      edges = np.linspace(0, n_all, blocks + 1).round().astype(int)
      for block in range(blocks):
        start, stop = edges[block], edges[block + 1]
        if stop <= start:
            continue
        s = {key: value[start:stop + 1] for key, value in s_all.items()}
        t = {key: value[start:stop + 1] for key, value in t_all.items()}
        n = stop - start
        conditions = {name: (s["state"][:n] == int(name)) if name != "R" else (s["reporters"][:n] > 0) for name in SOURCES}
        for rate, a, b in RATES:
            if a is None:                                 # initiation: state 3 with the first position free
                risk = (t["state"][:n] == 3) & (t["first"][:n] == 0)
                event = risk & (t["first"][1:] == 1)
            else:
                risk = t["state"][:n] == a
                event = risk & (t["state"][1:] == b)
            for name, holds in conditions.items():
                rows.append(dict(cell=c, block=block, rate=rate, source=name, k1=(event & holds).sum(), T1=(risk & holds).sum(),
                                 k0=(event & ~holds).sum(), T0=(risk & ~holds).sum()))
    return pd.DataFrame(rows)


def rate_ratios(cells, source, target, weights, dt, blocks):
    counts = rate_counts(cells, source, target, blocks)
    out = []
    for (rate, name), part in counts.groupby(["rate", "source"], sort=False):
        part = part.sort_values(["cell", "block"])
        k1, T1, k0, T0 = (part[c].to_numpy(dtype=float) for c in ("k1", "T1", "k0", "T0"))
        T = np.where(T1 + T0 > 0, T1 + T0, 1)
        member = part["cell"].to_numpy()
        W = weights[:, member]                            # a stratum is drawn whenever its cell is
        def mh(w):
            return (w @ (k1 * T0 / T)) / (w @ (k0 * T1 / T))
        def pooled(w):
            return ((w @ k1) / (w @ T1)) / ((w @ k0) / (w @ T0))
        with np.errstate(divide="ignore", invalid="ignore"):
            one = np.ones(len(k1))
            boot_mh, boot_pooled = mh(W), pooled(W)
            out.append(dict(regulator=GENES[source], regulated=GENES[target], source_state=name, rate=rate,
                            architecture=f"{source}{name} -> {target} {rate}",
                            ratio_within_cells=mh(one), low_within=interval(boot_mh)[0], high_within=interval(boot_mh)[1],
                            ratio_pooled=pooled(one), low_pooled=interval(boot_pooled)[0], high_pooled=interval(boot_pooled)[1],
                            rate_when_on_per_min=k1.sum() / T1.sum() / dt, rate_when_off_per_min=k0.sum() / T0.sum() / dt,
                            events_on=int(k1.sum()), events_off=int(k0.sum()),
                            cells_informative=int(len(set(member[(T1 > 0) & (T0 > 0)])))))
    return pd.DataFrame(out)

# ---------------------------------------------------------------------------
# 3. Between cells: per-cell statistics and the multi-tau correlation, trace-centred and global
# ---------------------------------------------------------------------------

def cell_statistics(cells, dt):
    rows = []
    for c, cell in enumerate(cells):
        row = dict(cell=c)
        for gene, trace in cell.items():
            name = GENES[gene]
            row[f"{name}_ON_fraction"] = (trace["reporters"] > 0).mean()
            row[f"{name}_mean_reporters"] = trace["reporters"].mean()
            row[f"{name}_bursts_per_hour"] = onsets(trace["reporters"]).sum() / (len(trace["reporters"]) * dt / 60)
            for s in (1, 2, 3):
                row[f"{name}_fraction_state{s}"] = (trace["state"] == s).mean()
        rows.append(row)
    return pd.DataFrame(rows)


def activity_over_time(cells, dt, bin_frames, weights):
    rows = []
    longest = max(len(cell["E"]["state"]) for cell in cells)
    for gene in ("E", "G"):
        for kind in ("ON", "in state 3", "bursts per hour"):
            values = np.full((len(cells), longest), np.nan)
            for c, cell in enumerate(cells):
                t = cell[gene]
                v = {"ON": t["reporters"] > 0, "in state 3": t["state"] == 3,
                     "bursts per hour": onsets(t["reporters"]) * 60 / dt}[kind]
                values[c, :len(v)] = v
            for start in range(0, longest, bin_frames):
                block = values[:, start:start + bin_frames]
                total, count = np.nansum(block, axis=1), np.sum(np.isfinite(block), axis=1)
                boot = (weights @ total) / (weights @ count)
                rows.append(dict(gene=GENES[gene], measure=kind, minute=(start + bin_frames / 2) * dt,
                                 value=total.sum() / count.sum(), low=np.percentile(boot, 2.5), high=np.percentile(boot, 97.5),
                                 cells=int((count > 0).sum())))
    return pd.DataFrame(rows)


def plot_activity(table):
    measures = list(dict.fromkeys(table["measure"]))
    figure, axes = plt.subplots(1, len(measures), figsize=(13, 3.4))
    for axis, measure in zip(axes, measures):
        for gene, color in (("CCAT1", RED), ("MYC", BLUE)):
            part = table[(table["gene"] == gene) & (table["measure"] == measure)]
            axis.fill_between(part["minute"], part["low"], part["high"], color=color, alpha=0.15, linewidth=0)
            axis.plot(part["minute"], part["value"], color=color, linewidth=1.3, label=gene)
        axis.set_title(measure, fontsize=9)
        axis.set_xlabel("minute of the recording", fontsize=8)
        axis.set_ylim(bottom=0)
        style(axis)
    axes[0].legend(fontsize=7, frameon=False)
    figure.suptitle("Activity over the recording (all cells, 95% interval)", x=0.01, ha="left", fontsize=10, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    return figure


def spearman(x, y):
    rx, ry = pd.Series(x).rank().to_numpy(), pd.Series(y).rank().to_numpy()
    return np.corrcoef(rx, ry)[0, 1]


def cell_correlations(stats, rng, reps):
    rows = []
    for e in [c for c in stats.columns if c.startswith("CCAT1_")]:
        for g in [c for c in stats.columns if c.startswith("MYC_")]:
            x, y = stats[e].to_numpy(), stats[g].to_numpy()
            if x.std() == 0 or y.std() == 0:
                continue
            value = spearman(x, y)
            boot = [spearman(x[i], y[i]) for i in (rng.integers(0, len(x), len(x)) for _ in range(reps))]
            null = np.array([spearman(x, rng.permutation(y)) for _ in range(reps)])
            rows.append(dict(CCAT1=e[6:], MYC=g[4:], spearman=value, low=np.nanpercentile(boot, 2.5),
                             high=np.nanpercentile(boot, 97.5), p_permutation=(np.abs(null) >= abs(value)).mean()))
    return pd.DataFrame(rows)


def idl_round(x):
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def multi_tau(left, right, M=4, sign=1):
    """Port of multi_tau_corr in fcs_multitau_xcor_final.R (bg_correct = FALSE): lags (frames), Go, Md, Me, n."""
    selected = len(left) - 1
    N = int(idl_round(np.log2((selected + 1) / M)))
    L, R = left[:selected].astype(float), right[:selected].astype(float)
    rows = []
    for lag in range(1, 2 * M + 1):
        j = selected - lag
        rows.append((sign * lag, (L[:j] * R[lag:lag + j]).sum(), L[:j].sum(), R[lag:lag + j].sum(), j) if j >= 1 else (sign * lag, 0, 0, 0, 0))
    points = selected
    for p in range(1, N - 1):
        points -= points % 2
        L, R = L[0:points:2] + L[1:points:2], R[0:points:2] + R[1:points:2]
        points = len(L)
        for i in range(1, M + 1):
            lag = M + i
            j = points - lag
            rows.append((sign * lag * 2 ** p, (L[:j] * R[lag:lag + j]).sum(), L[:j].sum(), R[lag:lag + j].sum(), j) if j >= 1
                        else (sign * lag * 2 ** p, 0, 0, 0, 0))
    return np.array(rows, dtype=float)


def correlation_sums(cells, observable, pair, epsilon):
    """Per cell: the multi-tau sums of the forward and backward correlation, on the common lag grid."""
    a, b = {"EG": ("E", "G"), "EE": ("E", "E"), "GG": ("G", "G")}[pair]
    out = []
    for cell in cells:
        x, y = (cell[g]["reporters"].astype(float) for g in (a, b))
        if observable == "binary":
            x, y = (x > 0).astype(float), (y > 0).astype(float)
        x, y = x + epsilon, y + epsilon
        out.append(np.vstack([multi_tau(x, y), multi_tau(y, x, sign=-1)]))
    common = sorted(set.intersection(*(set(o[:, 0]) for o in out)))
    return np.array(common), np.stack([o[[list(o[:, 0]).index(l) for l in common]] for o in out])   # [cell, lag, column]


def coulon(sums, w):
    """Coulon average of fcs_multitau_xcor_final.R (each trace around its own means), with cell weights w.
    As in the R code, each curve's weight and the global denominator use its first channel (forward lag 1)."""
    Go, Md, Me, n = (sums[:, :, k] for k in (1, 2, 3, 4))
    f = int(np.flatnonzero(sums[0, :, 0] == 1)[0])
    with np.errstate(divide="ignore", invalid="ignore"):
        Gn = n * Go / (Md * Me) - 1
        first = Md[:, f:f + 1] * Me[:, f:f + 1] / n[:, f:f + 1] ** 2
        compiled = np.where(np.isfinite(Gn), Gn * first, 0)
        counts = np.isfinite(Gn)
        denominator = (w @ Md[:, f]) / (w @ n[:, f]) * (w @ Me[:, f]) / (w @ n[:, f])
        return (w @ compiled) / (w @ counts) / denominator


def global_average(sums, w):
    """The same correlation around the global means: pooled <xy> / (pooled <x> pooled <y>) - 1, per lag."""
    Go, Md, Me, n = (sums[:, :, k] for k in (1, 2, 3, 4))
    return (w @ Go) * (w @ n) / ((w @ Md) * (w @ Me)) - 1

# ---------------------------------------------------------------------------
# Plots and report
# ---------------------------------------------------------------------------

def style(axis):
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    axis.tick_params(labelsize=7, colors=MUTED, labelcolor=INK)


def plot_aligned(results, dt):
    pairs = list(results)
    kinds = list(results[pairs[0]])
    figure, axes = plt.subplots(len(pairs), len(kinds), figsize=(13, 2.6 * len(pairs)), squeeze=False, sharex=True)
    for i, pair in enumerate(pairs):
        for j, kind in enumerate(kinds):
            r, axis = results[pair][kind], axes[i, j]
            minutes = r["lags"] * dt
            axis.fill_between(minutes, 100 * r["low"], 100 * r["high"], color=BLUE, alpha=0.2, linewidth=0)
            axis.plot(minutes, 100 * r["change"], color=BLUE, linewidth=1.2)
            axis.axhline(0, color=MUTED, linewidth=0.6)
            axis.axvline(0, color=MUTED, linewidth=0.6)
            axis.axhline(100 * r["level"], color=RED, linewidth=0.9, linestyle="--")
            axis.set_title(f"{pair}: {kind}  ({r['events']:,} reference bursts)", fontsize=8.5)
            style(axis)
            if j == 0:
                axis.set_ylabel("change vs expected (%)", fontsize=8)
    for axis in axes[-1]:
        axis.set_xlabel("minutes after the reference burst onset", fontsize=8)
    figure.suptitle("Burst-aligned averages. Blue: change at each lag within cells (95% interval, bootstrap over cells). "
                    "Red dashed: between-cell level", x=0.01, ha="left", fontsize=10, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    return figure


def plot_windows(results, dt):
    pairs = list(results)
    kinds = list(results[pairs[0]])
    colors = {"before": "#eb6834", "after": BLUE, "around": INK}
    figure, axes = plt.subplots(len(pairs), len(kinds), figsize=(13, 2.6 * len(pairs)), squeeze=False, sharex=True)
    for i, pair in enumerate(pairs):
        for j, kind in enumerate(kinds):
            table, axis = results[pair][kind]["windows"], axes[i, j]
            for side, part in table.groupby("side", sort=False):
                minutes = part["window_frames"] * dt
                axis.fill_between(minutes, 100 * part["low"], 100 * part["high"], color=colors[side], alpha=0.15, linewidth=0)
                axis.plot(minutes, 100 * part["change"], color=colors[side], linewidth=1.2, label=side)
            axis.axhline(0, color=MUTED, linewidth=0.6)
            axis.set_title(f"{pair}: {kind}", fontsize=8.5)
            style(axis)
            if j == 0:
                axis.set_ylabel("average change\nwithin the window (%)", fontsize=8)
    axes[0, 0].legend(fontsize=7, frameon=False)
    for axis in axes[-1]:
        axis.set_xlabel("window width (minutes): before = [-w, 0), after = (0, w], around = [-w, w]", fontsize=8)
    figure.suptitle("Average change within every window, no cutoff (95% interval)", x=0.01, ha="left", fontsize=10, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    return figure


def plot_ratios(table):
    figure, axes = plt.subplots(1, 2, figsize=(13, 7), sharex=True)
    for axis, ((regulator, regulated), part) in zip(axes, table.groupby(["regulator", "regulated"], sort=False)):
        labels = [f"{regulator} {s}  ->  {regulated} {rate}" for rate in [r[0] for r in RATES] for s in SOURCES]
        part = part.set_index(["rate", "source_state"]).loc[[(r[0], s) for r in RATES for s in SOURCES]]
        y = np.arange(len(part))[::-1]
        for offset, kind, color in ((0.15, "within", BLUE), (-0.15, "pooled", MUTED)):
            value, low, high = part[f"ratio_{kind}_cells" if kind == "within" else "ratio_pooled"], part[f"low_{kind}"], part[f"high_{kind}"]
            axis.errorbar(value, y + offset, xerr=[value - low, high - value], fmt="o", color=color, markersize=3.5, elinewidth=1,
                          label="within cells" if kind == "within" else "pooled (includes between-cell)")
        axis.axvline(1, color=INK, linewidth=0.7)
        axis.set_xscale("log")
        axis.set_yticks(y, labels, fontsize=7)
        for k in range(1, len(RATES)):
            axis.axhline(y[0] - 4 * k + 0.5, color="#e1e0d9", linewidth=0.6)
        axis.set_xlabel("rate while the source condition holds / rate while it does not (95% interval)", fontsize=8)
        axis.set_title(f"{regulator} regulating {regulated}", fontsize=10)
        style(axis)
    axes[0].legend(fontsize=7, frameon=False, loc="lower right")
    figure.suptitle("Model-free fold of every architecture. 1 = no effect; source R = the regulator has reporters",
                    x=0.01, ha="left", fontsize=10, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.96))
    return figure


def plot_cells(stats, correlations):
    figure, axes = plt.subplots(1, 3, figsize=(14, 4.2), gridspec_kw=dict(width_ratios=[1, 1, 1.4]))
    for axis, (x, y) in zip(axes[:2], [("CCAT1_ON_fraction", "MYC_ON_fraction"), ("CCAT1_fraction_state2", "MYC_ON_fraction")]):
        axis.scatter(stats[x], stats[y], s=14, color=BLUE, alpha=0.6, linewidths=0)
        row = correlations[(correlations["CCAT1"] == x[6:]) & (correlations["MYC"] == y[4:])]
        if len(row):
            r = row.iloc[0]
            axis.set_title(f"Spearman {r['spearman']:+.2f} [{r['low']:+.2f}, {r['high']:+.2f}], p = {r['p_permutation']:.3f}", fontsize=8.5)
        axis.set_xlabel(x.replace("_", " "), fontsize=8)
        axis.set_ylabel(y.replace("_", " "), fontsize=8)
        style(axis)
    matrix = correlations.pivot(index="CCAT1", columns="MYC", values="spearman")
    image = axes[2].imshow(matrix.to_numpy(), cmap="RdBu_r", vmin=-0.6, vmax=0.6, aspect="auto")
    axes[2].set_xticks(range(matrix.shape[1]), [c.replace("_", " ") for c in matrix.columns], rotation=40, ha="right", fontsize=7)
    axes[2].set_yticks(range(matrix.shape[0]), [c.replace("_", " ") for c in matrix.index], fontsize=7)
    for (i, j), value in np.ndenumerate(matrix.to_numpy()):
        axes[2].text(j, i, f"{value:+.2f}", ha="center", va="center", fontsize=6)
    axes[2].set_xlabel("MYC, per cell", fontsize=8)
    axes[2].set_ylabel("CCAT1, per cell", fontsize=8)
    figure.colorbar(image, ax=axes[2], label="Spearman across cells")
    figure.suptitle("Between cells: one point per cell", x=0.01, ha="left", fontsize=10, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.94))
    return figure


def plot_correlations(table, max_lag_min):
    observables = list(dict.fromkeys(table["observable"]))
    figure, axes = plt.subplots(len(observables), 3, figsize=(14, 3.4 * len(observables)), squeeze=False)
    names = {"EG": "CCAT1 x MYC cross-correlation", "EE": "CCAT1 autocorrelation", "GG": "MYC autocorrelation"}
    for i, observable in enumerate(observables):
        for j, pair in enumerate(("EG", "EE", "GG")):
            part = table[(table["observable"] == observable) & (table["pair"] == pair) & (table["lag_min"].abs() <= max_lag_min)]
            axis = axes[i, j]
            for column, color, label in (("trace_centered", BLUE, "trace-centred (Coulon, as in your R script)"),
                                         ("global", RED, "global means (within + between cells)"),
                                         ("between", MUTED, "difference = between cells")):
                axis.fill_between(part["lag_min"], part[column] - part[f"{column}_se"], part[column] + part[f"{column}_se"],
                                  color=color, alpha=0.15, linewidth=0)
                axis.plot(part["lag_min"], part[column], color=color, linewidth=1.2, label=label,
                          linestyle="--" if column == "between" else "-")
            axis.axhline(0, color=MUTED, linewidth=0.6)
            axis.set_title(f"{names[pair]}, {observable}", fontsize=9)
            axis.set_xlabel("lag (min)" + ("; positive = MYC after CCAT1" if pair == "EG" else ""), fontsize=8)
            style(axis)
    axes[0, 0].legend(fontsize=7, frameon=False)
    figure.suptitle("Multi-tau correlation around each trace's own means and around the global means (± bootstrap SE over cells)",
                    x=0.01, ha="left", fontsize=10, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.95))
    return figure


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ccat1")
    parser.add_argument("myc")
    parser.add_argument("output")
    parser.add_argument("--interval", type=float, default=5 / 3, help="minutes per frame")
    parser.add_argument("--bootstrap", type=int, default=1000)
    parser.add_argument("--max-lag-min", type=float, default=120.0)
    parser.add_argument("--epsilon", type=float, default=1e-7, help="added to traces before correlating, as in step 0")
    parser.add_argument("--trend-window", type=int, default=25,
                        help="frames over which the population's time course is smoothed (0 = assume no drift)")
    parser.add_argument("--time-blocks", type=int, default=4, help="parts of each trace that rate ratios are compared within")
    parser.add_argument("--seed", type=int, default=20261117)
    arguments = parser.parse_args()
    dt, rng = arguments.interval, np.random.default_rng(arguments.seed)
    out = Path(arguments.output)
    (out / "plots").mkdir(parents=True, exist_ok=True)

    E, G = read_traces(Path(arguments.ccat1)), read_traces(Path(arguments.myc))
    if sorted(E) != sorted(G):
        raise SystemExit("the two files do not have the same trace numbers")
    cells = [dict(E=E[i], G=G[i]) for i in sorted(E) if len(E[i]["state"]) == len(G[i]["state"]) and len(E[i]["state"]) > 1]
    weights = resample(len(cells), arguments.bootstrap, rng)
    K = int(round(arguments.max_lag_min / dt))
    print(f"{len(cells)} cells, {sum(len(c['E']['state']) for c in cells):,} frames")

    # 1. burst-aligned averages
    activity = activity_over_time(cells, dt, 12, weights)
    activity.to_csv(out / "activity_over_time.csv", index=False)
    aligned = {f"{GENES[a]} burst -> {GENES[b]}": burst_aligned(cells, a, b, K, weights, arguments.trend_window) for a, b in (("E", "G"), ("G", "E"), ("E", "E"), ("G", "G"))}
    rows, windows = [], []
    for pair, kinds in aligned.items():
        for kind, r in kinds.items():
            rows.append(pd.DataFrame(dict(reference_to_response=pair, response=kind, lag_frames=r["lags"], lag_min=r["lags"] * dt,
                                          change=r["change"], low=r["low"], high=r["high"], level_between_cells=r["level"],
                                          level_low=r["level_low"], level_high=r["level_high"], reference_bursts=r["events"])))
            windows.append(r["windows"].assign(reference_to_response=pair, response=kind, window_min=r["windows"]["window_frames"] * dt))
    pd.concat(rows).to_csv(out / "burst_aligned.csv", index=False)
    pd.concat(windows).to_csv(out / "burst_aligned_windows.csv", index=False)

    # 2. rate ratios
    ratios = pd.concat([rate_ratios(cells, a, b, weights, dt, arguments.time_blocks) for a, b in (("E", "G"), ("G", "E"))])
    ratios.to_csv(out / "rate_ratios.csv", index=False)

    # 3. between cells
    stats = cell_statistics(cells, dt)
    stats.to_csv(out / "cell_statistics.csv", index=False)
    cell_corr = cell_correlations(stats, rng, min(arguments.bootstrap, 1000))
    cell_corr.to_csv(out / "cell_correlations.csv", index=False)
    parts = []
    for observable in ("reporter", "binary"):
        for pair in ("EG", "EE", "GG"):
            lags, sums = correlation_sums(cells, observable, pair, arguments.epsilon)
            one = np.ones(len(cells))
            centred, overall = coulon(sums, one), global_average(sums, one)
            boot_c = np.array([coulon(sums, w) for w in weights])
            boot_g = np.array([global_average(sums, w) for w in weights])
            parts.append(pd.DataFrame(dict(observable=observable, pair=pair, lag_frames=lags, lag_min=lags * dt,
                                           trace_centered=centred, trace_centered_se=boot_c.std(axis=0),
                                           **{"global": overall, "global_se": boot_g.std(axis=0)},
                                           between=overall - centred, between_se=(boot_g - boot_c).std(axis=0))))
    correlations = pd.concat(parts)
    correlations.to_csv(out / "correlations.csv", index=False)

    figures = [("Activity over the recording", plot_activity(activity), "activity_over_time",
                "Fraction of cells ON, in promoter state 3, and burst onsets per hour, in 20-minute bins. A drift means the "
                "recordings are not stationary: stationary models (the sweep, StochasticGene fits) cannot reproduce it, and it "
                "shapes the correlations. The analyses below compare against this time course, not a flat average."),
               ("Burst-aligned averages", plot_aligned(aligned, dt), "burst_aligned",
                "Each row: bursts of the first gene are the reference; the response is the second gene. Blue: the response at "
                "each lag relative to its expected value (the same cell's level times the population's time course, so only "
                "timing within cells counts, and drifts over the recording do not); for example +20% at +10 min "
                "= MYC is 20% more likely to be firing 10 min after a CCAT1 burst began than at a random time in that cell. "
                "Red dashed: how much more active the response gene is, on average, in cells weighted by their number of reference "
                "bursts, versus all cells (between-cell effect; not timing). The same-gene rows show each gene's refractory period."),
               ("Average change within every window", plot_windows(aligned, dt), "burst_aligned_windows",
                "The average change over every window width, so no cutoff has to be chosen: before the reference onset, after it, "
                "and around it. Where a curve crosses zero shows the time scale of the effect."),
               ("Model-free fold of every architecture", plot_ratios(ratios), "rate_ratios",
                "Rate of each transition of the regulated gene while the regulator is in a given state (or has reporters, R), divided "
                "by the rate while it is not: the fold the sweep's coupling would need. Blue: compared within each cell "
                "and within the same part of the recording (Mantel-Haenszel), which removes cells that are simply more active and "
                "drifts over time. Grey: all cells and times pooled. "
                "States come from Viterbi decoding with the uncoupled model; decoding errors pull ratios towards 1."),
               ("Between cells", plot_cells(stats, cell_corr), "between_cells",
                "Per-cell activity of CCAT1 against MYC. A coupling acting from a long-lived state (for example strong inhibition "
                "from CCAT1 state 2) predicts a clear correlation here; timing effects within cells do not."),
               ("Correlation within versus between cells", plot_correlations(correlations, arguments.max_lag_min), "correlations",
                "Blue: your multi-tau correlation (Coulon average: each trace around its own means, divided by the global means). "
                "Red: the same lags around the global means. Grey dashed: their difference, the part that comes from differences "
                "between cells; it is nearly constant across lags.")]
    cards = []
    for title, figure, stem, caption in figures:
        figure.savefig(out / "plots" / f"{stem}.pdf")
        figure.savefig(out / "plots" / f"{stem}.png", dpi=300)
        plt.close(figure)
        cards.append(f'<h2>{html.escape(title)}</h2><p>{html.escape(caption)} <a href="plots/{stem}.pdf">PDF</a> · '
                     f'<a href="plots/{stem}.png">PNG</a></p><img src="plots/{stem}.png" alt="{html.escape(title)}">')
    best = ratios.assign(distance=np.abs(np.log(ratios["ratio_within_cells"]))).sort_values("distance", ascending=False).head(12)
    table = best.drop(columns="distance")[["architecture", "ratio_within_cells", "low_within", "high_within", "ratio_pooled",
                                           "events_on", "events_off", "cells_informative"]].to_html(index=False, float_format=lambda x: f"{x:.3g}")
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>CCAT1 and MYC in decoded traces</title>
<style>body{{font:15px/1.5 system-ui,sans-serif;max-width:1200px;margin:24px auto;padding:0 16px;color:#0b0b0b;background:#fcfcfb}}
img{{max-width:100%;border:1px solid #e1e0d9}}table{{border-collapse:collapse;font-size:13px}}td,th{{border-bottom:1px solid #e1e0d9;padding:3px 8px;text-align:right}}
p{{color:#52514e}}</style></head><body><h1>CCAT1 and MYC in decoded traces</h1>
<p>{len(cells)} cells, {dt:.3g} min per frame, {arguments.bootstrap} bootstrap replicates over cells. Inputs: {html.escape(Path(arguments.ccat1).name)},
{html.escape(Path(arguments.myc).name)}. Tables: burst_aligned.csv, burst_aligned_windows.csv, rate_ratios.csv, cell_statistics.csv,
cell_correlations.csv, correlations.csv.</p>
<h2>Architectures with the largest within-cell effect</h2>{table}
{''.join(cards)}</body></html>"""
    (out / "relationships_report.html").write_text(page)
    print(f"Wrote {out / 'relationships_report.html'}")


if __name__ == "__main__":
    main()
