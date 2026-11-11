#!/usr/bin/env python3
"""Step 4. Plot the simulation results against the experimental data and rank the simulations.

    python3 4_plot_results.py config.toml

Reads <run_dir>/simulations.csv and <run_dir>/results/, and writes into <run_dir>:
    plots.pdf           pages 1-6: one heatmap per correlation curve (cross-correlation and both
                                   autocorrelations, for ON/OFF state and for reporter count)
                        page 7:    heatmap of the ON fractions relative to experiment
                        page 8:    the ten best simulations as line plots against experiment
                        Heatmap rows are simulations, grouped by model group, then model, then parameters.
    ranked_scores.csv   every simulation with its parameters, metrics and score, best first
"""
import sys, tomllib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap

PLOT_MAX_LAG_MIN = 120      # lag range shown in every plot
TOP_N = 10                  # number of best simulations drawn as line plots


def score(metrics):
    """The ranking score of every simulation; LOWER IS BETTER. Edit this function to design your own score.

    `metrics` has one row per simulation with these columns:
        nrmse_<observable>_<pair>   0 = identical to experiment
        ccc_<observable>_<pair>     1 = identical to experiment
        on_E, on_G, reporters_E, reporters_G
    observable: binary (ON/OFF state) or reporter (reporter count)
    pair:       EG (CCAT1 x MYC cross-correlation), EE (CCAT1 autocorrelation), GG (MYC autocorrelation)
    Example of a combined score:  metrics["nrmse_reporter_EG"] + metrics["nrmse_binary_EG"]
    """
    return metrics["nrmse_reporter_EG"]


CURVES = [(observable, pair, f"{pair_name}, {observable_name}")
          for pair, pair_name in [("EG", "CCAT1 x MYC cross-correlation"), ("EE", "CCAT1 autocorrelation"),
                                  ("GG", "MYC autocorrelation")]
          for observable, observable_name in [("binary", "ON/OFF state"), ("reporter", "reporter count")]]
PARAMETERS = ["fold", "h12", "h21", "h23", "h32"]
INK, MUTED, BLUE = "#0b0b0b", "#898781", "#2a78d6"
DIVERGING = LinearSegmentedColormap.from_list("blue_gray_red", ["#104281", "#2a78d6", "#f0efec", "#e34948", "#8f1f1f"])
DIVERGING.set_bad("#c3c2b7")      # no data


def load(config_path):
    config = tomllib.loads(config_path.read_text())
    folder = config_path.parent
    run_dir = folder / config["output"]["run_dir"]
    read_all = lambda pattern: pd.concat([pd.read_csv(f) for f in sorted((run_dir / "results").glob(pattern))])
    simulations = pd.read_csv(run_dir / "simulations.csv")
    metrics, curves = read_all("metrics_*.csv"), read_all("curves_*.csv")
    experiment = pd.read_csv(folder / config["inputs"]["experimental_curves"])
    experiment_on = pd.read_csv(folder / config["inputs"]["experimental_metrics"]).set_index("metric")["mean"]
    print(f"{len(metrics)} of {len(simulations)} simulations have results")
    # One row per simulation: parameters, metrics, score and rank (1 = best).
    table = simulations.merge(metrics, on=["sim_id", "model_id"]).sort_values("sim_id")
    table["score"] = score(table)
    table["rank"] = table["score"].rank(method="first").astype("Int64")
    for frame in (curves, experiment):
        frame.drop(frame[frame["lag_min"].abs() > PLOT_MAX_LAG_MIN].index, inplace=True)
    return run_dir, table, curves.reset_index(drop=True), experiment, experiment_on


def bin_edges(centers):
    """Edges of heatmap cells centred on the (unevenly spaced) multi-tau lags."""
    middle = (centers[:-1] + centers[1:]) / 2
    return np.concatenate([[2 * centers[0] - middle[0]], middle, [2 * centers[-1] - middle[-1]]])


def heatmap_page(pdf, title, table, values, x_edges, limit, colorbar_label, top_row=None, x_label="", x_ticks=None):
    """One page: a panel per model group, rows = that group's simulations in sim_id order.

    values: DataFrame indexed by sim_id, one column per heatmap column. top_row: experimental values."""
    groups = list(dict.fromkeys(table["group"]))
    heights = ([0.35] if top_row is not None else []) + [1] * len(groups)
    figure, axes = plt.subplots(len(heights), 1, figsize=(8.5, 11), gridspec_kw=dict(height_ratios=heights, hspace=0.12),
                                squeeze=False)
    axes = list(axes[:, 0])
    style = dict(cmap=DIVERGING, vmin=-limit, vmax=limit, rasterized=True)
    if top_row is not None:
        mesh = axes[0].pcolormesh(x_edges, [0, 1], top_row[None, :], **style)
        axes[0].set_ylabel("experiment", rotation=0, ha="right", va="center", fontweight="bold")
    for axis, group in zip(axes[len(heights) - len(groups):], groups):
        rows = table[table["group"] == group]
        block = values.reindex(rows["sim_id"]).to_numpy()
        mesh = axis.pcolormesh(x_edges, np.arange(len(rows) + 1), block, **style)
        axis.invert_yaxis()                                     # first simulation of the group on top
        models = rows["model_id"].to_numpy()
        boundaries = np.flatnonzero(models[1:] != models[:-1]) + 1
        if len(boundaries) <= 12:                               # thin lines between the group's models
            axis.hlines(boundaries, x_edges[0], x_edges[-1], color="white", linewidth=0.6)
        axis.set_ylabel(f"{group}\n{len(rows)} simulation{'s' * (len(rows) > 1)}", rotation=0, ha="right", va="center")
    for axis in axes:
        axis.set_yticks([])
        axis.tick_params(labelbottom=axis is axes[-1], length=3, colors=MUTED, labelcolor=INK)
        if x_ticks is not None:
            axis.set_xticks(*x_ticks)
        for spine in axis.spines.values():
            spine.set_visible(False)
    axes[-1].set_xlabel(x_label)
    figure.suptitle(title, x=0.06, ha="left", fontsize=13, fontweight="bold")
    figure.subplots_adjust(left=0.27, right=0.86, top=0.94, bottom=0.06)
    figure.colorbar(mesh, cax=figure.add_axes([0.89, 0.3, 0.018, 0.4]), label=colorbar_label, extend="both")
    pdf.savefig(figure, dpi=200)
    plt.close(figure)


def curve_pages(pdf, table, curves, experiment):
    for observable, pair, name in CURVES:
        select = lambda frame: frame[(frame["observable"] == observable) & (frame["pair"] == pair)]
        values = select(curves).pivot(index="sim_id", columns="lag_min", values="mean")
        lags = values.columns.to_numpy()
        target = select(experiment).set_index("lag_frames")["mean"]
        frames = select(curves).drop_duplicates("lag_min").set_index("lag_min")["lag_frames"].reindex(lags)
        top_row = target.reindex(frames.to_numpy()).to_numpy()            # experiment on the simulated lags
        limit = 2 * np.nanmax(np.abs(top_row))        # colour scale set by the experiment; stronger simulations saturate
        heatmap_page(pdf, name, table, values, bin_edges(lags), limit, "correlation (blue negative, red positive)",
                     top_row=top_row, x_label="lag (min); positive = MYC after CCAT1" if pair == "EG" else "lag (min)")


def on_fraction_page(pdf, table, experiment_on):
    ratio = pd.DataFrame({"E": np.log2(table["on_E"] / experiment_on["on_CCAT1"]),
                          "G": np.log2(table["on_G"] / experiment_on["on_MYC"])}).set_index(table["sim_id"])
    limit = max(1.0, np.nanpercentile(np.abs(ratio.to_numpy()), 99))
    title = (f'ON fraction relative to experiment (experiment: CCAT1 {experiment_on["on_CCAT1"]:.3f}, '
             f'MYC {experiment_on["on_MYC"]:.3f})')
    heatmap_page(pdf, title, table, ratio, np.array([0.0, 1.0, 2.0]), limit,
                 "log2(simulated / experimental); 0 = equal, +1 = twice as often ON",
                 x_ticks=([0.5, 1.5], ["CCAT1 (E)", "MYC (G)"]))


def top_page(pdf, table, curves, experiment, experiment_on):
    """The TOP_N best simulations: every curve as a line plot against experiment, one row per simulation."""
    best = table.sort_values("rank").head(TOP_N)
    figure, axes = plt.subplots(len(best), len(CURVES), figsize=(22, 1.9 * len(best) + 1.6), sharex=True, sharey="col",
                                squeeze=False)
    for row, (_, sim) in enumerate(best.iterrows()):
        parameters = ", ".join(f"{name} {sim[name]:g}" for name in PARAMETERS if pd.notna(sim[name]))
        label = (f'#{sim["rank"]}   simulation {sim["sim_id"]}   score {sim["score"]:.3f}\n{sim["model_id"]}\n{parameters}\n'
                 f'ON: CCAT1 {sim["on_E"]:.3f} (exp. {experiment_on["on_CCAT1"]:.3f}), '
                 f'MYC {sim["on_G"]:.3f} (exp. {experiment_on["on_MYC"]:.3f})')
        axes[row, 0].set_ylabel(label, rotation=0, ha="right", va="center", fontsize=8.5)
        for column, (observable, pair, name) in enumerate(CURVES):
            axis = axes[row, column]
            select = lambda frame: frame[(frame["observable"] == observable) & (frame["pair"] == pair)].sort_values("lag_min")
            for data, color, text in [(select(experiment), INK, "experiment ± SE"),
                                      (select(curves[curves["sim_id"] == sim["sim_id"]]), BLUE, "simulation ± SE over replicates")]:
                axis.fill_between(data["lag_min"], data["mean"] - data["se"], data["mean"] + data["se"],
                                  color=color, alpha=0.18, linewidth=0)
                axis.plot(data["lag_min"], data["mean"], color=color, linewidth=1.6, label=text)
            axis.axhline(0, color=MUTED, linewidth=0.6)
            axis.tick_params(labelsize=8, colors=MUTED, labelcolor=INK)
            for side in ("top", "right"):
                axis.spines[side].set_visible(False)
            if row == 0:
                axis.set_title(name, fontsize=9.5)
            if row == len(best) - 1:
                axis.set_xlabel("lag (min)", fontsize=9)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper right", ncol=2, frameon=False)
    figure.suptitle(f"The {len(best)} best simulations (lowest score)", x=0.02, ha="left", fontsize=13, fontweight="bold")
    figure.subplots_adjust(left=0.19, right=0.99, top=1 - 1.0 / figure.get_figheight(), bottom=0.6 / figure.get_figheight(),
                           wspace=0.22, hspace=0.25)
    pdf.savefig(figure)
    plt.close(figure)


def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    run_dir, table, curves, experiment, experiment_on = load(Path(sys.argv[1]).resolve())
    table.sort_values("rank").to_csv(run_dir / "ranked_scores.csv", index=False)
    with PdfPages(run_dir / "plots.pdf") as pdf:
        curve_pages(pdf, table, curves, experiment)
        on_fraction_page(pdf, table, experiment_on)
        top_page(pdf, table, curves, experiment, experiment_on)
    print(f"Wrote {run_dir / 'plots.pdf'} and {run_dir / 'ranked_scores.csv'}")
    print(table.sort_values("rank").head(TOP_N)[["rank", "sim_id", "model_id", "score"]].to_string(index=False))


if __name__ == "__main__":
    main()
