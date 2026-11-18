#!/usr/bin/env python3
"""Report: ON/OFF dwell-time histograms against the fitted models, and which models fit best (step 2 of 2).

    python3 dwelltime_report.py RESULTFOLDER [RESULTFOLDER ...] [--out dwelltime_report.pdf] [--root SG_ROOT]

Reads <RESULTFOLDER>/dwelltime_report/fits.csv and distributions.csv written by dwelltime_report.jl (RESULTFOLDER under
<root>/results, or an absolute path). Fits are compared only with fits of the same histograms (same source, sample
and condition, i.e. the same data files); folders with fits of the same histograms (for example the expanded
topologies and the splicing models of the enhancer) are compared together.

For every set of histograms, two pages:
  1. the models ranked by AIC: delta AIC (AIC minus the best AIC; below 2 = as good as the best, above 10 = clearly
     worse), Akaike weight, Kolmogorov-Smirnov distance of ON and OFF (largest gap between the data and model CDFs,
     independent of the number of parameters), and the largest Rhat (above 1.1 = the chains have not converged)
  2. ON and OFF: data (points) and models (lines; the three best by AIC in colour), as probability per bin and as
     survival (fraction of dwell times in this bin or later), both on a log scale
Also writes model_ranking.csv (next to the PDF): every fit with its delta AIC and weight within its histograms.
"""
import argparse, os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

BEST = ["#d1495b", "#2a78d6", "#3c9d5d"]     # the three best models
OTHER = "#c8c8c8"


def read(folders, root):
    fits, dists = [], []
    for folder in folders:
        path = Path(folder) if os.path.isabs(folder) else Path(root) / "results" / folder
        report = path / "dwelltime_report"
        if not (report / "fits.csv").exists():
            raise SystemExit(f"{report / 'fits.csv'} not found: run dwelltime_report.jl {folder} first")
        f = pd.read_csv(report / "fits.csv")
        d = pd.read_csv(report / "distributions.csv")
        f["folder"] = d["folder"] = path.name
        fits.append(f)
        dists.append(d)
    return pd.concat(fits, ignore_index=True), pd.concat(dists, ignore_index=True)


def rank(fits):
    """Delta AIC and Akaike weight within each set of histograms; label = model, plus the folder where the same model
    name was fitted in two folders."""
    fits = fits.copy()
    fits["delta_AIC"] = fits["AIC"] - fits.groupby("data")["AIC"].transform("min")
    w = np.exp(-0.5 * fits["delta_AIC"])
    fits["weight"] = w / w.groupby(fits["data"]).transform("sum")
    repeated = fits.duplicated(["data", "model"], keep=False)
    fits["label"] = np.where(repeated, fits["model"] + " (" + fits["folder"] + ")", fits["model"])
    return fits.sort_values(["data", "AIC"]).reset_index(drop=True)


def title_of(group):
    """'model / 3Prime / enhancer' from the data files of a group."""
    first = str(group["data"].iloc[0]).split(";")[0]
    source = Path(first).parent.name
    name = Path(first).stem.split("_")                  # 3Prime_enhancer_ON
    return f"{source} histograms, {name[0]}, {name[1] if len(name) > 1 else ''}"


def ranking_page(group, title):
    n = len(group)
    fig, ax = plt.subplots(figsize=(13, max(3.0, 0.32 * n + 1.6)))
    y = np.arange(n)[::-1]
    delta = group["delta_AIC"].to_numpy()
    colours = [BEST[i] if i < 3 else "#7a7a7a" for i in range(n)]
    ax.barh(y, np.maximum(delta, 0) + 1, color=colours, height=0.7)
    ax.set_xscale("log")
    ax.set_yticks(y)
    ax.set_yticklabels(group["label"], fontsize=8)
    ax.axvline(3, color="black", lw=0.6, ls=":")      # delta AIC = 2
    ax.axvline(11, color="black", lw=0.6, ls="--")    # delta AIC = 10
    ax.set_xlabel("delta AIC + 1 (log scale; dotted = delta AIC 2, dashed = 10)")
    for yi, (_, row) in zip(y, group.iterrows()):
        flag = "  not converged" if row["max_rhat"] > 1.1 else ""
        ax.text(1.03, yi, f'dAIC {row["delta_AIC"]:8.1f}   w {row["weight"]:.2f}   KS ON {row["KS_ON"]:.3f}   '
                f'KS OFF {row["KS_OFF"]:.3f}   Rhat {row["max_rhat"]:.2f}   k {row["parameters"]}{flag}',
                va="center", fontsize=7, family="monospace", transform=ax.get_yaxis_transform())
    events = f'{group["events_ON"].iloc[0]:.0f} ON and {group["events_OFF"].iloc[0]:.0f} OFF events'
    ax.set_title(f"{title}: models ranked by AIC ({events})", fontsize=10, loc="left")
    fig.subplots_adjust(left=0.17, right=0.47)
    return fig


def shown(p, kind):
    """Probability per bin, or survival (fraction in this bin or later), with values the log scale cannot show
    (empty bins; values below 1e-5 of the largest, where the numerical solution of the model is noise) left out."""
    values = p if kind == "probability" else np.cumsum(p[::-1])[::-1]
    return np.where(values > 1e-5 * np.nanmax(values), values, np.nan)


def distribution_page(group, dists, title):
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    order = list(group.index)
    for col, dttype in enumerate(["ON", "OFF"]):
        for row, kind in enumerate(["probability", "survival"]):
            ax = axes[row, col]
            data_drawn = False
            for rank_i, i in reversed(list(enumerate(order))):      # best drawn last, on top
                fit = group.loc[i]
                d = dists[(dists["folder"] == fit["folder"]) & (dists["key"] == fit["key"]) & (dists["dttype"] == dttype)]
                if d.empty:
                    continue
                d = d.sort_values("time")
                if not data_drawn:
                    ax.plot(d["time"], shown(d["data"].to_numpy(), kind), "o", ms=3, color="black", label="data", zorder=10)
                    data_drawn = True
                best = rank_i < 3
                ax.plot(d["time"], shown(d["model"].to_numpy(), kind), color=BEST[rank_i] if best else OTHER,
                        lw=1.6 if best else 0.7, zorder=5 if best else 1,
                        label=f'{fit["label"]} (dAIC {fit["delta_AIC"]:.1f})' if best else None)
            ax.set_yscale("log")
            ax.set_xlabel("dwell time (as in the histogram files)")
            ax.set_ylabel("probability per bin" if kind == "probability" else "fraction in this bin or later")
            ax.set_title(f"{dttype} times, {kind}", fontsize=9)
    axes[0, 0].legend(fontsize=7, frameon=False)
    fig.suptitle(f"{title}: data and models (grey = the other models)", fontsize=10, x=0.01, ha="left")
    fig.tight_layout()
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("folders", nargs="+")
    parser.add_argument("--out", default="dwelltime_report.pdf")
    parser.add_argument("--root", default=os.environ.get("SG_ROOT", "/data/Larsonlab/Nadya/StochasticGene"))
    args = parser.parse_args()
    fits, dists = read(args.folders, args.root)
    fits = rank(fits)
    out = Path(args.out)
    fits.drop(columns="label").to_csv(out.with_name(out.stem + "_model_ranking.csv"), index=False)
    with PdfPages(out) as pdf:
        for _, group in fits.groupby("data", sort=False):
            title = title_of(group)
            for figure in (ranking_page(group, title), distribution_page(group, dists, title)):
                pdf.savefig(figure)
                plt.close(figure)
    print(f"wrote {out} and {out.with_name(out.stem + '_model_ranking.csv')}")


if __name__ == "__main__":
    main()
