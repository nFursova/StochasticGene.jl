#!/usr/bin/env python3
"""Optional. The simulations of one architecture whose coupling fold is nearest a given value (for example the
coupling a direct StochasticGene fit found), drawn against the experiment and the uncoupled control.

    python3 5_nearest_simulations.py config.toml E_to_G__state2__R23 --gamma -0.062
    python3 5_nearest_simulations.py config.toml E_to_G__state2__R23 --fold 0.938 --param fold_EG_R23 --n 5
    python3 5_nearest_simulations.py config_expanded.toml <architecture_id> --gamma -0.062 --profile full

gamma is StochasticGene's coupling parameter: fold = 1 + gamma. --param names the fold column when the architecture
has more than one (reciprocal models, or several regulated rates in an extensive or expanded sweep); without it
the script uses the only fold column of the architecture. Needs ranked_scores.csv from 4_plot_results.py.

Writes into the run folder: nearest_<architecture>_<param>_<fold>.csv (the nearest simulations, with all their
parameters and fits) and .pdf / .png (their curves; thin lines = every other simulation of the architecture,
coloured by fold).
"""
import argparse, tomllib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

PAIRS = [("EG", "CCAT1 x MYC cross-correlation"), ("EE", "CCAT1 autocorrelation"), ("GG", "MYC autocorrelation")]
MAX_LAG_MIN = 120


def run_folder(config, folder, profile):
    """The run folder: run_dir of a regular config, run_root/<profile> of an expanded one."""
    if "profiles" in config:
        return folder / config["output"]["run_root"] / profile
    return folder / config["output"]["run_dir"]


def read_curves(run_dir, sim_ids, observable):
    """Mean and SE of the curves of `sim_ids`: DataFrame sim_id, pair, lag_min, mean, se (only the needed files are read)."""
    wanted, parts = set(int(i) for i in sim_ids), []
    results = run_dir / "results"
    binary = sorted(results.glob("curves_*.f32"))
    layout = pd.read_csv(results / "curve_layout.csv") if binary else None
    for path in binary or sorted(results.glob("curves_*.csv")):
        first, last = (int(part) for part in path.stem.split("_")[1].split("-"))
        inside = sorted(i for i in wanted if first <= i <= last)
        if not inside:
            continue
        if binary:                                       # Float32 [simulation, layout row, (mean, se)]
            values = np.fromfile(path, dtype=np.float32).reshape(last - first + 1, len(layout), 2)
            for i in inside:
                parts.append(layout.assign(sim_id=i, mean=values[i - first, :, 0], se=values[i - first, :, 1]))
        else:
            frame = pd.read_csv(path, usecols=["sim_id", "observable", "pair", "lag_min", "mean", "se"])
            parts.append(frame[frame["sim_id"].isin(inside)])
    curves = pd.concat(parts)
    return curves[(curves["observable"] == observable) & (curves["lag_min"].abs() <= MAX_LAG_MIN)]


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config")
    parser.add_argument("architecture", help="model_id as in models.csv / ranked_scores.csv, e.g. E_to_G__state2__R23")
    value = parser.add_mutually_exclusive_group(required=True)
    value.add_argument("--gamma", type=float, help="StochasticGene coupling parameter; fold = 1 + gamma")
    value.add_argument("--fold", type=float)
    parser.add_argument("--param", help="fold column, e.g. fold_EG or fold_EG_R23 (needed if the architecture has several)")
    parser.add_argument("--n", type=int, default=5, help="number of nearest simulations (default 5)")
    parser.add_argument("--observable", default="reporter", choices=["reporter", "binary"])
    parser.add_argument("--profile", default="full", help="expanded runs only")
    arguments = parser.parse_args()
    fold = arguments.fold if arguments.fold is not None else 1 + arguments.gamma
    if fold <= 0:
        raise SystemExit("the fold must be positive (gamma above -1)")

    config_path = Path(arguments.config).resolve()
    config = tomllib.loads(config_path.read_text())
    run_dir = run_folder(config, config_path.parent, arguments.profile)
    table = pd.read_csv(run_dir / "ranked_scores.csv", low_memory=False)
    model = table[table["model_id"] == arguments.architecture]
    if model.empty:
        close = [m for m in table["model_id"].unique() if arguments.architecture.split("__")[0] in m][:20]
        raise SystemExit(f"no architecture {arguments.architecture} in ranked_scores.csv. Some of that group: {close}")
    folds = [c for c in table.columns if c.startswith("fold_") and model[c].notna().any() and (model[c] != 1).any()]
    param = arguments.param or (folds[0] if len(folds) == 1 else None)
    if param not in folds:
        raise SystemExit(f"choose the fold column with --param; this architecture has: {', '.join(folds)}")

    model = model.assign(distance=np.abs(np.log(model[param]) - np.log(fold))).sort_values("distance")
    nearest = model.head(arguments.n)
    others = [c for c in folds + ["h12", "h21", "h23", "h32"] if c != param and c in model and model[c].notna().any()]
    fits = [f"nrmse_{arguments.observable}_EG", f"ccc_{arguments.observable}_EG", f"nrmse3_{arguments.observable}",
            f"center_{arguments.observable}", f"left_{arguments.observable}", f"right_{arguments.observable}",
            f"hat_shape_{arguments.observable}", "on_E", "on_G"]
    shown = ["sim_id", param] + others + [c for c in fits if c in model]
    print(f"{arguments.architecture}: {len(model)} simulations; {param} from {model[param].min():.4g} to {model[param].max():.4g}")
    print(f"nearest to {param} = {fold:.4g} (gamma {fold - 1:+.4g}):")
    print(nearest[shown].to_string(index=False, float_format=lambda x: f"{x:.4g}"))

    stem = f"nearest_{arguments.architecture}_{param}_{fold:.4g}"
    nearest.drop(columns="distance").to_csv(run_dir / f"{stem}.csv", index=False)

    uncoupled = table.loc[table["model_id"] == "uncoupled", "sim_id"].tolist()
    curves = read_curves(run_dir, model["sim_id"].tolist() + uncoupled, arguments.observable)
    experiment = pd.read_csv(config_path.parent / config["inputs"]["experimental_curves"])
    experiment = experiment[(experiment["observable"] == arguments.observable) & (experiment["lag_min"].abs() <= MAX_LAG_MIN)]

    figure, axes = plt.subplots(1, 3, figsize=(14, 4.4), constrained_layout=True)
    norm = matplotlib.colors.LogNorm(model[param].min(), model[param].max()) if model[param].min() != model[param].max() else None
    colors = plt.get_cmap("coolwarm")
    for axis, (pair, title) in zip(axes, PAIRS):
        part = curves[curves["pair"] == pair]
        for sim_id, value in zip(model["sim_id"], model[param]):     # every simulation of the architecture, thin
            if sim_id in set(nearest["sim_id"]):
                continue
            line = part[part["sim_id"] == sim_id]
            axis.plot(line["lag_min"], line["mean"], color=colors(norm(value)) if norm else "#c3c2b7", linewidth=0.5, alpha=0.35)
        shades = plt.get_cmap("viridis")(np.linspace(0.05, 0.75, max(len(nearest), 2)))   # nearest = darkest
        for k, (sim_id, value) in enumerate(zip(nearest["sim_id"], nearest[param])):   # the nearest, bold with SE
            line = part[part["sim_id"] == sim_id]
            axis.fill_between(line["lag_min"], line["mean"] - line["se"], line["mean"] + line["se"], color=shades[k], alpha=0.12)
            axis.plot(line["lag_min"], line["mean"], color=shades[k], linewidth=1.6,
                      label=f"sim {sim_id}: {param} {value:.3g} (γ {value - 1:+.3g})")
        for sim_id in uncoupled[:1]:
            line = part[part["sim_id"] == sim_id]
            axis.plot(line["lag_min"], line["mean"], color="#52514e", linestyle="--", linewidth=1.2, label="uncoupled")
        target = experiment[experiment["pair"] == pair]
        axis.errorbar(target["lag_min"], target["mean"], yerr=target["se"], color="black", linewidth=1.2, elinewidth=0.6, label="experiment")
        axis.axhline(0, color="#898781", linewidth=0.6)
        axis.set_title(title, fontsize=10)
        axis.set_xlabel("lag (min)" + ("; positive = MYC after CCAT1" if pair == "EG" else ""))
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
    axes[0].legend(fontsize=7, frameon=False)
    if norm:
        bar = figure.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=colors), ax=axes, shrink=0.8)
        bar.set_label(f"{param} of the thin lines (other simulations)")
    figure.suptitle(f"{arguments.architecture}, {arguments.observable}: the {len(nearest)} simulations nearest "
                    f"{param} = {fold:.3g} (γ {fold - 1:+.3g})", fontsize=11, fontweight="bold")
    figure.savefig(run_dir / f"{stem}.pdf")
    figure.savefig(run_dir / f"{stem}.png", dpi=300)
    print(f"Wrote {run_dir / stem}.csv, .pdf and .png")


if __name__ == "__main__":
    main()
