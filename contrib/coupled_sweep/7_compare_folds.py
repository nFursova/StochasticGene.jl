#!/usr/bin/env python3
"""Optional. The folds the sweep needs to make the experimental cross-correlation shape, against the folds seen in
the decoded traces; and the list of architectures to fit as coupled models.

    python3 7_compare_folds.py config.toml RELATIONSHIPS_FOLDER [--selection shape] [--observable reporter]
            [--estimate simple] [--profile full]

RELATIONSHIPS_FOLDER is the output folder of 6_trace_relationships.py (its rate_ratios.csv). The run folder needs
<selection>_passing_<observable>.csv from 4_plot_results.py: every simulation inside the shape selection.

How it works, per connection component (for example "CCAT1 state 2 -> MYC R23"):
- Sweep: every simulation in the selection that contains the component (one-way or reciprocal), with the fold that
  simulation used for that rate (its own fold per rate; in runs with one fold per connection, that fold).
- Data: the rate of that MYC (or CCAT1) transition while CCAT1 (or MYC) is in that state, divided by the rate
  otherwise, from the decoded traces (--estimate simple = one formula over whole cells, with a bootstrap interval;
  --estimate regression = the drift-adjusted regression of 6_trace_relationships.py).
- A simulation is "consistent with the data" when every one of its gene-to-gene folds lies inside the data's 95%
  interval. Hidden-unit connections have no data counterpart (the hidden unit is not observed): those models are
  marked for validation with exported traces (3_run_simulations.jl --export-traces) instead.

Writes into the run folder: folds_vs_data.html, folds_vs_data.png/.pdf, folds_vs_data_components.csv (one row per
component) and fit_candidates.csv (one row per architecture with at least one simulation in the selection: how many
simulations pass, whether any is consistent with the data, the best fit, and the parameters of the best simulation
as starting values for a coupled fit, with StochasticGene's coupling gamma = fold - 1).

Caveats: decoded states come from the uncoupled model, so the data folds are pulled towards 1 (a true 1/30 can look
like 0.8); the data estimate for source R is "the regulator has any reporters", while the sweep's R multiplies the
effect by the number of occupied reporter positions; state subsets of expanded runs (1+2, ...) have no data estimate.
"""
import argparse, html, tomllib
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

GENE = {"E": "CCAT1", "G": "MYC"}
COLOURS = {"E_to_G": "#2a78d6", "G_to_E": "#eb6834", "reciprocal": "#1baf7a"}
INK, MUTED = "#0b0b0b", "#898781"


def run_folder(config, folder, profile):
    if "profiles" in config:
        return folder / config["output"]["run_root"] / profile
    return folder / config["output"]["run_dir"]


def components(sim):
    """The gene-to-gene components of one simulation: [(label, fold)], label like 'E2 -> G R23'."""
    out = []
    for edge, source, target in (("EG", "E", "G"), ("GE", "G", "E")):
        rates, state = sim.get(f"{edge}_rates"), sim.get(f"{edge}_state")
        if not isinstance(rates, str) or not rates:
            continue
        for rate in rates.split("+"):
            fold = sim.get(f"fold_{edge}_{rate}")
            fold = sim.get(f"fold_{edge}") if fold is None or not np.isfinite(fold) else fold
            out.append((f"{source}{state} -> {target} {rate}", float(fold)))
    return out


def has_hidden(sim):
    return any(isinstance(sim.get(f"{edge}_rates"), str) and sim.get(f"{edge}_rates") for edge in ("HE", "HG"))


def describe(label):
    source, target = label.split(" -> ")
    rate = target.split(" ")[1]
    state = source[1:]
    return f"{GENE[source[0]]} {'with reporters' if state == 'R' else 'state ' + state} -> {GENE[target[0]]} {rate}"


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config")
    parser.add_argument("relationships", help="output folder of 6_trace_relationships.py")
    parser.add_argument("--selection", default="shape", choices=["shape", "shape_nrmse"])
    parser.add_argument("--observable", default="reporter", choices=["reporter", "binary"])
    parser.add_argument("--estimate", default="simple", choices=["simple", "regression"])
    parser.add_argument("--profile", default="full", help="expanded runs only")
    arguments = parser.parse_args()
    config_path = Path(arguments.config).resolve()
    config = tomllib.loads(config_path.read_text())
    run_dir = run_folder(config, config_path.parent, arguments.profile)
    obs = arguments.observable
    passing = pd.read_csv(run_dir / f"{arguments.selection}_passing_{obs}.csv", low_memory=False,
                          dtype={c: str for c in ("EG_state", "GE_state", "HE_state", "HG_state")})
    data = pd.read_csv(Path(arguments.relationships) / "rate_ratios.csv")
    key = {"simple": ("ratio_simple", "low_simple", "high_simple"),
           "regression": ("ratio_within_cells", "low_within", "high_within")}[arguments.estimate]
    data = data.set_index("architecture")[list(key)].set_axis(["estimate", "low", "high"], axis=1)
    fit = f"nrmse_{obs}_EG"

    # every (simulation, component, fold)
    rows, consistent, hidden = [], {}, {}
    for _, sim in passing.iterrows():
        sim = sim.to_dict()
        parts = components(sim)
        hidden[sim["sim_id"]] = has_hidden(sim)
        ok = bool(parts)
        for label, fold in parts:
            rows.append(dict(sim_id=sim["sim_id"], model_id=sim["model_id"], group=sim["group"], component=label, fold=fold))
            if label in data.index and np.isfinite(data.loc[label, "low"]) and np.isfinite(data.loc[label, "high"]):
                ok &= data.loc[label, "low"] <= fold <= data.loc[label, "high"]
            else:
                ok = False                              # no usable data estimate for this component
        consistent[sim["sim_id"]] = ok
    folds = pd.DataFrame(rows, columns=["sim_id", "model_id", "group", "component", "fold"])

    # one row per component
    summary = []
    for label, part in folds.groupby("component"):
        d = data.loc[label] if label in data.index else pd.Series(dict(estimate=np.nan, low=np.nan, high=np.nan))
        f = part["fold"].to_numpy()
        inside = ((f >= d["low"]) & (f <= d["high"])).mean() if np.isfinite(d["low"]) and np.isfinite(d["high"]) else np.nan
        verdict = ("no data estimate" if not np.isfinite(d["estimate"]) else
                   "overlaps the data" if inside > 0 else
                   "sweep needs a smaller fold than the data show" if np.median(f) < d["low"] else
                   "sweep needs a larger fold than the data show")
        summary.append(dict(component=label, meaning=describe(label), simulations=len(part), architectures=part["model_id"].nunique(),
                            sweep_median=np.median(f), sweep_p10=np.percentile(f, 10), sweep_p90=np.percentile(f, 90),
                            data_estimate=d["estimate"], data_low=d["low"], data_high=d["high"],
                            fraction_inside_data_interval=inside, verdict=verdict))
    summary = pd.DataFrame(summary).sort_values("simulations", ascending=False)
    summary.to_csv(run_dir / "folds_vs_data_components.csv", index=False)

    # one row per architecture: the candidate list
    passing["consistent"] = passing["sim_id"].map(consistent)
    passing["hidden"] = passing["sim_id"].map(hidden)
    parameters = [c for c in passing.columns if c.startswith("fold_") or c in ("h12", "h21", "h23", "h32")]
    candidates = []
    for model_id, part in passing.groupby("model_id"):
        part = part.sort_values(fit)
        good = part[part["consistent"]]
        start = (good if len(good) else part).iloc[0]
        values = {c: start[c] for c in parameters if pd.notna(start[c])}
        gammas = "; ".join(f"{label}: gamma {fold - 1:+.3g}" for label, fold in components(start.to_dict()))
        if start["hidden"]:
            status = "hidden unit: validate with exported traces"
        elif len(good):
            status = "consistent with the data"
        elif not components(start.to_dict()):
            status = "no gene-to-gene connection"
        else:
            status = "needs folds outside the data's intervals"
        candidates.append(dict(model_id=model_id, group=start["group"], status=status, simulations_in_selection=len(part),
                               simulations_consistent=len(good), best_nrmse_EG=part[fit].iloc[0],
                               best_consistent_nrmse_EG=good[fit].iloc[0] if len(good) else np.nan,
                               start_sim_id=start["sim_id"], start_from="best consistent" if len(good) else "best in selection",
                               start_gammas=gammas, **{f"start_{c}": v for c, v in values.items()}))
    order = {"consistent with the data": 0, "hidden unit: validate with exported traces": 1,
             "needs folds outside the data's intervals": 2, "no gene-to-gene connection": 3}
    candidates = pd.DataFrame(candidates)
    if len(candidates):
        candidates = candidates.sort_values(["status", "best_nrmse_EG"], key=lambda c: c.map(order) if c.name == "status" else c)
    candidates.to_csv(run_dir / "fit_candidates.csv", index=False)

    # plot: per component, the sweep's folds (points, coloured by group) against the data interval (black)
    shown = summary.head(40)
    figure, axis = plt.subplots(figsize=(11, 0.32 * len(shown) + 1.6))
    rng = np.random.default_rng(0)
    for i, label in enumerate(shown["component"]):
        y = len(shown) - 1 - i
        part = folds[folds["component"] == label]
        part = part if len(part) <= 400 else part.sample(400, random_state=0)
        for group, sub in part.groupby("group"):
            axis.scatter(sub["fold"], y + rng.uniform(-0.25, 0.25, len(sub)), s=6, alpha=0.5, color=COLOURS.get(group, MUTED),
                         linewidths=0, label=group if i == 0 else None)
        if label in data.index and np.isfinite(data.loc[label, "estimate"]):
            d = data.loc[label]
            axis.plot([d["low"], d["high"]], [y, y], color=INK, linewidth=2.2, solid_capstyle="butt")
            axis.plot([d["estimate"]], [y], marker="|", color=INK, markersize=12, markeredgewidth=2)
    axis.axvline(1, color=MUTED, linewidth=0.7)
    axis.set_xscale("log")
    labels = [f"{describe(c)}  (n={n:,})" for c, n in zip(shown["component"], shown["simulations"])]
    axis.set_yticks(range(len(shown))[::-1], labels, fontsize=7)
    axis.set_xlabel("fold: below 1 inhibits, above 1 activates", fontsize=8)
    handles = [plt.Line2D([], [], marker="o", linestyle="", color=c, label=g) for g, c in COLOURS.items()]
    handles.append(plt.Line2D([], [], color=INK, linewidth=2.2, label=f"data: {arguments.estimate} estimate, 95% interval"))
    axis.legend(handles=handles, fontsize=7, frameon=False, loc="lower right")
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    figure.suptitle(f"Folds that make the experimental shape (points: every simulation in the {arguments.selection} selection) "
                    "against the folds in the decoded traces (black)", x=0.01, ha="left", fontsize=9.5, fontweight="bold")
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(run_dir / "folds_vs_data.pdf")
    figure.savefig(run_dir / "folds_vs_data.png", dpi=300)
    plt.close(figure)

    counts = candidates["status"].value_counts().to_dict() if len(candidates) else {}
    fmt = lambda x: f"{x:.3g}" if isinstance(x, float) else x
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8"><title>Sweep folds and data</title>
<style>body{{font:15px/1.5 system-ui,sans-serif;max-width:1200px;margin:24px auto;padding:0 16px;color:#0b0b0b;background:#fcfcfb}}
img{{max-width:100%;border:1px solid #e1e0d9}}table{{border-collapse:collapse;font-size:12.5px}}td,th{{border-bottom:1px solid #e1e0d9;
padding:3px 7px;text-align:left}}p{{color:#52514e}}</style></head><body><h1>Folds the sweep needs, and folds in the data</h1>
<p>{len(passing):,} simulations in the {arguments.selection} selection ({obs}), {passing['model_id'].nunique():,} architectures. Data:
{arguments.estimate} estimate from {html.escape(str(arguments.relationships))}. Architectures: {', '.join(f'{v} {k}' for k, v in counts.items())}.</p>
<p>Each row of the plot is one connection component. Points: the fold of that component in every simulation inside the shape
selection (coloured by model group). Black: the data's fold and 95% interval. Decoding with the uncoupled model pulls data
folds towards 1, so a sweep that needs a stronger fold than the data show is not ruled out by this alone.</p>
<img src="folds_vs_data.png" alt="folds against data">
<h2>Components</h2>{summary.to_html(index=False, float_format=lambda x: f"{x:.3g}")}
<h2>Architectures to fit (fit_candidates.csv)</h2>
<p>Start values: the best simulation that is consistent with the data, otherwise the best in the selection; gamma = fold - 1.</p>
{candidates.drop(columns=[c for c in candidates.columns if c.startswith("start_fold") or c.startswith("start_h")], errors="ignore").head(200).to_html(index=False, float_format=lambda x: f"{x:.3g}")}
</body></html>"""
    (run_dir / "folds_vs_data.html").write_text(page)
    print(f"{len(passing):,} simulations, {len(candidates)} architectures: {counts}")
    print(f"Wrote {run_dir / 'folds_vs_data.html'}, folds_vs_data_components.csv and fit_candidates.csv")


if __name__ == "__main__":
    main()
