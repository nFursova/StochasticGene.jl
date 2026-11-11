#!/usr/bin/env python3
"""Step 4. Report the simulation results against the experimental data and rank the simulations.

    python3 4_plot_results.py config.toml

Reads <run_dir>/simulations.csv and <run_dir>/results/, and writes into <run_dir>:
    results_report.html   the report: an interactive explorer (choose a model, vary one parameter, see the
                          curves change), summary tables, and every plot below
    plots/                every plot as a PDF and as a 600 dpi PNG
    ranked_scores.csv     every simulation with its parameters, metrics and score, best first

Heatmap rows are simulations, grouped by model group, then model, then parameters.
"""
import html, json, sys, tomllib
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

PLOT_MAX_LAG_MIN = 120      # lag range shown in every plot
TOP_N = 10                  # number of best simulations shown as line plots and in the table
PNG_DPI = 600
SHOW_REPLICATES = True      # also draw every replicate as a thin line in the top-simulation plots

SCORE_NAME = "NRMSE of the reporter-count cross-correlation"


def score(metrics):
    """The ranking score of every simulation; LOWER IS BETTER. Edit this function (and SCORE_NAME) to design your own.

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
RATES = ["R12", "R21", "R23", "R32", "Initiate"]                      # the rates a connection can regulate
BASE_COLUMNS = [f"base_{gene}_{rate}" for gene in "EG" for rate in RATES]
INK, MUTED, BLUE = "#0b0b0b", "#898781", "#2a78d6"
DIVERGING = LinearSegmentedColormap.from_list("blue_gray_red", ["#104281", "#2a78d6", "#f0efec", "#e34948", "#8f1f1f"])
DIVERGING.set_bad("#c3c2b7")      # no data

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

def load(config_path):
    config = tomllib.loads(config_path.read_text())
    folder = config_path.parent
    run_dir = folder / config["output"]["run_dir"]
    read_all = lambda pattern: pd.concat([pd.read_csv(f) for f in sorted((run_dir / "results").glob(pattern))])
    simulations = pd.read_csv(run_dir / "simulations.csv", dtype={"EG_state": str, "GE_state": str})
    metrics, curves = read_all("metrics_*.csv"), read_all("curves_*.csv")
    experiment = pd.read_csv(folder / config["inputs"]["experimental_curves"])
    experiment_on = pd.read_csv(folder / config["inputs"]["experimental_metrics"]).set_index("metric")["mean"]
    # One row per simulation: parameters, metrics, score and rank (1 = best).
    table = simulations.merge(metrics, on=["sim_id", "model_id"]).sort_values("sim_id")
    table["score"] = score(table)
    table["rank"] = table["score"].rank(method="first").astype("Int64")
    for frame in (curves, experiment):
        frame.drop(frame[frame["lag_min"].abs() > PLOT_MAX_LAG_MIN].index, inplace=True)
    done = f'{len(table):,} of {len(simulations):,} simulations have results'
    print(done)
    statistic = config["inputs"]["rate_statistic"]
    fitted = read_rates(folder / config["inputs"]["rates_E"], statistic) + read_rates(folder / config["inputs"]["rates_G"], statistic)
    return config, run_dir, table, curves.reset_index(drop=True), experiment, experiment_on, done, fitted


def read_rates(path, statistic):
    """The fitted R12, R21, R23, R32 and Initiate of one gene (rows of a rates_*.txt file: ml, mean, median, last)."""
    lines = [line for line in Path(path).read_text().splitlines() if line.strip()]
    try:
        float(lines[0].split(",")[0])
    except ValueError:
        lines = lines[1:]                               # drop the header line
    row = {"ml": 0, "mean": 1, "median": 2, "last": 3}[statistic]
    return [float(value) for value in lines[row].split(",")[:len(RATES)]]


def regulated_rates(sim):
    """Which rates a model regulates, and which unit and state regulates each."""
    hidden_states = ("H state 2", "H state 2") if sim["hidden"] == "H2" else ("H state 1", "H state 3")
    regulated = []
    for column, gene, source in [("GE_rates", "E", f'G state {sim["GE_state"]}'), ("HE_rates", "E", hidden_states[0]),
                                 ("EG_rates", "G", f'E state {sim["EG_state"]}'), ("HG_rates", "G", hidden_states[1])]:
        if pd.notna(sim[column]):
            regulated += [dict(gene=gene, rate=rate, by=source, index=BASE_COLUMNS.index(f"base_{gene}_{rate}"))
                          for rate in sim[column].split("+")]
    return regulated


def curve_matrix(curves, experiment, observable, pair, column="mean"):
    """One curve for every simulation (rows = sim_id, columns = lag in minutes) and the experiment on the same lags.
    column: "mean" (over the replicates of a simulation) or "se" (standard error over the replicates)."""
    select = lambda frame: frame[(frame["observable"] == observable) & (frame["pair"] == pair)]
    values = select(curves).pivot(index="sim_id", columns="lag_min", values=column)
    frames = select(curves).drop_duplicates("lag_min").set_index("lag_min")["lag_frames"].reindex(values.columns)
    target = select(experiment).set_index("lag_frames").reindex(frames.to_numpy())
    return values, target["mean"].to_numpy(), target["se"].to_numpy()


def free_parameters(rows):
    """The parameters that vary among these simulations, leaving out those tied to an earlier parameter."""
    free = []
    for name in PARAMETERS:
        if rows[name].nunique() > 1 and not any(rows[name].equals(rows[other]) for other in free):
            free.append(name)
    return free


def parameter_text(sim):
    return ", ".join(f"{name} {sim[name]:g}" for name in PARAMETERS if pd.notna(sim[name]))

# ---------------------------------------------------------------------------
# Plots (each returns a matplotlib figure)
# ---------------------------------------------------------------------------

def bin_edges(centers):
    """Edges of heatmap cells centred on the (unevenly spaced) multi-tau lags."""
    middle = (centers[:-1] + centers[1:]) / 2
    return np.concatenate([[2 * centers[0] - middle[0]], middle, [2 * centers[-1] - middle[-1]]])


def heatmap(title, table, values, x_edges, limit, colorbar_label, top_row=None, x_label="", x_ticks=None):
    """A panel per model group, rows = that group's simulations in sim_id order; top_row = experiment."""
    groups = list(dict.fromkeys(table["group"]))
    heights = ([0.35] if top_row is not None else []) + [1] * len(groups)
    figure, axes = plt.subplots(len(heights), 1, figsize=(7.5, 9.5), gridspec_kw=dict(height_ratios=heights, hspace=0.12),
                                squeeze=False)
    axes = list(axes[:, 0])
    style = dict(cmap=DIVERGING, vmin=-limit, vmax=limit, rasterized=True)
    if top_row is not None:
        mesh = axes[0].pcolormesh(x_edges, [0, 1], top_row[None, :], **style)
        axes[0].set_ylabel("experiment", rotation=0, ha="right", va="center", fontweight="bold")
    for axis, group in zip(axes[len(heights) - len(groups):], groups):
        rows = table[table["group"] == group]
        mesh = axis.pcolormesh(x_edges, np.arange(len(rows) + 1), values.reindex(rows["sim_id"]).to_numpy(), **style)
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
    figure.suptitle(title, x=0.06, ha="left", fontsize=12, fontweight="bold")
    figure.subplots_adjust(left=0.27, right=0.86, top=0.94, bottom=0.06)
    figure.colorbar(mesh, cax=figure.add_axes([0.89, 0.3, 0.018, 0.4]), label=colorbar_label, extend="both")
    return figure


def curve_heatmap(table, curves, experiment, observable, pair, name):
    values, target, _ = curve_matrix(curves, experiment, observable, pair)
    limit = 2 * np.nanmax(np.abs(target))        # colour scale set by the experiment; stronger simulations saturate
    return heatmap(name, table, values, bin_edges(values.columns.to_numpy()), limit,
                   "correlation (blue negative, red positive)", top_row=target,
                   x_label="lag (min); positive = MYC after CCAT1" if pair == "EG" else "lag (min)")


def on_fraction_heatmap(table, experiment_on):
    ratio = pd.DataFrame({"E": np.log2(table["on_E"] / experiment_on["on_CCAT1"]),
                          "G": np.log2(table["on_G"] / experiment_on["on_MYC"])}).set_index(table["sim_id"])
    limit = max(1.0, np.nanpercentile(np.abs(ratio.to_numpy()), 99))
    title = (f'ON fraction relative to experiment (experiment: CCAT1 {experiment_on["on_CCAT1"]:.3f}, '
             f'MYC {experiment_on["on_MYC"]:.3f})')
    return heatmap(title, table, ratio, np.array([0.0, 1.0, 2.0]), limit,
                   "log2(simulated / experimental); 0 = equal, +1 = twice as often ON",
                   x_ticks=([0.5, 1.5], ["CCAT1 (E)", "MYC (G)"]))


def top_curves(table, curves, experiment, experiment_on):
    """The TOP_N best simulations: every curve against experiment and the uncoupled control, one row per simulation."""
    best = table.sort_values("rank").head(TOP_N)
    control = table.loc[table["group"] == "uncoupled", "sim_id"]
    figure, axes = plt.subplots(len(best), len(CURVES), figsize=(15, 1.25 * len(best) + 1.4), sharex=True, sharey="col",
                                squeeze=False)
    for row, (_, sim) in enumerate(best.iterrows()):
        label = (f'#{sim["rank"]}  simulation {sim["sim_id"]}  score {sim["score"]:.3f}\n{sim["model_id"]}\n{parameter_text(sim)}\n'
                 f'ON: CCAT1 {sim["on_E"]:.3f} (exp. {experiment_on["on_CCAT1"]:.3f}), '
                 f'MYC {sim["on_G"]:.3f} (exp. {experiment_on["on_MYC"]:.3f})')
        axes[row, 0].set_ylabel(label, rotation=0, ha="right", va="center", fontsize=6.5)
        for column, (observable, pair, name) in enumerate(CURVES):
            axis = axes[row, column]
            select = lambda frame: frame[(frame["observable"] == observable) & (frame["pair"] == pair)].sort_values("lag_min")
            if len(control):
                data = select(curves[curves["sim_id"] == control.iloc[0]])
                axis.plot(data["lag_min"], data["mean"], color=MUTED, linewidth=0.9, linestyle="--", label="uncoupled control")
            if SHOW_REPLICATES:
                data = select(curves[curves["sim_id"] == sim["sim_id"]])
                for k, replicate in enumerate(c for c in curves.columns if c.startswith("rep")):
                    axis.plot(data["lag_min"], data[replicate], color=BLUE, linewidth=0.5, alpha=0.55,
                              label="individual replicates" if k == 0 else None)
            for data, color, text in [(select(experiment), INK, "experiment ± SE"),
                                      (select(curves[curves["sim_id"] == sim["sim_id"]]), BLUE, "simulation mean ± SE over replicates")]:
                axis.fill_between(data["lag_min"], data["mean"] - data["se"], data["mean"] + data["se"],
                                  color=color, alpha=0.18, linewidth=0)
                axis.plot(data["lag_min"], data["mean"], color=color, linewidth=1.3, label=text)
            axis.axhline(0, color=MUTED, linewidth=0.5)
            axis.tick_params(labelsize=6.5, colors=MUTED, labelcolor=INK)
            for side in ("top", "right"):
                axis.spines[side].set_visible(False)
            if row == 0:
                axis.set_title(name, fontsize=7.5)
            if row == len(best) - 1:
                axis.set_xlabel("lag (min)", fontsize=7.5)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper right", ncol=4, frameon=False, fontsize=8)
    figure.suptitle(f"The {len(best)} best simulations (lowest score)", x=0.02, ha="left", fontsize=11, fontweight="bold")
    height = figure.get_figheight()
    figure.subplots_adjust(left=0.2, right=0.99, top=1 - 0.8 / height, bottom=0.5 / height, wspace=0.24, hspace=0.28)
    return figure


def score_by_group(table):
    """Every simulation's score, one row of points per model group; the bar marks the group median."""
    groups = list(dict.fromkeys(table["group"]))
    figure, axis = plt.subplots(figsize=(7.5, 0.5 * len(groups) + 1.3))
    jitter = np.random.default_rng(0)
    for y, group in enumerate(groups):
        scores = table.loc[table["group"] == group, "score"].dropna()
        axis.scatter(scores, y + jitter.uniform(-0.25, 0.25, len(scores)), s=7, color=BLUE, alpha=0.35, linewidths=0)
        axis.plot([scores.median()] * 2, [y - 0.36, y + 0.36], color=INK, linewidth=1.6)
    axis.set_yticks(range(len(groups)), groups)
    axis.invert_yaxis()
    axis.set_xlabel(f"score: {SCORE_NAME} (lower is better); black bar = group median")
    axis.grid(axis="x", color="#e1e0d9", linewidth=0.6)
    axis.set_axisbelow(True)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    figure.tight_layout()
    return figure


def score_vs_parameters(table):
    """Score against each swept parameter, one row of panels per model group; the line joins the medians."""
    groups = [g for g in dict.fromkeys(table["group"]) if free_parameters(table[table["group"] == g])]
    figure, axes = plt.subplots(len(groups), len(PARAMETERS), figsize=(11, 1.55 * len(groups) + 0.9), sharey=True, squeeze=False)
    for row, group in enumerate(groups):
        rows = table[table["group"] == group]
        free = free_parameters(rows)
        for column, name in enumerate(PARAMETERS):
            axis = axes[row, column]
            if name not in free:
                axis.set_visible(False)
                continue
            axis.scatter(rows[name], rows["score"], s=6, color=BLUE, alpha=0.3, linewidths=0)
            if rows[name].nunique() <= 12:
                median = rows.groupby(name)["score"].median()
                axis.plot(median.index, median.to_numpy(), color=INK, linewidth=1.2, marker="o", markersize=3)
            axis.set_xscale("log")
            axis.set_title(f"{group}: {name}", fontsize=8)
            axis.tick_params(labelsize=7, colors=MUTED, labelcolor=INK)
            for side in ("top", "right"):
                axis.spines[side].set_visible(False)
        axes[row, 0].set_ylabel("score", fontsize=8)
    figure.tight_layout()
    return figure

# ---------------------------------------------------------------------------
# HTML report
# ---------------------------------------------------------------------------

def table_html(header, rows):
    row_html = lambda cells, tag: "<tr>" + "".join(f"<{tag}>{html.escape(str(cell))}</{tag}>" for cell in cells) + "</tr>"
    return '<div class="scroll"><table>' + row_html(header, "th") + "".join(row_html(row, "td") for row in rows) + "</table></div>"


def group_summary(table):
    rows = []
    for group in dict.fromkeys(table["group"]):
        part = table[table["group"] == group]
        best = part.loc[part["score"].idxmin()]
        rows.append([group, len(part), f'{best["score"]:.3f}', f'{part["score"].median():.3f}',
                     f'simulation {best["sim_id"]}: {best["model_id"]}; {parameter_text(best)}'])
    return table_html(["Group", "Simulations", "Best score", "Median score", "Best simulation"], rows)


def top_table(table, experiment_on):
    rows = [[sim["rank"], sim["sim_id"], sim["model_id"], parameter_text(sim), f'{sim["score"]:.3f}',
             f'{sim["nrmse_binary_EG"]:.3f} / {sim["ccc_binary_EG"]:.3f}', f'{sim["nrmse_reporter_EG"]:.3f} / {sim["ccc_reporter_EG"]:.3f}',
             f'{sim["on_E"]:.3f}', f'{sim["on_G"]:.3f}'] for _, sim in table.sort_values("rank").head(TOP_N).iterrows()]
    return table_html(["Rank", "Simulation", "Model", "Parameters", "Score", "x-corr state NRMSE / CCC",
                       "x-corr reporter NRMSE / CCC", f'ON CCAT1 (exp. {experiment_on["on_CCAT1"]:.3f})',
                       f'ON MYC (exp. {experiment_on["on_MYC"]:.3f})'], rows)


def explorer_data(table, curves, experiment, fitted, matching):
    """Everything the interactive explorer needs, as one JSON object."""
    clean = lambda values: [None if v != v else round(float(v), 4) for v in values]
    matrices = [curve_matrix(curves, experiment, observable, pair) for observable, pair, _ in CURVES]
    # Every replicate of every curve; the page computes mean and SE from them.
    replicates = [[curve_matrix(curves, experiment, observable, pair, column)[0].reindex(table["sim_id"]).to_numpy()
                   for column in curves.columns if column.startswith("rep")] for observable, pair, _ in CURVES]
    sims = [dict(id=int(sim["sim_id"]), model=sim["model_id"],
                 par={name: float(sim[name]) for name in PARAMETERS if pd.notna(sim[name])},
                 nrmse=clean(sim[f"nrmse_{observable}_{pair}"] for observable, pair, _ in CURVES),
                 base=[float(f"{sim[column]:.4g}") for column in BASE_COLUMNS],
                 reps=[[clean(values[i]) for values in one_curve] for one_curve in replicates])
            for i, (_, sim) in enumerate(table.iterrows())]
    models = {model: dict(group=rows["group"].iloc[0], free=free_parameters(rows), regulated=regulated_rates(rows.iloc[0]))
              for model, rows in table.groupby("model_id", sort=False)}
    control = table.loc[table["group"] == "uncoupled", "sim_id"]
    return dict(curve_names=[name for _, _, name in CURVES], lags=clean(matrices[0][0].columns),
                experiment=[clean(m[1]) for m in matrices], experiment_se=[clean(m[2]) for m in matrices],
                models=models, sims=sims, control=int(control.iloc[0]) if len(control) else None, fitted=fitted,
                state_average=matching == "state_average")


STYLE = """
:root{--bg:#f9f9f7;--card:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--line:#e1e0d9}
@media (prefers-color-scheme: dark){:root{--bg:#0d0d0d;--card:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--line:#2c2c2a}}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1180px;margin:auto;padding:24px 16px 48px} h1{font-size:26px;margin:0 0 4px} h2{font-size:19px;margin:30px 0 6px}
h3{font-size:15px;margin:0 0 2px} p{margin:6px 0} .sub{color:var(--ink2)} a{color:#2a78d6}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:12px 0}
.scroll{overflow-x:auto} table{border-collapse:collapse;font-size:13px;width:100%}
th,td{text-align:left;padding:5px 10px;border-bottom:1px solid var(--line);font-weight:400;font-variant-numeric:tabular-nums}
th{color:var(--ink2)} figure{margin:0} figure img{width:100%;height:auto;background:#fff;border-radius:4px}
.controls{display:flex;flex-wrap:wrap;gap:10px 14px;margin:8px 0} label{display:flex;flex-direction:column;font-size:12px;color:var(--ink2)}
select{font:inherit;font-size:14px;padding:5px 6px;max-width:340px;color:var(--ink);background:var(--card);border:1px solid var(--line);border-radius:6px}
#chart{width:100%;height:auto;display:block} #chart text{font:12px system-ui,sans-serif;fill:var(--ink2)}
#legend{display:flex;flex-wrap:wrap;gap:4px 18px;font-size:13px;margin-top:6px} #legend i{display:inline-block;width:22px;height:3px;vertical-align:middle;margin-right:6px}
"""

# The interactive explorer: choose a curve, a model and the parameter to vary; the other parameters are held fixed.
EXPLORER_JS = r"""
const $ = id => document.getElementById(id), NAMES = ["fold", "h12", "h21", "h23", "h32"];
const groups = [...new Set(Object.values(DATA.models).map(m => m.group))];
const options = (select, values, labels) => { select.innerHTML = values.map((v, i) => `<option value="${v}">${labels ? labels[i] : v}</option>`).join(""); };
const distinct = (sims, name) => [...new Set(sims.map(s => s.par[name]))].sort((a, b) => a - b);
const shade = t => { const a = [158, 197, 244], b = [13, 54, 107]; return `rgb(${a.map((v, i) => Math.round(v + (b[i] - v) * t)).join(",")})`; };

function modelChanged() {                       // which parameter to vary
  options($("vary"), DATA.models[$("model").value].free);
  varyChanged();
}
function varyChanged() {                        // one selector for every other free parameter (if it has few values)
  const sims = DATA.sims.filter(s => s.model === $("model").value), vary = $("vary").value;
  $("fixed").innerHTML = "";
  for (const name of DATA.models[$("model").value].free.filter(n => n !== vary)) {
    const values = distinct(sims, name);
    if (values.length > 12) continue;
    const label = document.createElement("label");
    label.innerHTML = `hold ${name} at<select data-name="${name}"></select>`;
    options(label.firstElementChild, values);
    label.firstElementChild.onchange = draw;
    $("fixed").appendChild(label);
  }
  draw();
}
function draw() {
  const curve = +$("curve").value, vary = $("vary").value, show = $("show").value;
  let sims = DATA.sims.filter(s => s.model === $("model").value);
  for (const select of $("fixed").querySelectorAll("select")) sims = sims.filter(s => s.par[select.dataset.name] == select.value);
  if (vary) sims.sort((a, b) => a.par[vary] - b.par[vary]);
  const lags = DATA.lags, target = DATA.experiment[curve], se = DATA.experiment_se[curve];
  const control = DATA.sims.find(s => s.id === DATA.control);
  const all = [...target.map((v, i) => v + (se[i] || 0)), ...target.map((v, i) => v - (se[i] || 0)), 0,
               ...sims.flatMap(s => show === "replicates" ? s.reps[curve].flat() : s.curves[curve]),
               ...(control ? control.curves[curve] : [])].filter(v => v !== null);
  const few = sims.length <= 12;                // few simulations: draw each one's mean ± SE over its replicates
  let lo = Math.min(...all), hi = Math.max(...all); const pad = (hi - lo || 1) * 0.06; lo -= pad; hi += pad;
  const W = 900, H = 380, L = 56, R = 12, T = 10, B = 40;
  const x = t => L + (W - L - R) * (t - lags[0]) / (lags[lags.length - 1] - lags[0]), y = v => T + (H - T - B) * (hi - v) / (hi - lo);
  const path = values => values.map((v, i) => v === null ? "" : `${x(lags[i])},${y(v)}`).filter(Boolean).join(" ");
  let svg = `<line x1="${L}" x2="${W - R}" y1="${y(0)}" y2="${y(0)}" stroke="var(--line)"/>`;
  for (let k = 0; k <= 4; k++) { const v = lo + (hi - lo) * k / 4; svg += `<text x="${L - 6}" y="${y(v) + 4}" text-anchor="end">${v.toPrecision(2)}</text>`; }
  for (const t of [-120, -90, -60, -30, 0, 30, 60, 90, 120].filter(t => t >= lags[0] && t <= lags[lags.length - 1]))
    svg += `<text x="${x(t)}" y="${H - B + 16}" text-anchor="middle">${t}</text>`;
  svg += `<text x="${(L + W - R) / 2}" y="${H - 6}" text-anchor="middle">lag (min)</text>`;
  const band = (mean, se, color, opacity) => `<polygon points="${mean.map((v, i) => `${x(lags[i])},${y(v + (se[i] || 0))}`)
    .concat(mean.map((v, i) => `${x(lags[i])},${y(v - (se[i] || 0))}`).reverse()).join(" ")}" fill="${color}" opacity="${opacity}"/>`;
  svg += band(target, se, "var(--ink)", 0.12);
  if (few && show === "se") sims.forEach((s, i) => { svg += band(s.curves[curve], s.se[curve], shade(sims.length > 1 ? i / (sims.length - 1) : 0.6), 0.16); });
  if (control) svg += `<polyline points="${path(control.curves[curve])}" fill="none" stroke="var(--muted)" stroke-width="1.5" stroke-dasharray="5 4"/>`;
  sims.forEach((s, i) => {
    const color = shade(sims.length > 1 ? i / (sims.length - 1) : 0.6);
    if (show === "replicates") for (const replicate of s.reps[curve]) svg += `<polyline points="${path(replicate)}" fill="none" stroke="${color}" stroke-width="1"/>`;
    else svg += `<polyline points="${path(s.curves[curve])}" fill="none" stroke="${color}" stroke-width="${few ? 2 : 1}"/>`;
  });
  svg += `<polyline points="${path(target)}" fill="none" stroke="var(--ink)" stroke-width="2.4"/>`;
  $("chart").innerHTML = svg;
  const item = (color, text) => `<span><i style="background:${color}"></i>${text}</span>`;
  let legend = item("var(--ink)", "experiment ± SE") + (control ? item("var(--muted)", "uncoupled control") : "");
  if (sims.length > 12) legend += item(shade(0), `${sims.length} simulations, ${vary} from ${sims[0].par[vary].toPrecision(3)} (light)`) + item(shade(1), `to ${sims[sims.length - 1].par[vary].toPrecision(3)} (dark)`);
  else sims.forEach((s, i) => { legend += item(shade(sims.length > 1 ? i / (sims.length - 1) : 0.6), `${vary ? vary + " " + s.par[vary] : "fold 1"} · simulation ${s.id} · NRMSE ${s.nrmse[curve] === null ? "n/a" : s.nrmse[curve].toFixed(3)}`); });
  $("legend").innerHTML = legend;
  rateTable(sims, vary);
}
function rateTable(sims, vary) {                // the regulated rates of the drawn simulations
  const regulated = DATA.models[$("model").value].regulated, f = v => v.toPrecision(3);
  if (!regulated.length || sims.length > 12) { $("rates").innerHTML = regulated.length ? "<p class='sub'>The rate table is shown when 12 or fewer simulations are drawn.</p>" : ""; return; }
  let rows = `<tr><th>Rate</th><th>Regulated by</th><th>Fitted (average)</th>` +
    sims.map(s => `<th>${vary ? vary + " " + s.par[vary] : "simulation " + s.id}: base → regulated (source active)</th>`).join("") + "</tr>";
  for (const r of regulated) {
    const fitted = DATA.fitted[r.index];
    rows += `<tr><td>${r.gene === "E" ? "CCAT1" : "MYC"} ${r.rate}</td><td>${r.by}</td><td>${f(fitted)}</td>` + sims.map(s => {
      const base = s.base[r.index], fold = s.par.fold, activity = (fitted / base - 1) / (fold - 1);
      const when = r.by.endsWith(" R") ? ` per occupied reporter position (mean ${activity.toFixed(2)} occupied)` : ` (${(100 * activity).toFixed(0)}% of the time)`;
      return `<td>${f(base)} → ${f(base * fold)}` + (fold > 1 ? when : "") + "</td>";
    }).join("") + "</tr>";
  }
  $("rates").innerHTML = `<table>${rows}</table><p class="sub">Rates per minute. Base = the rate while its source is inactive; regulated = base × fold. For source state R every occupied reporter position of the source gene adds (fold − 1) × base.
    ${DATA.state_average ? "The base is set so that the average rate, over the time when the gene is in the state where that rate can act, equals the fitted rate; the source activity shown is over that time too."
      : "The base is set so that the average over all time equals the fitted rate."} All other rates equal their fitted values.</p>`;
}
for (const s of DATA.sims) {                   // mean and SE over the replicates, for every curve and lag
  const stats = s.reps.map(reps => reps[0].map((_, i) => {
    const v = reps.map(r => r[i]).filter(x => x !== null), n = v.length, mean = n ? v.reduce((a, b) => a + b, 0) / n : null;
    return [mean, n > 1 ? Math.sqrt(v.reduce((a, b) => a + (b - mean) ** 2, 0) / (n - 1) / n) : 0];
  }));
  s.curves = stats.map(c => c.map(p => p[0])); s.se = stats.map(c => c.map(p => p[1]));
}
options($("curve"), DATA.curve_names.map((_, i) => i), DATA.curve_names);
options($("group"), groups);
$("group").onchange = () => { options($("model"), Object.keys(DATA.models).filter(m => DATA.models[m].group === $("group").value)); modelChanged(); };
$("model").onchange = modelChanged; $("vary").onchange = varyChanged; $("curve").onchange = draw; $("show").onchange = draw;
$("curve").value = 1; $("group").value = groups[Math.min(1, groups.length - 1)]; $("group").onchange();
"""


def write_report(run_dir, config, done, table, curves, experiment, experiment_on, fitted, plots):
    matching = config["simulation"].get("rate_matching", "time_average")
    data = json.dumps(explorer_data(table, curves, experiment, fitted, matching), separators=(",", ":")).replace("</", "<\\/")
    figures = "".join(
        f'<section class="card"><figure><h3>{html.escape(title)}</h3><p class="sub">{html.escape(caption)} '
        f'<a href="plots/{stem}.pdf">PDF</a> · <a href="plots/{stem}.png">PNG</a></p>'
        f'<img src="plots/{stem}.png" alt="{html.escape(title)}" loading="lazy"></figure></section>' for stem, title, caption in plots)
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Coupled-model sweep results</title><style>{STYLE}</style></head>
<body><main><h1>Coupled-model sweep: results</h1>
<p class="sub">{done} · {config["simulation"]["replicates"]} replicate datasets each · simulator: {config["simulation"]["simulator"]} ·
sweep mode: {config["sweep"]["mode"]} · rate matching: {config["simulation"].get("rate_matching", "time_average")} · score = {html.escape(SCORE_NAME)} (lower is better) · lags within ±{PLOT_MAX_LAG_MIN} min shown.</p>

<h2>Explore: how a parameter changes the curves</h2>
<section class="card"><div class="controls">
<label>curve<select id="curve"></select></label><label>model group<select id="group"></select></label>
<label>model<select id="model"></select></label><label>vary<select id="vary"></select></label>
<label>show<select id="show"><option value="se">mean ± SE</option><option value="mean">mean only</option><option value="replicates">individual replicates</option></select></label><span id="fixed" class="controls" style="margin:0"></span></div>
<svg id="chart" viewBox="0 0 900 380" role="img" aria-label="selected curves against experiment"></svg><div id="legend"></div>
<div id="rates" class="scroll" style="margin-top:10px"></div>
<p class="sub">"show" switches between the mean of a simulation's replicate datasets with a ± SE band (band drawn when 12 or
fewer simulations are shown), the mean alone, or every replicate as its own line. Light to dark blue = increasing value of the varied parameter. If a model was sampled with many random
parameter sets (extensive sweep), all of them are drawn and the other parameters are not held fixed.</p></section>

<h2>Summary by model group</h2><section class="card">{group_summary(table)}</section>
<h2>The {TOP_N} best simulations</h2><section class="card">{top_table(table, experiment_on)}
<p class="sub">All simulations: <a href="ranked_scores.csv">ranked_scores.csv</a></p></section>
<h2>Plots</h2>{figures}
<script>const DATA={data};</script><script>{EXPLORER_JS}</script></main></body></html>"""
    (run_dir / "results_report.html").write_text(page)

# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    config, run_dir, table, curves, experiment, experiment_on, done, fitted = load(Path(sys.argv[1]).resolve())
    table.sort_values("rank").to_csv(run_dir / "ranked_scores.csv", index=False)
    folder = run_dir / "plots"
    folder.mkdir(exist_ok=True)
    plots = []      # (file name, title, caption) of every plot, in report order

    def save(figure, stem, title, caption):
        figure.savefig(folder / f"{stem}.pdf", dpi=300)
        figure.savefig(folder / f"{stem}.png", dpi=PNG_DPI)
        plt.close(figure)
        plots.append((stem, title, caption))

    save(top_curves(table, curves, experiment, experiment_on), "01_top_simulations", f"The {TOP_N} best simulations",
         "Every curve of each top simulation (blue) against experiment (black) and the uncoupled control (grey, dashed).")
    save(score_by_group(table), "02_score_by_group", "Score by model group", "One point per simulation.")
    save(score_vs_parameters(table), "03_score_vs_parameters", "Score against each swept parameter",
         "One point per simulation; the black line joins the medians.")
    for number, (observable, pair, name) in enumerate(CURVES, start=4):
        save(curve_heatmap(table, curves, experiment, observable, pair, name), f"{number:02d}_heatmap_{pair}_{observable}", name,
             "Each row is one simulation (the mean of its replicate datasets), grouped by model group, then model, "
             "then parameters. Top strip: experiment.")
    save(on_fraction_heatmap(table, experiment_on), "10_heatmap_on_fraction", "ON fraction relative to experiment",
         "Rows as in the curve heatmaps.")
    write_report(run_dir, config, done, table, curves, experiment, experiment_on, fitted, plots)
    print(f"Wrote {run_dir / 'results_report.html'}, {folder}/ and {run_dir / 'ranked_scores.csv'}")


if __name__ == "__main__":
    main()
