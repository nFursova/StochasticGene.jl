#!/usr/bin/env python3
"""Step 4. Report the simulation results against the experimental data and rank the simulations.

    python3 4_plot_results.py config.toml
    python3 4_plot_results.py config_expanded.toml --profile full      (expanded runs; full is the default)

Reads the simulation list (simulations.csv, or expanded_models.csv and expanded_simulations.csv) and
<run_dir>/results/, and writes into <run_dir>:
    results_report.html   the report: an interactive explorer (choose a model, vary one parameter, see the
                          curves change), summary tables, and every plot below
    plots/                every plot as a PDF and as a 600 dpi PNG
    explorer/             the data the report's explorer loads when you pick a model (one small file per model)
    ranked_scores.csv     every simulation with its parameters, metrics and score, best first

Heatmap rows are simulations, grouped by model group, then model, then parameters.
"""
import html, inspect, json, os, sys, tomllib
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

PLOT_MAX_LAG_MIN = 120      # lag range shown in every plot
TOP_N = 10                  # number of best simulations shown as line plots and in the table
TOP_ARCHITECTURES = 100     # number of architectures on the "top architectures" pages
PNG_DPI = 600
A4 = (8.27, 11.69)          # page size of the heatmaps, inches
HAT_CENTER_MIN = 10         # Mexican-hat shape: the cross-correlation should be positive within +/- this many minutes
HAT_FLANKS_MIN = (20, 60)   # and negative between these lags (minutes) on both sides

SCORE_NAME = "NRMSE of the reporter-count cross-correlation"
MAX_POINTS = 50_000         # scatter plots draw at most this many randomly chosen points per row (medians use all)


def score(metrics):
    """The ranking score of every simulation; LOWER IS BETTER. Edit this function (and SCORE_NAME) to design your own.

    `metrics` has one row per simulation with these columns:
        nrmse_<observable>_<pair>   0 = identical to experiment
        ccc_<observable>_<pair>     1 = identical to experiment
        on_E, on_G, reporters_E, reporters_G
        center_<observable>, left_<observable>, right_<observable>   mean cross-correlation within +/- HAT_CENTER_MIN
                                    and in the left and right flanks (HAT_FLANKS_MIN)
        hat_<observable>            Mexican-hat depth: center minus the mean of the two flanks
        hat_shape_<observable>      True if center > 0 and both flanks < 0
    observable: binary (ON/OFF state) or reporter (reporter count)
    pair:       EG (CCAT1 x MYC cross-correlation), EE (CCAT1 autocorrelation), GG (MYC autocorrelation)
    Example of a combined score:  metrics["nrmse_reporter_EG"] + metrics["nrmse_binary_EG"]
    """
    return metrics["nrmse_reporter_EG"]


CURVES = [(observable, pair, f"{pair_name}, {observable_name}")
          for pair, pair_name in [("EG", "CCAT1 x MYC cross-correlation"), ("EE", "CCAT1 autocorrelation"),
                                  ("GG", "MYC autocorrelation")]
          for observable, observable_name in [("binary", "ON/OFF state"), ("reporter", "reporter count")]]
# The metrics drawn per model group (violins) and against every swept parameter: (column, title, axis scale).
METRICS = [(f"{name}_{observable}{pair}", f"{title}, {observable_name}", scale)
           for observable, observable_name in [("reporter", "reporter count"), ("binary", "ON/OFF state")]
           for name, pair, title, scale in [
               ("nrmse", "_EG", "NRMSE of the cross-correlation (0 = identical to experiment)", "log"),
               ("ccc", "_EG", "CCC of the cross-correlation (1 = identical to experiment)", "linear"),
               ("hat", "", "Mexican-hat depth of the cross-correlation (centre minus flanks)", "linear")]]
PARAMETERS = ["fold_EG", "fold_GE", "fold_HE", "fold_HG", "h12", "h21", "h23", "h32"]    # columns of simulations.csv
RATES = ["R12", "R21", "R23", "R32", "Initiate"]                      # the rates a connection can regulate
BASE_COLUMNS = [f"base_{gene}_{rate}" for gene in "EG" for rate in RATES]
INK, MUTED, BLUE = "#0b0b0b", "#898781", "#2a78d6"
DIVERGING = LinearSegmentedColormap.from_list("blue_gray_red", ["#104281", "#2a78d6", "#f0efec", "#e34948", "#8f1f1f"])
DIVERGING.set_bad("#c3c2b7")      # no data
TYPE_COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]     # sub-type bar: neighbouring sub-types get different colours
FIT_COLORS = LinearSegmentedColormap.from_list("fit", ["#e8f1fd", "#86b6ef", "#2a78d6", "#0d366b"])   # fit bar: darker = better

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

class Curves:
    """All correlation curves in one array: values[simulation, row within a simulation, column].

    Every simulation has the same rows (observable, pair, lag), listed in `layout`. `columns` names the last axis:
    mean, se, rep1, rep2, ... `position` gives the index of a sim_id along the first axis.
    A simulation that is in more than one file (a test run and the swarm, say) is taken from the newest file."""

    def __init__(self, files):
        cores = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else os.cpu_count()
        if files[0].suffix == ".f32":                   # expanded runs: mean and SE, rows listed in curve_layout.csv
            layout = pd.read_csv(files[0].parent / "curve_layout.csv")
            keep = (layout["lag_min"].abs() <= PLOT_MAX_LAG_MIN).to_numpy()    # only the plotted lags are kept in memory
            self.layout = layout[keep].reset_index(drop=True)
            self.columns = ["mean", "se"]
            with ProcessPoolExecutor(max_workers=cores) as pool:
                parts = list(pool.map(read_binary_curve_file, files, [keep] * len(files)))
            self.store(parts)
            return
        first = pd.read_csv(files[0])
        block = int((first["sim_id"] == first["sim_id"].iloc[0]).sum())           # rows per simulation
        self.layout = first.iloc[:block][["observable", "pair", "lag_frames", "lag_min"]].reset_index(drop=True)
        self.columns = [c for c in first.columns if c in ("mean", "se") or c.startswith("rep")]
        with ProcessPoolExecutor(max_workers=cores) as pool:                       # read the files in parallel
            parts = list(pool.map(read_curve_file, files, [block] * len(files), [self.columns] * len(files)))
        self.store(parts)

    def store(self, parts):
        """Keep the curves of every simulation, from (sim_ids, array) of each file."""
        ids = np.concatenate([part[0] for part in parts])
        keep = ~pd.Index(ids).duplicated(keep="last")                              # files are ordered oldest to newest
        self.values = np.concatenate([part[1] for part in parts])[keep]
        self.position = pd.Series(np.arange(keep.sum()), index=ids[keep])

    def rows(self, observable, pair):
        """Row numbers and lags (in frames and in minutes) of one curve, within the plotted lag range."""
        layout = self.layout
        part = layout[(layout["observable"] == observable) & (layout["pair"] == pair) &
                      (layout["lag_min"].abs() <= PLOT_MAX_LAG_MIN)].sort_values("lag_min")
        return part.index.to_numpy(), part["lag_frames"].to_numpy(), part["lag_min"].to_numpy()


def read_binary_curve_file(path, keep):
    """The sim_ids of one .f32 curves file (from its name, curves_FIRST-LAST.f32) and its mean and SE of the kept rows."""
    first, last = (int(part) for part in path.stem.split("_")[1].split("-"))
    values = np.fromfile(path, dtype=np.float32).reshape(last - first + 1, len(keep), 2)
    return np.arange(first, last + 1), np.ascontiguousarray(values[:, keep, :])


def read_curve_file(path, block, columns):
    """The sim_ids of one curves file and its numbers as an array [simulation, row, column]."""
    frame = pd.read_csv(path, usecols=["sim_id"] + columns)
    return frame["sim_id"].to_numpy()[::block], frame[columns].to_numpy(dtype=np.float32).reshape(-1, block, len(columns))


def run_folder(config, folder, profile):
    """The run folder, whether the run is expanded, and how many replicates each simulation has."""
    if "profiles" in config:                            # config_expanded.toml of 1_generate_models_expanded.py
        if profile not in config["profiles"]:
            raise SystemExit(f"no profile {profile}; available: {', '.join(config['profiles'])}")
        return folder / config["output"]["run_root"] / profile, True, config["profiles"][profile]["replicates"]
    return folder / config["output"]["run_dir"], False, config["simulation"]["replicates"]


def read_simulations(run_dir, expanded):
    """The simulation list, with the columns of a regular simulations.csv plus HE_state and HG_state.
    A state is one state ("3"), several states that each regulate ("1+2"), or "R" (Rsum)."""
    if expanded:
        models = pd.read_csv(run_dir / "expanded_models.csv", dtype=str)
        points = pd.read_csv(run_dir / "expanded_simulations.csv", usecols=["sim_id", "architecture_id"] + PARAMETERS)
        table = pd.DataFrame(dict(model_id=models["architecture_id"], group=models["group"],
                                  hidden=models["group"].map({"H2": "H2", "H3_symmetric": "H3", "H3_asymmetric": "H3"}).fillna("none")))
        for edge in ("EG", "GE", "HE", "HG"):
            table[f"{edge}_state"] = models[f"{edge}_states"].mask(models[f"{edge}_source_kind"] == "reporter_count", "R")
            table[f"{edge}_rates"] = models[f"{edge}_rates"]
        return points.rename(columns={"architecture_id": "model_id"}).merge(table, on="model_id")
    table = pd.read_csv(run_dir / "simulations.csv", dtype={"EG_state": str, "GE_state": str}, low_memory=False)
    if "fold" in table:                                 # results made before the folds were split per connection
        for edge in ("EG", "GE", "HE", "HG"):
            table[f"fold_{edge}"] = table["fold"].where(table[f"{edge}_rates"].notna())
    for name in PARAMETERS:
        table[name] = table.get(name, np.nan)
    default = {"H2": ("2", "2"), "H3": ("1", "3")}      # hidden states that regulate E and G
    for i, edge in enumerate(("HE", "HG")):
        table[f"{edge}_state"] = table["hidden"].map({h: states[i] for h, states in default.items()}).where(table[f"{edge}_rates"].notna())
    return table


def load(config_path, profile):
    config = tomllib.loads(config_path.read_text())
    folder = config_path.parent
    run_dir, expanded, replicates = run_folder(config, folder, profile)
    simulations = read_simulations(run_dir, expanded)
    oldest_first = lambda pattern: sorted((run_dir / "results").glob(pattern), key=lambda f: f.stat().st_mtime)
    metrics = pd.concat([pd.read_csv(f) for f in oldest_first("metrics_*.csv")])
    repeated = int(metrics["sim_id"].duplicated().sum())
    if repeated:
        print(f"{repeated:,} simulations are in more than one results file; the newest file is used for each")
    metrics = metrics.drop_duplicates("sim_id", keep="last")
    curves = Curves(oldest_first("curves_*.f32") or oldest_first("curves_*.csv"))
    experiment = pd.read_csv(folder / config["inputs"]["experimental_curves"])
    experiment_on = pd.read_csv(folder / config["inputs"]["experimental_metrics"]).set_index("metric")["mean"]
    # One row per simulation: parameters, metrics, score, rank (1 = best) and where its curves are in `curves`.
    table = simulations.merge(metrics, on=["sim_id", "model_id"]).sort_values("sim_id").reset_index(drop=True)
    table["curve_row"] = curves.position.reindex(table["sim_id"]).to_numpy()
    hat_experiment = add_hat_metrics(table, curves, experiment)
    table["score"] = score(table)
    table["rank"] = table["score"].rank(method="first").astype("Int64")
    table = display_order(table)
    done = f'{len(table):,} of {len(simulations):,} simulations have results'
    print(done)
    sweep = f"expanded, profile {profile}" if expanded else f'sweep mode: {config["sweep"]["mode"]}'
    done += f' · {replicates} replicate datasets each · simulator: {config["simulation"]["simulator"]} · {sweep}' 
    statistic = config["inputs"]["rate_statistic"]
    fitted = read_rates(folder / config["inputs"]["rates_E"], statistic) + read_rates(folder / config["inputs"]["rates_G"], statistic)
    return config, run_dir, table, curves, experiment, experiment_on, done, fitted, hat_experiment


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
    regulated = []
    for column, gene, source, fold in [("GE_rates", "E", f'G state {sim["GE_state"]}', "fold_GE"), ("HE_rates", "E", f'H state {sim["HE_state"]}', "fold_HE"),
                                       ("EG_rates", "G", f'E state {sim["EG_state"]}', "fold_EG"), ("HG_rates", "G", f'H state {sim["HG_state"]}', "fold_HG")]:
        if pd.notna(sim[column]):
            regulated += [dict(gene=gene, rate=rate, by=source, fold=fold, index=BASE_COLUMNS.index(f"base_{gene}_{rate}"))
                          for rate in sim[column].split("+")]
    return regulated


def experiment_curve(experiment, observable, pair, lag_frames):
    """Experimental mean and SE of one curve on the given lags (NaN where the experiment has no value)."""
    part = experiment[(experiment["observable"] == observable) & (experiment["pair"] == pair)]
    part = part.set_index("lag_frames").reindex(lag_frames)
    return part["mean"].to_numpy(), part["se"].to_numpy()


def mean_curves(curves, table, rows):
    """The mean curve (over replicates) of every simulation of `table`: array [simulation, lag]."""
    return curves.values[table["curve_row"].to_numpy()[:, None], rows[None, :], 0]


def add_hat_metrics(table, curves, experiment):
    """Add the Mexican-hat metrics of the cross-correlation to `table` (see score() for the column names)
    and return the same quantities for the experiment: {observable: {center, left, right, hat}}."""
    reference = {}
    for observable in ("binary", "reporter"):
        rows, lag_frames, lags = curves.rows(observable, "EG")
        windows = dict(center=np.abs(lags) <= HAT_CENTER_MIN,
                       left=(lags >= -HAT_FLANKS_MIN[1]) & (lags <= -HAT_FLANKS_MIN[0]),
                       right=(lags >= HAT_FLANKS_MIN[0]) & (lags <= HAT_FLANKS_MIN[1]))
        simulated, target = mean_curves(curves, table, rows), experiment_curve(experiment, observable, "EG", lag_frames)[0]
        for name, window in windows.items():
            table[f"{name}_{observable}"] = np.nanmean(simulated[:, window], axis=1)
        center, left, right = (table[f"{name}_{observable}"] for name in windows)
        table[f"hat_{observable}"] = center - (left + right) / 2
        table[f"hat_shape_{observable}"] = (center > 0) & (left < 0) & (right < 0)
        reference[observable] = {name: float(np.nanmean(target[window])) for name, window in windows.items()}
        reference[observable]["hat"] = reference[observable]["center"] - (reference[observable]["left"] + reference[observable]["right"]) / 2
    return reference


BLOCKS = ["block1", "block2", "block3"]      # names of the blocks a simulation belongs to, from coarse to fine


def display_order(table):
    """Sort the simulations for the heatmaps and name the blocks they form (columns block1, block2, block3).

    Order: model group; the fold-1 control first; source state; regulated rates; hidden-unit rates; folds
    (all from smallest to largest). Blocks, coarse to fine:
        E_to_G        E1                  E1_R12
        G_to_E        G1                  G1_R12
        reciprocal    E1/G1               E1_R12              E1_R12 · G1_R21
        hidden        H2 R12              h12 0.01 h21 0.1    (H3 asymmetric: also h23, h32)
    """
    text = lambda column: table[column].fillna("")
    number = lambda column: table[column].map(lambda value: f"{value:g}", na_action="ignore").fillna("")
    group = table["group"]
    e, g = "E" + text("EG_state"), "G" + text("GE_state")
    e_rates, g_rates = e + "_" + text("EG_rates"), g + "_" + text("GE_rates")
    # Hidden models: "H2 R12" when both genes have the same regulated rates and the usual hidden states, otherwise
    # "H2 E2:R12 G2:R21" (states and rates of each gene).
    usual = ((table["hidden"] == "H2") & (text("HE_state") == "2") & (text("HG_state") == "2")) | \
            ((table["hidden"] == "H3") & (text("HE_state") == "1") & (text("HG_state") == "3"))
    same = usual & (text("HE_rates") == text("HG_rates"))
    hidden = group.map({"H2": "H2", "H3_symmetric": "H3s", "H3_asymmetric": "H3a"}).fillna("") + " " + np.where(
        same, text("HE_rates"), "E" + text("HE_state") + ":" + text("HE_rates") + " G" + text("HG_state") + ":" + text("HG_rates"))
    hidden_first = "h12 " + number("h12") + " h21 " + number("h21")
    kinds = [group == "E_to_G", group == "G_to_E", group == "reciprocal", table["hidden"] != "none"]
    table = table.assign(
        block1=np.select(kinds, [e, g, e + "/" + g, hidden], default=group),
        block2=np.select(kinds, [e_rates, g_rates, e_rates, hidden_first], default=group))
    table["block3"] = np.select([kinds[2], group == "H3_asymmetric"],
                                [e_rates + " · " + g_rates, hidden_first + " h23 " + number("h23") + " h32 " + number("h32")],
                                default=table["block2"])
    control = table["model_id"].str.endswith("gamma0").to_numpy()
    table.loc[control, BLOCKS] = "control"
    # States and rates in the order the generator lists them (the controls are left out, as they come first anyway).
    first_seen = lambda column: pd.Categorical(table[column], categories=table.loc[~control, column].dropna().unique()).codes
    keys = [first_seen("group"), ~control] + [first_seen(column) for column in
            ("EG_state", "GE_state", "EG_rates", "GE_rates", "HE_state", "HG_state", "HE_rates", "HG_rates")] + \
           [table[name].fillna(-1) for name in PARAMETERS[4:] + PARAMETERS[:4]]
    return table.iloc[np.lexsort(keys[::-1])].reset_index(drop=True)


def sorted_within_groups(table, column, ascending=True):
    """`table` with the order of the model groups kept and the simulations of each group sorted by `column`."""
    order = {group: i for i, group in enumerate(dict.fromkeys(table["group"]))}
    return (table.assign(group_order=table["group"].map(order))
            .sort_values(["group_order", column], ascending=[True, ascending]).drop(columns="group_order"))


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


def blocks_at(levels, depth):
    """Start rows and names of the blocks at one level: a block ends where any name down to that level changes."""
    change = np.zeros(len(levels[0]), dtype=bool)
    change[0] = True
    for names in levels[:depth + 1]:
        change[1:] |= names[1:] != names[:-1]
    starts = np.flatnonzero(change)
    return starts, levels[depth][starts]


def panel_marks(levels, height):
    """What to draw on a heatmap panel that is `height` points high. levels: block names of its rows, coarse to fine.
    Returns the rows above which a white line goes (the finest level whose blocks are at least 1.5 points high),
    and (positions, names, font size, level) of the finest level whose names fit at 4 points or more, or None."""
    rows, lines, names = len(levels[0]), [], None
    for depth in range(len(levels)):
        starts, labels = blocks_at(levels, depth)
        if 1 < len(starts) <= min(height / 1.5, rows / 2):
            lines = starts[1:]
        if 1 < len(starts) <= height / 5:
            names = ((starts + np.r_[starts[1:], rows]) / 2, list(labels), min(6.0, 0.8 * height / len(starts)), depth)
    return lines, names


def heatmap_figure(title, panels, x_edges, limit, colorbar_label, column_titles=("",), x_label="", x_ticks=None):
    """An A4 page of heatmap panels stacked vertically, with one or more columns.

    panels: dictionaries with label, height (relative), arrays (one [row, heatmap column] array per column) and
    optionally levels (per column, the block names of the rows from coarse to fine; see panel_marks).
    Blocks are separated by thin white lines in every column and named beside the first column."""
    figure, axes = plt.subplots(len(panels), len(column_titles), figsize=A4, squeeze=False,
                                gridspec_kw=dict(height_ratios=[panel["height"] for panel in panels], hspace=0.07,
                                                 wspace=0.18 if len(column_titles) > 1 else 0.05))
    points = 72 * A4[1] * 0.9 * 0.93 / sum(panel["height"] for panel in panels)        # height of one unit, in points
    for i, panel in enumerate(panels):
        for j, values in enumerate(panel["arrays"]):
            axis = axes[i, j]
            image = axis.pcolorfast(x_edges, np.arange(len(values) + 1), np.ma.masked_invalid(values),
                                    cmap=DIVERGING, vmin=-limit, vmax=limit)
            axis.invert_yaxis()                                     # first row on top
            axis.set_yticks([])
            axis.tick_params(labelbottom=i == len(panels) - 1, length=3, colors=MUTED, labelcolor=INK, labelsize=7)
            if x_ticks is not None:
                axis.set_xticks(*x_ticks)
            for spine in axis.spines.values():
                spine.set_visible(False)
            if panel.get("levels") is None or panel["levels"][j] is None:
                continue
            lines, names = panel_marks(panel["levels"][j], panel["height"] * points)
            axis.hlines(lines, x_edges[0], x_edges[-1], color="white", linewidth=0.5)
            if j == 0 and names is not None:
                positions, texts, size, depth = names
                axis.set_yticks(positions, texts, fontsize=size)
                axis.tick_params(axis="y", length=0, pad=2)
                if depth > 0 and len(column_titles) > 1:            # also name the coarsest blocks, further left
                    starts, coarse = blocks_at(panel["levels"][j], 0)
                    for position, text in zip((starts + np.r_[starts[1:], len(values)]) / 2, coarse):
                        axis.annotate(text, (0, position), xycoords=("axes fraction", "data"), xytext=(-78, 0),
                                      textcoords="offset points", ha="right", va="center", fontsize=6, fontweight="bold")
        axes[i, 0].set_ylabel(panel["label"], rotation=0, ha="right", va="center", fontsize=8, labelpad=48)
    for axis, text in zip(axes[0], column_titles):
        axis.set_title(text, fontsize=6.5)
    for axis in axes[-1]:
        axis.set_xlabel(x_label, fontsize=8)
    figure.suptitle(title, x=0.04, y=0.985, ha="left", fontsize=11, fontweight="bold")
    figure.subplots_adjust(left=0.3, right=0.87, top=0.93 if len(column_titles) > 1 else 0.95, bottom=0.05)
    figure.colorbar(image, cax=figure.add_axes([0.9, 0.35, 0.015, 0.3]), label=colorbar_label, extend="both")
    return figure


def overview_heatmap(title, table, values, x_edges, limit, colorbar_label, top_row=None, blocks=True, **options):
    """Every model group on one page, rows in the order of `table` (values = array [simulation, heatmap column]).
    A panel of 1 simulation is as high as the experiment strip; bigger groups are higher, but not to scale.
    blocks: separate and name the model types (only meaningful when `table` is in display order)."""
    panels = [] if top_row is None else [dict(label="experiment", height=1, arrays=[top_row[None, :]])]
    for group in dict.fromkeys(table["group"]):
        inside = (table["group"] == group).to_numpy()
        count = int(inside.sum())
        panels.append(dict(label=f"{group}\n{count:,} simulation{'s' * (count > 1)}", height=1 + 2 * np.log10(count),
                           arrays=[values[inside]],
                           levels=[[table[name].to_numpy()[inside] for name in BLOCKS]] if blocks else None))
    return heatmap_figure(title, panels, x_edges, limit, colorbar_label, **options)


def curve_heatmap(table, curves, experiment, observable, pair, title, blocks=True):
    """One curve of every simulation, all model groups, rows in the order of `table`; top strip = experiment."""
    rows, lag_frames, lags = curves.rows(observable, pair)
    target, _ = experiment_curve(experiment, observable, pair, lag_frames)
    limit = 2 * np.nanmax(np.abs(target))        # colour scale set by the experiment; stronger simulations saturate
    return overview_heatmap(title, table, mean_curves(curves, table, rows), bin_edges(lags), limit,
                            "correlation (blue negative, red positive)", top_row=target, blocks=blocks,
                            x_label="lag (min); positive = MYC after CCAT1" if pair == "EG" else "lag (min)")


def group_heatmap(table, curves, experiment, group, observable, pair, title):
    """One curve of one model group on a full page, in three columns: in parameter order (model types separated and
    named), and sorted by NRMSE and by CCC within each model type (best on top)."""
    rows, lag_frames, lags = curves.rows(observable, pair)
    target, _ = experiment_curve(experiment, observable, pair, lag_frames)
    part = table[table["group"] == group]
    part = part.assign(model_order=pd.factorize(part["model_id"])[0])            # model types in display order
    by_nrmse = part.sort_values(["model_order", f"nrmse_{observable}_{pair}"])
    by_ccc = part.sort_values(["model_order", f"ccc_{observable}_{pair}"], ascending=[True, False])
    model_blocks = lambda ordered: [ordered["model_id"].to_numpy()]
    panels = [dict(label="experiment", height=1, arrays=[target[None, :]] * 3),
              dict(label="", height=45, arrays=[mean_curves(curves, ordered, rows) for ordered in (part, by_nrmse, by_ccc)],
                   levels=[[part[name].to_numpy() for name in BLOCKS], model_blocks(by_nrmse), model_blocks(by_ccc)])]
    return heatmap_figure(title, panels, bin_edges(lags), 2 * np.nanmax(np.abs(target)), "correlation (blue negative, red positive)",
                          column_titles=["in parameter order\n(hidden rates, then fold;\nsmallest on top)",
                                         "each model type\nsorted by NRMSE\n(best on top)",
                                         "each model type\nsorted by CCC\n(best on top)"], x_label="lag (min)")


def on_fraction_heatmap(table, experiment_on):
    ratio = np.column_stack([np.log2(table["on_E"] / experiment_on["on_CCAT1"]), np.log2(table["on_G"] / experiment_on["on_MYC"])])
    limit = max(1.0, np.nanpercentile(np.abs(ratio), 99))
    title = (f'ON fraction relative to experiment (experiment: CCAT1 {experiment_on["on_CCAT1"]:.3f}, '
             f'MYC {experiment_on["on_MYC"]:.3f})')
    return overview_heatmap(title, table, ratio, np.array([0.0, 1.0, 2.0]), limit,
                            "log2(simulated / experimental); 0 = equal, +1 = twice as often ON",
                            x_ticks=([0.5, 1.5], ["CCAT1 (E)", "MYC (G)"]))


def top_curves(table, curves, experiment, experiment_on):
    """The TOP_N best simulations: every curve against experiment and the uncoupled control, one row per simulation."""
    best = table.sort_values("rank").head(TOP_N)
    control = table.loc[table["group"] == "uncoupled", "curve_row"]
    figure, axes = plt.subplots(len(best), len(CURVES), figsize=(15, 1.25 * len(best) + 1.4), sharex=True, sharey="col",
                                squeeze=False)
    for row, (_, sim) in enumerate(best.iterrows()):
        label = (f'#{sim["rank"]}  simulation {sim["sim_id"]}  score {sim["score"]:.3f}\n{sim["model_id"]}\n{parameter_text(sim)}\n'
                 f'ON: CCAT1 {sim["on_E"]:.3f} (exp. {experiment_on["on_CCAT1"]:.3f}), '
                 f'MYC {sim["on_G"]:.3f} (exp. {experiment_on["on_MYC"]:.3f})')
        axes[row, 0].set_ylabel(label, rotation=0, ha="right", va="center", fontsize=6.5)
        for column, (observable, pair, name) in enumerate(CURVES):
            axis = axes[row, column]
            rows, lag_frames, lags = curves.rows(observable, pair)
            simulated = curves.values[sim["curve_row"], rows, :]            # [lag, column]: mean, se, rep1, rep2, ...
            if len(control):
                axis.plot(lags, curves.values[control.iloc[0], rows, 0], color=MUTED, linewidth=0.9, linestyle="--",
                          label="uncoupled control")
            target, target_se = experiment_curve(experiment, observable, pair, lag_frames)
            for mean, se, color, text in [(target, target_se, INK, "experiment ± SE"),
                                          (simulated[:, 0], simulated[:, 1], BLUE, "simulation mean ± SE over replicates")]:
                axis.fill_between(lags, mean - se, mean + se, color=color, alpha=0.18, linewidth=0)
                axis.plot(lags, mean, color=color, linewidth=1.3, label=text)
            axis.axhline(0, color=MUTED, linewidth=0.5)
            axis.tick_params(labelsize=6.5, colors=MUTED, labelcolor=INK)
            for side in ("top", "right"):
                axis.spines[side].set_visible(False)
            if row == 0:
                axis.set_title(name, fontsize=7.5)
            if row == len(best) - 1:
                axis.set_xlabel("lag (min)", fontsize=7.5)
    handles, labels = axes[0, 0].get_legend_handles_labels()
    figure.legend(handles, labels, loc="upper right", ncol=3, frameon=False, fontsize=8)
    figure.suptitle(f"The {len(best)} best simulations (lowest score)", x=0.02, ha="left", fontsize=11, fontweight="bold")
    height = figure.get_figheight()
    figure.subplots_adjust(left=0.2, right=0.99, top=1 - 0.8 / height, bottom=0.5 / height, wspace=0.24, hspace=0.28)
    return figure


def best_per_architecture(table, column, ascending, top=None):
    """The best parameter set of every architecture (model_id) by `column`, optionally only the `top` best
    architectures; ordered by group and sub-type as in `table` (display order), then by fit, best first."""
    best = table.sort_values(column, ascending=ascending, na_position="last").drop_duplicates("model_id")
    if top is not None:
        best = best.head(top)
    group_rank = {group: i for i, group in enumerate(dict.fromkeys(table["group"]))}
    type_rank = {key: i for i, key in enumerate(dict.fromkeys(zip(table["group"], table["block1"])))}
    return (best.assign(group_rank=best["group"].map(group_rank), type_rank=[type_rank[key] for key in zip(best["group"], best["block1"])])
            .sort_values(["group_rank", "type_rank", column], ascending=[True, True, ascending], na_position="last")
            .drop(columns=["group_rank", "type_rank"]))


def architecture_heatmap(title, rows, curves, experiment, observable, column, ascending):
    """One row per architecture (its best parameter set), all groups on one A4 page.
    Left bar: model sub-type (neighbouring sub-types differ in colour; names where they fit). Middle: the
    cross-correlation. Right bar: the fit of the row (darker = better)."""
    lag_rows, lag_frames, lags = curves.rows(observable, "EG")
    target, _ = experiment_curve(experiment, observable, "EG", lag_frames)
    values = mean_curves(curves, rows, lag_rows)
    fit = rows[column].to_numpy(dtype=float)
    low, high = np.nanpercentile(fit, [2, 98]) if np.isfinite(fit).any() else (0, 1)
    groups = list(dict.fromkeys(rows["group"]))
    counts = [int((rows["group"] == group).sum()) for group in groups]
    figure, axes = plt.subplots(len(groups) + 1, 3, figsize=A4, squeeze=False,
                                gridspec_kw=dict(width_ratios=[0.035, 1, 0.035], wspace=0.03, hspace=0.07,
                                                 height_ratios=[1] + [1 + 2 * np.log10(count) for count in counts]))
    points = 72 * A4[1] * 0.9 * 0.93 / (1 + sum(1 + 2 * np.log10(count) for count in counts))
    x_edges, limit = bin_edges(lags), 2 * np.nanmax(np.abs(target))
    for axis in axes.flat:                              # no ticks except the lag axis and the sub-type names set below
        axis.set_yticks([])
        axis.set_xticks([])
        for spine in axis.spines.values():
            spine.set_visible(False)
    axes[0, 1].pcolorfast(x_edges, [0, 1], np.ma.masked_invalid(target[None, :]), cmap=DIVERGING, vmin=-limit, vmax=limit)
    axes[0, 0].set_ylabel("experiment", rotation=0, ha="right", va="center", fontsize=8, fontweight="bold")
    for i, (group, count) in enumerate(zip(groups, counts), start=1):
        inside = (rows["group"] == group).to_numpy()
        types = rows["block1"].to_numpy()[inside]
        starts, names = blocks_at([types], 0)
        number = np.cumsum(np.r_[True, types[1:] != types[:-1]]) - 1                   # sub-type number of every row
        colors = matplotlib.colors.ListedColormap(TYPE_COLORS)
        axes[i, 0].pcolorfast([0, 1], np.arange(count + 1), (number % len(TYPE_COLORS))[:, None], cmap=colors,
                              vmin=-0.5, vmax=len(TYPE_COLORS) - 0.5)
        image = axes[i, 1].pcolorfast(x_edges, np.arange(count + 1), np.ma.masked_invalid(values[inside]),
                                      cmap=DIVERGING, vmin=-limit, vmax=limit)
        fit_image = axes[i, 2].pcolorfast([0, 1], np.arange(count + 1), np.ma.masked_invalid(fit[inside][:, None]),
                                          cmap=FIT_COLORS if not ascending else FIT_COLORS.reversed(), vmin=low, vmax=high)
        height = (1 + 2 * np.log10(count)) * points
        if 1 < len(starts) <= height / 1.5:
            axes[i, 1].hlines(starts[1:], x_edges[0], x_edges[-1], color="white", linewidth=0.5)
        if len(starts) <= height / 5:                    # sub-type names, where they fit
            axes[i, 0].set_yticks((starts + np.r_[starts[1:], count]) / 2, list(names), fontsize=min(6.0, 0.8 * height / len(starts)))
            axes[i, 0].tick_params(axis="y", length=0, pad=2)
        axes[i, 0].set_ylabel(f"{group}\n{count:,} architecture{'s' * (count > 1)}", rotation=0, ha="right", va="center",
                              fontsize=8, labelpad=48 if len(starts) <= height / 5 else 4)
    for axis in axes.flat:
        axis.invert_yaxis()
    axes[-1, 1].set_xticks([-100, -50, 0, 50, 100])
    axes[-1, 1].set_xlim(x_edges[0], x_edges[-1])
    axes[-1, 1].tick_params(axis="x", length=3, colors=MUTED, labelcolor=INK, labelsize=7)
    axes[-1, 1].set_xlabel("lag (min); positive = MYC after CCAT1", fontsize=8)
    axes[0, 2].set_title("fit", fontsize=7)
    axes[0, 0].set_title("type", fontsize=7)
    figure.suptitle(title, x=0.04, y=0.985, ha="left", fontsize=11, fontweight="bold")
    figure.subplots_adjust(left=0.3, right=0.84, top=0.94, bottom=0.05)
    figure.colorbar(image, cax=figure.add_axes([0.9, 0.55, 0.015, 0.25]), label="correlation", extend="both")
    figure.colorbar(fit_image, cax=figure.add_axes([0.9, 0.2, 0.015, 0.25]),
                    label=f"{column.split('_')[0].upper()} of the cross-correlation ({'lower' if ascending else 'higher'} is better)")
    return figure


def shape_map(table, observable, hat_experiment):
    """Centre of the cross-correlation against the mean of its two flanks, one panel per group: every simulation of
    the group in blue, all other simulations in grey; the shaded corner is the Mexican-hat shape
    (centre positive, flanks negative); the cross is the experiment."""
    x = table[f"center_{observable}"].to_numpy(dtype=float)
    y = ((table[f"left_{observable}"] + table[f"right_{observable}"]) / 2).to_numpy(dtype=float)
    reference = hat_experiment[observable]
    ex, ey = reference["center"], (reference["left"] + reference["right"]) / 2
    finite = np.isfinite(x) & np.isfinite(y)
    xl = np.r_[np.nanpercentile(x[finite], [5, 95]), ex, 0]           # axes cover the central 90% of simulations
    yl = np.r_[np.nanpercentile(y[finite], [5, 95]), ey, 0]           # and the experiment
    pad = lambda v: (v.min() - 0.1 * np.ptp(v), v.max() + 0.1 * np.ptp(v))
    groups = [g for g in dict.fromkeys(table["group"]) if g != "uncoupled"]
    ncol = 3
    figure, axes = plt.subplots(int(np.ceil(len(groups) / ncol)), ncol, figsize=(10, 3.1 * np.ceil(len(groups) / ncol) + 0.6),
                                sharex=True, sharey=True, squeeze=False)
    rng = np.random.default_rng(0)
    background = np.flatnonzero(finite)
    background = background if len(background) <= MAX_POINTS else rng.choice(background, MAX_POINTS, replace=False)
    uncoupled = table["group"].to_numpy() == "uncoupled"
    for axis, group in zip(axes.flat, groups):
        (xmin, xmax), (ymin, ymax) = pad(xl), pad(yl)
        axis.fill_between([0, xmax], ymin, 0, color="#e1e0d9", alpha=0.6, linewidth=0)
        axis.scatter(x[background], y[background], s=3, color="#c3c2b7", alpha=0.3, linewidths=0, rasterized=True)
        mine = np.flatnonzero((table["group"].to_numpy() == group) & finite)
        mine = mine if len(mine) <= MAX_POINTS else rng.choice(mine, MAX_POINTS, replace=False)
        axis.scatter(x[mine], y[mine], s=4, color=BLUE, alpha=0.25, linewidths=0, rasterized=True)
        if uncoupled.any():
            axis.scatter(x[uncoupled], y[uncoupled], s=30, marker="D", facecolor="white", edgecolor=INK, linewidths=1, zorder=4)
        axis.scatter([ex], [ey], s=90, marker="X", color=INK, zorder=5)
        axis.axhline(0, color=MUTED, linewidth=0.6)
        axis.axvline(0, color=MUTED, linewidth=0.6)
        axis.set_xlim(xmin, xmax)
        axis.set_ylim(ymin, ymax)
        axis.set_title(f"{group} ({len(np.flatnonzero(table['group'].to_numpy() == group)):,} simulations)", fontsize=9)
        axis.tick_params(labelsize=7, colors=MUTED, labelcolor=INK)
        for side in ("top", "right"):
            axis.spines[side].set_visible(False)
    for axis in axes.flat[len(groups):]:
        axis.set_visible(False)
    for axis in axes[-1]:
        axis.set_xlabel(f"centre: mean within ±{HAT_CENTER_MIN} min", fontsize=8)
    for axis in axes[:, 0]:
        axis.set_ylabel(f"flanks: mean from {HAT_FLANKS_MIN[0]} to {HAT_FLANKS_MIN[1]} min", fontsize=8)
    figure.suptitle("Shape of the cross-correlation: X = experiment, ◇ = uncoupled control, shaded = Mexican-hat shape "
                    "(axes show the central 90% of simulations)",
                    x=0.02, ha="left", fontsize=10, fontweight="bold")
    figure.tight_layout()
    return figure


def metric_by_group(table, column, title, scale, reference=None):
    """One metric per model group: a violin, every simulation as a faint point, and the group median
    (black bar, value at the right). reference: the experiment's value, drawn as a dashed line."""
    groups = list(dict.fromkeys(table["group"]))
    forward = np.log10 if scale == "log" else (lambda values: values)
    figure, axis = plt.subplots(figsize=(8, 0.62 * len(groups) + 1.3))
    jitter = np.random.default_rng(0)
    sideways = dict(orientation="horizontal") if "orientation" in inspect.signature(axis.violinplot).parameters else dict(vert=False)
    for y, group in enumerate(groups):
        values = table.loc[table["group"] == group, column].dropna().to_numpy(dtype=float)
        values = forward(values[values > 0] if scale == "log" else values)
        if len(values) == 0:
            continue
        if len(values) >= 3 and np.ptp(values) > 0:
            body = axis.violinplot(values, positions=[y], widths=0.85, showextrema=False, **sideways)["bodies"][0]
            body.set(facecolor=BLUE, edgecolor=BLUE, alpha=0.3)
        shown = values if len(values) <= MAX_POINTS else jitter.choice(values, MAX_POINTS, replace=False)
        axis.scatter(shown, y + jitter.uniform(-0.3, 0.3, len(shown)), s=4, color=BLUE, alpha=0.08, linewidths=0, rasterized=True)
        median = np.median(values)
        axis.plot([median] * 2, [y - 0.42, y + 0.42], color=INK, linewidth=1.8)
        axis.annotate(f"{10 ** median if scale == 'log' else median:.3g}", (1.01, y), xycoords=("axes fraction", "data"),
                      va="center", fontsize=8)
    if reference is not None:
        axis.axvline(reference, color=MUTED, linestyle="--", linewidth=1)
        axis.annotate("experiment", (reference, -0.62), ha="center", fontsize=7, color=MUTED, annotation_clip=False)
    if scale == "log":                                  # the axis is in log10 units: label it with 1, 2, 5, 10, ...
        low, high = axis.get_xlim()
        ticks = [m * 10.0 ** e for e in range(int(np.floor(low)) - 1, int(np.ceil(high)) + 1) for m in (1, 2, 5)]
        ticks = [tick for tick in ticks if low <= np.log10(tick) <= high]
        axis.set_xticks(np.log10(ticks), [f"{tick:g}" for tick in ticks])
    axis.set_yticks(range(len(groups)), groups)
    axis.set_ylim(len(groups) - 0.4, -0.75)
    axis.set_xlabel(f"{title}; bar and number = group median", fontsize=8)
    axis.grid(axis="x", color="#e1e0d9", linewidth=0.6)
    axis.set_axisbelow(True)
    for side in ("top", "right", "left"):
        axis.spines[side].set_visible(False)
    figure.subplots_adjust(left=0.19, right=0.91, top=0.97, bottom=0.14)
    return figure


def metric_vs_parameters(table, column, scale):
    """One metric against each swept parameter, one row of panels per model group; the line joins the medians.
    None if no group sweeps a parameter (one parameter set per model)."""
    groups = [g for g in dict.fromkeys(table["group"]) if free_parameters(table[table["group"] == g])]
    if not groups:
        return None
    shown = [name for name in PARAMETERS if any(name in free_parameters(table[table["group"] == g]) for g in groups)]
    figure, axes = plt.subplots(len(groups), len(shown), figsize=(2.2 * len(shown) + 0.8, 1.55 * len(groups) + 0.9),
                                sharey=True, squeeze=False)
    for row, group in enumerate(groups):
        rows = table[table["group"] == group]
        free = free_parameters(rows)
        for position, name in enumerate(shown):
            axis = axes[row, position]
            if name not in free:
                axis.set_visible(False)
                continue
            drawn = rows if len(rows) <= MAX_POINTS else rows.sample(MAX_POINTS, random_state=0)
            axis.scatter(drawn[name], drawn[column], s=5, color=BLUE, alpha=0.12, linewidths=0, rasterized=True)
            if rows[name].nunique() <= 12:
                median = rows.groupby(name)[column].median()
                axis.plot(median.index, median.to_numpy(), color=INK, linewidth=1.2, marker="o", markersize=3)
            axis.set_xscale("log")
            if scale == "log":
                axis.set_yscale("log")
            axis.set_title(f"{group}: {name}", fontsize=8)
            axis.tick_params(labelsize=7, colors=MUTED, labelcolor=INK)
            for side in ("top", "right"):
                axis.spines[side].set_visible(False)
        axes[row, 0].set_ylabel(column, fontsize=7)
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
        rows.append([group, f"{len(part):,}", f'{best["score"]:.3f}', f'{part["score"].median():.3f}',
                     f'{part["ccc_reporter_EG"].median():.3f}',
                     f'{int(part["hat_shape_reporter"].sum()):,} ({100 * part["hat_shape_reporter"].mean():.0f}%)',
                     f'simulation {best["sim_id"]}: {best["model_id"]}; {parameter_text(best)}'])
    return table_html(["Group", "Simulations", "Best score", "Median score", "Median CCC (x-corr, reporter)",
                       "With Mexican-hat shape (x-corr, reporter)", "Best simulation"], rows)


def top_table(table, experiment_on, hat_experiment):
    rows = [[sim["rank"], sim["sim_id"], sim["model_id"], parameter_text(sim), f'{sim["score"]:.3f}',
             f'{sim["nrmse_binary_EG"]:.3f} / {sim["ccc_binary_EG"]:.3f}', f'{sim["nrmse_reporter_EG"]:.3f} / {sim["ccc_reporter_EG"]:.3f}',
             f'{sim["hat_reporter"]:.3f} ({"yes" if sim["hat_shape_reporter"] else "no"})',
             f'{sim["on_E"]:.3f}', f'{sim["on_G"]:.3f}'] for _, sim in table.sort_values("rank").head(TOP_N).iterrows()]
    return table_html(["Rank", "Simulation", "Model", "Parameters", "Score", "x-corr state NRMSE / CCC",
                       "x-corr reporter NRMSE / CCC", f'Hat depth, reporter (exp. {hat_experiment["reporter"]["hat"]:.3f}); shape?',
                       f'ON CCAT1 (exp. {experiment_on["on_CCAT1"]:.3f})', f'ON MYC (exp. {experiment_on["on_MYC"]:.3f})'], rows)


def json_ready(array):
    """A numeric array as nested lists for JSON, rounded, with None where there is no value."""
    array = np.round(np.asarray(array, dtype=float), 4)
    return array.tolist() if not np.isnan(array).any() else np.where(np.isnan(array), None, array.astype(object)).tolist()


def write_explorer(run_dir, table, curves, experiment, fitted, matching):
    """Write one small data file per model into <run_dir>/explorer/ (the report loads a model's file when it is
    picked) and return the data the report page itself embeds."""
    folder = run_dir / "explorer"
    folder.mkdir(exist_ok=True)
    for old in folder.glob("model_*.js"):
        old.unlink()
    curve_rows = np.array([curves.rows(observable, pair)[0] for observable, pair, _ in CURVES])      # [curve, lag]

    def simulations(rows):
        """The simulations of `rows` with the mean and SE (over replicates) of every curve: curves[curve][lag], se[curve][lag]."""
        values = curves.values[rows["curve_row"].to_numpy()[:, None, None], curve_rows[None, :, :], :]
        means, errors = json_ready(values[..., 0]), json_ready(np.nan_to_num(values[..., 1]))
        return [dict(id=int(sim["sim_id"]), par={name: sim[name] for name in PARAMETERS if sim[name] == sim[name]},
                     nrmse=json_ready([sim[f"nrmse_{observable}_{pair}"] for observable, pair, _ in CURVES]),
                     base=[float(f"{sim[column]:.4g}") for column in BASE_COLUMNS], curves=mean, se=error)
                for sim, mean, error in zip(rows.to_dict("records"), means, errors)]

    models = {}
    for number, (model, rows) in enumerate(table.groupby("model_id", sort=False), start=1):
        file = f"model_{number:04d}.js"
        (folder / file).write_text(f"explorerLoaded({json.dumps(model)},{json.dumps(simulations(rows), separators=(',', ':'))});")
        models[model] = dict(group=rows["group"].iloc[0], free=free_parameters(rows), regulated=regulated_rates(rows.iloc[0]), file=file)
    control = table[table["group"] == "uncoupled"]
    targets = [experiment_curve(experiment, observable, pair, curves.rows(observable, pair)[1]) for observable, pair, _ in CURVES]
    return dict(curve_names=[name for _, _, name in CURVES], lags=json_ready(curves.rows(*CURVES[0][:2])[2]),
                experiment=[json_ready(mean) for mean, _ in targets], experiment_se=[json_ready(se) for _, se in targets],
                models=models, control=simulations(control)[0] if len(control) else None, fitted=fitted,
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
const $ = id => document.getElementById(id);
const groups = [...new Set(Object.values(DATA.models).map(m => m.group))];
const options = (select, values, labels) => { select.innerHTML = values.map((v, i) => `<option value="${v}">${labels ? labels[i] : v}</option>`).join(""); };
const distinct = (sims, name) => [...new Set(sims.map(s => s.par[name]))].sort((a, b) => a - b);
const shade = t => { const a = [158, 197, 244], b = [13, 54, 107]; return `rgb(${a.map((v, i) => Math.round(v + (b[i] - v) * t)).join(",")})`; };

const loaded = {};                              // model name => its simulations, once its data file has arrived
let current = [];                               // the simulations of the model that is shown
window.explorerLoaded = (model, sims) => {      // called by a model's data file (explorer/model_NNNN.js)
  loaded[model] = sims;
  if (model === $("model").value) modelChanged();
};
function modelChanged() {                       // fetch the model's data if needed, then list the parameters to vary
  const model = $("model").value;
  if (!loaded[model]) {
    $("legend").textContent = "loading…";
    const script = document.createElement("script");
    script.src = "explorer/" + DATA.models[model].file;
    script.onerror = () => { $("legend").textContent = "could not load explorer/" + DATA.models[model].file + " (keep the explorer folder next to this report)"; };
    document.head.appendChild(script);
    return;
  }
  current = loaded[model];
  options($("vary"), DATA.models[model].free);
  varyChanged();
}
function varyChanged() {                        // one selector for every other free parameter (if it has few values)
  const sims = current, vary = $("vary").value;
  $("fixed").innerHTML = "";
  for (const name of DATA.models[$("model").value].free.filter(n => n !== vary)) {
    const values = distinct(sims, name);
    if (values.length > 20) continue;
    const label = document.createElement("label");
    label.innerHTML = `hold ${name} at<select data-name="${name}"></select>`;
    options(label.firstElementChild, values);
    label.firstElementChild.onchange = draw;
    $("fixed").appendChild(label);
  }
  draw();
}
function draw() {
  const curve = +$("curve").value, vary = $("vary").value;
  let sims = current.slice();
  for (const select of $("fixed").querySelectorAll("select")) sims = sims.filter(s => s.par[select.dataset.name] == select.value);
  if (vary) sims.sort((a, b) => a.par[vary] - b.par[vary]);
  const lags = DATA.lags, target = DATA.experiment[curve], se = DATA.experiment_se[curve];
  const control = DATA.control;
  const all = [...target.map((v, i) => v + (se[i] || 0)), ...target.map((v, i) => v - (se[i] || 0)), 0,
               ...sims.flatMap(s => s.curves[curve]),
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
  if (few) sims.forEach((s, i) => { svg += band(s.curves[curve], s.se[curve], shade(sims.length > 1 ? i / (sims.length - 1) : 0.6), 0.16); });
  if (control) svg += `<polyline points="${path(control.curves[curve])}" fill="none" stroke="var(--muted)" stroke-width="1.5" stroke-dasharray="5 4"/>`;
  sims.forEach((s, i) => { svg += `<polyline points="${path(s.curves[curve])}" fill="none" stroke="${shade(sims.length > 1 ? i / (sims.length - 1) : 0.6)}" stroke-width="${few ? 2 : 1}"/>`; });
  svg += `<polyline points="${path(target)}" fill="none" stroke="var(--ink)" stroke-width="2.4"/>`;
  $("chart").innerHTML = svg;
  const item = (color, text) => `<span><i style="background:${color}"></i>${text}</span>`;
  let legend = item("var(--ink)", "experiment ± SE") + (control ? item("var(--muted)", "uncoupled control") : "");
  if (sims.length > 12) legend += item(shade(0), `${sims.length} simulations, ${vary} from ${sims[0].par[vary].toPrecision(3)} (light)`) + item(shade(1), `to ${sims[sims.length - 1].par[vary].toPrecision(3)} (dark)`);
  else sims.forEach((s, i) => { legend += item(shade(sims.length > 1 ? i / (sims.length - 1) : 0.6), `${vary ? vary + " " + s.par[vary] : "control (fold 1)"} · simulation ${s.id} · NRMSE ${s.nrmse[curve] === null ? "n/a" : s.nrmse[curve].toFixed(3)}`); });
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
      const base = s.base[r.index], fold = s.par[r.fold], activity = (fitted / base - 1) / (fold - 1);
      const when = r.by.endsWith(" R") ? ` per occupied reporter position (mean ${activity.toFixed(2)} occupied)` : ` (${(100 * activity).toFixed(0)}% of the time)`;
      return `<td>${f(base)} → ${f(base * fold)}` + (fold !== 1 ? when : "") + "</td>";
    }).join("") + "</tr>";
  }
  $("rates").innerHTML = `<table>${rows}</table><p class="sub">Rates per minute. Base = the rate while its source is inactive; regulated = base × fold. For source state R every occupied reporter position of the source gene adds (fold − 1) × base.
    ${DATA.state_average ? "The base is set so that the average rate, over the time when the gene is in the state where that rate can act, equals the fitted rate; the source activity shown is over that time too."
      : "The base is set so that the average over all time equals the fitted rate."} All other rates equal their fitted values.</p>`;
}
options($("curve"), DATA.curve_names.map((_, i) => i), DATA.curve_names);
options($("group"), groups);
$("group").onchange = () => { options($("model"), Object.keys(DATA.models).filter(m => DATA.models[m].group === $("group").value)); modelChanged(); };
$("model").onchange = modelChanged; $("vary").onchange = varyChanged; $("curve").onchange = draw;
$("curve").value = 1; $("group").value = groups[Math.min(1, groups.length - 1)]; $("group").onchange();
"""


def write_report(run_dir, config, done, table, curves, experiment, experiment_on, fitted, hat_experiment, plots):
    matching = config["simulation"].get("rate_matching", "time_average")
    data = json.dumps(write_explorer(run_dir, table, curves, experiment, fitted, matching), separators=(",", ":")).replace("</", "<\\/")
    figures = ""
    for section in dict.fromkeys(plot[0] for plot in plots):                    # one heading per group of plots
        figures += f"<h2>{html.escape(section)}</h2>" + "".join(
            f'<section class="card"><figure><h3>{html.escape(title)}</h3><p class="sub">{html.escape(caption)} '
            f'<a href="plots/{stem}.pdf">PDF</a> · <a href="plots/{stem}.png">PNG</a></p>'
            f'<img src="plots/{stem}.png" alt="{html.escape(title)}" loading="lazy"></figure></section>'
            for part, stem, title, caption in plots if part == section)
    hat = hat_experiment["reporter"]
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Coupled-model sweep results</title><style>{STYLE}</style></head>
<body><main><h1>Coupled-model sweep: results</h1>
<p class="sub">{done} · rate matching: {config["simulation"].get("rate_matching", "time_average")} · score = {html.escape(SCORE_NAME)} (lower is better) · lags within ±{PLOT_MAX_LAG_MIN} min shown.</p>

<h2>Explore: how a parameter changes the curves</h2>
<section class="card"><div class="controls">
<label>curve<select id="curve"></select></label><label>model group<select id="group"></select></label>
<label>model<select id="model"></select></label><label>vary<select id="vary"></select></label><span id="fixed" class="controls" style="margin:0"></span></div>
<svg id="chart" viewBox="0 0 900 380" role="img" aria-label="selected curves against experiment"></svg><div id="legend"></div>
<div id="rates" class="scroll" style="margin-top:10px"></div>
<p class="sub">Every simulation curve is the mean of its replicate datasets; the band is ± SE over the replicates (drawn when 12 or
fewer simulations are shown). The individual replicates are in the results files (columns rep1, rep2, ...). Light to dark blue = increasing value of the varied parameter. If a model was sampled with many random
parameter sets (extensive sweep), all of them are drawn and the other parameters are not held fixed.</p></section>

<h2>Summary by model group</h2><section class="card">{group_summary(table)}</section>
<p class="sub">Mexican-hat shape: mean cross-correlation positive within ±{HAT_CENTER_MIN} min and negative in both flanks
({HAT_FLANKS_MIN[0]} to {HAT_FLANKS_MIN[1]} min on each side). Depth = centre minus the mean of the flanks. Experiment (reporter count):
centre {hat["center"]:.3f}, left {hat["left"]:.3f}, right {hat["right"]:.3f}, depth {hat["hat"]:.3f}.</p>
<h2>The {TOP_N} best simulations</h2><section class="card">{top_table(table, experiment_on, hat_experiment)}
<p class="sub">All simulations: <a href="ranked_scores.csv">ranked_scores.csv</a></p></section>
{figures}
<script>const DATA={data};</script><script>{EXPLORER_JS}</script></main></body></html>"""
    (run_dir / "results_report.html").write_text(page)

# ---------------------------------------------------------------------------

def main():
    arguments = sys.argv[1:]
    profile = "full"
    if "--profile" in arguments:
        at = arguments.index("--profile")
        profile = arguments[at + 1]
        del arguments[at:at + 2]
    if len(arguments) != 1:
        raise SystemExit(__doc__)
    config, run_dir, table, curves, experiment, experiment_on, done, fitted, hat_experiment = load(Path(arguments[0]).resolve(), profile)
    # [report] in the config: which observables and curves to draw (default: all).
    global CURVES, METRICS
    observables = config.get("report", {}).get("observables", ["reporter", "binary"])
    pairs = config.get("report", {}).get("curves", ["EG", "EE", "GG"])
    CURVES = [curve for curve in CURVES if curve[0] in observables and curve[1] in pairs]
    METRICS = [metric for metric in METRICS if metric[0].split("_")[1] in observables]
    table.drop(columns=["curve_row"] + BLOCKS).sort_values("rank").to_csv(run_dir / "ranked_scores.csv", index=False)
    folder = run_dir / "plots"
    folder.mkdir(exist_ok=True)
    for old in [*folder.glob("*.pdf"), *folder.glob("*.png")]:
        old.unlink()
    plots = []      # (report section, file name, title, caption) of every plot, in report order

    def save(section, figure, name, title, caption):
        if figure is None:
            return
        stem = f"{len(plots) + 1:02d}_{name}"
        figure.savefig(folder / f"{stem}.pdf", dpi=300)
        figure.savefig(folder / f"{stem}.png", dpi=PNG_DPI)
        plt.close(figure)
        plots.append((section, stem, title, caption))

    save("Best simulations", top_curves(table, curves, experiment, experiment_on), "top_simulations", f"The {TOP_N} best simulations",
         "Every curve of each top simulation (blue, mean ± SE) against experiment (black) and the uncoupled control (grey, dashed).")

    names = {"reporter": "reporter count", "binary": "ON/OFF state"}
    for observable in [o for o in ("reporter", "binary") if o in observables]:
        for metric, ascending in [("nrmse", True), ("ccc", False)]:
            column = f"{metric}_{observable}_EG"
            for top, section, what in [(None, "Best parameter set of every architecture", "every architecture"),
                                       (TOP_ARCHITECTURES, f"Top {TOP_ARCHITECTURES} architectures", f"the {TOP_ARCHITECTURES} best architectures")]:
                rows = best_per_architecture(table, column, ascending, top)
                title = f"Cross-correlation, {names[observable]}: {what}, best parameter set by {metric.upper()}"
                save(section, architecture_heatmap(title, rows, curves, experiment, observable, column, ascending),
                     f"architectures_{'top' if top else 'all'}_{observable}_by_{metric}", title,
                     "One row per architecture, with its best parameter set. Rows are grouped by model group and sub-type "
                     "(left bar; neighbouring sub-types differ in colour) and sorted by fit within each sub-type, best on top. "
                     "Right bar: the fit (darker = better). Panel heights are not to scale.")
        save("Shape map", shape_map(table, observable, hat_experiment), f"shape_map_{observable}",
             f"Shape of the cross-correlation, {names[observable]}",
             "Every simulation as a point: the mean of the cross-correlation near zero lag against the mean of its two flanks. "
             "Simulations in the shaded corner have the Mexican-hat shape. At most 50,000 points drawn per panel.")

    for column, title, scale in METRICS:
        reference = hat_experiment[column.split("_")[1]]["hat"] if column.startswith("hat") else None
        save("Metrics by model group", metric_by_group(table, column, title, scale, reference), f"by_group_{column}", title,
             "Violin and one faint point per simulation (at most 50,000 drawn per group); the bar and the number at the right are the group median.")
    for column, title, scale in METRICS:
        save("Metrics against each swept parameter", metric_vs_parameters(table, column, scale), f"vs_parameters_{column}", title,
             "One point per simulation (at most 50,000 drawn per panel); the black line joins the medians.")

    rows_text = ("Each row is one simulation (the mean of its replicates). Top strip: experiment. White lines separate model types "
                 "(source state and regulated rates; for hidden models also the hidden rates), named where there is room. "
                 "Within a type, rows go from the smallest to the largest parameter values. Panel heights are not to scale.")
    for observable, pair, name in CURVES:
        save("Heatmaps of all model groups, by model type and parameters", curve_heatmap(table, curves, experiment, observable, pair, name),
             f"heatmap_{pair}_{observable}", name, rows_text)
    save("Heatmaps of all model groups, by model type and parameters", on_fraction_heatmap(table, experiment_on),
         "heatmap_on_fraction", "ON fraction relative to experiment", "Rows as in the curve heatmaps.")

    for observable, name in [(o, names[o]) for o in ("reporter", "binary") if o in observables]:
        for metric, ascending in [("nrmse", True), ("ccc", False)]:
            ordered = sorted_within_groups(table, f"{metric}_{observable}_EG", ascending)
            title = f"CCAT1 x MYC cross-correlation, {name}: sorted by {metric.upper()} within each group"
            save("Cross-correlation heatmaps, each whole group sorted by fit",
                 curve_heatmap(ordered, curves, experiment, observable, "EG", title, blocks=False),
                 f"heatmap_EG_{observable}_by_{metric}", title, "Best simulations on top of each group; model types are mixed.")

    for group in dict.fromkeys(table["group"]):
        if (table["group"] == group).sum() < 2:
            continue
        for observable, pair, name in CURVES:
            title = f"{group}: {name}"
            save(f"One model group per page: {group}", group_heatmap(table, curves, experiment, group, observable, pair, title),
                 f"group_{group}_{pair}_{observable}", title,
                 "Left: by model type, then parameters; the names beside it give the model types. Middle and right: the same "
                 "model types in the same order, each sorted by NRMSE and by CCC of this curve. The columns are independent: "
                 "a row in one column is not the same simulation as the row beside it.")

    write_report(run_dir, config, done, table, curves, experiment, experiment_on, fitted, hat_experiment, plots)
    print(f"Wrote {run_dir / 'results_report.html'} (needs plots/ and explorer/ next to it) and {run_dir / 'ranked_scores.csv'}")


if __name__ == "__main__":
    main()
