#!/usr/bin/env python3
"""Enumerate expanded coupled-model architectures and sampled parameter sets.

    python3 1_generate_models_expanded.py config_expanded.toml --profile full --count-only
    python3 1_generate_models_expanded.py config_expanded.toml --profile limited
    python3 1_generate_models_expanded.py config_expanded.toml --profile full

Generation only: this script does NOT run simulations. Steps 2 to 4 read its files when given
config_expanded.toml and the same --profile (e.g. julia 3_run_simulations.jl config_expanded.toml --profile full --swarm).
"""
import argparse
import csv
import hashlib
import html
import itertools
import json
import math
import random
import sys
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    try:
        import tomli as tomllib
    except ModuleNotFoundError:
        raise SystemExit("Use Python 3.11+, or install tomli for Python 3.10.") from None

SCHEMA = "expanded-v2"         # v2: one fold per regulated rate (fold_EG_R12, ...); v1 had one fold per connection
RATES = ("R12", "R21", "R23", "R32", "Initiate")
GROUPS = ("E_to_G", "G_to_E", "reciprocal", "H2", "H3_symmetric", "H3_asymmetric")
EDGES = ("EG", "GE", "HE", "HG")
AXES = ("fold_EG", "fold_GE", "fold_HE", "fold_HG", "h12", "h21", "h23", "h32")      # names used in the config
RATE_FOLDS = tuple(f"fold_{edge}_{rate}" for edge in EDGES for rate in RATES)          # sampled: one per regulated rate
HIDDEN_FREE = {"H2": ("h12", "h21"), "H3_symmetric": ("h12", "h21"),
               "H3_asymmetric": ("h12", "h21", "h23", "h32")}
ARCH_FIELDS = ["schema_version", "architecture_index", "architecture_id", "group", "hidden_nstates", "is_control"]
ARCH_FIELDS += [f"{edge}_{field}" for edge in EDGES for field in ("source_kind", "states", "rates", "always_active")]
ARCH_FIELDS += ["behavior_note", "sampling_axes", "n_parameter_sets", "first_sim_id", "last_sim_id"]
POINT_FIELDS = ["schema_version", "sim_id", "point_id", "architecture_id", "parameter_index", *AXES, *RATE_FOLDS,
                "replicates", "point_seed"]


def integer(value, name, minimum=1):
    if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return value


def positive(value, name):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be a finite positive number")
    return float(value)


def subset_sizes(values, maximum, name):
    if not isinstance(values, list) or not values:
        raise ValueError(f"{name} must be a nonempty list")
    if any(integer(v, name) > maximum for v in values) or len(set(values)) != len(values):
        raise ValueError(f"{name}: use distinct subset sizes from 1 to {maximum}")
    return sorted(values)


def subsets(items, sizes):
    return [combo for size in sizes for combo in itertools.combinations(items, size)]


def read_settings(path, profile_name):
    raw = path.read_bytes()
    config = tomllib.loads(raw.decode("utf-8"))
    design = config["design"]
    chosen = design.get("groups", list(GROUPS))
    if not isinstance(chosen, list) or not chosen or len(set(chosen)) != len(chosen) or set(chosen) - set(GROUPS):
        raise ValueError(f"design.groups must be a nonempty subset of {GROUPS}")
    groups = [g for g in GROUPS if g in chosen]
    sizes = subset_sizes(design.get("rate_subset_sizes", [1, 2, 3, 4]), 5, "rate_subset_sizes")
    state_sizes = subset_sizes(design.get("promoter_subset_sizes", [1, 2, 3]), 3, "promoter_subset_sizes")
    rate_sets = subsets(RATES, sizes)
    include_r = design.get("include_reporter_source", True)
    if not isinstance(include_r, bool):
        raise ValueError("include_reporter_source must be true or false")
    sources = [("promoter_any", s) for s in subsets((1, 2, 3), state_sizes)]
    if include_r:
        sources.append(("reporter_count", ()))
    hidden_mode = design.get("hidden_sources", "fixed")
    if hidden_mode not in ("fixed", "all_nonempty"):
        raise ValueError('hidden_sources must be "fixed" or "all_nonempty"')
    profiles = config["profiles"]
    if profile_name not in profiles:
        raise ValueError(f"No profile {profile_name!r}; available: {', '.join(profiles)}")
    profile = profiles[profile_name]
    npoints = integer(profile["parameter_sets"], "profile.parameter_sets")
    replicates = integer(profile["replicates"], "profile.replicates")
    levels = profile.get("values", {})
    if not isinstance(levels, dict) or set(levels) - set(AXES):
        raise ValueError(f"profile.values can contain only {AXES}")
    for key, value in levels.items():
        positive(value, f"profile.values.{key}")
    if levels and npoints != 1:
        raise ValueError("Explicit profile.values are supported only when parameter_sets = 1")
    sampling = config["sampling"]
    integer(sampling["seed"], "sampling.seed", 0)
    for key in ("coupling_fold", "hidden_rate"):
        checked_range(sampling[key], f"sampling.{key}")
    overrides = sampling.get("ranges", {})
    if not isinstance(overrides, dict) or set(overrides) - set(AXES):
        raise ValueError(f"sampling.ranges can contain only {AXES}")
    for key, values in overrides.items():
        checked_range(values, f"sampling.ranges.{key}")
    positive(sampling.get("control_hidden_rate", 1.0), "control_hidden_rate")
    sim = config["simulation"]
    for key in ("ntraces", "nframes"):
        integer(sim[key], f"simulation.{key}")
    positive(sim["interval_min"], "simulation.interval_min")
    integer(sim["seed"], "simulation.seed", 0)
    if sim.get("rate_matching") not in ("state_average", "time_average"):
        raise ValueError("simulation.rate_matching must be state_average or time_average")
    if config["output"]["run_root"] == "":
        raise ValueError("output.run_root cannot be empty")
    per_job = integer(config.get("swarm", {}).get("simulations_per_job", 96), "simulations_per_job")
    return dict(config=config, raw=raw, profile=profile, profile_name=profile_name, groups=groups,
                rate_sets=rate_sets, sources=sources, hidden_mode=hidden_mode, npoints=npoints,
                replicates=replicates, per_job=per_job, folder=path.parent)


def checked_range(values, name):
    if not isinstance(values, list) or len(values) < 2:
        raise ValueError(f"{name} must contain at least two positive range endpoints")
    checked = [positive(value, name) for value in values]
    return min(checked), max(checked)


def hidden_pairs(group, settings):
    if settings["hidden_mode"] == "fixed":
        return [((2,), (2,))] if group == "H2" else [((1,), (3,))]
    nstates = 2 if group == "H2" else 3
    choices = subsets(tuple(range(1, nstates + 1)), range(1, nstates + 1))
    return list(itertools.product(choices, repeat=2))


def model_count(group, settings):
    rate_count, source_count = len(settings["rate_sets"]), len(settings["sources"])
    if group in ("E_to_G", "G_to_E"):
        return rate_count * source_count
    if group == "reciprocal":
        return (rate_count * source_count) ** 2
    return rate_count ** 2 * len(hidden_pairs(group, settings))


def plan_counts(settings):
    rows = [dict(group="uncoupled", noncontrol_architectures=0, controls=1, parameter_sets=1)]
    for group in settings["groups"]:
        count = model_count(group, settings)
        rows.append(dict(group=group, noncontrol_architectures=count, controls=1,
                         parameter_sets=count * settings["npoints"] + 1))
    for row in rows:
        row["architecture_entries"] = row["noncontrol_architectures"] + row["controls"]
        row["independent_datasets"] = row["parameter_sets"] * settings["replicates"]
    totals = {key: sum(row[key] for row in rows) for key in rows[0] if key != "group"}
    totals["swarm_batches_if_supported"] = math.ceil(totals["parameter_sets"] / settings["per_job"])
    totals["paired_traces"] = totals["independent_datasets"] * settings["config"]["simulation"]["ntraces"]
    totals["paired_frames"] = totals["paired_traces"] * settings["config"]["simulation"]["nframes"]
    return dict(schema_version=SCHEMA, profile=settings["profile_name"], rate_subsets=len(settings["rate_sets"]),
                direct_source_options=len(settings["sources"]), parameter_sets_per_noncontrol=settings["npoints"],
                independent_datasets_per_set=settings["replicates"], groups=rows, totals=totals)


def edge(kind, states, rates):
    return dict(kind=kind, states=states, rates=rates)


def architecture(group, connections=None, control=False):
    connections = connections or {}
    nh = 2 if group == "H2" else 3 if group.startswith("H3") else 1
    labels = []
    for name in EDGES:
        if name not in connections:
            continue
        value = connections[name]
        source = "R" if value["kind"] == "reporter_count" else "s" + "-".join(map(str, value["states"]))
        labels.append(name + "_" + source + "_" + "+".join(value["rates"]))
    model_id = "uncoupled" if group == "uncoupled" else group + "__" + ("gamma0" if control else "__".join(labels))
    return dict(group=group, architecture_id=model_id, hidden_nstates=nh, connections=connections, is_control=control)


def all_architectures(settings):
    yield architecture("uncoupled", control=True)
    targets = settings["rate_sets"]
    options = [(kind, states, rates) for kind, states in settings["sources"] for rates in targets]
    for group in settings["groups"]:
        if group in ("E_to_G", "G_to_E", "reciprocal"):
            names = ("EG", "GE") if group == "reciprocal" else ("EG",) if group == "E_to_G" else ("GE",)
            yield architecture(group, {name: edge("promoter_any", (3,), ("R12",)) for name in names}, True)
            choices = itertools.product(options, repeat=2) if len(names) == 2 else ((option,) for option in options)
            for selected in choices:
                yield architecture(group, {name: edge(*option) for name, option in zip(names, selected)})
        else:
            fixed_states = ((2,), (2,)) if group == "H2" else ((1,), (3,))
            yield architecture(group, {name: edge("hidden_any", states, ("R12",))
                                       for name, states in zip(("HE", "HG"), fixed_states)}, True)
            for source_e, source_g in hidden_pairs(group, settings):
                for rates_e, rates_g in itertools.product(targets, repeat=2):
                    yield architecture(group, {"HE": edge("hidden_any", source_e, rates_e),
                                               "HG": edge("hidden_any", source_g, rates_g)})


def axes_for(model):
    """Sampled axes of a model: one fold per regulated rate of every connection, then the hidden rates."""
    return tuple(f"fold_{name}_{rate}" for name in EDGES if name in model["connections"]
                 for rate in model["connections"][name]["rates"]) + HIDDEN_FREE.get(model["group"], ())


def config_name(name):
    """The config name of an axis: fold_EG_R12 takes its range and limited-profile value from fold_EG."""
    return "_".join(name.split("_")[:2]) if name.startswith("fold_") else name


def axis_range(name, settings):
    sampling = settings["config"]["sampling"]
    key = "coupling_fold" if name.startswith("fold_") else "hidden_rate"
    return checked_range(sampling.get("ranges", {}).get(config_name(name), sampling[key]), name)


def log_levels(low, high, count):
    if count == 1:
        return [math.exp((math.log(low) + math.log(high)) / 2)]
    values = [math.exp(math.log(low) + (math.log(high) - math.log(low)) * i / (count - 1)) for i in range(count)]
    values[0], values[-1] = low, high
    return values


def parameter_sets(model, settings):
    axes = axes_for(model)
    count = 1 if model["is_control"] else settings["npoints"]
    columns = {}
    for name in axes:
        if model["is_control"]:
            columns[name] = [1.0 if name.startswith("fold_") else settings["config"]["sampling"].get("control_hidden_rate", 1.0)]
            continue
        low, high = axis_range(name, settings)
        if config_name(name) in settings["profile"].get("values", {}):
            value = float(settings["profile"]["values"][config_name(name)])
            if not low <= value <= high:
                raise ValueError(f"Limited-profile {name}={value} lies outside [{low}, {high}]")
            columns[name] = [value]
        else:
            values = log_levels(low, high, count)
            seed_text = f'{settings["config"]["sampling"]["seed"]}|{model["architecture_id"]}|{name}'
            seed = int.from_bytes(hashlib.sha256(seed_text.encode()).digest(), "big")
            random.Random(seed).shuffle(values)
            columns[name] = values
    for index in range(count):
        point = {name: column[index] for name, column in columns.items()}
        if model["group"] == "H3_symmetric":
            point.update(h23=point["h21"], h32=point["h12"])
        yield point


def architecture_row(model, index, first_id, count):
    row = dict(schema_version=SCHEMA, architecture_index=index, architecture_id=model["architecture_id"],
               group=model["group"], hidden_nstates=model["hidden_nstates"], is_control=str(model["is_control"]).lower(),
               sampling_axes=";".join(axes_for(model)), n_parameter_sets=count, first_sim_id=first_id,
               last_sim_id=first_id + count - 1)
    constant = []
    for name in EDGES:
        connection = model["connections"].get(name)
        if connection is None:
            continue
        nstates = model["hidden_nstates"] if name.startswith("H") else 3
        always = connection["kind"] != "reporter_count" and set(connection["states"]) == set(range(1, nstates + 1))
        row.update({f"{name}_source_kind": connection["kind"], f"{name}_states": "+".join(map(str, connection["states"])),
                    f"{name}_rates": "+".join(connection["rates"]), f"{name}_always_active": str(always).lower()})
        if always:
            constant.append(name)
    row["behavior_note"] = ("uncoupled control" if model["is_control"] else
                            "; ".join(f"{name}: constant multiplier is absorbed by basal-rate matching" for name in constant))
    return row


def print_counts(plan):
    print(f'Profile {plan["profile"]}: {plan["rate_subsets"]} rate subsets, {plan["direct_source_options"]} direct-source options')
    print(f'{"group":<17}{"architectures":>16}{"controls":>11}{"parameter sets":>18}{"datasets":>18}')
    for row in plan["groups"]:
        print(f'{row["group"]:<17}{row["noncontrol_architectures"]:>16,}{row["controls"]:>11,}'
              f'{row["parameter_sets"]:>18,}{row["independent_datasets"]:>18,}')
    total = plan["totals"]
    print(f'{"TOTAL":<17}{total["noncontrol_architectures"]:>16,}{total["controls"]:>11,}'
          f'{total["parameter_sets"]:>18,}{total["independent_datasets"]:>18,}')
    print("Next: 2_model_report.py, 3_run_simulations.jl and 4_plot_results.py with this config and --profile " + plan["profile"])


def report_html(plan, settings):
    rows = "".join("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in (
        r["group"], f'{r["noncontrol_architectures"]:,}', r["controls"], f'{r["parameter_sets"]:,}',
        f'{r["independent_datasets"]:,}')) + "</tr>" for r in plan["groups"])
    total = plan["totals"]
    rate_choices = ", ".join("+".join(r) for r in settings["rate_sets"])
    return f'''<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Expanded coupled-model design</title><style>body{{max-width:1100px;margin:40px auto;padding:0 20px;font:16px/1.6 system-ui;color:#163047}}
table{{border-collapse:collapse;width:100%}}td,th{{text-align:left;padding:8px;border-bottom:1px solid #ccd7df}}code{{background:#eef3f7}}</style>
<h1>Expanded coupled-model design</h1><p>Profile: <b>{html.escape(plan['profile'])}</b>. E = CCAT1; G = MYC.</p>
<p>{total['noncontrol_architectures']:,} noncontrol architectures + {total['controls']} controls;
{total['parameter_sets']:,} parameter sets; {total['independent_datasets']:,} independent datasets.</p>
<table><tr><th>Family</th><th>Architectures</th><th>Controls</th><th>Parameter sets</th><th>Datasets</th></tr>{rows}</table>
<p>Each noncontrol architecture has {settings['npoints']} parameter sets, each with {settings['replicates']} independent datasets.
Selected promoter states are an OR condition and share one rate subset. Every regulated rate has its own fold
(fold_EG_R12, fold_EG_R23, ...), sampled independently over the range of its connection.</p>
<p>When included, Rsum is a separate observed-source option: each occupied reporter position adds (fold − 1),
giving multiplier 1 + (fold − 1) × reporter_count. It is not a binary reporter-ON indicator.</p>
<p>Rate subsets ({len(settings['rate_sets'])}): {html.escape(rate_choices)}</p>
<p>H2: state 2 normally regulates E and G. H3: state 1 normally regulates E, state 3 regulates G, and state 2 is inactive.
Hidden-source selection mode: {html.escape(settings['hidden_mode'])}. Symmetric H3 ties h23=h21 and h32=h12;
asymmetric H3 samples all four independently.</p>
<p>Every free parameter gets the requested marginal log-spaced levels, independently shuffled into the parameter sets.
This sampling does not enumerate every numerical parameter combination. All-state sources are retained and flagged as constant:
their fold is absorbed by average-rate matching.</p><p>Simulate with 3_run_simulations.jl and report with 2_model_report.py and
4_plot_results.py, giving them this config and --profile {html.escape(plan['profile'])}.</p></html>'''


def write_manifests(settings, plan):
    config = settings["config"]
    folder = settings["folder"] / config["output"]["run_root"] / settings["profile_name"]
    folder.mkdir(parents=True, exist_ok=False)
    (folder / "INCOMPLETE").write_text("Generation has not completed. Do not consume these manifests.\n")
    (folder / "config_snapshot.toml").write_bytes(settings["raw"])
    (folder / "generator_snapshot.py").write_bytes(Path(__file__).read_bytes())
    fingerprint = hashlib.sha256(SCHEMA.encode() + settings["raw"] + settings["profile_name"].encode() + Path(__file__).read_bytes()).hexdigest()
    architecture_path, point_path = folder / "expanded_models.csv.part", folder / "expanded_simulations.csv.part"
    sim_id, architecture_index = 0, 0
    with architecture_path.open("w", newline="") as model_file, point_path.open("w", newline="") as point_file:
        models = csv.DictWriter(model_file, fieldnames=ARCH_FIELDS)
        points = csv.DictWriter(point_file, fieldnames=POINT_FIELDS)
        models.writeheader()
        points.writeheader()
        for architecture_index, model in enumerate(all_architectures(settings), 1):
            count = 1 if model["is_control"] else settings["npoints"]
            models.writerow(architecture_row(model, architecture_index, sim_id + 1, count))
            for parameter_index, values in enumerate(parameter_sets(model, settings), 1):
                sim_id += 1
                payload = json.dumps(values, sort_keys=True, separators=(",", ":"), allow_nan=False)
                key = f'{fingerprint}|{model["architecture_id"]}|{parameter_index}|{payload}'
                point_id = hashlib.sha256(key.encode()).hexdigest()[:32]
                seed_key = f'{config["simulation"]["seed"]}|{point_id}'
                point_seed = int.from_bytes(hashlib.sha256(seed_key.encode()).digest()[:8], "big") & ((1 << 63) - 1)
                points.writerow(dict(schema_version=SCHEMA, sim_id=sim_id, point_id=point_id,
                                     architecture_id=model["architecture_id"], parameter_index=parameter_index,
                                     **{name: format(value, ".17g") for name, value in values.items()},
                                     replicates=settings["replicates"], point_seed=point_seed))
    if architecture_index != plan["totals"]["architecture_entries"] or sim_id != plan["totals"]["parameter_sets"]:
        raise RuntimeError("Generated counts disagree with the combinatorial plan")
    architecture_path.rename(folder / "expanded_models.csv")
    point_path.rename(folder / "expanded_simulations.csv")
    (folder / "design_summary.html").write_text(report_html(plan, settings), encoding="utf-8")
    plan.update(run_fingerprint=fingerprint, requires_expanded_runner=True,
                rate_matching=config["simulation"]["rate_matching"], inputs=config.get("inputs", {}),
                input_base_directory=str(settings["folder"]), note="Inputs are declared, not read or validated by this generator.")
    (folder / "design_summary.json").write_text(json.dumps(plan, indent=2, allow_nan=False) + "\n")
    # Publishing this marker is the final step. There is deliberately no overwrite or resume mode.
    (folder / "COMPLETE.json").write_text(json.dumps(dict(schema_version=SCHEMA, run_fingerprint=fingerprint,
                                                         architectures=architecture_index, parameter_sets=sim_id)) + "\n")
    (folder / "INCOMPLETE").unlink()
    print(f"Wrote {folder}")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("config", type=Path)
    parser.add_argument("--profile", choices=("limited", "full"), default="full")
    parser.add_argument("--count-only", action="store_true", help="Print exact combinatorial counts without creating files")
    args = parser.parse_args()
    settings = read_settings(args.config.resolve(), args.profile)
    plan = plan_counts(settings)
    print_counts(plan)
    if not args.count_only:
        write_manifests(settings, plan)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, FileExistsError) as error:
        sys.exit(f"Generation stopped: {error}")
