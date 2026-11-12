#!/usr/bin/env python3
"""Step 1. List every coupled model and every parameter set to simulate.

    python3 1_generate_models.py config.toml

Writes two files into [output] run_dir:
    models.csv        one row per model: what regulates what, and which parameter values are tested
    simulations.csv   one row per simulation: a model plus one set of parameter values
                      (this is the work list that 3_run_simulations.jl reads)

Names used everywhere: E = CCAT1 (enhancer), G = MYC (gene), H = hidden unit.
"""
import csv, itertools, math, random, sys, tomllib
from pathlib import Path

# ---------------------------------------------------------------------------
# The models
# ---------------------------------------------------------------------------

# Rates of the target gene that a connection multiplies by the coupling fold.
RATE_SETS = ["R12", "R21", "R23", "R32", "Initiate",
             "R12+R23", "R21+R32", "R12+R23+Initiate", "R21+R32+Initiate"]

# States of the source gene that regulate, in the gene-to-gene models (E_to_G, G_to_E and reciprocal):
# promoter state 1, 2 or 3, or R = StochasticGene's "Rsum": every occupied reporter position of the source
# gene adds one coupling effect, so the rate is multiplied by 1 + (fold - 1) * (number of occupied positions).
SOURCE_STATES = ["1", "2", "3", "R"]

# The hidden groups.
#   H2: two hidden states 1 <-> 2; state 2 regulates both genes.
#   H3: three hidden states 1 <-> 2 <-> 3; state 1 regulates E, state 3 regulates G, state 2 regulates nothing.
#   "rates" are the hidden-unit rates the group sweeps; "tied" rates are set equal to another one.
HIDDEN_GROUPS = {
    "H2":            dict(hidden="H2", rates=["h12", "h21"], tied={}),
    "H3_symmetric":  dict(hidden="H3", rates=["h12", "h21"], tied={"h23": "h21", "h32": "h12"}),   # R32 = R12, R23 = R21
    "H3_asymmetric": dict(hidden="H3", rates=["h12", "h21", "h23", "h32"], tied={}),
}
# Every parameter of a simulation: the coupling fold of each kind of connection, then the hidden-unit rates.
#   fold_EG: E regulates G    fold_GE: G regulates E    fold_HE: H regulates E    fold_HG: H regulates G
# A regulated rate is multiplied by 1 + (fold - 1) for every active coupling effect: fold above 1 activates,
# fold below 1 inhibits, fold 1 does nothing.
PARAMETERS = ["fold_EG", "fold_GE", "fold_HE", "fold_HG", "h12", "h21", "h23", "h32"]


def model(group, name, EG="", GE="", HE="", HG="", EG_state="3", GE_state="3", control=False):
    """One model. EG = rates of G regulated by E, while E is in EG_state; GE = rates of E regulated by G, while G
    is in GE_state. HE and HG = rates of E and G regulated by the hidden unit."""
    return dict(model_id=f"{group}__{name}" if name else group, group=group,
                hidden=HIDDEN_GROUPS.get(group, {}).get("hidden", "none"),
                EG_state=EG_state if EG else "", EG_rates=EG, GE_state=GE_state if GE else "", GE_rates=GE,
                HE_rates=HE, HG_rates=HG, control=control)


def all_models(sweep):
    """Every model, group by group. Each group starts with a control whose coupling folds are 1 (gamma = 0)."""
    first = RATE_SETS[0]
    models = [model("uncoupled", "", control=True)]                      # no connections at all
    models.append(model("E_to_G", "gamma0", EG=first, control=True))
    models += [model("E_to_G", f"state{s}__{r}", EG=r, EG_state=s) for s in SOURCE_STATES for r in RATE_SETS]
    models.append(model("G_to_E", "gamma0", GE=first, control=True))
    models += [model("G_to_E", f"state{s}__{r}", GE=r, GE_state=s) for s in SOURCE_STATES for r in RATE_SETS]
    models.append(model("reciprocal", "gamma0", EG=first, GE=first, control=True))
    # Reciprocal: the two directions use the same source state, or every pair of source states.
    same = sweep.get("reciprocal_states", "same") == "same"
    state_pairs = [(s, t) for s in SOURCE_STATES for t in SOURCE_STATES if s == t or not same]
    models += [model("reciprocal", f"EG_state{s}_{r}__GE_state{t}_{q}", EG=r, GE=q, EG_state=s, GE_state=t)
               for s, t in state_pairs for r in RATE_SETS for q in RATE_SETS]
    for group in HIDDEN_GROUPS:                                           # the hidden unit regulates the same rates of E and G
        models.append(model(group, "gamma0", HE=first, HG=first, control=True))
        models += [model(group, r, HE=r, HG=r) for r in RATE_SETS]
    return models


# ---------------------------------------------------------------------------
# The parameter sets
# ---------------------------------------------------------------------------

def log_spaced(low, high, n):
    return [10 ** (math.log10(low) + (math.log10(high) - math.log10(low)) * i / (n - 1)) for i in range(n)]


def swept_parameters(m, sweep):
    """The parameters a model sweeps, as (free, tied).

    free: {parameter: its grid values (grid mode) or its [min, max] (extensive mode)}.
    tied: {parameter: the parameter it is set equal to}.
    - Gene-to-gene folds (fold_EG, fold_GE): the grid uses every listed fold and its inverse (1 / fold, inhibition).
      Reciprocal models use the same fold in both directions when reciprocal_states is "same", and every pair of
      folds when it is "all_pairs".
    - Hidden folds (fold_HE, fold_HG): activation only. The grid uses one fold for both connections; the extensive
      sweep samples the two independently.
    """
    grid = sweep["mode"] == "grid"
    gene_folds = sorted([1 / fold for fold in sweep["coupling_fold"]] + list(sweep["coupling_fold"]))
    extensive = extensive_ranges(sweep)
    values = lambda name, grid_values: grid_values if grid else extensive[name]
    free, tied = {}, {}
    if m["group"] == "E_to_G":
        free["fold_EG"] = values("fold_EG", gene_folds)
    elif m["group"] == "G_to_E":
        free["fold_GE"] = values("fold_GE", gene_folds)
    elif m["group"] == "reciprocal":
        free["fold_EG"] = values("fold_EG", gene_folds)
        if sweep.get("reciprocal_states", "same") == "same":
            tied["fold_GE"] = "fold_EG"
        else:
            free["fold_GE"] = values("fold_GE", gene_folds)
    elif m["group"] in HIDDEN_GROUPS:
        spec = HIDDEN_GROUPS[m["group"]]
        free["fold_HE"] = values("fold_HE", sweep["coupling_fold"])
        if grid:
            tied["fold_HG"] = "fold_HE"
        else:
            free["fold_HG"] = values("fold_HG", None)
        free.update({rate: values(rate, sweep["hidden_rate"]) for rate in spec["rates"]})
        tied.update(spec["tied"])
    return free, tied


def extensive_ranges(sweep):
    """[min, max] of every parameter in the extensive sweep: [sweep.extensive] where given, otherwise taken from
    the grid lists (gene-to-gene folds 1/largest to largest fold, hidden folds smallest to largest fold,
    hidden rates smallest to largest rate)."""
    folds, rates = sweep["coupling_fold"], sweep["hidden_rate"]
    ranges = dict(fold_EG=[1 / max(folds), max(folds)], fold_GE=[1 / max(folds), max(folds)],
                  fold_HE=[min(folds), max(folds)], fold_HG=[min(folds), max(folds)],
                  **{rate: [min(rates), max(rates)] for rate in ("h12", "h21", "h23", "h32")})
    ranges.update(sweep.get("extensive", {}))
    return ranges


def parameter_sets(m, sweep):
    """All parameter sets of one model: dictionaries of the parameters that apply to it."""
    free, tied = swept_parameters(m, sweep)
    if m["control"]:                                    # one simulation: every fold 1, a middle hidden rate
        middle = sorted(sweep["hidden_rate"])[len(sweep["hidden_rate"]) // 2]
        sets = [{name: 1.0 if name.startswith("fold") else middle for name in free}]
    elif sweep["mode"] == "grid":                       # every combination of the values
        sets = [dict(zip(free, combination)) for combination in itertools.product(*free.values())]
    else:                                               # n random combinations: each parameter takes n log-spaced
        n = sweep["n_extensive"]                        # values over its [min, max], each used once, in shuffled order
        rng = random.Random(f'{sweep["seed"]}|{m["model_id"]}')
        columns = {name: rng.sample(log_spaced(low, high, n), n) for name, (low, high) in free.items()}
        sets = [{name: column[i] for name, column in columns.items()} for i in range(n)]
    for s in sets:
        s.update({name: s[source] for name, source in tied.items()})
    return sets


def describe(m, name, sets, sweep):
    """Text for models.csv saying which values of one parameter a model is tested with."""
    tied = swept_parameters(m, sweep)[1]
    if name in tied:
        return f"= {tied[name]}"
    if name not in sets[0]:
        return ""
    values = sorted({s[name] for s in sets})
    if len(values) > 16:
        return f"{len(values)} values {values[0]:g} to {values[-1]:g}"
    return "; ".join(f"{v:g}" for v in values)


def check(sweep):
    """Stop with a clear message if the sweep settings cannot be used."""
    if sweep["mode"] not in ("grid", "extensive"):
        raise SystemExit('sweep.mode must be "grid" or "extensive"')
    if min(sweep["coupling_fold"]) <= 1:
        raise SystemExit("sweep.coupling_fold lists activation folds, all above 1; their inverses are added for gene-to-gene models")
    if sweep["mode"] == "extensive":
        ranges = extensive_ranges(sweep)
        if min(ranges["fold_HE"][0], ranges["fold_HG"][0]) < 1:
            raise SystemExit("[sweep.extensive] fold_HE and fold_HG must start at 1 or above: hidden models only activate")


# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    config_path = Path(sys.argv[1]).resolve()
    config = tomllib.loads(config_path.read_text())
    sweep = config["sweep"]
    check(sweep)
    run_dir = config_path.parent / config["output"]["run_dir"]
    if (run_dir / "results").exists():
        raise SystemExit(f"{run_dir} already holds results; choose a new [output] run_dir")
    run_dir.mkdir(parents=True, exist_ok=True)

    model_rows, simulation_rows = [], []
    for m in all_models(sweep):
        sets = parameter_sets(m, sweep)
        first_id = len(simulation_rows) + 1
        for s in sets:
            simulation_rows.append(dict(
                sim_id=len(simulation_rows) + 1, **{k: m[k] for k in m if k != "control"},
                **{name: f"{s[name]:.6g}" if name in s else "" for name in PARAMETERS}))
        model_rows.append(dict(
            **{k: m[k] for k in m if k != "control"}, control="yes" if m["control"] else "",
            **{name: describe(m, name, sets, sweep) for name in PARAMETERS},
            n_simulations=len(sets), first_sim_id=first_id, last_sim_id=len(simulation_rows)))

    for name, rows in [("models.csv", model_rows), ("simulations.csv", simulation_rows)]:
        with open(run_dir / name, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)

    print(f'Sweep mode: {sweep["mode"]}')
    print(f'{"group":<16}{"models":>8}{"simulations":>13}')
    for group in dict.fromkeys(m["group"] for m in model_rows):
        rows = [m for m in model_rows if m["group"] == group]
        print(f'{group:<16}{len(rows):>8}{sum(m["n_simulations"] for m in rows):>13}')
    print(f'{"total":<16}{len(model_rows):>8}{len(simulation_rows):>13}')
    print(f"Wrote {run_dir / 'models.csv'} and {run_dir / 'simulations.csv'}")


if __name__ == "__main__":
    main()
