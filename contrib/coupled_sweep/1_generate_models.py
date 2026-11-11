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

# States of the source gene that regulate, in the one-way models (E_to_G and G_to_E):
# promoter state 1, 2 or 3, or R = StochasticGene's "Rsum": every occupied reporter position of the source
# gene adds one coupling effect, so the rate is multiplied by 1 + (fold - 1) * (number of occupied positions).
SOURCE_STATES = ["1", "2", "3", "R"]

# Hidden-unit rates that each hidden group sweeps ("free") and the rates set equal to them ("tied").
#   H2: two hidden states 1 <-> 2; state 2 regulates both genes.
#   H3: three hidden states 1 <-> 2 <-> 3; state 1 regulates E, state 3 regulates G, state 2 regulates nothing.
HIDDEN_GROUPS = {
    "H2":            dict(hidden="H2", free=["h12", "h21"], tied={}),
    "H3_symmetric":  dict(hidden="H3", free=["h12", "h21"], tied={"h23": "h21", "h32": "h12"}),   # R32 = R12, R23 = R21
    "H3_asymmetric": dict(hidden="H3", free=["h12", "h21", "h23", "h32"], tied={}),
}
HIDDEN_RATES = ["h12", "h21", "h23", "h32"]


def model(group, name, EG="", GE="", HE="", HG="", state="3", control=False):
    """One model. EG = rates of G regulated by E; GE = rates of E regulated by G; `state` = the state of the
    regulating gene (EG_state for E, GE_state for G). HE and HG = rates of E and G regulated by the hidden unit."""
    return dict(model_id=f"{group}__{name}" if name else group, group=group,
                hidden=HIDDEN_GROUPS.get(group, {}).get("hidden", "none"),
                EG_state=state if EG else "", EG_rates=EG, GE_state=state if GE else "", GE_rates=GE,
                HE_rates=HE, HG_rates=HG, control=control)


def all_models():
    """Every model, group by group. Each group starts with a control whose coupling fold is 1 (gamma = 0)."""
    first = RATE_SETS[0]
    models = [model("uncoupled", "", control=True)]                      # no connections at all
    models.append(model("E_to_G", "gamma0", EG=first, control=True))
    models += [model("E_to_G", f"state{s}__{r}", EG=r, state=s) for s in SOURCE_STATES for r in RATE_SETS]
    models.append(model("G_to_E", "gamma0", GE=first, control=True))
    models += [model("G_to_E", f"state{s}__{r}", GE=r, state=s) for s in SOURCE_STATES for r in RATE_SETS]
    models.append(model("reciprocal", "gamma0", EG=first, GE=first, control=True))
    models += [model("reciprocal", f"EG_{r}__GE_{q}", EG=r, GE=q) for r in RATE_SETS for q in RATE_SETS]   # both from state 3
    for group in HIDDEN_GROUPS:                                           # the hidden unit regulates the same rates of E and G
        models.append(model(group, "gamma0", HE=first, HG=first, control=True))
        models += [model(group, r, HE=r, HG=r) for r in RATE_SETS]
    return models


# ---------------------------------------------------------------------------
# The parameter sets
# ---------------------------------------------------------------------------

def log_spaced(values, n):
    lo, hi = math.log10(min(values)), math.log10(max(values))
    return [10 ** (lo + (hi - lo) * i / (n - 1)) for i in range(n)]


def parameter_sets(m, sweep):
    """All parameter sets of one model, as dictionaries with fold and the four hidden rates."""
    spec = HIDDEN_GROUPS.get(m["group"], dict(free=[], tied={}))
    values = {"fold": sweep["coupling_fold"], **{h: sweep["hidden_rate"] for h in spec["free"]}}
    if m["control"]:                                    # one simulation: fold 1, a middle hidden rate
        middle = sorted(sweep["hidden_rate"])[len(sweep["hidden_rate"]) // 2]
        sets = [{"fold": 1.0, **{h: middle for h in spec["free"]}}]
    elif sweep["mode"] == "grid":                       # every combination
        sets = [dict(zip(values, combo)) for combo in itertools.product(*values.values())]
    elif sweep["mode"] == "extensive":                  # n random combinations: each parameter takes n
        n = sweep["n_extensive"]                        # log-spaced values, each used once, in shuffled order
        rng = random.Random(f'{sweep["seed"]}|{m["model_id"]}')
        columns = {name: rng.sample(log_spaced(v, n), n) for name, v in values.items()}
        sets = [{name: column[i] for name, column in columns.items()} for i in range(n)]
    else:
        raise SystemExit('sweep.mode must be "grid" or "extensive"')
    for s in sets:
        s.update({tied: s[source] for tied, source in spec["tied"].items()})
    return sets


def describe(m, name, sets, sweep):
    """Text for models.csv saying which values of one parameter a model is tested with."""
    spec = HIDDEN_GROUPS.get(m["group"], dict(free=[], tied={}))
    if name in spec["tied"]:
        return f'= {spec["tied"][name]}'
    if name not in sets[0]:
        return ""
    values = sorted({s[name] for s in sets})
    if len(values) > 12:
        return f"{len(values)} values {values[0]:g} to {values[-1]:g} (log-spaced)"
    return "; ".join(f"{v:g}" for v in values)


# ---------------------------------------------------------------------------

def main():
    if len(sys.argv) != 2:
        raise SystemExit(__doc__)
    config_path = Path(sys.argv[1]).resolve()
    config = tomllib.loads(config_path.read_text())
    sweep = config["sweep"]
    run_dir = config_path.parent / config["output"]["run_dir"]
    if (run_dir / "results").exists():
        raise SystemExit(f"{run_dir} already holds results; choose a new [output] run_dir")
    run_dir.mkdir(parents=True, exist_ok=True)

    model_rows, simulation_rows = [], []
    for m in all_models():
        sets = parameter_sets(m, sweep)
        first_id = len(simulation_rows) + 1
        for s in sets:
            simulation_rows.append(dict(
                sim_id=len(simulation_rows) + 1, **{k: m[k] for k in m if k != "control"},
                **{name: f"{s[name]:.6g}" if name in s else "" for name in ["fold"] + HIDDEN_RATES}))
        model_rows.append(dict(
            **{k: m[k] for k in m if k != "control"},
            **{name: describe(m, name, sets, sweep) for name in ["fold"] + HIDDEN_RATES},
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
