#!/usr/bin/env python3
"""Step 2. Draw what step 1 generated, so you can check it before simulating.

    python3 2_model_report.py config.toml
    python3 2_model_report.py config_expanded.toml --profile full      (expanded runs; full is the default)

Reads <run_dir>/models.csv and writes <run_dir>/models_report.html: for every model group a
schematic, the rates that are regulated, and the parameter values that are tested.
"""
import csv, html, sys, tomllib
from pathlib import Path

PARAMETERS = {"fold_EG": "fold of E regulating G", "fold_GE": "fold of G regulating E",
              "fold_HE": "fold of H regulating E", "fold_HG": "fold of H regulating G",
              "h12": "hidden rate h12 (1 → 2), per min", "h21": "hidden rate h21 (2 → 1), per min",
              "h23": "hidden rate h23 (2 → 3), per min", "h32": "hidden rate h32 (3 → 2), per min"}

GROUP_TEXT = {
    "uncoupled": "Two independent genes. The control for every other group.",
    "E_to_G": "CCAT1 (E) regulates MYC (G) from one source state: promoter state 1, 2 or 3 multiplies the chosen rates "
              "of G by the fold; R (StochasticGene's Rsum) adds one coupling effect per occupied reporter position.",
    "G_to_E": "MYC (G) regulates CCAT1 (E) from one source state: promoter state 1, 2 or 3 multiplies the chosen rates "
              "of E by the fold; R (StochasticGene's Rsum) adds one coupling effect per occupied reporter position.",
    "reciprocal": "Both directions at once. Each direction has its own choice of rates and regulates from source state "
                  "1, 2, 3 or R. The two directions use the same source state and fold, or every pair of them "
                  "(setting reciprocal_states); the parameters below show which.",
    "H2": "Two-state hidden unit. While it is in state 2, the chosen rates of both genes are multiplied by the fold "
          "(activation only).",
    "H3_symmetric": "Three-state hidden unit, mirror-symmetric rates (h23 = h21, h32 = h12). "
                    "State 1 multiplies the chosen rates of E, state 3 those of G; state 2 regulates nothing.",
    "H3_asymmetric": "Three-state hidden unit with four independent rates. "
                     "State 1 multiplies the chosen rates of E, state 3 those of G; state 2 regulates nothing.",
}

# ---------------------------------------------------------------------------
# Schematics (inline SVG)
# ---------------------------------------------------------------------------
E_BOX, G_BOX, BOX_Y = 10, 230, 118          # left edge of the E and G gene boxes, and their top edge


def gene(x, label):
    """A gene: promoter states 1 <-> 2 <-> 3, and RNA made from state 3."""
    parts = [f'<rect class="box" x="{x}" y="{BOX_Y}" width="150" height="62" rx="8"/>',
             f'<text class="label" x="{x + 10}" y="{BOX_Y + 16}">{label}</text>']
    for i in range(3):
        cx = x + 24 + 38 * i
        parts.append(f'<circle class="state" cx="{cx}" cy="{BOX_Y + 40}" r="11"/>'
                     f'<text class="num" x="{cx}" y="{BOX_Y + 44}">{i + 1}</text>')
        if i < 2:
            parts.append(f'<line class="thin" x1="{cx + 13}" y1="{BOX_Y + 37}" x2="{cx + 25}" y2="{BOX_Y + 37}" marker-end="url(#tip)"/>'
                         f'<line class="thin" x1="{cx + 25}" y1="{BOX_Y + 43}" x2="{cx + 13}" y2="{BOX_Y + 43}" marker-end="url(#tip)"/>')
    parts.append(f'<line class="thin" x1="{x + 113}" y1="{BOX_Y + 40}" x2="{x + 124}" y2="{BOX_Y + 40}" marker-end="url(#tip)"/>'
                 f'<text class="small" x="{x + 127}" y="{BOX_Y + 44}">RNA</text>')
    return "".join(parts)


def hidden_unit(n, labels):
    """Hidden states in a row above the genes; `labels` are the rate labels above and below each pair."""
    xs = [150, 240] if n == 2 else [105, 195, 285]
    parts = [f'<text class="label" x="{xs[0] - 12}" y="18">hidden unit H</text>']
    for i, cx in enumerate(xs):
        parts.append(f'<circle class="state" cx="{cx}" cy="48" r="12"/><text class="num" x="{cx}" y="52">{i + 1}</text>')
    for i, (above, below) in enumerate(labels):
        a, b = xs[i] + 14, xs[i + 1] - 16
        parts.append(f'<line class="thin" x1="{a}" y1="45" x2="{b}" y2="45" marker-end="url(#tip)"/>'
                     f'<line class="thin" x1="{b + 2}" y1="51" x2="{a + 2}" y2="51" marker-end="url(#tip)"/>'
                     f'<text class="small mid" x="{(a + b) / 2 + 1}" y="36">{above}</text>'
                     f'<text class="small mid" x="{(a + b) / 2 + 1}" y="68">{below}</text>')
    return "".join(parts), xs


def coupling(path, label_x, label_y):
    """A regulation arrow with its "× fold" label."""
    return (f'<path class="couple" d="{path}" marker-end="url(#ctip)"/>'
            f'<text class="fold" x="{label_x}" y="{label_y}">× fold</text>')


def schematic(group):
    parts = [gene(E_BOX, "E · CCAT1"), gene(G_BOX, "G · MYC")]
    # A gene regulates from one of its states (1, 2, 3 or R), so its arrow leaves the gene box, not one state.
    if group in ("E_to_G", "reciprocal"):       # from E, over the top, onto G
        parts.append(coupling(f"M{E_BOX + 75},{BOX_Y} C{E_BOX + 75},40 {G_BOX + 60},40 {G_BOX + 60},{BOX_Y - 4}", 178, 50))
    if group in ("G_to_E", "reciprocal"):       # from G, underneath, onto E
        parts.append(coupling(f"M{G_BOX + 75},{BOX_Y + 62} C{G_BOX + 75},245 {E_BOX + 60},245 {E_BOX + 60},{BOX_Y + 66}", 178, 244))
    if group.startswith("H"):
        tied = group.endswith("_symmetric")
        if group.startswith("H2"):
            drawing, xs = hidden_unit(2, [("h12", "h21")])
            sources = (xs[1], xs[1])            # state 2 regulates both genes
        else:
            drawing, xs = hidden_unit(3, [("h12", "h21"), ("h23 = h21" if tied else "h23", "h32 = h12" if tied else "h32")])
            sources = (xs[0], xs[2])            # state 1 regulates E, state 3 regulates G
        parts.append(drawing)
        e_top, g_top = E_BOX + 82, G_BOX + 68        # where the arrows meet the gene boxes
        parts.append(coupling(f"M{sources[0]},61 L{e_top},{BOX_Y - 4}", (sources[0] + e_top) / 2 - 30, 80 if group.startswith("H2") else 94))
        parts.append(coupling(f"M{sources[1]},61 L{g_top},{BOX_Y - 4}", (sources[1] + g_top) / 2 + 30, 94))
    # Crop the drawing to what this group shows.
    top = 0 if group.startswith("H") else 34 if group in ("E_to_G", "reciprocal") else BOX_Y - 8
    bottom = 250 if group in ("G_to_E", "reciprocal") else BOX_Y + 70
    return (f'<svg viewBox="0 {top} 390 {bottom - top}" role="img" aria-label="{group} schematic">'
            f'{"".join(parts)}</svg>')


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------

def chips(values):
    return "".join(f'<span class="chip">{html.escape(v)}</span>' for v in values)


def regulated_rates(rows):
    """Which rates the group's models regulate: a list, a source-state by rate matrix (one-way models),
    or an E→G by G→E matrix (reciprocal models)."""
    if rows[0]["group"] in ("E_to_G", "G_to_E"):
        side = "EG" if rows[0]["group"] == "E_to_G" else "GE"
        states = list(dict.fromkeys(m[f"{side}_state"] for m in rows))
        sets = list(dict.fromkeys(m[f"{side}_rates"] for m in rows))
        have = {(m[f"{side}_state"], m[f"{side}_rates"]) for m in rows}
        head = "".join(f'<th class="rot"><span>{html.escape(r)}</span></th>' for r in sets)
        body = "".join(f'<tr><th>source state {html.escape(s)}</th>' +
                       "".join(f'<td class="{"on" if (s, r) in have else "off"}"></td>' for r in sets) + "</tr>" for s in states)
        return (f'<p class="what">Rows = state of the regulating gene, columns = rates multiplied. '
                f'A filled cell is one model ({len(have)} models).</p>'
                f'<table class="matrix"><tr><th></th>{head}</tr>{body}</table>')
    if rows[0]["group"] == "reciprocal":
        eg = list(dict.fromkeys(m["EG_rates"] for m in rows))
        ge = list(dict.fromkeys(m["GE_rates"] for m in rows))
        have = {(m["EG_rates"], m["GE_rates"]) for m in rows}
        head = "".join(f'<th class="rot"><span>{html.escape(q)}</span></th>' for q in ge)
        body = "".join(f'<tr><th>{html.escape(r)}</th>' +
                       "".join(f'<td class="{"on" if (r, q) in have else "off"}"></td>' for q in ge) + "</tr>" for r in eg)
        states = list(dict.fromkeys(f'{m["EG_state"]} / {m["GE_state"]}' for m in rows))
        return (f'<p class="what">Rates multiplied: rows = rates of G regulated by E, columns = rates of E regulated by G. '
                f'Every filled cell is simulated with each of these {len(states)} source-state combinations '
                f'(state of E / state of G), {len(rows)} models in all:</p><p>{chips(states)}</p>'
                f'<table class="matrix"><tr><th></th>{head}</tr>{body}</table>')
    key = "EG_rates" if rows[0]["EG_rates"] else "GE_rates" if rows[0]["GE_rates"] else "HE_rates"
    return f'<p class="what">Rates multiplied, one model each ({len(rows)} models):</p><p>{chips(m[key] for m in rows)}</p>'


def parameters_of(row):
    """(column, label) of every parameter in a row of models.csv or of the simulation list, including the folds of
    single regulated rates (fold_EG_R12: "fold of E regulating G, rate R12") of the extensive sweep."""
    rate_folds = [(k, f'{PARAMETERS["_".join(k.split("_")[:2])]}, rate {k.split("_")[2]}')
                  for k in row if k.startswith("fold_") and k.count("_") == 2]
    return list(PARAMETERS.items()) + rate_folds


def parameter_grid(m, mode):
    """The parameter values one model of the group is tested with."""
    rows, counts = [], []
    for name, label in parameters_of(m):
        text = m.get(name) or ""
        if not text:
            continue
        if text.startswith("="):
            cell = f'<span class="tied">{html.escape(text)}</span>'
        elif "values" in text:
            cell = f'<span class="range">{html.escape(text)}</span>'
        else:
            values = text.split("; ")
            counts.append(str(len(values)))
            cell = chips(values)
        rows.append(f"<tr><th>{label}</th><td>{cell}</td></tr>")
    how = (f'{m["n_simulations"]} random combinations per model' if mode == "extensive"
           else f'{m["n_simulations"]} simulations per model' if len(counts) == 1
           else f'every combination: {" × ".join(counts)} = {m["n_simulations"]} simulations per model')
    return f'<table class="params">{"".join(rows)}</table><p class="what">{how}</p>'


# ---------------------------------------------------------------------------

STYLE = """
:root{--bg:#f9f9f7;--card:#fcfcfb;--ink:#0b0b0b;--ink2:#52514e;--muted:#898781;--line:#e1e0d9;--accent:#2a78d6;--fill:#cde2fb}
@media (prefers-color-scheme: dark){:root{--bg:#0d0d0d;--card:#1a1a19;--ink:#f0efec;--ink2:#c3c2b7;--line:#2c2c2a;--accent:#3987e5;--fill:#184f95}}
body{margin:0;background:var(--bg);color:var(--ink);font:15px/1.5 system-ui,-apple-system,"Segoe UI",sans-serif}
main{max-width:1080px;margin:auto;padding:24px 16px 48px}
h1{font-size:26px;margin:0 0 4px} h2{font-size:18px;margin:0 0 2px} p{margin:6px 0}
.sub,.what{color:var(--ink2)} .what{font-size:13px}
table{border-collapse:collapse} th,td{text-align:left;padding:5px 10px;border-bottom:1px solid var(--line);font-weight:400}
.overview{width:100%;margin:14px 0 22px;background:var(--card);border:1px solid var(--line);border-radius:8px}
.overview th{color:var(--ink2);font-size:13px} .overview td.n,.overview th.n{text-align:right;font-variant-numeric:tabular-nums}
.overview tr:last-child td{border-bottom:0;font-weight:600}
.card{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:14px 16px;margin:14px 0;
      display:grid;grid-template-columns:minmax(260px,390px) 1fr;gap:8px 22px;align-items:start}
.card header{grid-column:1/-1} .scroll{overflow-x:auto}
@media (max-width:760px){.card{grid-template-columns:1fr}}
.chip{display:inline-block;border:1px solid var(--line);border-radius:6px;padding:1px 8px;margin:2px 4px 2px 0;
      font-size:13px;font-variant-numeric:tabular-nums;white-space:nowrap}
.params th{color:var(--ink2);font-size:13px;white-space:nowrap} .params td,.params th{border-bottom:0;padding:3px 10px 3px 0}
.tied{color:var(--ink2);font-size:13px}
.range{display:inline-block;font-size:13px;padding:1px 10px;border-radius:6px;background:var(--fill)}
.matrix th{font-size:12px;color:var(--ink2);border:0;padding:1px 6px 1px 0;white-space:nowrap}
.matrix td{width:16px;height:16px;padding:0;border:2px solid var(--card)} .matrix td.on{background:var(--accent)} .matrix td.off{background:var(--line)}
.matrix th.rot{height:112px;vertical-align:bottom;padding:0} .matrix th.rot span{display:block;width:16px;writing-mode:vertical-rl;transform:rotate(180deg)}
svg{width:100%;height:auto} svg text{font:12px system-ui,sans-serif;fill:var(--ink)}
.box{fill:none;stroke:var(--muted);stroke-width:1.2} .state{fill:var(--card);stroke:var(--ink2);stroke-width:1.4}
.thin{stroke:var(--ink2);stroke-width:1.1} .label{font-weight:600} .num,.mid{text-anchor:middle}
svg text.small{font-size:10.5px;fill:var(--ink2)} .couple{fill:none;stroke:var(--accent);stroke-width:2.2}
svg text.fold{fill:var(--ink);font-weight:600;font-size:11.5px;text-anchor:middle}
"""

MARKERS = ('<svg width="0" height="0" style="position:absolute"><defs>'
           '<marker id="tip" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="6" markerHeight="6" orient="auto">'
           '<path d="M0,0 L8,4 L0,8 z" fill="#898781"/></marker>'
           '<marker id="ctip" viewBox="0 0 8 8" refX="6" refY="4" markerWidth="5" markerHeight="5" orient="auto">'
           '<path d="M0,0 L8,4 L0,8 z" fill="#2a78d6"/></marker></defs></svg>')


def expanded_page(config, run_dir, profile):
    """Summary page of an expanded run (1_generate_models_expanded.py): too many architectures to draw one by one."""
    import pandas as pd
    models = pd.read_csv(run_dir / "expanded_models.csv", dtype=str)
    points = pd.read_csv(run_dir / "expanded_simulations.csv",
                         usecols=lambda column: column == "architecture_id" or column in PARAMETERS or column.startswith("fold_"))
    points = points.merge(models[["architecture_id", "group", "is_control"]], on="architecture_id")
    replicates = config["profiles"][profile]["replicates"]
    source = lambda edge: models[f"{edge}_states"].mask(models[f"{edge}_source_kind"] == "reporter_count", "R")
    overview, cards = [], []
    for group, rows in models.groupby("group", sort=False):
        swept = rows[rows["is_control"] != "true"]
        sims = points[points["group"] == group]
        overview.append(f'<tr><td>{group}</td><td class="n">{len(swept):,}</td><td class="n">{len(rows) - len(swept)}</td>'
                        f'<td class="n">{len(sims):,}</td><td class="n">{len(sims) * replicates:,}</td></tr>')
        if group == "uncoupled":
            cards.append(f'<section class="card"><header><h2>{group}</h2><p class="sub">{GROUP_TEXT[group]}</p></header>'
                         f'<div>{schematic(group)}</div></section>')
            continue
        lines = []
        for edge, text in [("EG", "E regulates G"), ("GE", "G regulates E"), ("HE", "H regulates E"), ("HG", "H regulates G")]:
            used = swept[swept[f"{edge}_rates"].notna()]
            if len(used):
                lines.append(f'<p class="what">{text}, from {source(edge)[used.index].nunique()} source options:</p>'
                             f'<p>{chips(dict.fromkeys(source(edge)[used.index]))}</p>'
                             f'<p class="what">regulated rates, {used[f"{edge}_rates"].nunique()} subsets:</p>'
                             f'<p>{chips(dict.fromkeys(used[f"{edge}_rates"]))}</p>')
        noncontrol = sims[sims["architecture_id"].isin(swept["architecture_id"])]
        ranges = ""                                     # per connection fold (its single-rate folds pooled) and hidden rate
        for name, label in PARAMETERS.items():
            columns = [c for c in noncontrol if c == name or c.startswith(name + "_")]
            values = pd.concat([noncontrol[c] for c in columns]).dropna()
            if len(values):
                each = " (one per regulated rate)" if any(c != name and noncontrol[c].notna().any() for c in columns) else ""
                ranges += f"<tr><th>{label}{each}</th><td>{values.min():g} to {values.max():g}</td></tr>"
        combined = ("Every architecture combines one choice of each connection above: "
                    f"{len(swept):,} architectures, each with {int(noncontrol.groupby('architecture_id').size().max())} parameter sets.")
        cards.append(f'<section class="card"><header><h2>{group}</h2><p class="sub">{GROUP_TEXT[group]} '
                     f'A source with several states (e.g. 1+2) regulates while it is in any of them.</p></header>'
                     f'<div>{schematic(group)}</div><div class="scroll">{"".join(lines)}<p class="what">{combined}</p>'
                     f'<table class="params">{ranges}</table>'
                     f'<p class="what">Plus 1 control with fold 1 (gamma = 0), simulated once.</p></div></section>')
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Expanded coupled-model sweep</title><style>{STYLE}</style></head>
<body>{MARKERS}<main><h1>Expanded coupled-model sweep: what will be simulated</h1>
<p class="sub">Profile {html.escape(profile)} · {(models["is_control"] != "true").sum():,} architectures and {(models["is_control"] == "true").sum()} controls · {len(points):,} parameter sets × {replicates} replicate datasets ·
E = CCAT1 (enhancer), G = MYC (gene), H = hidden unit. "× fold" multiplies the chosen rates while the source is in a regulating state:
fold above 1 activates, fold below 1 inhibits. From R every occupied reporter position adds (fold − 1).</p>
<table class="overview"><tr><th>Group</th><th class="n">Architectures</th><th class="n">Controls</th><th class="n">Parameter sets</th>
<th class="n">Datasets</th></tr>{"".join(overview)}
<tr><td>Total</td><td class="n">{(models["is_control"] != "true").sum():,}</td><td class="n">{(models["is_control"] == "true").sum()}</td>
<td class="n">{len(points):,}</td><td class="n">{len(points) * replicates:,}</td></tr></table>
{"".join(cards)}</main></body></html>"""
    out = run_dir / "models_report.html"
    out.write_text(page)
    print(f"Wrote {out}")


def main():
    arguments = sys.argv[1:]
    profile = "full"
    if "--profile" in arguments:
        at = arguments.index("--profile")
        profile = arguments[at + 1]
        del arguments[at:at + 2]
    if len(arguments) != 1:
        raise SystemExit(__doc__)
    config_path = Path(arguments[0]).resolve()
    config = tomllib.loads(config_path.read_text())
    if "profiles" in config:                            # config_expanded.toml of 1_generate_models_expanded.py
        if profile not in config["profiles"]:
            raise SystemExit(f"no profile {profile}; available: {', '.join(config['profiles'])}")
        return expanded_page(config, config_path.parent / config["output"]["run_root"] / profile, profile)
    run_dir = config_path.parent / config["output"]["run_dir"]
    mode, replicates = config["sweep"]["mode"], config["simulation"]["replicates"]
    models = list(csv.DictReader(open(run_dir / "models.csv")))
    groups = {}
    for m in models:
        groups.setdefault(m["group"], []).append(m)

    overview, cards = [], []
    for group, rows in groups.items():
        swept = [m for m in rows if m["control"] != "yes" or group == "uncoupled"]            # without the gamma = 0 control
        example = swept[0]
        n_sims = sum(int(m["n_simulations"]) for m in rows)
        names = ", ".join(dict.fromkeys("_".join(k.split("_")[:2]) if k.startswith("fold_") else k for k, _ in parameters_of(example)
                                        if example.get(k) and not example[k].startswith("="))) if group != "uncoupled" else "none"
        controls = len(rows) - len(swept)
        overview.append(f'<tr><td>{group}</td><td class="n">{len(swept)}</td><td>{names}</td>'
                        f'<td class="n">{example["n_simulations"]}</td><td class="n">{controls}</td><td class="n">{n_sims:,}</td></tr>')
        detail = "" if group == "uncoupled" else (
            f'<div class="scroll">{regulated_rates(swept)}{parameter_grid(example, mode)}'
            f'<p class="what">Plus {controls} control with fold 1 (gamma = 0), simulated once.</p></div>')
        cards.append(f'<section class="card"><header><h2>{group}</h2><p class="sub">{GROUP_TEXT[group]}</p></header>'
                     f'<div>{schematic(group)}</div>{detail}</section>')

    total = sum(int(m["n_simulations"]) for m in models)
    page = f"""<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>Coupled-model sweep</title><style>{STYLE}</style></head>
<body>{MARKERS}<main><h1>Coupled-model sweep: what will be simulated</h1>
<p class="sub">{len(models)} models · {total:,} simulations × {replicates} replicate datasets · sweep mode: {mode} ·
E = CCAT1 (enhancer), G = MYC (gene), H = hidden unit. "× fold" multiplies the chosen rates while the source is in its regulating state:
fold above 1 activates, fold below 1 inhibits. Gene-to-gene models are tested with every fold and its inverse; hidden models only activate.</p>
<table class="overview"><tr><th>Group</th><th class="n">Models</th><th>Parameters swept</th><th class="n">Parameter sets per model</th>
<th class="n">Fold-1 controls</th><th class="n">Simulations</th></tr>{"".join(overview)}
<tr><td>Total</td><td class="n"></td><td></td><td class="n"></td><td class="n"></td><td class="n">{total:,}</td></tr></table>
{"".join(cards)}</main></body></html>"""
    out = run_dir / "models_report.html"
    out.write_text(page)
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
