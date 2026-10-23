import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="SAC explorer", css_file="sac_explorer.css",
                 html_head_file="sac_explorer_head.html")


@app.cell
def _():
    import os
    import re

    import altair as alt
    import duckdb
    import marimo as mo
    import numpy as np
    import pandas as pd

    # the website's metal colours (site/template.html), so a metal looks the same everywhere
    METAL_COLOURS = {"Fe": "#c4682f", "Co": "#3f6fc0", "Ni": "#3f9a6a", "Cu": "#1f8f8a", "Zn": "#7a8796",
                     "Mo": "#8a68b0", "W": "#5d5a93", "Mn": "#b0427a", "V": "#b28a12"}
    SUPPORT_COLOURS = {"graphene": "#2c6bd3", "MoS2": "#8a68b0", "SiO2": "#b28a12"}

    @alt.theme.register("sac", enable=True)
    def _sac_theme():
        # mid-tone text so the charts read on the light and on the dark marimo theme
        _ink, _muted, _grid = "#5d6b70", "#7d8b91", "#e1e7ea"
        return {"config": {
            "background": "transparent", "font": "Source Sans 3", "view": {"stroke": None},
            "axis": {"labelColor": _muted, "titleColor": _ink, "gridColor": _grid, "domainColor": "#c9d2d6",
                     "tickColor": "#c9d2d6", "labelFontSize": 12, "titleFontSize": 12, "titleFontWeight": 600},
            "legend": {"labelColor": _ink, "titleColor": _ink, "labelFontSize": 12, "symbolSize": 90},
            "title": {"font": "Newsreader", "fontSize": 16, "color": _ink},
            "range": {"category": list(METAL_COLOURS.values())}}}
    return METAL_COLOURS, SUPPORT_COLOURS, alt, duckdb, mo, np, os, pd, re


@app.cell
def _(ENZYME_METALS, dft, mo, n_entries, rx, sites, xtb):
    # header: what this is, and the numbers that matter (independent of the filters)
    _rep = sites[sites.is_representative & sites.metal.isin(ENZYME_METALS)]
    _d = dft.dropna(subset=["dG_H"])
    _near = int((_d.dG_H.round(2).abs() <= 0.2).sum())
    _rx_sites = int(rx.filter(like="U_").notna().any(axis=1).sum()) if not rx.empty else 0
    _facts = [
        (f"{len(sites):,}", f"metal sites in {n_entries:,} PDB entries", ""),
        (f"{_rep.groupby(['metal', 'formula']).ngroups:,}", "distinct first shells of the nine enzyme metals", ""),
        (f"{xtb[xtb.source != 'library'].motif.nunique():,}",
         f"graphene sites screened with GFN2-xTB, plus {xtb[xtb.source == 'library'].motif.nunique():,} library sites", ""),
        (f"{len(_d)}", "H positions computed in DFT (results as of 22 Oct, 23:30)", ""),
        (f"{_near}", "sites within 0.2 eV of the optimum in DFT", "good"),
        (f"{_rx_sites}", "sites with O₂, CO₂ and N₂ intermediates (xTB)", ""),
    ]
    _tiles = "".join(f'<div class="sac-fact {c}"><b>{v}</b><span>{t}</span></div>' for v, t, c in _facts)
    mo.Html(f"""
    <header class="sac-hero">
      <div class="sac-eyebrow">Hackathon project · draft results · DFT as of 22 Oct, 23:30; the rest GFN2-xTB</div>
      <h1>Enzyme sites as single-atom catalysts</h1>
      <p class="sac-lede">Every transition-metal site in the Protein Data Bank, reduced to its first coordination shell and
      rebuilt as a single-atom site on graphene, MoS₂ or silica; screened for hydrogen evolution, then for O₂, CO₂ and N₂
      electrochemistry. Filters in the sidebar apply to every tab.</p>
      <div class="sac-facts">{_tiles}</div>
    </header>""")
    return


@app.cell
def _(METAL_COLOURS, calib, dft_window, metals, min_entries, mo, se_cl, xtb_window):
    _legend = "".join(f'<span><i style="background:{c}"></i>{m}</span>' for m, c in METAL_COLOURS.items())
    mo.sidebar([
        mo.Html('<p class="sac-side-title">Filters</p>'),
        metals, calib, xtb_window, dft_window, min_entries, se_cl,
        mo.Html('<p class="sac-note">Windows in eV around ΔG<sub>H</sub> = 0; min. entries: PDB entries a donor '
                'environment needs to be screened.</p>'),
        mo.Html(f'<p class="sac-side-title" style="margin-top:10px">Metal colours</p><div class="sac-legend">{_legend}</div>'),
        mo.Html('<p class="sac-note">DFT as of the last pull from the cluster (22 Oct, 23:30); everything else is '
                'GFN2-xTB. All values are drafts that need expert checking.</p>'),
    ])
    return


@app.cell
def _(ask_view, browse_view, calib_view, candidates_view, detail_view, dft_view, funnel_view, get_tab, mo,
      reactions_view, rx_cal_view, rx_controls, set_tab, sql_view, sql_text, support_controls, supports_view,
      zoom_html):
    _zoom = (mo.iframe('<!doctype html><html data-theme="light"><meta charset="utf-8"><meta name="viewport" '
                       'content="width=device-width, initial-scale=1">' + zoom_html, height="860px") if zoom_html else
             mo.callout(mo.md("Build the 3D view first: `python site/zoom_viewer.py`."), kind="warn"))
    mo.ui.tabs({
        "From enzyme to catalyst (3D)": _zoom,
        "Hydrogen evolution": mo.vstack([funnel_view, candidates_view, detail_view,
                                         mo.hstack([calib_view, dft_view], widths="equal", wrap=True, align="start")],
                                        gap=1.5),
        "MoS₂ and silica": mo.vstack([support_controls, supports_view], gap=1),
        "O₂ · CO₂ · N₂": mo.vstack([rx_controls, reactions_view, rx_cal_view], gap=1.2),
        "Ask and browse": mo.vstack([ask_view, mo.accordion({"Write the SQL yourself": mo.vstack([sql_text, sql_view])}),
                                     browse_view], gap=1.2),
    }, value=get_tab(), on_change=set_tab)
    return


@app.cell
def _(mo):
    get_tab, set_tab = mo.state("From enzyme to catalyst (3D)")    # the open tab survives a re-render of the tabs
    return get_tab, set_tab


@app.cell
def _(mo, os, pd):
    # data: the release database (Parquet) and the result tables of the pipeline
    ROOT = os.path.dirname(str(mo.notebook_dir()))
    DATA = os.path.join(ROOT, "data")
    REL = os.path.join(DATA, "release", "parquet")
    sites = pd.read_parquet(os.path.join(REL, "sites.parquet"), columns=[
        "site_id", "pdb_id", "metal", "is_representative", "formula", "graphene_motif", "preferred_site",
        "dv_formula", "dv_isomer", "sv_formula", "reaction"])
    motifs = pd.read_parquet(os.path.join(REL, "motifs.parquet"))
    n_entries = len(pd.read_parquet(os.path.join(REL, "processing_log.parquet"), columns=["pdb_id"]))   # swept entries
    xtb = pd.read_csv(os.path.join(DATA, "xtb_scan.csv"))                  # scan and library (column source)
    offsets = pd.read_csv(os.path.join(DATA, "xtb_offsets.csv"))           # screen/calibrate.py
    pairs = pd.read_csv(os.path.join(DATA, "xtb_calibration.csv"))
    dft = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv"))
    novelty = pd.read_csv(os.path.join(DATA, "novelty.csv"))
    ENZYME_METALS = ["V", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W"]

    def _read(path, reader):
        return reader(path) if os.path.exists(path) else pd.DataFrame()

    # MoS2 and silica (bio/supports.py, screen/xtb_supports.py) and the O2/CO2/N2 screen (screen/xtb_reactions.py)
    support_motifs = _read(os.path.join(REL, "support_motifs.parquet"), pd.read_parquet)
    xtb_supports = _read(os.path.join(DATA, "xtb_supports.csv"), pd.read_csv)
    rx = _read(os.path.join(DATA, "xtb_reactions.csv"), pd.read_csv)
    rx_cal = _read(os.path.join(DATA, "xtb_orr_calibration.csv"), pd.read_csv)
    # the 3D zoom from enzyme to support (site/zoom_viewer.py)
    _zoom = os.path.join(ROOT, "site", "zoom_viewer.html")
    zoom_html = open(_zoom).read() if os.path.exists(_zoom) else ""
    return (DATA, ENZYME_METALS, REL, dft, motifs, n_entries, novelty, offsets, pairs, rx, rx_cal, sites, support_motifs, xtb,
            xtb_supports, zoom_html)


@app.cell
def _(re):
    # motif names: 'Fe-N2SO:DV:NSNO' = Fe with two N, one S, one O in a double vacancy, donor order N,S,N,O
    _SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
    _ORDER = ["N", "S", "O", "C", "Se", "Cl"]

    def pretty(motif):
        head, site = (motif.split(":") + [""])[:2]
        metal, donors = head.split("-")
        return f"{metal}–{donors.translate(_SUB)} ({site})"

    def canonical(motif):
        """Same site, whatever the order the donors were written in: 'Fe-C2S2:DV:cis5' -> 'Fe|S2C2|DV|cis5|'."""
        head, *opts = motif.split(":")
        metal, donors = head.split("-")
        n = {el: int(k or 1) for el, k in re.findall(r"([A-Z][a-z]?)(\d*)", donors)}
        formula = "".join(el + (str(n[el]) if n[el] > 1 else "") for el in _ORDER if el in n)
        iso = next((o for o in opts[1:] if not o.startswith("ax=")), "")
        ax = next((o for o in opts[1:] if o.startswith("ax=")), "")
        return "|".join([metal, formula, opts[0] if opts else "", iso, ax])

    def builder_isomer(formula, isomer):
        """Sweep isomer (cis/trans) -> the site builder's isomer, as in screen/xtb_scan.py."""
        n = {el: int(k or 1) for el, k in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}
        counts = sorted(n.values())
        if counts == [2, 2]:
            return "trans" if isomer == "trans" else "cis5"
        if counts == [1, 1, 2]:
            pair = next(el for el, k in n.items() if k == 2)
            a, b = [el for el, k in n.items() if k == 1]
            return pair + a + b + pair if isomer == "trans" else pair + a + pair + b
        if counts == [1, 1, 1, 1]:
            return "".join(n)
        return None

    def site_label(kind, site, metal):
        """Readable name of a site of the reaction screen."""
        if kind == "graphene":
            return pretty(site)
        return f"{metal}–S₃ (MoS₂)" if kind == "MoS2" else f"{metal} (silica)"
    return builder_isomer, canonical, pretty, site_label


@app.cell
def _(ENZYME_METALS, mo):
    min_entries = mo.ui.slider(1, 50, value=5, label="min. entries", show_value=True)
    se_cl = mo.ui.switch(value=False, label="Se and Cl donors")
    metals = mo.ui.multiselect(options=ENZYME_METALS, value=ENZYME_METALS, label="metals")
    calib = mo.ui.dropdown(options=["one offset per metal (used)", "one line for all metals", "raw xTB"],
                           value="one offset per metal (used)", label="H calibration")
    xtb_window = mo.ui.slider(0.05, 1.0, step=0.05, value=0.3, label="xTB |ΔG_H| ≤", show_value=True)
    dft_window = mo.ui.slider(0.05, 0.5, step=0.05, value=0.2, label="DFT |ΔG_H| ≤", show_value=True)
    return calib, dft_window, metals, min_entries, se_cl, xtb_window


@app.cell
def _(metals, min_entries, np, se_cl, sites):
    # graphene donor environments of the unique sites of the chosen metals, counted in PDB entries
    rep = sites[sites.is_representative & sites.metal.isin(metals.value)]
    _mapped = rep[rep.preferred_site.isin(["DV", "SV"])].copy()
    _mapped["env_formula"] = np.where(_mapped.preferred_site == "DV", _mapped.dv_formula, _mapped.sv_formula)
    _mapped["env_isomer"] = np.where(_mapped.preferred_site == "DV", _mapped.dv_isomer.fillna(""), "")
    env = (_mapped.groupby(["env_formula", "preferred_site", "env_isomer"]).pdb_id.nunique()
           .rename("entries").reset_index())
    _keep = env.entries >= min_entries.value
    if not se_cl.value:
        _keep &= ~env.env_formula.str.contains("Se|Cl")
    frequent = env[_keep]
    return env, frequent, rep


@app.cell
def _(calib, np, offsets, pairs):
    # xTB -> DFT calibration: the offsets of screen/calibrate.py (data/xtb_offsets.csv; a metal needs two DFT pairs)
    _p = pairs.copy()
    _p["metal"] = _p.motif_dft.str.split("-").str[0]
    _used = offsets[offsets.used]
    _off = _used[_used.state == "H_M"].set_index("metal").offset
    _off_s = _used[_used.state == "H_S"].offset.iloc[0]
    _a, _b = np.polyfit(_p.dG_H_xtb, _p.dG_H_dft, 1)

    def calibrate(metal, state, value):
        if calib.value == "raw xTB":
            return value
        if calib.value == "one line for all metals":
            return _a * value + _b
        if state == "H_S":
            return value + _off_s
        return value + _off[metal] if metal in _off.index else np.nan   # no DFT pair for this metal yet

    _p["dG_cal"] = [calibrate(m, s, v) for m, s, v in zip(_p.metal, _p.state_dft, _p.dG_H_xtb)]
    calib_pairs = _p
    if calib.value == "one offset per metal (used)":      # leave-one-out: each site predicted from the others
        calib_mae = float((_used.n_pairs * _used.loo_mae).sum() / _used.n_pairs.sum())
    else:                                                  # one line, or raw: in-sample
        calib_mae = float((_p.dG_cal - _p.dG_H_dft).abs().mean())
    return calib_mae, calib_pairs, calibrate


@app.cell
def _(builder_isomer, calibrate, canonical, dft, frequent, metals, novelty, np, pd, pretty, xtb, xtb_window):
    # every frequent environment with every chosen metal, its xTB result and the calibrated value
    _rows = []
    for _e in frequent.itertuples():
        _iso = builder_isomer(_e.env_formula, _e.env_isomer) if _e.preferred_site == "DV" else None
        for _m in metals.value:
            _rows.append({"motif": f"{_m}-{_e.env_formula}:{_e.preferred_site}" + (f":{_iso}" if _iso else ""),
                          "metal": _m, "entries": int(_e.entries)})
    to_screen = pd.DataFrame(_rows, columns=["motif", "metal", "entries"])
    scanned = to_screen[to_screen.motif.isin(xtb.motif)]
    _x = xtb.merge(to_screen, on="motif")
    _ok = _x.converged.astype(bool) & _x.bare_converged.astype(bool) & _x.dG_H.notna()
    _intended = np.where(_x.state == "H_M", _x.H_partner == "M", _x.H_partner.astype(str).str.endswith("(donor)"))
    _n_c = (_x.motif.str.split("-").str[1].str.split(":").str[0].str.extract(r"C(\d*)(?![a-z])")[0]
            .map(lambda v: 0 if pd.isna(v) else int(v or 1)))
    _x = _x[_ok & _intended & (_x.max_M_X_bare <= 2.6) & (_n_c < 2)].copy()   # quality filters of the scan
    _x["dG_cal"] = [calibrate(m, s, v) for m, s, v in zip(_x.metal, _x.state, _x.dG_H)]
    _hand = {canonical(m) for m, o in zip(dft.motif, dft.enzyme) if not str(o).startswith("xTB scan")}
    _c = _x[(_x.dG_cal.abs() <= xtb_window.value) & ~_x.motif.map(canonical).isin(_hand)].copy()
    _c["env"] = _c.motif.str.split("-").str[1].str.split(":").str[:2].str.join(":")
    _c = _c.sort_values("dG_cal", key=abs).drop_duplicates(["metal", "env"])   # as in screen/calibrate.py
    _d = dft.dropna(subset=["dG_H"])
    _c["dG_DFT"] = (_c.motif.map(canonical) + "#" + _c.state).map(
        dict(zip(_d.motif.map(canonical) + "#" + _d.state, _d.dG_H)))
    _lit = dict(zip(novelty.motif.map(lambda m: canonical(m).split("|")[1]) + "|" + novelty.motif.str.split("-").str[0],
                    novelty.openalex_all))
    _c["papers"] = (_c.motif.map(lambda m: canonical(m).split("|")[1]) + "|" + _c.metal).map(_lit)
    _c["site"] = _c.motif.map(pretty)
    _c["H on"] = np.where(_c.state == "H_S", "S", _c.metal)
    candidates = _c.sort_values("dG_cal", key=abs)
    return candidates, scanned, to_screen


@app.cell
def _(alt, candidates, dft, dft_window, env, frequent, metals, mo, n_entries, pd, rep, scanned, sites, to_screen):
    _d = dft[dft.motif.str.split("-").str[0].isin(metals.value)].dropna(subset=["dG_H"])
    _near = _d[_d.dG_H.round(2).abs() <= dft_window.value]
    _f = pd.DataFrame([
        (n_entries, "PDB entries with a transition metal", "selection"),
        (len(sites), "metal sites, every copy", "more per entry or site"),
        (int(sites.is_representative.sum()), "unique sites", "merging"),
        (len(rep), "unique sites of the chosen metals", "selection"),
        (rep.groupby(["metal", "formula"]).ngroups, "distinct first shells", "merging"),
        (len(env), "distinct graphene sites", "merging"),
        (len(frequent), "frequent enough, donor filter", "selection"),
        (len(to_screen), "each with every chosen metal", "more per entry or site"),
        (len(scanned), "of these, scanned with xTB", "selection"),
        (len(candidates), "xTB within the window (new sites)", "selection"),
        (_d.motif.nunique(), "sites with DFT values (incl. hand-picked)", "selection"),
        (_near.motif.nunique(), "DFT within the window", "selection"),
    ], columns=["count", "step", "kind"])
    _f["bar"] = _f["count"].clip(lower=1)
    _order = list(_f.step)
    _bars = alt.Chart(_f).mark_bar(cornerRadiusEnd=3, height=14).encode(
        x=alt.X("bar:Q", scale=alt.Scale(type="log", domain=[0.9, 1e6], nice=False), title="count (log scale)"),
        x2=alt.datum(0.9),
        y=alt.Y("step:N", sort=_order, title=None),
        color=alt.Color("kind:N", title=None, legend=alt.Legend(orient="bottom"),
                        scale=alt.Scale(domain=["merging", "selection", "more per entry or site"],
                                        range=["#c7cfd3", "#2c6bd3", "#8f9ca2"])),
        tooltip=["step", alt.Tooltip("count:Q", format=",")])
    _labels = alt.Chart(_f).mark_text(align="left", dx=4, font="IBM Plex Mono", fontSize=11, color="#5d6b70").encode(
        x="bar:Q", y=alt.Y("step:N", sort=_order), text=alt.Text("count:Q", format=","))
    _missing = len(to_screen) - len(scanned)
    funnel_view = mo.vstack([
        mo.Html('<h2 class="sac-h2">From the PDB to candidates</h2><p class="sac-intro">Most of the drop is merging, not rejecting: '
                'the PDB holds the same few coordination motifs thousands of times. The real filters are frequency, '
                'the xTB screen and the window around ΔG<sub>H</sub> = 0.</p>'),
        (_bars + _labels).properties(width=520, height=360),
        mo.callout(mo.md(f"**{_missing} site(s)** pass the filters but have no xTB result yet: the scan covered "
                         "environments seen in at least 5 entries, without Se and Cl."), kind="warn")
        if _missing else mo.md(""),
    ])
    return (funnel_view,)


@app.cell
def _(candidates, mo):
    _show = candidates[["site", "H on", "dG_H", "dG_cal", "dG_DFT", "entries", "papers", "motif"]].rename(columns={
        "dG_H": "xTB raw (eV)", "dG_cal": "xTB calibrated (eV)", "dG_DFT": "DFT (eV)", "entries": "PDB entries"})
    table = mo.ui.table(_show.round(2).reset_index(drop=True), selection="single", page_size=8,
                        label=f"{len(_show)} new xTB candidates in the window (best isomer and H position per "
                              "metal and site; hand-picked DFT sites left out); select one")
    candidates_view = mo.vstack([mo.Html("<h2 class='sac-h2'>Candidates</h2>"), table])
    return candidates_view, table


@app.cell
def _(mo, motifs, re, table):
    # where the selected site's donor environment occurs in the PDB, with any metal
    if table.value is None or len(table.value) == 0:
        detail_view = mo.Html('<p class="sac-intro"><i>Select a candidate to see where its first shell occurs in the '
                              'PDB.</i></p>')
    else:
        _r = table.value.iloc[0]
        _formula, _site = _r["motif"].split("-", 1)[1].split(":")[:2]
        _m = motifs[(motifs.set_name == "complete") & (motifs.kind == "graphene")
                    & motifs.motif.str.match(r"^\w+-%s\(%s\)" % (re.escape(_formula), _site))]
        _m = _m.sort_values("n_entries_rep", ascending=False)
        detail_view = mo.vstack([
            mo.md(f"### {_r['site']}, H on {_r['H on']}\n"
                  f"xTB {_r['xTB calibrated (eV)']:+.2f} eV calibrated"
                  + (f", DFT {_r['DFT (eV)']:+.2f} eV" if _r['DFT (eV)'] == _r['DFT (eV)'] else ", no DFT yet")
                  + ". This donor environment occurs in the PDB with these metals:"),
            mo.ui.table(_m[["motif", "n_entries_rep", "n_sites_rep", "reactions", "example_pdb"]].rename(columns={
                "n_entries_rep": "entries", "n_sites_rep": "unique sites"}).reset_index(drop=True),
                selection=None, page_size=6)])
    return (detail_view,)


@app.cell
def _(METAL_COLOURS, alt, calib, calib_mae, calib_pairs, mo, pd):
    _lo, _hi = -2.5, 2.5
    _pts = alt.Chart(calib_pairs.dropna(subset=["dG_cal"])).mark_circle(size=90, opacity=0.9).encode(
        x=alt.X("dG_cal:Q", title=f"xTB, {calib.value} (eV)", scale=alt.Scale(domain=[_lo, _hi])),
        y=alt.Y("dG_H_dft:Q", title="DFT (eV)", scale=alt.Scale(domain=[_lo, _hi])),
        color=alt.Color("metal:N", title=None, scale=alt.Scale(domain=list(METAL_COLOURS), range=list(METAL_COLOURS.values()))),
        tooltip=["motif_dft", "state_dft", alt.Tooltip("dG_H_xtb:Q", format="+.2f"),
                 alt.Tooltip("dG_cal:Q", format="+.2f"), alt.Tooltip("dG_H_dft:Q", format="+.2f")])
    _diag = alt.Chart(pd.DataFrame({"x": [_lo, _hi], "y": [_lo, _hi]})).mark_line(color="#8f9ca2", strokeDash=[4, 4]).encode(
        x="x:Q", y="y:Q")
    calib_view = mo.vstack([
        mo.Html(f'<h3 class="sac-h3">xTB against DFT</h3><p class="sac-intro">Mean absolute error <b>{calib_mae:.2f} eV</b> over '
                f'{len(calib_pairs)} sites that have both ' + ('(leave-one-out: each site predicted from the offset of '
                'the others; a metal needs two DFT pairs, so Mo and W stay raw).</p>' if calib.value.startswith("one offset")
                else '(in-sample).</p>')),
        (_diag + _pts).properties(width=380, height=340)])
    return (calib_view,)


@app.cell
def _(METAL_COLOURS, alt, dft, dft_window, metals, mo, pd, pretty):
    _d = dft.dropna(subset=["dG_H"]).copy()
    _d = _d[_d.motif.str.split("-").str[0].isin(metals.value)]
    _d["site"] = _d.motif.map(pretty) + _d.state.map({"H_M": "", "H_S": ", H on S"})
    _d["metal"] = _d.motif.str.split("-").str[0]
    _d["near"] = _d.dG_H.round(2).abs() <= dft_window.value
    _band = alt.Chart(pd.DataFrame({"lo": [-dft_window.value], "hi": [dft_window.value]})).mark_rect(
        color="#2f6b4f", opacity=0.13).encode(x="lo:Q", x2="hi:Q")
    _pts = alt.Chart(_d).mark_circle(size=80).encode(
        x=alt.X("dG_H:Q", title="ΔG_H from DFT (eV)"), y=alt.Y("site:N", sort="x", title=None),
        color=alt.Color("metal:N", legend=None, scale=alt.Scale(domain=list(METAL_COLOURS), range=list(METAL_COLOURS.values()))),
        opacity=alt.condition("datum.near", alt.value(1.0), alt.value(0.45)),
        tooltip=["site", alt.Tooltip("dG_H:Q", format="+.2f")])
    dft_view = mo.vstack([
        mo.Html(f'<h3 class="sac-h3">DFT results</h3><p class="sac-intro"><b>{int(_d.near.sum())}</b> of {len(_d)} computed H '
                f'positions lie within ±{dft_window.value:.2f} eV of the optimum (shaded band).</p>'),
        (_band + _pts).properties(width=380, height=max(120, 15 * len(_d)))])
    return (dft_view,)


@app.cell
def _(mo):
    browse_reaction = mo.ui.dropdown(options=["all", "HER", "ORR", "CO2RR", "NRR"], value="all",
                                     label="reaction class (enzyme set only)")
    browse_site = mo.ui.radio(options=["both", "DV", "SV"], value="both", label="vacancy", inline=True)
    return browse_reaction, browse_site


@app.cell
def _(browse_reaction, browse_site, mo, rep):
    _r = rep[rep.preferred_site.isin(["DV", "SV"] if browse_site.value == "both" else [browse_site.value])]
    if browse_reaction.value != "all":
        _r = _r[_r.reaction.fillna("").str.contains(browse_reaction.value)]
    _top = (_r.groupby("graphene_motif").agg(entries=("pdb_id", "nunique"), sites=("site_id", "size"),
                                             example=("pdb_id", "first"))
            .sort_values("entries", ascending=False).reset_index())
    browse_view = mo.vstack([
        mo.Html('<h2 class="sac-h2">Browse the sites</h2><p class="sac-intro">Unique sites of the chosen metals, grouped by the '
                'graphene site they map to (DV double vacancy, SV single vacancy).</p>'),
        mo.hstack([browse_reaction, browse_site], justify="start"),
        mo.ui.table(_top, selection=None, page_size=10,
                    label=f"{len(_r):,} unique sites in {len(_top):,} graphene motifs (metal, donors, vacancy, axial ligand)")])
    return (browse_view,)


@app.cell
def _(REL, duckdb, mo, os):
    # every release table as a view, so a query can just say FROM sites
    db = duckdb.connect()
    TABLES = sorted(f[:-8] for f in os.listdir(REL) if f.endswith(".parquet"))
    for _t in TABLES:
        db.execute(f"CREATE VIEW {_t} AS SELECT * FROM read_parquet('{os.path.join(REL, _t + '.parquet')}')")
    sql_text = mo.ui.text_area(
        value="SELECT metal, COUNT(*) AS unique_sites, COUNT(DISTINCT pdb_id) AS entries,\n"
              "       ROUND(AVG(cn), 2) AS mean_coordination\n"
              "FROM sites\n"
              "WHERE is_representative AND metal IN ('V', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Mo', 'W')\n"
              "GROUP BY metal ORDER BY unique_sites DESC",
        rows=6, full_width=True, label="DuckDB SQL; tables: " + ", ".join(TABLES))
    return db, sql_text


@app.cell
def _(db, mo, sql_text):
    try:
        sql_view = mo.ui.table(db.execute(sql_text.value).df(), selection=None, page_size=10)
    except Exception as _e:              # a typo in the query: say what DuckDB says
        sql_view = mo.callout(mo.md(f"**The query did not run:** {_e}"), kind="danger")
    return (sql_view,)


@app.cell
def _(mo, os):
    # ask in plain language: Claude writes the SQL and runs it read-only on the release tables (db/ask.py); the key
    # comes from notebooks/.env (a line ANTHROPIC_API_KEY=..., git-ignored), which marimo loads at start
    import sys as _sys
    _root = os.path.dirname(str(mo.notebook_dir()))
    if _root not in _sys.path:
        _sys.path.insert(0, _root)
    from db.ask import EXAMPLES as _EXAMPLES, Agent as _Agent, has_credentials as _has_credentials
    ask_chat = mo.ui.chat(_Agent().chat, prompts=_EXAMPLES, max_height=560)
    _status = (mo.Html('<p class="sac-note">Claude Opus 5.5 writes each query; only reading SELECTs run, on the '
                       'release tables. Every query and its first rows stay visible above the answer.</p>')
               if _has_credentials() else
               mo.callout(mo.md("**No Anthropic API key yet.** Add a line `ANTHROPIC_API_KEY=...` to `notebooks/.env` "
                                "(git-ignored) and restart the explorer; until then, write the SQL yourself below."),
                          kind="warn"))
    ask_view = mo.vstack([
        mo.Html('<h2 class="sac-h2">Ask the database</h2><p class="sac-intro">Ask in plain language, in English or '
                'German. Claude writes the SQL, runs it on the release database and answers from the result.</p>'),
        _status, ask_chat], gap=0.6)
    return ask_chat, ask_view


@app.cell
def _(mo):
    support = mo.ui.radio(options={"MoS₂, metal on three S": "MoS2", "MoS₂, metal in an S vacancy": "MoS2vS",
                                   "silica, siloxide O": "SiO2"},
                          value="MoS₂, metal on three S", label="support", inline=True)
    support_set = mo.ui.dropdown(options={"enzymes (HER, ORR, CO2RR, NRR)": "enzyme", "whole PDB": "complete"},
                                 value="enzymes (HER, ORR, CO2RR, NRR)", label="sites")
    fits_only = mo.ui.switch(value=False, label="only where the geometry fits the support")
    support_controls = mo.vstack([
        mo.Html('<h2 class="sac-h2">MoS₂ and silica</h2><p class="sac-intro">Here the support supplies the donors: three lattice S, '
                '3.16 Å apart, on MoS₂ (or three Mo for a metal in an S vacancy); one to three siloxide O on silica. '
                'A motif reads <b>face(support)+ax:donors</b>, where the axial donors would have to come from a thiol '
                'or a silane linker.</p>'),
        mo.hstack([support, support_set, fits_only], justify="start", wrap=True)])
    return fits_only, support, support_controls, support_set


@app.cell
def _(METAL_COLOURS, alt, fits_only, metals, mo, pd, support, support_motifs, support_set, xtb_supports):
    if support_motifs.empty:
        _tab = mo.callout(mo.md("Rebuild the database (`python db/build.py`) for the support tables."), kind="warn")
    else:
        _m = support_motifs[(support_motifs.support == support.value) & (support_motifs.set_name == support_set.value)
                            & support_motifs.metal.isin(metals.value)]
        if fits_only.value:
            _m = _m[(_m.n_sites_hq >= 3) & (_m.frac_fit >= 0.5)]
        _m = _m.sort_values(["n_clusters30_rep", "n_entries_rep"], ascending=False)
        _tab = mo.ui.table(_m[["motif", "n_entries_rep", "n_clusters30_rep", "n_sites_hq", "median_xx", "frac_fit",
                               "graphene_motif_top", "reactions", "example_descriptions"]].rename(columns={
            "n_entries_rep": "entries", "n_clusters30_rep": "families (30 %)", "n_sites_hq": "X-ray ≤ 2 Å",
            "median_xx": "face spacing (Å)", "frac_fit": "share that fits", "graphene_motif_top": "on graphene",
            "example_descriptions": "examples"}).round(2).reset_index(drop=True), selection=None, page_size=8,
            label=f"{len(_m)} motifs")
    _sup = "MoS2" if support.value.startswith("MoS2") else "SiO2"
    _x = (xtb_supports[(xtb_supports.support == _sup) & (xtb_supports.state == "H_M") & xtb_supports.metal.isin(metals.value)]
          if not xtb_supports.empty else pd.DataFrame())
    if _x.empty:
        _chart = mo.md("")
    else:
        _long = pd.concat([_x.assign(kind="xTB", value=_x.dG_H),
                           _x.dropna(subset=["dG_H_offset"]).assign(kind="xTB + graphene offset",
                                                                    value=lambda d: d.dG_H_offset)])
        _band = alt.Chart(pd.DataFrame({"lo": [-0.2], "hi": [0.2]})).mark_rect(color="#2f6b4f", opacity=0.13).encode(
            x="lo:Q", x2="hi:Q")
        _rule = alt.Chart(_x).mark_rule(color="#c9d2d6").encode(x="dG_H:Q", x2="dG_H_offset:Q", y=alt.Y("metal:N", sort="x"))
        _pts = alt.Chart(_long).mark_point(size=110, strokeWidth=2).encode(
            x=alt.X("value:Q", title="ΔG_H, H on the metal (eV)"), y=alt.Y("metal:N", title=None, sort="x"),
            color=alt.Color("metal:N", legend=None, scale=alt.Scale(domain=list(METAL_COLOURS), range=list(METAL_COLOURS.values()))),
            fill=alt.condition("datum.kind == 'xTB'", alt.value("transparent"),
                               alt.Color("metal:N", legend=None, scale=alt.Scale(domain=list(METAL_COLOURS),
                                                                                 range=list(METAL_COLOURS.values())))),
            shape=alt.Shape("kind:N", title=None, legend=alt.Legend(orient="bottom"),
                            scale=alt.Scale(domain=["xTB", "xTB + graphene offset"], range=["circle", "diamond"])),
            tooltip=["metal", "kind", alt.Tooltip("value:Q", format="+.2f"), "site_bare"])
        _chart = mo.vstack([
            mo.Html(f"<h3 class='sac-h3'>Hydrogen on the bare {'MoS₂' if _sup == 'MoS2' else 'silica'} sites (GFN2-xTB)</h3>"
                    '<p class="sac-intro">Open circles: raw xTB; filled diamonds: plus the per-metal offset fitted on graphene against '
                    'DFT (none for Mn, Mo and W). Shaded: ±0.2 eV.</p>'),
            (_band + _rule + _pts).properties(width=460, height=30 * _x.metal.nunique() + 40)])
    supports_view = mo.vstack([_tab, _chart], gap=1.2)
    return (supports_view,)


@app.cell
def _(mo):
    REACTIONS = {"O₂ → H₂O (ORR, 4 e⁻)": ("U_ORR", "max", 1.23), "O₂ → H₂O₂ (2 e⁻)": ("U_H2O2", "max", 0.70),
                 "H₂O → O₂ (OER)": ("U_OER", "min", 1.23), "CO₂ → CO": ("U_CO", "max", -0.10),
                 "CO₂ → HCOOH": ("U_HCOOH", "max", -0.12), "N₂ → NH₃": ("U_NRR", "max", 0.06)}
    rx_reaction = mo.ui.dropdown(options=list(REACTIONS), value="O₂ → H₂O (ORR, 4 e⁻)", label="reaction")
    rx_support = mo.ui.multiselect(options=["graphene", "MoS2", "SiO2"], value=["graphene", "MoS2", "SiO2"],
                                   label="supports")
    beats_her = mo.ui.switch(value=False, label="only sites that beat H₂ evolution (CO₂, N₂)")
    rx_controls = mo.vstack([
        mo.Html('<h2 class="sac-h2">The same sites for O₂, CO₂ and N₂</h2><p class="sac-intro">Limiting potential (V vs RHE, pH 0, '
                'computational hydrogen electrode, GFN2-xTB): the most positive potential at which every step of a '
                'reduction runs downhill; for O₂ evolution the least positive. Read the values as a ranking: against the '
                'published PBE-D3 values of Cho et al. (Front. Chem. 2022, 18 sites) xTB orders the sites like DFT '
                '(Spearman 0.66) but is about 0.7 V off, and a linear correction does not help, so all values are raw '
                'xTB; CO₂ and N₂ have no DFT reference.</p>'),
        mo.hstack([rx_reaction, rx_support, beats_her], justify="start", wrap=True)])
    return REACTIONS, beats_her, rx_controls, rx_reaction, rx_support


@app.cell
def _(REACTIONS, SUPPORT_COLOURS, alt, beats_her, metals, mo, pd, rx, rx_reaction, rx_support, site_label):
    _col, _sense, _eq = REACTIONS[rx_reaction.value]
    _cal = _col + "_cal" if _col in ("U_ORR", "U_H2O2", "U_OER") else None
    _use = _col                         # raw xTB: the O2 line does not beat raw xTB once a site is left out
    _calibrated = False
    _chip = '<span class="sac-chip raw">raw xTB' + (": a ranking, not volts" if _cal else "") + "</span>"
    _scale = alt.Scale(domain=list(SUPPORT_COLOURS), range=list(SUPPORT_COLOURS.values()))
    if rx.empty or _use not in rx:
        reactions_view = mo.callout(mo.md("No values for this reaction yet: the screen is still running. They appear "
                                          "after `python screen/xtb_reactions.py table` and a reload."), kind="info")
    else:
        _all = rx[~rx.calibration_only.astype(bool)]
        _d = _all[_all.kind.isin(rx_support.value) & _all.metal.isin(metals.value)].dropna(subset=[_use]).copy()
        if beats_her.value and _col in ("U_CO", "U_HCOOH", "U_NRR"):
            _d = _d[_d[_use] > -_d.dG_H.abs()]            # needs less overpotential than H2 evolution
        _d["site_label"] = [site_label(k, s, m) for k, s, m in zip(_d.kind, _d.site, _d.metal)]
        _d["U"] = _d[_use]
        _d = _d.sort_values("U", ascending=(_sense == "min"))
        _dom = [min(float(_d.U.min()), _eq) - 0.15, max(float(_d.U.max()), _eq) + 0.15] if len(_d) else [_eq - 1, _eq + 1]
        _rule = alt.Chart(pd.DataFrame({"eq": [_eq]})).mark_rule(color="#2f6b4f", strokeDash=[4, 3], strokeWidth=1.5).encode(
            x=alt.X("eq:Q", scale=alt.Scale(domain=_dom)))
        _pts = alt.Chart(_d).mark_circle(size=80, opacity=0.95).encode(
            x=alt.X("U:Q", title="limiting potential (V vs RHE)", scale=alt.Scale(domain=_dom)),
            y=alt.Y("site_label:N", sort=list(_d.site_label), title=None),
            color=alt.Color("kind:N", title=None, legend=alt.Legend(orient="bottom"), scale=_scale),
            tooltip=["site_label", "kind", alt.Tooltip("U:Q", format="+.2f", title="U_L (V)"),
                     alt.Tooltip("dG_H:Q", format="+.2f", title="ΔG_H (eV)"),
                     alt.Tooltip("dG_H_source:N", title="ΔG_H from")])
        _n = int(_all[_use].notna().sum())
        _parts = [mo.Html(f'<p class="sac-intro">{_chip} &nbsp;<b>{_n}</b> of {len(_all)} sites have values for '
                          f'{rx_reaction.value}. Dashed line: equilibrium potential {_eq:+.2f} V. Best first.</p>'),
                  (_rule + _pts).properties(width=520, height=max(120, 15 * len(_d)))]
        if _col in ("U_ORR", "U_OER", "U_H2O2"):        # the activity volcano against *OH binding
            _x = "dG_OH_cal" if _calibrated and "dG_OH_cal" in _d else "dG_OH"
            _v = alt.Chart(_d.dropna(subset=[_x])).mark_circle(size=80).encode(
                x=alt.X(f"{_x}:Q", title="ΔG(*OH) (eV)"), y=alt.Y("U:Q", title="limiting potential (V)"),
                color=alt.Color("kind:N", legend=None, scale=_scale),
                tooltip=["site_label", alt.Tooltip("U:Q", format="+.2f"), alt.Tooltip(f"{_x}:Q", format="+.2f")])
            _parts += [mo.Html("<h3 class='sac-h3'>Volcano: activity against *OH binding</h3>"), _v.properties(width=440, height=300)]
        reactions_view = mo.vstack(_parts)
    return (reactions_view,)


@app.cell
def _(alt, mo, np, pd, rx_cal):
    if rx_cal.empty:
        rx_cal_view = mo.md("")
    else:
        _fits = []
        for _a, _q in rx_cal.groupby("adsorbate"):
            if len(_q) >= 2:
                _k, _c = np.polyfit(_q.dE_xtb, _q.dE_dft, 1)
                _fits.append(f"*{_a}: {len(_q)} sites, mean error {(_q.dE_xtb - _q.dE_dft).abs().mean():.2f} eV raw, "
                             f"{(_k * _q.dE_xtb + _c - _q.dE_dft).abs().mean():.2f} eV after the line")
        _lo = float(min(rx_cal.dE_xtb.min(), rx_cal.dE_dft.min())) - 0.3
        _hi = float(max(rx_cal.dE_xtb.max(), rx_cal.dE_dft.max())) + 0.3
        _diag = alt.Chart(pd.DataFrame({"x": [_lo, _hi], "y": [_lo, _hi]})).mark_line(
            color="#8f9ca2", strokeDash=[4, 4]).encode(x="x:Q", y="y:Q")
        _pts = alt.Chart(rx_cal).mark_circle(size=80).encode(
            x=alt.X("dE_xtb:Q", title="xTB reaction energy (eV)", scale=alt.Scale(domain=[_lo, _hi])),
            y=alt.Y("dE_dft:Q", title="DFT, Cho et al. 2022 (eV)", scale=alt.Scale(domain=[_lo, _hi])),
            color=alt.Color("adsorbate:N", title=None, legend=alt.Legend(orient="bottom"),
                            scale=alt.Scale(domain=["OH", "O", "OOH"], range=["#2c6bd3", "#c4682f", "#3f9a6a"])),
            tooltip=["site", "adsorbate", alt.Tooltip("dE_xtb:Q", format="+.2f"), alt.Tooltip("dE_dft:Q", format="+.2f")])
        rx_cal_view = mo.vstack([
            mo.Html('<h3 class="sac-h3">O₂ check: xTB against published DFT</h3><p class="sac-intro">' + "; ".join(_fits) + "</p>"),
            (_diag + _pts).properties(width=400, height=340)])
    return (rx_cal_view,)


if __name__ == "__main__":
    app.run()
