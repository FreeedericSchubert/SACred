import marimo

__generated_with = "0.25.2"
app = marimo.App(width="full", app_title="SAC synthesis check", css_file="sac_explorer.css",
                 html_head_file="sac_explorer_head.html")


@app.cell
def _():
    import glob
    import json
    import os

    import altair as alt
    import marimo as mo
    import numpy as np
    import pandas as pd

    # the website's metal colours (site/template.html), as in the explorer
    METAL_COLOURS = {"Fe": "#c4682f", "Co": "#3f6fc0", "Ni": "#3f9a6a", "Cu": "#1f8f8a", "Zn": "#7a8796",
                     "Mo": "#8a68b0", "W": "#5d5a93", "Mn": "#b0427a", "V": "#b28a12", "Cr": "#a0522d",
                     "Ru": "#4b6b8a", "Pt": "#6b6f73"}
    ENZYME_METALS = ["V", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W"]

    @alt.theme.register("sac_synthesis", enable=True)
    def _sac_theme():
        _ink, _muted, _grid = "#5d6b70", "#7d8b91", "#e1e7ea"
        return {"config": {
            "background": "transparent", "font": "Source Sans 3", "view": {"stroke": None},
            "axis": {"labelColor": _muted, "titleColor": _ink, "gridColor": _grid, "domainColor": "#c9d2d6",
                     "tickColor": "#c9d2d6", "labelFontSize": 12, "titleFontSize": 12, "titleFontWeight": 600},
            "legend": {"labelColor": _ink, "titleColor": _ink, "labelFontSize": 12, "symbolSize": 90},
            "title": {"font": "Newsreader", "fontSize": 16, "color": _ink}}}
    return ENZYME_METALS, METAL_COLOURS, alt, glob, json, mo, np, os, pd


@app.cell
def _(glob, json, mo, os, pd):
    # data: screen/geometry_mismatch.py, screen/xtb_stability.py, screen/stability_descriptors.py, dft/analyse.py
    ROOT = os.path.dirname(str(mo.notebook_dir()))
    DATA = os.path.join(ROOT, "data")

    def _read(path, reader):
        return reader(path) if os.path.exists(path) else pd.DataFrame()

    geo = _read(os.path.join(DATA, "complete", "geometry_mismatch_sites.parquet"), pd.read_parquet)
    geo_models = _read(os.path.join(DATA, "geometry_mismatch.csv"), pd.read_csv)
    dft_stab = _read(os.path.join(DATA, "stability_graphene.csv"), pd.read_csv)
    dft_stab = dft_stab.dropna(subset=["E_f"]) if len(dft_stab) else dft_stab
    stab = _read(os.path.join(DATA, "xtb_stability.csv"), pd.read_csv)
    neb = pd.DataFrame([json.load(open(f)) for f in sorted(glob.glob(os.path.join(DATA, "screen", "xtb_neb", "*.json")))])
    _sites_done = len(glob.glob(os.path.join(DATA, "screen", "xtb_sites", "*.json")))
    progress = {"empty_sites": _sites_done}
    _desc = os.path.join(DATA, "xtb_stability_descriptors.md")
    descriptors_md = open(_desc).read() if os.path.exists(_desc) else ""
    three_way = _read(os.path.join(DATA, "geometry_three_way.csv"), pd.read_csv)
    return descriptors_md, dft_stab, geo, geo_models, neb, progress, stab, three_way


@app.cell
def _(dft_stab, geo, mo):
    # header: the question and the numbers that answer it
    _kept = geo[geo.lost == 0]
    _facts = [
        (f"{len(geo):,}", "enzyme sites (X-ray ≤ 2.0 Å, mononuclear) compared with their relaxed model on graphene", ""),
        (f"{geo.motif.nunique()}", "graphene models (DFT where finished, else GFN2-xTB)", ""),
        (f"{(geo.lost > 0).mean():.0%}", "of the sites: the model lets go of a donor (all in single vacancies)", "bad"),
        (f"{_kept[(_kept.site == 'DV') & (_kept.holder == 'side chains')].angle_rms.median():.0f}°",
         "angles off in double vacancies (side-chain sites, RMS)", ""),
        (f"{_kept[(_kept.site == 'SV') & (_kept.holder == 'side chains')].dh.median():+.1f} Å",
         "metal further out of the donor plane in single vacancies", ""),
        (f"{int((dft_stab.E_f < 0).sum())} of {len(dft_stab)}", "sites hold their metal better than the bulk metal (DFT E_f < 0)",
         "bad" if len(dft_stab) and (dft_stab.E_f < 0).mean() < 0.5 else "good"),
    ]
    _tiles = "".join(f'<div class="sac-fact {c}"><b>{v}</b><span>{t}</span></div>' for v, t, c in _facts)
    mo.Html(f"""
    <header class="sac-hero">
      <div class="sac-eyebrow">Hackathon project · draft results · synthesis check</div>
      <h1>Can carbon hold the enzyme's metal?</h1>
      <p class="sac-lede">An enzyme holds its metal with flexible side chains; a graphene vacancy is rigid. Here the same
      first shell is compared in the enzyme (PDB) and in its relaxed single-atom model, and the model's stability against
      clustering into metal particles is tested. Filters in the sidebar apply to every tab.</p>
      <div class="sac-facts">{_tiles}</div>
    </header>""")
    return


@app.cell
def _(ENZYME_METALS, mo):
    metals = mo.ui.multiselect(options=ENZYME_METALS, value=ENZYME_METALS, label="metals")
    vacancy = mo.ui.radio(options=["both", "DV", "SV"], value="both", label="vacancy", inline=True)
    holder = mo.ui.radio(options=["all", "side chains", "cofactor"], value="all", label="enzyme donors from", inline=True)
    keep_lost = mo.ui.switch(value=False, label="include models that lose a donor")
    return holder, keep_lost, metals, vacancy


@app.cell
def _(METAL_COLOURS, holder, keep_lost, metals, mo, vacancy):
    _legend = "".join(f'<span><i style="background:{c}"></i>{m}</span>' for m, c in METAL_COLOURS.items()
                      if m in metals.value)
    mo.sidebar([
        mo.Html('<p class="sac-side-title">Filters</p>'),
        metals, vacancy, holder, keep_lost,
        mo.Html('<p class="sac-note">DV: double vacancy, four donors in the sheet. SV: single vacancy, three donors. '
                '"Cofactor": heme N, Fe–S cluster S or dithiolene S hold the metal; "side chains": amino acids, '
                'backbone and water.</p>'),
        mo.Html(f'<p class="sac-side-title" style="margin-top:10px">Metal colours</p><div class="sac-legend">{_legend}</div>'),
        mo.Html('<p class="sac-note">Model minus enzyme throughout. All values are drafts that need expert checking.</p>'),
    ])
    return


@app.cell
def _(geo, geo_models, holder, keep_lost, metals, vacancy):
    # the filtered enzyme sites and models
    g = geo[geo.metal.isin(metals.value)]
    if vacancy.value != "both":
        g = g[g.site == vacancy.value]
    if holder.value != "all":
        g = g[g.holder == holder.value]
    if not keep_lost.value:
        g = g[g.lost == 0]
    gm = geo_models[geo_models.motif.isin(g.motif.unique())].copy()
    gm["metal"] = gm.motif.str.split("-").str[0]
    return g, gm


@app.cell
def _(METAL_COLOURS, alt, pd):
    def box_chart(df, value, by, x_title, domain=None):
        """Box plot from statistics computed here (35,000 sites are too many rows to hand to the browser):
        median tick, 25-75 % bar, 5-95 % line, per metal (and the facet column `by`)."""
        keys = ["metal"] + ([by] if by else [])
        q = df.dropna(subset=[value]).groupby(keys)[value].quantile([0.05, 0.25, 0.5, 0.75, 0.95]).unstack()
        q.columns = ["p05", "q1", "median", "q3", "p95"]
        q = q.reset_index().merge(df.groupby(keys).size().rename("n").reset_index(), on=keys)
        scale = alt.Scale(domain=list(METAL_COLOURS), range=list(METAL_COLOURS.values()))
        xs = alt.Scale(domain=domain, clamp=True) if domain else alt.Undefined
        y = alt.Y("metal:N", title=None, sort=list(METAL_COLOURS))
        tip = ["metal", "n"] + ([by] if by else []) + [alt.Tooltip(c + ":Q", format="+.2f") for c in ("q1", "median", "q3")]
        base = alt.Chart()
        layers = [base.mark_rule(color="#2f6b4f", strokeDash=[4, 3]).encode(x=alt.datum(0)),
                  base.mark_rule(strokeWidth=1.5).encode(x=alt.X("p05:Q", title=x_title, scale=xs), x2="p95:Q", y=y,
                                                         color=alt.Color("metal:N", scale=scale, legend=None)),
                  base.mark_bar(size=12, opacity=0.55).encode(x="q1:Q", x2="q3:Q", y=y, tooltip=tip,
                                                              color=alt.Color("metal:N", scale=scale, legend=None)),
                  base.mark_tick(thickness=3, size=14, color="#1c2427").encode(x="median:Q", y=y, tooltip=tip)]
        chart = alt.layer(*layers, data=q).properties(width=200, height=24 * q.metal.nunique() + 20)
        return chart.facet(column=alt.Column(by + ":N", title=None, header=alt.Header(labelFontSize=13))) if by else chart
    return (box_chart,)


@app.cell
def _(box_chart, g, mo):
    # bond lengths per donor element
    _long = g.melt(id_vars=["metal", "site"], value_vars=["dd_N", "dd_S", "dd_O"], var_name="donor",
                   value_name="dd").dropna()
    _long["donor"] = _long.donor.str[-1]
    bonds_view = mo.vstack([
        mo.Html('<h3 class="sac-h3">Bond lengths</h3><p class="sac-intro">Graphene shortens most metal–donor bonds by '
                '0.1–0.4 Å against the enzyme: the vacancy edge is tighter than a protein side chain. Zn is the exception, '
                'its bonds get longer. Bar 25–75 %, line 5–95 %, tick the median.</p>'),
        box_chart(_long, "dd", "donor", "metal–donor bond, model − enzyme (Å)", [-0.8, 1.0])
        if len(_long) else mo.md("_No sites in this selection._")])
    return (bonds_view,)


@app.cell
def _(box_chart, g, mo):
    # angles and the metal's position relative to its donors
    _ang = box_chart(g, "angle_rms", "site", "donor–metal–donor angles, RMS change (°)")
    _dh = box_chart(g, "dh", "site", "metal out of the donor plane, model − enzyme (Å)")
    angles_view = mo.vstack([
        mo.Html('<h3 class="sac-h3">Angles and planarity</h3><p class="sac-intro">In a double vacancy the support '
                'flattens the site and pulls the metal into the plane; in a single vacancy it pushes the metal out of '
                'the plane of its three donors. Heme-like cofactors fit the double vacancy almost unchanged (about 3°).</p>'),
        mo.vstack([_ang, _dh], gap=1) if len(g) else mo.md("_No sites in this selection._")])
    return (angles_view,)


@app.cell
def _(METAL_COLOURS, alt, gm, mo, pd):
    # per model: metal out of the donor plane, enzyme against model
    _scale = alt.Scale(domain=list(METAL_COLOURS), range=list(METAL_COLOURS.values()))
    _diag = alt.Chart(pd.DataFrame({"x": [0, 3.5], "y": [0, 3.5]})).mark_line(color="#8f9ca2", strokeDash=[4, 4]).encode(
        x="x:Q", y="y:Q")
    _pts = alt.Chart(gm).mark_circle(opacity=0.85, stroke="white", strokeWidth=0.5).encode(
        x=alt.X("h_enzyme:Q", title="metal out of the donor plane, enzyme median (Å)", scale=alt.Scale(domain=[0, 3.5])),
        y=alt.Y("h_model:Q", title="… in the relaxed model (Å)", scale=alt.Scale(domain=[0, 3.5])),
        size=alt.Size("enzyme_sites:Q", title="enzyme sites", scale=alt.Scale(range=[30, 600])),
        color=alt.Color("metal:N", scale=_scale, legend=None),
        shape=alt.Shape("site:N", title=None),
        tooltip=["motif", "model", "enzyme_sites", "donors_lost", alt.Tooltip("angle_rms:Q", format=".0f", title="angles RMS (°)"),
                 alt.Tooltip("h_enzyme:Q", format=".2f"), alt.Tooltip("h_model:Q", format=".2f")])
    models_view = mo.vstack([
        mo.Html('<h3 class="sac-h3">Model by model</h3><p class="sac-intro">Each point is one graphene model, sized by '
                'the enzyme sites behind it. On the diagonal the support keeps the enzyme\'s geometry; above it the metal '
                'sits further out of its donor plane than in the enzyme.</p>'),
        (_diag + _pts).properties(width=480, height=380)])
    return (models_view,)


@app.cell
def _(METAL_COLOURS, alt, mo, pd, three_way):
    # enzyme -> free complex -> support (screen/xtb_free_complexes.py)
    if three_way.empty:
        three_way_view = mo.callout(mo.md("**Running:** the free complexes (screen/xtb_free_complexes.py) fill this part."),
                                    kind="info")
    else:
        _t = three_way.melt(id_vars=["motif", "metal", "site"], value_vars=["h_enzyme", "h_free", "h_support"],
                            var_name="where", value_name="h")
        _t["where"] = _t["where"].map({"h_enzyme": "1 enzyme", "h_free": "2 free complex", "h_support": "3 support"})
        _scale = alt.Scale(domain=list(METAL_COLOURS), range=list(METAL_COLOURS.values()))
        _line = alt.Chart(_t).mark_line(opacity=0.35).encode(
            x=alt.X("where:N", title=None), y=alt.Y("median(h):Q", title="metal out of the donor plane (Å)"),
            detail="motif:N", color=alt.Color("metal:N", scale=_scale, legend=None),
            tooltip=["motif", "where", alt.Tooltip("median(h):Q", format=".2f")])
        _pts = _line.mark_circle(size=50)
        _share = (three_way.angles_enzyme_free / (three_way.angles_enzyme_free + three_way.angles_free_support)).median()
        _ang = pd.DataFrame({"step": ["enzyme → free complex", "free complex → support"],
                             "deg": [three_way.angles_enzyme_free.median(), three_way.angles_free_support.median()]})
        _bar = alt.Chart(_ang).mark_bar(size=22, color="#59666b").encode(
            x=alt.X("deg:Q", title="angles changed, median RMS (°)"), y=alt.Y("step:N", title=None),
            tooltip=["step", alt.Tooltip("deg:Q", format=".0f")])
        three_way_view = mo.vstack([
            mo.Html('<h3 class="sac-h3">Who changes the geometry: protein or support?</h3><p class="sac-intro">The same '
                    f'first shells as free complexes (M(II) with model ligands, relaxed in vacuum): {len(three_way)} '
                    f'sites of {three_way.motif.nunique()} cores. The protein accounts for a median '
                    f'{100 * _share:.0f} % of the angle changes from enzyme over free complex to support, the support '
                    'for the rest.</p>'),
            mo.hstack([_bar.properties(width=320, height=90), (_line + _pts).properties(width=360, height=260)],
                      justify="start", gap=2)])
    return (three_way_view,)


@app.cell
def _(gm, mo):
    _cols = {"motif": "model", "site": "vacancy", "model": "method", "enzyme_sites": "enzyme sites",
             "entries": "PDB entries", "donors_lost": "donors lost", "angle_rms": "angles RMS (°)",
             "h_enzyme": "out of plane, enzyme (Å)", "h_model": "out of plane, model (Å)",
             "tau4_enzyme": "τ4 enzyme", "tau4_model": "τ4 model", "share_cofactor": "share cofactor"}
    _t = (gm[[c for c in _cols if c in gm]].rename(columns=_cols).sort_values("enzyme sites", ascending=False)
          .round(2).reset_index(drop=True))
    table_view = mo.vstack([mo.Html('<h3 class="sac-h3">All models</h3>'),
                            mo.ui.table(_t, selection=None, page_size=10)])
    return (table_view,)


@app.cell
def _(METAL_COLOURS, alt, dft_stab, mo, pd):
    # DFT formation energies of the metal in its site
    _d = dft_stab.copy()
    _d["metal"] = _d.motif.str.split("-").str[0]
    _scale = alt.Scale(domain=list(METAL_COLOURS), range=list(METAL_COLOURS.values()))
    _zero = alt.Chart(pd.DataFrame({"x": [0]})).mark_rule(color="#2f6b4f", strokeDash=[4, 3]).encode(x="x:Q")
    _pts = alt.Chart(_d).mark_circle(size=140).encode(
        x=alt.X("E_f:Q", title="E_f = E(metal in site) − E(empty site) − E(bulk metal per atom) (eV)"),
        y=alt.Y("motif:N", title=None, sort="x"),
        color=alt.Color("metal:N", scale=_scale, legend=None),
        tooltip=["motif", alt.Tooltip("E_f:Q", format="+.2f"), alt.Tooltip("U_diss:Q", format="+.2f", title="U_diss (V)")])
    dft_view = mo.vstack([
        mo.Html('<h3 class="sac-h3">DFT: single atom or particle?</h3><p class="sac-intro">Left of the dashed line the '
                'metal prefers its site to the bulk metal; to the right it would rather cluster. Only Zn–N₃ is on the left. '
                'N donors hold better than S (Cu: N₃ < N₂S < NS₂), and across metals the order follows the bulk cohesive '
                'energy (Zn < Cu < Ni < Fe for NS₂). Against the bulk metal this is the strictest test.</p>'),
        (_zero + _pts).properties(width=520, height=max(160, 24 * len(_d)))])
    return (dft_view,)


@app.cell
def _(METAL_COLOURS, alt, geo_models, mo, pd, progress, stab):
    # GFN2-xTB: insertion and escape energy of every site, and the link to the geometry
    if stab.empty:
        xtb_view = mo.callout(mo.md(f"**Running:** the 96 empty donor arrangements ({progress['empty_sites']}/96 done); "
                                    "then `python screen/xtb_stability.py table` fills this part."), kind="info")
    else:
        _s = stab[stab.converged == True].copy()   # noqa: E712
        _s["dE"] = _s.E_ins - _s.groupby("metal").E_ins.transform("mean")
        _scale = alt.Scale(domain=list(METAL_COLOURS), range=list(METAL_COLOURS.values()))
        _strip = alt.Chart(_s).mark_tick(thickness=2, size=16, opacity=0.6).encode(
            x=alt.X("dE:Q", title="E_ins relative to the metal's mean (eV; lower holds more strongly)"),
            y=alt.Y("metal:N", title=None, sort=list(METAL_COLOURS)),
            color=alt.Color("metal:N", scale=_scale, legend=None),
            tooltip=["motif", alt.Tooltip("E_ins:Q", format="+.2f"), alt.Tooltip("dE:Q", format="+.2f")])
        _parts = [mo.Html('<h3 class="sac-h3">GFN2-xTB: every donor arrangement</h3><p class="sac-intro">E_ins = '
                          'E(metal in the site) − E(empty site) for every scanned and library site; within one metal it '
                          'ranks the arrangements directly.</p>'), _strip.properties(width=520, height=260)]
        _m = _s.merge(geo_models[["motif", "angle_rms", "h_model", "h_enzyme", "donors_lost"]], on="motif")
        if len(_m):
            _sc = alt.Chart(_m).mark_circle(size=70, opacity=0.8).encode(
                x=alt.X("angle_rms:Q", title="angles changed against the enzyme, RMS (°)"),
                y=alt.Y("dE:Q", title="E_ins relative to the metal's mean (eV)"),
                color=alt.Color("metal:N", scale=_scale, legend=None),
                tooltip=["motif", alt.Tooltip("angle_rms:Q", format=".0f"), alt.Tooltip("dE:Q", format="+.2f")])
            _parts += [mo.Html('<h3 class="sac-h3">Does a forced geometry cost stability?</h3>'),
                       _sc.properties(width=520, height=300)]
        if "E_escape" in _s and _s.E_escape.notna().any():
            _rules = alt.Chart(pd.DataFrame({"x": [3.5, 4.2], "label": ["1 h at 800 °C", "1 h at 1000 °C"]})).mark_rule(
                color="#b0427a", strokeDash=[4, 3]).encode(x="x:Q", tooltip=["label"])
            _esc = alt.Chart(_s.dropna(subset=["E_escape"])).mark_tick(thickness=2, size=16, opacity=0.6).encode(
                x=alt.X("E_escape:Q", title="escape energy: site → adatom on clean carbon (eV)"),
                y=alt.Y("metal:N", title=None, sort=list(METAL_COLOURS)),
                color=alt.Color("metal:N", scale=_scale, legend=None), tooltip=["motif", alt.Tooltip("E_escape:Q", format=".2f")])
            _parts += [mo.Html('<h3 class="sac-h3">Kinetics: can the metal leave during pyrolysis?</h3><p class="sac-intro">'
                               'The escape energy is a lower bound of the barrier out of the site. An hour of pyrolysis gets '
                               'an atom over about 3.5 eV at 800 °C and 4.2 eV at 1000 °C (dashed lines); sites to the left '
                               'empty during synthesis.</p>'), (_rules + _esc).properties(width=520, height=260)]
        xtb_view = mo.vstack(_parts)
    return (xtb_view,)


@app.cell
def _(descriptors_md, mo, neb):
    _parts = []
    if len(neb) and "barrier" in neb:
        _n = neb.dropna(subset=["barrier"])[["motif", "barrier", "reaction", "converged"]].round(2)
        _parts += [mo.Html('<h3 class="sac-h3">Escape barriers (NEB, pilot)</h3><p class="sac-intro">Climbing-image NEB '
                           'from the site to an intact hexagon 5 Å away; if the barrier sits at the escape energy, the '
                           'cheap number is enough.</p>'), mo.ui.table(_n, selection=None)]
    else:
        _parts.append(mo.callout(mo.md("**NEB pilot** (12 sites with DFT stability values) starts after the empty sites "
                                       "and fills this part when it finishes."), kind="info"))
    if descriptors_md:
        _parts += [mo.Html('<h3 class="sac-h3">Descriptor analysis</h3>'), mo.md(descriptors_md)]
    kinetics_view = mo.vstack(_parts, gap=1)
    return (kinetics_view,)


@app.cell
def _(mo):
    method_view = mo.md("""
    ### How this is computed

    - **Enzyme side:** representative X-ray sites at ≤ 2.0 Å of the nine enzyme metals, without a metal neighbour. Of
      each site only the donors that go into the vacancy are compared (the four of the double vacancy or the three of
      the single vacancy, chosen by the same rules as the mapping; reproduced for 372/372 DV and 399/400 SV sites).
    - **Model side:** the relaxed bare site, DFT where it finished, else GFN2-xTB. A donor counts as lost when the
      model's bond is more than 0.3 Å beyond the donor cutoff of the PDB analysis.
    - **Mismatch:** model minus enzyme; angles compared as sorted lists, so no donor-by-donor assignment is needed.
    - **Stability:** DFT E_f against the bulk metal (dft/analyse.py); xTB E_ins = E(metal in the site) − E(empty site),
      the empty sites relaxed once per arrangement; the escape energy adds the metal adsorbed on clean graphene.
    - **Reference:** the enzyme stands in for a freely coordinated complex. Harding's work on protein metal sites
      (Acta Cryst. D, 1999–2006) uses distances from small-molecule complexes (CSD) as the target, which well-refined
      protein sites follow; a free-complex calculation would make the three-way comparison explicit.
    - Code: `screen/geometry_mismatch.py`, `screen/xtb_stability.py`, `screen/stability_descriptors.py`.
    """)
    return (method_view,)


@app.cell
def _(angles_view, bonds_view, dft_view, kinetics_view, method_view, mo, models_view, table_view, three_way_view,
      xtb_view):
    mo.ui.tabs({
        "Geometry: enzyme → support": mo.vstack([bonds_view, angles_view, models_view, three_way_view, table_view],
                                                gap=1.5),
        "Stability and kinetics": mo.vstack([dft_view, xtb_view, kinetics_view], gap=1.5),
        "Method": method_view,
    })
    return


if __name__ == "__main__":
    app.run()
