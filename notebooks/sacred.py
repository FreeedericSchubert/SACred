import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium", app_title="SACred", css_file="sac_pitch.css",
                 html_head_file="sac_explorer_head.html")


@app.cell
def _():
    import base64
    import io
    import os

    import altair as alt
    import marimo as mo
    import pandas as pd

    @alt.theme.register("sac", enable=True)
    def _sac_theme():
        # mid-tone text so the charts read on the light and on the dark marimo theme
        _ink, _muted, _grid = "#5d6b70", "#7d8b91", "#e1e7ea"
        return {"config": {
            "background": "transparent", "font": "Source Sans 3", "view": {"stroke": None},
            "axis": {"labelColor": _muted, "titleColor": _ink, "gridColor": _grid, "domainColor": "#c9d2d6",
                     "tickColor": "#c9d2d6", "labelFontSize": 12, "titleFontSize": 12, "titleFontWeight": 600},
            "legend": {"labelColor": _muted, "titleColor": _ink, "labelFontSize": 12, "titleFontSize": 12},
            "title": {"font": "Newsreader", "fontSize": 16, "color": _ink}}}

    def pretty(motif):
        """Fe-NS2:SV -> Fe–NS₂; arrangement labels dropped."""
        return motif.split(":")[0].replace("-", "–", 1).translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))

    METALS = ["Fe", "Co", "Ni", "Cu", "Zn", "Mn", "Mo", "W", "V", "Cr", "Ru", "Pt"]
    COLORS = ["#c4682f", "#3f6fc0", "#3f9a6a", "#1f8f8a", "#7a8796", "#b0427a", "#8a68b0", "#5d5a93", "#b28a12",
              "#6b7f2a", "#a5523f", "#6f7c86"]
    metal_color = alt.Scale(domain=METALS, range=COLORS)
    return alt, base64, io, metal_color, mo, os, pd, pretty


@app.cell
def _(base64, io, mo, os, pd):
    # On GitHub Pages (WebAssembly) the files sit in notebooks/public/ next to the page (copied there by the Pages
    # workflow); locally they are read from data/, notebooks/assets/ and site/.
    _loc = mo.notebook_location()
    WEB = str(_loc).startswith("http")
    _root = os.path.dirname(str(mo.notebook_dir())) if not WEB else ""
    _dirs = [os.path.join(_root, d) for d in ("data", "notebooks/assets", "site")]

    def src(name):
        if WEB:
            return str(_loc / "public" / name)
        return next(os.path.join(d, name) for d in _dirs if os.path.exists(os.path.join(d, name)))

    def embed(name, mime):
        """URL on the web, data URI locally (so a static export carries the file too)."""
        if WEB:
            return "public/" + name
        return "data:%s;base64,%s" % (mime, base64.b64encode(open(src(name), "rb").read()).decode())

    def text(name):
        # fetched by hand on the web: GitHub Pages serves gzip, the browser unpacks it, and pandas would try again
        if WEB:
            from pyodide.http import open_url
            return io.StringIO(open_url(src(name)).read())
        return src(name)

    dgh = pd.read_csv(text("dgh_graphene.csv"))
    stab = pd.read_csv(text("stability_graphene.csv"))
    xtb = pd.read_csv(text("xtb_scan.csv"))
    offsets = pd.read_csv(text("xtb_offsets.csv"))
    supports = pd.read_csv(text("xtb_supports.csv"))
    reactions = pd.read_csv(text("xtb_reactions.csv"))
    literature = pd.read_json(text("literature_check_final.json"))
    gallery = pd.read_json(text("gallery.json"), typ="series").to_dict()
    assets = {"logo": embed("sacred_logo.webp", "image/webp"), "md": embed("md_pyrolysis_1500K.mp4", "video/mp4"),
              "zoom": src("zoom_viewer.html") if WEB else open(src("zoom_viewer.html")).read()}
    return WEB, assets, dgh, gallery, literature, offsets, reactions, stab, supports, xtb


@app.cell
def _(dgh, pd, pretty, stab):
    # one DFT value per motif: the H position closest to the optimum, with the metal's formation energy if computed
    _d = (dgh.dropna(subset=["dG_H"]).assign(a=lambda d: d.dG_H.abs()).sort_values("a").drop_duplicates("motif")
          .drop(columns=["E_f", "U_diss"], errors="ignore"))
    dft = _d.merge(stab[["motif", "E_f", "U_diss"]], on="motif", how="left")
    dft = dft.assign(site=dft.motif.map(pretty), metal=dft.motif.str.split("-").str[0],
                     H_on=dft.state.map({"H_M": "metal", "H_S": "S donor"}).fillna(dft.state),
                     stability=pd.Series(pd.NA, index=dft.index, dtype="object"))
    dft.loc[dft.E_f.notna(), "stability"] = "prefers metal particles"
    dft.loc[dft.E_f < 0, "stability"] = "holds its metal (E_f < 0)"
    dft.loc[(dft.E_f < 0) & (dft.U_diss > 0), "stability"] = "holds its metal, also in water (U_diss > 0)"
    dft["stability"] = dft.stability.fillna("not computed")
    n_usable = int((dft.dG_H.round(2).abs() <= 0.2).sum())
    return dft, n_usable


@app.cell
def _(assets, dft, mo, n_usable, xtb):
    _facts = [("291,977", "metal sites from the Protein Data Bank"),
              (f"{xtb.motif.nunique():,}", "single-atom sites in the xTB funnel"),
              (str(len(dft)), "sites refined with DFT"),
              (str(n_usable), "usable motifs for H₂ evolution (preliminary)")]
    _tiles = "".join(f'<div class="sac-fact{" good" if k == 3 else ""}"><b>{v}</b><span>{t}</span></div>'
                     for k, (v, t) in enumerate(_facts))
    mo.Html(f"""
    <header class="sac-hero">
      <div class="pt-brand"><img src="{assets['logo']}" alt="SACred logo"><div>
      <div class="sac-eyebrow">SACred · bioinspired single-atom catalysts from the PDB</div>
      <h1>Turning enzyme active sites into usable single-atom catalysts</h1></div></div>
      <p class="sac-lede">An agentic AI pipeline and open dataset: every metal site in the Protein Data Bank is reduced
      to its first coordination shell, rebuilt as a single-atom site in graphene, MoS₂ or silica, screened with
      GFN2-xTB and refined with DFT, for activity and for survival in synthesis. This notebook lets you browse the
      results. All values are hackathon drafts that need expert checking.</p>
      <div class="sac-facts">{_tiles}</div>
    </header>""")
    return


@app.cell
def _(WEB, assets, mo):
    _zoom = (mo.Html(f'<iframe src="{assets["zoom"]}" style="width:100%;height:820px;border:0"></iframe>') if WEB else
             mo.iframe('<!doctype html><html data-theme="light"><meta charset="utf-8"><meta name="viewport" '
                       'content="width=device-width, initial-scale=1">' + assets["zoom"], height="820px"))
    mo.vstack([
        mo.Html('<section class="pt-section"><div class="pt-kicker">1 · From enzyme to catalyst</div>'
                '<h2>Every metal site of the PDB as a template</h2><p class="pt-lede">Press Play: from the whole enzyme '
                'to its active site, to the first coordination shell, and to the same metal as a single-atom site on '
                'graphene, MoS₂ and silica.</p></section>'),
        _zoom], gap=0.4)
    return


@app.cell
def _(mo):
    mo.Html("""
    <section class="pt-section">
      <div class="pt-kicker">2 · The pipeline</div>
      <h2>From protein structures to a shortlist</h2>
      <div class="pt-flow pt-steps">
        <div class="pt-card"><h3>Protein Data Bank</h3><p>every metal site, reduced to its first coordination shell</p></div>
        <span class="pt-arrow">→</span>
        <div class="pt-card"><h3>Descriptor funnel</h3><p>each shell rebuilt in a support; GFN2-xTB and cheap
          descriptors for activity and stability</p></div>
        <span class="pt-arrow">→</span>
        <div class="pt-card"><h3>DFT</h3><p>the best candidates refined with Quantum ESPRESSO (PBE-D3)</p></div>
        <span class="pt-arrow">→</span>
        <div class="pt-card"><h3>Learning from DFT</h3><p>the DFT data calibrate the xTB screen and train a MACE
          force field</p></div>
        <span class="pt-arrow">→</span>
        <div class="pt-card"><h3>Literature agents</h3><p>Valency agents check what has been reported and benchmark
          the models</p></div>
      </div>
    </section>""")
    return


@app.cell
def _(dft, mo):
    window = mo.ui.slider(0.1, 1.0, step=0.05, value=0.2, label="Target window |ΔG_H| ≤", show_value=True)
    metals = mo.ui.multiselect(sorted(dft.metal.unique()), value=sorted(dft.metal.unique()), label="Metals")
    mo.vstack([
        mo.Html('<section class="pt-section"><div class="pt-kicker">3 · Hydrogen evolution on graphene (DFT)</div>'
                '<h2>Which motifs bind hydrogen near the optimum, and keep their metal?</h2><p class="pt-lede">'
                'ΔG_H close to 0 eV makes a good H₂ evolution site. The colour shows the stability test: a negative '
                'formation energy E_f means the metal prefers the site over a metal particle; U_diss &gt; 0 means it '
                'also resists dissolution. Hover a point for the motif.</p></section>'),
        mo.hstack([metals, window], justify="start", gap=2)], gap=0.6)
    return metals, window


@app.cell
def _(alt, dft, metals, mo, window):
    _d = dft[dft.metal.isin(metals.value)]
    _w = window.value
    _band = alt.Chart().mark_rect(color="#dcefe4", opacity=0.6).encode(y=alt.datum(-_w), y2=alt.datum(_w))
    _stab = alt.Scale(domain=["holds its metal, also in water (U_diss > 0)", "holds its metal (E_f < 0)",
                              "prefers metal particles", "not computed"],
                      range=["#2f6b4f", "#7fb89a", "#c4682f", "#b5bfc3"])
    _pts = alt.Chart(_d).mark_circle(size=90, opacity=0.9, stroke="white", strokeWidth=1.5).encode(
        x=alt.X("metal:N", title=None, sort=sorted(dft.metal.unique()), axis=alt.Axis(labelAngle=0),
                scale=alt.Scale(paddingInner=0.4)),
        xOffset=alt.XOffset("jitter:Q"),
        y=alt.Y("dG_H:Q", title="ΔG_H (eV)", scale=alt.Scale(domain=[-1.5, 1.5], clamp=True, nice=False)),
        color=alt.Color("stability:N", scale=_stab, title="Stability",
                        legend=alt.Legend(orient="bottom", columns=2)),
        tooltip=[alt.Tooltip("site:N", title="motif"), alt.Tooltip("dG_H:Q", format="+.2f", title="ΔG_H (eV)"),
                 alt.Tooltip("H_on:N", title="H on"), alt.Tooltip("E_f:Q", format="+.2f", title="E_f (eV)"),
                 alt.Tooltip("U_diss:Q", format="+.2f", title="U_diss (V)"),
                 alt.Tooltip("enzyme:N", title="origin")],
    ).transform_calculate(jitter="random()")
    _chart = alt.layer(_band, _pts, data=_d).properties(width="container", height=360)
    _near = (_d[_d.dG_H.round(2).abs() <= _w + 1e-9].sort_values("dG_H", key=abs)
             [["site", "dG_H", "H_on", "E_f", "U_diss", "stability", "enzyme"]]
             .rename(columns={"site": "motif", "dG_H": "ΔG_H (eV)", "H_on": "H on", "E_f": "E_f (eV)",
                              "U_diss": "U_diss (V)", "enzyme": "origin"}).round(2))
    mo.vstack([mo.ui.altair_chart(_chart, chart_selection=False, legend_selection=False),
               mo.md(f"**{len(_near)} motifs** within ±{_w:.2f} eV (values rounded to 0.01 eV):"),
               mo.ui.table(_near, selection=None, page_size=10, show_column_summaries=False)], gap=0.6)
    return


@app.cell
def _(mo, offsets, pretty, xtb):
    # calibrated xTB: raw GFN2-xTB ΔG_H plus a per-metal offset fitted to the DFT pairs (H on S: one offset)
    _used = offsets[offsets.used]
    _off_m = dict(zip(_used[_used.state == "H_M"].metal, _used[_used.state == "H_M"].offset))
    _off_s = float(_used[_used.state == "H_S"].offset.iloc[0])
    _x = xtb[xtb.converged.astype(bool)].copy()
    _x["metal"] = _x.motif.str.split("-").str[0]
    _x["offset"] = [(_off_m.get(m) if s == "H_M" else _off_s) for m, s in zip(_x.metal, _x.state)]
    _x["dG_H_cal"] = _x.dG_H + _x.offset
    funnel = (_x.dropna(subset=["dG_H_cal"]).assign(a=lambda d: d.dG_H_cal.abs()).sort_values("a")
              .drop_duplicates("motif"))
    funnel = funnel.assign(site=funnel.motif.map(pretty),
                           arrangement=funnel.motif.str.split(":").str[1:].str.join(" "),
                           origin=funnel.source.map({"scan": "enzyme first shell", "library": "systematic library"}),
                           H_on=funnel.state.map({"H_M": "metal", "H_S": "S donor"}))
    n_xtb = xtb.motif.nunique()
    _loo = _used[_used.state == "H_M"].set_index("metal").loo_mae.round(2).to_dict()
    mo.Html(f"""<section class="pt-section"><div class="pt-kicker">4 · The xTB funnel</div>
      <h2>{n_xtb:,} graphene sites, screened in minutes each</h2>
      <p class="pt-lede">GFN2-xTB relaxes every site with and without H. A per-metal offset fitted to the DFT values
      calibrates it; the leave-one-out error is about 0.1–0.5 eV (Co {_loo.get('Co')}, Fe {_loo.get('Fe')},
      Cu {_loo.get('Cu')}, Ni {_loo.get('Ni')}, Zn {_loo.get('Zn')} eV), so the funnel ranks, the DFT decides. Below:
      the {len(funnel):,} converged sites whose metal has a calibration. Sites come from enzyme first shells or from a
      systematic library of all N/S/O/C arrangements.</p></section>""")
    return (funnel,)


@app.cell
def _(alt, funnel, metal_color, mo):
    _counts = (funnel.assign(near=funnel.dG_H_cal.abs() <= 0.3)
               .groupby(["metal", "origin"]).near.mean().mul(100).round(0).reset_index(name="share"))
    _bars = alt.Chart(_counts).mark_bar(cornerRadiusEnd=4, height=12).encode(
        y=alt.Y("metal:N", title=None, sort="-x"),
        yOffset=alt.YOffset("origin:N"),
        x=alt.X("share:Q", title="share of sites within ±0.3 eV (%)"),
        color=alt.Color("origin:N", title="Sites from", scale=alt.Scale(range=["#2c6bd3", "#b5bfc3"]),
                        legend=alt.Legend(orient="bottom")),
        tooltip=["metal", "origin", alt.Tooltip("share:Q", title="% near optimum")],
    ).properties(width="container", height=alt.Step(30))
    _tbl = (funnel[["site", "arrangement", "metal", "origin", "H_on", "dG_H_cal"]]
            .rename(columns={"site": "motif", "H_on": "H on", "dG_H_cal": "ΔG_H calibrated (eV)"})
            .round(2).sort_values("ΔG_H calibrated (eV)", key=abs))
    mo.vstack([
        mo.ui.tabs({"Share near the optimum": mo.ui.altair_chart(_bars, chart_selection=False,
                                                                  legend_selection=False),
                    "All sites (search and filter)": mo.ui.table(_tbl, selection=None, page_size=12)}),
        mo.Html('<p class="pt-note">Do enzyme shells do better than arbitrary arrangements? For Cu, Ni and Zn they '
                'do; on partial data this is a hint, not a result.</p>')], gap=0.4)
    _ = metal_color
    return


@app.cell
def _(alt, mo, supports):
    # calibrated xTB on the other supports, only where a per-metal calibration exists
    _s = supports[(supports.state == "H_M") & supports.dG_H_offset.notna()].copy()
    _s["support"] = _s.support.map({"MoS2": "MoS₂", "SiO2": "silica"})
    _base = alt.Chart(_s).encode(x=alt.X("metal:N", title=None, axis=alt.Axis(labelAngle=0)),
                                 y=alt.Y("support:N", title=None, axis=alt.Axis(labelFontSize=13, labelPadding=8)))
    _heat = _base.mark_rect(cornerRadius=4, stroke="white", strokeWidth=2).encode(
        color=alt.Color("dG_H_offset:Q", title="ΔG_H, eV", legend=alt.Legend(gradientLength=100),
                        scale=alt.Scale(scheme="blueorange", domainMid=0, domain=[-1.6, 1.6], clamp=True)),
        tooltip=["support", "metal", alt.Tooltip("dG_H_offset:Q", format="+.2f", title="ΔG_H calibrated (eV)")])
    _txt = _base.mark_text(fontSize=12, font="IBM Plex Mono").encode(
        text=alt.Text("dG_H_offset:Q", format="+.2f"),
        color=alt.condition("abs(datum.dG_H_offset) > 0.9", alt.value("white"), alt.value("#1c2427")))
    mo.vstack([
        mo.Html('<section class="pt-section"><div class="pt-kicker">5 · Beyond graphene</div>'
                '<h2>The same metals on MoS₂ and silica</h2><p class="pt-lede">Enzyme sites mapped onto MoS₂ (a metal '
                'on three S of the sheet, like the Fe of an Fe–S cluster) and silica (siloxide O, like the water and '
                'carboxylate O of carbonic anhydrase or ferritin). Calibrated xTB with H on the metal; white is the '
                'optimum. On MoS₂, Fe is confirmed by DFT.</p></section>'),
        mo.ui.altair_chart((_heat + _txt).properties(width="container", height=130), chart_selection=False,
                           legend_selection=False)], gap=0.4)
    return


@app.cell
def _(mo):
    reaction = mo.ui.dropdown(
        {"O₂ reduction (limiting potential, higher is better)": "U_ORR_cal",
         "CO₂ → CO": "U_CO", "CO₂ → formate": "U_HCOOH", "N₂ → NH₃": "U_NRR"},
        value="O₂ reduction (limiting potential, higher is better)", label="Reaction")
    mo.vstack([
        mo.Html('<section class="pt-section"><div class="pt-kicker">6 · Other reactions</div>'
                '<h2>O₂, CO₂ and N₂ on the same sites</h2><p class="pt-lede">The intermediates of each reaction on '
                'every graphene candidate, with GFN2-xTB. These are trends: for O₂ reduction, our DFT check on Co–N₄ '
                'reproduces the published 0.93 V.</p></section>'),
        reaction], gap=0.6)
    return (reaction,)


@app.cell
def _(alt, metal_color, mo, pretty, reaction, reactions):
    _c = reaction.value
    _r = (reactions[(reactions.kind == "graphene") & ~reactions.calibration_only.astype(bool)]
          .dropna(subset=[_c]).nlargest(15, _c))
    _r = _r.assign(motif=_r.site.map(pretty) + "  " + _r.site.str.split(":").str[2].fillna(""))
    _chart = alt.Chart(_r).mark_bar(cornerRadiusEnd=4, height=14).encode(
        y=alt.Y("motif:N", sort="-x", title=None),
        x=alt.X(f"{_c}:Q", title="limiting potential (V vs RHE)"),
        color=alt.Color("metal:N", scale=metal_color, legend=None),
        tooltip=[alt.Tooltip("motif:N"), alt.Tooltip(f"{_c}:Q", format="+.2f", title="U_L (V)"),
                 alt.Tooltip("dG_H:Q", format="+.2f", title="ΔG_H (eV)")],
    ).properties(width="container", height=alt.Step(22))
    mo.ui.altair_chart(_chart, chart_selection=False, legend_selection=False)
    return


@app.cell
def _(gallery, mo):
    # six motifs from enzyme to catalyst as small 3D models (site/gallery.py); the marks are preliminary
    import json as _json
    _cards = [
        ("fe_s2c", "Nitrogenase FeMo-cofactor", "Fe–S₂C", "graphene · H₂ evolution", "yes", "no", "yes", "Fe"),
        ("fe_ns2", "Fe–S cluster with a His ligand", "Fe–NS₂", "graphene · H₂ evolution", "yes", "no", "yes", "Fe"),
        ("co_silica", "Carbonic anhydrase, Co form", "Co–O₃", "silica · H₂ evolution", "close", "open", "open", "Co"),
        ("fe_mos2", "Fe of an Fe–S cluster", "Fe–S₃", "MoS₂ · H₂ evolution", "yes", "open", "open", "Fe"),
        ("cu_n3", "Cytochrome c oxidase CuB, laccase", "Cu–N₃", "graphene · H₂ evolution", "yes", "no", "no", "Cu"),
        ("co_n2oc", "Metal environment in the PDB", "Co–N₂OC", "graphene · H₂ evolution", "yes", "yes", "yes", "Co"),
    ]
    _mark = {"yes": ("✓", "ok"), "no": ("✗", "no"), "open": ("–", "open"), "close": ("≈", "open")}

    def _chip(label, state):
        sym, cls = _mark[state]
        return f'<span class="chip {cls}">{sym} {label}</span>'

    _html = "".join(
        f'<div class="card"><div class="mol" id="m-{k}"></div><div class="from">{s}</div><div class="arrow">↓</div>'
        f'<div class="to">{m}</div><div class="where">{w}</div>'
        f'<div class="chips">{_chip("active", a)}{_chip("holds its metal", h)}{_chip("new", n)}</div></div>'
        for k, s, m, w, a, h, n, _ in _cards)
    _metals = {c[0]: c[-1] for c in _cards}
    _page = f"""<!doctype html><html><head><meta charset="utf-8">
    <link rel="stylesheet" href="https://fonts.googleapis.com/css2?family=Newsreader:opsz,wght@6..72,600&family=Source+Sans+3:wght@400;600&display=swap">
    <style>
      body {{ margin: 0; font-family: "Source Sans 3", system-ui, sans-serif; color: #1c2427; background: transparent; }}
      .grid {{ display: grid; grid-template-columns: repeat(3, minmax(0, 1fr)); gap: 12px; }}
      @media (max-width: 640px) {{ .grid {{ grid-template-columns: 1fr; }} }}
      .card {{ background: #fff; border: 1px solid #d3dbde; border-radius: 10px; padding: 10px 12px 12px; text-align: center; }}
      .mol {{ position: relative; width: 100%; height: 170px; }}
      .from {{ font-size: 0.86rem; color: #59666b; margin-top: 4px; }}
      .arrow {{ color: #9aa5ab; line-height: 1; }}
      .to {{ font-family: "Newsreader", Georgia, serif; font-size: 1.3rem; font-weight: 600; }}
      .where {{ font-size: 0.82rem; color: #59666b; margin-bottom: 8px; }}
      .chips {{ display: flex; flex-wrap: wrap; gap: 4px; justify-content: center; }}
      .chip {{ font-size: 0.74rem; padding: 1px 8px; border-radius: 999px; font-weight: 600; }}
      .chip.ok {{ background: #dcefe4; color: #2f6b4f; }}
      .chip.no {{ background: #f6e1de; color: #9c3a2c; }}
      .chip.open {{ background: #eef1f2; color: #59666b; }}
    </style>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.4.2/3Dmol-min.js"></script></head>
    <body><div class="grid">{_html}</div>
    <script>
      const XYZ = {_json.dumps(gallery)}, METAL = {_json.dumps(_metals)};
      const COL = {{Fe: "#c4682f", Co: "#3f6fc0", Cu: "#1f8f8a", Ni: "#3f9a6a"}};
      for (const k in XYZ) {{
        const v = $3Dmol.createViewer(document.getElementById("m-" + k), {{backgroundColor: "white", antialias: true}});
        v.addModel(XYZ[k], "xyz");
        v.setStyle({{}}, {{stick: {{radius: 0.13, colorscheme: "Jmol"}}, sphere: {{scale: 0.22, colorscheme: "Jmol"}}}});
        v.setStyle({{elem: "C"}}, {{stick: {{radius: 0.13, color: "#5d6b70"}}, sphere: {{scale: 0.2, color: "#5d6b70"}}}});
        v.setStyle({{elem: METAL[k]}}, {{stick: {{radius: 0.13, color: COL[METAL[k]]}}, sphere: {{scale: 0.42, color: COL[METAL[k]]}}}});
        v.setStyle({{elem: "Mo"}}, {{stick: {{radius: 0.13, color: "#8a68b0"}}, sphere: {{scale: 0.26, color: "#8a68b0"}}}});
        v.zoomTo(); v.rotate(-58, "x"); v.rotate(18, "y"); v.zoom(1.2); v.render();
      }}
    </script></body></html>"""
    mo.vstack([
        mo.Html('<section class="pt-section"><div class="pt-kicker">7 · Gallery</div><h2>Six motifs from enzyme to '
                'catalyst</h2><p class="pt-lede">Each card goes from the enzyme site to the catalyst motif we computed '
                '(relaxed structure, H on the metal; drag to rotate). The marks are preliminary.</p></section>'),
        mo.iframe(_page, height="640px")], gap=0.4)
    return


@app.cell
def _(assets, mo):
    mo.Html(f"""
    <section class="pt-section pt-md">
      <div class="pt-kicker">8 · Will it survive synthesis?</div>
      <h2>Pyrolysis with a machine-learned force field</h2>
      <p class="pt-lede">Many motifs bind hydrogen well, but most would lose their metal during pyrolysis. A MACE force
      field fine-tuned on our DFT frames simulates it at 1500 K: the weak site (Zn–NS₂) loses its metal within a
      fraction of a picosecond, while Co–N₂OC keeps all four donors over 5 ps, as the escape barriers predicted.</p>
      <video src="{assets['md']}" autoplay loop muted playsinline controls></video>
    </section>""")
    return


@app.cell
def _(literature, mo, pd, pretty):
    _label = {"reported_experimental": "made before", "related_only": "only related sites reported",
              "not_found": "no report found"}
    _lit = pd.DataFrame({
        "motif": literature.motif.map(pretty),
        "ΔG_H (eV)": literature.dg,
        "verdict": [_label.get(v["verdict"], v["verdict"]) for v in literature.verdict],
        "closest report": [next((f'{c["title"]} ({c["year"]})' for c in v["closest"]), "") for v in literature.verdict],
        "H₂ evolution reported?": [v["her_reported"] for v in literature.verdict],
    })
    mo.vstack([
        mo.Html('<section class="pt-section"><div class="pt-kicker">9 · Literature check</div><h2>What is new?</h2>'
                '<p class="pt-lede">Multi-agent literature screen on Valency: Claude Haiku triages titles, Claude Sonnet '
                'judges abstracts and tries to refute novelty claims, Claude Opus reads the flagged papers in full. '
                'Final check of the near-optimal DFT motifs:</p></section>'),
        mo.ui.table(_lit, selection=None, show_column_summaries=False, wrapped_columns=["closest report",
                                                                                        "H₂ evolution reported?"])],
        gap=0.4)
    return


@app.cell
def _(mo):
    mo.Html("""
    <section class="pt-section">
      <div class="pt-kicker">10 · Data and code</div>
      <h2>Open dataset, open pipeline</h2>
      <div class="pt-grid">
        <div class="pt-card"><h3>Dataset (Zenodo)</h3><p>Every PDB metal site with its first shell, the support
          mappings, all xTB and DFT results, and about 11,000 DFT and 5,000 xTB structures with energies and forces for
          force-field training. Database and tables CC BY 4.0, structures CC0.</p></div>
        <div class="pt-card"><h3>Code (this repository)</h3><p>PDB sweep, structure builders, xTB screens, Quantum
          ESPRESSO inputs, analysis, database and these notebooks. MIT licence. Run this notebook locally with
          <code>marimo run notebooks/sacred.py</code>.</p></div>
        <div class="pt-card"><h3>Built with AI</h3><p>Claude Opus in Claude Code wrote and ran the pipeline; Valency
          literature agents; a MACE-MP-0 force field fine-tuned on our DFT. The scientific decisions and checks are
          ours.</p></div>
      </div>
    </section>
    <p class="pt-foot">SACred · hackathon release, October 2026. Every number is a draft and needs expert checking.
    Structures from the RCSB PDB (CC0); quantum chemistry with GFN2-xTB (tblite) and Quantum ESPRESSO.</p>""")
    return


if __name__ == "__main__":
    app.run()
