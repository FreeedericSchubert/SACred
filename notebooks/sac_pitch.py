import marimo

__generated_with = "0.25.2"
app = marimo.App(width="medium", app_title="SACred", css_file="sac_pitch.css",
                 html_head_file="sac_explorer_head.html")


@app.cell
def _():
    import base64
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
            "title": {"font": "Newsreader", "fontSize": 16, "color": _ink}}}

    def pretty(motif):
        """Fe-NS2:SV -> Fe–NS₂; arrangement labels dropped."""
        return motif.split(":")[0].replace("-", "–", 1).translate(str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉"))

    def fmt(n):
        return f"{n:,}"
    return alt, base64, fmt, mo, os, pd, pretty


@app.cell
def _(base64, mo, os, pd):
    # the few numbers this page shows are read from the pipeline's files, so a refresh updates the pitch
    ROOT = os.path.dirname(str(mo.notebook_dir()))
    DATA = os.path.join(ROOT, "data")
    sites = pd.read_parquet(os.path.join(DATA, "release", "parquet", "sites.parquet"), columns=["pdb_id"])
    xtb = pd.read_csv(os.path.join(DATA, "xtb_scan.csv"))
    dgh = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv"))
    stab = pd.read_csv(os.path.join(DATA, "stability_graphene.csv"))
    _zoom = os.path.join(ROOT, "site", "zoom_viewer.html")
    zoom_html = open(_zoom).read() if os.path.exists(_zoom) else ""

    def _uri(name, mime):
        # logo and MD movie inlined, so the static export carries them too
        path = os.path.join(str(mo.notebook_dir()), "assets", name)
        return "data:%s;base64,%s" % (mime, base64.b64encode(open(path, "rb").read()).decode()) if os.path.exists(path) else ""

    assets = {"logo": _uri("sacred_logo.webp", "image/webp"), "md": _uri("md_pyrolysis_1500K.mp4", "video/mp4")}
    return assets, dgh, sites, stab, xtb, zoom_html


@app.cell
def _(dgh, stab):
    # DFT sites with both a hydrogen binding energy and a formation energy; the candidates pass both tests
    _d = dgh.dropna(subset=["dG_H"]).assign(a=lambda d: d.dG_H.abs()).sort_values("a").drop_duplicates("motif")
    both = _d[["motif", "dG_H"]].merge(stab.dropna(subset=["E_f"])[["motif", "E_f", "U_diss"]], on="motif")
    passing = both[(both.dG_H.abs() <= 0.2) & (both.E_f < 0) & (both.U_diss > 0)]
    n_dft = int(_d.motif.nunique())
    n_usable = int((_d.dG_H.round(2).abs() <= 0.2).sum())       # DFT motifs near the hydrogen optimum
    return both, n_dft, n_usable, passing


@app.cell
def _(assets, fmt, mo, n_dft, n_usable, sites, xtb):
    _facts = [(fmt(len(sites)), "metal sites from the Protein Data Bank screened"),
              (fmt(int(xtb.motif.nunique())), "single-atom sites in the xTB funnel"),
              (str(n_dft), "candidates refined with DFT"),
              (str(n_usable), "usable motifs found (preliminary)")]
    _tiles = "".join(f'<div class="sac-fact{" good" if k == 3 else ""}"><b>{v}</b><span>{t}</span></div>'
                     for k, (v, t) in enumerate(_facts))
    mo.Html(f"""
    <header class="sac-hero">
      <div class="pt-brand"><img src="{assets['logo']}" alt="SACred logo"><div>
      <div class="sac-eyebrow">SACred · bioinspired single-atom catalysts from the PDB</div>
      <h1>Turning enzyme active sites into usable single-atom catalysts</h1></div></div>
      <p class="sac-lede">An agentic AI pipeline and open dataset: every metal site in the Protein Data Bank becomes a
      single-atom catalyst candidate, and the pipeline finds the few that are both active and can actually be made,
      opening up new synthesis strategies.</p>
      <div class="sac-facts">{_tiles}</div>
    </header>""")
    return


@app.cell
def _(mo):
    mo.Html("""
    <section class="pt-section">
      <div class="pt-kicker">1 · The problem</div>
      <h2>Single-atom catalysts are designed by trial and error</h2>
      <div class="pt-grid">
        <div class="pt-card"><h3>Promising, but hard to make</h3><p>Single metal atoms in carbon could replace
          platinum for hydrogen, fuel cells and CO₂ conversion. But during synthesis most metals clump into particles,
          and whether a site survives only shows up in the lab.</p></div>
        <div class="pt-card"><h3>Nature already solved it</h3><p>Enzymes hold single metals in pockets of N, S and O
          and run the same reactions: hydrogenases make H₂, nitrogenases fix N₂. Bio-inspired catalysts so far copy a
          handful of hand-picked enzymes.</p></div>
        <div class="pt-card"><h3>Our idea</h3><p>Use <b>all of them</b>: every metal site in the Protein Data Bank
          as a template, and test each one in silico for activity <b>and</b> for survival in synthesis before anyone
          goes to the lab.</p></div>
      </div>
    </section>""")
    return


@app.cell
def _(mo):
    _steps = [("Protein Data Bank", "every metal site, reduced to its first coordination shell"),
              ("Descriptor funnel", "each shell rebuilt in graphene, MoS₂ or silica; GFN2-xTB and cheap descriptors "
               "for activity and stability"),
              ("DFT", "the best candidates refined with density functional theory"),
              ("Learning from DFT", "the DFT data calibrate the xTB screen and train a MACE force field for "
               "dynamics"),
              ("Candidates", "sites that are active, stable and survive pyrolysis")]
    _flow = '<span class="pt-arrow">→</span>'.join(
        f'<div class="pt-card"><h3>{t}</h3><p>{d}</p></div>' for t, d in _steps)
    mo.Html(f"""
    <section class="pt-section">
      <div class="pt-kicker">2 · What we built</div>
      <h2>An agentic AI pipeline from protein structures to catalysts</h2>
      <p class="pt-lede">We screened the PDB and used a funnel of descriptors and semi-empirical methods to find good
      candidates, which we treated with DFT. The DFT data in turn calibrate the semi-empirical screen and train a
      machine-learned force field (MACE) that simulates the synthesis.</p>
      <div class="pt-flow pt-steps">{_flow}</div>
    </section>""")
    return


@app.cell
def _(mo, zoom_html):
    _zoom = (mo.iframe('<!doctype html><html data-theme="light"><meta charset="utf-8"><meta name="viewport" '
                       'content="width=device-width, initial-scale=1">' + zoom_html, height="820px") if zoom_html else
             mo.callout(mo.md("Build the 3D view first: `python site/zoom_viewer.py`."), kind="warn"))
    mo.vstack([
        mo.Html('<div class="pt-kicker" style="margin-top:6px">From enzyme to catalyst</div><p class="pt-lede">Press '
                'Play: from the whole enzyme to its active site, to the first shell, and to the same metal as a '
                'single-atom site on graphene, MoS₂ and silica.</p>'),
        _zoom], gap=0.4)
    return


@app.cell
def _(mo):
    mo.Html("""
    <section class="pt-section">
      <div class="pt-kicker">3 · Agents and tools</div>
      <h2>A multi-agent literature screen and interactive results</h2>
      <div class="pt-grid">
        <div class="pt-card"><h3>Valency: a multi-agent literature screen</h3><p>Teams of agents searched the
          literature for simulated and experimental data on single-atom catalysts. We used it to validate and
          benchmark our models, to check which candidates are new, and to find the biggest problems in synthesising
          single-atom catalysts, each of which became a calculation in the pipeline.</p></div>
        <div class="pt-card"><h3>marimo: datasets and results made visible</h3><p>We used marimo to visualise the
          datasets and results: an explorer of all sites with the 3D view and a chat agent over the database, a
          synthesis dashboard, and this pitch.</p></div>
      </div>
    </section>""")
    return


@app.cell
def _(mo):
    mo.Html("""<section class="pt-section">
          <div class="pt-kicker">4 · What it finds</div>
          <h2>A shortlist instead of trial and error</h2>
          <p class="pt-lede">The same pipeline ranks motifs for several reactions and supports. Preliminary results:</p>
          <table class="pt-table">
            <thead><tr><th>Reaction · support</th><th>Most promising motifs</th><th>Basis</th></tr></thead>
            <tbody>
              <tr><td>H₂ evolution · graphene</td><td>Co–N₂OC, Co–N₃C, Fe–NS₂, Fe–S₂C; Ni–N₃ has already been made
                and tested</td><td>DFT</td></tr>
              <tr><td>H₂ evolution · MoS₂</td><td>Fe on three S of the sheet, like the Fe of an Fe–S cluster: brings
                MoS₂ from far off to near the optimum</td><td>DFT</td></tr>
              <tr><td>H₂ evolution · silica</td><td>Co held by siloxide O, as in carbonic anhydrase</td>
                <td>xTB</td></tr>
              <tr><td>O₂ reduction · graphene</td><td>Co sites: Co–N₄ (our DFT reproduces the literature), Co–N₂OC,
                Co–NS₂C</td><td>DFT check, xTB ranking</td></tr>
              <tr><td>CO₂ → CO · graphene</td><td>Mo and W sites (Mo–O₂C₂, Mo–NS₂O, W–S₄)</td><td>xTB trends</td></tr>
              <tr><td>CO₂ → formate · graphene</td><td>Fe and Cu sites (Fe–O₂C₂, Fe–NS₂O, Cu–NOC), which also make
                H₂</td><td>xTB trends</td></tr>
            </tbody>
          </table>
          </section>""")
    return


@app.cell
def _(mo, os):
    # gallery: six motifs from enzyme to catalyst, each a small 3D model (site/gallery.py); marks are preliminary
    import json as _json
    _xyz = _json.load(open(os.path.join(str(mo.notebook_dir()), "assets", "gallery.json")))
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
        f'<div class="card"><div class="mol" id="m-{k}"></div><div class="from">{src}</div><div class="arrow">↓</div>'
        f'<div class="to">{motif}</div><div class="where">{where}</div>'
        f'<div class="chips">{_chip("active", a)}{_chip("holds its metal", h)}{_chip("new", n)}</div></div>'
        for k, src, motif, where, a, h, n, _m in _cards)
    _metals = {k: m for k, *_rest, m in [(c[0],) + c[1:] for c in _cards]}
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
      const XYZ = {_json.dumps(_xyz)}, METAL = {_json.dumps(_metals)};
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
        mo.Html('<div class="pt-kicker" style="margin-top:6px">From enzyme to catalyst: six motifs</div><p class="pt-lede">'
                'Each card goes from the enzyme site to the catalyst motif we computed (relaxed structure, H on the '
                'metal; drag to rotate). The marks are preliminary.</p>'),
        mo.iframe(_page, height="640px")], gap=0.4)
    return


@app.cell
def _(assets, mo):
    mo.Html(f"""
    <div class="pt-md">
      <div class="pt-kicker" style="margin-top:6px">Will it survive synthesis?</div>
      <p class="pt-lede">Many motifs bind hydrogen well, but most would lose their metal during pyrolysis. Our force
      field simulates it, here at 1500 K: the weak site loses its metal within a fraction of a picosecond, the
      promising motif keeps it, as the screen predicted.</p>
      <video src="{assets['md']}" autoplay loop muted playsinline controls></video>
    </div>""") if assets["md"] else None
    return


@app.cell
def _(mo):
    mo.Html("""
    <section class="pt-section">
      <div class="pt-kicker">5 · What comes next</div>
      <h2>From a hackathon pipeline to a design tool for chemists</h2>
      <div class="pt-grid">
        <div class="pt-card"><h3>Pick a site before the lab</h3><p>For any metal and reaction: which first shells
          nature uses, how active they are, and whether they survive synthesis, including the temperature window for
          pyrolysis.</p></div>
        <div class="pt-card"><h3>Scale it up</h3><p>All transition metals, more supports and reactions; the pipeline
          reruns with every PDB release, and the force field takes over more of the expensive DFT.</p></div>
        <div class="pt-card"><h3>Closer to reality</h3><p>The enzyme's second coordination sphere, pH and potential,
          and competing particles and sulfides, so that the predictions match the conditions in the reactor.</p></div>
        <div class="pt-card"><h3>Make it</h3><p>Synthesis of the best candidates with experimental partners, and an
          open dataset to train the next generation of machine-learning models for catalyst stability.</p></div>
      </div>
    </section>""")
    return


@app.cell
def _(mo):
    mo.Html("""
    <p class="pt-foot">SACred · drafts from a hackathon: every number needs expert checking. Structures from the RCSB
    PDB (CC0); quantum chemistry with GFN2-xTB and Quantum ESPRESSO; force field MACE; literature agents on Valency;
    built with Claude.</p>""")
    return


if __name__ == "__main__":
    app.run()
