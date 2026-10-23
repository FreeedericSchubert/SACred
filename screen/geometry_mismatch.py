"""How much does the support change the metal's first shell? Enzyme geometry (PDB) against the relaxed SAC model.

    python screen/geometry_mismatch.py      # -> data/geometry_mismatch.csv (per model), data/geometry_mismatch.md

Enzyme side: representative X-ray sites at <= 2.0 A of the nine enzyme metals without a metal neighbour (no clusters),
with the donors of their graphene mapping: the four of the double vacancy or the three of the single vacancy, picked
by bio/first_shell.py (dv_quad, sv_tri) from the donor vectors of the contacts table. Model side: the relaxed bare
site of the same motif, DFT where it finished (dft/runs/graphene/<motif>/bare), else GFN2-xTB (data/screen/xtb), with
its metal and donors from the builder (materials/graphene/sites.py).
Per site: metal-donor distance per donor element, the donor-metal-donor angles, tau4 for the four DV donors, and the
metal's distance from the plane of its donors. The enzyme sites are split by what holds the donors: amino-acid side
chains and backbone (flexible), or a rigid cofactor (porphyrin, Fe-S cluster S, dithiolene). Mismatch = model minus
enzyme; the angles are compared as sorted lists, so no donor-by-donor assignment is needed.
"""
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path[:0] = [os.path.join(ROOT, "bio"), os.path.join(ROOT, "screen")]
import first_shell as fs                 # noqa: E402
import xtb_scan as xs                    # noqa: E402

DATA = os.path.join(ROOT, "data")
REL = os.path.join(DATA, "release", "parquet")
RUNS = os.path.join(ROOT, "dft", "runs", "graphene")
MAX_RES = 2.0
RIGID = ("N_porphyrin", "S_inorganic", "S_dithiolene")
LOST = 0.3            # A beyond the donor cutoff of bio/first_shell.py: the model has let go of that donor


def geometry(vec, els):
    """Distances per element, sorted angles, tau4 (4 donors) and the metal's distance from the donors' plane."""
    vec = np.asarray(vec, float)
    ang = fs.angle_matrix(vec)
    n = len(vec)
    iu = np.triu_indices(n, 1)
    centred = vec - vec.mean(0)
    normal = np.linalg.svd(centred)[2][-1]
    out = {"angles": np.sort(ang[iu]), "h_plane": float(abs(-vec.mean(0) @ normal)),
           "tau4": float(fs.tau4(ang)) if n == 4 else np.nan}
    d = np.linalg.norm(vec, axis=1)
    for el in set(els):
        out["d_" + el] = float(np.mean([x for x, e in zip(d, els) if e == el]))
    return out


def model_geometry(motif):
    """Geometry of the relaxed bare model: DFT if finished, else xTB; None without a structure."""
    from ase.io import read
    d = motif.replace(":", "_").replace("=", "")
    _, info = xs.build(motif)
    atoms, source = None, None
    out = os.path.join(RUNS, d, "bare", "relax.out")
    if os.path.exists(out) and "End final coordinates" in open(out).read():
        atoms, source = read(out, format="espresso-out", index=-1), "DFT"
    elif os.path.exists(os.path.join(xs.OUT, d, "bare.xyz")):
        atoms, source = read(os.path.join(xs.OUT, d, "bare.xyz")), "xTB"
    if atoms is None:
        return None
    vec = atoms.get_distances(info["imetal"], info["idonors"], mic=True, vector=True)
    els = [atoms[i].symbol for i in info["idonors"]]
    g = geometry(vec, els)
    g["source"] = source
    g["lost"] = int(sum(np.linalg.norm(v) > fs.CUTOFF.get(el, 2.6) + LOST for v, el in zip(vec, els)))
    return g


def enzyme_sites():
    """The enzyme-side subsets: one row per site with its model motif name and geometry."""
    s = pd.read_parquet(os.path.join(REL, "sites.parquet"), columns=[
        "site_id", "pdb_id", "metal", "is_representative", "preferred_site", "dv_formula", "dv_isomer", "sv_formula",
        "resolution", "method", "n_metal_neighbors", "cn"])
    s = s[s.is_representative & s.metal.isin(xs.METALS) & s.preferred_site.isin(["DV", "SV"])
          & (s.n_metal_neighbors == 0) & s.method.str.contains("X-RAY", case=False, na=False) & (s.resolution <= MAX_RES)]
    s = s[~(s.dv_formula.fillna("") + s.sv_formula.fillna("")).str.contains("Se|Cl")]
    c = pd.read_parquet(os.path.join(REL, "contacts.parquet"), filters=[("kind", "==", "donor")],
                        columns=["site_id", "rank", "element", "contact_class", "distance", "dx", "dy", "dz"])
    c = c[c.site_id.isin(s.site_id)].sort_values(["site_id", "rank"])
    rows = []
    meta = s.set_index("site_id")
    for sid, g in c.groupby("site_id", sort=False):
        r = meta.loc[sid]
        donors = [{"el": e, "cls": k, "dist": d} for e, k, d in zip(g.element, g.contact_class, g.distance)]
        vec = g[["dx", "dy", "dz"]].to_numpy(float)
        if len(donors) < (4 if r.preferred_site == "DV" else 3):
            continue
        ang = fs.angle_matrix(vec)
        if r.preferred_site == "DV":
            _, sub = fs.dv_quad(donors, ang)
            iso = xs.builder_isomer(r.dv_formula, r.dv_isomer if isinstance(r.dv_isomer, str) else "")
            motif = "%s-%s:DV%s" % (r.metal, r.dv_formula, ":" + iso if iso else "")
        else:
            sub, _ = fs.sv_tri(donors, ang)
            motif = "%s-%s:SV" % (r.metal, r.sv_formula)
        sub = list(sub)
        geo = geometry(vec[sub], [donors[i]["el"] for i in sub])
        holder = "cofactor" if any(donors[i]["cls"] in RIGID for i in sub) else "side chains"
        rows.append(dict(geo, site_id=sid, pdb_id=r.pdb_id, metal=r.metal, site=r.preferred_site, motif=motif,
                         holder=holder, cn=int(r.cn), n_in_support=len(sub)))
    return pd.DataFrame(rows)


def main():
    e = enzyme_sites()
    models = {m: model_geometry(m) for m in sorted(e.motif.unique())}
    models = {m: g for m, g in models.items() if g is not None}
    e = e[e.motif.isin(models)].copy()
    els = ["N", "S", "O", "C"]
    for col in ["lost", "angle_rms", "dh", "dtau4"] + ["dd_" + el for el in els]:
        e[col] = np.nan
    for i, r in e.iterrows():
        g = models[r.motif]
        e.at[i, "lost"] = g["lost"]
        e.at[i, "angle_rms"] = float(np.sqrt(np.mean((g["angles"] - r.angles) ** 2)))
        e.at[i, "dh"] = g["h_plane"] - r.h_plane
        e.at[i, "dtau4"] = g["tau4"] - r.tau4 if r.site == "DV" else np.nan
        for el in els:
            if "d_" + el in g and pd.notna(r.get("d_" + el)):
                e.at[i, "dd_" + el] = g["d_" + el] - r["d_" + el]
    e.drop(columns="angles").to_parquet(os.path.join(DATA, "complete", "geometry_mismatch_sites.parquet"))

    per = []
    for m, g in e.groupby("motif"):
        mod = models[m]
        row = {"motif": m, "site": g.site.iloc[0], "model": mod["source"], "donors_lost": mod["lost"],
               "enzyme_sites": len(g),
               "entries": g.pdb_id.nunique(), "share_cofactor": (g.holder == "cofactor").mean(),
               "enzyme_cn": g.cn.median(), "angle_rms": g.angle_rms.median(),
               "h_enzyme": g.h_plane.median(), "h_model": mod["h_plane"],
               "tau4_enzyme": g.tau4.median(), "tau4_model": mod["tau4"]}
        for el in els:
            if "d_" + el in mod and g["d_" + el].notna().any():
                row["d_%s_enzyme" % el] = g["d_" + el].median()
                row["d_%s_model" % el] = mod["d_" + el]
        per.append(row)
    per = pd.DataFrame(per).sort_values("enzyme_sites", ascending=False)
    per.to_csv(os.path.join(DATA, "geometry_mismatch.csv"), index=False, float_format="%.3f")

    lines = ["# How much the support changes the first shell", "",
             "%d enzyme sites (X-ray <= %.1f A, mononuclear, nine enzyme metals) on %d graphene models (%d DFT, %d xTB). "
             "Model minus enzyme; medians with the 25-75 %% range." % (
                 len(e), MAX_RES, len(models), sum(g["source"] == "DFT" for g in models.values()),
                 sum(g["source"] == "xTB" for g in models.values())), "",
             "In %.0f %% of the enzyme sites the model lets go of at least one donor (longer than the donor cutoff + "
             "%.1f A): %s. The table below covers the sites whose model keeps all its donors." % (
                 100 * (e.lost > 0).mean(), LOST, ", ".join("%s %.0f %%" % (k, 100 * v) for k, v in
                                                            e.groupby("site").lost.apply(lambda x: (x > 0).mean()).items())),
             "",
             "| | sites | M-N (A) | M-S (A) | M-O (A) | angles, RMS (deg) | metal out of donor plane (A) | tau4 (DV) |",
             "|---|---|---|---|---|---|---|---|"]

    def q(x):
        x = x.dropna()
        return "%+.2f (%+.2f to %+.2f)" % (x.median(), x.quantile(0.25), x.quantile(0.75)) if len(x) else "-"

    for (site, holder), g in e[e.lost == 0].groupby(["site", "holder"]):
        lines.append("| %s, %s | %d | %s | %s | %s | %.0f (%.0f-%.0f) | %s | %s |" % (
            site, holder, len(g), q(g.dd_N), q(g.dd_S), q(g.dd_O), g.angle_rms.median(), g.angle_rms.quantile(0.25),
            g.angle_rms.quantile(0.75), q(g.dh), q(g.dtau4) if site == "DV" else "-"))
    lost = per[per.donors_lost > 0].sort_values("enzyme_sites", ascending=False)
    lines += ["", "Models that let go of a donor (by enzyme sites behind them): " + "; ".join(
        "%s (%s, %d sites, %d lost)" % (r.motif, r.model, r.enzyme_sites, r.donors_lost) for r in lost.head(10).itertuples())
              + ".", "", "Largest angular mismatch among the models that keep their donors (>= 10 enzyme sites):", ""]
    top = per[(per.enzyme_sites >= 10) & (per.donors_lost == 0)].sort_values("angle_rms", ascending=False).head(8)
    lines += ["- %s (%s, %d sites): %.0f deg; metal %.2f A out of the donor plane in enzymes, %.2f A in the model" % (
        r.motif, r.model, r.enzyme_sites, r.angle_rms, r.h_enzyme, r.h_model) for r in top.itertuples()]
    open(os.path.join(DATA, "geometry_mismatch.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
