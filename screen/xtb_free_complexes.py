"""The third geometry: the enzyme's first shell as a free molecule, between the protein and the support.

    python screen/xtb_free_complexes.py run      # build and relax the complexes (GFN2-xTB, vacuum); resumes
    python screen/xtb_free_complexes.py table    # -> data/geometry_three_way.csv and data/geometry_three_way.md

For the most common side-chain cores of Mn, Fe, Co, Ni, Cu and Zn (screen/geometry_mismatch.py) the best-resolved
enzyme sites are rebuilt as complexes: the metal as M(II) with its whole first shell as model ligands placed on the
enzyme's donor vectors (His -> imidazole, Cys -> methanethiolate, Asp/Glu/acetate -> acetate, water, Met ->
methanethiol, Tyr -> methanolate, Ser/Thr -> methanol, backbone or Asn/Gln O -> formaldehyde, Lys or backbone N -> NH3;
sites with other donors are left out), then relaxed without any constraint. Of each complex the donors that go into the
support are compared, as in geometry_mismatch.py: enzyme -> free complex is what the protein imposes, free complex ->
support model what the support imposes. Settings of the scan (spin-restricted, 1000 K, BFGS to 0.05 eV/A).
"""
import json
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from ase import Atoms                    # noqa: E402
from ase.io import write                 # noqa: E402
from ase.optimize import BFGS            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import geometry_mismatch as gm           # noqa: E402

ROOT, DATA = gm.ROOT, gm.DATA
OUT = os.path.join(DATA, "screen", "xtb_complexes")
METALS = ["Mn", "Fe", "Co", "Ni", "Cu", "Zn"]
N_MOTIFS, N_SITES = 30, 3
WORKERS = int(os.environ.get("COMPLEX_WORKERS", 4))
SCF_RETRY = {"max_iterations": 1000, "mixer_damping": 0.1}
LIGAND = {"N_His": "imidazole", "S_Cys": "thiolate", "O_Asp": "acetate", "O_Glu": "acetate", "O_ACT": "acetate",
          "O_water": "water", "S_Met": "thiol", "O_Tyr": "methanolate", "O_Ser/Thr": "methanol",
          "O_backbone": "formaldehyde", "O_ASN": "formaldehyde", "O_GLN": "formaldehyde", "N_Lys": "ammonia",
          "N_backbone": "ammonia"}
CHARGE = {"thiolate": -1, "acetate": -1, "methanolate": -1}


# ---------------------------------------------------------------- ligands (donor at the origin, away from the metal = +z)

def _dir(angle_from_z, phi=0.0):
    a, p = np.radians(angle_from_z), np.radians(phi)
    return np.array([np.sin(a) * np.cos(p), np.sin(a) * np.sin(p), np.cos(a)])


def _methyl(c, towards):
    """Three H on carbon c, tetrahedral about the bond from c in direction `towards`."""
    b = towards / np.linalg.norm(towards)
    p = np.cross(b, [0.3, 1.0, 0.2])
    p /= np.linalg.norm(p)
    q = np.cross(b, p)
    return [c + 1.09 * (np.cos(np.radians(109.5)) * b + np.sin(np.radians(109.5)) * (np.cos(t) * p + np.sin(t) * q))
            for t in np.radians([0, 120, 240])]


def ligand(kind):
    """(symbols, positions): the donor first."""
    if kind == "water":
        return ["O", "H", "H"], [np.zeros(3), [0.757, 0, 0.586], [-0.757, 0, 0.586]]
    if kind == "ammonia":
        return ["N"] + ["H"] * 3, [np.zeros(3)] + [1.01 * _dir(70.5, p) for p in (0, 120, 240)]
    if kind == "imidazole":
        c = np.array([0, 0, 1.13])
        ring = [c + 1.13 * np.array([np.sin(t), 0, -np.cos(t)]) for t in np.radians([0, 72, 144, 216, 288])]
        hs = [r + 1.05 * (r - c) / np.linalg.norm(r - c) for r in ring[1:]]
        return ["N", "C", "C", "N", "C", "H", "H", "H", "H"], ring + hs
    if kind in ("thiolate", "thiol"):
        cpos = 1.82 * _dir(70)
        sym, pos = ["S", "C"], [np.zeros(3), cpos]
        sym += ["H"] * 3
        pos += _methyl(cpos, -cpos)
        if kind == "thiol":
            sym.append("H")
            pos.append(1.34 * _dir(70, 180))
        return sym, pos
    if kind in ("methanol", "methanolate"):
        cpos = 1.42 * _dir(55)
        sym, pos = ["O", "C", "H", "H", "H"], [np.zeros(3), cpos] + _methyl(cpos, -cpos)
        if kind == "methanol":
            sym.append("H")
            pos.append(0.96 * _dir(55, 180))
        return sym, pos
    if kind == "acetate":
        cc = 1.27 * _dir(55)
        o2 = cc + 1.25 * np.array([0, 0, 1.0])
        cm = cc + 1.52 * np.array([0.887, 0, -0.462])
        return ["O", "C", "O", "C", "H", "H", "H"], [np.zeros(3), cc, o2, cm] + _methyl(cm, cc - cm)
    if kind == "formaldehyde":
        cpos = 1.21 * _dir(40)
        return ["O", "C", "H", "H"], [np.zeros(3), cpos, cpos + 1.10 * np.array([-0.326, 0, 0.946]),
                                      cpos + 1.10 * np.array([0.988, 0, -0.156])]
    raise KeyError(kind)


def place(sym, pos, donor, spin):
    """Ligand atoms with the donor at `donor` and +z along the metal-donor direction, turned by `spin` about it."""
    u = donor / np.linalg.norm(donor)
    e1 = np.cross(u, [0.31, 0.95, 0.07])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(u, e1)
    c, s = np.cos(spin), np.sin(spin)
    rot = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])
    frame = np.c_[e1, e2, u]
    return [donor + frame @ rot @ np.asarray(p, float) for p in pos]


def build(metal, donors):
    """The complex: metal at the origin, each ligand turned to keep clear of the atoms placed before it."""
    atoms = Atoms(metal, positions=[np.zeros(3)])
    idonor, charge = [], 2
    for d in donors:
        kind = LIGAND[d["cls"]]
        sym, pos = ligand(kind)
        best = None
        for spin in np.radians(np.arange(0, 360, 30)):
            p = place(sym, pos, d["vec"], spin)
            clear = min(np.linalg.norm(np.asarray(p)[1:, None, :] - atoms.positions[None, :, :], axis=2).min(), 9.0) \
                if len(p) > 1 else 9.0
            if best is None or clear > best[0]:
                best = (clear, p)
        idonor.append(len(atoms))
        atoms += Atoms(sym, positions=best[1])
        charge += CHARGE.get(kind, 0)
    return atoms, idonor, charge


# ---------------------------------------------------------------- sites and runs

def sites():
    """The enzyme sites to rebuild: N_SITES best-resolved per common side-chain core, all donors mappable."""
    g = pd.read_parquet(os.path.join(DATA, "complete", "geometry_mismatch_sites.parquet"))
    g = g[(g.holder == "side chains") & (g.lost == 0) & g.metal.isin(METALS)]
    res = pd.read_parquet(os.path.join(gm.REL, "sites.parquet"), columns=["site_id", "resolution"])
    g = g.merge(res, on="site_id")
    c = pd.read_parquet(os.path.join(gm.REL, "contacts.parquet"), filters=[("kind", "==", "donor")],
                        columns=["site_id", "rank", "element", "contact_class", "distance", "dx", "dy", "dz"])
    c = c[c.site_id.isin(g.site_id)]
    ok = c.groupby("site_id").contact_class.apply(lambda s: s.isin(list(LIGAND)).all())
    g = g[g.site_id.map(ok).fillna(False)]
    top = g.motif.value_counts().head(N_MOTIFS).index
    chosen = g[g.motif.isin(top)].sort_values("resolution").groupby("motif").head(N_SITES)
    out = []
    for r in chosen.itertuples():
        d = c[c.site_id == r.site_id].sort_values("rank")
        donors = [{"el": e, "cls": k, "dist": x, "vec": [dx, dy, dz]}
                  for e, k, x, dx, dy, dz in zip(d.element, d.contact_class, d.distance, d.dx, d.dy, d.dz)]
        out.append({"site_id": r.site_id, "motif": r.motif, "metal": r.metal, "site": r.site, "donors": donors})
    return out


def subset(job):
    vec = np.array([d["vec"] for d in job["donors"]], float)
    ang = gm.fs.angle_matrix(vec)
    if job["site"] == "DV":
        return list(gm.fs.dv_quad(job["donors"], ang)[1])
    return list(gm.fs.sv_tri(job["donors"], ang)[0])


def relax(job):
    from tblite.ase import TBLite
    res_file = os.path.join(OUT, job["site_id"] + ".json")
    if os.path.exists(res_file):
        return json.load(open(res_file))
    os.makedirs(OUT, exist_ok=True)
    for d in job["donors"]:
        d["vec"] = np.asarray(d["vec"], float)
    atoms, idonor, charge = build(job["metal"], job["donors"])
    res = {"site_id": job["site_id"], "motif": job["motif"], "metal": job["metal"], "site": job["site"],
           "charge": charge, "n_donors": len(idonor), "n_atoms": len(atoms), "idonor": idonor,
           "subset": subset(job), "classes": [d["cls"] for d in job["donors"]]}
    t = time.time()
    for scf in ({}, SCF_RETRY):
        try:
            a = atoms.copy()
            a.calc = TBLite(method="GFN2-xTB", charge=charge, electronic_temperature=1000.0, verbosity=0, **scf)
            opt = BFGS(a, logfile=None, maxstep=0.15)
            conv = opt.run(fmax=0.05, steps=500)
            res.update(energy=float(a.get_potential_energy()), steps=opt.nsteps, converged=bool(conv), retried=bool(scf))
            atoms = a
            break
        except Exception as exc:
            res.update(energy=None, error="%s: %s" % (type(exc).__name__, exc))
    res["seconds"] = round(time.time() - t, 1)
    if res.get("energy") is not None:
        vec = atoms.positions[idonor] - atoms.positions[0]
        res["vec"] = vec.round(4).tolist()
        dist = np.linalg.norm(vec, axis=1)
        res["lost"] = int(sum(x > gm.fs.CUTOFF.get(atoms[i].symbol, 2.6) + gm.LOST for x, i in zip(dist, idonor)))
    write(os.path.join(OUT, job["site_id"] + ".xyz"), atoms)
    json.dump(res, open(res_file, "w"), indent=1)
    return res


def run():
    jobs = sites()
    print("%d sites of %d cores" % (len(jobs), len({j["motif"] for j in jobs})), flush=True)
    t = time.time()
    with Pool(WORKERS) as pool:
        for k, r in enumerate(pool.imap_unordered(relax, jobs), 1):
            print("  %d/%d %s %s: %s, %.0f s (%.0f min)" % (
                k, len(jobs), r["site_id"], r["motif"], "E=%.2f, %d donors lost" % (r["energy"], r.get("lost", 0))
                if r.get("energy") is not None else r.get("error"), r["seconds"], (time.time() - t) / 60), flush=True)


def table():
    """Enzyme, free complex and support for the donors that go into the support."""
    rows = []
    enzyme = pd.read_parquet(os.path.join(DATA, "complete", "geometry_mismatch_sites.parquet")).set_index("site_id")
    contacts = pd.read_parquet(os.path.join(gm.REL, "contacts.parquet"), filters=[("kind", "==", "donor")],
                               columns=["site_id", "rank", "element", "dx", "dy", "dz"])
    models = {}
    for f in sorted(os.listdir(OUT)):
        if not f.endswith(".json"):
            continue
        r = json.load(open(os.path.join(OUT, f)))
        if r.get("energy") is None:
            continue
        sub = r["subset"]
        d = contacts[contacts.site_id == r["site_id"]].sort_values("rank")
        els = list(d.element)
        e_vec = d[["dx", "dy", "dz"]].to_numpy(float)[sub]
        f_vec = np.asarray(r["vec"], float)[sub]
        sub_els = [els[i] for i in sub]
        ge, gf = gm.geometry(e_vec, sub_els), gm.geometry(f_vec, sub_els)
        if r["motif"] not in models:
            models[r["motif"]] = gm.model_geometry(r["motif"])
        gs = models[r["motif"]]
        rms = lambda a, b: float(np.sqrt(np.mean((a - b) ** 2)))   # noqa: E731
        row = {"site_id": r["site_id"], "motif": r["motif"], "metal": r["metal"], "site": r["site"],
               "free_lost": r["lost"], "support_lost": gs["lost"] if gs else np.nan,
               "angles_enzyme_free": rms(ge["angles"], gf["angles"]),
               "angles_free_support": rms(gf["angles"], gs["angles"]) if gs else np.nan,
               "angles_enzyme_support": rms(ge["angles"], gs["angles"]) if gs else np.nan,
               "h_enzyme": ge["h_plane"], "h_free": gf["h_plane"], "h_support": gs["h_plane"] if gs else np.nan}
        for el in set(sub_els):
            row["d_%s_enzyme" % el], row["d_%s_free" % el] = ge["d_" + el], gf["d_" + el]
            row["d_%s_support" % el] = gs.get("d_" + el, np.nan) if gs else np.nan
        rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(os.path.join(DATA, "geometry_three_way.csv"), index=False, float_format="%.3f")
    lines = ["# Enzyme, free complex and support: who changes the geometry?", "",
             "%d enzyme sites of %d cores (Mn-Zn, side-chain donors) rebuilt as free M(II) complexes with model ligands "
             "and relaxed in vacuum (GFN2-xTB); compared on the donors that go into the support. Medians." % (
                 len(t), t.motif.nunique()), "",
             "| | sites | angles: enzyme -> free (deg) | free -> support (deg) | enzyme -> support (deg) | metal out of donor "
             "plane: enzyme / free / support (A) | free complex lost donors |", "|---|---|---|---|---|---|---|"]
    for (site, metal), g in list(t.groupby(["site", "metal"])) + [(("all", ""), t)]:
        lines.append("| %s %s | %d | %.0f | %.0f | %.0f | %.2f / %.2f / %.2f | %d |" % (
            site, metal, len(g), g.angles_enzyme_free.median(), g.angles_free_support.median(),
            g.angles_enzyme_support.median(), g.h_enzyme.median(), g.h_free.median(), g.h_support.median(),
            int((g.free_lost > 0).sum())))
    share = (t.angles_enzyme_free / (t.angles_enzyme_free + t.angles_free_support)).median()
    lines += ["", "The protein accounts for a median %.0f %% of the angle changes along enzyme -> free -> support, the "
              "support for the rest." % (100 * share), ""]
    open(os.path.join(DATA, "geometry_three_way.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "run":
        run()
    elif mode == "table":
        table()
    else:
        sys.exit(__doc__)
