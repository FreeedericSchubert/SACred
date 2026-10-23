"""Write the pw.x run directories for a batch of graphene single-atom sites.

    python dft/make_inputs.py graphene1      # dft/runs/graphene/<motif>/{bare,H_M,H_S}/relax.in + meta.json
    python dft/make_inputs.py graphene2      # batch 2 (motifs from the PDB sweep); existing directories are kept
    python dft/make_inputs.py graphene3      # batch 3: the xTB candidates (data/xtb_selected.csv), started from the
                                             # xTB geometries, bare + the H position xTB found near-optimal
    python dft/make_inputs.py library        # the best new sites of the xTB library (data/xtb_selected_library.csv,
                                             # screen/library_next.py), started from the xTB geometries
    python dft/make_inputs.py anchors        # calibration anchors for metals with fewer than two xTB/DFT pairs
                                             # (data/xtb_selected_anchors.csv), from the xTB geometries
    python dft/make_inputs.py stability      # stability references of the sites with |dG_H| <= 0.3 eV so far
                                             # (data/dgh_graphene.csv): the doped vacancy without the metal and
                                             # the bulk metals (dft/runs/bulk/<M>)
    python dft/make_inputs.py conditions     # Co-N2OC under H2-evolution conditions (H on donors, H-capped sites)
    python dft/make_inputs.py neb_points     # Gamma single points on the xTB escape paths (force-field check)
    python dft/make_inputs.py points <file.extxyz> <label> [gamma|331]   # single points on external structures
    python dft/make_inputs.py orr            # O2 cross-check: OH, O, OOH on two Co sites (ORR_SITES) + H2O reference

Each structure gets its own directory (one SLURM job each, see dft/curta/submit.sh). The H2 reference for
the hydrogen adsorption free energy goes to dft/runs/graphene/_ref/H2, in the same cell and settings.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
from ase import Atoms
from ase.build import bulk
from ase.io import read

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "materials", "graphene"))
from batches import GRAPHENE_BATCH1, GRAPHENE_BATCH2
from qe import pw_input, anchor
from sites import build, with_h, h_sites, C_CELL, Z_SHEET

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs")
BATCHES = {"graphene1": ("graphene", GRAPHENE_BATCH1), "graphene2": ("graphene", GRAPHENE_BATCH2)}
# bulk metal references: structure and starting lattice constants (A; relaxed with vc-relax), k mesh. Mn is left
# out until a Mn site needs it (alpha-Mn: 29 atoms per cell, non-collinear magnetism)
BULK = {"V": ("bcc", 3.00, None), "Fe": ("bcc", 2.83, None), "Co": ("hcp", 2.49, 4.02), "Ni": ("fcc", 3.52, None),
        "Cu": ("fcc", 3.63, None), "Zn": ("hcp", 2.66, 4.95), "Mo": ("bcc", 3.15, None), "W": ("bcc", 3.17, None),
        "Cr": ("bcc", 2.88, None), "Ru": ("hcp", 2.71, 4.28), "Pt": ("fcc", 3.92, None)}
KPTS_BULK = {"bcc": (16, 16, 16), "fcc": (16, 16, 16), "hcp": (16, 16, 10)}
WINDOW = 0.3            # eV: sites with |dG_H| up to this get stability references
ORR_SITES = {"Co-N4:DV": "O2 cross-check: the same site in Cho et al. 2022 (PBE-D3, U_ORR 0.95 V)",
             "Co-N2OC:DV:NONC": "O2 cross-check: best new site (xTB U_ORR 0.58 V, DFT dG_H +0.03 eV)"}
ORR_STATES = ("OH", "O", "OOH")


def dirname(motif):
    return motif.replace(":", "_").replace("=", "")


def write_run(path, atoms, meta, fixed=(), calc="relax", **pw):
    """Input + metadata of one structure; an existing directory is never overwritten (it may be running)."""
    if os.path.exists(os.path.join(path, "relax.in")):
        return
    os.makedirs(path, exist_ok=True)
    with open(os.path.join(path, "relax.in"), "w") as f:
        f.write(pw_input(atoms, calc, fixed=fixed, **pw))
    with open(os.path.join(path, "meta.json"), "w") as f:
        json.dump(meta, f, indent=1)


def h2_reference(base, cell):
    h2 = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.75]], cell=cell, pbc=True)
    h2.positions += h2.cell.sum(axis=0) / 2 - [0, 0, C_CELL / 2 - Z_SHEET - 0.4]
    write_run(os.path.join(base, "_ref", "H2"), h2, {"motif": "H2", "state": "ref"})


def from_xtb(base, table="xtb_selected.csv", label="xTB scan"):
    """Batch 3 (and the library): the selected xTB candidates; geometries from data/screen/xtb/<motif>/<state>.xyz."""
    root = os.path.join(os.path.dirname(ROOT), "..", "data")
    sel = pd.read_csv(os.path.join(root, table))
    xtb = os.path.join(root, "screen", "xtb")
    for r in sel.itertuples():
        _, info = build(r.motif)
        d = dirname(r.motif)
        for state in ("bare", r.state):
            atoms = read(os.path.join(xtb, d, state + ".xyz"))
            fixed = [anchor(atoms, info["imetal"])]
            meta = {"motif": r.motif, "enzyme": "%s: %s dG_H %+.2f eV (%s)" % (
                        label, "raw" if pd.isna(r.dG_H_cal) else "calibrated",
                        r.dG_H if pd.isna(r.dG_H_cal) else r.dG_H_cal, r.state),
                    "imetal": info["imetal"], "idonors": info["idonors"], "fixed": fixed, "state": state,
                    "start": "GFN2-xTB geometry"}
            if state != "bare":
                meta["ih"] = len(atoms) - 1
            write_run(os.path.join(base, d, state), atoms, meta, fixed)
    print("%d xTB candidates" % len(sel))


def stability(base):
    """Doped vacancy without the metal (_ref/site_<donors>: same cell, settings and fixed carbon, spin-polarised
    from the donors) for each site near dG_H = 0, and the bulk metals (vc-relax, energy from the relaxation,
    4 cores, 2 h). Metal formation energy and dissolution potential: dft/analyse.py."""
    dgh = pd.read_csv(os.path.join(os.path.dirname(ROOT), "..", "data", "dgh_graphene.csv"))
    new_metal = dgh.motif.str.split("-").str[0].isin(["Cr", "Ru"])    # first E_f of the proof-of-concept metals
    near = sorted({m for m in dgh[(dgh.dG_H.abs() <= WINDOW) | (new_metal & (dgh.max_M_X < 2.6))].motif
                   if ":ax=" not in m})   # axial: own ligand reference
    for motif in near:
        atoms, info = build(motif)
        fixed = [anchor(atoms, info["imetal"])]
        del atoms[info["imetal"]]          # the last atom: the other indices stay
        env = motif.split("-", 1)[1]
        moments = {atoms[i].symbol: 1 for i in info["idonors"] if atoms[i].symbol != "C"}
        meta = {"motif": env, "state": "ref", "from": motif, "idonors": info["idonors"], "fixed": fixed,
                "start": "builder geometry without the metal"}
        path = os.path.join(base, "_ref", "site_" + dirname(env))
        if not os.path.exists(path):
            print("  site_%s (from %s)" % (dirname(env), motif))
        write_run(path, atoms, meta, fixed, moments=moments)
    for metal, (struct, a, c) in BULK.items():
        path = os.path.join(ROOT, "bulk", metal)
        atoms = bulk(metal, struct, a=a, c=c) if c else bulk(metal, struct, a=a)
        write_run(path, atoms, {"motif": metal, "state": "bulk", "structure": struct}, calc="vc-relax",
                  kpts=KPTS_BULK[struct])
        open(os.path.join(path, "gamma_only"), "w").close()       # no separate SCF: energy from the vc-relax
        with open(os.path.join(path, "resources"), "w") as f:
            f.write("NTASKS=4\nMEM=4G\nHOURS=2\n")
    envs = {m.split("-", 1)[1] for m in near}
    print("%d metal-free sites for %d sites with |dG_H| <= %.1f eV, %d bulk metals" % (
        len(envs), len(near), WINDOW, len(BULK)))


def orr(base):
    """O2 cross-check: OH, O and OOH on the metal (dft/runs/graphene/<motif>/ORR_<adsorbate>), started from the
    DFT-relaxed bare site with the adsorbate where GFN2-xTB relaxed it (data/screen/xtb_reactions/graphene/<motif>,
    relative to the metal), and the H2O reference (_ref/H2O) in the same cell. Analysis: dft/orr.py."""
    xtb = os.path.join(os.path.dirname(ROOT), "..", "data", "screen", "xtb_reactions", "graphene")
    cell = None
    for motif, why in ORR_SITES.items():
        d = dirname(motif)
        bare = os.path.join(base, d, "bare")
        meta = json.load(open(os.path.join(bare, "meta.json")))
        atoms = read(os.path.join(bare, "relax.out"), format="espresso-out", index=-1)
        atoms.set_constraint()
        cell = atoms.cell
        im = meta["imetal"]
        for state in ORR_STATES:
            ax = read(os.path.join(xtb, d, state + ".xyz"))
            ads = ax[len(atoms):]
            ads.positions += atoms.positions[im] - ax.positions[im]
            a = atoms + ads
            m = dict(meta, enzyme=why, state="ORR_" + state, iads=list(range(len(atoms), len(a))),
                     start="DFT bare geometry, %s from GFN2-xTB" % state)
            write_run(os.path.join(base, d, "ORR_" + state), a, m, meta["fixed"])
    h2o = Atoms("OH2", positions=[[0, 0, 0], [0.76, 0, 0.59], [-0.76, 0, 0.59]], cell=cell, pbc=True)
    h2o.positions += h2o.cell.sum(axis=0) / 2 - [0, 0, C_CELL / 2 - Z_SHEET - 0.4]
    write_run(os.path.join(base, "_ref", "H2O"), h2o, {"motif": "H2O", "state": "ref"})
    print("%d sites x %d O intermediates + the H2O reference" % (len(ORR_SITES), len(ORR_STATES)))


def final(run):
    atoms = read(os.path.join(run, "relax.out"), format="espresso-out", index=-1)
    atoms.set_constraint()
    return atoms


def with_donor_h(atoms, donors, centre, sides, lean):
    """H on donor atoms (N-H 1.02, O-H 0.98, C-H 1.10 A), out of plane on the given side (+1 above, -1 below the
    sheet) and tilted in plane by `lean` (> 0 away from the site centre, < 0 into the vacancy, where the lone pairs
    of a metal-free site point), so H on neighbouring donors do not clash."""
    bond = {"N": 1.02, "O": 0.98, "C": 1.10, "S": 1.35}
    out = atoms.copy()
    for i, side in zip(donors, sides):
        p = atoms.positions[i]
        away = p - centre
        away[2] = 0.0
        v = lean * away / np.linalg.norm(away) + np.array([0.0, 0.0, side])
        out += Atoms("H", positions=[p + bond[atoms[i].symbol] * v / np.linalg.norm(v)])
    return out


def conditions(base):
    """Co-N2OC under H2-evolution conditions (dft/conditions.py): H on its C and on its O donor (Hdonor_C,
    Hdonor_O), H on the metal with the C donor hydrogenated (CH_HM), and the metal-free site with H on two or on all
    four donors (_ref/site_<env>_H2, _H4), the product of acid demetalation; the same with two H for Co-N4.
    All started from the relaxed DFT structures. The names do not start with H_, so dft/analyse.py leaves them out of
    the dG_H table."""
    site = os.path.join(base, "Co-N2OC_DV_NONC")
    meta = json.load(open(os.path.join(site, "bare", "meta.json")))
    im, don = meta["imetal"], meta["idonors"]                    # N, O, N, C
    bare, hm = final(os.path.join(site, "bare")), final(os.path.join(site, "H_M"))
    centre = bare.positions[im]
    common = {"motif": meta["motif"], "imetal": im, "idonors": don, "fixed": meta["fixed"]}
    why = "operating conditions: %s (Li et al. arXiv 2602.04475; Meng & Zhong, Angew. 2026)"
    for name, atoms, note in (
            ("Hdonor_C", with_donor_h(bare, [don[3]], centre, [1], 0.2), "H on the C donor"),
            ("Hdonor_O", with_donor_h(bare, [don[1]], centre, [1], 0.2), "H on the O donor"),
            ("CH_HM", with_donor_h(hm, [don[3]], centre, [-1], 0.2), "H on Co with the C donor hydrogenated (below)")):
        write_run(os.path.join(site, name), atoms, dict(common, state=name, enzyme=why % note,
                                                        ih=list(range(len(bare), len(atoms))),
                                                        start="relaxed DFT structure + H"), meta["fixed"])
    for env, donors in (("N2OC_DV_NONC", don), ("N4_DV", None)):
        ref = os.path.join(base, "_ref", "site_" + env)
        rmeta = json.load(open(os.path.join(ref, "meta.json")))
        donors = donors or rmeta["idonors"]
        out = sorted(f for f in os.listdir(ref) if f.startswith("relax.out"))
        last = "relax.out" if "relax.out" in out and "End final coordinates" in open(os.path.join(ref, "relax.out")).read() \
            else out[-1]
        empty = read(os.path.join(ref, last), format="espresso-out", index=-1)
        empty.set_constraint()
        centre = empty.positions[donors].mean(axis=0)
        n_ = [i for i in donors if empty[i].symbol == "N"]
        sets = {"H2": ([n_[0], n_[-1]], [1, -1])}
        if env == "N2OC_DV_NONC":
            sets["H4"] = (donors, [1, -1, 1, -1])
        for tag, (idx, sides) in sets.items():
            atoms = with_donor_h(empty, idx, centre, sides, -0.8)
            moments = {atoms[i].symbol: 1 for i in donors if atoms[i].symbol != "C"}
            write_run(os.path.join(base, "_ref", "site_%s_%s" % (env, tag)), atoms,
                      dict(rmeta, state="ref", motif="%s+%s" % (rmeta["motif"], tag), ih=list(range(len(empty), len(atoms))),
                           start="relaxed DFT site (%s) + H on %s" % (last, ", ".join(empty[i].symbol for i in idx))),
                      rmeta["fixed"], moments=moments)
    print("Co-N2OC: 3 runs with H on donors; metal-free sites with H: N2OC (H2, H4), N4 (H2)")


NEB_SITES = ["Zn-NS2:SV", "Ni-NS2:SV"]     # the two escapes that MD with the force field can reach (ps, ns)


def single_point(path, atoms, meta, kpts="gamma", resources="NTASKS=8\nMEM=10G\nHOURS=1\n"):
    """One pw.x single point (no relaxation, no separate SCF step: a gamma_only marker tells step.py the run is
    done once relax.out converged). kpts "gamma" as the relaxations, "331" as the energies behind dG_H."""
    atoms = atoms.copy()
    atoms.set_constraint()
    write_run(path, atoms, meta, (), calc="scf")
    inp = os.path.join(path, "relax.in")
    if kpts == "gamma":
        text = open(inp).read()
        open(inp, "w").write(text[:text.index("K_POINTS automatic")] + "K_POINTS gamma\n")
    open(os.path.join(path, "gamma_only"), "w").close()
    with open(os.path.join(path, "resources"), "w") as f:
        f.write(resources)


def neb_points(base):
    """DFT single points (Gamma point, energies and forces) on every image of the GFN2-xTB escape paths
    (data/screen/xtb_neb/<site>.xyz): <motif>/NEB_<k>. They check the force field and xTB along the path the metal
    takes out of the site, and add off-minimum configurations to the training data. 8 cores, 1 h."""
    xtb = os.path.join(os.path.dirname(ROOT), "..", "data", "screen", "xtb_neb")
    for motif in NEB_SITES:
        d = dirname(motif)
        meta = json.load(open(os.path.join(base, d, "bare", "meta.json")))
        for k, atoms in enumerate(read(os.path.join(xtb, d + ".xyz"), index=":")):
            single_point(os.path.join(base, d, "NEB_%d" % k), atoms,
                         dict(meta, state="NEB_%d" % k, image=k, start="GFN2-xTB NEB image"))
    print("NEB single points for %s" % ", ".join(NEB_SITES))


ESCAPE_SITES = ["Zn-NS2:SV", "Ni-NS2:SV", "Cu-N3:SV"]


def escape_end(base):
    """DFT escape energy: the end state of the GFN2-xTB escape path (metal as an adatom on clean graphene about 5 A
    from its emptied site, the last NEB image) relaxed in DFT -> <motif>/ESC_end. The metal may only move up and
    down (if_pos 0 0 1), so it stays on its hollow site instead of sliding back; the far carbon is fixed as in every
    run. E(ESC_end) - E(bare) is the DFT escape energy, the lower bound of the barrier (dft/conditions.py)."""
    xtb = os.path.join(os.path.dirname(ROOT), "..", "data", "screen", "xtb_neb")
    for motif in ESCAPE_SITES:
        d = dirname(motif)
        meta = json.load(open(os.path.join(base, d, "bare", "meta.json")))
        atoms = read(os.path.join(xtb, d + ".xyz"), index=-1)
        atoms.set_constraint()
        path = os.path.join(base, d, "ESC_end")
        write_run(path, atoms, dict(meta, state="ESC_end", start="last image of the GFN2-xTB escape path"),
                  meta["fixed"])
        inp = os.path.join(path, "relax.in")
        lines = open(inp).read().split("\n")
        start = lines.index("ATOMIC_POSITIONS angstrom") + 1
        k = start + meta["imetal"]
        lines[k] = lines[k].rstrip() + "  0 0 1"
        open(inp, "w").write("\n".join(lines))
        with open(os.path.join(path, "resources"), "w") as f:
            f.write("NTASKS=24\nMEM=24G\nHOURS=6\n")
    print("escape end states for %s" % ", ".join(ESCAPE_SITES))


def points(base, path, label, kpts):
    """DFT single points on structures made elsewhere, e.g. the force-field optimisation with a second sphere:
    every frame of <path> (extended XYZ; info keys motif, e.g. Co-N2OC:DV:NONC, and state, e.g. bare or H_M)
    -> dft/runs/graphene/<motif>/<label>_<state>. Same settings as every other run (dft/qe.py); 24 cores, 3 h.
    Analysis: python dft/points.py (energies, forces left on the structure, dG_H where bare and H_M exist)."""
    n = 0
    for k, atoms in enumerate(read(path, index=":")):
        motif = atoms.info.get("motif")
        state = str(atoms.info.get("state", k))
        if not motif:
            raise SystemExit("frame %d of %s has no 'motif' in its info line" % (k, path))
        meta = {"motif": motif, "state": "%s_%s" % (label, state), "label": label, "kind": state,
                "source": os.path.basename(path), "frame": k, "kpts": kpts,
                "start": "external structure (single point)",
                "info": {key: v for key, v in atoms.info.items() if isinstance(v, (str, int, float))}}
        single_point(os.path.join(base, dirname(motif), "%s_%s" % (label, state)), atoms, meta, kpts,
                     "NTASKS=24\nMEM=24G\nHOURS=3\n")
        n += 1
    print("%d single points (%s, k points %s) from %s" % (n, label, kpts, path))


def main(batch):
    if batch == "escape_end":
        escape_end(os.path.join(ROOT, "graphene"))
        return
    if batch == "points":                    # python dft/make_inputs.py points <file.extxyz> <label> [gamma|331]
        points(os.path.join(ROOT, "graphene"), sys.argv[2], sys.argv[3], sys.argv[4] if len(sys.argv) > 4 else "331")
        return
    if batch == "neb_points":
        neb_points(os.path.join(ROOT, "graphene"))
        return
    if batch == "conditions":
        conditions(os.path.join(ROOT, "graphene"))
        return
    if batch == "orr":
        orr(os.path.join(ROOT, "graphene"))
        return
    if batch == "graphene3":
        from_xtb(os.path.join(ROOT, "graphene"))
        return
    if batch == "library":
        from_xtb(os.path.join(ROOT, "graphene"), "xtb_selected_library.csv", "xTB library")
        return
    if batch == "anchors":
        from_xtb(os.path.join(ROOT, "graphene"), "xtb_selected_anchors.csv", "xTB calibration anchor")
        return
    if batch == "stability":
        stability(os.path.join(ROOT, "graphene"))
        return
    sub, motifs = BATCHES[batch]
    base = os.path.join(ROOT, sub)
    cell = None
    for motif, why in motifs:
        atoms, info = build(motif)
        cell = atoms.cell
        fixed = [anchor(atoms, info["imetal"])]
        common = {"motif": motif, "enzyme": why, "imetal": info["imetal"], "idonors": info["idonors"],
                  "fixed": fixed}
        write_run(os.path.join(base, dirname(motif), "bare"), atoms, dict(common, state="bare"), fixed)
        for site in h_sites(info, atoms):
            ah = with_h(atoms, info, site)
            write_run(os.path.join(base, dirname(motif), "H_" + site), ah,
                      dict(common, state="H_" + site, ih=len(ah) - 1), fixed)
    h2_reference(base, cell)
    runs = sorted(os.path.relpath(r, base) for r, _, fs in os.walk(base) if "relax.in" in fs)
    print("%d run directories under %s" % (len(runs), base))
    for r in runs:
        if r.split(os.sep)[0] in {dirname(m) for m, _ in motifs}:
            print("  " + r)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "graphene1")
