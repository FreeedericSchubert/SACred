"""GFN2-xTB hydrogen binding on the MoS2 and silica sites of the support sweep (bio/supports.py).

    python screen/xtb_supports.py run      # relax every structure (resumes where it stopped), then the table
    python screen/xtb_supports.py table    # -> data/xtb_supports.csv
    python screen/xtb_supports.py todo     # the unfinished structures -> data/xtb_supports_todo.json
    python screen/xtb_supports.py shard <k> <n>   # structures k, k+n, ... of that list (one SLURM job on Curta)

MoS2: the metal on three lattice S of the relaxed MoS2 of the DFT runs (dft/runs/mos2/33_mos2_fe without its Fe),
cut from 6x6 to 3x3 for speed, over a Mo (top) and over the empty hexagon (hollow), placed as the Fe of
dft/curta/mos2_fe.py; H on the metal. The Mo layer is fixed, as in DFT. Clean MoS2 with H on a basal S is the
check against DFT (31_mos2_h, 6x6).
Silica: the metal in the open corner of the trisilanol silsesquioxane Me7Si7O9(OH)3, the molecular model of
tripodal grafting on silica: two siloxide O and the remaining silanol O bind the metal, so the cluster is
neutral for M(II). H on the metal, and H on a siloxide O (the support's own proton site) for comparison.
A molecule: no periodic cell, H2 in vacuum.
Settings of screen/xtb_scan.py: spin-restricted, Fermi smearing at 1000 K, BFGS to 0.05 eV/A (<= 400 steps).
dG_H = E(*H) - E(*) - E(H2)/2 + 0.24 eV; on MoS2 the lower of the two bare sites is the reference, as for the DFT
runs. dG_H_offset adds the per-metal offset of H on the metal fitted on graphene (data/xtb_offsets.csv,
screen/calibrate.py; none for Mn, Mo, W): that assumes xTB's error belongs to the metal, not to the support.
"""
import json
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from ase import Atom, Atoms              # noqa: E402
from ase.constraints import FixAtoms     # noqa: E402
from ase.io import read, write           # noqa: E402
from ase.optimize import BFGS            # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "data", "screen", "xtb_supports")
TABLE = os.path.join(ROOT, "data", "xtb_supports.csv")
MOS2_IN = os.path.join(ROOT, "dft", "runs", "mos2", "33_mos2_fe", "relax.in")
TRANSITION_METALS = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
                     "Rh", "Pd", "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]
# first the metals with S3 faces (Fe, Co, Ni, Cu, Zn, Mo, W) or siloxide faces (Mn, Fe, Co, Ni, Cu, Zn) in the sweep,
# now the nine enzyme metals and Cr, Ru, Pt as for the graphene library (screen/xtb_scan.py library); the other
# transition metals are the roadmap
LIBRARY_METALS = ["V", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W", "Cr", "Ru", "Pt"]
MOS2_METALS = LIBRARY_METALS
SILICA_METALS = LIBRARY_METALS
TODO = os.path.join(ROOT, "data", "xtb_supports_todo.json")
CELL_CUT = 2                     # 6x6 DFT sheet -> 3x3 for xTB
VACUUM = 8.0                     # A above and below the sheet (the DFT cell has 35 A for the gate runs)
D_MS = {"Mo": 2.35, "W": 2.35}   # starting M-S (A), else D_MS_DEFAULT
D_MS_DEFAULT = 2.25
D_MO, D_SIO, D_SIC = 1.95, 1.62, 1.86   # starting M-O, Si-O, Si-C (A)
D_MH, D_SH, D_OH, D_CH = 1.55, 1.36, 0.97, 1.09
T8_HALF = 1.55                   # half the Si...Si edge of the T8 cube (A)
WORKERS = int(os.environ.get("SLURM_NTASKS", 13))   # 13 of the 14 local cores, else the job's cores
FMAX, STEPS, ETEMP = 0.05, 400, 1000.0
DG_CORR = 0.24
SCF_RETRY = {"max_iterations": 1000, "mixer_damping": 0.1}   # second attempt after an SCF failure (Mo on MoS2)
DFT_CLEAN_S = 2.0888             # dG_H (eV) of H on a basal S of clean MoS2, DFT step 31 (data/dgh_mos2.csv)


# ---------------------------------------------------------------- MoS2

def mos2():
    """3x3 cut of the relaxed 6x6 sheet: a 6x6 SCF costs ~60 s per xTB step (graphene ~5 s), and 3x3 still folds
    MoS2's K point onto Gamma (4x4 would not); sites 9.5 A apart."""
    a = read(MOS2_IN, format="espresso-in")
    del a[[x.index for x in a if x.symbol not in ("Mo", "S")]]
    frac = (a.get_scaled_positions(wrap=False) + 0.01) % 1.0    # atoms sit on multiples of 1/18: no ties at 0.5
    keep = (frac[:, 0] < 0.5) & (frac[:, 1] < 0.5)
    a = a[keep]
    cell = a.cell.array.copy()
    cell[:2] /= CELL_CUT
    z = a.positions[:, 2]
    cell[2, 2] = z.max() - z.min() + 2 * VACUUM
    a.set_cell(cell)
    a.positions[:, 2] += VACUUM - z.min()
    a.set_initial_magnetic_moments(None)
    a.set_constraint(FixAtoms([x.index for x in a if x.symbol == "Mo"]))
    return a


def mos2_site(metal, site):
    """metal over the Mo next to the top-layer S closest to the cell centre (top) or across that Mo (hollow)."""
    a = mos2()
    z_mo = a.positions[a.symbols == "Mo", 2].mean()
    top_s = [x.index for x in a if x.symbol == "S" and x.position[2] > z_mo]
    centre = a.cell[:2, :2].sum(axis=0) / 2
    probe = Atoms("X", positions=[[*centre, a.positions[top_s[0], 2]]], cell=a.cell, pbc=True)
    both = a + probe
    s = min(top_s, key=lambda i: both.get_distance(len(a), i, mic=True))
    mo = min((x.index for x in a if x.symbol == "Mo"),
             key=lambda i: np.linalg.norm(a.get_distance(i, s, mic=True, vector=True)[:2]))
    v = a.get_distance(mo, s, mic=True, vector=True)[:2]          # Mo -> S, ~1.83 A in plane
    d = D_MS.get(metal, D_MS_DEFAULT)
    z = a.positions[s, 2] + np.sqrt(d ** 2 - v @ v)
    xy = a.positions[mo, :2] if site == "top" else a.positions[mo, :2] - v
    a.append(Atom(metal, (*xy, z)))
    return a, {"imetal": len(a) - 1, "iS": s}


# ---------------------------------------------------------------- silica

def unit(v):
    return v / np.linalg.norm(v)


def methyl(si, out):
    """C and three H of a methyl on Si, pointing along out."""
    c = si + D_SIC * out
    e1 = unit(np.cross(out, [1.0, 0.3, 0.1]))
    e2 = np.cross(out, e1)
    tilt = np.radians(180 - 109.47)
    hs = [c + D_CH * (np.cos(tilt) * out + np.sin(tilt) * (np.cos(p) * e1 + np.sin(p) * e2))
          for p in np.radians([0, 120, 240])]
    return [("C", c)] + [("H", h) for h in hs]


def silica_site(metal, state="bare"):
    """Me7Si7O9(O)2(OH)M: the T8 cube Si8O12 without its (+,+,+) Si, methyls outside, the metal on the three O of
    the open corner (two siloxides and one silanol). state H_M: H on the metal; H_O: H on a siloxide O."""
    s = T8_HALF
    corners = [np.array(c) * s for c in np.array(np.meshgrid([-1, 1], [-1, 1], [-1, 1])).T.reshape(-1, 3)]
    gone = np.array([s, s, s])
    atoms = []
    for c in corners:
        if np.allclose(c, gone):
            continue
        atoms.append(("Si", c))
        atoms += methyl(c, unit(c))
    push = np.sqrt(D_SIO ** 2 - s ** 2)       # bridging O sits outside the edge midpoint
    corner_o = []
    for i, a in enumerate(corners):
        for b in corners[i + 1:]:
            if np.sum(~np.isclose(a, b)) != 1:
                continue
            mid = (a + b) / 2
            out = unit(np.where(np.isclose(a, b), mid, 0.0))
            o = mid + push * out
            atoms.append(("O", o))
            if np.allclose(a, gone) or np.allclose(b, gone):
                corner_o.append(len(atoms) - 1)
    centroid = np.mean([atoms[i][1] for i in corner_o], axis=0)
    r = np.linalg.norm(atoms[corner_o[0]][1] - centroid)
    m = centroid + np.sqrt(D_MO ** 2 - r ** 2) * unit(gone)
    atoms.append((metal, m))
    # the silanol: H on the first corner O, pointing away from the metal and the cage
    o = atoms[corner_o[0]][1]
    atoms.append(("H", o + D_OH * unit(unit(o - m) + unit(o))))
    if state == "H_M":
        atoms.append(("H", m + D_MH * unit(gone)))
    elif state == "H_O":
        o2 = atoms[corner_o[1]][1]
        atoms.append(("H", o2 + D_OH * unit(unit(o2 - m) + unit(o2))))
    mol = Atoms([x for x, _ in atoms], positions=[p for _, p in atoms])
    mol.center(vacuum=8.0)
    imetal = [x.index for x in mol if x.symbol == metal and x.index >= len(mol) - 4][0]
    return mol, {"imetal": imetal, "iO": corner_o}


# ---------------------------------------------------------------- jobs

def jobs():
    j = [("MoS2", "H2", "", "ref"), ("MoS2", "clean", "", "bare"), ("MoS2", "clean", "", "H_S")]
    for m in MOS2_METALS:
        for site in ("top", "hollow"):
            j += [("MoS2", m, site, "bare"), ("MoS2", m, site, "H_M")]
    j += [("SiO2", "H2", "", "ref")]
    for m in SILICA_METALS:
        j += [("SiO2", m, "T7", st) for st in ("bare", "H_M", "H_O")]
    return j


def structure(support, metal, site, state):
    if support == "MoS2":
        if metal == "H2":
            a = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.75]], cell=mos2().cell, pbc=True)
            a.positions += a.cell.sum(axis=0) / 2
            return a, None
        if metal == "clean":
            a = mos2()
            info = {"imetal": None, "iS": None}
            if state == "H_S":
                z_mo = a.positions[a.symbols == "Mo", 2].mean()
                s = max((x.index for x in a if x.symbol == "S" and x.position[2] > z_mo),
                        key=lambda i: -np.linalg.norm(a.positions[i, :2] - a.cell[:2, :2].sum(axis=0) / 2))
                a.append(Atom("H", a.positions[s] + [0, 0, D_SH]))
                info["iS"] = s
            return a, info
        a, info = mos2_site(metal, site)
        if state == "H_M":
            a.append(Atom("H", a.positions[info["imetal"]] + [0, 0, D_MH]))
        return a, info
    if metal == "H2":
        return Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.75]]), None
    return silica_site(metal, state)


def relax(job):
    from tblite.ase import TBLite
    support, metal, site, state = job
    name = "_".join(x for x in (metal, site, state) if x)
    d = os.path.join(OUT, support)
    res_file = os.path.join(d, name + ".json")
    retry = False
    if os.path.exists(res_file):
        res = json.load(open(res_file))
        if res.get("energy") is not None or res.get("retried"):
            return res
        retry = True                         # failed once: again with a damped, longer SCF
    os.makedirs(d, exist_ok=True)
    atoms, info = structure(*job)
    atoms.calc = TBLite(method="GFN2-xTB", electronic_temperature=ETEMP, verbosity=0, **(SCF_RETRY if retry else {}))
    t = time.time()
    res = {"support": support, "metal": metal, "site": site, "state": state, "n_atoms": len(atoms), "retried": retry}
    try:
        opt = BFGS(atoms, logfile=None, maxstep=0.15)
        conv = opt.run(fmax=FMAX, steps=STEPS)
        res.update(energy=float(atoms.get_potential_energy()), steps=opt.nsteps, converged=bool(conv))
    except Exception as exc:                 # SCC failure etc.: record it, the rest goes on
        res.update(energy=None, steps=None, converged=False, error="%s: %s" % (type(exc).__name__, exc))
    res["seconds"] = round(time.time() - t, 1)
    if info is not None and info.get("imetal") is not None and res.get("energy") is not None:
        im = info["imetal"]
        el = "S" if support == "MoS2" else "O"
        near = [x.index for x in atoms if x.symbol == el]
        dist = np.sort(atoms.get_distances(im, near, mic=atoms.pbc.any()))
        res.update(M_X=[round(float(x), 3) for x in dist[:4]])
        if support == "MoS2":
            n_s = int((atoms.symbols == "S").sum())
            z_top = np.sort(atoms.positions[atoms.symbols == "S", 2])[-(n_s // 2):].mean()
            res["M_height"] = round(float(atoms.positions[im, 2] - z_top), 3)
    if state.startswith("H_") and res.get("energy") is not None:
        ih = len(atoms) - 1
        dh = atoms.get_distances(ih, range(ih), mic=atoms.pbc.any())
        j = int(np.argmin(dh))
        res.update(H_partner="M" if info and j == info.get("imetal") else atoms[j].symbol, H_dist=round(float(dh[j]), 3))
    write(os.path.join(d, name + ".xyz"), atoms)
    json.dump(res, open(res_file, "w"), indent=1)
    return res


def todo_list():
    """The structures without a result yet (a failed one counts as done once it was retried)."""
    todo = []
    for job in jobs():
        support, metal, site, state = job
        res = os.path.join(OUT, support, "_".join(x for x in (metal, site, state) if x) + ".json")
        if os.path.exists(res):
            r = json.load(open(res))
            if r.get("energy") is not None or r.get("retried"):
                continue
        todo.append(job)
    return todo


def run(todo=None):
    todo = jobs() if todo is None else todo
    t = time.time()
    # the MoS2 cells are the slow ones: start them first
    with Pool(WORKERS) as pool:
        for k, r in enumerate(pool.imap_unordered(relax, todo), 1):
            print("  %d/%d %s %s %s %s: %s, %.0f s (%.0f min total)" % (
                k, len(todo), r["support"], r["metal"], r["site"], r["state"],
                "E=%.3f" % r["energy"] if r.get("energy") is not None else r.get("error", "failed"),
                r["seconds"], (time.time() - t) / 60), flush=True)


def offsets():
    """Per-metal offsets for H on the metal from screen/calibrate.py (only metals with enough DFT pairs)."""
    o = pd.read_csv(os.path.join(ROOT, "data", "xtb_offsets.csv"))
    return o[o.used & (o.state == "H_M")].set_index("metal").offset


def table():
    rows = [json.load(open(os.path.join(r, f))) for r, _, fs in os.walk(OUT) for f in fs if f.endswith(".json")]
    df = pd.DataFrame(rows)
    off = offsets()
    out = []
    for sup, g in df.groupby("support"):
        e_h2 = g[g.metal == "H2"].energy.iloc[0]
        for metal, m in g[g.metal != "H2"].groupby("metal"):
            bare = m[(m.state == "bare") & m.energy.notna()].sort_values("energy")
            if bare.empty:
                continue
            ref = bare.iloc[0]
            for state, h in m[m.state.str.startswith("H_") & m.energy.notna()].groupby("state"):
                best = h.sort_values("energy").iloc[0]
                de = best.energy - ref.energy - e_h2 / 2
                row = {"support": sup, "metal": metal, "state": state, "site_bare": ref.site, "site_H": best.site,
                       "dE_H": de, "dG_H": de + DG_CORR,
                       "dG_H_offset": de + DG_CORR + off[metal] if state == "H_M" and metal in off else np.nan,
                       "H_partner": best.get("H_partner"), "H_dist": best.get("H_dist"),
                       "M_X_bare": ref.get("M_X"), "M_X_H": best.get("M_X"), "M_height_bare": ref.get("M_height"),
                       "converged": bool(best.converged and ref.converged),
                       "seconds": float(h.seconds.sum() + bare.seconds.sum())}
                if metal == "clean":
                    row["dG_H_dft"] = DFT_CLEAN_S
                out.append(row)
    t = pd.DataFrame(out).sort_values(["support", "state", "metal"])
    t.to_csv(TABLE, index=False, float_format="%.3f")
    print(t.drop(columns=["seconds"]).to_string(index=False, float_format="%.2f"))
    print("-> " + TABLE)
    return t


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "run":
        run()
        table()
    elif mode == "todo":
        json.dump(todo_list(), open(TODO, "w"))
        print("%d unfinished structures -> %s" % (len(json.load(open(TODO))), TODO))
    elif mode == "shard":
        k, n = int(sys.argv[2]), int(sys.argv[3])
        run([tuple(j) for j in json.load(open(TODO))][k::n])
    elif mode == "table":
        table()
    else:
        sys.exit(__doc__)
