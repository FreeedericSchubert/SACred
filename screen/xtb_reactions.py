"""GFN2-xTB intermediates of O2, CO2 and N2 electrochemistry on the sites we already have (hydrogen evolution is
screen/xtb_scan.py and screen/xtb_supports.py), plus an O2 calibration against published DFT.

    python screen/xtb_reactions.py run       # relax every structure (resumes where it stopped), then the tables
    python screen/xtb_reactions.py table     # -> data/xtb_reactions.csv, data/xtb_orr_calibration.csv
    python screen/xtb_reactions.py retry     # every structure whose SCF failed and was not retried yet
    python screen/xtb_reactions.py todo      # the unfinished structures -> data/xtb_reactions_todo.json
    python screen/xtb_reactions.py shard <k> <n> [list]   # structures k, k+n, ... of that list (one SLURM job on
                                             # Curta); list: another job list, e.g. data/xtb_reactions_library_todo.json

Sites: the 58 graphene SAC models with DFT dG_H (data/dgh_graphene.csv), M-S3 on MoS2 (the lower bare site of
screen/xtb_supports.py) and M on the silsesquioxane silica model. Calibration: the double-vacancy sites of Cho et al.,
Front. Chem. 10, 873609 (2022) that the graphene builder can make (N/O/S donors), with their PBE-D3 reaction energies
from Catalysis-Hub (data/external/catalysis_hub_cho2022_orr.csv).
Adsorbates on the metal, on the side away from the support, starting from the relaxed bare site: OH, O, OOH
(O2 reduction to water or H2O2, O2 evolution), COOH, CO, OCHO (CO2 reduction to CO or formate), N2 end-on, N2H, NH2,
NH3 (N2 reduction). Settings of screen/xtb_scan.py (spin-restricted, Fermi smearing 1000 K, BFGS to 0.05 eV/A).
Free energies, computational hydrogen electrode at 0 V vs RHE: G = E(xTB) + ZPE - TS with typical literature values
(ZPE_TS, assumed, not computed). Gas molecules from xTB: H2, H2O, CO2 and N2, each in the support's own cell; O2 enters
as 2 H2O - 2 H2 + 4.92 eV, CO, HCOOH and NH3 through their experimental reaction free energies (GAS_ANCHOR), the usual
way to keep gas-phase errors out of the comparison. dG_*_cal: O-intermediates after the linear xTB -> DFT fit on the
Cho sites; CO2 and N2 intermediates have no DFT reference and stay uncalibrated.
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

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path[:0] = [os.path.join(ROOT, "materials", "graphene"), os.path.join(ROOT, "dft"), HERE]
from qe import anchor                    # noqa: E402
from sites import build                  # noqa: E402

DATA = os.path.join(ROOT, "data")
OUT = os.path.join(DATA, "screen", "xtb_reactions")
SCAN = os.path.join(DATA, "screen", "xtb")
SUPPORTS = os.path.join(DATA, "screen", "xtb_supports")
WORKERS = int(os.environ.get("SLURM_NTASKS", 13))   # 13 of the 14 local cores, else the job's cores
FMAX, STEPS, ETEMP = 0.05, 300, 1000.0
SCF_RETRY = {"max_iterations": 1000, "mixer_damping": 0.1}
ADSORBATES = ["OH", "O", "OOH", "COOH", "CO", "OCHO", "N2", "N2H", "NH2", "NH3"]
O_STATES = ["OH", "O", "OOH"]
CALIBRATION_METALS = ["Co", "Cr", "Fe", "Mn", "Ni", "V", "Zn"]   # Cho's metals the builder knows
# ZPE - TS (eV, 298 K), typical values used in SAC screening papers (adsorbates: vibrations only; gases: standard
# entropies; H2O at 0.035 bar, the vapour pressure that makes it equivalent to liquid water)
ZPE_TS = {"H": 0.16, "OH": 0.29, "O": 0.04, "OOH": 0.34, "COOH": 0.42, "CO": 0.10, "OCHO": 0.44,
          "N2": 0.12, "N2H": 0.37, "NH2": 0.58, "NH3": 0.83,
          "gas_H2": -0.13, "gas_H2O": -0.11, "gas_CO2": -0.35, "gas_N2": -0.44}
GAS_ANCHOR = {"O2_4e": 4.92,         # 2 H2O -> O2 + 2 H2: 4 x 1.23 eV
              "CO": 0.20,            # CO2 + H2 -> CO + H2O(l), CO2/CO at -0.10 V vs RHE
              "HCOOH": 0.24,         # CO2 + H2 -> HCOOH(aq), -0.12 V vs RHE
              "NH3": -0.17}          # 1/2 N2 + 3/2 H2 -> NH3(g)
BOND = {"O": 1.85, "C": 1.90, "N": 1.90}
MIN_CALIBRATION = 5             # sites per O intermediate before the xTB -> DFT line is used


# ---------------------------------------------------------------- adsorbates

def frame(up):
    up = up / np.linalg.norm(up)
    e1 = np.cross(up, [0.3, 1.0, 0.2])
    e1 /= np.linalg.norm(e1)
    return up, e1, np.cross(up, e1)


def toward(up, e1, deg):
    """Unit vector at `deg` from up, towards e1 (negative: towards -e1)."""
    a = np.radians(deg)
    return np.cos(a) * up + np.sin(a) * e1


def adsorbate(name, m, up):
    """Atoms of one adsorbate above metal position m, with `up` pointing away from the support."""
    up, e1, e2 = frame(up)
    u = lambda deg: toward(up, e1, deg)     # noqa: E731
    if name == "OH":
        o = m + BOND["O"] * up
        return [("O", o), ("H", o + 0.98 * u(70))]
    if name == "O":
        return [("O", m + 1.65 * up)]
    if name == "OOH":
        o1 = m + 1.90 * up
        o2 = o1 + 1.45 * u(60)
        return [("O", o1), ("O", o2), ("H", o2 + 0.98 * toward(up, -e1, 30))]
    if name == "COOH":
        c = m + 1.95 * up
        o1, o2 = c + 1.22 * u(60), c + 1.35 * u(-60)
        return [("C", c), ("O", o1), ("O", o2), ("H", o2 + 0.98 * u(-10))]
    if name == "CO":
        c = m + 1.85 * up
        return [("C", c), ("O", c + 1.15 * up)]
    if name == "OCHO":
        o1 = m + 2.00 * up
        c = o1 + 1.27 * u(50)
        return [("O", o1), ("C", c), ("O", c + 1.24 * u(-10)), ("H", c + 1.10 * u(110))]
    if name == "N2":
        n1 = m + 1.90 * up
        return [("N", n1), ("N", n1 + 1.12 * up)]
    if name == "N2H":
        n1 = m + 1.85 * up
        n2 = n1 + 1.22 * u(20)
        return [("N", n1), ("N", n2), ("H", n2 + 1.03 * u(-60))]
    if name == "NH2":
        n = m + 1.95 * up
        return [("N", n), ("H", n + 1.02 * (0.45 * up + 0.89 * e1)), ("H", n + 1.02 * (0.45 * up - 0.45 * e1 + 0.77 * e2))]
    if name == "NH3":
        n = m + 2.05 * up
        return [("N", n)] + [("H", n + 1.02 * (0.33 * up + 0.94 * (np.cos(p) * e1 + np.sin(p) * e2)))
                             for p in np.radians([0, 120, 240])]
    raise ValueError(name)


# ---------------------------------------------------------------- sites

def graphene_motifs():
    d = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv"))
    return sorted(d.motif.unique())


def calibration_motifs():
    ext = pd.read_csv(os.path.join(DATA, "external", "catalysis_hub_cho2022_orr.csv"), comment="#")
    keep = ext[ext.metal.isin(CALIBRATION_METALS) & ~ext.donors.str.contains("B|P")]
    return sorted({"%s-%s:DV" % (m, d) for m, d in zip(keep.metal, keep.donors)})


def mos2_sites():
    t = pd.read_csv(os.path.join(DATA, "xtb_supports.csv"))
    t = t[(t.support == "MoS2") & (t.state == "H_M") & (t.metal != "clean")]
    return {m: s for m, s in zip(t.metal, t.site_bare)}


def bare_site(kind, name):
    """(atoms with constraints, metal index, up vector, fixed description) of a bare site."""
    if kind == "graphene":
        path = os.path.join(SCAN, name.replace(":", "_").replace("=", ""), "bare.xyz")
        atoms, info = build(name)
        if os.path.exists(path):                       # the relaxed bare site of the hydrogen scan
            relaxed = read(path)
            if len(relaxed) == len(atoms) and (relaxed.symbols == atoms.symbols).all():
                atoms.positions = relaxed.positions
        im = info["imetal"]
        atoms.set_constraint(FixAtoms([anchor(atoms, im)]))
        return atoms, im, np.array([0, 0, 1.0])
    if kind == "MoS2":
        atoms = read(os.path.join(SUPPORTS, "MoS2", "%s_%s_bare.xyz" % (name, mos2_sites()[name])))
        atoms.set_constraint(FixAtoms([a.index for a in atoms if a.symbol == "Mo"]))
        im = [a.index for a in atoms if a.symbol == name][-1]
        return atoms, im, np.array([0, 0, 1.0])
    atoms = read(os.path.join(SUPPORTS, "SiO2", "%s_T7_bare.xyz" % name))
    atoms.pbc = False
    im = [a.index for a in atoms if a.symbol == name][-1]
    cage = atoms.positions[[a.index for a in atoms if a.symbol == "Si"]].mean(axis=0)
    return atoms, im, atoms.positions[im] - cage


def gas(kind, molecule):
    """A gas molecule in the support's own cell (graphene, MoS2) or in vacuum (silica)."""
    shapes = {"H2": [("H", (0, 0, 0)), ("H", (0, 0, 0.75))],
              "H2O": [("O", (0, 0, 0)), ("H", (0.76, 0, 0.59)), ("H", (-0.76, 0, 0.59))],
              "CO2": [("O", (0, 0, -1.16)), ("C", (0, 0, 0)), ("O", (0, 0, 1.16))],
              "N2": [("N", (0, 0, 0)), ("N", (0, 0, 1.10))]}
    mol = Atoms([s for s, _ in shapes[molecule]], positions=[p for _, p in shapes[molecule]])
    if kind == "SiO2":
        mol.center(vacuum=8.0)
        return mol
    ref, _, _ = bare_site(kind, "Co-N4:DV" if kind == "graphene" else "Co")
    mol.set_cell(ref.cell)
    mol.pbc = True
    mol.positions += ref.cell.sum(axis=0) / 2
    return mol


def jobs():
    j = []
    for kind in ("graphene", "MoS2", "SiO2"):
        j += [(kind, "gas", m) for m in ("H2", "H2O", "CO2", "N2")]
    for m in graphene_motifs():
        j += [("graphene", m, "bare")] + [("graphene", m, a) for a in ADSORBATES]
    for m in calibration_motifs():
        if m not in graphene_motifs():
            j += [("graphene", m, "bare")] + [("graphene", m, a) for a in O_STATES]
    for m in mos2_sites():
        j += [("MoS2", m, "bare")] + [("MoS2", m, a) for a in ADSORBATES]
    for m in ("Mn", "Fe", "Co", "Ni", "Cu", "Zn"):
        j += [("SiO2", m, "bare")] + [("SiO2", m, a) for a in ADSORBATES]
    return j


def structure(kind, name, state):
    if name == "gas":
        return gas(kind, state), None
    atoms, im, up = bare_site(kind, name)
    if state != "bare":
        for sym, pos in adsorbate(state, atoms.positions[im], up):
            atoms.append(Atom(sym, pos))
    return atoms, im


def relax(job):
    from tblite.ase import TBLite
    kind, name, state = job
    d = os.path.join(OUT, kind, name.replace(":", "_").replace("=", ""))
    res_file = os.path.join(d, state + ".json")
    retry = False
    if os.path.exists(res_file):
        res = json.load(open(res_file))
        if res.get("energy") is not None or res.get("retried"):
            return res
        retry = True
    os.makedirs(d, exist_ok=True)
    try:
        atoms, im = structure(kind, name, state)
    except Exception as exc:
        res = {"kind": kind, "site": name, "state": state, "energy": None, "error": "build: %s" % exc, "retried": True}
        json.dump(res, open(res_file, "w"))
        return res
    n_bare = len(atoms) - (0 if state in ("bare", "gas") or name == "gas" else len(adsorbate(state, np.zeros(3), np.array([0, 0, 1.0]))))
    atoms.calc = TBLite(method="GFN2-xTB", electronic_temperature=ETEMP, verbosity=0, **(SCF_RETRY if retry else {}))
    t = time.time()
    res = {"kind": kind, "site": name, "state": state, "n_atoms": len(atoms), "retried": retry}
    try:
        opt = BFGS(atoms, logfile=None, maxstep=0.15)
        conv = opt.run(fmax=FMAX, steps=STEPS)
        res.update(energy=float(atoms.get_potential_energy()), steps=opt.nsteps, converged=bool(conv))
    except Exception as exc:
        res.update(energy=None, steps=None, converged=False, error="%s: %s" % (type(exc).__name__, exc))
    res["seconds"] = round(time.time() - t, 1)
    if im is not None and state not in ("bare",) and res.get("energy") is not None:
        ads = list(range(n_bare, len(atoms)))
        heavy = [i for i in ads if atoms[i].symbol != "H"]
        dist = atoms.get_distances(im, heavy, mic=atoms.pbc.any())
        res["M_X"] = round(float(dist.min()), 3)            # metal to the nearest adsorbate atom
        res["bound"] = bool(dist.min() < 2.6)               # still on the metal (not desorbed or moved off)
    write(os.path.join(d, state + ".xyz"), atoms)
    json.dump(res, open(res_file, "w"), indent=1)
    return res


def todo_list():
    """Every structure without a result yet (a failed one counts as done once it was retried)."""
    order = {"gas": 0, "bare": 1, "OH": 2, "O": 2, "OOH": 2, "COOH": 3, "CO": 3, "OCHO": 3}
    df = load()
    done = set() if df.empty else {(r.kind, r.site, r.state) for r in df.itertuples()
                                   if pd.notna(r.energy) or r.retried == True}   # noqa: E712 (NaN-safe)
    todo = [j for j in jobs() if j not in done]
    return sorted(todo, key=lambda j: (0 if j[1] == "gas" else order.get(j[2], 4), j[0] != "SiO2", j[0] != "MoS2"))


def run(todo=None):
    todo = todo_list() if todo is None else todo
    t = time.time()
    with Pool(WORKERS) as pool:
        for k, r in enumerate(pool.imap_unordered(relax, todo), 1):
            if k % 20 == 0 or k == len(todo):
                print("  %d/%d structures, %.0f min" % (k, len(todo), (time.time() - t) / 60), flush=True)


# ---------------------------------------------------------------- free energies

def load():
    """Every finished structure; a file still being written by a running worker is skipped."""
    rows = []
    for r, _, fs in os.walk(OUT):
        for f in fs:
            if f.endswith(".json"):
                try:
                    rows.append(json.load(open(os.path.join(r, f))))
                except json.JSONDecodeError:
                    pass
    return pd.DataFrame(rows)


def gases(df):
    """Per support kind: electronic energies E of H2, H2O, CO2, N2 and free energies G including the anchored
    molecules (O2 enters only through the 4.92 eV of the O2 steps)."""
    out = {}
    for kind, g in df[df.site == "gas"].groupby("kind"):
        E = g.set_index("state").energy.to_dict()
        G = {m: E[m] + ZPE_TS["gas_" + m] for m in ("H2", "H2O", "CO2", "N2")}
        G["CO"] = G["CO2"] + G["H2"] - G["H2O"] + GAS_ANCHOR["CO"]
        G["HCOOH"] = G["CO2"] + G["H2"] + GAS_ANCHOR["HCOOH"]
        G["NH3"] = 0.5 * G["N2"] + 1.5 * G["H2"] + GAS_ANCHOR["NH3"]
        out[kind] = (E, G)
    return out


def references(X):
    """What each adsorbate takes from the gas phase (X = E or G of the gas molecules)."""
    return {"OH": X["H2O"] - 0.5 * X["H2"], "O": X["H2O"] - X["H2"], "OOH": 2 * X["H2O"] - 1.5 * X["H2"],
            "COOH": X["CO2"] + 0.5 * X["H2"], "CO": X.get("CO"), "OCHO": X["CO2"] + 0.5 * X["H2"],
            "N2": X["N2"], "N2H": X["N2"] + 0.5 * X["H2"], "NH2": 0.5 * X["N2"] + X["H2"],
            "NH3": 0.5 * X["N2"] + 1.5 * X["H2"]}


def binding(df):
    """One row per site: dE (O species, on the Catalysis-Hub scale: vs H2O and H2) and dG of every adsorbate."""
    gas = gases(df)
    rows = []
    for (kind, site), s in df[df.site != "gas"].groupby(["kind", "site"]):
        e = s.set_index("state")
        if pd.isna(e.energy.get("bare")):
            continue
        E, G = gas[kind]
        refE, refG = references(E), references(G)
        row = {"kind": kind, "site": site}
        for a in ADSORBATES:
            if pd.isna(e.energy.get(a)):
                continue
            ea = e.energy[a] - e.energy["bare"]
            if a in O_STATES:
                row["dE_" + a] = ea - refE[a]
            row["dG_" + a] = ea + ZPE_TS[a] - refG[a]
            row["ok_" + a] = bool(e.converged.get(a)) and bool(e.get("bound", pd.Series(dtype=bool)).get(a, True))
        rows.append(row)
    return pd.DataFrame(rows)


def calibrate(b):
    """dE(DFT) = k dE(xTB) + c per O intermediate, on the double-vacancy sites of Cho et al. 2022 that we built."""
    ext = pd.read_csv(os.path.join(DATA, "external", "catalysis_hub_cho2022_orr.csv"), comment="#")
    ext["site"] = ext.metal + "-" + ext.donors + ":DV"
    g = b[b.kind == "graphene"].set_index("site")
    pairs = []
    for r in ext.itertuples():
        col = "dE_" + r.adsorbate
        if r.site in g.index and col in g and pd.notna(g.loc[r.site, col]) and g.loc[r.site, "ok_" + r.adsorbate]:
            pairs.append({"site": r.site, "adsorbate": r.adsorbate, "dE_dft": r.dE_eV, "dE_xtb": g.loc[r.site, col]})
    p = pd.DataFrame(pairs)
    fit = {}
    for a, q in (p.groupby("adsorbate") if len(p) else []):
        if len(q) < MIN_CALIBRATION:      # a line through a handful of points calibrates nothing: stay raw
            continue
        k, c = np.polyfit(q.dE_xtb, q.dE_dft, 1)
        fit[a] = {"k": k, "c": c, "r": float(np.corrcoef(q.dE_xtb, q.dE_dft)[0, 1]), "n": len(q),
                  "mae_raw": float((q.dE_xtb - q.dE_dft).abs().mean()),
                  "mae_fit": float((k * q.dE_xtb + c - q.dE_dft).abs().mean())}
    return p, fit


def limiting(row, tag=""):
    """Limiting potentials (V vs RHE) at pH 0 from the free energies of one site. An adsorbate that left the metal
    or did not converge (ok_* False) counts as missing."""
    out = {}
    bad = {k for k, v in row.items() if k.startswith("ok_") and v is not None and v == False}   # noqa: E712 (np.bool_, NaN)
    row = {k: (None if k.startswith("dG_") and "ok_" + k[3:].replace("_cal", "") in bad else v) for k, v in row.items()}
    oh, o, ooh = (row.get("dG_%s%s" % (a, tag)) for a in ("OH", "O", "OOH"))
    if all(v is not None and pd.notna(v) for v in (oh, o, ooh)):
        steps = [ooh - 4.92, o - ooh, oh - o, -oh]                  # O2 -> *OOH -> *O -> *OH -> H2O
        out["U_ORR" + tag] = -max(steps)
        out["U_OER" + tag] = max(-x for x in steps)                 # the same steps backwards
        out["U_H2O2" + tag] = min(4.92 - ooh, ooh - 3.52)           # O2 -> *OOH -> H2O2 (equilibrium 0.70 V)
    if tag:
        return out
    cooh, co, ocho = row.get("dG_COOH"), row.get("dG_CO"), row.get("dG_OCHO")
    if all(v is not None and pd.notna(v) for v in (cooh, co)):
        out["U_CO"] = -max(cooh, co - cooh + GAS_ANCHOR["CO"])      # CO2 -> *COOH -> *CO + H2O
        out["dG_CO_release"] = -co                                  # *CO -> CO(g); > ~0.5 eV: CO poisons the site
    if ocho is not None and pd.notna(ocho):
        out["U_HCOOH"] = -max(ocho, GAS_ANCHOR["HCOOH"] - ocho)     # CO2 -> *OCHO -> HCOOH
    n2, n2h, nh2, nh3 = (row.get("dG_" + a) for a in ("N2", "N2H", "NH2", "NH3"))
    if all(v is not None and pd.notna(v) for v in (n2, n2h, nh2, nh3)):
        out["dG_N2_N2H"] = n2h - n2                                 # first protonation
        out["dG_NH2_NH3"] = nh3 - nh2                               # last protonation
        u = -max(n2h - n2, nh3 - nh2)
        # only the first and the last step are computed: above the equilibrium potential (NH3: +0.06 V) one of
        # the middle steps must limit, so such a value says nothing about the site
        out["U_NRR"] = u if u <= -GAS_ANCHOR["NH3"] / 3 else np.nan
    return out


def her_values():
    """dG_H (H on the metal) of the same sites, for selectivity against H2 evolution: DFT where there is a value,
    else GFN2-xTB with the per-metal offset of screen/calibrate.py, else raw xTB -> {(kind, site or metal): (dG_H,
    source)}."""
    out = {}
    off = pd.read_csv(os.path.join(DATA, "xtb_offsets.csv"))
    off = off[off.used & (off.state == "H_M")].set_index("metal").offset
    x = pd.read_csv(os.path.join(DATA, "xtb_scan.csv"))
    for m, v in x[x.state == "H_M"].groupby("motif").dG_H.min().items():
        o = off.get(m.split("-")[0])
        out[("graphene", m)] = (v + o, "xTB calibrated") if pd.notna(o) else (v, "xTB raw")
    d = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv"))
    for m, v in d[d.state == "H_M"].dropna(subset=["dG_H"]).groupby("motif").dG_H.min().items():
        out[("graphene", m)] = (v, "DFT")
    s = pd.read_csv(os.path.join(DATA, "xtb_supports.csv"))
    for r in s[s.state == "H_M"].itertuples():
        out[(r.support, r.metal)] = ((r.dG_H_offset, "xTB calibrated") if pd.notna(r.dG_H_offset) else
                                     (r.dG_H, "xTB raw"))
    fe = pd.read_csv(os.path.join(DATA, "dgh_mos2.csv")).set_index("run").dG_H.get("34_mos2_fe_h")   # Fe-S3 on MoS2
    if fe is not None and pd.notna(fe):
        out[("MoS2", "Fe")] = (fe, "DFT")
    return out


def table():
    df = load()
    b = binding(df)
    pairs, fit = calibrate(b)
    pairs.to_csv(os.path.join(DATA, "xtb_orr_calibration.csv"), index=False, float_format="%.3f")
    for a, f in fit.items():
        print("O calibration %-3s: dE_DFT = %.2f dE_xTB %+.2f eV, r = %.2f, MAE %.2f eV raw -> %.2f eV fitted "
              "(%d sites)" % (a, f["k"], f["c"], f["r"], f["mae_raw"], f["mae_fit"], f["n"]))
    gas = gases(df)
    her = her_values()
    calibration_only = set(calibration_motifs()) - set(graphene_motifs())   # Cho's sites, only for the O fit
    rows = []
    for r in b.to_dict("records"):
        E, G = gas[r["kind"]]
        for a in O_STATES:          # calibrated O species: fitted dE, then the same ZPE - TS as the raw value
            if pd.notna(r.get("dE_" + a)) and a in fit:
                f = fit[a]
                r["dG_%s_cal" % a] = r["dG_" + a] + (f["k"] * r["dE_" + a] + f["c"]) - r["dE_" + a]
        r.update(limiting(r))
        r.update(limiting(r, "_cal"))
        metal = r["site"].split("-")[0]
        r["metal"] = metal
        r["dG_H"], r["dG_H_source"] = her.get((r["kind"], r["site"] if r["kind"] == "graphene" else metal), (None, None))
        r["calibration_only"] = r["kind"] == "graphene" and r["site"] in calibration_only
        rows.append(r)
    out = pd.DataFrame(rows)
    lead = ["kind", "site", "metal", "calibration_only", "dG_H", "dG_H_source", "U_ORR_cal", "U_OER_cal", "U_H2O2_cal", "U_CO",
            "dG_CO_release", "U_HCOOH", "U_NRR", "dG_N2_N2H", "dG_NH2_NH3", "U_ORR", "U_OER", "U_H2O2"]
    rest = sorted(c for c in out.columns if c not in lead)
    out = out[[c for c in lead if c in out.columns] + rest]
    out.to_csv(os.path.join(DATA, "xtb_reactions.csv"), index=False, float_format="%.3f")
    json.dump({a: {k: round(v, 4) if isinstance(v, float) else v for k, v in f.items()} for a, f in fit.items()},
              open(os.path.join(DATA, "xtb_orr_calibration_fit.json"), "w"), indent=1)
    show = out[~out.calibration_only]
    for col, best in (("U_ORR_cal", "max"), ("U_OER_cal", "min"), ("U_CO", "max"), ("U_HCOOH", "max"), ("U_NRR", "max")):
        if col not in show:
            continue
        top = show.dropna(subset=[col]).sort_values(col, ascending=(best == "min")).head(6)
        print("\n%s (%s is best):" % (col, best))
        print(top[["kind", "site", col, "dG_H"]].round(2).to_string(index=False))
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "run":
        run()
    elif mode == "todo":
        todo = todo_list()
        json.dump(todo, open(os.path.join(DATA, "xtb_reactions_todo.json"), "w"))
        print("%d unfinished structures -> data/xtb_reactions_todo.json" % len(todo))
    elif mode == "shard":
        k, n = int(sys.argv[2]), int(sys.argv[3])
        todo = sys.argv[4] if len(sys.argv) > 4 else os.path.join(DATA, "xtb_reactions_todo.json")
        run([tuple(j) for j in json.load(open(todo))][k::n])
    elif mode == "retry":
        failed = [(r["kind"], r["site"], r["state"]) for r in load().to_dict("records")
                  if pd.isna(r.get("energy")) and r.get("retried") != True]   # noqa: E712 (NaN-safe)
        print("%d structures with an SCF failure, once more with the damped SCF" % len(failed))
        run(failed)
    elif mode == "table":
        table()
    else:
        sys.exit(__doc__)
