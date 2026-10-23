"""How strongly each donor arrangement holds each metal, with GFN2-xTB: the thermodynamic side of getting single
atoms into carbon instead of metal particles.

    python screen/xtb_stability.py sites     # relax the 96 empty donor arrangements (20 SV, 76 DV); resumes
    python screen/xtb_stability.py adatoms   # clean graphene, and each metal on it (hollow and on top of a C)
    python screen/xtb_stability.py table     # -> data/xtb_stability.csv: insertion and escape energy of every site
    python screen/xtb_stability.py neb [motif ...]   # escape barrier by NEB (default: the sites with DFT E_f)

E_ins = E(metal in the site) - E(empty site), both GFN2-xTB with the settings of screen/xtb_scan.py (same cell, same
fixed carbon). The empty site is the builder geometry without the metal, as the DFT references of
dft/make_inputs.py stability, and it is the same for every metal, so 96 relaxations cover the whole library.
Within one metal E_ins ranks the arrangements directly; across metals it needs the metal's own reference (E_f =
E_ins - E per atom of the bulk metal in DFT): one constant per metal, fitted on the DFT formation energies of
data/stability_graphene.csv where a metal has at least two, else that metal is compared only within itself.

Kinetics, first step: E_escape = [E(metal on clean graphene) - E(clean graphene)] - E_ins, the energy for the metal to
leave its site and sit on intact carbon as an adatom, the start of every route to a particle. It is a lower bound
of the escape barrier. With a 1e13 /s attempt frequency an atom crosses about kT ln(1e13 s^-1 x 3600 s) = 38 kT in an
hour: 3.5 eV at 800 C, 4.2 eV at 1000 C, so sites below that empty during a typical pyrolysis.
Kinetics, second step (a pilot): the barrier itself, by a climbing-image NEB from the metal in its site to the metal
on an intact carbon hexagon 4.3-6 A away (7 images between them, IDPP start, FIRE to 0.1 eV/A). If the barriers sit
at the escape energy, that cheap number is the descriptor; if not, the path has a barrier of its own.
"""
import json
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from ase.constraints import FixAtoms     # noqa: E402
from ase.io import write                 # noqa: E402
from ase.optimize import BFGS            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import xtb_scan as xs                    # noqa: E402

ROOT = xs.ROOT
OUT = os.path.join(ROOT, "data", "screen", "xtb_sites")
ADATOMS = os.path.join(ROOT, "data", "screen", "xtb_adatoms")
NEB_OUT = os.path.join(ROOT, "data", "screen", "xtb_neb")
N_IMAGES, NEB_FMAX, NEB_STEPS = 7, 0.1, 300
TABLE = os.path.join(ROOT, "data", "xtb_stability.csv")
WORKERS = int(os.environ.get("SLURM_NTASKS", 8))      # 8 of the 14 local cores, so the laptop stays usable
MIN_PAIRS = 2


def arrangements():
    """Every arrangement as (site, donor order): 20 sorted triples (SV), 76 orders of d1..d4 up to symmetry (DV)."""
    out = []
    for m in xs.library_motifs():
        if m.startswith(xs.LIBRARY_METALS[0] + "-"):
            out.append(xs.arrangement(m)[1:])
    return out


def name(arr):
    return "%s_%s" % (arr[0], "".join(arr[1]))


def relax(arr):
    """Relax one empty donor arrangement -> <OUT>/<site>_<order>.json/.xyz; skip if done."""
    from tblite.ase import TBLite
    res_file = os.path.join(OUT, name(arr) + ".json")
    if os.path.exists(res_file):
        res = json.load(open(res_file))
        if res.get("energy") is not None:
            return res
    os.makedirs(OUT, exist_ok=True)
    site, order = arr
    donors = list(order)
    motif = "Fe-%s:%s" % (xs.formula(donors), site) + ("" if site == "SV" else ":" + "".join(order))
    atoms, info = xs.build(motif)                         # any metal: the donor positions do not depend on it
    fixed = xs.anchor(atoms, info["imetal"])
    del atoms[info["imetal"]]                             # the last atom: the other indices stay
    atoms.set_constraint(FixAtoms([fixed]))
    atoms.calc = TBLite(method="GFN2-xTB", electronic_temperature=xs.ETEMP, verbosity=0)
    res = {"site": site, "order": "".join(order), "formula": xs.formula(donors), "n_atoms": len(atoms)}
    t = time.time()
    try:
        opt = BFGS(atoms, logfile=None, maxstep=0.15)
        conv = opt.run(fmax=xs.FMAX, steps=xs.STEPS)
        res.update(energy=float(atoms.get_potential_energy()), steps=opt.nsteps, converged=bool(conv))
    except Exception as exc:
        res.update(energy=None, converged=False, error="%s: %s" % (type(exc).__name__, exc))
    res["seconds"] = round(time.time() - t, 1)
    write(os.path.join(OUT, name(arr) + ".xyz"), atoms)
    json.dump(res, open(res_file, "w"), indent=1)
    return res


def run():
    todo = arrangements()
    t = time.time()
    with Pool(WORKERS) as pool:
        for k, r in enumerate(pool.imap_unordered(relax, todo), 1):
            print("  %d/%d %s_%s: %s, %.0f s (%.0f min total)" % (
                k, len(todo), r["site"], r["order"], "E=%.3f" % r["energy"] if r.get("energy") is not None else
                r.get("error", "failed"), r["seconds"], (time.time() - t) / 60), flush=True)


def relax_adatom(job):
    """Clean graphene (metal None) or one metal atom on it, starting over a hexagon (hollow) or a carbon (top)."""
    from tblite.ase import TBLite
    from ase import Atom
    from sites import graphene, neighbours
    metal, start = job
    label = "graphene" if metal is None else "%s_%s" % (metal, start)
    res_file = os.path.join(ADATOMS, label + ".json")
    if os.path.exists(res_file):
        res = json.load(open(res_file))
        if res.get("energy") is not None:
            return res
    os.makedirs(ADATOMS, exist_ok=True)
    atoms = graphene()
    centre = atoms.cell.sum(axis=0) / 2
    i = int(np.argmin(np.linalg.norm(atoms.positions[:, :2] - centre[:2], axis=1)))
    fixed = 0
    if metal is not None:
        j, k = neighbours(atoms, i)[:2]
        p = atoms.positions
        pos = p[i] + [0, 0, 2.0] if start == "top" else p[i] + (p[j] - p[i]) + (p[k] - p[i]) + [0, 0, 1.8]
        atoms.append(Atom(metal, pos))
        fixed = xs.anchor(atoms, len(atoms) - 1)
    atoms.set_constraint(FixAtoms([fixed]))
    atoms.calc = TBLite(method="GFN2-xTB", electronic_temperature=xs.ETEMP, verbosity=0)
    res = {"metal": metal, "start": start, "n_atoms": len(atoms)}
    t = time.time()
    try:
        opt = BFGS(atoms, logfile=None, maxstep=0.15)
        conv = opt.run(fmax=xs.FMAX, steps=xs.STEPS)
        res.update(energy=float(atoms.get_potential_energy()), steps=opt.nsteps, converged=bool(conv))
        if metal is not None:
            res["height"] = round(float(atoms.positions[-1, 2] - np.median(atoms.positions[:-1, 2])), 3)
    except Exception as exc:
        res.update(energy=None, converged=False, error="%s: %s" % (type(exc).__name__, exc))
    res["seconds"] = round(time.time() - t, 1)
    write(os.path.join(ADATOMS, label + ".xyz"), atoms)
    json.dump(res, open(res_file, "w"), indent=1)
    return res


def adatoms():
    todo = [(None, "")] + [(m, s) for m in xs.LIBRARY_METALS for s in ("hollow", "top")]
    with Pool(int(os.environ.get("ADATOM_WORKERS", 3))) as pool:
        for k, r in enumerate(pool.imap_unordered(relax_adatom, todo), 1):
            print("  %d/%d %s %s: %s, %.0f s" % (k, len(todo), r["metal"], r["start"], "E=%.3f" % r["energy"]
                                                if r.get("energy") is not None else r.get("error", "failed"),
                                                r["seconds"]), flush=True)


def away_hollow(atoms, im):
    """In-plane position of an intact carbon hexagon (six C within 1.7 A) 4.3-6 A from the metal, closest to 5 A."""
    from ase import Atom
    from sites import graphene, N_CELL
    g = graphene()
    a1, a2 = g.cell[0] / N_CELL, g.cell[1] / N_CELL                       # primitive lattice vectors
    z = np.median([a.position[2] for a in atoms if a.symbol == "C"])
    best = None
    for i in range(N_CELL):
        for j in range(N_CELL):
            h = i * a1 + j * a2 - (a1 + a2) / 3                            # hexagon centres of the lattice
            probe = atoms.copy()
            probe.set_constraint()
            probe.append(Atom("He", [h[0], h[1], z]))
            vec = probe.get_distances(len(probe) - 1, [im], mic=True, vector=True)[0]
            r = float(np.linalg.norm(vec[:2]))
            dist = probe.get_distances(len(probe) - 1, range(len(atoms)), mic=True)
            ring = np.argsort(dist)[:6]
            intact = dist[ring].max() < 1.7 and all(atoms[k].symbol == "C" for k in ring)
            if 4.3 <= r <= 6.0 and intact and (best is None or abs(r - 5.0) < abs(best[0] - 5.0)):
                best = (r, atoms.positions[im] - vec)
    if best is None:
        return None
    pos = best[1].copy()
    pos[2] = z
    return pos


def escape_neb(motif):
    """escape_neb_once, and again with a damped, longer SCF if the first attempt's SCF fails; a site that fails
    twice is recorded with its error, so one bad site does not stop the others."""
    try:
        return escape_neb_once(motif)
    except Exception:
        try:
            return escape_neb_once(motif, scf=xs_scf_retry())
        except Exception as exc:
            res = {"motif": motif, "metal": motif.split("-")[0], "error": "%s: %s" % (type(exc).__name__, exc)}
            os.makedirs(NEB_OUT, exist_ok=True)
            json.dump(res, open(os.path.join(NEB_OUT, motif.replace(":", "_").replace("=", "") + ".json"), "w"), indent=1)
            return res


def xs_scf_retry():
    return {"max_iterations": 1000, "mixer_damping": 0.1}       # as SCF_RETRY of the reaction and support screens


def escape_neb_once(motif, scf=None):
    """Climbing-image NEB of the metal out of its site -> <NEB_OUT>/<motif>.json (barrier, reaction energy)."""
    from tblite.ase import TBLite
    from ase.io import read
    from ase.mep import NEB
    from ase.optimize import FIRE
    d = motif.replace(":", "_").replace("=", "")
    res_file = os.path.join(NEB_OUT, d + ".json")
    if os.path.exists(res_file):
        return json.load(open(res_file))
    os.makedirs(NEB_OUT, exist_ok=True)
    calc = lambda: TBLite(method="GFN2-xTB", electronic_temperature=xs.ETEMP, verbosity=0, **(scf or {}))   # noqa: E731
    metal = motif.split("-")[0]
    initial = read(os.path.join(xs.OUT, d, "bare.xyz"))
    im = [a.index for a in initial if a.symbol == metal][-1]
    fix = FixAtoms([xs.anchor(initial, im)])
    initial.set_constraint(fix)
    initial.calc = calc()
    e0 = initial.get_potential_energy()
    res = {"motif": motif, "metal": metal}
    t = time.time()
    target = away_hollow(initial, im)
    if target is None:
        res["error"] = "no intact hexagon 4.3-6 A away"
        json.dump(res, open(res_file, "w"), indent=1)
        return res
    final = initial.copy()
    final.positions[im] = target + [0, 0, 1.8]
    final.set_constraint(fix)
    final.calc = calc()
    BFGS(final, logfile=None, maxstep=0.15).run(fmax=xs.FMAX, steps=xs.STEPS)
    back = np.linalg.norm((final.positions[im] - initial.positions[im])[:2])
    res.update(E_final_rel=float(final.get_potential_energy() - e0), final_distance=round(float(back), 2))
    if back < 3.0:                                         # the adatom slid back into the site: no separate minimum
        res["error"] = "the metal returns to the site (%.1f A)" % back
        json.dump(res, open(res_file, "w"), indent=1)
        return res
    images = [initial] + [initial.copy() for _ in range(N_IMAGES)] + [final]
    for image in images[1:-1]:
        image.set_constraint(fix)
        image.calc = calc()
    neb = NEB(images, climb=True)
    neb.interpolate("idpp", mic=True)
    opt = FIRE(neb, logfile=None)
    try:
        conv = opt.run(fmax=NEB_FMAX, steps=NEB_STEPS)
        energies = [image.get_potential_energy() - e0 for image in images]
        res.update(barrier=float(max(energies)), reaction=float(energies[-1]), energies=[round(e, 3) for e in energies],
                   steps=opt.nsteps, converged=bool(conv))
    except Exception as exc:
        res.update(error="%s: %s" % (type(exc).__name__, exc))
    res["seconds"] = round(time.time() - t, 1)
    write(os.path.join(NEB_OUT, d + ".xyz"), images)
    json.dump(res, open(res_file, "w"), indent=1)
    return res


def nebs(motifs):
    with Pool(int(os.environ.get("NEB_WORKERS", WORKERS))) as pool:
        for k, r in enumerate(pool.imap_unordered(escape_neb, motifs), 1):
            print("  %d/%d %s: %s" % (k, len(motifs), r["motif"], "barrier %.2f eV, reaction %+.2f eV (%s)" % (
                r["barrier"], r["reaction"], "converged" if r.get("converged") else "not converged")
                if "barrier" in r else r.get("error")), flush=True)


def adsorption():
    """E(metal on clean graphene) - E(clean graphene) per metal, the lower of the two starts (eV)."""
    rows = [json.load(open(os.path.join(ADATOMS, f))) for f in os.listdir(ADATOMS) if f.endswith(".json")] \
        if os.path.isdir(ADATOMS) else []
    rows = [r for r in rows if r.get("energy") is not None and r.get("converged")]
    clean = [r["energy"] for r in rows if r["metal"] is None]
    if not clean:
        return pd.Series(dtype=float)
    ads = pd.DataFrame([r for r in rows if r["metal"] is not None])
    return (ads.groupby("metal").energy.min() - clean[0]) if len(ads) else pd.Series(dtype=float)


def table():
    """E_ins of every xTB site with a relaxed bare structure (scan and library), and E_f where DFT has one."""
    empty = {}
    for f in os.listdir(OUT):
        if f.endswith(".json"):
            r = json.load(open(os.path.join(OUT, f)))
            if r.get("energy") is not None and r.get("converged"):
                empty[(r["site"], tuple(r["order"]))] = r["energy"]
    rows = []
    for d in os.listdir(xs.OUT):
        res = os.path.join(xs.OUT, d, "bare.json")
        if not os.path.exists(res):
            continue
        r = json.load(open(res))
        if r.get("energy") is None or r.get("motif") is None:
            continue
        try:
            key = xs.arrangement(r["motif"])
        except (ValueError, KeyError):
            continue
        if key is None or key[1:] not in empty:
            continue
        rows.append({"motif": r["motif"], "metal": key[0], "site": key[1], "order": "".join(key[2]),
                     "formula": xs.formula(list(key[2])), "E_ins": r["energy"] - empty[key[1:]],
                     "converged": r.get("converged"), "max_M_X": r.get("max_M_X"), "M_height": r.get("M_height")})
    t = pd.DataFrame(rows).sort_values("motif").drop_duplicates(["metal", "site", "order"])   # one name per site
    for el in xs.DONORS:
        t["n_" + el] = t.order.str.count(el)
    # the metal's reference from DFT: E_f(DFT) - E_ins(xTB) per metal, used from MIN_PAIRS pairs on
    dft = pd.read_csv(os.path.join(ROOT, "data", "stability_graphene.csv")).dropna(subset=["E_f"])
    dft["key"] = dft.motif.map(lambda m: xs.arrangement(m))
    t["key"] = [(m, s, tuple(o)) for m, s, o in zip(t.metal, t.site, t.order)]
    pairs = dft.merge(t[["key", "E_ins"]], on="key")
    pairs["metal"] = pairs.motif.str.split("-").str[0]
    ref = (pairs.E_f - pairs.E_ins).groupby(pairs.metal).agg(["mean", "size"])
    ref = ref[ref["size"] >= MIN_PAIRS]["mean"]
    t["E_f_xtb"] = t.E_ins + t.metal.map(ref)
    t["E_escape"] = t.metal.map(adsorption()) - t.E_ins          # lower bound of the barrier out of the site
    t = t.drop(columns="key").sort_values(["metal", "E_ins"])
    t.to_csv(TABLE, index=False, float_format="%.4f")
    print("%d sites with E_ins (%d metals); metals with a DFT reference: %s -> %s" % (
        len(t), t.metal.nunique(), ", ".join("%s (%d pairs)" % (m, n) for m, n in
                                             pairs.groupby("metal").size().items() if m in ref), TABLE))
    print(pairs[["motif", "E_f", "E_ins"]].assign(E_f_xtb=lambda p: p.E_ins + p.motif.str.split("-").str[0].map(ref))
          .round(2).to_string(index=False))
    return t, pairs


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "sites":
        run()
    elif mode == "adatoms":
        adatoms()
    elif mode == "neb":
        motifs = sys.argv[2:] or sorted(pd.read_csv(os.path.join(ROOT, "data", "stability_graphene.csv"))
                                        .dropna(subset=["E_f"]).motif)
        nebs(motifs)
    elif mode == "table":
        table()
    else:
        sys.exit(__doc__)
