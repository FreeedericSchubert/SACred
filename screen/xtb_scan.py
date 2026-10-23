"""Fast first scan of graphene single-atom sites with GFN2-xTB (tblite, periodic, Gamma point), before DFT.

    python screen/xtb_scan.py bench           # the 28 DFT motifs of dft/batches.py (calibration against DFT)
    python screen/xtb_scan.py scan            # enzyme donor environments x 9 enzyme metals
    python screen/xtb_scan.py table           # -> data/xtb_scan.csv (dE_H, dG_H per site and H position)
    python screen/xtb_scan.py jobs            # the scan's structure list -> data/xtb_scan_jobs.json
    python screen/xtb_scan.py shard <k> <n>   # structures k, k+n, ... of that list (one SLURM job on Curta)
    python screen/xtb_scan.py library         # every donor arrangement x the 9 enzyme metals + Cr, Ru, Pt, the sites
                                              # not yet computed -> data/xtb_library_jobs.json + starting structures
    python screen/xtb_scan.py library-shard <k> <n>   # structures k, k+n, ... of the library list
    python screen/xtb_scan.py retry           # every structure whose SCF failed, once more with a damped SCF

Each site is relaxed bare, with H on the metal and (if it has S donors) with H on S, same cell, same fixed
carbon and same starting geometries as the DFT runs (materials/graphene/sites.py). dG_H = E(*H) - E(*) - E(H2)/2
+ 0.24 eV with H2 relaxed in the same cell. Speed over polish: spin-restricted, Fermi smearing at 1000 K,
BFGS to 0.05 eV/A (at most 400 steps). One JSON + one .xyz per structure in data/screen/xtb/, so the scan
resumes where it stopped and the relaxed geometries can start the DFT runs of the best candidates.
"""
import glob
import itertools
import json
import os
import sys
import time
from multiprocessing import Pool

os.environ.setdefault("OMP_NUM_THREADS", "1")
import numpy as np                       # noqa: E402
import pandas as pd                      # noqa: E402
from ase import Atoms                    # noqa: E402
from ase.constraints import FixAtoms     # noqa: E402
from ase.io import write                 # noqa: E402
from ase.optimize import BFGS            # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path[:0] = [os.path.join(ROOT, "materials", "graphene"), os.path.join(ROOT, "dft")]
from batches import GRAPHENE_BATCH1, GRAPHENE_BATCH2   # noqa: E402
from qe import anchor                                   # noqa: E402
from sites import build, with_h, h_sites, parse, _assign_dv, C_CELL, Z_SHEET   # noqa: E402

OUT = os.path.join(ROOT, "data", "screen", "xtb")
TABLE = os.path.join(ROOT, "data", "xtb_scan.csv")
METALS = ["V", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W"]     # the enzyme metals of the sweep
WORKERS = int(os.environ.get("SLURM_NTASKS", 13))   # 13 of the 14 local cores, else the job's cores
FMAX, STEPS, ETEMP = 0.05, 400, 1000.0
SCF_RETRY = {"max_iterations": 1000, "mixer_damping": 0.1}
DG_CORR = 0.24
MIN_ENTRIES = 5                  # a donor environment enters the scan if >= 5 representative entries have it
TRANSITION_METALS = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",
                     "Rh", "Pd", "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]
PROOF_METALS = ["Cr", "Ru", "Pt"]   # beyond the enzyme metals, a proof of concept; the other 17 are the roadmap
LIBRARY_METALS = METALS + PROOF_METALS
DONORS = ["N", "S", "O", "C"]
DV_SYMMETRY = [(0, 1, 2, 3), (2, 3, 0, 1), (1, 0, 3, 2), (3, 2, 1, 0)]   # graphene symmetries acting on d1..d4
LIBRARY_JOBS = os.path.join(ROOT, "data", "xtb_library_jobs.json")
LIBRARY_START = os.path.join(ROOT, "data", "screen", "xtb_library_start.extxyz")


def relax(job):
    """job = (motif, state): relax one structure, write <OUT>/<motif dir>/<state>.json/.xyz; skip if done."""
    from tblite.ase import TBLite
    motif, state = job
    d = os.path.join(OUT, motif.replace(":", "_").replace("=", ""))
    res_file = os.path.join(d, state + ".json")
    retry = False
    if os.path.exists(res_file):
        res = json.load(open(res_file))
        if res.get("energy") is not None or res.get("retried"):
            return res
        retry = True                         # SCF failure: once more with the damped SCF of xtb_reactions.py
    os.makedirs(d, exist_ok=True)
    if motif == "H2":
        atoms = Atoms("H2", positions=[[0, 0, 0], [0, 0, 0.75]], cell=build("Co-N4:DV")[0].cell, pbc=True)
        atoms.positions += atoms.cell.sum(axis=0) / 2
        info = None
    else:
        atoms, info = build(motif)
        atoms.set_constraint(FixAtoms([anchor(atoms, info["imetal"])]))
        if state != "bare":
            atoms = with_h(atoms, info, state[2:])
    atoms.calc = TBLite(method="GFN2-xTB", electronic_temperature=ETEMP, verbosity=0, **(SCF_RETRY if retry else {}))
    t = time.time()
    res = {"motif": motif, "state": state, "n_atoms": len(atoms), "retried": retry}
    try:
        opt = BFGS(atoms, logfile=None, maxstep=0.15)
        conv = opt.run(fmax=FMAX, steps=STEPS)
        res.update(energy=float(atoms.get_potential_energy()), steps=opt.nsteps, converged=bool(conv))
    except Exception as exc:                 # SCC failure etc.: record it, the scan goes on
        res.update(energy=None, steps=None, converged=False, error="%s: %s" % (type(exc).__name__, exc))
    res["seconds"] = round(time.time() - t, 1)
    if info is not None and res.get("energy") is not None:
        im = info["imetal"]
        dist = atoms.get_distances(im, info["idonors"], mic=True)
        zc = np.median(atoms.positions[[a.index for a in atoms if a.symbol == "C"], 2])
        res.update(max_M_X=float(dist.max()), M_height=float(atoms.positions[im, 2] - zc))
        if state != "bare":
            dh = atoms.get_distances(len(atoms) - 1, range(len(atoms) - 1), mic=True)
            j = int(np.argmin(dh))
            res.update(H_partner="M" if j == im else atoms[j].symbol + ("(donor)" if j in info["idonors"] else ""),
                       H_dist=float(dh[j]))
    write(os.path.join(d, state + ".xyz"), atoms)
    json.dump(res, open(res_file, "w"), indent=1)
    return res


def jobs_for(motifs):
    jobs = [("H2", "ref")]
    for m in motifs:
        atoms, info = build(m)
        jobs += [(m, "bare")] + [(m, "H_" + s) for s in h_sites(info, atoms)]
    return jobs


def environments(min_entries=MIN_ENTRIES):
    """In-plane donor environments of the PDB sweep (biological metals, representative copies), as builder
    motif suffixes: ('N2S', 'SV', None), ('N2S2', 'DV', 'cis5'), ... with the number of entries."""
    s = pd.read_parquet(os.path.join(ROOT, "data", "complete", "sites.parquet"),
                        columns=["entry_id", "metal", "is_representative", "preferred_site",
                                 "dv_formula", "dv_isomer", "sv_formula"])
    s = s[s.is_representative & s.metal.isin(METALS) & s.preferred_site.isin(["DV", "SV"])]
    s["formula"] = np.where(s.preferred_site == "DV", s.dv_formula, s.sv_formula)
    s["isomer"] = np.where(s.preferred_site == "DV", s.dv_isomer.fillna(""), "")
    s = s[~s.formula.str.contains("Se|Cl", regex=True)]
    env = s.groupby(["formula", "preferred_site", "isomer"]).entry_id.nunique().reset_index(name="entries")
    return env[env.entries >= min_entries].sort_values("entries", ascending=False)


def builder_isomer(formula, isomer):
    """Sweep isomer (cis/trans) -> sites.py isomer: cis5/trans for 2+2, explicit d1..d4 order for 2+1+1."""
    import re
    n = {el: int(k or 1) for el, k in re.findall(r"([A-Z][a-z]?)(\d*)", formula)}
    counts = sorted(n.values())
    if counts == [2, 2]:
        return "trans" if isomer == "trans" else "cis5"
    if counts == [1, 1, 2]:
        pair = next(el for el, k in n.items() if k == 2)
        a, b = [el for el, k in n.items() if k == 1]
        return pair + a + b + pair if isomer == "trans" else pair + a + pair + b    # trans: d1,d4; cis: d1,d3
    if counts == [1, 1, 1, 1]:
        return "".join(n)
    return None


def scan_motifs():
    env = environments()
    motifs = []
    for r in env.itertuples():
        iso = builder_isomer(r.formula, r.isomer) if r.preferred_site == "DV" else None
        for m in METALS:
            motifs.append("%s-%s:%s%s" % (m, r.formula, r.preferred_site, ":" + iso if iso else ""))
    return motifs, env


def formula(donors):
    return "".join(el + (str(donors.count(el)) if donors.count(el) > 1 else "") for el in DONORS if el in donors)


def arrangement(motif):
    """(metal, site, donor order up to symmetry) of a motif: SV as a sorted triple, DV as the smallest of d1..d4
    under the graphene symmetries; motifs with axial ligands or other donors give None."""
    m = parse(motif)
    if m["axial"] or not set(m["donors"]) <= set(DONORS):
        return None
    if m["site"] == "SV":
        return m["metal"], "SV", tuple(sorted(m["donors"]))
    order = _assign_dv(m["donors"], m["isomer"])
    return m["metal"], "DV", min(tuple(order[g[i]] for i in range(4)) for g in DV_SYMMETRY)


def library_motifs():
    """Every arrangement of N/S/O/C donors (20 in a single vacancy, 76 in a double vacancy) for every metal of
    LIBRARY_METALS (TRANSITION_METALS for the full design space), named as the scan names them (cis5/cis6/trans for 2+2, else the explicit order d1..d4)."""
    sv = sorted({tuple(sorted(c)) for c in itertools.product(DONORS, repeat=3)})
    dv = sorted({min(tuple(c[g[i]] for i in range(4)) for g in DV_SYMMETRY)
                 for c in itertools.product(DONORS, repeat=4)})
    names = ["%s:SV" % formula(list(t)) for t in sv]
    for t in dv:
        counts = sorted(t.count(el) for el in set(t))
        if counts in ([4], [1, 3]):
            names.append("%s:DV" % formula(list(t)))
        elif counts == [2, 2]:
            iso = next(i for i in ("cis5", "cis6", "trans") if arrangement("Fe-%s:DV:%s" % (formula(list(t)), i))[2] == t)
            names.append("%s:DV:%s" % (formula(list(t)), iso))
        else:
            names.append("%s:DV:%s" % (formula(list(t)), "".join(t)))
    return ["%s-%s" % (m, n) for m in LIBRARY_METALS for n in names]


def computed_arrangements():
    """Arrangements that already have a relaxed bare site in data/screen/xtb (scan, bench, library)."""
    done = set()
    for d in os.listdir(OUT) if os.path.isdir(OUT) else []:
        res = os.path.join(OUT, d, "bare.json")
        if os.path.exists(res):
            r = json.load(open(res))
            if r.get("energy") is not None:
                try:
                    done.add(arrangement(r["motif"]))
                except (ValueError, KeyError):
                    pass
    return done


def library():
    """The library sites not computed yet -> LIBRARY_JOBS (bare + H positions) and their starting structures."""
    from ase.io import write as write_xyz
    motifs = library_motifs()
    done = computed_arrangements()
    todo = [m for m in motifs if arrangement(m) not in done]
    jobs = jobs_for(todo)
    json.dump(jobs, open(LIBRARY_JOBS, "w"))
    seen = {(r.preferred_site, r.formula) for r in environments(min_entries=1).itertuples()}   # sweep: N, S, O, C order
    frames = []
    for m in todo:
        atoms, info = build(m)
        atoms.info.update(motif=m, metal=info["metal"], site=info["site"],
                          donor_set_in_pdb=(info["site"], formula(sorted(info["donors"]))) in seen)
        frames.append(atoms)
    os.makedirs(os.path.dirname(LIBRARY_START), exist_ok=True)
    write_xyz(LIBRARY_START, frames, format="extxyz")
    print("%d arrangements x %d metals = %d sites; %d already computed; %d sites -> %d structures -> %s" % (
        len(motifs) // len(LIBRARY_METALS), len(LIBRARY_METALS), len(motifs), len(motifs) - len(todo),
        len(todo), len(jobs), LIBRARY_JOBS))
    return jobs


def run(jobs):
    jobs = [j for j in dict.fromkeys(jobs)]
    t = time.time()
    with Pool(WORKERS) as pool:
        for k, r in enumerate(pool.imap_unordered(relax, jobs), 1):
            if k % 10 == 0 or k == len(jobs):
                print("  %d/%d structures, %.0f min" % (k, len(jobs), (time.time() - t) / 60), flush=True)


def table():
    rows = []
    for r, _, fs in os.walk(OUT):
        for f in fs:
            if f.endswith(".json"):
                try:
                    rows.append(json.load(open(os.path.join(r, f))))
                except json.JSONDecodeError:         # still being written by a running worker
                    pass
    df = pd.DataFrame(rows)
    e_h2 = df[(df.motif == "H2")].energy.iloc[0]
    bare = df[df.state == "bare"].set_index("motif")
    h = df[df.state.str.startswith("H_")].copy()
    h["E_bare"] = h.motif.map(bare.energy)
    h["bare_converged"] = h.motif.map(bare.converged)
    h["max_M_X_bare"] = h.motif.map(bare.max_M_X)
    h["dE_H"] = h.energy - h.E_bare - e_h2 / 2
    h["dG_H"] = h.dE_H + DG_CORR
    library = ({j[0] for j in json.load(open(LIBRARY_JOBS))} if os.path.exists(LIBRARY_JOBS) else set()) - {"H2"}
    h["source"] = np.where(h.motif.isin(library), "library", "scan")   # scan: enzyme environments and the DFT bench
    cols = ["motif", "source", "state", "dG_H", "dE_H", "H_partner", "H_dist", "max_M_X", "max_M_X_bare", "M_height",
            "converged", "bare_converged", "steps", "seconds", "error"]
    h = h.reindex(columns=cols).sort_values("motif")
    h.to_csv(TABLE, index=False, float_format="%.4f")
    print("%d H structures, %d with dG_H -> %s" % (len(h), h.dG_H.notna().sum(), TABLE))
    return h


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "bench":
        run(jobs_for([m for m, _ in GRAPHENE_BATCH1 + GRAPHENE_BATCH2]))
        table()
    elif mode == "scan":
        motifs, env = scan_motifs()
        print("%d environments x %d metals = %d sites" % (len(env), len(METALS), len(motifs)))
        print(env.to_string(index=False))
        run(jobs_for(motifs))
        table()
    elif mode == "jobs":
        motifs, env = scan_motifs()
        jobs = jobs_for(motifs)
        json.dump(jobs, open(os.path.join(ROOT, "data", "xtb_scan_jobs.json"), "w"))
        print("%d structures of %d sites -> data/xtb_scan_jobs.json" % (len(jobs), len(motifs)))
    elif mode == "shard":
        k, n = int(sys.argv[2]), int(sys.argv[3])
        jobs = [tuple(j) for j in json.load(open(os.path.join(ROOT, "data", "xtb_scan_jobs.json")))]
        run(jobs[k::n])
    elif mode == "library":
        library()
    elif mode == "library-shard":
        k, n = int(sys.argv[2]), int(sys.argv[3])
        run([tuple(j) for j in json.load(open(LIBRARY_JOBS))][k::n])
    elif mode == "retry":
        failed = []
        for f in sorted(glob.glob(os.path.join(OUT, "*", "*.json"))):
            r = json.load(open(f))
            if r.get("energy") is None and not r.get("retried"):
                failed.append((r["motif"], r["state"]))
        print("%d structures with an SCF failure, once more with the damped SCF" % len(failed))
        run(failed)
    elif mode == "table":
        table()
    else:
        sys.exit(__doc__)
