"""Does the metal leave its site at pyrolysis temperature? Langevin MD with any ASE calculator (the trained force
field on a GPU), counting how long each metal stays.

    python screen/md_escape.py --calc mypkg.potential:make_calc --kwargs '{"model": "ff.model", "device": "cuda"}' \
        --sites Zn-NS2:SV Ni-NS2:SV Co-N2OC:DV:NONC Fe-NO3:DV --T 1000 --ps 200 --replicas 4
    python screen/md_escape.py --calc ... --sites Ni-NS2:SV --T 1200 1400 1600 1800 --ps 200 --replicas 4  # Arrhenius
    python screen/md_escape.py --calc xtb --sites Zn-NS2:SV --T 1000 --ps 0.05      # CPU smoke test (GFN2-xTB)
    python screen/md_escape.py --summary                                           # tables only

--calc module:function returns an ASE calculator (called with --kwargs); "xtb" uses tblite on the CPU.
Start: the relaxed bare site in the 6x6 cell of the training data (DFT relax.out where finished, else GFN2-xTB), one
far carbon fixed as in the DFT runs so the sheet does not drift. Langevin thermostat (friction 0.01 /fs), velocities
from Maxwell-Boltzmann at T, time step 1 fs (no H in the bare sites).
Escape: fewer than two donors within their bond cutoff (bio/first_shell.py) + 0.5 A for 1 ps in a row; the time it
started is the first-passage time. A run stops 2 ps after the escape or at --ps.
Out: data/md/<site>_<T>K_r<k>.json (time series every 10 fs, escape time) and .extxyz (a frame every 50 fs);
data/md/summary.csv. Escape rate per temperature k = n_escaped / total time (maximum likelihood for first passages,
runs without escape count with their full length); with escapes at two or more temperatures an Arrhenius fit
ln k = ln A - Ea / kT gives the apparent barrier, set against the xTB NEB barrier (data/screen/xtb_neb). At one
temperature, k = 1e13 /s exp(-E/kT) gives the expected mean time: 8 ps for 0.38 eV, 7 ns for 0.96 eV at 1000 K.
"""
import argparse
import glob
import importlib
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path[:0] = [os.path.join(ROOT, "bio"), os.path.join(ROOT, "dft"), HERE]
import first_shell as fs                 # noqa: E402
import xtb_scan as xs                    # noqa: E402

OUT = os.path.join(ROOT, "data", "md")
RUNS = os.path.join(ROOT, "dft", "runs", "graphene")
NEB = os.path.join(ROOT, "data", "screen", "xtb_neb")
KB = 8.617333e-5                         # eV/K
LOOSE = 0.5                              # A beyond the bond cutoff: the donor has let go
HOLD = 1.0                               # ps below two donors before it counts as an escape


def calculator(spec, kwargs):
    if spec == "xtb":
        from tblite.ase import TBLite
        return TBLite(method="GFN2-xTB", electronic_temperature=1000.0, verbosity=0)
    module, func = spec.split(":")
    return getattr(importlib.import_module(module), func)(**kwargs)


def start(motif):
    """Relaxed bare site (DFT if finished, else xTB) with metal and donor indices and the fixed carbon."""
    from ase.constraints import FixAtoms
    from ase.io import read
    from qe import anchor
    d = xs.dirname(motif) if hasattr(xs, "dirname") else motif.replace(":", "_").replace("=", "")
    _, info = xs.build(motif)
    out = os.path.join(RUNS, d, "bare", "relax.out")
    if os.path.exists(out) and "End final coordinates" in open(out).read():
        atoms, source = read(out, format="espresso-out", index=-1), "DFT"
    else:
        atoms, source = read(os.path.join(xs.OUT, d, "bare.xyz")), "xTB"
    atoms.set_constraint(FixAtoms([anchor(atoms, info["imetal"])]))
    return atoms, info["imetal"], info["idonors"], source, d


def neb_barrier(d):
    path = os.path.join(NEB, d + ".json")
    return json.load(open(path)).get("barrier") if os.path.exists(path) else None


def run(motif, temp, ps, replica, calc, dt=1.0, seed0=0):
    from ase import units
    from ase.io import write
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
    atoms, im, donors, source, d = start(motif)
    name = "%s_%dK_r%d" % (d, temp, replica)
    res_file = os.path.join(OUT, name + ".json")
    if os.path.exists(res_file):
        return json.load(open(res_file))
    atoms.calc = calc
    rng = np.random.default_rng(seed0 + 1000 * replica + temp)
    MaxwellBoltzmannDistribution(atoms, temperature_K=temp, rng=rng)
    Stationary(atoms)
    dyn = Langevin(atoms, timestep=dt * units.fs, temperature_K=temp, friction=0.01 / units.fs, rng=rng)
    cut = np.array([fs.CUTOFF.get(atoms[i].symbol, 2.6) + LOOSE for i in donors])
    every, frame = max(1, int(round(10 / dt))), max(1, int(round(50 / dt)))
    series, frames, t_low, t_escape = [], [], None, None
    zc = lambda a: np.median(a.positions[[x.index for x in a if x.symbol == "C"], 2])     # noqa: E731
    t0, n_steps = time.time(), int(round(ps * 1000 / dt))
    for step in range(0, n_steps + 1, every):
        if step:
            dyn.run(every)
        t = step * dt / 1000.0                                            # ps
        dist = atoms.get_distances(im, donors, mic=True)
        bound = int((dist < cut).sum())
        series.append([round(t, 3), bound, round(float(dist.max()), 3),
                       round(float(atoms.positions[im, 2] - zc(atoms)), 3), round(float(atoms.get_temperature()), 1)])
        if step % frame == 0:
            frames.append(atoms.copy())
        if t_escape is None:
            t_low = (t_low if t_low is not None else t) if bound < 2 else None
            if t_low is not None and t - t_low >= HOLD:
                t_escape = t_low
        elif t - t_escape >= HOLD + 2.0:
            break
    os.makedirs(OUT, exist_ok=True)
    write(os.path.join(OUT, name + ".extxyz"), frames, format="extxyz")
    res = {"site": motif, "T": temp, "replica": replica, "start": source, "dt_fs": dt, "ps": series[-1][0],
           "escaped": t_escape is not None, "t_escape_ps": t_escape, "neb_barrier": neb_barrier(d),
           "seconds": round(time.time() - t0, 1), "columns": ["t_ps", "donors_bound", "max_M_X", "M_height", "T_K"],
           "series": series}
    json.dump(res, open(res_file, "w"))
    return res


def summary():
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(OUT, "*.json")))]
    if not rows:
        print("no MD runs in %s" % OUT)
        return
    df = pd.DataFrame([{k: r[k] for k in ("site", "T", "replica", "ps", "escaped", "t_escape_ps", "neb_barrier")}
                       for r in rows])
    df["time_ps"] = df.t_escape_ps.where(df.escaped, df.ps)
    per = df.groupby(["site", "T"]).agg(runs=("replica", "size"), escaped=("escaped", "sum"),
                                        total_ps=("time_ps", "sum"), neb_barrier=("neb_barrier", "first")).reset_index()
    per["rate_per_ps"] = per.escaped / per.total_ps                      # first passages: n / total time
    per["expected_ps_from_neb"] = [1 / (1e13 * np.exp(-b / (KB * t))) * 1e12 if pd.notna(b) else np.nan
                                   for b, t in zip(per.neb_barrier, per["T"])]
    fits = []
    for site, g in per[per.escaped > 0].groupby("site"):
        if g["T"].nunique() >= 2:
            slope, ln_a = np.polyfit(1 / (KB * g["T"]), np.log(g.rate_per_ps * 1e12), 1)
            fits.append({"site": site, "Ea_md_eV": -slope, "prefactor_per_s": float(np.exp(ln_a)),
                         "neb_barrier_eV": g.neb_barrier.iloc[0]})
    per.to_csv(os.path.join(OUT, "summary.csv"), index=False, float_format="%.4g")
    print(per.to_string(index=False, float_format=lambda x: "%.3g" % x))
    if fits:
        pd.DataFrame(fits).to_csv(os.path.join(OUT, "arrhenius.csv"), index=False, float_format="%.4g")
        print(pd.DataFrame(fits).to_string(index=False, float_format=lambda x: "%.3g" % x))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--calc", default="xtb")
    p.add_argument("--kwargs", default="{}", help="JSON keyword arguments for the calculator factory")
    p.add_argument("--sites", nargs="+", default=["Zn-NS2:SV", "Ni-NS2:SV", "Co-N2OC:DV:NONC", "Fe-NO3:DV"])
    p.add_argument("--T", nargs="+", type=int, default=[1000])
    p.add_argument("--ps", type=float, default=100.0)
    p.add_argument("--replicas", type=int, default=1)
    p.add_argument("--dt", type=float, default=1.0)
    p.add_argument("--summary", action="store_true")
    a = p.parse_args()
    if not a.summary:
        calc = calculator(a.calc, json.loads(a.kwargs))
        for motif in a.sites:
            for temp in a.T:
                for k in range(a.replicas):
                    r = run(motif, temp, a.ps, k, calc, dt=a.dt)
                    print("%s %d K r%d: %s after %.2f ps (%.0f s)" % (
                        motif, temp, k, "escaped at %.2f ps" % r["t_escape_ps"] if r["escaped"] else "stayed",
                        r["ps"], r["seconds"]), flush=True)
    summary()


if __name__ == "__main__":
    main()
