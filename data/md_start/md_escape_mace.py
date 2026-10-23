"""Escape MD of single metal atoms in graphene at pyrolysis temperature, standalone (ASE + MACE), e.g. on Modal.

    python md_escape_mace.py --model sac_mace.model --device cuda --T 1000 --ps 100 --replicas 3
    python md_escape_mace.py --model sac_mace.model --sites Zn-NS2:SV --T 1000 --ps 20       # one site
    python md_escape_mace.py --check --model sac_mace.model    # force field against the DFT escape-path points

Input: md_start_sites.extxyz next to this script (DFT-relaxed sites; info: motif, metal_index, donor_indices,
fixed_index). Langevin dynamics (friction 0.01 /fs, 1 fs steps, velocities from Maxwell-Boltzmann), the far carbon
of the DFT runs fixed so the sheet does not drift. Escape: fewer than two donors within their bond cutoff + 0.5 A
for 1 ps in a row; the time it started is the escape time; a run ends 2 ps later or at --ps.
Output per run: <site>_<T>K_r<k>.extxyz (a frame every 50 fs) and .json (time series every 10 fs: donors still
bound, largest metal-donor distance, metal height above the sheet, temperature; escape time); summary.csv.
--check: MACE energies and forces on escape_paths/dft_single_points_on_xtb_paths.extxyz against DFT (relative
energy along each path, force RMSE), the configurations the metal passes on its way out.
"""
import argparse
import json
import os
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CUTOFF = {"N": 2.6, "O": 2.6, "S": 2.8, "C": 2.3}      # metal-donor bond cutoffs (A) of the PDB sweep
LOOSE, HOLD = 0.5, 1.0


def mace(model, device):
    from mace.calculators import MACECalculator
    return MACECalculator(model_paths=model, device=device, default_dtype="float32")


def run(atoms0, temp, ps, replica, calc, out, dt=1.0):
    from ase import units
    from ase.constraints import FixAtoms
    from ase.io import write
    from ase.md.langevin import Langevin
    from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
    atoms = atoms0.copy()
    info = atoms0.info
    im = int(info["metal_index"])
    donors = [int(x) for x in str(info["donor_indices"]).split()]
    name = "%s_%dK_r%d" % (info["motif"].replace(":", "_"), temp, replica)
    if os.path.exists(os.path.join(out, name + ".json")):
        return json.load(open(os.path.join(out, name + ".json")))
    atoms.set_constraint(FixAtoms([int(info["fixed_index"])]))
    atoms.calc = calc
    rng = np.random.default_rng(1000 * replica + temp)
    MaxwellBoltzmannDistribution(atoms, temperature_K=temp, rng=rng)
    Stationary(atoms)
    dyn = Langevin(atoms, timestep=dt * units.fs, temperature_K=temp, friction=0.01 / units.fs, rng=rng)
    cut = np.array([CUTOFF.get(atoms[i].symbol, 2.6) + LOOSE for i in donors])
    carbons = [a.index for a in atoms if a.symbol == "C"]
    every, frame = int(round(10 / dt)), int(round(50 / dt))
    series, frames, t_low, t_escape, t0 = [], [], None, None, time.time()
    for step in range(0, int(round(ps * 1000 / dt)) + 1, every):
        if step:
            dyn.run(every)
        t = step * dt / 1000.0
        dist = atoms.get_distances(im, donors, mic=True)
        bound = int((dist < cut).sum())
        h = float(atoms.positions[im, 2] - np.median(atoms.positions[carbons, 2]))
        series.append([round(t, 3), bound, round(float(dist.max()), 3), round(h, 3), round(float(atoms.get_temperature()), 1)])
        if step % frame == 0:
            frames.append(atoms.copy())
        if t_escape is None:
            t_low = (t_low if t_low is not None else t) if bound < 2 else None
            if t_low is not None and t - t_low >= HOLD:
                t_escape = t_low
        elif t - t_escape >= HOLD + 2.0:
            break
    write(os.path.join(out, name + ".extxyz"), frames, format="extxyz")
    res = {"site": info["motif"], "T": temp, "replica": replica, "ps": series[-1][0], "escaped": t_escape is not None,
           "t_escape_ps": t_escape, "seconds": round(time.time() - t0, 1),
           "columns": ["t_ps", "donors_bound", "max_M_X", "M_height", "T_K"], "series": series}
    json.dump(res, open(os.path.join(out, name + ".json"), "w"))
    return res


def check(calc):
    from ase.io import read
    path = os.path.join(HERE, "escape_paths", "dft_single_points_on_xtb_paths.extxyz")
    frames = read(path, index=":")
    by_site = {}
    for a in frames:
        e_dft, f_dft = a.get_potential_energy(), a.get_forces()
        b = a.copy()
        b.calc = calc
        by_site.setdefault(a.info["motif"], []).append((a.info["image"], e_dft, b.get_potential_energy(),
                                                        float(np.sqrt(((b.get_forces() - f_dft) ** 2).mean()))))
    for site, rows in by_site.items():
        rows.sort()
        e0d, e0m = rows[0][1], rows[0][2]
        print("%s: image, E_DFT - E_DFT(0), E_MACE - E_MACE(0) (eV), force RMSE (eV/A)" % site)
        for k, ed, em, fr in rows:
            print("  %d  %+.3f  %+.3f  %.3f" % (k, ed - e0d, em - e0m, fr))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--model", required=True)
    p.add_argument("--device", default="cuda")
    p.add_argument("--sites", nargs="+", default=None, help="motifs, default: every frame of md_start_sites.extxyz")
    p.add_argument("--T", nargs="+", type=int, default=[1000])
    p.add_argument("--ps", type=float, default=100.0)
    p.add_argument("--replicas", type=int, default=1)
    p.add_argument("--out", default=os.path.join(HERE, "md_runs"))
    p.add_argument("--check", action="store_true")
    a = p.parse_args()
    calc = mace(a.model, a.device)
    if a.check:
        check(calc)
        return
    from ase.io import read
    os.makedirs(a.out, exist_ok=True)
    rows = []
    for atoms in read(os.path.join(HERE, "md_start_sites.extxyz"), index=":"):
        if a.sites and atoms.info["motif"] not in a.sites:
            continue
        for temp in a.T:
            for k in range(a.replicas):
                r = run(atoms, temp, a.ps, k, calc, a.out)
                print("%s %d K r%d: %s (%.1f ps, %.0f s)" % (r["site"], temp, k, "escaped at %.2f ps" % r["t_escape_ps"]
                                                            if r["escaped"] else "stayed", r["ps"], r["seconds"]), flush=True)
                rows.append({k2: r[k2] for k2 in ("site", "T", "replica", "ps", "escaped", "t_escape_ps", "seconds")})
    import csv
    with open(os.path.join(a.out, "summary.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
