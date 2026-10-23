"""DFT single points on external structures (dft/make_inputs.py points), e.g. force-field optimisations with a second
coordination sphere: energy, the largest force left on the structure, magnetization, and dG_H where a label has
both a bare and an H_M structure of the same site.

    bash dft/sync.sh pull && python dft/points.py      # -> data/dft_points.csv, data/dft_points.md

dG_H = E(H_M) - E(bare) - E(H2)/2 + 0.24 eV with both single points of the same label and k points; H2 from the
3x3x1 SCF (k points 331) or the Gamma relaxation (gamma) of _ref/H2. The largest force says how far the force-field
geometry is from a DFT minimum: below about 0.05 eV/A the single point is as good as a relaxation; above about
0.2 eV/A, relax in DFT before quoting energies.
"""
import glob
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analyse import DG_CORR, RUNS, RY, relax_energy, scf_energy   # noqa: E402

DATA = os.path.join(HERE, "..", "data")
BOHR = 0.529177


def single(run):
    """(energy eV, max force eV/A, magnetization) of a converged single point, else NaNs."""
    out = os.path.join(run, "relax.out")
    if not os.path.exists(out):
        return np.nan, np.nan, np.nan
    text = open(out).read()
    if "convergence has been achieved" not in text or "JOB DONE" not in text:
        return np.nan, np.nan, np.nan
    e = float(re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)", text)[-1]) * RY
    forces = re.findall(r"atom\s+\d+\s+type\s+\d+\s+force\s+=\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)", text)
    fmax = float(np.linalg.norm(np.array(forces, float), axis=1).max()) * RY / BOHR if forces else np.nan
    mag = re.findall(r"total magnetization\s+=\s+(-?\d+\.\d+)", text)
    return e, fmax, float(mag[-1]) if mag else 0.0


def main():
    rows = []
    for meta_file in sorted(glob.glob(os.path.join(RUNS, "*", "*", "meta.json"))):
        meta = json.load(open(meta_file))
        if "label" not in meta:
            continue
        e, fmax, mag = single(os.path.dirname(meta_file))
        rows.append({"motif": meta["motif"], "label": meta["label"], "kind": meta["kind"], "kpts": meta["kpts"],
                     "source": meta["source"], "frame": meta["frame"], "E": e, "fmax_eV_A": fmax, "mag": mag})
    if not rows:
        print("no external single points under %s" % RUNS)
        return
    df = pd.DataFrame(rows)
    h2 = {"331": scf_energy(os.path.join(RUNS, "_ref", "H2", "scf.out"))[0],
          "gamma": relax_energy(os.path.join(RUNS, "_ref", "H2", "relax.out"))[0]}
    dg = []
    for (motif, label, kpts), g in df.groupby(["motif", "label", "kpts"]):
        e = g.set_index("kind").E
        if "bare" in e and "H_M" in e:
            dg.append({"motif": motif, "label": label, "kpts": kpts,
                       "dG_H": e["H_M"] - e["bare"] - h2[kpts] / 2 + DG_CORR})
    df.to_csv(os.path.join(DATA, "dft_points.csv"), index=False, float_format="%.4f")
    dref = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv"))
    lines = ["# DFT single points on external structures", "",
             "| site | label | structure | k points | max force (eV/A) | magnetization |", "|---|---|---|---|---|---|"]
    lines += ["| %s | %s | %s | %s | %s | %s |" % (r.motif, r.label, r.kind, r.kpts,
                                                 "%.3f" % r.fmax_eV_A if np.isfinite(r.fmax_eV_A) else "running",
                                                 "%.2f" % r.mag if np.isfinite(r.mag) else "") for r in df.itertuples()]
    if dg:
        lines += ["", "| site | label | dG_H (eV) | dG_H without the second sphere (DFT, relaxed) |", "|---|---|---|---|"]
        for r in dg:
            base = dref[(dref.motif == r["motif"]) & (dref.state == "H_M")].dG_H
            lines.append("| %s | %s | %s | %s |" % (r["motif"], r["label"],
                                                    "%+.2f" % r["dG_H"] if np.isfinite(r["dG_H"]) else "running",
                                                    "%+.2f" % base.iloc[0] if len(base) else "-"))
        pd.DataFrame(dg).to_csv(os.path.join(DATA, "dft_points_dgh.csv"), index=False, float_format="%.4f")
    open(os.path.join(DATA, "dft_points.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
