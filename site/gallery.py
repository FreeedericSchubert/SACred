"""Structures for the pitch gallery: each catalyst motif as relaxed (DFT where there is one, else GFN2-xTB), with H
on the metal, cut to a round patch around the metal, support normal along +z.

    python site/gallery.py      # -> notebooks/assets/gallery.json (XYZ per motif, drawn with 3Dmol.js in the pitch)
"""
import json
import os

import numpy as np
from ase.io import read

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "notebooks", "assets")
RUNS = os.path.join(ROOT, "dft", "runs")
SITES = {   # key: (structure file, format, metal, patch radius in A, view)
    "fe_s2c": (os.path.join(RUNS, "graphene", "Fe-S2C_SV", "H_M", "relax.out"), "espresso-out", "Fe", 6.0, "-62x"),
    "fe_ns2": (os.path.join(RUNS, "graphene", "Fe-NS2_SV", "H_M", "relax.out"), "espresso-out", "Fe", 6.0, "-62x"),
    "cu_n3": (os.path.join(RUNS, "graphene", "Cu-N3_SV", "H_M", "relax.out"), "espresso-out", "Cu", 6.0, "-62x"),
    "co_n2oc": (os.path.join(RUNS, "graphene", "Co-N2OC_DV_NONC", "H_M", "relax.out"), "espresso-out", "Co", 6.5, "-62x"),
    "fe_mos2": (os.path.join(RUNS, "mos2", "34_mos2_fe_h", "relax.out"), "espresso-out", "Fe", 5.2, "-62x"),
    "co_silica": (os.path.join(ROOT, "data", "screen", "xtb_supports", "SiO2", "Co_T7_H_M.xyz"), "extxyz", "Co", 9.0,
                  "-70x"),
}


def patch(atoms, metal, radius):
    """Atoms within `radius` of the (last) metal atom, unwrapped around it, metal at the origin."""
    im = [a.index for a in atoms if a.symbol == metal][-1]
    if atoms.pbc.any():
        d = atoms.get_distances(im, range(len(atoms)), mic=True, vector=True)
    else:
        d = atoms.positions - atoms.positions[im]
    keep = np.linalg.norm(d[:, :2] if atoms.pbc.any() else d, axis=1) <= radius
    out = atoms[keep]
    out.positions = d[keep]
    out.set_cell(np.eye(3) * 0.0)
    out.pbc = False
    return out


def main():
    os.makedirs(OUT, exist_ok=True)
    out = {}
    for key, (path, fmt, metal, radius, _) in SITES.items():
        p = patch(read(path, format=fmt, index=-1), metal, radius)
        if key == "co_silica":                       # turn the cage so the metal sits on top
            cage = p.positions[[a.index for a in p if a.symbol == "Si"]].mean(axis=0)
            p.rotate(-cage, (0, 0, 1), center=(0, 0, 0))
        out[key] = "%d\n\n" % len(p) + "\n".join("%s %.3f %.3f %.3f" % (a.symbol, *a.position) for a in p) + "\n"
        print(key, len(p), "atoms")
    json.dump(out, open(os.path.join(OUT, "gallery.json"), "w"))


if __name__ == "__main__":
    main()
