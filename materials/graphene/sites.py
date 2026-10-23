"""Single-atom sites in graphene: a metal in a single (SV) or double vacancy (DV) with N/S/O/C donors.

    python materials/graphene/sites.py                      # build every motif of dft/batches.py, print a summary
    python materials/graphene/sites.py Ni-N2S2:DV:cis5 out/  # build one motif, write bare and H* as .xyz

A motif is "<M>-<donors>:<SV|DV>[:<isomer>][:ax=<ligand>]", e.g. "Fe-N4:DV:ax=NH3" or "Cu-N2S:SV".
Donors are the dangling atoms of the vacancy; "C" means an undoped edge carbon (the graphene analogue of
the CO/CN carbons of hydrogenases). The vacancy sits at the centre of a 6x6 graphene cell.

DV geometry: the two removed atoms lie along x, the four dangling atoms d1..d4 at (-x,+y), (-x,-y), (+x,+y),
(+x,-y). d1-d3 close a five-membered ring with the metal, d1-d2 a six-membered ring, d1-d4 sit trans.
Isomers of X2Y2: cis5 (the two Y on a pentagon pair, default), cis6 (hexagon pair), trans; any other
arrangement can be given explicitly as the d1..d4 elements, e.g. "SNSN".
SV: the three dangling atoms are equivalent by the C3 axis, so a composition fixes the site; the metal
starts 1.25 A above the plane (it cannot fit in a single vacancy).
Axial ligands sit below the sheet, so the top face stays free for H: ax=NH3 (proximal His of a heme),
ax=OH, ax=H2O.
"""
import os
import re
import sys

import numpy as np
from ase import Atoms
from ase.io import write

A_GRAPHENE = 2.467      # PBE lattice constant (A)
N_CELL = 6              # 6x6 cell: 14.8 A between periodic images of the site
C_CELL = 20.0           # cell height (A); sheet at z = Z_SHEET, H and axial ligands stay >= 6 A from the dipole layer
Z_SHEET = 8.0
METALS = {"Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
          "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"}   # 3d, 4d, 5d transition metals
DV_ISOMERS = {"cis5": (0, 2), "cis6": (0, 1), "trans": (0, 3)}   # positions of the minority/second donor pair
BOND = {"N": 1.90, "C": 1.90, "O": 1.95, "S": 2.20, "Se": 2.35}   # target M-X (A) for the starting metal height
S_LIFT = {"DV": 1.0, "SV": 0.5}   # S/Se start this far above the sheet: in plane a C-S bond would be 1.42 A
H_RANGE = {"DV": (0.0, 1.0), "SV": (0.8, 2.5)}   # allowed starting heights of the metal above the sheet (A)
D_MH = 1.55             # starting M-H distance (A)
D_SH = 1.36             # starting S-H distance (A)
RATTLE = 0.02           # random displacement (A) so relaxations do not start on a symmetric saddle


def parse(motif):
    """'Ni-N2S2:DV:cis5' -> dict(metal, donors list, site, isomer, axial)."""
    head, *opts = motif.split(":")
    metal, formula = head.split("-")
    donors = []
    for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        donors += [el] * int(n or 1)
    site = opts[0] if opts else ("DV" if len(donors) == 4 else "SV")
    isomer = next((o for o in opts[1:] if not o.startswith("ax=")), None)
    axial = next((o[3:] for o in opts[1:] if o.startswith("ax=")), None)
    if metal not in METALS or site not in ("SV", "DV") or len(donors) != {"SV": 3, "DV": 4}[site]:
        raise ValueError("bad motif %r" % motif)
    return {"metal": metal, "donors": donors, "site": site, "isomer": isomer, "axial": axial}


def graphene():
    a1 = A_GRAPHENE * np.array([1.0, 0.0, 0.0])
    a2 = A_GRAPHENE * np.array([0.5, np.sqrt(3) / 2, 0.0])
    basis = [np.zeros(3), (a1 + a2) / 3]
    pos = [i * a1 + j * a2 + b for i in range(N_CELL) for j in range(N_CELL) for b in basis]
    cell = [N_CELL * a1, N_CELL * a2, [0, 0, C_CELL]]
    atoms = Atoms("C%d" % len(pos), positions=pos, cell=cell, pbc=True)
    atoms.positions[:, 2] = Z_SHEET
    return atoms


def neighbours(atoms, i, rmax=1.6):
    d = atoms.get_distances(i, range(len(atoms)), mic=True)
    return [j for j in np.argsort(d) if 0 < d[j] < rmax]


def _assign_dv(donors, isomer):
    """Elements for d1..d4."""
    if isomer and isomer not in DV_ISOMERS:          # explicit order, e.g. "SNSN"
        order = re.findall(r"[A-Z][a-z]?", isomer)
        if sorted(order) != sorted(donors):
            raise ValueError("isomer %s does not match donors %s" % (isomer, donors))
        return order
    counts = {el: donors.count(el) for el in dict.fromkeys(donors)}
    if len(counts) == 1:
        return list(donors)
    els = sorted(counts, key=lambda el: (counts[el], el))     # minority first
    if sorted(counts.values()) == [1, 3]:                     # all four positions are equivalent
        return [els[0]] + [els[1]] * 3
    if sorted(counts.values()) == [2, 2]:
        pair = DV_ISOMERS[isomer or "cis5"]
        return [els[0] if k in pair else els[1] for k in range(4)]
    raise ValueError("give the d1..d4 order explicitly for %s" % donors)


def build(motif, seed=0):
    """Bare site: (atoms, info); info = the parsed motif plus imetal and idonors (atom indices)."""
    m = parse(motif)
    atoms = graphene()
    centre = atoms.cell.sum(axis=0) / 2
    centre[2] = Z_SHEET
    i = int(np.argmin(np.linalg.norm(atoms.positions - centre, axis=1)))
    if m["site"] == "SV":
        removed = [i]
        hole = atoms.positions[i].copy()
        dangling = neighbours(atoms, i)[:3]
    else:
        j = neighbours(atoms, i)[0]
        removed = [i, j]
        hole = (atoms.positions[i] + atoms.positions[j]) / 2
        ex = atoms.positions[j] - atoms.positions[i]
        ex /= np.linalg.norm(ex)
        ey = np.cross([0, 0, 1.0], ex)
        dangling = [k for k in set(neighbours(atoms, i) + neighbours(atoms, j)) if k not in removed]
        local = {k: ((atoms.positions[k] - hole) @ ex, (atoms.positions[k] - hole) @ ey) for k in dangling}
        dangling = sorted(dangling, key=lambda k: (local[k][0] > 0, local[k][1] < 0))   # d1 (-,+) d2 (-,-) d3 (+,+) d4 (+,-)
    elements = m["donors"] if m["site"] == "SV" else _assign_dv(m["donors"], m["isomer"])
    for k, el in zip(dangling, elements):
        atoms[k].symbol = el
        if el in ("S", "Se"):        # large donors pucker towards the metal side
            atoms.positions[k, 2] += S_LIFT[m["site"]]
    keep = [k for k in range(len(atoms)) if k not in removed]
    donor_pos = atoms.positions[dangling].copy()
    atoms = atoms[keep]
    metal_pos = hole + [0, 0, metal_height(hole, donor_pos, elements, m["site"])]
    atoms += Atoms(m["metal"], positions=[metal_pos])
    if m["axial"]:
        atoms += axial_ligand(m["axial"], metal_pos)
    rng = np.random.default_rng(seed)
    atoms.positions += rng.normal(0, RATTLE, atoms.positions.shape)
    imetal = [a.index for a in atoms if a.symbol == m["metal"]][0]
    idonor = [int(np.argmin(np.linalg.norm(atoms.positions - p, axis=1))) for p in donor_pos]
    return atoms, dict(m, imetal=imetal, idonors=idonor)


def metal_height(hole, donor_pos, elements, site):
    """Height above the sheet that best matches the target M-X bond lengths (least squares on a grid)."""
    hs = np.arange(*H_RANGE[site], 0.01)
    target = np.array([BOND[el] for el in elements])
    d = np.linalg.norm(hole[None, None, :2] - donor_pos[None, :, :2], axis=2)
    dz = hs[:, None] - (donor_pos[None, :, 2] - hole[2])
    err = ((np.sqrt(d ** 2 + dz ** 2) - target) ** 2).sum(axis=1)
    return hs[int(np.argmin(err))]


def axial_ligand(name, metal_pos):
    """Ligand below the metal (the H side is +z)."""
    z = metal_pos + [0, 0, 0]
    if name == "NH3":
        n = z - [0, 0, 2.05]
        hs = [n + [0.94 * np.cos(t), 0.94 * np.sin(t), -0.33] for t in np.radians([0, 120, 240])]
        return Atoms("NH3", positions=[n] + hs)
    if name == "OH":
        o = z - [0, 0, 1.85]
        return Atoms("OH", positions=[o, o - [0, 0, 0.97]])
    if name == "H2O":
        o = z - [0, 0, 2.15]
        return Atoms("OH2", positions=[o, o + [0.76, 0, -0.59], o + [-0.76, 0, -0.59]])
    raise ValueError("unknown axial ligand %s" % name)


def with_h(atoms, info, where="M"):
    """Copy with one H on the metal (where='M') or on the first S donor (where='S')."""
    atoms = atoms.copy()
    pm = atoms.positions[info["imetal"]]
    if where == "M":
        h = pm + [0, 0, D_MH]
    else:
        k = next(k for k in info["idonors"] if atoms[k].symbol == where)
        ps = atoms.positions[k]
        out = (ps - pm) * [1, 1, 0]
        out = out / np.linalg.norm(out)
        h = ps + D_SH * (0.55 * out + 0.835 * np.array([0, 0, 1.0]))   # ~57 deg from the normal, away from M
    atoms += Atoms("H", positions=[h])
    return atoms


def h_sites(info, atoms):
    """Adsorption sites worth computing: always the metal, plus S when the site has S donors."""
    return ["M"] + (["S"] if any(atoms[k].symbol == "S" for k in info["idonors"]) else [])


if __name__ == "__main__":
    if len(sys.argv) > 1:
        motif = sys.argv[1]
        out = sys.argv[2] if len(sys.argv) > 2 else "."
        os.makedirs(out, exist_ok=True)
        at, inf = build(motif)
        write(os.path.join(out, "bare.xyz"), at)
        for w in h_sites(inf, at):
            write(os.path.join(out, "H_%s.xyz" % w), with_h(at, inf, w))
        print(motif, at.get_chemical_formula(), "metal", inf["imetal"], "donors", inf["idonors"])
    else:
        sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", "dft"))
        from batches import GRAPHENE_BATCH1
        for motif, why in GRAPHENE_BATCH1:
            at, inf = build(motif)
            d = at.get_distances(inf["imetal"], inf["idonors"], mic=True)
            print("%-22s %-10s M-X %s  H on %s  # %s" % (motif, at.get_chemical_formula(),
                  " ".join("%.2f" % x for x in d), "+".join(h_sites(inf, at)), why))
