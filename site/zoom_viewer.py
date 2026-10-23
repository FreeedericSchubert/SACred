"""Enzyme -> active site -> first shell -> the same metal on graphene, MoS2 and silica, as one interactive 3D page.

    python site/zoom_viewer.py      # -> site/zoom_viewer.html (self-contained: 3Dmol.js from cdnjs, structures inline)

One slider zooms from the whole protein (cartoon) into the active-site pocket (residues within POCKET A as sticks)
down to the first shell (metal and donors, labelled); its last stop shows the metal as a single-atom site with H on
it: on graphene (relaxed DFT structure where there is one, else xTB), on MoS2 and on the silica model (xTB), each with
its dG_H. The protein is turned so the metal faces the viewer; the supports so their surface faces up.
"""
import json
import os

import gemmi
import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
CUTOFF = {"N": 2.6, "O": 2.6, "S": 2.8}      # metal-donor (A), as bio/first_shell.py
POCKET = 7.0                                  # A around the metal, whole residues
PATCH = 7.2                                   # A: radius of the graphene and MoS2 patches around the metal
TILT = 28.0                                   # deg: supports seen a little from above

EXAMPLES = [   # one example per support; "support" is where the zoom ends
    dict(key="graphene", support="graphene", name="Carbonic anhydrase II, Co form", pdb="6LV1",
         note="CO₂ hydration · PDB 6LV1, 1.20 Å", site="6LV1_A_CO_302_CO", assembly=False, chains="one chain",
         shell="Co bound by three His, a water and bicarbonate; two His N and two O lie nearly in one plane.",
         graphene=("xtb", "Co-N2O2:DV:cis5", "H_M",
                   "Co–N₂O₂ in a double vacancy: the in-plane His N and O become N and O dopants; the third His "
                   "and the bicarbonate O would sit axially and are left out, so H binds there (xTB)")),
    dict(key="mos2", support="MoS2", name="[NiFe] hydrogenase", pdb="4U9H", note="H₂ evolution · PDB 4U9H, 0.89 Å",
         site="4U9H_L_EJ2_601_NI", assembly=False, chains="large subunit and small subunit",
         shell="Ni bound by four Cys; two of them also bind the Fe next to it.",
         graphene=("dft", "Ni-S3:SV", "H_M", "Ni–S₃ in a single vacancy")),
    dict(key="silica", support="SiO2", name="γ-Carbonic anhydrase, Co form", pdb="1QRE",
         note="CO₂ hydration · PDB 1QRE, 1.46 Å", site="1QRE_A_CO_214_CO", assembly=True,
         chains="homotrimer; the site sits between two subunits",
         shell="Co bound by three His (one from the next subunit), bicarbonate and two waters: the three O are the "
               "side silica can take over.",
         graphene=("xtb", "Co-N2O2:DV:cis5", "H_M", "Co–N₂O₂ in a double vacancy (xTB)")),
]


# ---------------------------------------------------------------- protein

def rotation_to(v, target):
    """Rotation matrix taking unit vector v onto target."""
    v, target = v / np.linalg.norm(v), target / np.linalg.norm(target)
    c = float(v @ target)
    axis = np.cross(v, target)
    s = np.linalg.norm(axis)
    if s < 1e-8:
        return np.eye(3) if c > 0 else -np.eye(3)
    k = axis / s
    K = np.array([[0, -k[2], k[1]], [k[2], 0, -k[0]], [-k[1], k[0], 0]])
    return np.eye(3) + s * K + (1 - c) * K @ K


def protein(ex):
    """PDB text of the protein (metal facing +z), the metal's serial, its donors and the metal residue name."""
    st = gemmi.read_structure(os.path.join(DATA, "pdb", ex["pdb"] + ".cif.gz"))
    st.setup_entities()
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    if ex["assembly"]:
        model = gemmi.make_assembly(st.assemblies[0], st[0], gemmi.HowToNameCopiedChain.Short)
        new = gemmi.Structure()
        new.add_model(model)
        new.setup_entities()
        st = new
    pdb_id, chain, resname, seqid, atom = ex["site"].split("_")
    model = st[0]
    metal = next(a for ch in model if ch.name == chain for r in ch if r.name == resname and str(r.seqid.num) == seqid
                 for a in r if a.name == atom)
    mpos = np.array(metal.pos.tolist())
    # waters: keep only those bound to the metal
    for ch in model:
        for i in reversed(range(len(ch))):
            r = ch[i]
            if r.is_water() and not any(np.linalg.norm(np.array(a.pos.tolist()) - mpos) < CUTOFF["O"] for a in r):
                del ch[i]
    # turn the protein so the metal faces the viewer (+z)
    ca = np.array([a.pos.tolist() for ch in model for r in ch for a in r if a.name == "CA"])
    centre = ca.mean(axis=0)
    R = rotation_to(mpos - centre, np.array([0, 0, 1.0]))
    for ch in model:
        for r in ch:
            for a in r:
                p = R @ (np.array(a.pos.tolist()) - centre)
                a.pos = gemmi.Position(*p)
    mpos = R @ (mpos - centre)
    # atoms are picked in the page by chain, residue number and atom name (PDB serials shift at every TER)
    spec = lambda ch, r, a: {"chain": ch.name, "resi": r.seqid.num, "atom": a.name}   # noqa: E731
    metal_spec, donors = None, []
    for ch in model:
        for r in ch:
            for a in r:
                d = np.linalg.norm(np.array(a.pos.tolist()) - mpos)
                if r.name == resname and a.name == atom and ch.name == chain and d < 0.01:
                    metal_spec = spec(ch, r, a)
                elif a.element.name in CUTOFF and d <= CUTOFF[a.element.name]:
                    label = "water" if r.is_water() else ("bicarbonate" if r.name == "BCT" else
                                                          "%s%d" % (r.name.capitalize(), r.seqid.num))
                    donors.append(dict(spec(ch, r, a), label="%s %s" % (a.element.name, label), dist=round(d, 2),
                                       pos=[round(float(x), 3) for x in a.pos.tolist()]))
    return st.make_pdb_string(), metal_spec, donors, resname, mpos


# ---------------------------------------------------------------- supports

def oriented(atoms, im, normal):
    """Metal at the origin, support normal turned up (+y) and tilted towards the viewer."""
    atoms = atoms.copy()
    atoms.positions -= atoms.positions[im]
    atoms.rotate(normal, (0, 1, 0), center=(0, 0, 0))
    atoms.rotate(TILT, "x", center=(0, 0, 0))
    return atoms


def patch_2d(atoms, im, tiles):
    """Round patch of a periodic sheet around the metal (tiles: images per direction; the metal's own images left
    out) -> (patch, index of the metal in it)."""
    a1, a2 = atoms.cell[0], atoms.cell[1]
    c = atoms.positions[im, :2]
    sym, pos, new_im = [], [], None
    for i in range(-tiles, tiles + 1):
        for j in range(-tiles, tiles + 1):
            for at in atoms:
                p = at.position + i * a1 + j * a2
                if np.linalg.norm(p[:2] - c) <= PATCH and (at.index != im or (i, j) == (0, 0)):
                    if at.index == im:
                        new_im = len(sym)
                    sym.append(at.symbol)
                    pos.append(p)
    patch = Atoms(sym, positions=pos)
    patch.info["im"] = new_im
    if "C" in sym:                       # graphene: drop edge carbons left with one neighbour
        for _ in range(3):
            d = patch.get_all_distances()
            nn = (d < 1.6).sum(axis=1) - 1
            keep = [k for k in range(len(patch)) if nn[k] >= 2 or patch[k].symbol != "C"]
            new_im = keep.index(new_im)
            patch = patch[keep]
    return patch, new_im


def to_xyz(atoms):
    return "%d\n\n" % len(atoms) + "\n".join("%s %.4f %.4f %.4f" % (a.symbol, *a.position) for a in atoms) + "\n"


def metal_index(atoms, metal):
    return [a.index for a in atoms if a.symbol == metal][-1]


def kabsch(P, Q):
    """Proper rotation R minimising |R P_i - Q_i| (rows are vectors, both centred on the metal)."""
    H = P.T @ Q
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1, 1, d]) @ U.T


def quaternion(R):
    """(w, x, y, z) of a rotation matrix."""
    w = np.sqrt(max(0.0, 1 + np.trace(R))) / 2
    if w > 1e-6:
        return [w, (R[2, 1] - R[1, 2]) / (4 * w), (R[0, 2] - R[2, 0]) / (4 * w), (R[1, 0] - R[0, 1]) / (4 * w)]
    i = int(np.argmax(np.diag(R)))
    j, k = (i + 1) % 3, (i + 2) % 3
    x = np.sqrt(max(0.0, 1 + R[i, i] - R[j, j] - R[k, k])) / 2
    q = [0.0, 0.0, 0.0, 0.0]
    q[1 + i] = x
    q[1 + j] = (R[j, i] + R[i, j]) / (4 * x)
    q[1 + k] = (R[k, i] + R[i, k]) / (4 * x)
    q[0] = (R[k, j] - R[j, k]) / (4 * x)
    return q


def align(sup, metal, mpos, donors):
    """Rigid placement of a support model (metal at the origin) onto the enzyme's first shell: the rotation that puts
    its donors onto the enzyme donors (best assignment over all pairings, elements ignored) and the metal onto the
    metal -> {"q": quaternion, "t": metal position in the protein frame, "rmsd"}. The page morphs from this
    placement to the support's own view."""
    import itertools
    im = metal_index(sup, metal)
    sd = [a.index for a in sup if a.symbol not in ("H", metal) and np.linalg.norm(a.position - sup.positions[im]) < 2.9]
    S = sup.positions[sd] - sup.positions[im]
    E = np.array([d["pos"] for d in donors]) - mpos
    n = min(len(S), len(E), 4)
    S = S[np.argsort(np.linalg.norm(S, axis=1))[:n]]          # the n closest support donors, in a fixed order
    P = S / np.linalg.norm(S, axis=1)[:, None]
    best = None
    for ei in itertools.permutations(range(len(E)), n):      # every assignment of enzyme donors to them
        Q = E[list(ei)] / np.linalg.norm(E[list(ei)], axis=1)[:, None]
        R = kabsch(P, Q)
        rmsd = float(np.sqrt(((P @ R.T - Q) ** 2).sum(axis=1).mean()))
        if best is None or rmsd < best[0]:
            best = (rmsd, R)
    rmsd, R = best
    return {"q": [round(float(x), 6) for x in quaternion(R)], "t": [round(float(x), 3) for x in mpos],
            "rmsd": round(rmsd, 3)}


def supports(ex, metal, mpos, donors):
    """[{support, label, xyz, text, align}] for graphene, MoS2 and silica."""
    out = []
    kind, motif, state, text = ex["graphene"]
    d = motif.replace(":", "_").replace("=", "")
    if kind == "dft":
        atoms = read(os.path.join(ROOT, "dft", "runs", "graphene", d, state, "relax.out"), format="espresso-out", index=-1)
        dgh = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv"))
        val = dgh[(dgh.motif == motif) & (dgh.state == state)].dG_H.dropna()
        value = "ΔG_H = %+.2f eV (DFT)" % val.iloc[0] if len(val) else ""
    else:
        atoms = read(os.path.join(DATA, "screen", "xtb", d, state + ".xyz"))
        x = pd.read_csv(os.path.join(DATA, "xtb_scan.csv"))
        o = pd.read_csv(os.path.join(DATA, "xtb_offsets.csv"))
        off = o[o.used & (o.state == "H_M")].set_index("metal").offset.get(metal)
        val = x[(x.motif == motif) & (x.state == state)].dG_H
        value = ("ΔG_H = %+.2f eV (xTB)" % val.iloc[0] + ("" if off is None else ", %+.2f eV corrected" % (val.iloc[0] + off))
                 if len(val) else "")
    p, pm = patch_2d(atoms, metal_index(atoms, metal), 0)
    g = oriented(p, pm, np.array([0, 0, 1.0]))
    out.append({"support": "graphene", "xyz": to_xyz(g), "text": text, "value": value,
                "align": align(g, metal, mpos, donors)})
    t = pd.read_csv(os.path.join(DATA, "xtb_supports.csv"))
    for sup, label in (("MoS2", "on three S of MoS₂ (xTB, 3×3 cut of the sheet)"),
                       ("SiO2", "on the silsesquioxane silica model, two siloxide O and one silanol O (xTB)")):
        r = t[(t.support == sup) & (t.metal == metal) & (t.state == "H_M")]
        if not len(r):
            continue
        r = r.iloc[0]
        site = r.site_H if sup == "MoS2" else "T7"
        atoms = read(os.path.join(DATA, "screen", "xtb_supports", sup, "%s_%s_H_M.xyz" % (metal, site)))
        im = metal_index(atoms, metal)
        if sup == "MoS2":
            p, pm = patch_2d(atoms, im, 1)
            p = oriented(p, pm, np.array([0, 0, 1.0]))
        else:
            cage = atoms.positions[[a.index for a in atoms if a.symbol == "Si"]].mean(axis=0)
            p = oriented(atoms, im, atoms.positions[im] - cage)
        value = "ΔG_H = %+.2f eV (xTB)" % r.dG_H + ("" if pd.isna(r.dG_H_offset) else ", %+.2f eV corrected" % r.dG_H_offset)
        out.append({"support": sup, "xyz": to_xyz(p), "text": "%s %s" % (metal, label), "value": value,
                    "align": align(p, metal, mpos, donors)})
    return out


def main():
    data = []
    for ex in EXAMPLES:
        pdb, metal_spec, donors, resname, mpos = protein(ex)
        metal = ex["site"].split("_")[-1].capitalize()
        data.append({"key": ex["key"], "support": ex["support"], "name": ex["name"], "note": ex["note"], "chains": ex["chains"],
                     "shell": ex["shell"], "metal": metal, "metal_sel": metal_spec, "metal_res": resname,
                     "donors": donors, "pdb": pdb, "supports": supports(ex, metal, mpos, donors)})
        print("%-8s %s: %d atoms, metal %s, %d donors, supports %s" % (
            ex["key"], ex["pdb"], pdb.count("\nATOM") + pdb.count("\nHETATM"), metal_spec, len(donors),
            [(s["support"], s["align"]["rmsd"]) for s in data[-1]["supports"]]))
    blob = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    page = open(os.path.join(HERE, "zoom_template.html")).read().replace("__DATA__", blob)
    out = os.path.join(HERE, "zoom_viewer.html")
    with open(out, "w") as f:
        f.write(page)
    print("wrote %s (%.1f MB)" % (out, os.path.getsize(out) / 1e6))


if __name__ == "__main__":
    main()
