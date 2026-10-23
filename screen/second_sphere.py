"""Second coordination sphere from the enzyme on the graphene site: start structures for a force-field optimisation
followed by DFT single points (dft/make_inputs.py points).

    python screen/second_sphere.py      # -> data/second_sphere/second_sphere_start.extxyz, .csv, README.md

For each site the DFT-relaxed model (bare and with H on the metal) gets the polar groups that surround the metal's
first shell in its source enzyme: the enzyme is laid onto the model by its metal and first-shell donors (rotation
from the best donor assignment, as in site/zoom_viewer.py), and every residue or water with an N, O or S atom
3.0-6.5 A from the metal that is not itself a donor becomes a small neutral molecule on that atom (water; Asp/Glu ->
acetic acid; His -> imidazole; Ser/Thr/Tyr -> methanol; Cys/Met -> methanethiol; backbone, Asn/Gln O ->
formaldehyde; Lys, Arg, backbone, Asn/Gln N -> ammonia), turned to keep clear of the sheet and of each other.
Groups that would sit inside the sheet (the enzyme's donor residues continue where graphene is) or clash are left
out, and the bare and the H_M structure carry the same groups, so dG_H with the second sphere is a clean
difference. Atoms are tagged: 0 = site, k = k-th second-sphere molecule. Starting geometries, not relaxed.
"""
import itertools
import os
import sys

import gemmi
import numpy as np
import pandas as pd
from ase import Atoms
from ase.io import read, write

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
from xtb_free_complexes import ligand, place     # noqa: E402

DATA = os.path.join(ROOT, "data")
OUT = os.path.join(DATA, "second_sphere")
RUNS = os.path.join(ROOT, "dft", "runs", "graphene")
SITES = [("Zn-NS2:SV", "6LAC_A_ZN_201_ZN", "beta-carbonic anhydrase"),
         ("Ni-NS2:SV", "8EZH_B_NI_514_NI", "lactate racemase (Ni pincer cofactor)"),
         ("Co-N2OC:DV:NONC", "3HTW_A_CO_2_CO", "organophosphate hydrolase")]
R_MIN, R_MAX, MAX_GROUPS = 3.0, 6.5, 6
CLEAR_HEAVY, CLEAR_H = 2.6, 1.9          # A from the site's atoms
GROUP = {("HOH", "O"): "water", ("ASP", "OD1"): "acetic acid", ("ASP", "OD2"): "acetic acid",
         ("GLU", "OE1"): "acetic acid", ("GLU", "OE2"): "acetic acid", ("HIS", "ND1"): "imidazole",
         ("HIS", "NE2"): "imidazole", ("SER", "OG"): "methanol", ("THR", "OG1"): "methanol", ("TYR", "OH"): "methanol",
         ("CYS", "SG"): "thiol", ("MET", "SD"): "thiol", ("ASN", "OD1"): "formaldehyde", ("GLN", "OE1"): "formaldehyde",
         ("LYS", "NZ"): "ammonia", ("ARG", "NE"): "ammonia", ("ARG", "NH1"): "ammonia", ("ARG", "NH2"): "ammonia",
         ("ASN", "ND2"): "ammonia", ("GLN", "NE2"): "ammonia", ("TRP", "NE1"): "ammonia"}
AMINO = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO",
         "SER", "THR", "TRP", "TYR", "VAL"}


def molecule(kind):
    if kind == "acetic acid":                    # acetate with H on the second O: neutral
        sym, pos = ligand("acetate")
        o2, c = np.asarray(pos[2]), np.asarray(pos[1])
        return sym + ["H"], pos + [o2 + 0.97 * (o2 - c) / np.linalg.norm(o2 - c)]
    return ligand(kind)


def kabsch(P, Q):
    U, _, Vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    return Vt.T @ np.diag([1, 1, d]) @ U.T


def rotation(model_vecs, enzyme_vecs):
    """Rotation taking enzyme donor directions onto the model's, best assignment (elements ignored)."""
    P = enzyme_vecs / np.linalg.norm(enzyme_vecs, axis=1)[:, None]
    Q = model_vecs / np.linalg.norm(model_vecs, axis=1)[:, None]
    n = min(len(P), len(Q))
    best = None
    for pi in itertools.permutations(range(len(P)), n):
        R = kabsch(P[list(pi)], Q[:n])
        rmsd = float(np.sqrt(((P[list(pi)] @ R.T - Q[:n]) ** 2).sum(axis=1).mean()))
        if best is None or rmsd < best[0]:
            best = (rmsd, R)
    return best[1], best[0]


def enzyme(site_id):
    """Metal position, donor positions and the candidate second-sphere atoms of an enzyme site."""
    pdb, chain, resname, seq, atom = site_id.split("_")
    st = gemmi.read_structure(os.path.join(DATA, "pdb", pdb + ".cif.gz"))
    st.remove_hydrogens()
    st.remove_alternative_conformations()
    model = st[0]
    metal = next(a for ch in model if ch.name == chain for r in ch if r.name == resname and str(r.seqid.num) == seq
                 for a in r if a.name == atom)
    m = np.array(metal.pos.tolist())
    con = pd.read_parquet(os.path.join(DATA, "release", "parquet", "contacts.parquet"),
                          filters=[("site_id", "==", site_id), ("kind", "==", "donor")])
    donors = m + con[["dx", "dy", "dz"]].to_numpy(float)
    donor_res = {(c, int(s)) for c, s in zip(con.auth_chain, con.seqid)}
    cand = {}
    ns = gemmi.NeighborSearch(model, st.cell, 7).populate()
    for mark in ns.find_atoms(metal.pos, "\0", radius=R_MAX):
        cra = mark.to_cra(model)
        r, a = cra.residue, cra.atom
        if (cra.chain.name, r.seqid.num) in donor_res or a.element.name not in ("N", "O", "S"):
            continue
        name = "HOH" if r.is_water() else r.name
        kind = GROUP.get((name, a.name)) or ("formaldehyde" if (name in AMINO and a.name == "O") else
                                            "ammonia" if (name in AMINO and a.name == "N") else None)
        p = np.array(st.cell.find_nearest_pbc_position(metal.pos, mark.pos, mark.image_idx).tolist())   # the image next to the metal
        d = float(np.linalg.norm(p - m))
        if kind is None or d < R_MIN:
            continue
        key = (cra.chain.name, r.seqid.num, name)
        if key not in cand or d < cand[key]["d"]:
            cand[key] = {"pos": p, "d": d, "kind": kind, "label": "%s%d %s" % (name.capitalize(), r.seqid.num, a.name)}
    return m, donors, sorted(cand.values(), key=lambda x: x["d"])


def mic_min(a, b, cell):
    """Smallest periodic distance between point sets a and b (graphene cell, periodic in x and y)."""
    inv = np.linalg.inv(cell)
    d = a[:, None, :] - b[None, :, :]
    f = d @ inv
    f[..., :2] -= np.round(f[..., :2])
    return np.linalg.norm(f @ cell, axis=2).min() if len(a) and len(b) else 9.0


def dress(site, im, groups, R, m_enz):
    """Second-sphere molecules on one model structure -> (atoms, tags, used group labels)."""
    out = site.copy()
    tags = [0] * len(site)
    m_mod = site.positions[im]
    used = []
    for k, g in enumerate(groups, 1):
        anchor = m_mod + R @ (g["pos"] - m_enz)
        sym, pos = molecule(g["kind"])
        best = None
        for spin in np.radians(np.arange(0, 360, 20)):
            p = np.asarray(place(sym, pos, anchor - m_mod, spin)) + m_mod     # template +z points away from the metal
            p += anchor - p[0]
            heavy = np.array([s != "H" for s in sym])
            others = out.positions
            oh = np.array([a.symbol != "H" for a in out])
            clear_h = mic_min(p[heavy], others[oh], out.cell) - CLEAR_HEAVY
            clear_l = mic_min(p[~heavy], others, out.cell) - CLEAR_H if (~heavy).any() else 9.0
            score = min(clear_h, clear_l)
            if best is None or score > best[0]:
                best = (score, p)
        if best[0] < 0:
            continue
        out += Atoms(sym, positions=best[1])
        tags += [k] * len(sym)
        used.append("%s: %s" % (g["label"], g["kind"]))
        if len(used) >= MAX_GROUPS:
            break
    out.set_tags(tags)
    return out, used


def main():
    os.makedirs(OUT, exist_ok=True)
    frames, rows = [], []
    for motif, site_id, enzyme_name in SITES:
        d = motif.replace(":", "_").replace("=", "")
        import json
        meta = json.load(open(os.path.join(RUNS, d, "bare", "meta.json")))
        im, idon = meta["imetal"], meta["idonors"]
        bare = read(os.path.join(RUNS, d, "bare", "relax.out"), format="espresso-out", index=-1)
        hm = read(os.path.join(RUNS, d, "H_M", "relax.out"), format="espresso-out", index=-1)
        for a in (bare, hm):
            a.set_constraint()
        m_enz, donors, cand = enzyme(site_id)
        R, rmsd = rotation(bare.positions[idon] - bare.positions[im], donors - m_enz)
        # the same groups on both structures: keep only those that fit next to H on the metal as well
        b1, used_b = dress(bare, im, cand, R, m_enz)
        h1, used_h = dress(hm, im, cand, R, m_enz)
        common = [g for g in cand if any(u.startswith(g["label"] + ":") for u in used_b)
                  and any(u.startswith(g["label"] + ":") for u in used_h)][:MAX_GROUPS]
        b1, used = dress(bare, im, common, R, m_enz)
        h1, used2 = dress(hm, im, common, R, m_enz)
        if used != used2:                           # a group placed in one but not the other: drop to the common set
            keep = [g for g in common if any(u.startswith(g["label"] + ":") for u in used2)]
            b1, used = dress(bare, im, keep, R, m_enz)
            h1, used2 = dress(hm, im, keep, R, m_enz)
        for state, atoms in (("bare", b1), ("H_M", h1)):
            atoms.info = {"motif": motif, "state": state, "source_site": site_id, "enzyme": enzyme_name,
                          "second_sphere": "; ".join(used), "n_groups": len(used), "metal_index": im,
                          "donor_indices": " ".join(map(str, idon)), "fixed_index": meta["fixed"][0],
                          "align_rmsd_unit": round(rmsd, 3)}
            frames.append(atoms)
        rows.append({"motif": motif, "source_site": site_id, "enzyme": enzyme_name, "candidates": len(cand),
                     "groups": len(used), "second_sphere": "; ".join(used), "atoms_bare": len(b1), "atoms_H_M": len(h1),
                     "align_rmsd_unit": round(rmsd, 3)})
        print("%-16s %s: %d of %d groups (%s), %d atoms" % (motif, site_id, len(used), len(cand), "; ".join(used), len(b1)))
    write(os.path.join(OUT, "second_sphere_start.extxyz"), frames, format="extxyz")
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "second_sphere_start.csv"), index=False)


if __name__ == "__main__":
    main()
