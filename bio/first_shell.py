"""Extract the first coordination shell of each transition-metal atom and map it to a graphene SAC site.

Usage:
    python bio/first_shell.py enzymes     # hits_enzymes.csv entries -> data/enzyme_sites.csv, data/motifs_enzyme.csv
    python bio/first_shell.py all         # hits_all.csv entries     -> data/pdb_sites_all.csv, data/motifs_pdb.csv
    python bio/first_shell.py check       # literature sanity checks -> data/sanity_checks.csv
    python bio/first_shell.py summary     # data/sweep_summary.md (also printed)
    python bio/first_shell.py show 3U7Q   # print the sites of one entry
    python bio/first_shell.py complete    # all d-block sites of data/complete/entries.parquet -> data/complete/
    python bio/first_shell.py complete retry    # process the entries whose last attempt failed once more
    python bio/first_shell.py complete finish   # only re-join labels, motif table, provenance from the parts

First shell: hydrogens removed, per residue the alternative conformer with the highest mean occupancy
is kept (ties: first altloc), neighbours are searched with crystal symmetry (sites on crystal contacts and
special positions are complete). A neighbour is a donor if it lies within CUTOFF[element] of the metal;
other metals within METAL_CUTOFF are listed as metal_neighbors (cluster nuclearity), not donors.
Copies of the same site in several chains of one entry (same metal atom, same donor residues/atoms)
are counted once (column n_copies). Column donors lists "class chain:residue:seqid:atom distance" per
donor (class = <element>_<type>, * = atom of a symmetry mate).

Mapping to a graphene site (the key output):
  DV  double vacancy: four in-plane donors (the four pyridinic/C sites of the vacancy), metal in the plane.
      Among all 4-subsets of the donors take the most square-planar one, score =
      RMS(4 cis angles - 90 deg) + RMS(2 trans angles - 180 deg) over the best trans pairing. A tetrapyrrole
      (4 N_porphyrin) always defines the plane; near-ties (< DV_TIE deg) leave the more labile donors axial.
      The left-over donors are axial (+ax:..., e.g. the proximal His of a heme). 2+2 and 2+1+1
      compositions get an isomer: cis (like donors ~90 deg apart) or trans (~180 deg).
  SV  single vacancy: three donors, metal raised out of the plane, apex open for substrate/H.
      Take the 3-subset forming a face (no angle > FACE_MAX_ANGLE) whose complement holds the most labile
      donors (water, then other O, then the longest bond); score = RMS(3 angles - 90 deg).
  preferred_site = DV if CN >= 4 and dv_score <= DV_SCORE_MAX (square planar, square pyramidal,
      octahedral), else SV if CN >= 3 (tetrahedral, trigonal), else none (CN <= 2).
  graphene_motif = "Ni-N2S2(DV,cis)", "Fe-N4(DV)+ax:N", "Cu-N3(SV)".
"""
import itertools
import json
import multiprocessing as mp
import os
import signal
import sys
import time
from collections import Counter, defaultdict
from datetime import date
from multiprocessing import Pool, get_context
from pathlib import Path

import gemmi
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PDB_DIR = DATA / "pdb"

METALS = ["Fe", "Ni", "Co", "Cu", "Mo", "W", "Mn", "V", "Zn", "Cr"]
CUTOFF = {"N": 2.6, "O": 2.6, "S": 2.8, "Se": 2.9, "C": 2.3, "Cl": 2.8}  # metal-donor, A
METAL_CUTOFF = 3.0     # metal-metal contacts recorded as metal_neighbors
WEAK_EXTRA = 0.5       # cutoff..cutoff+0.5 A listed as weak_contacts (e.g. Cu-S(Met) ~3.1 A in T1 Cu)
MIN_DIST = 1.4         # closer than any real M-L bond (M=O ~1.6 A): self image or clash
METAL_ALT = 2.0        # two metals this close are alternative positions, not a cluster
MINOR_OCC = 0.5        # contacts with occupancy < 0.5 x metal occupancy are minority states -> weak_contacts
OVERLAP_ANGLE = 25.0   # deg; L-M-L below this (side-on O2/N2 >= 32 deg) means two copies of one ligand
ORDER = ["N", "S", "O", "C", "Se", "Cl"]
DV_SCORE_MAX = 55.0    # deg; minimum between the two modes of the dv_score distribution (see summary), ~ tau4 0.5
DV_TIE = 5.0           # deg; DV subsets this close count as equally planar
FACE_MAX_ANGLE = 140.0 # deg; a 3-subset with a wider angle contains a trans pair -> no facial (SV) site
ISOMER_TRANS = 135.0   # deg; like donors further apart than this are trans
N_WORKERS = 12

SIDECHAIN = {("HIS", "ND1"): "N_His", ("HIS", "NE2"): "N_His", ("LYS", "NZ"): "N_Lys",
             ("CYS", "SG"): "S_Cys", ("MET", "SD"): "S_Met", ("SEC", "SE"): "Se_Sec", ("CSE", "SE"): "Se_Sec",
             ("ASP", "OD1"): "O_Asp", ("ASP", "OD2"): "O_Asp", ("GLU", "OE1"): "O_Glu", ("GLU", "OE2"): "O_Glu",
             ("TYR", "OH"): "O_Tyr", ("SER", "OG"): "O_Ser/Thr", ("THR", "OG1"): "O_Ser/Thr"}
# tetrapyrroles whose ring N count as N_porphyrin even if fewer than 4 of them bind (4 always do)
TETRAPYRROLES = {"HEM", "HEC", "HEA", "HEB", "HEO", "HDD", "HAS", "HDM", "HEV", "DHE", "SRM", "F43", "M43",
                 "B12", "COB", "CNC", "BCB", "CLA", "CHL", "BCL", "VER", "HNI", "COH", "ZNH", "MNH", "MNR"}
LABILE_RANK = {"O_water": 3}  # then any other O: 2, everything else: 1; ties -> longer bond is more labile


# ---------------------------------------------------------------- structure preparation

def keep_best_altloc(model):
    """Keep one conformer: per residue the altloc with the highest mean occupancy (ties: first letter);
    ligands modelled as separate residues with different altlocs (e.g. both P-cluster states in 3U7Q)
    are resolved the same way where they overlap."""
    for chain in model:
        groups = defaultdict(list)  # microheterogeneity: several residues share one seqid
        for i, res in enumerate(chain):
            groups[(res.subchain, res.seqid.num, res.seqid.icode)].append(i)
        empty = []
        for idx in groups.values():
            occ = defaultdict(list)
            for i in idx:
                for a in chain[i]:
                    if a.altloc != "\0":
                        occ[a.altloc].append(a.occ)
            if not occ:
                continue
            best = max(sorted(occ), key=lambda k: sum(occ[k]) / len(occ[k]))
            for i in idx:
                res = chain[i]
                for j in reversed(range(len(res))):
                    if res[j].altloc not in ("\0", best):
                        del res[j]
                if len(res) == 0:
                    empty.append(i)
        for i in sorted(empty, reverse=True):
            del chain[i]
    alt_res = [(ci, ri, res[0].altloc, np.mean([a.occ for a in res]), np.array([a.pos.tolist() for a in res]))
               for ci, chain in enumerate(model) for ri, res in enumerate(chain)
               if res.entity_type == gemmi.EntityType.NonPolymer and len({a.altloc for a in res}) == 1
               and res[0].altloc != "\0"]
    kept, drop = [], []
    for r in sorted(alt_res, key=lambda r: (-r[3], r[2])):
        overlaps = any(r[2] != k[2] and np.linalg.norm(r[4][:, None] - k[4][None], axis=2).min() < 1.0 for k in kept)
        (drop if overlaps else kept).append(r)
    for ci, ri, *_ in sorted(drop, key=lambda r: (r[0], r[1]), reverse=True):
        del model[ci][ri]


def bonded(atom, res, max_extra=0.35):
    """Non-metal heavy atoms of the same residue covalently bonded to atom (distance < r1 + r2 + 0.35)."""
    return [b for b in res if b.name != atom.name and not b.element.is_metal and b.element.name not in ("H", "D")
            and atom.pos.dist(b.pos) < atom.element.covalent_r + b.element.covalent_r + max_extra]


def is_dithiolene(s_atom, res):
    """S-C-C-S with two thiolate S (bonded to C only) within the residue: ene-1,2-dithiolate of MGD/MPT/MCD.
    The thiolate test keeps CoM (HS-CH2-CH2-SO3-) out; C=C is not required (some pterins are modelled
    with a long C-C bond)."""
    def thiolate(s):
        return s.element.name == "S" and all(x.element.name == "C" for x in bonded(s, res))
    if not thiolate(s_atom):
        return False
    for c1 in bonded(s_atom, res):
        for c2 in bonded(c1, res):
            if c1.element.name == c2.element.name == "C" and any(
                    s2.name != s_atom.name and thiolate(s2) for s2 in bonded(c2, res)):
                return True
    return False


def donor_class(atom, res):
    el, name = atom.element.name, res.name
    if res.is_water():
        return "O_water"
    if (name, atom.name) in SIDECHAIN:
        return SIDECHAIN[(name, atom.name)]
    if res.entity_type == gemmi.EntityType.Polymer and atom.name in ("N", "O", "OXT"):
        return "N_backbone" if atom.name == "N" else "O_backbone"
    partners = bonded(atom, res)
    if el == "S":
        if not partners:
            return "S_inorganic"
        if is_dithiolene(atom, res):
            return "S_dithiolene"
    if el == "C":
        if not partners and sum(a.element.is_metal for a in res) >= 2:
            return "C_carbide"
        if len(partners) == 1 and partners[0].element.name in ("N", "O") and len(bonded(partners[0], res)) == 1:
            return "C_CO/CN"
    return f"{el}_{name}"


# ---------------------------------------------------------------- geometry

def angle_matrix(vecs):
    u = vecs / np.linalg.norm(vecs, axis=1)[:, None]
    return np.degrees(np.arccos(np.clip(u @ u.T, -1, 1)))


def tau4(ang):
    """Yang, Powell, Houser 2007: 0 square planar, 1 tetrahedral."""
    a = sorted(ang[i, j] for i, j in itertools.combinations(range(4), 2))
    return (360 - a[-1] - a[-2]) / 141


def tau5(ang):
    """Addison 1984: 0 square pyramidal, 1 trigonal bipyramidal."""
    a = sorted(ang[i, j] for i, j in itertools.combinations(range(5), 2))
    return (a[-1] - a[-2]) / 60


def square_planar_score(ang, quad):
    """Best (score, trans pairs) of 4 donors: RMS(cis - 90) + RMS(trans - 180) over the 3 pairings."""
    a, b, c, d = quad
    best = None
    for t1, t2 in [((a, b), (c, d)), ((a, c), (b, d)), ((a, d), (b, c))]:
        trans = np.array([ang[t1], ang[t2]])
        cis = np.array([ang[i, j] for i, j in itertools.combinations(quad, 2) if {i, j} not in ({*t1}, {*t2})])
        s = np.sqrt(np.mean((cis - 90) ** 2)) + np.sqrt(np.mean((trans - 180) ** 2))
        if best is None or s < best[0]:
            best = (s, (t1, t2))
    return best


def formula(elements):
    c = Counter(elements)
    keys = [e for e in ORDER if e in c] + sorted(e for e in c if e not in ORDER)
    return "".join(f"{e}{c[e] if c[e] > 1 else ''}" for e in keys)


def lability(d):
    return (LABILE_RANK.get(d["cls"], 2 if d["el"] == "O" else 1), d["dist"])


def labile_key(donors):
    return sorted((lability(d) for d in donors), reverse=True)


def dv_quad(donors, ang):
    """(score, quad): the most square-planar 4-subset (see module docstring), as donor indices."""
    n = len(donors)
    por = [i for i, d in enumerate(donors) if d["cls"] == "N_porphyrin"]
    quads = [tuple(por[:4])] if len(por) >= 4 else list(itertools.combinations(range(n), 4))
    scored = [(square_planar_score(ang, q), q) for q in quads]
    best = min(s[0][0] for s in scored)
    near = [s for s in scored if s[0][0] <= best + DV_TIE]
    (score, _), quad = max(near, key=lambda s: (labile_key([donors[i] for i in range(n) if i not in s[1]]),
                                                 -s[0][0]))
    return score, quad


def dv_site(donors, ang):
    """Most square-planar 4-subset (see module docstring)."""
    n = len(donors)
    score, quad = dv_quad(donors, ang)
    els = [donors[i]["el"] for i in quad]
    isomer = ""
    twice = [e for e, k in Counter(els).items() if k == 2]
    if twice and len(set(els)) <= 3:  # 2+2 or 2+1+1 composition
        i, j = [q for q in quad if donors[q]["el"] == twice[0]]
        isomer = "trans" if ang[i, j] > ISOMER_TRANS else "cis"
    axial = [donors[i] for i in range(n) if i not in quad]
    return {"dv_formula": formula(els), "dv_isomer": isomer, "dv_score": round(score, 1),
            "dv_tau4": round(tau4(ang[np.ix_(quad, quad)]), 3),
            "dv_classes": " ".join(sorted(donors[i]["cls"] for i in quad)),
            "axial_formula": formula(d["el"] for d in axial),
            "axial_classes": " ".join(sorted(d["cls"] for d in axial))}


def sv_face_score(ang, t):
    return float(np.sqrt(np.mean([(ang[i, j] - 90) ** 2 for i, j in itertools.combinations(t, 2)])))


def sv_tri(donors, ang):
    """(tri, is a face): the 3-subset whose complement holds the most labile donors, as donor indices."""
    n = len(donors)
    triples = list(itertools.combinations(range(n), 3))
    faces = [t for t in triples if max(ang[i, j] for i, j in itertools.combinations(t, 2)) <= FACE_MAX_ANGLE]
    tri = max(faces or triples, key=lambda t: (labile_key([donors[i] for i in range(n) if i not in t]),
                                               -sv_face_score(ang, t)))
    return tri, bool(faces)


def sv_site(donors, ang):
    """3-subset forming a face whose complement holds the most labile donors (see module docstring)."""
    n = len(donors)
    tri, is_face = sv_tri(donors, ang)
    score = lambda t: sv_face_score(ang, t)   # noqa: E731
    faces = [tri] if is_face else []
    removed = [donors[i] for i in range(n) if i not in tri]
    return {"sv_formula": formula(donors[i]["el"] for i in tri), "sv_score": round(score(tri), 1),
            "sv_face": bool(faces), "sv_classes": " ".join(sorted(donors[i]["cls"] for i in tri)),
            "sv_removed": " ".join(sorted(d["cls"] for d in removed))}


# ---------------------------------------------------------------- per entry

def entry_sites(entry_id):
    """All transition-metal sites of one entry (copies in several chains merged)."""
    return merge_copies(raw_sites(entry_id))


def raw_sites(entry_id, metals=METALS, detail=False, stats=None):
    """One site per atom of `metals` in the first model (copies not merged, copy_key still set). detail: add the
    per-site fields and the full shell ("_detail") of the complete sweep; stats (dict) receives the file
    revision, model and atom counts."""
    path = str(PDB_DIR / f"{entry_id}.cif.gz")
    if stats is None:
        st = gemmi.read_structure(path)
    else:
        doc = gemmi.cif.Document()
        st = gemmi.read_structure(path, save_doc=doc)
        stats.update(file_revision(doc), n_models=len(st))
        del doc
    st.remove_hydrogens()
    model = st[0]
    if stats is not None:  # all conformers, before keep_best_altloc
        stats["n_metal_atoms_file"] = sum(a.element.name in metals for ch in model for r in ch for a in r)
    keep_best_altloc(model)
    st.setup_cell_images()
    crystal = st.cell.is_crystal()  # EM/NMR models have a dummy 1x1x1 cell: no symmetry
    cell = st.cell if crystal else gemmi.UnitCell()
    back = None  # complete sweep: rotation and shift from the standard back to the deposited crystal frame
    try:
        ns = gemmi.NeighborSearch(model, cell, 5).populate(include_h=False)
    except RuntimeError as exc:  # a few old entries have a non-standard crystal frame (SCALEn)
        if not detail or "standard orientation" not in str(exc):
            raise
        back = standardize_frame(st)
        model, cell = st[0], st.cell
        ns = gemmi.NeighborSearch(model, cell, 5).populate(include_h=False)
    sub2ent = {sub: ent.name for ent in st.entities for sub in ent.subchains}
    chain_poly = {ch.name: sub2ent.get(ch.get_polymer().subchain_id(), "") for ch in model if len(ch.get_polymer())}
    b_median = float(np.median([a.b_iso for ch in model for r in ch for a in r])) or 1.0
    sites, skipped = [], Counter()
    for ci, chain in enumerate(model):
        for ri, res in enumerate(chain):
            for ai, metal in enumerate(res):
                if metal.element.name not in metals:
                    continue
                if metal.occ <= 0 or is_minor_alternative(model, ns, cell, crystal, metal, (ci, ri, ai)):
                    skipped["zero_occupancy" if metal.occ <= 0 else "minor_alternative"] += 1
                    continue
                site = shell(model, ns, cell, crystal, chain, res, metal, sub2ent, chain_poly, detail)
                site.update(entry_id=entry_id, resolution=st.resolution or float("nan"),
                            metal_b=round(metal.b_iso, 1), metal_b_rel=round(metal.b_iso / b_median, 2),
                            metal_occ=round(metal.occ, 2))
                if detail:
                    xyz = np.array(metal.pos.tolist())
                    if back is not None:
                        xyz = back[0] @ (xyz - back[1])
                    site.update(model=model.num, label_asym_id=res.subchain, auth_seq_id=res.seqid.num,
                                icode=res.seqid.icode.strip(), altloc=metal.altloc.strip("\0"),
                                metal_entity=sub2ent.get(res.subchain, ""), x=round(float(xyz[0]), 3),
                                y=round(float(xyz[1]), 3), z=round(float(xyz[2]), 3))
                sites.append(site)
    if back is not None:
        to_deposited_frame(sites, back[0])
    if stats is not None:
        n_kept = len(sites) + sum(skipped.values())
        stats.update(crystal_symmetry=crystal, crystal_frame_standardized=back is not None,
                     n_atoms_model1=model.count_atom_sites(), n_metal_atoms=n_kept,
                     n_skipped_altloc=stats["n_metal_atoms_file"] - n_kept,
                     n_skipped_zero_occupancy=skipped["zero_occupancy"],
                     n_skipped_minor_alternative=skipped["minor_alternative"])
    return sites


def standardize_frame(st):
    """Move a structure with a non-standard crystal frame (SCALEn of some old entries, which gemmi's neighbour
    search cannot use) to the standard frame; returns (inverse rotation, shift) back to the deposited frame."""
    orth = np.array(gemmi.UnitCell(*st.cell.parameters).orth.mat.tolist())
    rot = orth @ np.array(st.cell.frac.mat.tolist())  # x_std = rot x_dep + shift
    shift = orth @ np.array(st.cell.frac.vec.tolist())
    st.standardize_crystal_frame()
    st.setup_cell_images()
    return np.linalg.inv(rot), shift


def to_deposited_frame(sites, rot_inv):
    """Contact vectors of the shell detail back in the deposited frame (rotated after rounding: <= 0.001 A off)."""
    for s in sites:
        for c in [c for k in ("donors", "weak", "other", "metals") for c in s["_detail"][k]]:
            c["dx"], c["dy"], c["dz"] = (round(float(v), 3) for v in rot_inv @ np.array([c["dx"], c["dy"], c["dz"]]))


def file_revision(doc):
    """Latest 'Structure model' record of _pdbx_audit_revision_history (other content types such as structure
    factors or EM maps have their own numbering): revision "major.minor" and its date."""
    t = doc.sole_block().find("_pdbx_audit_revision_history.",
                              ["ordinal", "major_revision", "minor_revision", "revision_date", "?data_content_type"])
    rows = sorted((int(r[0]), gemmi.cif.as_string(r[1]), gemmi.cif.as_string(r[2]), gemmi.cif.as_string(r[3]))
                  for r in t if not r.has(4) or gemmi.cif.as_string(r[4]) == "Structure model")
    return {"file_revision": f"{rows[-1][1]}.{rows[-1][2]}" if rows else "",
            "file_revision_date": rows[-1][3] if rows else ""}


def is_minor_alternative(model, ns, cell, crystal, metal, idx):
    """Two metals closer than METAL_ALT are alternative positions modelled without altloc (e.g. the C-cluster
    Fe1 of 1OAO, 0.7/0.3): keep the higher occupancy (ties: the first). Own symmetry images do not count."""
    for mark in ns.find_atoms(metal.pos, "\0", radius=METAL_ALT):
        other = mark.to_cra(model).atom
        oidx = (mark.chain_idx, mark.residue_idx, mark.atom_idx)
        if not other.element.is_metal or oidx == idx:
            continue
        pos = cell.find_nearest_pbc_position(metal.pos, other.pos, mark.image_idx) if crystal else other.pos
        if pos.dist(metal.pos) < METAL_ALT and (other.occ > metal.occ or (other.occ == metal.occ and oidx < idx)):
            return True
    return False


def contact_ref(cra, mark, pos, metal, dist, sub2ent):
    """Identity and metal-relative position (symmetry applied) of a contact atom, for shell_json."""
    r, v = cra.residue, pos - metal.pos
    return {"el": cra.atom.element.name, "auth_chain": cra.chain.name, "label_asym": r.subchain, "resname": r.name,
            "seqid": r.seqid.num, "icode": r.seqid.icode.strip(), "atom": cra.atom.name, "dist": round(dist, 3),
            "dx": round(v.x, 3), "dy": round(v.y, 3), "dz": round(v.z, 3), "occ": round(cra.atom.occ, 2),
            "b": round(cra.atom.b_iso, 2), "sym": mark.image_idx,
            "polymer": r.entity_type == gemmi.EntityType.Polymer, "entity": sub2ent.get(r.subchain, "")}


def shell(model, ns, cell, crystal, chain, res, metal, sub2ent, chain_poly, detail=False):
    donors, metals, weak, seen = [], [], [], set()
    other, metal_refs = [], []  # detail only: atoms of elements without a cutoff, metal neighbour atoms
    # a half-occupied metal off any special position belongs to a ligand disordered across a symmetry
    # element (e.g. decavanadate): half-occupied atoms of symmetry images are then the other orientation
    disordered = crystal and metal.occ <= 0.5 and not cell.is_special_position(metal.pos, 0.5)
    for mark in ns.find_atoms(metal.pos, "\0", radius=max(CUTOFF.values()) + WEAK_EXTRA):
        cra = mark.to_cra(model)
        atom = cra.atom
        pos = cell.find_nearest_pbc_position(metal.pos, atom.pos, mark.image_idx) if crystal else atom.pos
        dist = pos.dist(metal.pos)
        key = (cra.chain.name, str(cra.residue.seqid), atom.name, round(pos.x, 1), round(pos.y, 1), round(pos.z, 1))
        if dist < MIN_DIST or atom.occ <= 0 or key in seen:
            continue
        seen.add(key)
        el = atom.element.name
        minor = atom.occ < MINOR_OCC * metal.occ  # minority state (e.g. a 0.3-occupancy water)
        if atom.element.is_metal:
            if METAL_ALT <= dist <= METAL_CUTOFF and not minor:
                metals.append((el, dist))
                if detail:
                    metal_refs.append(contact_ref(cra, mark, pos, metal, dist, sub2ent))
            continue
        if el not in CUTOFF or dist > CUTOFF[el] + WEAK_EXTRA:
            if detail and el not in CUTOFF:
                other.append({**contact_ref(cra, mark, pos, metal, dist, sub2ent), "reason": "element_not_counted"})
            continue
        d = {"el": el, "res": cra.residue.name, "atom": atom.name, "chain": cra.chain.name,
             "seqid": str(cra.residue.seqid), "dist": dist, "vec": np.array((pos - metal.pos).tolist()), "pos": pos,
             "cls": donor_class(atom, cra.residue), "sym": mark.image_idx != 0, "occ": atom.occ, "minor": minor,
             "cov_r": atom.element.covalent_r, "polymer": cra.residue.entity_type == gemmi.EntityType.Polymer,
             "entity": sub2ent.get(cra.residue.subchain, ""),
             "ring_n": el == "N" and [b.element.name for b in bonded(atom, cra.residue)] == ["C", "C"],
             "rkey": (cra.chain.name, str(cra.residue.seqid), cra.residue.name, mark.image_idx)}
        d["alt"] = disordered and d["sym"] and atom.occ <= 0.5
        # protein carbons closer than 2.3 A are model strain/clashes, never CO/CN/carbide ligands
        bond = dist <= CUTOFF[el] and not minor and not d["alt"] and not (el == "C" and d["polymer"])
        if detail:
            d["ref"] = contact_ref(cra, mark, pos, metal, dist, sub2ent)
            d["why"] = [r for r, hit in (("beyond_cutoff", dist > CUTOFF[el]), ("minor_occupancy", minor),
                                         ("alt_disorder", d["alt"]), ("protein_carbon", el == "C" and d["polymer"]))
                        if hit]
        (donors if bond else weak).append(d)
    # two donors < OVERLAP_ANGLE apart are alternatives (e.g. a ligand disordered across a special position
    # and its symmetry copy): keep the original, then the shorter bond
    kept = []
    for d in sorted(donors, key=lambda d: (d["sym"], d["dist"])):
        if any(np.degrees(np.arccos(np.clip(d["vec"] @ k["vec"] / d["dist"] / k["dist"], -1, 1))) < OVERLAP_ANGLE
               for k in kept):
            d["alt"] = True
            d["why"] = ["alt_overlap"]
            weak.append(d)
        else:
            kept.append(d)
    donors = kept
    # a His binds through one ring N only: keep the closer of ND1/NE2 (both within cutoff = misplaced metal)
    his = defaultdict(list)
    for d in donors:
        if d["cls"] == "N_His":
            his[d["rkey"]].append(d)
    for pair in his.values():
        for d in sorted(pair, key=lambda d: d["dist"])[1:]:
            donors.remove(d)
            d["why"] = ["his_second_ring_n"]
            weak.append(d)
    # the far end of a bound ligand (O of CO, N of CN) is not a contact of its own
    weak = [w for w in weak if not any(w["rkey"] == d["rkey"] and d["dist"] < w["dist"] and
                                       w["pos"].dist(d["pos"]) < w["cov_r"] + d["cov_r"] + 0.35 for d in donors)]
    mark_porphyrin(donors)
    donors.sort(key=lambda d: d["dist"])
    n = len(donors)
    ang = angle_matrix(np.array([d["vec"] for d in donors])) if n else None
    ents = Counter(d["entity"] for d in donors if d["polymer"])
    site = {"metal": metal.element.name, "chain": chain.name, "metal_res": res.name, "seqid": str(res.seqid),
            "metal_atom": metal.name, "cn": n, "formula": formula(d["el"] for d in donors),
            "donor_classes": " ".join(sorted(d["cls"] for d in donors)),
            "donors": "; ".join(f"{d['cls']} {d['chain']}:{d['res']}:{d['seqid']}:{d['atom']}{'*' if d['sym'] else ''} "
                                f"{d['dist']:.2f}" for d in donors),
            "mean_dist": round(float(np.mean([d["dist"] for d in donors])), 3) if n else float("nan"),
            "n_metal_neighbors": len(metals),
            "metal_neighbors": " ".join(f"{e}{x:.2f}" for e, x in sorted(metals, key=lambda m: m[1])),
            "weak_contacts": " ".join(f"{d['cls']}:{d['dist']:.2f}" + (f"(occ{d['occ']:.2f})" if d["minor"] else "")
                                      + ("(alt)" if d.get("alt") else "")
                                      for d in sorted(weak, key=lambda d: d["dist"])),
            "sym_donors": sum(d["sym"] for d in donors),
            "tau4": round(tau4(ang), 3) if n == 4 else float("nan"),
            "tau5": round(tau5(ang), 3) if n == 5 else float("nan"),
            "host_entity": ents.most_common(1)[0][0] if ents else chain_poly.get(chain.name, ""),
            "copy_key": (metal.element.name, res.name, metal.name,
                         tuple(sorted((d["cls"], d["res"], d["atom"], d["seqid"] if d["polymer"] else "")
                                      for d in donors)),
                         "" if ents else anchor(model, ns, cell, crystal, res, metal))}
    site.update(dv_site(donors, ang) if n >= 4 else {})
    site.update(sv_site(donors, ang) if n >= 3 else {})
    if n >= 4 and site["dv_score"] <= DV_SCORE_MAX:
        site["preferred_site"] = "DV"
        iso = f",{site['dv_isomer']}" if site["dv_isomer"] else ""
        ax = f"+ax:{site['axial_formula']}" if site["axial_formula"] else ""
        site["graphene_motif"] = f"{site['metal']}-{site['dv_formula']}(DV{iso}){ax}"
    elif n >= 3:
        site["preferred_site"] = "SV"
        site["graphene_motif"] = f"{site['metal']}-{site['sv_formula']}(SV)"
    else:
        site["preferred_site"] = "none"
        site["graphene_motif"] = f"{site['metal']}-{site['formula'] or '0'}(none)"
    site["raw_motif"] = f"{site['metal']}-{site['formula'] or '0'}"
    if detail:
        def ref(d):  # the donor class is final only now (mark_porphyrin)
            r = dict(d["ref"])
            return {"el": r.pop("el"), "cls": d["cls"], **r}
        site["_detail"] = {"donors": [ref(d) for d in donors],
                           "weak": [{**ref(d), "reason": "+".join(d["why"])}
                                    for d in sorted(weak, key=lambda d: d["dist"])],
                           "other": sorted(other, key=lambda o: o["dist"]),
                           "metals": sorted(metal_refs, key=lambda m: m["dist"])}
    return site


def anchor(model, ns, cell, crystal, res, metal):
    """Lowest-numbered protein residue within 0.5 A of the nearest one (<= 5 A), else the metal's own residue:
    tells apart sites without protein donors (e.g. two Mn(H2O)6) when copies are merged; the 0.5 A window
    keeps near-ties from splitting copies."""
    near = {}
    for mark in ns.find_atoms(metal.pos, "\0", radius=5.0):
        cra = mark.to_cra(model)
        if cra.residue.entity_type == gemmi.EntityType.Polymer:
            pos = cell.find_nearest_pbc_position(metal.pos, cra.atom.pos, mark.image_idx) if crystal else cra.atom.pos
            key = (cra.residue.seqid.num, cra.residue.seqid.icode, cra.residue.name)
            near[key] = min(near.get(key, 9.9), pos.dist(metal.pos))
    if not near:
        return (res.seqid.num, res.seqid.icode, res.name)
    d_min = min(near.values())
    return min(k for k, d in near.items() if d <= d_min + 0.5)


def mark_porphyrin(donors):
    """Ring N (bonded to two C) of a tetrapyrrole: >= 4 from one ligand (most square-planar 4) or from a known
    tetrapyrrole."""
    groups = defaultdict(list)
    for i, d in enumerate(donors):
        if d["el"] == "N" and not d["polymer"] and d["ring_n"]:
            groups[d["rkey"]].append(i)
    for idx in groups.values():
        if len(idx) >= 4:
            ang = angle_matrix(np.array([donors[i]["vec"] for i in idx]))
            quad = min(itertools.combinations(range(len(idx)), 4), key=lambda q: square_planar_score(ang, q)[0])
            idx = [idx[i] for i in quad]
        elif donors[idx[0]]["res"] not in TETRAPYRROLES:
            continue
        for i in idx:
            donors[i]["cls"] = "N_porphyrin"


def merge_copies(sites):
    """Count identical sites from several copies of a chain once; keep the copy with the lowest metal B."""
    groups = defaultdict(list)
    for s in sites:
        groups[s.pop("copy_key")].append(s)
    out = []
    for copies in groups.values():
        best = min(copies, key=lambda s: s["metal_b"])
        best["n_copies"] = len(copies)
        best["chains"] = " ".join(sorted({s["chain"] for s in copies}))
        out.append(best)
    return out


def safe_sites(entry_id):
    try:
        return entry_sites(entry_id)
    except Exception as exc:  # one broken file must not stop the sweep
        print(f"  {entry_id}: {type(exc).__name__}: {exc}", flush=True)
        return []


# ---------------------------------------------------------------- tables

def read_csv(path):
    return pd.read_csv(path, dtype={"entry_id": str, "entity_id": str, "host_entity": str, "cluster30": str},
                       keep_default_na=False, na_values=[""])


def label_sites(df, hits):
    """Join enzyme/EC/reaction labels: host entity first, else the union over the entry (label_source)."""
    ents = read_csv(DATA / "entities.csv").fillna("")
    ent = ents.set_index(["entry_id", "entity_id"])
    keys = ["entry_id", "reaction", "enzyme", "query_ec"]
    labelled = pd.concat([ents[ents.reaction != ""][keys], hits[keys] if "reaction" in hits else None]).fillna("")
    entry_lab = {e: {k: ";".join(sorted({x for v in g[k] if v for x in v.split(";")})) for k in keys[1:]}
                 for e, g in labelled.groupby("entry_id")}
    rows = []
    for s in df.itertuples():
        key = (s.entry_id, s.host_entity)
        r = ent.loc[key] if key in ent.index else None
        lab = {"description": r.description if r is not None else "", "ec": r.ec if r is not None else "",
               "organism": r.organism if r is not None else "",
               "cluster30": str(r.cluster30) if r is not None else f"{s.entry_id}_{s.host_entity}"}
        if r is not None and r.reaction:
            lab.update(reaction=r.reaction, enzyme=r.enzyme, query_ec=r.query_ec, label_source="entity")
        elif s.entry_id in entry_lab and entry_lab[s.entry_id]["reaction"]:
            lab.update(entry_lab[s.entry_id], label_source="entry")
        else:
            lab.update(reaction="", enzyme="", query_ec="", label_source="")
        rows.append(lab)
    return pd.concat([df.reset_index(drop=True), pd.DataFrame(rows)], axis=1)


def motif_table(df):
    """Aggregate per graphene motif and per raw first-shell motif."""
    out = []
    for kind, col in [("graphene", "graphene_motif"), ("first_shell", "raw_motif")]:
        for motif, g in df.groupby(col):
            g = g.sort_values("resolution")
            best = g.drop_duplicates("entry_id")
            out.append({
                "kind": kind, "motif": motif, "metal": g.metal.iloc[0], "n_sites": len(g),
                "n_entries": g.entry_id.nunique(), "n_sequence_clusters": g.cluster30.nunique(),
                "reactions": ";".join(sorted({x for v in g.reaction if v for x in v.split(";")})),
                "example_enzymes": " | ".join(dict.fromkeys(d for d in best.description if d))[:200]
                if "description" in g else "",
                "example_pdb": " ".join(f"{e}({r:.2f})" for e, r in zip(best.entry_id[:6], best.resolution[:6])),
                "best_resolution": g.resolution.min(),
                "median_tau4": round(g.tau4.median(), 3), "median_dv_tau4": round(g.dv_tau4.median(), 3),
                "median_dv_score": round(g.dv_score.median(), 1), "median_sv_score": round(g.sv_score.median(), 1),
                "frac_mononuclear": round((g.n_metal_neighbors == 0).mean(), 2),
                "donor_classes_top": g.donor_classes.value_counts().index[0]})
    m = pd.DataFrame(out)
    return m.sort_values(["kind", "n_entries", "n_sequence_clusters", "n_sites"], ascending=False)  # graphene first


SITE_COLS = ["entry_id", "metal", "chain", "metal_res", "seqid", "metal_atom", "cn", "formula", "raw_motif",
             "graphene_motif", "preferred_site", "dv_formula", "dv_isomer", "dv_score", "dv_tau4", "dv_classes",
             "axial_formula", "axial_classes", "sv_formula", "sv_score", "sv_face", "sv_classes", "sv_removed",
             "tau4", "tau5", "mean_dist", "donor_classes", "donors", "n_metal_neighbors", "metal_neighbors",
             "weak_contacts", "sym_donors", "metal_b", "metal_b_rel", "metal_occ", "resolution", "n_copies",
             "chains", "host_entity", "cluster30", "reaction", "enzyme", "query_ec", "ec", "description",
             "organism", "label_source"]


def run_set(hits_file, sites_file, motifs_file):
    hits = read_csv(DATA / hits_file)
    if "selected" in hits:
        hits = hits[hits.selected.astype(str) == "True"]
    entries = sorted(e for e in set(hits.entry_id) if (PDB_DIR / f"{e}.cif.gz").exists())
    print(f"{len(entries)} entries ({hits.entry_id.nunique() - len(entries)} missing on disk)", flush=True)
    with Pool(N_WORKERS) as pool:
        sites = [s for chunk in pool.imap(safe_sites, entries, chunksize=4) for s in chunk]
    df = label_sites(pd.DataFrame(sites), hits)[SITE_COLS]
    df.to_csv(DATA / sites_file, index=False)
    motif_table(df).to_csv(DATA / motifs_file, index=False)
    print(f"{len(df)} unique sites in {df.entry_id.nunique()} entries -> data/{sites_file}, data/{motifs_file}")
    print(df.metal.value_counts().to_string())
    return df


# ---------------------------------------------------------------- sanity checks against the literature

def n_cls(site, cls):
    return site["donor_classes"].split().count(cls)


def n_el(site, el):
    return sum(c.startswith(el + "_") for c in site["donor_classes"].split())


def is_t1(s):
    return s["metal"] == "Cu" and n_cls(s, "S_Cys") == 1


def t1_ok(s):
    return n_cls(s, "N_His") == 2 and set(s["donor_classes"].split()) <= {"N_His", "S_Cys", "S_Met"}


def acs_ok(s):
    return (s["dv_formula"] == "N2S2" and n_cls(s, "N_backbone") == 2 and n_cls(s, "S_Cys") == 2
            and s["dv_isomer"] == "cis" and s["preferred_site"] == "DV")


def mcr_ok(s):
    return s["dv_classes"] == " ".join(["N_porphyrin"] * 4) and "O" in s["axial_formula"]


def fdh_ok(s):
    return n_cls(s, "S_dithiolene") == 4 and n_cls(s, "Se_Sec") == 1


# (check, entry, expectation, site selector, per-site test, test on the number of selected sites)
CHECKS = [
    ("[NiFe]-H2ase DvMF Ni (Ni-R)", "4U9H", "Ni: S4 = 4 Cys (2 bridging)", lambda s: s["metal"] == "Ni",
     lambda s: s["formula"] == "S4" and n_cls(s, "S_Cys") == 4, None),
    ("[NiFe]-H2ase DvMF Fe (Ni-R)", "4U9H", "Fe: C3S2 = 2 CN + CO + 2 bridging Cys",
     lambda s: s["metal"] == "Fe" and n_cls(s, "C_CO/CN") > 0,
     lambda s: s["formula"] == "S2C3" and n_cls(s, "C_CO/CN") == 3 and n_cls(s, "S_Cys") == 2, None),
    ("[NiFe]-H2ase DvMF Ni (Ni-A, alt)", "1WUI", "Ni: S4 = 4 Cys (2 bridging)", lambda s: s["metal"] == "Ni",
     lambda s: s["formula"] == "S4" and n_cls(s, "S_Cys") == 4, None),
    ("[FeFe]-H2ase CpI H-cluster Fe", "3C8Y", "Fe_p/Fe_d: CO/CN carbon + S donors",
     lambda s: s["metal"] == "Fe" and n_cls(s, "C_CO/CN") > 0,
     lambda s: n_cls(s, "C_CO/CN") >= 2 and n_el(s, "S") >= 2, lambda k: k == 2),
    ("Nitrogenase MoFe Mo", "3U7Q", "Mo: S3 N O2 (His, homocitrate)", lambda s: s["metal"] == "Mo",
     lambda s: s["formula"] == "NS3O2" and n_cls(s, "N_His") == 1 and n_cls(s, "S_inorganic") == 3, None),
    ("Nitrogenase FeMo-co belt Fe", "3U7Q", "6 Fe: S3 C (inorganic S, carbide)",
     lambda s: s["metal"] == "Fe" and n_cls(s, "C_carbide") > 0,
     lambda s: s["formula"] == "S3C" and n_cls(s, "S_inorganic") == 3, lambda k: k % 6 == 0),
    ("ACS distal Ni", "1RU3", "Ni: N2S2 (2 backbone N, 2 Cys), cis", lambda s: s["metal"] == "Ni"
     and n_cls(s, "N_backbone") > 0, acs_ok, None),
    ("ACS distal Ni (alt)", "1OAO", "Ni: N2S2 (2 backbone N, 2 Cys), cis", lambda s: s["metal"] == "Ni"
     and n_cls(s, "N_backbone") > 0, acs_ok, None),
    ("CcO heme a3 Fe", "1V54", "Fe: N4 porphyrin + N(His) axial",
     lambda s: s["metal"] == "Fe" and n_cls(s, "N_porphyrin") == 4 and n_cls(s, "N_His") == 1,
     lambda s: s["dv_formula"] == "N4" and "N_His" in s["axial_classes"] and s["preferred_site"] == "DV", None),
    ("CcO CuB", "1V54", "CuB: N3 (3 His)", lambda s: s["metal"] == "Cu" and n_cls(s, "N_His") >= 3,
     lambda s: n_cls(s, "N_His") == 3 and n_el(s, "S") == 0, None),
    ("CcO heme a3 Fe (alt)", "5B1A", "Fe: N4 porphyrin + N(His) axial",
     lambda s: s["metal"] == "Fe" and n_cls(s, "N_porphyrin") == 4 and n_cls(s, "N_His") == 1,
     lambda s: s["dv_formula"] == "N4" and "N_His" in s["axial_classes"] and s["preferred_site"] == "DV", None),
    ("CcO CuB (alt)", "5B1A", "CuB: N3 (3 His)", lambda s: s["metal"] == "Cu" and n_cls(s, "N_His") >= 3,
     lambda s: n_cls(s, "N_His") == 3 and n_el(s, "S") == 0, None),
    ("Laccase T1 Cu", "LACCASE", "T1 Cu: N2S (2 His, Cys) (+S_Met)", is_t1, t1_ok, None),
    ("Laccase T3 Cu", "LACCASE", "T3 Cu: N3 (3 His)", lambda s: s["metal"] == "Cu" and n_cls(s, "N_His") >= 3,
     lambda s: n_cls(s, "N_His") == 3, lambda k: k >= 2),
    ("Ascorbate oxidase T1 Cu (Met)", "1AOZ", "T1 Cu: N2S + S_Met", is_t1,
     lambda s: t1_ok(s) and "S_Met" in s["donor_classes"] + s["weak_contacts"], None),
    ("MCR Ni F430", "1MRO", "Ni: N4 (F430) + O axial", lambda s: s["metal"] == "Ni", mcr_ok, None),
    ("MCR Ni F430 (alt)", "1HBN", "Ni: N4 (F430) + O axial", lambda s: s["metal"] == "Ni", mcr_ok, None),
    ("Carbonic anhydrase II Zn", "3KS3", "Zn: N3 (His) + O_water", lambda s: s["metal"] == "Zn",
     lambda s: s["formula"] == "N3O" and n_cls(s, "N_His") == 3 and n_cls(s, "O_water") == 1, None),
    ("FDH-H Mo (oxidized)", "1FDO", "Mo: S4 dithiolene + Se (SeCys) (+S/O)", lambda s: s["metal"] == "Mo",
     fdh_ok, None),
    ("FDH-H Mo (reduced, alt)", "1AA6", "Mo: S4 dithiolene + Se (SeCys) (+S/O)", lambda s: s["metal"] == "Mo",
     fdh_ok, None),
]


CHECK_NOTES = {
    "[NiFe]-H2ase DvMF Ni (Ni-A, alt)": "Ni-A (oxidised, unready) model: Cys84/Cys546 partly oxidised to "
    "sulfenates (CSO) and bridging O species superposed at partial occupancy; the structure differs, S4 holds "
    "for the active Ni-R state (4U9H)",
}


def best_laccase():
    """Best-resolution EC 1.10.3.2 hit with a T1 copper."""
    hits = read_csv(DATA / "hits_enzymes.csv")
    for e in hits[hits.query_ec.str.contains("1.10.3.2", regex=False)].sort_values("resolution").entry_id:
        if any(is_t1(s) for s in entry_sites(e)):
            return e


def run_checks():
    from pdb_sweep import fetch_structures  # test entries not among the hits are downloaded on demand
    fetch_structures([c[1] for c in CHECKS if c[1] != "LACCASE"])
    laccase, cache, rows = best_laccase(), {}, []
    for name, entry, expected, select, test, count_ok in CHECKS:
        entry = laccase if entry == "LACCASE" else entry
        sites = cache.setdefault(entry, entry_sites(entry))
        sel = [s for s in sites if select(s)]
        ok = bool(sel) and all(test(s) for s in sel) and (count_ok is None or count_ok(len(sel)))
        observed = "; ".join(
            f"{s['metal']} {s['metal_res']} {s['seqid']}/{s['metal_atom']}: {s['formula']} [{s['donor_classes']}]"
            f" -> {s['graphene_motif']}" + (f" (weak: {s['weak_contacts']})" if s["weak_contacts"] else "")
            for s in sel) or "no matching site"
        rows.append({"check": name, "entry": entry, "expected": expected, "observed": observed,
                     "n_sites": len(sel), "passed": ok, "note": CHECK_NOTES.get(name, "")})
        print(f"{'PASS' if ok else 'FAIL'}  {name:34s} {entry}  {observed[:220]}")
    pd.DataFrame(rows).to_csv(DATA / "sanity_checks.csv", index=False)


# ---------------------------------------------------------------- summary

def md_table(df):
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    lines += ["| " + " | ".join(str(v) for v in row) + " |" for row in df.itertuples(index=False)]
    return "\n".join(lines)


def split_count(series_by_entry, sites):
    """Entries and sites per reaction label (labels joined by ';' count for each class)."""
    out = defaultdict(lambda: [set(), 0, set()])
    for e, lab in series_by_entry:
        for r in str(lab).split(";"):
            out[r][0].add(e)
    for s in sites.itertuples():
        for r in str(s.reaction).split(";") if s.reaction else []:
            out[r][1] += 1
            out[r][2].add(s.entry_id)
    return out


def dv_threshold_report(sites):
    """KDE minimum of the bimodal dv_score distribution (CN >= 4) and agreement with tau4(best 4) <= 0.5."""
    from scipy.stats import gaussian_kde
    d = sites[sites.cn >= 4].drop_duplicates(["entry_id", "chain", "metal_res", "seqid", "metal_atom"])
    grid = np.arange(20, 80, 0.5)
    dens = gaussian_kde(d.dv_score, bw_method=0.08)(grid)
    agree = ((d.dv_score <= DV_SCORE_MAX) == (d.dv_tau4 <= 0.5)).mean()
    return len(d), float(grid[dens.argmin()]), float(agree)


def write_summary():
    counts = json.loads((DATA / "sweep_counts.json").read_text())
    hits = read_csv(DATA / "hits_enzymes.csv")
    es, ps = read_csv(DATA / "enzyme_sites.csv"), read_csv(DATA / "pdb_sites_all.csv")
    es["reaction"] = es.reaction.fillna("")
    motifs, checks = read_csv(DATA / "motifs_enzyme.csv"), read_csv(DATA / "sanity_checks.csv").fillna("")
    pdb_motifs = read_csv(DATA / "motifs_pdb.csv").fillna("")
    per_rxn = split_count(zip(hits.entry_id, hits.reaction), es)
    rxn_tab = pd.DataFrame([{"reaction": r, "hit entries": len(v[0]), "entries with TM sites": len(v[2]),
                             "unique TM sites": v[1]} for r, v in sorted(per_rxn.items())])
    ec_tab = pd.DataFrame([{"reaction": c["reaction"], "EC": ec, "enzyme": c["enzyme"],
                            "entities <=2.5 A": c["entities_total"], "95% representatives": c["representatives"],
                            "entries": c["entries"]} for ec, c in counts["enzymes"]["per_ec"].items()])
    empty = [f"{ec} ({counts['enzymes']['per_ec'][ec].get('entities_any_resolution', '?')} entities at any "
             f"resolution)" for ec in counts["enzymes"]["empty_ec"]]
    metal_tab = pd.DataFrame([{"metal": m, "CCD components": c["components"],
                               "protein entities <=2.0 A": c["protein_entities"], "30% clusters": c["clusters_30"],
                               "representative entries": c["entries"],
                               "sites (metal) in sweep": int((ps.metal == m).sum())}
                              for m, c in counts["all"]["per_metal"].items()])
    graphene = motifs[motifs.kind == "graphene"]
    top = graphene[~graphene.motif.str.endswith("(none)")].head(25)[
        ["motif", "n_sites", "n_entries", "n_sequence_clusters", "reactions", "frac_mononuclear", "example_pdb"]]
    pdb_graphene = pdb_motifs[(pdb_motifs.kind == "graphene") & ~pdb_motifs.motif.str.endswith("(none)")]
    top_pdb = pdb_graphene.head(15)[["motif", "n_sites", "n_entries", "n_sequence_clusters", "frac_mononuclear",
                                     "example_pdb"]]
    chk = checks[["check", "entry", "expected", "passed"]].assign(
        observed=checks.observed.str.slice(0, 160), note=checks.note)
    n_dv, n_sv, n_none = ((es.preferred_site == k).sum() for k in ("DV", "SV", "none"))
    n_kde, kde_min, agree = dv_threshold_report(pd.concat([es, ps], ignore_index=True))
    text = f"""# PDB sweep: transition-metal first shells for graphene SAC motifs

Generated {date.today()} by `bio/pdb_sweep.py` (search, download) and `bio/first_shell.py` (extraction, motifs,
checks) from the RCSB PDB state of that day. All outputs are drafts that need expert review.

## Searches
* Enzyme search (RCSB Search API v2): EC lineage per class, resolution <= {2.5} A, polymer entities grouped at
  95 % sequence identity, best resolution as representative: **{counts['enzymes']['entities']} representative
  entities in {counts['enzymes']['entries']} entries**.
* EC classes without any hit at <= 2.5 A: {', '.join(empty) if empty else 'none'}.
* Whole-PDB sweep: entries with a bound component containing Fe/Ni/Co/Cu/Mo/W/Mn/V, resolution <= 2.0 A,
  protein entities grouped at 30 % identity: **{counts['all']['entries']} unique representative entries**
  ({'below the 6000 cap, all' if counts['all']['entries'] <= 6000 else 'capped at 6000 best-resolution'}
  selected and downloaded). The chemical-service formula search (match_subset) returned only 130 Fe
  components and missed SF4/FES/F3S/ICS/CLF, so the metal component ids were taken from the wwPDB Chemical
  Component Dictionary; rcsb_entry_info.nonpolymer_bound_components is not searchable, the filter uses
  rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id.

{md_table(ec_tab)}

{md_table(metal_tab)}

## Sites
Enzyme set: {len(es)} unique transition-metal sites in {es.entry_id.nunique()} entries ({n_dv} DV, {n_sv} SV,
{n_none} none).
Whole-PDB set: {len(ps)} unique sites in {ps.entry_id.nunique()} entries.

{md_table(rxn_tab)}

## Top 25 graphene motifs, enzyme set (sorted by n_entries; the {n_none} "none" sites, CN <= 2, not listed)
{md_table(top)}

## Top 15 graphene motifs, whole-PDB sweep ({(ps.preferred_site == 'none').sum()} "none" sites not listed)
{md_table(top_pdb)}

## Sanity checks ({int(checks.passed.astype(str).eq('True').sum())}/{len(checks)} passed)
{md_table(chk)}

## Assumptions and thresholds
* Donor cutoffs (A): {', '.join(f'{k} {v}' for k, v in CUTOFF.items())}; metal-metal <= {METAL_CUTOFF} A recorded
  as metal_neighbors; F/Br/I and other elements are not counted as donors. Contacts up to cutoff + {WEAK_EXTRA} A
  are listed as weak_contacts (e.g. Cu-S(Met) 2.9-3.3 A of T1 Cu), not as donors.
* Hydrogens removed; per residue the altloc with the highest mean occupancy is kept (ties: first letter); ligands
  modelled as separate residues with different altlocs are resolved the same way where they overlap (< 1 A).
* Metals closer than {METAL_ALT} A to another metal are alternative positions: only the higher-occupancy one is
  kept. Contacts with occupancy < {MINOR_OCC} x metal occupancy are minority states (weak_contacts).
* Protein carbon atoms are never donors; of the two His ring N only the closer one can bind; two donors less
  than {OVERLAP_ANGLE} deg apart (L-M-L) are copies of one ligand (the original is kept); for a half-occupied
  metal off special positions, half-occupied atoms of symmetry mates are the other orientation of a ligand
  disordered across a symmetry element (both go to weak_contacts as "alt").
* Symmetry: neighbour search with all crystal symmetry operators (sym_donors counts donors from images);
  EM entries without a crystal cell are searched without symmetry.
* Donor classes from residue/atom names; S_inorganic = S without non-metal neighbour, S_dithiolene =
  S-C-C-S with two thiolate S, C_CO/CN = C bonded to one terminal N/O, C_carbide = C without non-metal
  neighbour in a multi-metal cluster, N_porphyrin = >= 4 ring N (bonded to two C) of one ligand around the
  metal (or of a known tetrapyrrole); everything else <element>_<residue>.
* Copies of a site in several chains (same metal atom, same donor residues/atoms; without protein donors also
  the same nearest protein residue) are counted once (n_copies).
* DV: most square-planar 4-subset, score = RMS(cis - 90) + RMS(trans - 180) deg; a tetrapyrrole always defines
  the plane; subsets within {DV_TIE} deg of the best leave the more labile donors axial. DV is preferred when
  CN >= 4 and score <= {DV_SCORE_MAX} deg, set at the minimum between the two modes (square planar/octahedral
  ~12 deg, tetrahedral ~87 deg) of the score distribution: the KDE minimum of the {n_kde} CN >= 4 sites of both
  sets lies at {kde_min:.1f} deg, and the split agrees with tau4(best 4-subset) <= 0.5 for {agree:.1%} of them.
* SV (CN >= 3): 3-subset forming a face (all angles <= {FACE_MAX_ANGLE} deg) whose complement holds the most
  labile donors (water > other O > longest bond); score = RMS(angles - 90) deg. CN <= 2: none.
* Isomer for 2+2 and 2+1+1 DV compositions: trans if the like donors are > {ISOMER_TRANS} deg apart.
* Reaction/enzyme labels: EC lineage of the polymer entity that provides most protein donors (label_source =
  entity), else the union over the entry (label_source = entry). n_sequence_clusters uses RCSB 30 % clusters.

## Limitations
* Representatives only: 95 % (enzymes) and 30 % (sweep) clusters, chosen by best resolution; the sweep
  representative of a cluster need not be the chain that binds the metal (entry-level component filter).
* Formate dehydrogenase classes 1.17.1.9/1.2.1.2 contain metal-free NAD-dependent enzymes; ferroxidase 1.16.3.1
  is dominated by ferritins; these entries add few or no sites.
* Zn and Cr sites in the sweep are incidental (Zn/Cr were not sweep metals) and therefore incomplete.
* No functional filter: the sweep contains adventitious ions from crystallisation (e.g. Mn/Co/Ni with water and
  carboxylate shells), polyoxometalates and metallocene ligands (each Cp ring counts as 5 C donors).
* Very recent entries may not be in the RCSB cluster files yet and then form their own cluster (e.g. 30OH/30OR,
  two refinements of one rubredoxin data set).
* Fe-S cluster iron dominates site counts (frac_mononuclear shows cluster sites); geometric DV/SV mapping does not
  consider oxidation state, spin or electronic structure.
"""
    (DATA / "sweep_summary.md").write_text(text)
    print(text)


def show(entry_id):
    cols = ["metal", "chain", "metal_res", "seqid", "metal_atom", "cn", "formula", "graphene_motif", "dv_score",
            "tau4", "n_copies", "donors", "metal_neighbors", "weak_contacts"]
    for s in entry_sites(entry_id):
        print(" | ".join(str(s.get(c, "")) for c in cols))


# ---------------------------------------------------------------- complete d-block sweep (dataset release)

COMPLETE = DATA / "complete"
PARTS = COMPLETE / "sites_parts"
DBLOCK = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
          "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]  # = pdb_sweep.DBLOCK
BATCH_SIZE = 2000        # entries per part file
LARGE_ATOMS = 150_000    # entries with more deposited atoms run last, in batches of LARGE_BATCH with fewer workers
LARGE_BATCH = 100
LARGE_WORKERS = 4
ENTRY_TIMEOUT = 1800     # s; an entry still running after this is logged as failed (timeout)
WATCHDOG_EXTRA = 600     # s; no result for ENTRY_TIMEOUT + this: a worker died, the rest runs one by one
CLUSTER_LEVELS = (30, 50, 95, 100)

# column order and types of data/complete/sites.parquet (labels are joined in finish_complete)
COMPLETE_TYPES = {
    "site_id": "string", "entry_id": "string", "model": "Int64", "metal": "string", "chain": "string",
    "label_asym_id": "string", "metal_res": "string", "seqid": "string", "auth_seq_id": "Int64", "icode": "string",
    "metal_atom": "string", "altloc": "string", "metal_entity": "string", "x": "float64", "y": "float64",
    "z": "float64", "metal_occ": "float64", "metal_b": "float64", "metal_b_rel": "float64",
    "copy_group": "string", "is_representative": "boolean", "n_copies": "Int64", "chains": "string",
    "cn": "Int64", "formula": "string", "raw_motif": "string", "graphene_motif": "string",
    "preferred_site": "string", "dv_formula": "string", "dv_isomer": "string", "dv_score": "float64",
    "dv_tau4": "float64", "dv_classes": "string", "axial_formula": "string", "axial_classes": "string",
    "sv_formula": "string", "sv_score": "float64", "sv_face": "boolean", "sv_classes": "string",
    "sv_removed": "string", "tau4": "float64", "tau5": "float64", "mean_dist": "float64",
    "donor_classes": "string", "donors": "string", "n_metal_neighbors": "Int64", "metal_neighbors": "string",
    "weak_contacts": "string", "sym_donors": "Int64", "resolution": "float64", "host_entity": "string",
    "shell_json": "string"}
LABEL_TYPES = {"method": "string", "cluster30": "Int64", "cluster50": "Int64", "cluster95": "Int64",
               "cluster100": "Int64", "reaction": "string", "enzyme": "string", "query_ec": "string", "ec": "string",
               "description": "string", "organism": "string", "taxonomy_id": "string", "uniprot": "string",
               "label_source": "string"}
LOG_TYPES = {"entry_id": "string", "status": "string", "error": "string", "seconds": "float64",
             "n_sites": "Int64", "n_representative": "Int64", "file_revision": "string",
             "file_revision_date": "string", "n_models": "Int64", "crystal_symmetry": "boolean",
             "n_atoms_model1": "Int64", "n_metal_atoms_file": "Int64", "n_skipped_altloc": "Int64",
             "n_metal_atoms": "Int64", "n_skipped_zero_occupancy": "Int64", "n_skipped_minor_alternative": "Int64",
             "crystal_frame_standardized": "boolean", "part": "string", "worker_pid": "Int64",
             "batch_seconds": "float64", "batch_workers": "Int64"}


class EntryTimeout(Exception):
    pass


def _alarm(signum, frame):
    raise EntryTimeout(f"no result after {ENTRY_TIMEOUT} s")


def _init_worker():
    signal.signal(signal.SIGALRM, _alarm)


def group_copies(sites, entry_id):
    """merge_copies without dropping: every copy keeps its row; copy_group = same key as merge_copies, the copy
    with the lowest metal B (first one on ties, as merge_copies) is the representative."""
    groups = defaultdict(list)
    for s in sites:
        groups[s.pop("copy_key")].append(s)
    for k, copies in enumerate(groups.values(), 1):
        best = min(copies, key=lambda s: s["metal_b"])
        chains = " ".join(sorted({s["chain"] for s in copies}))
        for s in copies:
            s.update(copy_group=f"{entry_id}_{k}", n_copies=len(copies), chains=chains, is_representative=s is best)
    return sites


def assign_site_ids(sites, entry_id):
    """site_id = <entry>_<auth chain>_<resname>_<seqid><icode>_<atom>; where that is not unique within the entry
    the label_asym_id is appended (then a running number). Returns {(base id, label_asym_id): site_id}."""
    base = [f"{entry_id}_{s['chain']}_{s['metal_res']}_{s['seqid']}_{s['metal_atom']}" for s in sites]
    n = Counter(base)
    ids = [b if n[b] == 1 else f"{b}_{s['label_asym_id']}" for b, s in zip(base, sites)]
    n_ids, seen = Counter(ids), Counter()
    for s, i in zip(sites, ids):
        seen[i] += 1
        s["site_id"] = i if n_ids[i] == 1 else f"{i}_{seen[i]}"
    keys = Counter((b, s["label_asym_id"]) for b, s in zip(base, sites))
    return {(b, s["label_asym_id"]): s["site_id"] for b, s in zip(base, sites) if keys[(b, s["label_asym_id"])] == 1}


def complete_sites(entry_id):
    """All d-block sites of one entry (every copy) with site ids, copy groups and shell_json."""
    stats = {}
    sites = group_copies(raw_sites(entry_id, DBLOCK, detail=True, stats=stats), entry_id)
    by_atom = assign_site_ids(sites, entry_id)
    for s in sites:
        det = s.pop("_detail")
        for m in det["metals"]:  # neighbour metal -> its own site (asymmetric-unit copy), if it is one
            m["site_id"] = by_atom.get((f"{entry_id}_{m['auth_chain']}_{m['resname']}_{m['seqid']}{m['icode']}_"
                                        f"{m['atom']}", m["label_asym"]))
        s["shell_json"] = json.dumps(det, separators=(",", ":"))
    stats.update(n_sites=len(sites), n_representative=sum(s["is_representative"] for s in sites))
    return sites, stats


def complete_entry(entry_id):
    """Pool worker: (sites, log row) of one entry; never raises (errors and timeouts go to the log)."""
    t0 = time.perf_counter()
    log = {"entry_id": entry_id, "status": "ok", "error": "", "worker_pid": os.getpid()}
    try:
        signal.alarm(ENTRY_TIMEOUT)
        sites, stats = complete_sites(entry_id)
        log.update(stats)
    except Exception as exc:  # one broken file must not stop the sweep
        sites = []
        log.update(status="timeout" if isinstance(exc, EntryTimeout) else "failed",
                   error=f"{type(exc).__name__}: {exc}"[:2000])
    finally:
        signal.alarm(0)
    log["seconds"] = round(time.perf_counter() - t0, 3)
    return sites, log


def typed(df, types):
    df = df.reindex(columns=list(types))
    for c, t in types.items():
        if t == "boolean" and df[c].dtype == object:
            df[c] = df[c].map(lambda v: v if isinstance(v, (bool, np.bool_)) else pd.NA)
        df[c] = df[c].astype(t)
    return df


def isolated_entry(entry_id):
    """complete_entry in a fresh single-worker pool: a worker that dies (e.g. out of memory) is logged as failed."""
    with get_context("fork").Pool(1, initializer=_init_worker) as pool:
        try:
            return pool.apply_async(complete_entry, (entry_id,)).get(ENTRY_TIMEOUT + WATCHDOG_EXTRA)
        except mp.TimeoutError:
            return [], {"entry_id": entry_id, "status": "failed", "seconds": float("nan"),
                        "error": f"no result after {ENTRY_TIMEOUT + WATCHDOG_EXTRA} s: worker process died (e.g. "
                                 "killed for memory) or hung outside Python"}


def run_batch(entries, name, workers):
    """Process one batch, write sites_parts/<name>.parquet and then log_<name>.parquet (= batch done)."""
    t0, sites, logs = time.time(), [], []

    def collect(s, log):
        sites.extend(s)
        log["part"] = name
        logs.append(log)
        if len(logs) % 500 == 0 or log["seconds"] > 120:
            print(f"    {name}: {len(logs)}/{len(entries)} entries, {len(sites)} sites, {time.time() - t0:.0f} s"
                  + (f" (slow: {log['entry_id']} {log['seconds']:.0f} s)" if log["seconds"] > 120 else ""),
                  flush=True)

    try:
        with get_context("fork").Pool(workers, initializer=_init_worker, maxtasksperchild=200) as pool:
            results = pool.imap_unordered(complete_entry, entries)
            for _ in entries:  # a task whose worker died never returns: wait at most for one entry timeout
                collect(*results.next(timeout=ENTRY_TIMEOUT + WATCHDOG_EXTRA))
    except mp.TimeoutError:
        finished = {log["entry_id"] for log in logs}
        left = [e for e in entries if e not in finished]
        print(f"    {name}: no result for {ENTRY_TIMEOUT + WATCHDOG_EXTRA} s (a worker was lost); "
              f"{len(left)} entries one by one", flush=True)
        for e in left:
            collect(*isolated_entry(e))
    sites_df = typed(pd.DataFrame(sites), COMPLETE_TYPES)
    log_df = typed(pd.DataFrame(logs).assign(batch_seconds=round(time.time() - t0, 1), batch_workers=workers),
                   LOG_TYPES)
    for df, f in ((sites_df, PARTS / f"{name}.parquet"), (log_df, PARTS / f"log_{name}.parquet")):
        tmp = f.with_name(f.name + ".tmp")
        df.to_parquet(tmp, index=False, compression="zstd")
        tmp.rename(f)
    bad = (log_df.status != "ok").sum()
    print(f"  {name}: {len(entries)} entries, {len(sites_df)} sites, {bad} failed, {time.time() - t0:.0f} s, "
          f"slowest {log_df.seconds.max():.1f} s ({log_df.entry_id[log_df.seconds.idxmax()]})", flush=True)


def latest_log():
    """Entry log of all finished batches; an entry processed again (complete retry) keeps its last record."""
    files = sorted(PARTS.glob("log_part_*.parquet"))
    if not files:
        return typed(pd.DataFrame(), LOG_TYPES).assign(attempts=pd.Series(dtype="Int64"))
    log = typed(pd.concat([pd.read_parquet(f) for f in files], ignore_index=True), LOG_TYPES)
    attempts = log.entry_id.value_counts()
    log = log.sort_values("part", kind="stable").drop_duplicates("entry_id", keep="last", ignore_index=True)
    log["attempts"] = log.entry_id.map(attempts).astype("Int64")
    ok = log.status == "ok"  # batches written before this flag existed: no entry needed the standard frame
    log.loc[ok & log.crystal_frame_standardized.isna(), "crystal_frame_standardized"] = False
    return log


def run_complete(retry=False):
    """data/complete/entries.parquet (downloaded entries) -> sites_parts/ -> sites.parquet, motifs.parquet.
    Finished batches are skipped; retry: process the entries whose last attempt failed once more."""
    t0 = time.time()
    PARTS.mkdir(parents=True, exist_ok=True)
    entries = pd.read_parquet(COMPLETE / "entries.parquet")
    entries = entries[entries.downloaded.astype(bool)]
    log = latest_log()
    done = set(log.entry_id[log.status == "ok"] if retry else log.entry_id)
    todo = entries[~entries.entry_id.isin(done)]
    large = todo.deposited_atom_count.fillna(0) > LARGE_ATOMS
    normal = sorted(todo.entry_id[~large])
    big = list(todo[large].sort_values("deposited_atom_count").entry_id)
    n_done = log.part.nunique()
    batches = [(normal[i:i + BATCH_SIZE], N_WORKERS) for i in range(0, len(normal), BATCH_SIZE)]
    batches += [(big[i:i + LARGE_BATCH], LARGE_WORKERS) for i in range(0, len(big), LARGE_BATCH)]
    print(f"{len(entries)} downloaded entries: {len(done)} done in {n_done} parts, {len(normal)} to do in "
          f"{-(-len(normal) // BATCH_SIZE)} batches, {len(big)} large (> {LARGE_ATOMS} atoms) last", flush=True)
    numbers = [int(f.stem.split("_")[-1]) for f in PARTS.glob("log_part_*.parquet")]
    for k, (batch, workers) in enumerate(batches, max(numbers, default=-1) + 1):
        run_batch(batch, f"part_{k:04d}", workers)
    finish_complete(time.time() - t0)


def label_complete(sites):
    """Host-entity labels first (label_source = entity), else the union over the entry (entry), like label_sites;
    plus UniProt, taxonomy and sequence clusters of the host entity and the entry's experimental method."""
    ents = pd.read_parquet(COMPLETE / "entities.parquet")
    keys = ["reaction", "enzyme", "query_ec"]
    host = ents.rename(columns={"entity_id": "host_entity"})[
        ["entry_id", "host_entity", "description", "ec", "organism", "taxonomy_id", "uniprot", *keys]
        + [f"cluster{n}" for n in CLUSTER_LEVELS]]
    df = sites.merge(host, on=["entry_id", "host_entity"], how="left")
    entry_lab = (ents[ents.reaction != ""].groupby("entry_id")[keys]
                 .agg(lambda col: ";".join(sorted({x for v in col if v for x in v.split(";")}))))
    own = df.reaction.fillna("") != ""
    use_entry = ~own & df.entry_id.isin(entry_lab.index)
    for k in keys:
        df[k] = df[k].fillna("").astype(object)
        df.loc[use_entry, k] = df.loc[use_entry, "entry_id"].map(entry_lab[k])
    df["label_source"] = np.where(own, "entity", np.where(use_entry, "entry", ""))
    for c in ("description", "ec", "organism", "taxonomy_id", "uniprot"):
        df[c] = df[c].fillna("")
    meth = pd.read_parquet(COMPLETE / "entries.parquet", columns=["entry_id", "method"])
    df = df.merge(meth, on="entry_id", how="left")
    return typed(df, {**COMPLETE_TYPES, **LABEL_TYPES})


def motif_table_complete(df):
    """Per graphene motif and per raw first-shell motif: sites, entries and 30 % clusters, over all copies and
    over representative copies only; descriptive columns from the representative copies."""
    clu = df.cluster30.astype("string").fillna(df.entry_id + "_" + df.host_entity.fillna(""))
    df = df.assign(clu30=clu)
    out = []
    for kind, col in [("graphene", "graphene_motif"), ("first_shell", "raw_motif")]:
        rep_groups = dict(list(df[df.is_representative].groupby(col)))
        for motif, g in df.groupby(col):
            r = rep_groups.get(motif, g.iloc[:0])
            d = r if len(r) else g  # a motif seen only on non-representative copies: describe those
            best = d.sort_values(["resolution", "entry_id"]).drop_duplicates("entry_id")
            out.append({
                "kind": kind, "motif": motif, "metal": g.metal.iloc[0],
                "n_sites": len(g), "n_entries": g.entry_id.nunique(), "n_clusters30": g.clu30.nunique(),
                "n_sites_rep": len(r), "n_entries_rep": r.entry_id.nunique(), "n_clusters30_rep": r.clu30.nunique(),
                "n_xray_entries": g.entry_id[g.method == "X-ray"].nunique(),
                "n_em_entries": g.entry_id[g.method == "EM"].nunique(),
                "reactions": ";".join(sorted({x for v in g.reaction.fillna("") if v for x in v.split(";")})),
                "example_descriptions": " | ".join(dict.fromkeys(x for x in best.description.fillna("") if x))[:200],
                "example_pdb": " ".join(f"{e}({x:.2f})" if pd.notna(x) else e
                                        for e, x in zip(best.entry_id[:6], best.resolution[:6])),
                "best_resolution": d.resolution.min(),
                "median_tau4": d.tau4.median(), "median_dv_tau4": d.dv_tau4.median(),
                "median_dv_score": d.dv_score.median(), "median_sv_score": d.sv_score.median(),
                "frac_mononuclear": round(float((d.n_metal_neighbors == 0).mean()), 3),
                "donor_classes_top": d.donor_classes.fillna("").value_counts().index[0]})
    m = pd.DataFrame(out)
    for c in ("median_tau4", "median_dv_tau4", "median_dv_score", "median_sv_score"):
        m[c] = m[c].astype("float64").round(3)
    return m.sort_values(["kind", "n_entries", "n_clusters30", "n_sites"], ascending=[False, False, False, False],
                         ignore_index=True)


def sha256(path):
    import hashlib
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def finish_complete(seconds):
    """Concatenate the parts -> sites.parquet (+ labels), motifs.parquet, entry_log.parquet, failures.csv,
    provenance."""
    t0 = time.time()
    log = latest_log()
    frames = []
    for part, ids in log.groupby("part").entry_id:  # sites of each entry from the batch of its last attempt
        df = pd.read_parquet(PARTS / f"{part}.parquet")
        frames.append(df[df.entry_id.isin(set(ids))])
    sites = pd.concat(frames, ignore_index=True)
    sites = sites.sort_values("entry_id", kind="stable", ignore_index=True)  # keeps the model order per entry
    sites = label_complete(sites)
    sites.to_parquet(COMPLETE / "sites.parquet", index=False, compression="zstd")
    log = log.sort_values("entry_id", ignore_index=True)
    log.to_parquet(COMPLETE / "entry_log.parquet", index=False, compression="zstd")
    fail = log[log.status != "ok"]
    entries = pd.read_parquet(COMPLETE / "entries.parquet")
    lost = entries[~entries.downloaded.astype(bool)]
    pd.concat([pd.DataFrame({"entry_id": lost.entry_id, "stage": "download", "status": "failed",
                             "error": lost.download_status, "seconds": np.nan}),
               fail[["entry_id", "status", "error", "seconds"]].assign(stage="extraction")],
              ignore_index=True)[["entry_id", "stage", "status", "error", "seconds"]].to_csv(
        COMPLETE / "failures.csv", index=False)
    motifs = motif_table_complete(sites)
    motifs.to_parquet(COMPLETE / "motifs.parquet", index=False)
    rep = sites[sites.is_representative]
    from pdb_sweep import code_version, software_versions, update_provenance, utcnow
    update_provenance("extraction", {
        "finished_utc": utcnow(), "seconds_last_run": round(seconds + time.time() - t0),
        "entries_found": len(entries), "entries_downloaded": int(entries.downloaded.sum()),
        "entries_processed": int((log.status == "ok").sum()), "entries_failed": int((log.status != "ok").sum()),
        "entries_retried": log.entry_id[log.attempts > 1].tolist(),
        "entries_crystal_frame_standardized": log.entry_id[log.crystal_frame_standardized.fillna(False)].tolist(),
        "entries_with_sites": int(sites.entry_id.nunique()),
        "entries_processed_without_site": int(((log.status == "ok") & (log.n_sites == 0)).sum()),
        "sites": len(sites), "representative_sites": len(rep),
        "metal_atoms_model1": {"in_file_all_conformers": int(log.n_metal_atoms_file.sum()),
                               "not_selected_conformer": int(log.n_skipped_altloc.sum()),
                               "zero_occupancy": int(log.n_skipped_zero_occupancy.sum()),
                               "minor_alternative_position": int(log.n_skipped_minor_alternative.sum()),
                               "sites": len(sites)},
        "entry_seconds": {"sum": round(float(log.seconds.sum())), "median": float(log.seconds.median()),
                          "max": float(log.seconds.max()), "max_entry": log.entry_id[log.seconds.idxmax()]},
        "batches": int(log.part.nunique()),
        "batch_wall_seconds": round(float(log.drop_duplicates("part").batch_seconds.sum())),
        "metals": DBLOCK, "model": "first model only (model number in column model)",
        "thresholds": {"donor_cutoff_A": CUTOFF, "metal_metal_cutoff_A": METAL_CUTOFF,
                       "weak_extra_A": WEAK_EXTRA, "min_dist_A": MIN_DIST, "metal_alternative_A": METAL_ALT,
                       "minor_occupancy_fraction": MINOR_OCC, "overlap_angle_deg": OVERLAP_ANGLE,
                       "dv_score_max_deg": DV_SCORE_MAX, "dv_tie_deg": DV_TIE, "face_max_angle_deg": FACE_MAX_ANGLE,
                       "isomer_trans_deg": ISOMER_TRANS, "neighbour_search_radius_A": max(CUTOFF.values()) + WEAK_EXTRA,
                       "large_entry_atoms": LARGE_ATOMS, "entry_timeout_s": ENTRY_TIMEOUT},
        "copies": "every copy is a row; copy_group = same metal element, metal residue name and atom name, same "
                  "donor classes/residues/atoms (residue numbers of polymer donors only), and without polymer donors "
                  "the same nearest polymer residue (number, name); is_representative = lowest metal B in the group "
                  "(first in model order on ties), as merge_copies of the older modes",
        "shell_json": {
            "donors": "first-shell donors (sorted by distance): el, cls (donor class), auth_chain, label_asym, "
                      "resname, seqid, icode, atom, dist, dx/dy/dz (A, contact minus metal, symmetry applied), occ, "
                      "b, sym, polymer, entity",
            "weak": "contacts that are not donors, same fields + reason ('+'-joined): beyond_cutoff (cutoff..cutoff + "
                    f"{WEAK_EXTRA} A), minor_occupancy (occ < {MINOR_OCC} x metal occ), alt_disorder (half-occupied "
                    "symmetry mate of a ligand disordered across a symmetry element), alt_overlap (L-M-L < "
                    f"{OVERLAP_ANGLE} deg: second copy of one ligand), protein_carbon (protein C is never a donor), "
                    "his_second_ring_n (the farther of ND1/NE2)",
            "other": f"atoms of elements without a donor cutoff (e.g. P, F, Br, I, As, Te) within "
                     f"{max(CUTOFF.values()) + WEAK_EXTRA} A, reason element_not_counted; never donors",
            "metals": f"other metal atoms (any metal, incl. alkali/alkaline earth) at {METAL_ALT}-{METAL_CUTOFF} A, "
                      "not minority states; site_id of that atom when it is a site of this table (symmetry mate: "
                      "the asymmetric-unit copy), else null",
            "sym": "0 = atom of the deposited model; k > 0 = gemmi cell image k (crystal symmetry operators, plus NCS "
                   "operators that the file lists as not applied); EM/NMR entries are searched without symmetry"},
        "crystal_frame": "entries with a non-standard crystal frame (SCALEn) are moved to the standard frame for the "
                         "neighbour search; x/y/z and dx/dy/dz are reported in the deposited frame",
        "software": software_versions(), "code": code_version()})
    update_provenance("files", {f.name: {"bytes": f.stat().st_size, "sha256": sha256(f)} for f in
                                [COMPLETE / n for n in ("entries.parquet", "entities.parquet", "sites.parquet",
                                                        "motifs.parquet", "entry_log.parquet", "failures.csv")]})
    print(f"{len(sites)} sites ({len(rep)} representative) in {sites.entry_id.nunique()} entries; "
          f"{len(fail)} entries failed -> data/complete/sites.parquet, motifs.parquet, failures.csv", flush=True)
    print(rep.metal.value_counts().to_string(), flush=True)


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "enzymes":
        run_set("hits_enzymes.csv", "enzyme_sites.csv", "motifs_enzyme.csv")
    elif mode == "all":
        run_set("hits_all.csv", "pdb_sites_all.csv", "motifs_pdb.csv")
    elif mode == "check":
        run_checks()
    elif mode == "summary":
        write_summary()
    elif mode == "show":
        for e in sys.argv[2:]:
            show(e.upper())
    elif mode == "complete":
        if sys.argv[2:] == ["finish"]:
            finish_complete(0)
        else:
            run_complete(retry=sys.argv[2:] == ["retry"])
    else:
        sys.exit(__doc__)
