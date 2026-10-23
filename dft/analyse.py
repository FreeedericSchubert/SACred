"""Hydrogen adsorption free energies of the graphene single-atom sites from the pw.x outputs.

    bash dft/sync.sh pull && python dft/analyse.py      # -> data/dgh_graphene.csv, data/dgh_mos2.csv, ranked tables

dG_H = E(site+H) - E(site) - E(H2)/2 + 0.24 eV, all energies from the 3x3x1 SCF at the relaxed geometries.
0.24 eV is the zero-point and entropy correction of Norskov et al., J. Electrochem. Soc. 152, J23 (2005).
Per H structure the table also says where the H ended up (its nearest atom), so an H that left the metal
or the S during the relaxation is not mistaken for the intended site, and whether the metal stayed in the
vacancy (largest M-donor distance).

Stability of the metal in the site (data/stability_graphene.csv): E_f = E(site with M) - E(site without M)
- E(bulk M per atom), < 0 if the single atom is more stable in the site than in its bulk metal (no driving
force to cluster); U_diss = U0(M -> M^n+) - E_f / n (Greeley & Norskov, Electrochim. Acta 52, 5829 (2007)),
the potential above which the metal dissolves; > 0 V vs SHE is the usual criterion for stability under HER.
"""
import glob
import json
import os
import re

import numpy as np
import pandas as pd
from ase.io import read

RY = 13.605693123     # eV
DG_CORR = 0.24        # eV
HERE = os.path.dirname(os.path.abspath(__file__))
RUNS = os.path.join(HERE, "runs", "graphene")
OUT = os.path.join(HERE, "..", "data", "dgh_graphene.csv")
OUT_MOS2 = os.path.join(HERE, "..", "data", "dgh_mos2.csv")
OUT_STAB = os.path.join(HERE, "..", "data", "stability_graphene.csv")
MOS2 = os.path.join(HERE, "runs", "mos2")
BULK = os.path.join(HERE, "runs", "bulk")
BOHR = 0.529177       # A
ZVAL_MOS2 = {"Mo": 14, "S": 6, "C": 4, "H": 1, "N": 5, "Fe": 16}   # valence electrons, MoS2 pseudopotential set
# standard potentials (V vs SHE) and electrons of M^n+ + n e- -> M (CRC Handbook) for the dissolution potential
DISSOLUTION = {"V": (-1.175, 2), "Mn": (-1.185, 2), "Fe": (-0.447, 2), "Co": (-0.28, 2), "Ni": (-0.257, 2),
               "Cu": (0.342, 2), "Zn": (-0.762, 2), "Mo": (-0.200, 3), "W": (0.1, 3), "Cr": (-0.744, 3),
               "Ru": (0.455, 2), "Pt": (1.18, 2)}


def scf_energy(path):
    """(energy in eV, total magnetization) of a converged scf.out, else (nan, nan)."""
    if not os.path.exists(path):
        return np.nan, np.nan
    text = open(path).read()
    if "convergence has been achieved" not in text or "JOB DONE" not in text:
        return np.nan, np.nan
    e = float(re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)", text)[-1]) * RY
    mag = re.findall(r"total magnetization\s+=\s+(-?\d+\.\d+)", text)
    return e, float(mag[-1]) if mag else 0.0


def relax_energy(path):
    """Final energy (eV) of a finished Gamma-point relaxation (the MoS2 runs take energies from it), else nan."""
    if not os.path.exists(path):
        return np.nan, np.nan
    text = open(path).read()
    if "End final coordinates" not in text:
        return np.nan, np.nan
    mag = re.findall(r"total magnetization\s+=\s+(-?\d+\.\d+)", text)
    return float(re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)", text)[-1]) * RY, float(mag[-1]) if mag else 0.0


def gate_energy(path):
    """(energy in eV, magnetization) of a finished gate relaxation without QE 7.1's plate-constant artefact, as in
    the earlier MoS2 project (QE_MOS2 in dft/curta/mos2_h.py): E - 2 N_val gateamp L / 6 Ry, N_val = valence
    electrons of the neutral cell, gateamp = 'prefactor of the potential' [Ha a.u.], L = c in bohr."""
    e, mag = relax_energy(path)
    if not np.isfinite(e):
        return e, mag
    gateamp = float(re.findall(r"prefactor of the potential in \[Ha a\.u\.\]:\s+(-?\d+\.\d+)", open(path).read())[-1])
    atoms = read(path, format="espresso-out", index=-1)
    n_val = sum(ZVAL_MOS2[el] for el in atoms.get_chemical_symbols())
    return e - 2 * n_val * gateamp * atoms.cell[2][2] / BOHR / 6 * RY, mag


def mos2():
    """H on MoS2: on a basal S of clean MoS2 (host from the earlier MoS2 project), and on an Fe-S3
    site (host computed here, the lower of the top-Mo and hollow sites), also with the sheet charged against a
    plate (charge_e extra electrons per cell); see mos2_refs.json."""
    refs = json.load(open(os.path.join(HERE, "mos2_refs.json")))
    e_h2, _ = relax_energy(os.path.join(MOS2, refs["h2"], "relax.out"))
    rows = []
    hosts = {run: (ref["host"], ref["E_Ry"] * RY) for run, ref in refs["hosts"].items()}
    for run, candidates in refs.get("computed", {}).items():     # e.g. Fe on MoS2: lowest finished bare site
        done = {c: relax_energy(os.path.join(MOS2, c, "relax.out"))[0] for c in candidates}
        done = {c: e for c, e in done.items() if np.isfinite(e)}
        hosts[run] = min(done.items(), key=lambda kv: kv[1]) if done else (" or ".join(candidates), np.nan)
        if len(done) > 1:
            print("MoS2 %s: %s" % (run, ", ".join("%s %+.3f eV" % (c, e - min(done.values())) for c, e in done.items())))
    for run, (host, e_host) in hosts.items():
        e, mag = relax_energy(os.path.join(MOS2, run, "relax.out"))
        de = e - e_host - e_h2 / 2
        rows.append({"run": run, "host": host, "E_eV": e, "magnetization": mag, "dE_H": de, "dG_H": de + DG_CORR})
    for run, host in refs.get("gate", {}).items():          # charged against a plate: host at the same charge
        qdirs = sorted(os.listdir(os.path.join(MOS2, run))) if os.path.isdir(os.path.join(MOS2, run)) else []
        for q in qdirs:
            e, mag = gate_energy(os.path.join(MOS2, run, q, "relax.out"))
            e_host, _ = gate_energy(os.path.join(MOS2, host, q, "relax.out"))
            de = e - e_host - e_h2 / 2
            rows.append({"run": "%s/%s" % (run, q), "host": "%s/%s" % (host, q), "charge_e": float(q[2:]), "E_eV": e,
                         "magnetization": mag, "dE_H": de, "dG_H": de + DG_CORR})
    df = pd.DataFrame(rows)
    df.to_csv(OUT_MOS2, index=False, float_format="%.4f")
    print("\nMoS2 (H on a basal S, or on Fe):")
    print(df[["run", "dG_H", "magnetization"]].to_string(index=False, float_format="%.2f", na_rep="not finished"))


def bulk_energy(metal):
    """Energy per atom (eV) of the relaxed bulk metal (vc-relax), else nan. A vc-relax that ran into nstep while its
    energy no longer changed (Cu: the cell bounces between two volumes 0.2 % apart, stress noise at 40 Ry) counts:
    its lowest energy of the last 10 steps, if they span less than 1e-4 Ry."""
    out = os.path.join(BULK, metal, "relax.out")
    e, _ = relax_energy(out)
    if not np.isfinite(e) and os.path.exists(out):
        text = open(out).read()
        last = [float(x) for x in re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)", text)][-10:]
        if "maximum number of steps has been reached" in text and len(last) == 10 and max(last) - min(last) < 1e-4:
            e = min(last) * RY
    if not np.isfinite(e):
        return np.nan
    return e / int(re.search(r"number of atoms/cell\s+=\s+(\d+)", open(out).read()).group(1))


def stability(df):
    """Metal formation energy E_f and dissolution potential U_diss per site (bare reference as for dG_H)."""
    st = df.drop_duplicates("motif")[["motif", "bare_ref", "E_bare", "mag_bare"]].copy()
    st["metal"] = st.motif.str.split("-").str[0]
    env = st.motif.str.split("-", n=1).str[1]
    st["site_ref"] = "site_" + env.str.replace(":", "_").str.replace("=", "")
    site = {r: scf_energy(os.path.join(RUNS, "_ref", r, "scf.out")) for r in st.site_ref.unique()}
    st["E_site"] = st.site_ref.map(lambda r: site[r][0])
    st["mag_site"] = st.site_ref.map(lambda r: site[r][1])
    st["E_bulk_atom"] = st.metal.map({m: bulk_energy(m) for m in st.metal.unique()})
    st["E_f"] = st.E_bare - st.E_site - st.E_bulk_atom
    st["U_diss"] = st.metal.map(lambda m: DISSOLUTION.get(m, (np.nan,))[0]) - st.E_f / st.metal.map(
        lambda m: DISSOLUTION.get(m, (np.nan, np.nan))[1])
    st.to_csv(OUT_STAB, index=False, float_format="%.4f")
    return st.set_index("motif")[["E_f", "U_diss"]]


def final_atoms(run):
    out = os.path.join(run, "relax.out")
    if not os.path.exists(out) or "End final coordinates" not in open(out).read():
        return None
    return read(out, format="espresso-out", index=-1)


def describe(atoms, meta):
    """Largest metal-donor distance, metal height over the sheet, and the H partner atom."""
    im, idon = meta["imetal"], meta["idonors"]
    d = atoms.get_distances(im, idon, mic=True)
    carbons = [a.index for a in atoms if a.symbol == "C"]
    height = atoms.positions[im, 2] - np.median(atoms.positions[carbons, 2])
    row = {"max_M_X": d.max(), "M_height": height}
    if "ih" in meta:
        dh = atoms.get_distances(meta["ih"], range(len(atoms)), mic=True)
        dh[meta["ih"]] = np.inf
        j = int(np.argmin(dh))
        row["H_partner"] = "M" if j == im else atoms[j].symbol + ("(donor)" if j in idon else "")
        row["H_dist"] = dh[j]
    return row


def main():
    e_h2, _ = scf_energy(os.path.join(RUNS, "_ref", "H2", "scf.out"))
    rows = []
    for motif_dir in sorted(glob.glob(os.path.join(RUNS, "*"))):
        if os.path.basename(motif_dir).startswith("_"):
            continue
        # reference: the lowest bare configuration found (e.g. bare_Sdown, started from a relaxed H geometry,
        # when a donor moved during H adsorption); bare_ref says which one
        bares = sorted(glob.glob(os.path.join(motif_dir, "bare*")))
        energies = {b: scf_energy(os.path.join(b, "scf.out")) for b in bares}
        done = [b for b in bares if np.isfinite(energies[b][0])]
        bare = min(done, key=lambda b: energies[b][0]) if done else os.path.join(motif_dir, "bare")
        meta = json.load(open(os.path.join(motif_dir, "bare", "meta.json")))
        e0, m0 = energies.get(bare, (np.nan, np.nan))
        at0 = final_atoms(bare)
        base = {"motif": meta["motif"], "enzyme": meta["enzyme"], "bare_ref": os.path.basename(bare),
                "E_bare": e0, "mag_bare": m0,
                "bare_states_pending": " ".join(os.path.basename(b) for b in bares if b not in done)}
        if at0 is not None:
            base.update({k + "_bare": v for k, v in describe(at0, meta).items()})
        for hrun in sorted(glob.glob(os.path.join(motif_dir, "H_*"))):
            hmeta = json.load(open(os.path.join(hrun, "meta.json")))
            eh, mh = scf_energy(os.path.join(hrun, "scf.out"))
            row = dict(base, state=hmeta["state"], E_H=eh, mag_H=mh)
            ath = final_atoms(hrun)
            if ath is not None:
                row.update(describe(ath, hmeta))
            row["dE_H"] = eh - e0 - e_h2 / 2
            row["dG_H"] = row["dE_H"] + DG_CORR
            rows.append(row)
    df = pd.DataFrame(rows)
    df = df.join(stability(df), on="motif")
    df.to_csv(OUT, index=False, float_format="%.4f")
    done = df.dropna(subset=["dG_H"]).copy()
    print("H2 reference: %s" % ("%.4f eV" % e_h2 if np.isfinite(e_h2) else "not finished"))
    print("%d of %d H structures finished; metal formation energy for %d of %d sites\n" % (
        len(done), len(df), df.drop_duplicates("motif").E_f.notna().sum(), df.motif.nunique()))
    if len(done):
        done["absdG"] = done["dG_H"].abs()
        cols = ["motif", "state", "dG_H", "E_f", "U_diss", "bare_ref", "H_partner", "H_dist", "mag_bare", "mag_H",
                "max_M_X", "bare_states_pending"]
        cols = [c for c in cols if c in done]
        print(done.sort_values("absdG")[cols].to_string(index=False, float_format="%.2f"))
    mos2()


if __name__ == "__main__":
    main()
