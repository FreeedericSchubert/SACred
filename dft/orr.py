"""O2 cross-check: DFT (pw.x) O2-reduction intermediates on two Co sites against GFN2-xTB and Cho et al. 2022.

    bash dft/sync.sh pull && python dft/orr.py      # -> data/dft_orr.csv, data/dft_orr.md

Runs: dft/runs/graphene/<motif>/ORR_{OH,O,OOH} (dft/make_inputs.py orr) with the bare site and the H2 and H2O
references of the same cell, all energies from the 3x3x1 SCF. dE on the Catalysis-Hub scale (vs H2O and H2, the
numbers Cho et al. publish), dG = dE + the ZPE - TS terms of screen/xtb_reactions.py, limiting potentials with its
limiting(), so the DFT, xTB and Cho numbers go through the same free-energy steps. Per adsorbate the table also
says whether it stayed on the metal (M-O distance).
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path[:0] = [HERE, os.path.join(ROOT, "screen")]
from analyse import RUNS, scf_energy, final_atoms     # noqa: E402
from make_inputs import ORR_SITES, ORR_STATES, dirname   # noqa: E402
import xtb_reactions as xr                           # noqa: E402

DATA = os.path.join(ROOT, "data")
BOUND = 2.3           # A: an O farther than this from the metal has left it


def gas_shift(X):
    """refG - refE of each O intermediate from the gas ZPE - TS terms."""
    z = {m: xr.ZPE_TS["gas_" + m] for m in ("H2", "H2O")}
    return {a: v for a, v in xr.references(dict(z, CO2=0, N2=0)).items() if a in ORR_STATES}


def main():
    e_h2, _ = scf_energy(os.path.join(RUNS, "_ref", "H2", "scf.out"))
    e_h2o, _ = scf_energy(os.path.join(RUNS, "_ref", "H2O", "scf.out"))
    ref = {a: v for a, v in xr.references({"H2": e_h2, "H2O": e_h2o, "CO2": 0, "N2": 0}).items() if a in ORR_STATES}
    shift = gas_shift(None)
    ext = pd.read_csv(os.path.join(DATA, "external", "catalysis_hub_cho2022_orr.csv"), comment="#")
    ext["site"] = ext.metal + "-" + ext.donors + ":DV"
    cho = ext.set_index(["site", "adsorbate"]).dE_eV
    xtb = pd.read_csv(os.path.join(DATA, "xtb_reactions.csv")).set_index("site")
    rows, sites = [], []
    for motif in ORR_SITES:
        d = os.path.join(RUNS, dirname(motif))
        e0, _ = scf_energy(os.path.join(d, "bare", "scf.out"))
        lim = {"dft": {}, "xtb": {}, "cho": {}}
        for a in ORR_STATES:
            run = os.path.join(d, "ORR_" + a)
            meta = json.load(open(os.path.join(run, "meta.json")))
            e, mag = scf_energy(os.path.join(run, "scf.out"))
            atoms = final_atoms(run)
            d_mo = (float(atoms.get_distance(meta["imetal"], meta["iads"][0], mic=True))
                    if atoms is not None else np.nan)
            de = e - e0 - ref[a]
            row = {"site": motif, "adsorbate": a, "dE_dft": de, "dG_dft": de + xr.ZPE_TS[a] - shift[a],
                   "M_O": d_mo, "on_metal": bool(d_mo <= BOUND) if np.isfinite(d_mo) else None, "mag": mag,
                   "dE_xtb": xtb.loc[motif, "dE_" + a] if motif in xtb.index else np.nan,
                   "dG_xtb": xtb.loc[motif, "dG_" + a] if motif in xtb.index else np.nan,
                   "dE_cho": cho.get((motif, a), np.nan)}
            row["dG_cho"] = row["dE_cho"] + xr.ZPE_TS[a] - shift[a]
            rows.append(row)
            for k in lim:
                lim[k]["dG_" + a] = row["dG_" + k]
        s = {"site": motif}
        for k, v in lim.items():
            if all(pd.notna(x) for x in v.values()):
                s["U_ORR_" + k] = xr.limiting(v)["U_ORR"]
        sites.append(s)
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA, "dft_orr.csv"), index=False, float_format="%.3f")
    st = pd.DataFrame(sites)
    st.to_csv(os.path.join(DATA, "dft_orr_sites.csv"), index=False, float_format="%.3f")

    def f(x, fmt="%+.2f"):
        return fmt % x if pd.notna(x) else "-"

    lines = ["# O2 cross-check: DFT against GFN2-xTB and Cho et al. 2022", "",
             "H2O reference: %s; H2 reference: %s. %d of %d intermediates finished." % (
                 "done" if np.isfinite(e_h2o) else "not finished", "done" if np.isfinite(e_h2) else "not finished",
                 df.dE_dft.notna().sum(), len(df)), "",
             "| site | adsorbate | dE DFT (eV) | dE xTB | dE Cho | M-O (A) | on metal |", "|---|---|---|---|---|---|---|"]
    lines += ["| %s | %s | %s | %s | %s | %s | %s |" % (r.site, r.adsorbate, f(r.dE_dft), f(r.dE_xtb), f(r.dE_cho),
                                                       f(r.M_O, "%.2f"), "-" if r.on_metal is None else
                                                       ("yes" if r.on_metal else "no")) for r in df.itertuples()]
    lines += ["", "| site | U_ORR DFT (V) | U_ORR xTB | U_ORR Cho |", "|---|---|---|---|"]
    lines += ["| %s | %s | %s | %s |" % (r["site"], f(r.get("U_ORR_dft")), f(r.get("U_ORR_xtb")),
                                         f(r.get("U_ORR_cho"))) for r in sites]
    open(os.path.join(DATA, "dft_orr.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
