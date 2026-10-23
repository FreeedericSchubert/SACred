"""What follows the xTB library (screen/xtb_scan.py library); runs on Curta as soon as the library jobs end.

    python screen/library_next.py     # -> data/xtb_reactions_library_todo.json, data/xtb_selected_library.csv

Rebuilds data/xtb_scan.csv with the library sites and calibrates it as screen/calibrate.py does: one offset per
metal for H on the metal and one for H on S, fitted on every xTB/DFT pair so far; metals with fewer than two DFT pairs
(V, Cr, Ru, Pt, ...) keep the raw xTB value. A site counts with its H position closest to 0, if the H stayed where it was
put, both relaxations converged and the metal stayed in its site (bare metal-donor distances <= 2.6 A).
- O2, CO2 and N2 intermediates (screen/xtb_reactions.py shard ... data/xtb_reactions_library_todo.json): the library
  sites with |dG_H| <= WINDOW, closest to 0 first, at most MAX_REACTION_SITES.
- DFT (dft/make_inputs.py library): the best N_PROOF sites of each proof-of-concept metal (Cr, Ru, Pt: raw xTB, so
  DFT is the check) and the best calibrated site of N_ENZYME enzyme metals; at most one C donor, as calibrate.py.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(os.path.dirname(HERE), "data")
sys.path.insert(0, HERE)
import xtb_reactions                     # noqa: E402
import xtb_scan                          # noqa: E402
from calibrate import canonical, offsets  # noqa: E402

WINDOW = 0.3                 # eV, as screen/calibrate.py
MAX_REACTION_SITES = 100     # x 11 structures each
N_PROOF, N_ENZYME = 2, 4
TODO = os.path.join(DATA, "xtb_reactions_library_todo.json")
SELECTED = os.path.join(DATA, "xtb_selected_library.csv")


def calibrate(x):
    """dG_H_cal as in screen/calibrate.py, with the DFT values of data/dgh_graphene.csv."""
    d = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv")).dropna(subset=["dG_H"])
    x["key"] = x.motif.map(canonical) + "#" + x.state
    d["key"] = d.motif.map(canonical) + "#" + d.state
    pairs = d.merge(x[x.ok].drop_duplicates("key"), on="key", suffixes=("_dft", "_xtb"))
    pairs["metal"] = pairs.motif_dft.str.split("-").str[0]
    o = offsets(pairs)                       # as screen/calibrate.py: a metal needs two DFT pairs
    o = o[o.used]
    off = o[o.state == "H_M"].set_index("metal").offset
    off_s = o[o.state == "H_S"].offset.iloc[0]
    x["dG_H_cal"] = np.where(x.state == "H_S", x.dG_H + off_s, x.dG_H + x.metal.map(off))
    print("offsets from %d pairs: %s; H on S %+.2f eV" % (len(pairs), ", ".join("%s %+.2f" % kv for kv in off.items()),
                                                           off_s))
    return x, off


def main():
    x = xtb_scan.table()
    library = {j[0] for j in json.load(open(xtb_scan.LIBRARY_JOBS)) if j[0] != "H2"}
    x["metal"] = x.motif.str.split("-").str[0]
    x["ok"] = (x.converged == True) & (x.bare_converged == True) & x.dG_H.notna()   # noqa: E712 (NaN-safe)
    x, off = calibrate(x)
    intended = np.where(x.state == "H_M", x.H_partner == "M", x.H_partner.astype(str).str.endswith("(donor)"))
    good = x.ok & intended & (x.max_M_X_bare <= 2.6) & x.motif.isin(library)
    x["dG"] = x.dG_H_cal.fillna(x.dG_H)
    best = x[good].assign(a=lambda t: t.dG.abs()).sort_values("a").drop_duplicates("motif")

    rx = best[best.a <= WINDOW].head(MAX_REACTION_SITES)
    jobs = [("graphene", m, s) for m in rx.motif for s in ["bare"] + xtb_reactions.ADSORBATES]
    json.dump(jobs, open(TODO, "w"))
    print("%d library sites within %.1f eV (%d usable of %d); %d go to the O2/CO2/N2 screen -> %d structures" % (
        (best.a <= WINDOW).sum(), WINDOW, len(best), len(library), len(rx), len(jobs)))

    # DFT as in screen/calibrate.py: at most one C donor (all-carbon vacancies are neither enzyme-like nor new),
    # and one site per enzyme metal so the check covers several metals
    n_c = best.motif.str.split("-").str[1].str.split(":").str[0].str.extract(r"C(\d*)(?![a-z])")[0]
    dft = best[n_c.map(lambda v: 0 if pd.isna(v) else int(v or 1)) < 2]
    proof = dft[dft.metal.isin(xtb_scan.PROOF_METALS)].groupby("metal").head(N_PROOF)
    enzyme = dft[dft.metal.isin(off.index) & dft.dG_H_cal.notna() & (dft.a <= WINDOW)]
    enzyme = enzyme.drop_duplicates("metal").head(N_ENZYME)
    sel = pd.concat([proof, enzyme])[["motif", "state", "dG_H", "dG_H_cal"]]
    sel.to_csv(SELECTED, index=False, float_format="%.3f")
    print("DFT for %d sites -> %s\n%s" % (len(sel), SELECTED, sel.to_string(index=False, float_format="%.2f")))


if __name__ == "__main__":
    main()
