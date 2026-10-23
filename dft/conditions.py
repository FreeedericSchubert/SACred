"""Does Co-N2OC keep its metal under H2-evolution conditions? DFT check of the two caveats from the literature.

    bash dft/sync.sh pull && python dft/conditions.py      # -> data/dft_conditions.csv, data/dft_conditions.md

Runs (dft/make_inputs.py conditions), energies from the 3x3x1 SCF like dft/analyse.py:
- H on the donors of Co-N2OC (Hdonor_C, Hdonor_O): dG = E(site+H) - E(site) - E(H2)/2 + 0.24 eV, as dG_H.
  A C donor that takes up H near 0 V vs RHE is hydrogenated under H2 evolution (Li et al., arXiv 2602.04475).
- H on Co with the C donor hydrogenated (CH_HM): dG_H on the metal of the hydrogenated site, against +0.03 eV.
- Acid demetalation with H capping of the emptied site (Meng & Zhong, Angew. Chem. 2026), at pH 0:
  M@site + n H+ -> site-H_n + M^n+ (n = 2 for Co(II): no net electrons, so independent of the potential),
  dG_dem = -E_f + n e U0(M/M^n+) + dG_cap(n),  dG_cap(n) = E(site-H_n) - E(site) - n/2 E(H2) + n 0.24 eV.
  With four H the site also takes up two electrons: dG_dem(4, U) = -E_f + 2 e U0 + dG_cap(4) + 2 e U (U vs RHE).
  dG_dem > 0: the metal stays. The plain criterion U_diss > 0 (dft/analyse.py) leaves the H capping out.
The 0.24 eV per H is the zero-point and entropy correction of dG_H, an approximation for N-H, O-H and C-H.
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from analyse import DG_CORR, DISSOLUTION, RUNS, scf_energy, final_atoms   # noqa: E402

DATA = os.path.join(HERE, "..", "data")
U_HER = (0.0, -0.2)            # V vs RHE: equilibrium and a typical H2-evolution potential


def energy(path):
    return scf_energy(os.path.join(RUNS, path, "scf.out"))[0]


def main():
    e_h2 = energy("_ref/H2")
    stab = pd.read_csv(os.path.join(DATA, "stability_graphene.csv")).set_index("motif")
    rows = []
    site = "Co-N2OC_DV_NONC"
    e_bare, e_hm = energy(site + "/bare"), energy(site + "/H_M")
    e_hc = energy(site + "/Hdonor_C")
    for name, ref, what in (("Hdonor_C", e_bare, "H on the C donor"), ("Hdonor_O", e_bare, "H on the O donor"),
                            ("CH_HM", e_hc, "H on Co, C donor hydrogenated")):
        e = energy("%s/%s" % (site, name))
        rows.append({"site": "Co-N2OC:DV:NONC", "quantity": "dG_H", "state": name, "what": what,
                     "value": e - ref - e_h2 / 2 + DG_CORR})
    rows.append({"site": "Co-N2OC:DV:NONC", "quantity": "dG_H", "state": "H_M", "what": "H on Co (reference)",
                 "value": e_hm - e_bare - e_h2 / 2 + DG_CORR})
    for motif, env in (("Co-N2OC:DV:NONC", "N2OC_DV_NONC"), ("Co-N4:DV", "N4_DV")):
        e_site = energy("_ref/site_" + env)
        e_f = stab.E_f.get(motif, np.nan)
        u0, n = DISSOLUTION["Co"]
        for k in (2, 4):
            e_cap = energy("_ref/site_%s_H%d" % (env, k)) if os.path.isdir(os.path.join(RUNS, "_ref", "site_%s_H%d" % (env, k))) else None
            if e_cap is None:
                continue
            g_cap = e_cap - e_site - k / 2 * e_h2 + k * DG_CORR
            rows.append({"site": motif, "quantity": "dG_cap", "state": "H%d" % k,
                         "what": "H capping of the emptied site (%d H)" % k, "value": g_cap})
            for u in (U_HER if k == 4 else (0.0,)):
                rows.append({"site": motif, "quantity": "dG_dem", "state": "H%d" % k, "U_RHE": u,
                             "what": "acid demetalation, pH 0%s" % ("" if k == 2 else ", %+.1f V vs RHE" % u),
                             "value": -e_f + n * u0 + g_cap + (k - n) * u})
    # escape energy in DFT: the end state of the xTB escape path relaxed in DFT (dft/make_inputs.py escape_end)
    import json
    for motif in ("Zn-NS2:SV", "Ni-NS2:SV", "Cu-N3:SV"):
        d = motif.replace(":", "_")
        if not os.path.isdir(os.path.join(RUNS, d, "ESC_end")):
            continue
        neb = os.path.join(DATA, "screen", "xtb_neb", d + ".json")
        neb = json.load(open(neb)) if os.path.exists(neb) else {}
        rows.append({"site": motif, "quantity": "E_escape", "state": "ESC_end",
                     "what": "escape energy, DFT (xTB: %+.2f eV, NEB barrier %.2f eV)" % (
                         neb.get("reaction", np.nan), neb.get("barrier", np.nan)),
                     "value": energy(d + "/ESC_end") - energy(d + "/bare")})
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(DATA, "dft_conditions.csv"), index=False, float_format="%.3f")
    lines = ["# Co-N2OC under H2-evolution conditions, and DFT escape energies", "",
             "| site | quantity | state | what | value (eV) |", "|---|---|---|---|---|"]
    lines += ["| %s | %s | %s | %s | %s |" % (r.site, r.quantity, r.state, r.what,
                                             "%+.2f" % r.value if np.isfinite(r.value) else "running")
              for r in df.itertuples()]
    lines += ["", "dG_dem > 0: the metal stays in the site; dG_H on a donor below about 0 eV: the donor is "
              "hydrogenated under H2 evolution. E_escape: metal from its site to an adatom on clean graphene about 5 A away "
              "(metal kept on its hollow site laterally), the lower bound of the escape barrier."]
    open(os.path.join(DATA, "dft_conditions.md"), "w").write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
