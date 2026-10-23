"""Calibrate the GFN2-xTB scan against DFT and pick the candidates for full DFT treatment.

    python screen/calibrate.py          # -> data/xtb_calibration.csv, xtb_offsets.csv, xtb_candidates.csv

Pairs are matched by structure, not by name: the same metal, in-plane donors (canonical order N, S, O, C),
site, isomer and axial ligand, and the same H position. dG_H(DFT) = a dG_H(xTB) + b is fitted on those pairs;
the global line is printed for reference, but the ranking uses one offset per metal for H on the metal (xTB's
error is metal-dependent: about -1.5 eV on Cu, about 0 on Fe) and one offset for H on S. An offset is used only
from MIN_PAIRS pairs on: a single pair cannot be checked (left out, it predicts nothing), so such metals stay raw,
like metals without a pair. The offsets, their pair counts and leave-one-out errors go to data/xtb_offsets.csv, which
every later step reads (supports, reactions, database, website, explorer). Candidates: |calibrated dG_H| <= WINDOW, H stayed on the intended atom, both
relaxations converged, the metal still in its site (bare metal-donor distances <= 2.6 A), at most one carbon
donor (all-carbon vacancies are neither enzyme-like nor new), the best isomer per metal and site, no DFT yet.
"""
import os
import re
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
sys.path.insert(0, os.path.join(ROOT, "dft"))
from batches import GRAPHENE_BATCH1, GRAPHENE_BATCH2   # noqa: E402

WINDOW = 0.3            # eV, calibrated |dG_H| for a candidate (the near-optimal band plus the scan's error)
MIN_PAIRS = 2           # DFT pairs a metal needs before its offset is used
N_DFT = 24              # candidates that go to DFT
ORDER = ["N", "S", "O", "C", "Se", "Cl"]


def canonical(motif):
    """'Fe-C2S2:DV:cis5' and 'Fe-S2C2:DV:cis5' -> 'Fe|S2C2|DV|cis5|'."""
    head, *opts = motif.split(":")
    metal, donors = head.split("-")
    n = {el: int(k or 1) for el, k in re.findall(r"([A-Z][a-z]?)(\d*)", donors)}
    formula = "".join(el + (str(n[el]) if n[el] > 1 else "") for el in ORDER if el in n)
    site = opts[0] if opts else ""
    iso = next((o for o in opts[1:] if not o.startswith("ax=")), "")
    ax = next((o for o in opts[1:] if o.startswith("ax=")), "")
    return "|".join([metal, formula, site, iso, ax])


def offsets(pairs):
    """One offset per metal for H on the metal and one for H on S, each with its pair count and its leave-one-out
    error (every pair predicted from an offset fitted without it); used only from MIN_PAIRS pairs on."""
    rows = []
    groups = [("H_M", m, g) for m, g in pairs[pairs.state_dft == "H_M"].groupby("metal")]
    groups.append(("H_S", "", pairs[pairs.state_dft == "H_S"]))
    for state, metal, g in groups:
        diff = (g.dG_H_dft - g.dG_H_xtb).to_numpy()
        n = len(diff)
        loo = np.abs(diff - (diff.sum() - diff) / (n - 1)).mean() if n >= 2 else np.nan
        rows.append({"state": state, "metal": metal, "n_pairs": n, "offset": diff.mean(), "loo_mae": loo,
                     "used": n >= MIN_PAIRS})
    return pd.DataFrame(rows)


def main():
    x = pd.read_csv(os.path.join(DATA, "xtb_scan.csv"))
    d = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv")).dropna(subset=["dG_H"])
    x["key"] = x.motif.map(canonical) + "#" + x.state
    d["key"] = d.motif.map(canonical) + "#" + d.state
    ok = x.converged.astype(bool) & x.bare_converged.astype(bool) & x.dG_H.notna()
    pairs = d.merge(x[ok].drop_duplicates("key"), on="key", suffixes=("_dft", "_xtb"))
    a, b = np.polyfit(pairs.dG_H_xtb, pairs.dG_H_dft, 1)
    pairs["dG_H_cal"] = a * pairs.dG_H_xtb + b
    r = np.corrcoef(pairs.dG_H_xtb, pairs.dG_H_dft)[0, 1]
    rho = pairs[["dG_H_xtb", "dG_H_dft"]].rank().corr().iloc[0, 1]
    mae_raw = (pairs.dG_H_xtb - pairs.dG_H_dft).abs().mean()
    mae_cal = (pairs.dG_H_cal - pairs.dG_H_dft).abs().mean()
    print("%d pairs: dG_DFT = %.2f dG_xTB %+.2f eV; Pearson r %.2f, Spearman %.2f; MAE %.2f eV raw, %.2f eV calibrated"
          % (len(pairs), a, b, r, rho, mae_raw, mae_cal))
    cols = ["motif_dft", "state_dft", "dG_H_dft", "dG_H_xtb", "dG_H_cal", "H_partner_dft", "H_partner_xtb"]
    print(pairs[cols].sort_values("dG_H_dft").to_string(index=False, float_format="%.2f"))
    pairs[cols].to_csv(os.path.join(DATA, "xtb_calibration.csv"), index=False, float_format="%.3f")

    # the global line hides a metal-dependent error (xTB binds H ~1.2 eV too strongly on Ni and Cu, about right
    # on Fe): calibrate H on the metal with one offset per metal, H on S with one offset for all S sites
    pairs["metal"] = pairs.motif_dft.str.split("-").str[0]
    offs = offsets(pairs)
    offs.to_csv(os.path.join(DATA, "xtb_offsets.csv"), index=False, float_format="%.3f")
    used = offs[offs.used]
    off = used[used.state == "H_M"].set_index("metal").offset
    off_s = used[used.state == "H_S"].offset.iloc[0]
    print("offsets (eV; pairs; leave-one-out MAE):")
    print(offs.to_string(index=False, float_format="%.2f"))
    loo = (offs.n_pairs * offs.loo_mae)[offs.used].sum() / offs.n_pairs[offs.used].sum()
    print("leave-one-out MAE over the calibrated pairs: %.2f eV" % loo)
    x["metal"] = x.motif.str.split("-").str[0]
    x["dG_H_cal"] = np.where(x.state == "H_S", x.dG_H + off_s, x.dG_H + x.metal.map(off))
    done = {canonical(m) for m, _ in GRAPHENE_BATCH1 + GRAPHENE_BATCH2} | set(d.motif.map(canonical))   # DFT exists
    intended = np.where(x.state == "H_M", x.H_partner == "M", x.H_partner.astype(str).str.endswith("(donor)"))
    formula = x.motif.str.split("-").str[1].str.split(":").str[0]
    n_c = formula.str.extract(r"C(\d*)(?![a-z])")[0].map(lambda v: 0 if pd.isna(v) else int(v or 1))
    cand = x[ok & intended & (x.max_M_X_bare <= 2.6) & x.dG_H_cal.notna() & (x.dG_H_cal.abs() <= WINDOW)
             & (n_c < 2) & ~x.motif.map(canonical).isin(done)].copy()
    cand["env"] = formula + ":" + x.motif.str.split(":").str[1]
    cand["abs_cal"] = cand.dG_H_cal.abs()
    cand = cand.sort_values("abs_cal").drop_duplicates(["metal", "env"])     # best isomer per metal and site
    cand[["motif", "state", "dG_H", "dG_H_cal", "H_partner", "max_M_X_bare", "max_M_X"]].to_csv(
        os.path.join(DATA, "xtb_candidates.csv"), index=False, float_format="%.3f")
    selected = os.path.join(DATA, "xtb_selected.csv")     # the record of DFT batch 3: written once, then kept
    if not os.path.exists(selected):
        cand.head(N_DFT)[["motif", "state", "dG_H", "dG_H_cal"]].to_csv(selected, index=False, float_format="%.3f")
    print("\n%d candidate sites with |calibrated dG_H| <= %.2f eV and no DFT yet; the best %d:" % (len(cand), WINDOW, N_DFT))
    print(cand.head(N_DFT)[["motif", "state", "dG_H", "dG_H_cal"]].to_string(index=False, float_format="%.2f"))
    compare_library(x[ok & intended & (x.max_M_X_bare <= 2.6) & x.dG_H_cal.notna()], n_c)
    return a, b, r, rho, len(pairs)


def compare_library(x, n_c):
    """Do the donor arrangements that enzymes use (the scan) sit closer to the optimum than the other arrangements of
    the library? Calibrated metals only, at most one C donor (as in enzymes), each site at its best H position
    -> data/xtb_library_comparison.csv (per metal and all) and the Mann-Whitney p of |dG_H|."""
    from scipy.stats import mannwhitneyu
    if "source" not in x or not (x.source == "library").any():
        return None
    x = x[n_c.reindex(x.index) <= 1].assign(a=lambda t: t.dG_H_cal.abs())
    best = x.sort_values("a").drop_duplicates("motif")
    best = best[best.state.eq("H_M") | best.metal.isin(best[best.state == "H_M"].metal.unique())]
    rows = []
    for metal, g in list(best.groupby("metal")) + [("all", best)]:
        for source, h in g.groupby("source"):
            rows.append({"metal": metal, "source": source, "sites": len(h), "within_0.3": (h.a <= 0.3).mean(),
                         "within_0.2": (h.a <= 0.2).mean(), "median_abs_dG": h.a.median()})
    out = pd.DataFrame(rows)
    s, l = best[best.source == "scan"].a, best[best.source == "library"].a
    out["p_mannwhitney"] = np.where(out.metal == "all", mannwhitneyu(s, l).pvalue, np.nan)
    out.to_csv(os.path.join(DATA, "xtb_library_comparison.csv"), index=False, float_format="%.3f")
    print("\nenzyme arrangements (scan) against the other arrangements (library), calibrated metals, <= 1 C donor:")
    print(out.to_string(index=False, float_format="%.2f"))
    return out


if __name__ == "__main__":
    main()
