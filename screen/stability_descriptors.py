"""Which properties of a site predict how strongly the carbon holds its metal, and how easily the metal escapes?

    python screen/stability_descriptors.py    # -> data/xtb_stability_descriptors.md and .csv

Inputs: data/xtb_stability.csv (screen/xtb_stability.py table), data/stability_graphene.csv (DFT E_f), the PDB sweep
(data/complete/sites.parquet).
1. xTB against DFT: E_ins(xTB) against E_f(DFT) on the sites that have both: per metal, and with one constant per metal.
2. Within one metal (the metal's reference cancels): E_ins minus the metal's mean, fitted on the donor counts (relative
   to C), the vacancy type and the strain (metal height above the sheet, longest metal-donor bond); 5-fold
   cross-validation over the donor arrangements, R2 against predicting the mean, and each descriptor's weight.
3. Across metals: the mean E_f of a metal against its bulk cohesive energy (metals with a DFT reference only).
4. The enzyme prior: do donor sets that enzymes use with a metal hold it more strongly? Spearman correlation of E_ins
   (within the metal) with log(1 + PDB entries with that metal, vacancy type and donor set).
5. Escape energies against what survives an hour of pyrolysis (3.5 eV at 800 C, 4.2 eV at 1000 C).
6. Does a forced geometry cost stability? The enzyme-to-support mismatch of each model (screen/geometry_mismatch.py:
   angle change, metal shift out of the donor plane, bond-length change, donors lost) against E_ins within the
   metal, alone (Spearman) and on top of donors and vacancy type (cross-validated R2 on the same sites).
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
OUT_MD = os.path.join(DATA, "xtb_stability_descriptors.md")
OUT_CSV = os.path.join(DATA, "xtb_stability_descriptors.csv")
COHESIVE = {"V": 5.31, "Cr": 4.10, "Mn": 2.92, "Fe": 4.28, "Co": 4.39, "Ni": 4.44, "Cu": 3.49, "Zn": 1.35,
            "Mo": 6.82, "Ru": 6.74, "W": 8.90, "Pt": 5.84}        # eV/atom, experimental (Kittel)
FEATURES = ["n_N", "n_S", "n_O", "DV", "M_height", "max_M_X"]     # donor counts relative to C
KT = {800: 8.617e-5 * 1073.15, 1000: 8.617e-5 * 1273.15}
THRESHOLD = {t: kt * np.log(1e13 * 3600) for t, kt in KT.items()}  # barrier crossed within 1 h, eV


def cv_r2(X, y, groups, k=5, seed=0):
    """R2 of a least-squares fit in k-fold cross-validation, folds drawn over groups (donor arrangements)."""
    rng = np.random.default_rng(seed)
    uniq = np.array(sorted(set(groups)))
    fold = dict(zip(uniq, rng.permutation(len(uniq)) % k))
    pred = np.full(len(y), np.nan)
    for f in range(k):
        test = np.array([fold[g] == f for g in groups])
        A = np.c_[np.ones(int((~test).sum())), X[~test]]
        coef = np.linalg.lstsq(A, y[~test], rcond=None)[0]
        pred[test] = np.c_[np.ones(int(test.sum())), X[test]] @ coef
    return 1 - np.nansum((y - pred) ** 2) / np.sum((y - y.mean()) ** 2), pred


def enzyme_counts():
    """PDB entries per (metal, vacancy type, donor set) among the representative sites of the sweep."""
    s = pd.read_parquet(os.path.join(DATA, "complete", "sites.parquet"),
                        columns=["entry_id", "metal", "is_representative", "preferred_site", "dv_formula", "sv_formula"])
    s = s[s.is_representative & s.preferred_site.isin(["DV", "SV"])]
    s["formula"] = np.where(s.preferred_site == "DV", s.dv_formula, s.sv_formula)
    return s.groupby(["metal", "preferred_site", "formula"]).entry_id.nunique()


def geometry_part(t, rows):
    """Section 6: the enzyme-to-support mismatch of each model as a stability descriptor."""
    path = os.path.join(DATA, "geometry_mismatch.csv")
    if not os.path.exists(path):
        return []
    g = pd.read_csv(path)
    g["dh"] = (g.h_model - g.h_enzyme).abs()
    dd = [c[2:-6] for c in g if c.startswith("d_") and c.endswith("_model")]
    g["dd"] = pd.concat([(g["d_%s_model" % e] - g["d_%s_enzyme" % e]).abs() for e in dd], axis=1).mean(axis=1)
    m = t.merge(g[["motif", "angle_rms", "dh", "dd", "donors_lost", "enzyme_sites"]], on="motif")
    feats = {"angle_rms": "angles changed (RMS, deg)", "dh": "metal shift out of the donor plane (A)",
             "dd": "bond-length change (mean |model - enzyme|, A)", "donors_lost": "donors lost in the model"}
    lines = ["## 6. Does a forced geometry cost stability?", "",
             "%d sites have both an xTB E_ins and enzyme sites to compare with (screen/geometry_mismatch.py). Spearman "
             "with E_ins relative to the metal's mean (positive: the more the support changes the enzyme's geometry, "
             "the more weakly it holds the metal):" % len(m), ""]
    for f, label in feats.items():
        if m[f].nunique() > 1:
            r = spearmanr(m[f], m.dE)
            lines.append("- %s: %.2f (p = %.3f)" % (label, r.statistic, r.pvalue))
            rows.append({"part": "geometry", "descriptor": f, "weight": r.statistic})
    base = ["n_N", "n_S", "n_O", "DV"]
    mm = m.dropna(subset=base + list(feats))
    grp = (mm.site + "_" + mm.order).to_numpy()
    r2_base, _ = cv_r2(mm[base].to_numpy(float), mm.dE.to_numpy(), grp)
    r2_geo, _ = cv_r2(mm[base + list(feats)].to_numpy(float), mm.dE.to_numpy(), grp)
    lines += ["", "Cross-validated R2 on these %d sites: %.2f with donors and vacancy type, %.2f with the mismatch "
              "added." % (len(mm), r2_base, r2_geo), ""]
    rows += [{"part": "geometry", "descriptor": "cv_r2_base", "weight": r2_base},
             {"part": "geometry", "descriptor": "cv_r2_with_mismatch", "weight": r2_geo}]
    d = pd.read_csv(os.path.join(DATA, "stability_graphene.csv")).dropna(subset=["E_f"]).merge(
        g[["motif", "angle_rms", "dh", "donors_lost"]], on="motif")
    if len(d) >= 4:
        r = spearmanr(d.angle_rms, d.E_f)
        lines += ["DFT check on the %d sites with E_f and enzyme sites: Spearman of E_f with the angle change %.2f "
                  "(p = %.2f); %s." % (len(d), r.statistic, r.pvalue, ", ".join(
                      "%s E_f %+.2f eV, angles %.0f deg%s" % (x.motif, x.E_f, x.angle_rms,
                                                             ", %d donors lost" % x.donors_lost if x.donors_lost else "")
                      for x in d.sort_values("E_f").itertuples())), ""]
    return lines


def main():
    t = pd.read_csv(os.path.join(DATA, "xtb_stability.csv"))
    t = t[t.converged == True].copy()   # noqa: E712 (NaN-safe)
    t["DV"] = (t.site == "DV").astype(int)
    lines = ["# What holds a single metal atom in carbon: descriptors from GFN2-xTB", "",
             "E_ins = E(metal in the site) - E(empty site); more negative holds the metal more strongly. "
             "%d sites, %d metals, %d donor arrangements. Drafts: xTB, checked against %s." % (
                 len(t), t.metal.nunique(), t.groupby(["site", "order"]).ngroups, "the DFT formation energies"), ""]
    rows = []

    # 1. xTB against DFT
    dft = pd.read_csv(os.path.join(DATA, "stability_graphene.csv")).dropna(subset=["E_f"])
    p = dft.merge(t[["motif", "E_ins", "E_f_xtb"]], on="motif", how="inner")
    p["metal"] = p.motif.str.split("-").str[0]
    lines += ["## 1. xTB against DFT", "", "| site | E_f DFT (eV) | E_ins xTB (eV) | E_f from xTB (eV) |", "|---|---|---|---|"]
    lines += ["| %s | %+.2f | %+.2f | %s |" % (r.motif, r.E_f, r.E_ins, "%+.2f" % r.E_f_xtb if pd.notna(r.E_f_xtb) else "-")
              for r in p.sort_values(["metal", "E_f"]).itertuples()]
    within = []
    for m, g in p.groupby("metal"):
        if len(g) >= 2:
            within.append((m, len(g), spearmanr(g.E_f, g.E_ins).statistic, (g.E_f - g.E_ins).std()))
    ok = p.dropna(subset=["E_f_xtb"])
    lines += ["", "Within each metal (sites ordered the same way?): " + "; ".join(
        "%s %d sites, Spearman %.2f, spread of E_f - E_ins %.2f eV" % w for w in within) + ".",
        "With one constant per metal (in-sample) the mean error is %.2f eV over %d sites." % (
            (ok.E_f - ok.E_f_xtb).abs().mean(), len(ok)) if len(ok) else "", ""]

    # 2. within one metal
    t["dE"] = t.E_ins - t.groupby("metal").E_ins.transform("mean")
    d = t.dropna(subset=FEATURES)
    X = d[FEATURES].to_numpy(float)
    groups = (d.site + "_" + d.order).to_numpy()
    r2_all, _ = cv_r2(X, d.dE.to_numpy(), groups)
    r2_donor, _ = cv_r2(d[["n_N", "n_S", "n_O", "DV"]].to_numpy(float), d.dE.to_numpy(), groups)
    A = np.c_[np.ones(len(d)), (X - X.mean(0)) / X.std(0)]
    coef = np.linalg.lstsq(A, d.dE.to_numpy(), rcond=None)[0]
    lines += ["## 2. Within one metal", "",
              "E_ins relative to the metal's mean, %d sites, one model for all metals. Cross-validated R2 (5 folds over "
              "the arrangements): %.2f with donors, vacancy type and strain; %.2f with donors and vacancy type alone." % (
                  len(d), r2_all, r2_donor), "",
              "| descriptor | weight (eV per standard deviation) |", "|---|---|"]
    lines += ["| %s | %+.2f |" % (f, c) for f, c in zip(FEATURES, coef[1:])]
    rows += [{"part": "within metal", "descriptor": f, "weight": c} for f, c in zip(FEATURES, coef[1:])]
    rows += [{"part": "within metal", "descriptor": "cv_r2_all", "weight": r2_all},
             {"part": "within metal", "descriptor": "cv_r2_donors", "weight": r2_donor}]
    per = []
    for m, g in d.groupby("metal"):
        if len(g) >= 20:
            r2, _ = cv_r2(g[FEATURES].to_numpy(float), g.E_ins.to_numpy(), (g.site + "_" + g.order).to_numpy())
            best = g.nsmallest(3, "E_ins").motif.tolist()
            per.append("%s (%d sites): CV R2 %.2f, strongest %s" % (m, len(g), r2, ", ".join(best)))
    lines += ["", "Per metal: " + "; ".join(per) + ".", ""]

    # 3. across metals
    ref = p.groupby("metal").E_f.mean()
    if len(ref) >= 3:
        rho = spearmanr(ref.values, [COHESIVE[m] for m in ref.index]).statistic
        lines += ["## 3. Across metals", "", "Mean DFT E_f per metal against the bulk cohesive energy: " + ", ".join(
            "%s %.2f eV (E_coh %.2f)" % (m, v, COHESIVE[m]) for m, v in ref.items()) +
            "; Spearman %.2f over %d metals." % (rho, len(ref)), ""]

    # 4. the enzyme prior
    counts = enzyme_counts()
    t["pdb_entries"] = [counts.get((m, s, f), 0) for m, s, f in zip(t.metal, t.site, t.formula)]
    rho_all = spearmanr(np.log1p(t.pdb_entries), t.dE)
    used = t[t.pdb_entries > 0]
    lines += ["## 4. The enzyme prior", "",
              "Spearman of E_ins (within the metal) with log(1 + PDB entries of that metal with that donor set): "
              "%.2f (p = %.3f, %d sites; %d of them with the donor set seen in enzymes). Negative means donor sets "
              "that enzymes use with a metal hold it more strongly." % (rho_all.statistic, rho_all.pvalue, len(t), len(used)),
              ""]
    rows.append({"part": "enzyme prior", "descriptor": "spearman", "weight": rho_all.statistic})

    # 5. escape energies
    if "E_escape" in t and t.E_escape.notna().any():
        e = t.dropna(subset=["E_escape"])
        lines += ["## 5. Escape energy (lower bound of the barrier out of the site)", "",
                  "An hour of pyrolysis gets an atom over %.1f eV at 800 C and %.1f eV at 1000 C (1e13 /s)." % (
                      THRESHOLD[800], THRESHOLD[1000]), "",
                  "| metal | sites | median (eV) | max (eV) | above %.1f eV |" % THRESHOLD[800], "|---|---|---|---|---|"]
        for m, g in e.groupby("metal"):
            lines.append("| %s | %d | %.2f | %.2f | %d |" % (m, len(g), g.E_escape.median(), g.E_escape.max(),
                                                           int((g.E_escape >= THRESHOLD[800]).sum())))
        lines.append("")
    lines += geometry_part(t, rows)
    open(OUT_MD, "w").write("\n".join(lines) + "\n")
    pd.DataFrame(rows).to_csv(OUT_CSV, index=False, float_format="%.3f")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
