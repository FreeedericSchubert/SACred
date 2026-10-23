"""Results website: collects the current results into one JSON block and writes site/index.html.

    python dft/analyse.py && python site/build.py     # then publish site/index.html (images: figures/*.png)

The page (site/template.html) renders everything from that block, so a rebuild after new results is all it
takes. Values of MoS2 runs that are still relaxing come from their last step and are marked as such.
"""
import datetime
import json
import os
import re

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
MOS2 = os.path.join(ROOT, "dft", "runs", "mos2")
RY = 13.605693123
FROZEN = "22 Oct, 23:30"        # last pull from Curta
ENZYME_METALS = ["V", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W"]
MOS2_RUNS = [   # run, what sits on/under the sheet, H on
    ("31_mos2_h", "clean MoS₂", "S"), ("34_mos2_fe_h", "Fe–S₃ site, clean", "Fe"),
    ("36_gate_mos2_fe_h/q_-0.30", "Fe–S₃ site, plate −0.3 e (acceptor-like)", "Fe"),
    ("36_gate_mos2_fe_h/q_+0.30", "Fe–S₃ site, plate +0.3 e (donor-like)", "Fe")]


def clean(v):
    """JSON-safe value: NaN -> None, numpy scalars -> Python."""
    if isinstance(v, (float, np.floating)):
        return None if not np.isfinite(v) else round(float(v), 4)
    if isinstance(v, np.integer):
        return int(v)
    return v


ORDER = ["N", "S", "O", "C", "Se", "Cl"]


def donor_key(head):
    """'Fe-C2S2' and 'Fe-S2C2' -> 'Fe-S2C2' (donors in the order N, S, O, C)."""
    metal, donors = head.split("-")
    n = {el: int(k or 1) for el, k in re.findall(r"([A-Z][a-z]?)(\d*)", donors)}
    return metal + "-" + "".join(el + (str(n[el]) if n[el] > 1 else "") for el in ORDER if el in n)


def literature():
    """Verdicts of the literature check (data/literature_check.json), keyed by donor_key."""
    path = os.path.join(DATA, "literature_check.json")
    if not os.path.exists(path):
        return {}
    out = {}
    for r in json.load(open(path)):
        v = r["verdict"]
        out[donor_key(r["motif"])] = {
            "verdict": v["verdict"], "confidence": v["confidence"],
            "refuted": bool(r.get("challenge") and r["challenge"].get("refuted")),
            "closest": [{k: c.get(k) for k in ("title", "year", "doi", "coordination", "match")} for c in v["closest"][:2]]}
    return out


def graphene():
    lit_check = literature()
    d = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv"))
    nov = pd.read_csv(os.path.join(DATA, "novelty.csv")).set_index("motif")
    rows = []
    for r in d.dropna(subset=["dG_H"]).itertuples():
        head = r.motif.split(":")[0]
        notes = []
        if isinstance(r.bare_states_pending, str) and r.bare_states_pending.strip():
            notes.append("another bare configuration still running")
        if r.state == "H_M" and r.H_partner != "M":
            notes.append("H moved to %s" % r.H_partner)
        if r.max_M_X > 2.5:
            notes.append("a metal–donor bond opened to %.2f Å" % r.max_M_X)
        origin = r.enzyme
        if origin.startswith("xTB scan"):
            origin = "picked by the xTB scan"
        elif origin.startswith("benchmark"):
            origin = "benchmark: the most studied HER single-atom site"
        rows.append({k: clean(v) for k, v in {
            "motif": r.motif, "state": r.state, "dG_H": r.dG_H, "origin": origin,
            "lit": int(nov.openalex_all[head]) if head in nov.index else None,
            "E_f": getattr(r, "E_f", np.nan), "U_diss": getattr(r, "U_diss", np.nan),
            "provisional": bool(notes and notes[0].startswith("another")), "notes": "; ".join(notes),
            "lit_check": (lit_check.get(donor_key(head)) or {}).get("verdict")}.items()})
    rows.sort(key=lambda x: abs(x["dG_H"]))
    return rows, int(d.dG_H.isna().sum()), int(d.motif.nunique())


def last_energy(path):
    """(energy in eV of the last SCF, finished?) of a relax.out, or (None, False)."""
    if not os.path.exists(path):
        return None, False
    text = open(path).read()
    e = re.findall(r"!\s+total energy\s+=\s+(-?\d+\.\d+)", text)
    return (float(e[-1]) * RY if e else None), "End final coordinates" in text


def mos2():
    """dG_H on MoS2: finished values from data/dgh_mos2.csv; running ones from their last step."""
    refs = json.load(open(os.path.join(ROOT, "dft", "mos2_refs.json")))
    done = pd.read_csv(os.path.join(DATA, "dgh_mos2.csv")).set_index("run")
    e_h2, _ = last_energy(os.path.join(MOS2, refs["h2"], "relax.out"))
    rows = []
    for run, label, on in MOS2_RUNS:
        value, status = None, "not submitted"
        if run in done.index and pd.notna(done.dG_H.get(run)):
            value, status = done.dG_H[run], "relaxed"
        elif os.path.exists(os.path.join(MOS2, run, "relax.out")):
            e, finished = last_energy(os.path.join(MOS2, run, "relax.out"))
            host = refs["hosts"].get(run)
            status = "stopped before convergence"
            if e is not None and host and e_h2:
                value, status = e - host["E_Ry"] * RY - e_h2 / 2 + 0.24, "stopped before convergence, value from the last step"
        elif os.path.exists(os.path.join(MOS2, run, "relax.in")) and not run.startswith("36_"):
            status = "not run"
        rows.append({"run": run, "label": label, "on": on, "dG_H": clean(value), "status": status})
    return rows


def sweep():
    s = pd.read_parquet(os.path.join(DATA, "complete", "sites.parquet"), columns=["entry_id", "metal", "is_representative"])
    m = pd.read_parquet(os.path.join(DATA, "complete", "motifs.parquet"))
    top = m[(m.kind == "graphene") & m.metal.isin(ENZYME_METALS)].nlargest(10, "n_entries_rep")
    val = open(os.path.join(DATA, "validation", "validation.md")).read()
    grab = lambda label: re.search(re.escape(label) + r"[^|]*\|\s*([\d.,]+ ?%?)", val).group(1).strip()
    return {
        "sites": int(len(s)), "entries": int(s.entry_id.nunique()), "representative": int(s.is_representative.sum()),
        "entries_swept": int(len(pd.read_parquet(os.path.join(DATA, "complete", "entries.parquet"), columns=["entry_id"]))),
        "enzyme_sites": int(len(pd.read_csv(os.path.join(DATA, "enzyme_sites.csv"), usecols=[0]))),
        "exact": grab("exact agreement"), "precision": grab("precision"), "recall": grab("recall"),
        "top": [{"motif": t.motif, "entries": int(t.n_entries_rep), "example": str(t.example_pdb).split()[0] if t.example_pdb else ""}
                for t in top.itertuples()]}


def funnel():
    """From the PDB to the near-optimal sites: (count, label, kind); kind = merge (nothing lost), select, add, compute."""
    s = pd.read_parquet(os.path.join(DATA, "complete", "sites.parquet"),
                        columns=["entry_id", "metal", "is_representative", "preferred_site", "dv_formula", "dv_isomer",
                                 "sv_formula", "formula"])
    rep = s[s.is_representative & s.metal.isin(ENZYME_METALS)]
    mapped = rep[rep.preferred_site.isin(["DV", "SV"])]
    env = np.where(mapped.preferred_site == "DV", mapped.dv_formula + ":DV:" + mapped.dv_isomer.fillna(""),
                   mapped.sv_formula + ":SV")
    n_env = pd.Series(mapped.entry_id.values, index=env).groupby(level=0).nunique()
    jobs = json.load(open(os.path.join(DATA, "xtb_scan_jobs.json")))
    d = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv"))
    n_xtb_dft = len(pd.read_csv(os.path.join(DATA, "xtb_selected.csv")))
    steps = [
        (len(pd.read_parquet(os.path.join(DATA, "complete", "entries.parquet"), columns=["entry_id"])),
         "PDB entries with a transition metal", "select"),
        (len(s), "metal sites, every copy", "add"),
        (int(s.is_representative.sum()), "unique sites (one per set of identical copies)", "merge"),
        (len(rep), "unique sites of the nine enzyme metals", "select"),
        (rep.groupby(["metal", "formula"]).ngroups, "distinct first shells (metal and donor set)", "merge"),
        (len(n_env), "distinct graphene sites (donor set and vacancy type)", "merge"),
        (int(((n_env >= 5) & ~n_env.index.str.contains("Se|Cl")).sum()), "seen in at least 5 entries, N/S/O/C donors", "select"),
        (len({m for m, _ in jobs if m != "H2"}), "xTB scan: each of these with all nine metals", "add"),
        (len(pd.read_csv(os.path.join(DATA, "xtb_candidates.csv"))), "within ±0.3 eV after the xTB calibration", "select"),
        (n_xtb_dft, "of these, sent to DFT (compute time)", "select"),
        (int(d.motif.nunique()), "sites in DFT, with %d picked by hand from enzymes" % (d.motif.nunique() - n_xtb_dft), "add"),
        (int(d.dG_H.notna().sum()), "H positions finished in DFT, of %d" % len(d), "compute"),
        (int((d.dG_H.round(2).abs() <= 0.2).sum()), "within 0.2 eV of the optimum", "select")]
    return [{"n": int(n), "label": label, "kind": kind} for n, label, kind in steps]


def screen():
    x = pd.read_csv(os.path.join(DATA, "xtb_scan.csv"))
    lib = x[x.source == "library"]
    x = x[x.source != "library"]
    pairs = pd.read_csv(os.path.join(DATA, "xtb_calibration.csv"))
    o = pd.read_csv(os.path.join(DATA, "xtb_offsets.csv"))
    hm = o[o.used & (o.state == "H_M")]          # leave-one-out error of the per-metal offsets (screen/calibrate.py)
    jobs = json.load(open(os.path.join(DATA, "xtb_scan_jobs.json")))         # the scan: [motif, state] pairs
    return {"sites": len({m for m, _ in jobs if m != "H2"}), "structures": len(jobs),
            "with_dG": int(x.dG_H.notna().sum()), "pairs": int(len(pairs)),
            "mae_metal": round(float((hm.n_pairs * hm.loo_mae).sum() / hm.n_pairs.sum()), 2),
            "metals_calibrated": ", ".join(hm.metal), "library_sites": int(lib.motif.nunique()),
            "selected": int(len(pd.read_csv(os.path.join(DATA, "xtb_selected.csv")))), "compare": library_compare()}


def library_compare():
    """Enzyme arrangements against the rest of the library (screen/calibrate.py), for one sentence on the page."""
    path = os.path.join(DATA, "xtb_library_comparison.csv")
    if not os.path.exists(path):
        return None
    c = pd.read_csv(path)
    a = c[c.metal == "all"].set_index("source")
    per = c[c.metal != "all"].pivot(index="metal", columns="source", values="within_0.3").dropna()
    gap = per.scan - per.library
    return {"scan": round(float(a.loc["scan", "within_0.3"]), 2), "library": round(float(a.loc["library", "within_0.3"]), 2),
            "scan_median": round(float(a.loc["scan", "median_abs_dG"]), 2),
            "library_median": round(float(a.loc["library", "median_abs_dG"]), 2),
            "p": round(float(a.loc["scan", "p_mannwhitney"]), 3), "n_scan": int(a.loc["scan", "sites"]),
            "n_library": int(a.loc["library", "sites"]),
            "favoured": [m for m in gap.index if gap[m] >= 0.15], "reverse": [m for m in gap.index if gap[m] <= -0.15]}


SUPPORT_XX = {"MoS2": 3.16, "MoS2vS": 3.16, "SiO2": 2.9}   # what the support's face has (A), bio/supports.py


def supports():
    """MoS2 and silica: enzyme sites that fit each support (bio/supports.py) and xTB hydrogen binding on the bare
    single-atom sites (screen/xtb_supports.py)."""
    path = os.path.join(DATA, "complete", "motifs_supports.parquet")
    if not os.path.exists(path):
        return None
    m = pd.read_parquet(path)
    # three-donor faces, except for the enzyme sites on silica: grafting through one or two Si-O is common there
    m = m[m.metal.isin(ENZYME_METALS) & (m.site.str.contains(r"3\(") | ((m.kind == "SiO2") & (m.set_name == "enzyme")))]

    def rows(df, n, why):
        return [{"support": r.kind, "motif": r.motif, "why": why, "entries": int(r.n_entries_rep),
                 "clusters": int(r.n_clusters30_rep), "xx": clean(r.median_xx), "ideal": SUPPORT_XX[r.kind],
                 "fit": clean(r.frac_fit), "examples": [x for x in str(r.example_descriptions).split(" | ") if x][:2]}
                for r in df.head(n).itertuples()]

    comp, enz = m[m.set_name == "complete"], m[m.set_name == "enzyme"]
    rank = ["n_clusters30_rep", "n_entries_rep"]
    fits = comp[(comp.kind == "MoS2") & (comp.n_sites_hq >= 3) & (comp.frac_fit >= 0.5)].sort_values(rank, ascending=False)
    tables = {
        "MoS2": rows(fits, 4, "fits the S triangle")
                + rows(enz[enz.kind == "MoS2"].sort_values(rank, ascending=False), 3, "most common in the enzymes"),
        "MoS2vS": rows(enz[enz.kind == "MoS2vS"].sort_values(rank, ascending=False), 2, "most common in the enzymes"),
        "SiO2": rows(enz[enz.kind == "SiO2"].sort_values(rank, ascending=False), 3, "most common in the enzymes")
                + rows(comp[(comp.kind == "SiO2") & (comp.frac_fit >= 0.5)].sort_values(rank, ascending=False), 2,
                       "most common in the whole PDB")}
    x = pd.read_csv(os.path.join(DATA, "xtb_supports.csv")) if os.path.exists(os.path.join(DATA, "xtb_supports.csv")) else None
    xtb = [] if x is None else [{k: clean(v) for k, v in r.items()} for r in x.to_dict("records")]
    s = pd.read_parquet(os.path.join(DATA, "complete", "supports.parquet"),
                        columns=["site_id", "MoS2_k", "MoS2vS_k", "SiO2_k"])
    return {"tables": tables, "xtb": xtb,
            "faces": {k: int((s[k + "_k"] == 3).sum()) for k in ("MoS2", "MoS2vS", "SiO2")}}


REACTIONS = [   # column, label, "max" or "min" is best, equilibrium potential (V vs RHE), note
    ("U_ORR_cal", "O₂ → H₂O (ORR, 4 e⁻)", "max", 1.23, "calibrated against DFT"),
    ("U_H2O2_cal", "O₂ → H₂O₂ (2 e⁻)", "max", 0.70, "calibrated against DFT"),
    ("U_OER_cal", "H₂O → O₂ (OER)", "min", 1.23, "calibrated against DFT"),
    ("U_CO", "CO₂ → CO", "max", -0.10, "xTB only"),
    ("U_HCOOH", "CO₂ → HCOOH", "max", -0.12, "xTB only"),
    ("U_NRR", "N₂ → NH₃", "max", 0.06, "xTB only")]


def reactions():
    """Best sites per reaction from screen/xtb_reactions.py (limiting potentials, V vs RHE)."""
    path = os.path.join(DATA, "xtb_reactions.csv")
    if not os.path.exists(path):
        return None
    x = pd.read_csv(path)
    x = x[~x.calibration_only.astype(bool)]
    fit = json.load(open(os.path.join(DATA, "xtb_orr_calibration_fit.json")))
    origin = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv")).drop_duplicates("motif").set_index("motif").get("origin")
    best = []
    for col, label, sense, eq, note in REACTIONS:
        if col not in x:
            continue
        d = x.dropna(subset=[col]).sort_values(col, ascending=(sense == "min")).head(4)
        best.append({"label": label, "eq": eq, "note": note, "sense": sense, "n": int(x[col].notna().sum()),
                     "sites": [{"kind": r.kind, "site": r.site, "metal": r.metal, "U": clean(getattr(r, col)),
                                "dG_H": clean(r.dG_H), "dG_H_source": r.dG_H_source if pd.notna(r.dG_H_source) else "",
                                "origin": (origin.get(r.site) if origin is not None and r.kind == "graphene" else "")}
                               for r in d.itertuples()]})
    return {"best": best, "fit": fit, "sites": int(len(x))}


def main():
    rows, pending, n_sites = graphene()
    lit = literature()
    rows_by_key = {}
    for r in rows:
        rows_by_key.setdefault(donor_key(r["motif"].split(":")[0]), r)
    lit_rows = [dict(v, key=k, motif=rows_by_key[k]["motif"], dG_H=rows_by_key[k]["dG_H"]) for k, v in lit.items() if k in rows_by_key]
    lit_rows.sort(key=lambda r: abs(r["dG_H"]))
    data = {"updated": datetime.datetime.now().strftime("%d %b %Y, %H:%M"), "frozen": FROZEN, "graphene": rows, "pending": pending, "literature": lit_rows,
            "dft_sites": n_sites, "mos2": mos2(), "sweep": sweep(), "screen": screen(), "funnel": funnel(),
            "supports": supports(), "reactions": reactions()}
    blob = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    page = open(os.path.join(HERE, "template.html")).read().replace("__DATA__", blob)
    with open(os.path.join(HERE, "index.html"), "w") as f:
        f.write(page)
    print("wrote site/index.html: %d DFT values, %d pending, MoS2 %s" % (
        len(rows), pending, ", ".join("%s %s" % (r["run"].split("_")[0], r["status"]) for r in data["mos2"])))


if __name__ == "__main__":
    main()
