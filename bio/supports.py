"""Map every first shell of the complete PDB sweep onto two more supports: the MoS2 basal plane and silica.

    python bio/supports.py      # data/complete/sites.parquet -> data/complete/supports.parquet (per site),
                                #   data/complete/motifs_supports.parquet, data/supports_summary.md

On graphene the vacancy supplies the donors (N/O/S/C dopants, bio/first_shell.py). Here the support does:
  MoS2    sulfur. A metal over a Mo (or a hollow) of the basal plane sits on a face of three lattice S,
          3.16 A apart (2H-MoS2).
  MoS2vS  the metal in a sulfur vacancy: it binds the three Mo below (3.16 A apart), with the six S around the
          vacancy at about its own height. In the PDB only cluster corners look like that, e.g. an Fe of a
          [4Fe-4S] cube: three Fe on one side, its three mu3-S around it, the Cys on top. Face = three metal
          neighbours; donors less than RING_LIFT above their plane take the place of the S around the
          vacancy (fine if they are S, removed otherwise).
  SiO2    oxygen. Si-OH + M-L -> Si-O-M + HL: up to three siloxide O on one face (tripodal grafting), about
          2.6-3.2 A apart like vicinal and geminal silanols. Amorphous silica adapts, so the window is wide.
Per site and support the face is a set of support-element donors (3, else 2, else 1) whose pairwise angles
stay below FACE_MAX_ANGLE; among faces of one size the one that leaves fewest donors in the surface wins,
then the lowest score. score = RMS deviation of the face's X...X distances from the support's (A; 0 for one
donor). The other donors stay on the metal (apical: a thiol, a silane linker or the solvent has to supply
them); donors less than MIN_LIFT above the plane of the face atoms would sit in the surface and are removed.
motif = "Fe-S3(MoS2)+ax:S": face from the support, then the enzyme's apical donors; the bare single-atom
site is the part before "+ax".
"""
import itertools
import json
import os
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from first_shell import FACE_MAX_ANGLE, formula  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
COMPLETE = ROOT / "data" / "complete"
BIOLOGICAL = {"V", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W"}
SUPPORTS = {
    # element the support supplies, ideal X...X on the face (A), half-width of the flat window (A), and the
    # shortest X...X two support atoms can have (bonded pairs: persulfide, peroxide, both O of a carboxylate)
    "MoS2": {"el": "S", "xx": 3.16, "window": 0.0, "min_xx": 2.6},
    "SiO2": {"el": "O", "xx": 2.90, "window": 0.30, "min_xx": 2.4},
}
VACANCY = {"xx": 3.16, "window": 0.0}   # Mo...Mo around an S vacancy (A)
MIN_LIFT = 1.0     # A; an apical donor must sit at least this far above the plane of the face atoms
RING_LIFT = 2.2    # A; S vacancy: donors lower than this above the Mo3 plane stand in for the S ring
# A; a score up to this counts as a fit: MoS2 at the dip of the bimodal S...S distribution of X-ray sites
# <= 2.0 A (octahedral/prismatic S faces ~3.2-3.3 A vs tetrahedral ones ~3.6-3.8 A); silica: inside the window
GOOD = {"MoS2": 0.25, "SiO2": 0.10, "MoS2vS": 0.25}
HQ_RESOLUTION = 2.0   # geometry statistics only from X-ray structures at this resolution or better
COLS = ["site_id", "entry_id", "metal", "cn", "is_representative", "resolution", "method", "cluster30",
        "host_entity", "reaction", "description", "graphene_motif", "dv_score", "sv_score", "n_metal_neighbors",
        "metal_neighbors", "shell_json"]


def face_fit(donors, support):
    """Best face of one support for one site (module docstring), or None without a donor of its element."""
    spec = SUPPORTS[support]
    vec = np.array([(d["dx"], d["dy"], d["dz"]) for d in donors], dtype=float)
    dist = np.linalg.norm(vec, axis=1)
    unit = vec / dist[:, None]
    ang = np.degrees(np.arccos(np.clip(unit @ unit.T, -1, 1)))
    cand = [i for i, d in enumerate(donors) if d["el"] == spec["el"]]
    best = None
    for k in (3, 2, 1):
        for face in itertools.combinations(cand, k):
            pairs = list(itertools.combinations(face, 2))
            if any(ang[i, j] > FACE_MAX_ANGLE or np.linalg.norm(vec[i] - vec[j]) < spec["min_xx"] for i, j in pairs):
                continue
            up = -unit[list(face)].sum(axis=0)
            norm = np.linalg.norm(up)
            if norm < 1e-6:
                continue
            up /= norm
            height = vec @ up
            floor = height[list(face)].mean()
            rest = [i for i in range(len(donors)) if i not in face]
            removed = [i for i in rest if height[i] - floor < MIN_LIFT]
            dev = [max(0.0, abs(np.linalg.norm(vec[i] - vec[j]) - spec["xx"]) - spec["window"]) for i, j in pairs]
            score = float(np.sqrt(np.mean(np.square(dev)))) if dev else 0.0
            key = (len(removed), score)
            if best is None or key < best[0]:
                best = (key, face, rest, removed, score)
        if best:
            break
    if best is None:
        return None
    _, face, rest, removed, score = best
    apical = [i for i in rest if i not in removed]
    els = lambda idx: formula(donors[i]["el"] for i in idx)  # noqa: E731
    cls = lambda idx: " ".join(sorted(donors[i]["cls"] for i in idx))  # noqa: E731
    xx = [np.linalg.norm(vec[i] - vec[j]) for i, j in itertools.combinations(face, 2)]
    return {"k": len(face), "face_classes": cls(face), "score": round(score, 3),
            "xx_mean": round(float(np.mean(xx)), 3) if xx else np.nan,
            "apical": els(apical), "apical_classes": cls(apical), "removed": cls(removed)}


def vacancy_fit(donors, metals):
    """Three metal neighbours on one face as the Mo3 of an S vacancy (module docstring), or None."""
    if len(metals) < 3:
        return None
    mvec = np.array([(m["dx"], m["dy"], m["dz"]) for m in metals], dtype=float)
    munit = mvec / np.linalg.norm(mvec, axis=1)[:, None]
    best = None
    for face in itertools.combinations(range(len(metals)), 3):
        cos = [munit[i] @ munit[j] for i, j in itertools.combinations(face, 2)]
        if min(cos) < np.cos(np.radians(FACE_MAX_ANGLE)):
            continue
        dev = [np.linalg.norm(mvec[i] - mvec[j]) - VACANCY["xx"] for i, j in itertools.combinations(face, 2)]
        score = float(np.sqrt(np.mean(np.square(dev))))
        if best is None or score < best[0]:
            best = (score, face)
    if best is None:
        return None
    score, face = best
    up = -munit[list(face)].sum(axis=0)
    up /= np.linalg.norm(up)
    floor = (mvec[list(face)] @ up).mean()
    vec = np.array([(d["dx"], d["dy"], d["dz"]) for d in donors], dtype=float).reshape(-1, 3)
    lift = vec @ up - floor
    ring = [i for i in range(len(donors)) if lift[i] < RING_LIFT]
    apical = [i for i in range(len(donors)) if lift[i] >= RING_LIFT]
    xx = [np.linalg.norm(mvec[i] - mvec[j]) for i, j in itertools.combinations(face, 2)]
    return {"k": 3, "face_classes": " ".join(sorted(metals[i]["el"] for i in face)), "score": round(score, 3),
            "xx_mean": round(float(np.mean(xx)), 3),
            "ring": formula(donors[i]["el"] for i in ring),
            "apical": formula(donors[i]["el"] for i in apical),
            "apical_classes": " ".join(sorted(donors[i]["cls"] for i in apical)),
            "removed": " ".join(sorted(donors[i]["cls"] for i in ring if donors[i]["el"] != "S"))}


def map_site(args):
    sid, metal, js = args
    rec = json.loads(js)
    donors, metals = rec["donors"], rec.get("metals") or []
    out = {"site_id": sid}
    fits = {sup: face_fit(donors, sup) if donors else None for sup in SUPPORTS}
    fits["MoS2vS"] = vacancy_fit(donors, metals)
    for sup, fit in fits.items():
        if fit is None:
            continue
        if sup == "MoS2vS":
            motif = f"{metal}-Mo3(MoS2vS)"
        else:
            motif = f"{metal}-{SUPPORTS[sup]['el']}{fit['k'] if fit['k'] > 1 else ''}({sup})"
        out.update({f"{sup}_{k}": v for k, v in fit.items()})
        out[f"{sup}_site"] = motif
        out[f"{sup}_motif"] = motif + (f"+ax:{fit['apical']}" if fit["apical"] else "")
    return out


def motif_table(df, sup, set_name="complete"):
    """Per support motif: sites, entries and 30 % clusters over all copies and over representative copies.
    set_name "enzyme": only sites labelled with a reaction (HER/ORR/CO2RR/NRR enzymes of the sweep)."""
    col = f"{sup}_motif"
    if set_name == "enzyme":
        df = df[df.reaction.fillna("") != ""]
    d = df[df[col].notna()].assign(clu30=lambda x: x.cluster30.astype("string").fillna(
        x.entry_id + "_" + x.host_entity.astype("string").fillna("")))
    rep = dict(list(d[d.is_representative].groupby(col)))
    out = []
    for motif, g in d.groupby(col):
        r = rep.get(motif, g.iloc[:0])
        e = r if len(r) else g
        hq = e[(e.method == "X-ray") & (e.resolution <= HQ_RESOLUTION)]
        best = e.sort_values(["resolution", "entry_id"]).drop_duplicates("entry_id")
        out.append({
            "set_name": set_name, "kind": sup, "motif": motif, "site": g[f"{sup}_site"].iloc[0],
            "metal": g.metal.iloc[0],
            "n_sites": len(g), "n_entries": g.entry_id.nunique(), "n_clusters30": g.clu30.nunique(),
            "n_sites_rep": len(r), "n_entries_rep": r.entry_id.nunique(), "n_clusters30_rep": r.clu30.nunique(),
            "reactions": ";".join(sorted({x for v in g.reaction.fillna("") if v for x in v.split(";")})),
            "example_descriptions": " | ".join(dict.fromkeys(x for x in best.description.fillna("") if x))[:200],
            "example_pdb": " ".join(f"{a}({x:.2f})" if pd.notna(x) else a
                                    for a, x in zip(best.entry_id[:6], best.resolution[:6])),
            "best_resolution": e.resolution.min(),
            "n_sites_hq": len(hq),   # geometry columns below: X-ray <= HQ_RESOLUTION only (NaN without such sites)
            "median_score": round(float(hq[f"{sup}_score"].median()), 3) if len(hq) else np.nan,
            "median_xx": round(float(hq[f"{sup}_xx_mean"].median()), 3) if len(hq) else np.nan,
            "frac_fit": round(float((hq[f"{sup}_score"] <= GOOD[sup]).mean()), 3) if len(hq) else np.nan,
            "face_classes_top": e[f"{sup}_face_classes"].value_counts().index[0],
            "removed_top": e[f"{sup}_removed"].value_counts().index[0],
            "graphene_motif_top": e.graphene_motif.fillna("").value_counts().index[0],
            "frac_mononuclear": round(float((e.n_metal_neighbors == 0).mean()), 3)})
    m = pd.DataFrame(out)
    return m.sort_values(["n_entries_rep", "n_clusters30_rep", "n_sites"], ascending=False, ignore_index=True)


def md(df, cols):
    rows = ["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)]
    rows += ["| " + " | ".join("" if pd.isna(v) else str(v) for v in r) + " |" for r in df[cols].itertuples(index=False)]
    return "\n".join(rows)


NAMES = {"MoS2": "MoS2, metal on three lattice S", "MoS2vS": "MoS2, metal in an S vacancy (on three Mo)",
         "SiO2": "silica, metal on siloxide O"}


def summary(sites, motifs):
    rep = sites[sites.is_representative]
    bio = rep[rep.metal.isin(BIOLOGICAL)]
    enz = bio[bio.reaction.fillna("") != ""]
    hq = lambda x: x[(x.method == "X-ray") & (x.resolution <= HQ_RESOLUTION)]  # noqa: E731
    lines = ["# PDB sweep: first shells for MoS2 and silica single-atom sites", "",
             "Generated by `bio/supports.py` from `data/complete/sites.parquet` (every d-block site of the PDB, "
             "first shells of `bio/first_shell.py`). Counts are representative copies (one per group of identical "
             "sites in an entry); biological metals: " + ", ".join(sorted(BIOLOGICAL)) + "; enzymes: the HER, "
             "ORR, CO2RR and NRR enzyme labels of the sweep. Drafts that need expert review.", "",
             "Face = what the support supplies: three lattice S 3.16 A apart (MoS2), three Mo 3.16 A apart with "
             "the S ring of the vacancy around the metal (MoS2 S vacancy), up to three siloxide O 2.6-3.2 A apart "
             "(silica, Si-OH condensation). Score = RMS deviation of the face's X...X distances (A); a fit is "
             "score <= %.2f A (MoS2, S vacancy) or <= %.2f A (silica), counted on X-ray sites <= %.1f A only, "
             "because low-resolution and EM models of Fe4S4 cubes come out compressed. Apical donors (+ax) would "
             "have to come from a thiol, a silane linker or the solvent; the bare single-atom site keeps only the "
             "face." % (GOOD["MoS2"], GOOD["SiO2"], HQ_RESOLUTION), "",
             "| | all sites | biological metals | in HER/ORR/CO2RR/NRR enzymes |", "|---|---|---|---|",
             f"| representative sites | {len(rep)} | {len(bio)} | {len(enz)} |",
             f"| of these X-ray <= {HQ_RESOLUTION} A | {len(hq(rep))} | {len(hq(bio))} | {len(hq(enz))} |"]
    for sup in ("MoS2", "MoS2vS", "SiO2"):
        k, s = f"{sup}_k", f"{sup}_score"
        sel = lambda x: x[x[k] == 3]  # noqa: E731
        fit = lambda x: hq(x)[(hq(x)[k] == 3) & (hq(x)[s] <= GOOD[sup])]  # noqa: E731
        lines.append(f"| {NAMES[sup]}: three-donor face | {len(sel(rep))} | {len(sel(bio))} | {len(sel(enz))} |")
        lines.append(f"| ... X-ray <= {HQ_RESOLUTION} A and fit | {len(fit(rep))} | {len(fit(bio))} | "
                     f"{len(fit(enz))} |")
    cols = ["motif", "n_entries_rep", "n_clusters30_rep", "n_sites_hq", "median_xx", "frac_fit",
            "graphene_motif_top", "reactions", "example_descriptions"]
    for sup in ("MoS2", "MoS2vS", "SiO2"):
        three = (motifs.kind == sup) & motifs.metal.isin(BIOLOGICAL) & motifs.site.str.contains(r"3\(")
        m = motifs[three & (motifs.set_name == "complete")]
        lines += ["", f"## {NAMES[sup]}: three-donor faces, biological metals, most entries first", "",
                  md(m.head(15), cols)]
        e = motifs[three & (motifs.set_name == "enzyme")].sort_values(["n_clusters30_rep", "n_entries_rep"],
                                                                     ascending=False)
        lines += ["", "### the same, counting only sites of HER/ORR/CO2RR/NRR enzymes", "", md(e.head(12), cols)]
        g = m[(m.n_sites_hq >= 3) & (m.frac_fit >= 0.5)].sort_values(["n_clusters30_rep", "n_entries_rep"],
                                                                      ascending=False)
        lines += ["", f"### the same, geometry fits the support (>= 3 X-ray sites <= {HQ_RESOLUTION} A, at least "
                      "half of them fit)", "", md(g.head(12), cols)]
    return "\n".join(lines) + "\n"


def main():
    sites = pd.read_parquet(COMPLETE / "sites.parquet", columns=[c for c in COLS if c != "method"])
    sites = sites.merge(pd.read_parquet(COMPLETE / "entries.parquet", columns=["entry_id", "method"]),
                        on="entry_id", how="left")
    jobs = list(zip(sites.site_id, sites.metal, sites.shell_json))
    with Pool(max(1, (os.cpu_count() or 2) - 2)) as pool:
        mapped = pd.DataFrame(pool.map(map_site, jobs, chunksize=2000))
    out = sites.drop(columns="shell_json").merge(mapped, on="site_id", how="left")
    keep = ["site_id"] + [c for c in mapped.columns if c != "site_id"]
    out[keep].to_parquet(COMPLETE / "supports.parquet", index=False)
    motifs = pd.concat([motif_table(out, sup, s) for s in ("complete", "enzyme") for sup in [*SUPPORTS, "MoS2vS"]],
                       ignore_index=True)
    motifs.to_parquet(COMPLETE / "motifs_supports.parquet", index=False)
    text = summary(out, motifs)
    (ROOT / "data" / "supports_summary.md").write_text(text)
    print(text)


if __name__ == "__main__":
    main()
