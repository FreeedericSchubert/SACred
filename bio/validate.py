"""Technical validation of the complete site set (data/complete/sites.parquet) for the dataset release.

    python bio/validate.py metalc      # our donors vs the metal links the PDB annotates (struct_conn metalc)
    python bio/validate.py cutoffs     # how many motifs change when donor cutoffs or the DV threshold move
    python bio/validate.py summary     # -> data/validation/validation.md (+ the CSV tables next to it)

metalc: every mmCIF lists metal-ligand links in _struct_conn (conn_type_id metalc), set by the depositor or
by wwPDB annotation, independently of our distance rules. For each of our sites whose metal atom appears in
such links, the annotated partners are compared with our donors by identity (auth chain, residue number +
insertion code, residue name, atom name; altlocs and symmetry codes ignored for the match). Links between
two metals are left out, as are metals we do not treat as sites (Na, Mg, Ca, minor alternative positions).
Donors inside the metal's own component (heme N, the S of an Fe-S cluster) are never in _struct_conn,
whose links join two residues; they are counted separately (n_intra), not as disagreements.
An annotated partner we do not count as a donor gets the reason from our own records (weak contact and why,
an element without a donor cutoff, or not within 3.4 A); our matched distances are compared with the
annotated ones (pdbx_dist_value), which checks the symmetry handling.
cutoffs: donor formulas recomputed from the stored contacts with every element cutoff moved by -0.2 to
+0.2 A, and the share of CN >= 4 sites whose DV/SV choice sits within 10 deg of the 55 deg threshold.
Annotation is not ground truth: it is incomplete for many entries and sometimes generous, so the numbers
below measure agreement, and the disagreements are listed by kind.
"""
import json
import os
import sys
from collections import Counter
from multiprocessing import Pool

import gemmi
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
PDB_DIR = os.path.join(DATA, "pdb")
OUT = os.path.join(DATA, "validation")
SITES = os.path.join(DATA, "complete", "sites.parquet")
WORKERS = 8                     # of 14 local cores, at low priority (os.nice), so the machine stays usable
CUTOFF = {"N": 2.6, "O": 2.6, "S": 2.8, "Se": 2.9, "C": 2.3, "Cl": 2.8}     # as bio/first_shell.py
SHIFTS = (-0.2, -0.1, 0.1, 0.2)
DV_MAX, DV_BAND = 55.0, 10.0
COLS = ["conn_type_id", "ptnr1_auth_asym_id", "ptnr1_auth_seq_id", "pdbx_ptnr1_PDB_ins_code", "ptnr1_label_comp_id",
        "ptnr1_label_atom_id", "ptnr1_symmetry", "ptnr2_auth_asym_id", "ptnr2_auth_seq_id", "pdbx_ptnr2_PDB_ins_code",
        "ptnr2_label_comp_id", "ptnr2_label_atom_id", "ptnr2_symmetry", "pdbx_dist_value"]


def key(chain, seq, icode, res, atom):
    icode = "" if icode in ("?", ".", None) else icode
    return "%s|%s%s|%s|%s" % (chain, seq, icode, res, atom)


def links(entry):
    """All metalc rows of one entry: (entry, partner-1 key, partner-2 key, symmetry flag, distance)."""
    try:
        block = gemmi.cif.read(os.path.join(PDB_DIR, entry + ".cif.gz")).sole_block()
        tab = block.find("_struct_conn.", COLS)
    except Exception as exc:
        return [(entry, None, None, None, None, "read: %s" % exc)]
    rows = []
    for r in tab:
        if r[0] != "metalc":
            continue
        k1 = key(r[1], r[2], r[3], r[4], gemmi.cif.as_string(r[5]))
        k2 = key(r[7], r[8], r[9], r[10], gemmi.cif.as_string(r[11]))
        sym = r[6] not in ("1_555", "?", ".") or r[12] not in ("1_555", "?", ".")
        dist = float(r[13]) if r[13] not in ("?", ".") else np.nan
        rows.append((entry, k1, k2, sym, dist, ""))
    return rows


def our_sites():
    s = pd.read_parquet(SITES, columns=["site_id", "entry_id", "metal", "chain", "auth_seq_id", "icode", "metal_res",
                                        "metal_atom", "is_representative", "resolution", "method", "cn", "formula",
                                        "graphene_motif", "preferred_site", "dv_score", "shell_json"])
    s["mkey"] = [key(c, q, i, r, a) for c, q, i, r, a in zip(s.chain, s.auth_seq_id, s.icode, s.metal_res, s.metal_atom)]
    return s


def contacts_of(js):
    rec = json.loads(js)
    out = {}
    for kind in ("donors", "weak", "other"):
        for c in rec.get(kind) or []:
            k = key(c["auth_chain"], c["seqid"], c["icode"], c["resname"], c["atom"])
            if k not in out or kind == "donors":
                out[k] = (kind, c.get("reason", ""), c["dist"], c["el"], c.get("cls", ""))
    return out


def mode_metalc():
    os.makedirs(OUT, exist_ok=True)
    s = our_sites()
    entries = sorted(s.entry_id.unique())
    cache = os.path.join(OUT, "metalc_links.parquet")       # reading 26 GB of mmCIF takes minutes: once
    if os.path.exists(cache):
        L = pd.read_parquet(cache)
    else:
        with Pool(WORKERS, initializer=os.nice, initargs=(19,)) as pool:
            rows = [r for chunk in pool.imap_unordered(links, entries, chunksize=64) for r in chunk]
        L = pd.DataFrame(rows, columns=["entry_id", "k1", "k2", "symmetry", "dist", "error"])
        L.to_parquet(cache, index=False)
    errors = L[L.error != ""]
    L = L[L.error == ""].drop(columns="error")
    site_of = dict(zip(s.entry_id + "#" + s.mkey, s.site_id))
    # orient every link as (site, partner); metal-metal links (both ends sites) are left out
    a = (L.entry_id + "#" + L.k1).map(site_of)
    b = (L.entry_id + "#" + L.k2).map(site_of)
    one = a.notna() ^ b.notna()
    L = L[one].assign(site_id=np.where(a[one].notna(), a[one], b[one]),
                      partner=np.where(a[one].notna(), L.k2[one], L.k1[one]))
    ann = L.groupby("site_id").agg(partners=("partner", lambda p: sorted(set(p))),
                                   dists=("dist", list), pkeys=("partner", list), sym=("symmetry", "any"))
    out, reasons, ddiff = [], Counter(), []
    for r in s.itertuples():
        mine = contacts_of(r.shell_json)
        own = r.mkey.rsplit("|", 1)[0]                      # chain|seq|residue of the metal itself
        donors = {k for k, v in mine.items() if v[0] == "donors"}
        ours = {k for k in donors if k.rsplit("|", 1)[0] != own}
        a_set = set(ann.partners[r.site_id]) if r.site_id in ann.index else None
        row = {"site_id": r.site_id, "metal": r.metal, "resolution": r.resolution, "method": r.method,
               "is_representative": r.is_representative, "cn": r.cn, "n_intra": len(donors) - len(ours),
               "annotated": a_set is not None}
        if a_set is not None:
            both = ours & a_set
            missed = a_set - ours
            extra = ours - a_set
            why = []
            for k in missed:
                v = mine.get(k)
                tag = "not within 3.4 A" if v is None else ("%s: %s" % (v[0], v[1] or "element without cutoff"))
                why.append(tag)
                reasons[tag] += 1
            for k, d in zip(ann.pkeys[r.site_id], ann.dists[r.site_id]):
                if k in both and np.isfinite(d):
                    ddiff.append(mine[k][2] - d)
            row.update(n_ours=len(ours), n_annotated=len(a_set), n_both=len(both), n_missed=len(missed),
                       n_extra=len(extra), exact=len(missed) == 0 and len(extra) == 0, missed_reasons="; ".join(why),
                       extra_classes=" ".join(sorted(mine[k][4] for k in extra)))
        out.append(row)
    V = pd.DataFrame(out)
    V.to_parquet(os.path.join(OUT, "metalc_sites.parquet"), index=False)
    pd.Series(reasons).sort_values(ascending=False).rename("count").to_csv(os.path.join(OUT, "metalc_missed_reasons.csv"))
    np.save(os.path.join(OUT, "metalc_distance_diff.npy"), np.array(ddiff))
    errors.to_csv(os.path.join(OUT, "metalc_read_errors.csv"), index=False)
    print("%d entries read (%d errors), %d metalc links on our sites" % (len(entries), len(errors), len(L)))
    print_metalc(V, np.array(ddiff))


def print_metalc(V, dd):
    A = V[V.annotated]
    print("sites with annotated links: %d of %d (%.1f %%)" % (len(A), len(V), 100 * len(A) / len(V)))
    print("exact donor-set agreement: %.1f %%" % (100 * A.exact.mean()))
    print("donor level: precision %.3f (our donors annotated), recall %.3f (annotated partners we count)" % (
        A.n_both.sum() / max(A.n_ours.sum(), 1), A.n_both.sum() / max(A.n_annotated.sum(), 1)))
    if len(dd):
        print("our distance - annotated distance: median %+.3f A, 99 %% within %.3f A" % (
            np.median(dd), np.quantile(np.abs(dd), 0.99)))


def mode_cutoffs():
    os.makedirs(OUT, exist_ok=True)
    s = our_sites()
    base, changed = {}, Counter()
    for r in s.itertuples():
        rec = json.loads(r.shell_json)
        cand = [(c["el"], c["dist"], True) for c in rec.get("donors") or []]
        cand += [(c["el"], c["dist"], False) for c in rec.get("weak") or [] if c.get("reason") == "beyond_cutoff"]
        def formula(shift):
            n = Counter(el for el, d, is_donor in cand if el in CUTOFF and d <= CUTOFF[el] + shift
                        and (is_donor or shift > 0))
            return "".join("%s%s" % (el, n[el] if n[el] > 1 else "") for el in ["N", "S", "O", "C", "Se", "Cl"] if n[el])
        f0 = formula(0.0)
        for sh in SHIFTS:
            if formula(sh) != f0:
                changed[sh] += 1
    n = len(s)
    rows = [{"cutoff_shift_A": sh, "sites_with_other_formula": changed[sh], "share": changed[sh] / n} for sh in SHIFTS]
    near = s[(s.cn >= 4) & s.dv_score.notna()]
    band = ((near.dv_score - DV_MAX).abs() <= DV_BAND).mean()
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "cutoff_sensitivity.csv"), index=False)
    json.dump({"cn4_sites": int(len(near)), "share_within_10deg_of_threshold": float(band)},
              open(os.path.join(OUT, "dv_threshold_sensitivity.json"), "w"), indent=1)
    print(pd.DataFrame(rows).to_string(index=False))
    print("CN >= 4 sites within %g deg of the %g deg DV threshold: %.1f %%" % (DV_BAND, DV_MAX, 100 * band))


def mode_summary():
    V = pd.read_parquet(os.path.join(OUT, "metalc_sites.parquet"))
    A = V[V.annotated]
    dd = np.load(os.path.join(OUT, "metalc_distance_diff.npy"))
    miss = pd.read_csv(os.path.join(OUT, "metalc_missed_reasons.csv"), index_col=0)["count"]
    extra = A.extra_classes.str.split().explode().dropna().value_counts()
    bio = A[A.is_representative & A.metal.isin(["V", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W"])]
    per_metal = bio.groupby("metal").agg(sites=("exact", "size"), exact=("exact", "mean"), both=("n_both", "sum"),
                                         annotated=("n_annotated", "sum"))
    cut = pd.read_csv(os.path.join(OUT, "cutoff_sensitivity.csv"))
    dv = json.load(open(os.path.join(OUT, "dv_threshold_sensitivity.json")))
    md = ["# Technical validation", "",
          "Generated by `bio/validate.py` from `data/complete/sites.parquet` and the `_struct_conn` (metalc) records of "
          "the same mmCIF files. Annotation is not ground truth: it is incomplete for many entries and sometimes "
          "generous, so these numbers measure agreement.", "",
          "## Agreement with the PDB's metal-link annotation", "",
          "| measure | value |", "| --- | --- |",
          "| sites with annotated links | %d of %d (%.1f %%) |" % (len(A), len(V), 100 * len(A) / len(V)),
          "| exact agreement of the inter-residue donor set | %.1f %% |" % (100 * A.exact.mean()),
          "| precision (our donors that are annotated) | %.3f |" % (A.n_both.sum() / A.n_ours.sum()),
          "| recall (annotated partners we count as donors) | %.3f |" % (A.n_both.sum() / A.n_annotated.sum()),
          "| our minus annotated distance, median | %+.3f A |" % np.median(dd),
          "| 99 %% of matched distances within | %.3f A |" % np.quantile(np.abs(dd), 0.99),
          "| donors inside the metal's own component (not in struct_conn) | %d |" % V.n_intra.sum(), "",
          "Annotated partners we do not count as donors, by our reason:", "", "| reason | links |", "| --- | --- |"]
    md += ["| %s | %d |" % (k, v) for k, v in miss.items()]
    md += ["", "Our donors without an annotated link, by class (top 10):", "", "| class | donors |", "| --- | --- |"]
    md += ["| %s | %d |" % (k, v) for k, v in extra.head(10).items()]
    md += ["", "By metal (biological metals, representative copies):", "", "| metal | sites | exact | recall |",
           "| --- | --- | --- | --- |"]
    md += ["| %s | %d | %.3f | %.3f |" % (m, r.sites, r.exact, r.both / r.annotated) for m, r in per_metal.iterrows()]
    md += ["", "## Sensitivity to the thresholds", "", "| donor cutoffs moved by | sites with another donor formula |",
           "| --- | --- |"]
    md += ["| %+.1f A | %d (%.1f %%) |" % (r.cutoff_shift_A, r.sites_with_other_formula, 100 * r.share)
           for r in cut.itertuples()]
    md += ["", "CN >= 4 sites within 10 deg of the 55 deg double-vacancy threshold: %.1f %% of %d." % (
        100 * dv["share_within_10deg_of_threshold"], dv["cn4_sites"]), ""]
    open(os.path.join(OUT, "validation.md"), "w").write("\n".join(md))
    print("\n".join(md))


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "metalc":
        mode_metalc()
    elif mode == "cutoffs":
        mode_cutoffs()
    elif mode == "summary":
        mode_summary()
    else:
        sys.exit(__doc__)
