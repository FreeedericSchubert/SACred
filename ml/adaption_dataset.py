"""Fine-tuning data for an "enzyme site -> single-atom catalyst" model (Adaption AutoScientist or any SFT tool).

    python ml/adaption_dataset.py build                  # -> data/adaption/{train,test}.jsonl, *_meta.csv, README.md
    python ml/adaption_dataset.py score predictions.jsonl  # exact-match scores against data/adaption/test.jsonl

Every label comes from the pipeline, not from text a model wrote: the first shell (bio/first_shell.py), its graphene
site, and its MoS2 and silica faces (bio/supports.py). Two tasks, both answered as one JSON object:
  geometry   metal, donor atoms with distances and donor-metal-donor angles -> the mappings. The answer depends on
             the angles (square planar -> double vacancy, which face, which donor is left out).
  enzyme     enzyme name, organism and metal -> its first shell and the mappings; only (enzyme, metal) pairs whose
             sites all share one first shell, so the answer is unique.
Sites: representative copies of the nine enzyme metals, resolution <= 3.0 A, at least two donors, with a polymer
description. Split by 30 % sequence cluster (one hash per cluster): no test protein family is in training. Both
sets are capped per graphene motif so Fe4S4 cubes and zinc fingers do not drown the rest.
"""
import hashlib
import itertools
import json
import os
import sys

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COMPLETE = os.path.join(ROOT, "data", "complete")
OUT = os.path.join(ROOT, "data", "adaption")
METALS = ["V", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W"]
MAX_RESOLUTION = 3.0
TEST_EVERY = 10                     # one cluster in ten goes to the test set
CAP = {("geometry", "train"): 120, ("geometry", "test"): 12, ("enzyme", "train"): 80, ("enzyme", "test"): 8}
FIT = {"MoS2": 0.25, "SiO2": 0.10}  # score limits of bio/supports.py
FIELDS = ["first_shell", "graphene", "MoS2", "MoS2_fit", "MoS2_vacancy", "silica", "silica_fit"]
INSTRUCTION = ("Map this first shell onto single-atom catalyst sites: a graphene vacancy (SV = single, DV = double, "
               "+ax = axial donors left on the metal), the MoS2 basal plane (face of lattice S, or a metal in an S "
               "vacancy on three Mo) and silica (siloxide O). Answer with one JSON object.")


def label(d):
    """'S Cys81', 'O water', 'N HEM'."""
    if d["resname"] == "HOH":
        return "%s water" % d["el"]
    if d.get("polymer"):
        return "%s %s%s" % (d["el"], d["resname"].capitalize(), d["seqid"])
    return "%s %s" % (d["el"], d["resname"])


def geometry_prompt(metal, donors):
    vec = np.array([(d["dx"], d["dy"], d["dz"]) for d in donors], dtype=float)
    unit = vec / np.linalg.norm(vec, axis=1)[:, None]
    names = [label(d) for d in donors]
    seen = {}
    for i, n in enumerate(names):                       # two waters -> 'O water 1', 'O water 2'
        if names.count(n) > 1:
            seen[n] = seen.get(n, 0) + 1
            names[i] = "%s %d" % (n, seen[n])
    dist = "; ".join("%s %.2f Å" % (n, d["dist"]) for n, d in zip(names, donors))
    ang = "; ".join("%s / %s %.0f" % (names[i], names[j], np.degrees(np.arccos(np.clip(unit[i] @ unit[j], -1, 1))))
                    for i, j in itertools.combinations(range(len(donors)), 2))
    return ("Metal site in a protein structure. Metal: %s. Donor atoms (distance to the metal): %s. "
            "Donor-metal-donor angles (degrees): %s. %s" % (metal, dist, ang, INSTRUCTION))


def enzyme_prompt(desc, organism, metal):
    org = " from %s" % organism if organism else ""
    return ("Enzyme: %s%s. Which first coordination shell does its %s site have, and which single-atom catalyst "
            "sites does it suggest? %s" % (desc, org, metal, INSTRUCTION))


def answer(r):
    def fit(sup):
        k, score = r.get(sup + "_k"), r.get(sup + "_score")
        return None if pd.isna(k) or k != 3 or pd.isna(score) else bool(score <= FIT[sup])
    val = lambda x: None if pd.isna(x) else x       # noqa: E731
    return {"first_shell": r["raw_motif"], "graphene": val(r["graphene_motif"]), "MoS2": val(r.get("MoS2_motif")),
            "MoS2_fit": fit("MoS2"), "MoS2_vacancy": val(r.get("MoS2vS_motif")), "silica": val(r.get("SiO2_motif")),
            "silica_fit": fit("SiO2")}


def split_of(cluster):
    return "test" if int(hashlib.md5(cluster.encode()).hexdigest(), 16) % TEST_EVERY == 0 else "train"


def capped(df, cap, seed=7):
    return (df.sample(frac=1.0, random_state=seed).groupby("graphene_motif", dropna=False, group_keys=False)
            .head(cap).sort_values("site_id"))


def build():
    cols = ["site_id", "entry_id", "metal", "is_representative", "resolution", "cn", "raw_motif", "graphene_motif",
            "description", "organism", "cluster30", "host_entity", "reaction", "shell_json"]
    s = pd.read_parquet(os.path.join(COMPLETE, "sites.parquet"), columns=cols)
    s = s[s.is_representative & s.metal.isin(METALS) & (s.resolution <= MAX_RESOLUTION) & (s.cn >= 2)
          & s.description.fillna("").ne("")]
    s = s.merge(pd.read_parquet(os.path.join(COMPLETE, "supports.parquet")), on="site_id", how="left")
    s["cluster"] = s.cluster30.astype("string").fillna(s.entry_id + "_" + s.host_entity.astype("string").fillna(""))
    s["split"] = s.cluster.map(split_of)
    rows = []
    for split in ("train", "test"):
        part = s[s.split == split]
        for r in capped(part, CAP[("geometry", split)]).to_dict("records"):
            donors = json.loads(r["shell_json"])["donors"]
            rows.append({"split": split, "task": "geometry", "site_id": r["site_id"], "cluster": r["cluster"],
                         "prompt": geometry_prompt(r["metal"], donors), "completion": json.dumps(answer(r))})
        # enzyme task: (description, metal) pairs whose sites agree on the first shell and graphene site
        agree = part.groupby(["description", "metal"]).filter(
            lambda g: g.raw_motif.nunique() == 1 and g.graphene_motif.nunique(dropna=False) == 1)
        best = agree.sort_values("resolution").drop_duplicates(["description", "metal"])
        for r in capped(best, CAP[("enzyme", split)]).to_dict("records"):
            org = r["organism"] if isinstance(r["organism"], str) else ""
            rows.append({"split": split, "task": "enzyme", "site_id": r["site_id"], "cluster": r["cluster"],
                         "prompt": enzyme_prompt(r["description"], org, r["metal"]),
                         "completion": json.dumps(answer(r))})
    df = pd.DataFrame(rows)
    assert not set(df[df.split == "train"].cluster) & set(df[df.split == "test"].cluster)
    os.makedirs(OUT, exist_ok=True)
    for split in ("train", "test"):
        part = df[df.split == split].sample(frac=1.0, random_state=11).reset_index(drop=True)
        part["id"] = ["%s-%05d" % (split, i) for i in range(len(part))]
        part[["id", "prompt", "completion"]].to_json(os.path.join(OUT, split + ".jsonl"), orient="records",
                                                     lines=True, force_ascii=False)
        part[["id", "task", "site_id", "cluster"]].to_csv(os.path.join(OUT, split + "_meta.csv"), index=False)
    counts = df.groupby(["split", "task"]).size().unstack(fill_value=0)
    write_card(counts, df)
    print(counts.to_string())
    print("clusters: train %d, test %d" % (df[df.split == "train"].cluster.nunique(), df[df.split == "test"].cluster.nunique()))
    return df


def write_card(counts, df):
    ex = df[(df.split == "train") & (df.task == "geometry") & df.prompt.str.contains("Ni. ")].head(1)
    lines = [
        "# Enzyme metal sites as single-atom catalyst templates: fine-tuning data", "",
        "Prompt/completion pairs (JSONL, columns `id`, `prompt`, `completion`) built from every transition-metal "
        "site of the Protein Data Bank (RCSB data version 2026_41). Draft data from a hackathon; generated by "
        "`ml/adaption_dataset.py`, labels from the project's deterministic pipeline (`bio/first_shell.py`, "
        "`bio/supports.py`), not from a language model.", "",
        "## Tasks", "",
        "- **geometry**: metal, donor atoms with distances and donor-metal-donor angles -> first shell and its "
        "single-atom sites on graphene, MoS2 and silica. Needs the angles: e.g. four nearly square-planar donors go "
        "into a double vacancy, otherwise the best three-donor face into a single vacancy.",
        "- **enzyme**: enzyme name, organism and metal -> the same answer, for enzymes whose sites of that metal all "
        "share one first shell.", "",
        "## Answer fields", "",
        "| field | meaning |", "| --- | --- |",
        "| first_shell | metal and donor atoms, e.g. Ni-S4 |",
        "| graphene | SV/DV site, e.g. Ni-S3(SV) or Fe-N4(DV)+ax:N |",
        "| MoS2 | metal on a face of lattice S (+ax: donors that would need a thiol), e.g. Ni-S3(MoS2)+ax:S |",
        "| MoS2_fit | the site's S...S spacing within 0.25 A (RMS) of MoS2's 3.16 A (null without a three-S face) |",
        "| MoS2_vacancy | metal in an S vacancy on three Mo: only metal-cluster corners map, e.g. Fe-Mo3(MoS2vS)+ax:S |",
        "| silica | metal on siloxide O, e.g. Co-O3(SiO2)+ax:N3 |",
        "| silica_fit | O...O spacing within the 2.6-3.2 A window (null without a three-O face) |", "",
        "## Size and split", "",
        "Split by 30 % sequence-identity cluster: no protein family of the test set is in training. Capped per "
        "graphene motif (geometry: 120 train / 12 test; enzyme: 80 / 8) so the common motifs (Fe4S4 cubes, zinc "
        "fingers) do not dominate.", "",
        "| split | " + " | ".join(counts.columns) + " |", "| --- |" + " --- |" * len(counts.columns),
        *["| %s | %s |" % (i, " | ".join(str(v) for v in r)) for i, r in counts.iterrows()], "",
        "## Scoring", "",
        "Exact match per field and for the whole answer: `python ml/adaption_dataset.py score predictions.jsonl` "
        "(one line per test item: `{\"id\": ..., \"prediction\": ...}`). Compare the base model with the "
        "fine-tuned one on the same test file.", "",
        "## Example (geometry)", "",
        "```", ex.prompt.iloc[0] if len(ex) else "", "```", "",
        "```json", ex.completion.iloc[0] if len(ex) else "", "```", "",
        "## Licence and source", "",
        "Source structures: wwPDB archive, CC0 1.0. Derived data: licence to be decided by the authors before "
        "publication (proposal CC BY 4.0). No personal data.", ""]
    with open(os.path.join(OUT, "README.md"), "w") as f:
        f.write("\n".join(lines))


def score(pred_path):
    ref = pd.read_json(os.path.join(OUT, "test.jsonl"), lines=True).set_index("id")
    meta = pd.read_csv(os.path.join(OUT, "test_meta.csv")).set_index("id")
    pred = pd.read_json(pred_path, lines=True).set_index("id")
    rows = []
    for i, r in ref.iterrows():
        truth = json.loads(r.completion)
        raw = pred.prediction.get(i)
        try:
            p = json.loads(raw) if isinstance(raw, str) else (raw or {})
        except json.JSONDecodeError:
            p = {}
        hit = {f: p.get(f) == truth.get(f) for f in FIELDS}
        rows.append(dict(id=i, task=meta.task[i], parsed=bool(p), all=all(hit.values()), **hit))
    s = pd.DataFrame(rows)
    out = s.groupby("task")[["parsed", "all"] + FIELDS].mean().round(3)
    out.loc["overall"] = s[["parsed", "all"] + FIELDS].mean().round(3)
    print("%d test items, %d with a prediction" % (len(ref), ref.index.isin(pred.index).sum()))
    print(out.T.to_string())
    return out


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else ""
    if mode == "build":
        build()
    elif mode == "score" and len(sys.argv) > 2:
        score(sys.argv[2])
    else:
        sys.exit(__doc__)
