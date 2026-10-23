"""Search the RCSB PDB for transition-metal sites relevant to HER/ORR/CO2RR/NRR and download the structures.

Usage:
    python bio/pdb_sweep.py enzymes            # EC-class search      -> data/hits_enzymes.csv
    python bio/pdb_sweep.py all                # metal-component sweep -> data/hits_all.csv
    python bio/pdb_sweep.py download [IDS...]  # mmCIF of all hits (+ extra IDs) -> data/pdb/, data/entities.csv
    python bio/pdb_sweep.py complete [STEPS]   # every entry with a d-block metal component, no filter or grouping
                                               # -> data/complete/entries.parquet, entities.parquet, data/pdb/
                                               # (steps: search meta download tables; default all, each resumable)

`enzymes`: one RCSB Search API v2 query per EC class (rcsb_ec_lineage, so sub-classes count),
resolution <= 2.5 A, polymer entities grouped at 95 % sequence identity, best resolution as representative.
`all`: entries with a bound chemical component that contains Fe/Ni/Co/Cu/Mo/W/Mn/V, resolution <= 2.0 A,
protein entities grouped at 30 % identity. The metal-containing component ids come from the wwPDB
Chemical Component Dictionary, because the chemical formula search misses the Fe-S clusters (SF4, FES,
F3S), FeMo-co (ICS) and P-cluster (CLF); they are matched via
rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id, since rcsb_entry_info.nonpolymer_bound_components
is not searchable. Every HTTP response is cached in data/cache/, so reruns are free.
`complete`: the same component match for all d-block metals (Sc-Zn, Y-Cd, Hf-Hg), every experimental entry, no
resolution or method filter, no grouping; GraphQL metadata per entry and polymer entity, RCSB sequence clusters at
30/50/95/100 % and all structures. Its caches are in data/cache/complete/ (the older modes' caches are not touched),
query dates, data version and counts in data/complete/provenance.json.
"""
import fcntl
import gzip
import hashlib
import json
import re
import subprocess
import sys
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
PDB_DIR = DATA / "pdb"
CACHE = DATA / "cache"
SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
GRAPHQL_URL = "https://data.rcsb.org/graphql"
FILE_URL = "https://files.rcsb.org/download/{}.cif.gz"
CCD_URL = "https://files.wwpdb.org/pub/pdb/data/monomers/components.cif.gz"
CLUSTER_URL = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt"

MAX_WORKERS = 8          # concurrent downloads (be polite to RCSB)
PAGE_ROWS = 1000         # search pagination
META_BATCH = 100         # entries per GraphQL request
RES_ENZYMES, ID_ENZYMES = 2.5, 95
RES_ALL, ID_ALL = 2.0, 30
MAX_ALL = 6000
SWEEP_METALS = ["Fe", "Ni", "Co", "Cu", "Mo", "W", "Mn", "V"]

# reaction -> [(EC lineage id, enzyme label)]
EC_CLASSES = {
    "HER": [("1.12", "hydrogenase")],
    "ORR": [("1.10.3.2", "laccase"), ("1.10.3.3", "ascorbate oxidase"), ("1.3.3.5", "bilirubin oxidase"),
            ("1.16.3.1", "ferroxidase"), ("7.1.1.9", "cytochrome c oxidase"), ("1.9.3.1", "cytochrome c oxidase"),
            ("7.1.1.7", "cytochrome bd oxidase"), ("7.1.1.3", "bo3 quinol oxidase"),
            ("1.10.3.11", "alternative oxidase")],
    "CO2RR": [("1.2.7.4", "Ni-Fe CO dehydrogenase"), ("1.2.5.3", "Mo-Cu CO dehydrogenase"),
              ("1.17.1.9", "formate dehydrogenase"), ("1.17.1.10", "formate dehydrogenase"),
              ("1.17.1.11", "formate dehydrogenase"), ("1.17.5.3", "formate dehydrogenase"),
              ("1.17.2.3", "formate dehydrogenase"), ("1.17.98.3", "formate dehydrogenase"),
              ("1.17.99.7", "formate dehydrogenase"), ("1.2.1.2", "formate dehydrogenase"),
              ("2.3.1.169", "acetyl-CoA synthase"), ("1.2.7.12", "formylmethanofuran dehydrogenase"),
              ("2.8.4.1", "methyl-coenzyme M reductase"), ("4.2.1.1", "carbonic anhydrase")],
    "NRR": [("1.18.6.1", "nitrogenase"), ("1.19.6.1", "nitrogenase")],
}

META_QUERY = """query($ids: [String!]!) { entries(entry_ids: $ids) {
  rcsb_id rcsb_entry_info { resolution_combined experimental_method } struct { title }
  polymer_entities { rcsb_id
    rcsb_polymer_entity_container_identifiers { entity_id auth_asym_ids asym_ids }
    rcsb_polymer_entity { pdbx_description pdbx_ec rcsb_ec_lineage { id } }
    entity_poly { rcsb_entity_polymer_type } rcsb_entity_source_organism { scientific_name } } } }"""


def make_session():
    s = requests.Session()
    retry = Retry(total=6, backoff_factor=2, status_forcelist=(429, 500, 502, 503, 504), allowed_methods=None)
    s.mount("https://", HTTPAdapter(max_retries=retry, pool_maxsize=MAX_WORKERS))
    s.headers["User-Agent"] = "sac-motif-sweep/0.1 (FU Berlin hackathon)"
    return s


SESSION = make_session()


def text_node(attribute, operator, value):
    return {"type": "terminal", "service": "text",
            "parameters": {"attribute": attribute, "operator": operator, "value": value}}


def rcsb_search(nodes, identity):
    """AND-combine nodes, group polymer entities at `identity` %, return representatives (all pages)."""
    query = {"query": {"type": "group", "logical_operator": "and", "nodes": nodes},
             "return_type": "polymer_entity",
             "request_options": {
                 "group_by": {"aggregation_method": "sequence_identity", "similarity_cutoff": identity,
                              "ranking_criteria_type": {"sort_by": "rcsb_entry_info.resolution_combined",
                                                        "direction": "asc"}},
                 "group_by_return_type": "representatives",
                 "results_content_type": ["experimental"]}}
    key = hashlib.sha1(json.dumps(query, sort_keys=True).encode()).hexdigest()[:16]
    cached = CACHE / "search" / f"{key}.json"
    if cached.exists():
        return json.loads(cached.read_text())
    ids, total, groups, start = [], 0, 0, 0
    while True:
        query["request_options"]["paginate"] = {"start": start, "rows": PAGE_ROWS}
        r = SESSION.post(SEARCH_URL, json=query, timeout=300)
        if r.status_code == 204 or not r.text:  # 204 = no hits
            break
        r.raise_for_status()
        d = r.json()
        total, groups = d.get("total_count", 0), d.get("group_by_count", 0)
        ids += [x["identifier"] for x in d["result_set"]]
        start += PAGE_ROWS
        if start >= groups or not d["result_set"]:
            break
    res = {"entities_total": total, "groups": groups, "ids": ids}
    cached.parent.mkdir(parents=True, exist_ok=True)
    cached.write_text(json.dumps(res))
    return res


def entry_meta(entry_ids):
    """GraphQL metadata per entry (resolution, title, polymer entities with EC/description/organism)."""
    f = CACHE / "entry_meta.json"
    cache = json.loads(f.read_text()) if f.exists() else {}
    todo = sorted(set(entry_ids) - set(cache))
    for i in range(0, len(todo), META_BATCH):
        batch = todo[i:i + META_BATCH]
        r = SESSION.post(GRAPHQL_URL, json={"query": META_QUERY, "variables": {"ids": batch}}, timeout=300)
        r.raise_for_status()
        for e in r.json()["data"]["entries"] or []:
            cache[e["rcsb_id"]] = e
        for eid in batch:  # obsolete/unknown ids: remember the miss
            cache.setdefault(eid, None)
        f.write_text(json.dumps(cache))
        print(f"  metadata {min(i + META_BATCH, len(todo))}/{len(todo)}", flush=True)
    return {e: cache.get(e) for e in entry_ids}


def resolution(meta):
    res = (meta or {}).get("rcsb_entry_info", {}).get("resolution_combined") or []
    return min(res) if res else float("nan")


def ec_reactions(lineage_ids):
    """Our (reaction, EC class, label) tuples matched by an entity's EC lineage."""
    return [(rxn, ec, label) for rxn, ecs in EC_CLASSES.items() for ec, label in ecs if ec in lineage_ids]


def entity_rows(meta):
    """One dict per polymer entity of an entry."""
    rows = []
    for pe in (meta or {}).get("polymer_entities") or []:
        ids = pe["rcsb_polymer_entity_container_identifiers"]
        ent = pe.get("rcsb_polymer_entity") or {}
        lineage = [x["id"] for x in ent.get("rcsb_ec_lineage") or []]
        matched = ec_reactions(lineage)
        rows.append({
            "entry_id": meta["rcsb_id"], "entity_id": ids["entity_id"], "rcsb_entity": pe["rcsb_id"],
            "auth_asym_ids": " ".join(ids.get("auth_asym_ids") or []),
            "asym_ids": " ".join(ids.get("asym_ids") or []),
            "polymer_type": (pe.get("entity_poly") or {}).get("rcsb_entity_polymer_type", ""),
            "description": ent.get("pdbx_description") or "",
            "ec": ent.get("pdbx_ec") or "",
            "reaction": ";".join(sorted({m[0] for m in matched})),
            "query_ec": ";".join(sorted({m[1] for m in matched})),
            "enzyme": ";".join(sorted({m[2] for m in matched})),
            "organism": ";".join(sorted({o["scientific_name"] for o in pe.get("rcsb_entity_source_organism") or []
                                        if o.get("scientific_name")})),
        })
    return rows


def mode_enzymes():
    counts, hits = {}, defaultdict(set)  # rcsb_entity -> {(rxn, ec, label)}
    for rxn, ecs in EC_CLASSES.items():
        for ec, label in ecs:
            res = rcsb_search([text_node("rcsb_polymer_entity.rcsb_ec_lineage.id", "exact_match", ec),
                               text_node("rcsb_entry_info.resolution_combined", "less_or_equal", RES_ENZYMES)],
                              ID_ENZYMES)
            counts[ec] = {"reaction": rxn, "enzyme": label, "entities_total": res["entities_total"],
                          "representatives": len(res["ids"]),
                          "entries": len({i.split("_")[0] for i in res["ids"]})}
            for i in res["ids"]:
                hits[i].add((rxn, ec, label))
            print(f"{rxn:6s} {ec:10s} {label:34s} entities={res['entities_total']:6d} "
                  f"representatives={len(res['ids']):5d}", flush=True)
    empty = [ec for ec, c in counts.items() if c["representatives"] == 0]
    for ec in empty:  # tell "no structure at all" apart from "only low-resolution structures"
        n = rcsb_search([text_node("rcsb_polymer_entity.rcsb_ec_lineage.id", "exact_match", ec)], ID_ENZYMES)
        counts[ec]["entities_any_resolution"] = n["entities_total"]
        print(f"EC {ec} has no hit at <= {RES_ENZYMES} A ({n['entities_total']} entities at any resolution)")

    meta = entry_meta(sorted({i.split("_")[0] for i in hits}))
    rows = []
    for rid, matched in hits.items():
        entry, ent = rid.split("_")
        m = meta.get(entry)
        info = next((r for r in entity_rows(m) if r["entity_id"] == ent), {})
        rows.append({"entry_id": entry, "entity_id": ent, "reaction": ";".join(sorted({x[0] for x in matched})),
                     "query_ec": ";".join(sorted({x[1] for x in matched})),
                     "enzyme": ";".join(sorted({x[2] for x in matched})), "ec": info.get("ec", ""),
                     "description": info.get("description", ""), "organism": info.get("organism", ""),
                     "resolution": resolution(m),
                     "method": (m or {}).get("rcsb_entry_info", {}).get("experimental_method", ""),
                     "title": ((m or {}).get("struct") or {}).get("title", "")})
    df = pd.DataFrame(rows).sort_values(["reaction", "enzyme", "resolution", "entry_id", "entity_id"])
    df.to_csv(DATA / "hits_enzymes.csv", index=False)
    print(f"{len(df)} representative entities in {df.entry_id.nunique()} entries -> data/hits_enzymes.csv")
    print(df.groupby("reaction").agg(entities=("entity_id", "size"), entries=("entry_id", "nunique")).to_string())
    save_counts("enzymes", {"per_ec": counts, "empty_ec": empty, "entities": len(df),
                            "entries": int(df.entry_id.nunique())})


def metal_components():
    """{metal: [CCD comp ids whose formula contains that element]} from the Chemical Component Dictionary."""
    return ccd_components(SWEEP_METALS, CACHE / "metal_components.json")


def ccd_components(metals, out_file):
    """{metal: [CCD comp ids whose formula contains that element]} for `metals`, cached in out_file."""
    if out_file.exists():
        return json.loads(out_file.read_text())
    ccd = CACHE / "components.cif.gz"
    if not ccd.exists():
        download(CCD_URL, ccd)
    token = re.compile(r"([A-Z][a-z]?)\d*")
    comps, cid = defaultdict(list), None
    with gzip.open(ccd, "rt") as fh:  # stream: the dictionary is ~0.5 GB uncompressed
        for line in fh:
            if line.startswith("_chem_comp.id "):
                cid = line.split()[1]
            elif line.startswith("_chem_comp.formula "):
                elements = set(token.findall(line.split(None, 1)[1]))
                for metal in metals:
                    if metal in elements:
                        comps[metal].append(cid)
    out_file.write_text(json.dumps(comps))
    return comps


def mode_all():
    comps = metal_components()
    reps, counts = defaultdict(set), {}  # rcsb_entity -> {metal}
    for metal in SWEEP_METALS:
        res = rcsb_search([text_node("rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id", "in",
                                     comps[metal]),
                           text_node("rcsb_entry_info.resolution_combined", "less_or_equal", RES_ALL),
                           text_node("entity_poly.rcsb_entity_polymer_type", "exact_match", "Protein")], ID_ALL)
        for i in res["ids"]:
            reps[i].add(metal)
        counts[metal] = {"components": len(comps[metal]), "protein_entities": res["entities_total"],
                         "clusters_30": len(res["ids"]), "entries": len({i.split("_")[0] for i in res["ids"]})}
        print(f"{metal:3s} components={len(comps[metal]):4d} protein entities={res['entities_total']:6d} "
              f"30%-clusters={len(res['ids']):5d} representative entries={counts[metal]['entries']:5d}", flush=True)
    by_entry = defaultdict(lambda: {"entities": set(), "metals": set()})
    for rid, metals in reps.items():
        by_entry[rid.split("_")[0]]["entities"].add(rid.split("_")[1])
        by_entry[rid.split("_")[0]]["metals"] |= metals
    meta = entry_meta(sorted(by_entry))
    rows = [{"entry_id": e, "rep_entities": " ".join(sorted(v["entities"], key=int)),
             "metals": " ".join(m for m in SWEEP_METALS if m in v["metals"]), "resolution": resolution(meta[e]),
             "method": (meta[e] or {}).get("rcsb_entry_info", {}).get("experimental_method", ""),
             "title": ((meta[e] or {}).get("struct") or {}).get("title", "")} for e, v in by_entry.items()]
    df = pd.DataFrame(rows).sort_values(["resolution", "entry_id"]).reset_index(drop=True)
    df["selected"] = df.index < MAX_ALL
    df.to_csv(DATA / "hits_all.csv", index=False)
    print(f"{len(df)} unique representative entries over all metals"
          + (f"; > {MAX_ALL}, kept the {MAX_ALL} best-resolution ones (cut at {df.resolution[MAX_ALL - 1]:.2f} A)"
             if len(df) > MAX_ALL else "; all selected") + " -> data/hits_all.csv")
    save_counts("all", {"per_metal": counts, "entries": len(df), "selected": int(df.selected.sum())})


def save_counts(key, value):
    f = DATA / "sweep_counts.json"
    d = json.loads(f.read_text()) if f.exists() else {}
    d[key] = value
    f.write_text(json.dumps(d, indent=1))


def download(url, dest):
    """GET url -> dest (atomic rename), retried with backoff by the session."""
    r = SESSION.get(url, timeout=300)
    if r.status_code == 404:
        return False
    r.raise_for_status()
    tmp = dest.with_name(dest.name + ".part")
    tmp.write_bytes(r.content)
    tmp.rename(dest)
    return True


def fetch_structures(entry_ids):
    """Download <ID>.cif.gz for every id not yet in data/pdb/; returns the ids that failed."""
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    todo = [e for e in sorted(set(entry_ids)) if not (PDB_DIR / f"{e}.cif.gz").exists()]
    failed = []

    def get(e):
        try:
            if not download(FILE_URL.format(e), PDB_DIR / f"{e}.cif.gz"):
                failed.append(e)
        except requests.RequestException as exc:
            print(f"  {e}: {exc}", flush=True)
            failed.append(e)

    with ThreadPoolExecutor(MAX_WORKERS) as pool:
        for n, _ in enumerate(pool.map(get, todo), 1):
            if n % 250 == 0:
                print(f"  downloaded {n}/{len(todo)}", flush=True)
    print(f"{len(set(entry_ids))} entries requested, {len(todo)} fetched now, {len(failed)} failed "
          f"{' '.join(failed)}", flush=True)
    return failed


def cluster30():
    """{rcsb_entity: 30 %-identity cluster index} from the RCSB sequence cluster file."""
    f = CACHE / "clusters-by-entity-30.txt"
    if not f.exists():
        download(CLUSTER_URL, f)
    return {ent: i for i, line in enumerate(f.read_text().splitlines()) for ent in line.split()}


def mode_download(extra):
    ids = {e.upper() for e in extra}
    for name in ["hits_enzymes.csv", "hits_all.csv"]:
        if (DATA / name).exists():
            df = pd.read_csv(DATA / name, dtype={"entry_id": str})  # ids like 1E10 are not floats
            ids |= set(df[df.selected].entry_id if "selected" in df else df.entry_id)
    fetch_structures(ids)
    # data/pdb/ also holds the structures of the complete d-block sweep (data/complete/entities.parquet): list only
    # this sweep's entries, i.e. hits, extra ids and entries listed before (e.g. the sanity-check structures)
    if (DATA / "entities.csv").exists():
        ids |= set(pd.read_csv(DATA / "entities.csv", dtype={"entry_id": str}).entry_id)
    write_entities(sorted(e for e in ids if (PDB_DIR / f"{e}.cif.gz").exists()))


def write_entities(entry_ids):
    """data/entities.csv: polymer entities of every downloaded entry (labels + 30 % sequence cluster)."""
    meta, clusters = entry_meta(entry_ids), cluster30()
    rows = [r for e in entry_ids for r in entity_rows(meta[e])]
    for r in rows:  # entities missing from the cluster file (very new) form their own cluster
        r["cluster30"] = clusters.get(r["rcsb_entity"], r["rcsb_entity"])
        r["resolution"] = resolution(meta[r["entry_id"]])
    pd.DataFrame(rows).to_csv(DATA / "entities.csv", index=False)
    print(f"{len(rows)} polymer entities of {len(entry_ids)} entries -> data/entities.csv")


# ---------------------------------------------------------------- complete d-block sweep (dataset release)

DBLOCK = ["Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd",
          "Ag", "Cd", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg"]
COMPLETE = DATA / "complete"
CCACHE = CACHE / "complete"  # caches of the complete sweep, apart from those the older modes read
CLUSTER_LEVELS = (30, 50, 95, 100)
CLUSTER_FILE_URL = "https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-{}.txt"
HOLDINGS_URL = "https://data.rcsb.org/rest/v1/holdings/current/entry_ids"
STATUS_URL = "https://data.rcsb.org/rest/v1/holdings/status/{}"
SCHEMA_URL = "https://data.rcsb.org/rest/v1/schema/entry"
COMPLETE_PAGE_ROWS = 10000
COMPLETE_STEPS = ("search", "meta", "download", "tables")

COMPLETE_QUERY = """query($ids: [String!]!) { entries(entry_ids: $ids) {
  rcsb_id struct { title } exptl { method }
  rcsb_entry_info { resolution_combined experimental_method deposited_atom_count polymer_entity_count
                    nonpolymer_entity_count }
  rcsb_accession_info { initial_release_date deposit_date revision_date major_revision minor_revision }
  refine { pdbx_refine_id ls_R_factor_R_free ls_R_factor_R_work }
  rcsb_primary_citation { pdbx_database_id_DOI pdbx_database_id_PubMed year rcsb_journal_abbrev journal_abbrev }
  polymer_entities { rcsb_id
    rcsb_polymer_entity_container_identifiers { entity_id auth_asym_ids asym_ids uniprot_ids
      reference_sequence_identifiers { database_name database_accession } }
    rcsb_polymer_entity { pdbx_description pdbx_ec rcsb_ec_lineage { id } }
    entity_poly { rcsb_entity_polymer_type }
    rcsb_entity_source_organism { scientific_name ncbi_taxonomy_id } }
  nonpolymer_entities { rcsb_nonpolymer_entity_container_identifiers { entity_id nonpolymer_comp_id } } } }"""


def utcnow():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def metal_components_dblock():
    """{metal: [CCD comp ids]} for all d-block metals (Sc-Zn, Y-Cd, Hf-Hg) from the Chemical Component Dictionary."""
    return ccd_components(DBLOCK, CACHE / "metal_components_dblock.json")


def update_provenance(section, value):
    """Set one section of data/complete/provenance.json (file lock: several steps may run in parallel)."""
    COMPLETE.mkdir(parents=True, exist_ok=True)
    f = COMPLETE / "provenance.json"
    with open(COMPLETE / ".provenance.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        d = json.loads(f.read_text()) if f.exists() else {}
        d[section] = value
        tmp = f.with_name(f.name + ".tmp")
        tmp.write_text(json.dumps(d, indent=1, default=str))
        tmp.rename(f)


def code_version():
    """Git commit of the code, whether bio/ has uncommitted changes, and the sha256 of both scripts."""
    def git(*args):
        try:
            return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True,
                                  timeout=30).stdout.strip()
        except (OSError, subprocess.SubprocessError):
            return ""
    return {"git_commit": git("rev-parse", "HEAD"),
            "uncommitted_changes_in_bio": git("status", "--porcelain", "--", "bio").splitlines(),
            "sha256": {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in sorted((ROOT / "bio").glob("*.py"))}}


def software_versions():
    import platform
    out = {"python": platform.python_version(), "pandas": pd.__version__, "requests": requests.__version__}
    for mod in ("gemmi", "numpy", "pyarrow", "scipy"):
        try:
            out[mod] = __import__(mod).__version__
        except ImportError:
            pass
    return out


def rcsb_data_version():
    """What the RCSB APIs tell about the data state: weekly update id, newest release date, holdings count."""
    newest = {"query": text_node("rcsb_accession_info.initial_release_date", "exists", None),
              "return_type": "entry",
              "request_options": {"paginate": {"start": 0, "rows": 1}, "results_content_type": ["experimental"],
                                  "sort": [{"sort_by": "rcsb_accession_info.initial_release_date",
                                            "direction": "desc"}]}}
    del newest["query"]["parameters"]["value"]
    d = SESSION.post(SEARCH_URL, json=newest, timeout=120).json()
    entry = d["result_set"][0]["identifier"]
    q = '{ entry(entry_id: "%s") { rcsb_accession_info { initial_release_date } } }' % entry
    released = SESSION.post(GRAPHQL_URL, json={"query": q}, timeout=120).json()["data"]["entry"]
    status = SESSION.get(STATUS_URL.format(entry), timeout=120).json()
    return {"checked_utc": utcnow(),
            "update_id": status["rcsb_repository_holdings_combined_entry_container_identifiers"]["update_id"],
            "newest_initial_release_date": released["rcsb_accession_info"]["initial_release_date"][:10],
            "experimental_entries": d["total_count"],
            "current_holdings": len(SESSION.get(HOLDINGS_URL, timeout=300).json())}


def search_all_entries(comp_ids, cache_file):
    """Every experimental entry with a non-polymer entity whose component is in comp_ids: return_type entry,
    no grouping, no resolution/method filter, sorted by id, all pages; cached with the query time."""
    if cache_file.exists():
        return json.loads(cache_file.read_text())
    query = {"query": text_node("rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id", "in",
                                sorted(comp_ids)),
             "return_type": "entry",
             "request_options": {"results_content_type": ["experimental"],
                                 "sort": [{"sort_by": "rcsb_entry_container_identifiers.entry_id",
                                           "direction": "asc"}]}}
    started, ids, total, start = utcnow(), [], 0, 0
    while True:
        query["request_options"]["paginate"] = {"start": start, "rows": COMPLETE_PAGE_ROWS}
        r = SESSION.post(SEARCH_URL, json=query, timeout=300)
        if r.status_code == 204 or not r.text:  # 204 = no hits
            break
        r.raise_for_status()
        d = r.json()
        total = d["total_count"]
        ids += [x["identifier"] for x in d["result_set"]]
        start += COMPLETE_PAGE_ROWS
        if start >= total or not d["result_set"]:
            break
    if len(set(ids)) != total:
        raise RuntimeError(f"pagination returned {len(set(ids))} distinct ids for total_count {total}; rerun")
    del query["request_options"]["paginate"]
    res = {"url": SEARCH_URL, "query": query, "started_utc": started, "finished_utc": utcnow(),
           "total_count": total, "ids": sorted(set(ids))}
    cache_file.parent.mkdir(parents=True, exist_ok=True)
    cache_file.write_text(json.dumps(res))
    return res


def complete_search():
    """Per-metal entry lists -> {entry: [metals]} (d-block order)."""
    comps = metal_components_dblock()
    vfile = CCACHE / "search" / "rcsb_data_version.json"  # data state of the cached searches, not of today
    cached = vfile.exists() and all((CCACHE / "search" / f"{m}.json").exists() for m in DBLOCK if comps.get(m))
    before = json.loads(vfile.read_text())["before"] if cached else rcsb_data_version()
    per_metal = {}
    for metal in DBLOCK:
        if not comps.get(metal):
            per_metal[metal] = {"components": 0, "entries": 0}
            print(f"{metal:3s} no CCD component", flush=True)
            continue
        res = search_all_entries(comps[metal], CCACHE / "search" / f"{metal}.json")
        per_metal[metal] = {"components": len(comps[metal]), "entries": res["total_count"],
                            "query_utc": res["started_utc"]}
        print(f"{metal:3s} components={len(comps[metal]):4d} entries={res['total_count']:6d}", flush=True)
    if cached:
        after = json.loads(vfile.read_text())["after"]
    else:
        after = rcsb_data_version()
        vfile.write_text(json.dumps({"before": before, "after": after}, indent=1))
    if (before["update_id"], before["newest_initial_release_date"]) != (after["update_id"],
                                                                        after["newest_initial_release_date"]):
        print(f"WARNING: the RCSB data changed during the search ({before} -> {after})", flush=True)
    metals = entry_metals()
    print(f"{len(metals)} distinct entries over {sum(v['entries'] > 0 for v in per_metal.values())} metals",
          flush=True)
    update_provenance("search", {
        "endpoint": SEARCH_URL, "return_type": "entry", "grouping": "none",
        "results_content_type": ["experimental"], "filters": "none (no resolution or method filter)",
        "attribute": "rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id (operator 'in')",
        "component_source": f"wwPDB Chemical Component Dictionary ({CCD_URL}): every component whose "
                            "_chem_comp.formula contains the element; cached in "
                            "data/cache/metal_components_dblock.json",
        "ccd_file_mtime_utc": datetime.fromtimestamp((CACHE / "components.cif.gz").stat().st_mtime,
                                                     timezone.utc).isoformat(timespec="seconds")
        if (CACHE / "components.cif.gz").exists() else None,
        "metals": DBLOCK, "rcsb_data_version_before": before, "rcsb_data_version_after": after,
        "per_metal": per_metal, "entries_found": len(metals), "components_per_metal": comps})
    return metals


def entry_metals():
    """{entry: [metals]} from the cached per-metal searches."""
    metals = defaultdict(list)
    for metal in DBLOCK:
        f = CCACHE / "search" / f"{metal}.json"
        if f.exists():
            for e in json.loads(f.read_text())["ids"]:
                metals[e].append(metal)
    return dict(sorted(metals.items()))


def complete_meta(entry_ids):
    """GraphQL metadata per entry (batches of META_BATCH), appended to a JSON-lines cache: resumable."""
    f = CCACHE / "entry_meta.jsonl"
    f.parent.mkdir(parents=True, exist_ok=True)
    cache = {}
    if f.exists():
        for line in f.read_text().splitlines():
            try:
                d = json.loads(line)
            except json.JSONDecodeError:  # torn last line of an interrupted run
                continue
            cache[d["id"]] = d["meta"]
    todo = [e for e in sorted(set(entry_ids)) if e not in cache]
    print(f"metadata: {len(cache)} cached, {len(todo)} to fetch", flush=True)
    with f.open("a") as out:
        for i in range(0, len(todo), META_BATCH):
            batch = todo[i:i + META_BATCH]
            got = graphql_entries(batch)
            for eid in batch:  # obsolete/unknown ids are stored as null
                cache[eid] = got.get(eid)
                out.write(json.dumps({"id": eid, "fetched_utc": utcnow(), "meta": cache[eid]}) + "\n")
            out.flush()
            if (i // META_BATCH) % 20 == 0 or i + META_BATCH >= len(todo):
                print(f"  metadata {min(i + META_BATCH, len(todo))}/{len(todo)}", flush=True)
    return {e: cache.get(e) for e in entry_ids}


def graphql_entries(batch):
    """{entry: metadata} for a batch; a failing batch is split in halves (one bad entry must not block the rest)."""
    try:
        r = SESSION.post(GRAPHQL_URL, json={"query": COMPLETE_QUERY, "variables": {"ids": batch}}, timeout=600)
        r.raise_for_status()
        d = r.json()
        if d.get("errors"):
            raise RuntimeError(str(d["errors"])[:500])
        return {e["rcsb_id"]: e for e in d["data"]["entries"] or []}
    except (requests.RequestException, RuntimeError, ValueError) as exc:
        if len(batch) == 1:
            print(f"  metadata {batch[0]} failed: {exc}", flush=True)
            return {}
        print(f"  metadata batch of {len(batch)} failed ({type(exc).__name__}), splitting", flush=True)
        time.sleep(10)
        half = len(batch) // 2
        return {**graphql_entries(batch[:half]), **graphql_entries(batch[half:])}


def cluster_maps(entities):
    """{level: {rcsb_entity: 0-based line index of the RCSB cluster file}} for the given entities; the files
    are cached in data/cache/complete/ together with their download time and Last-Modified header."""
    info_file = CCACHE / "cluster_files.json"
    info = json.loads(info_file.read_text()) if info_file.exists() else {}
    wanted, out = set(entities), {}
    for n in CLUSTER_LEVELS:
        f = CCACHE / f"clusters-by-entity-{n}.txt"
        if not f.exists():
            r = SESSION.get(CLUSTER_FILE_URL.format(n), timeout=600)
            r.raise_for_status()
            tmp = f.with_name(f.name + ".part")
            tmp.write_bytes(r.content)
            tmp.rename(f)
            info[str(n)] = {"url": CLUSTER_FILE_URL.format(n), "downloaded_utc": utcnow(),
                            "last_modified": r.headers.get("Last-Modified"),
                            "sha256": hashlib.sha256(r.content).hexdigest(), "bytes": len(r.content)}
            info_file.write_text(json.dumps(info, indent=1))
        lines = f.read_text().splitlines()
        out[n] = {ent: i for i, line in enumerate(lines) for ent in line.split() if ent in wanted}
        info.setdefault(str(n), {})["clusters"] = len(lines)
    info_file.write_text(json.dumps(info, indent=1))
    return out, info


def first_refine(refine, key):
    """Value of `key` from the X-ray refine record, else from the first record that has it."""
    recs = sorted(refine or [], key=lambda r: r.get("pdbx_refine_id") != "X-RAY DIFFRACTION")
    return next((r[key] for r in recs if r.get(key) is not None), None)


def complete_entity_rows(meta, clusters):
    rows = []
    for pe in sorted((meta or {}).get("polymer_entities") or [],
                     key=lambda p: int(p["rcsb_polymer_entity_container_identifiers"]["entity_id"])):
        ids = pe["rcsb_polymer_entity_container_identifiers"]
        ent = pe.get("rcsb_polymer_entity") or {}
        lineage = [x["id"] for x in ent.get("rcsb_ec_lineage") or []]
        matched = ec_reactions(lineage)
        orgs = list(dict.fromkeys((o.get("scientific_name") or "", str(o.get("ncbi_taxonomy_id") or ""))
                                  for o in pe.get("rcsb_entity_source_organism") or []))
        uniprot = ids.get("uniprot_ids") or [r["database_accession"] for r in
                                             ids.get("reference_sequence_identifiers") or []
                                             if r.get("database_name") == "UniProt"]
        row = {"entry_id": meta["rcsb_id"], "entity_id": ids["entity_id"], "rcsb_entity": pe["rcsb_id"],
               "auth_asym_ids": " ".join(ids.get("auth_asym_ids") or []),
               "asym_ids": " ".join(ids.get("asym_ids") or []),
               "polymer_type": (pe.get("entity_poly") or {}).get("rcsb_entity_polymer_type") or "",
               "description": ent.get("pdbx_description") or "", "ec": ent.get("pdbx_ec") or "",
               "ec_lineage": " ".join(lineage),
               "reaction": ";".join(sorted({m[0] for m in matched})),
               "query_ec": ";".join(sorted({m[1] for m in matched})),
               "enzyme": ";".join(sorted({m[2] for m in matched})),
               "organism": ";".join(o[0] for o in orgs), "taxonomy_id": ";".join(o[1] for o in orgs),
               "uniprot": " ".join(dict.fromkeys(uniprot))}
        for n in CLUSTER_LEVELS:
            row[f"cluster{n}"] = clusters[n].get(pe["rcsb_id"])
        rows.append(row)
    return rows


def complete_entry_row(entry, metals, meta, comps_by_id, status):
    m = meta or {}
    info, acc = m.get("rcsb_entry_info") or {}, m.get("rcsb_accession_info") or {}
    cit = m.get("rcsb_primary_citation") or {}
    comp_ids = sorted({x["rcsb_nonpolymer_entity_container_identifiers"]["nonpolymer_comp_id"]
                       for x in m.get("nonpolymer_entities") or []} & set(comps_by_id))
    res = info.get("resolution_combined") or []
    pubmed = cit.get("pdbx_database_id_PubMed")
    return {"entry_id": entry, "metals": " ".join(metals), "metal_comp_ids": " ".join(comp_ids),
            "title": (m.get("struct") or {}).get("title") or "",
            "method": info.get("experimental_method") or "",
            "exptl_methods": ";".join(x["method"] for x in m.get("exptl") or [] if x.get("method")),
            "resolution": min(res) if res else float("nan"),
            "resolution_combined": ";".join(f"{x:g}" for x in res),
            "r_free": first_refine(m.get("refine"), "ls_R_factor_R_free"),
            "r_work": first_refine(m.get("refine"), "ls_R_factor_R_work"),
            "deposit_date": (acc.get("deposit_date") or "")[:10],
            "initial_release_date": (acc.get("initial_release_date") or "")[:10],
            "revision_date": (acc.get("revision_date") or "")[:10],
            "revision": f"{acc['major_revision']}.{acc['minor_revision']}" if acc.get("major_revision") else "",
            "citation_doi": cit.get("pdbx_database_id_DOI") or "",
            "citation_pubmed": pubmed if pubmed and pubmed > 0 else None,  # the archive writes -1 for "none"
            "citation_year": cit.get("year"),
            "citation_journal": cit.get("rcsb_journal_abbrev") or cit.get("journal_abbrev") or "",
            "deposited_atom_count": info.get("deposited_atom_count"),
            "polymer_entity_count": info.get("polymer_entity_count"),
            "nonpolymer_entity_count": info.get("nonpolymer_entity_count"),
            "metadata_found": meta is not None, **status}


def read_provenance():
    f = COMPLETE / "provenance.json"
    return json.loads(f.read_text()) if f.exists() else {}


def complete_download(entry_ids):
    """fetch_structures for all entries; failed ids are classified via the RCSB holdings status."""
    t0, started = time.time(), utcnow()
    on_disk = sum((PDB_DIR / f"{e}.cif.gz").exists() for e in entry_ids)
    failed = fetch_structures(entry_ids)
    status = {}
    for e in failed:
        try:
            h = SESSION.get(STATUS_URL.format(e), timeout=120).json().get("rcsb_repository_holdings_combined", {})
            status[e] = f"{h.get('status_code', '?')}" + (f" replaced by {h['id_code_replaced_by_latest']}"
                                                          if h.get("id_code_replaced_by_latest") else "")
        except (requests.RequestException, ValueError) as exc:
            status[e] = f"status lookup failed: {exc}"
    (COMPLETE / "download_failures.json").write_text(json.dumps(status, indent=1))
    run = {"started_utc": started, "finished_utc": utcnow(), "seconds": round(time.time() - t0),
           "requested": len(set(entry_ids)), "already_on_disk": on_disk,
           "fetched": len(set(entry_ids)) - on_disk - len(failed), "failed": len(failed)}
    update_provenance("download", {
        "endpoint": FILE_URL, "workers": MAX_WORKERS,
        "runs": read_provenance().get("download", {}).get("runs", []) + [run], "failed": status,
        "note": "data/pdb/ is shared with the earlier sweep: files already there were not fetched again; the "
                "revision of every processed file is recorded per entry (entry_log.parquet)"})
    return status


def complete_tables(entry_metals_map):
    """data/complete/entries.parquet and entities.parquet."""
    meta = complete_meta(list(entry_metals_map))
    comps = metal_components_dblock()
    comp_metal = defaultdict(list)
    for metal, ids in comps.items():
        for c in ids:
            comp_metal[c].append(metal)
    failures_file = COMPLETE / "download_failures.json"
    failures = json.loads(failures_file.read_text()) if failures_file.exists() else {}
    wanted = [pe["rcsb_id"] for m in meta.values() if m for pe in m.get("polymer_entities") or []]
    clusters, cluster_info = cluster_maps(wanted)
    entries, entities = [], []
    for e, metals in entry_metals_map.items():
        f = PDB_DIR / f"{e}.cif.gz"
        status = {"downloaded": f.exists(), "file_bytes": f.stat().st_size if f.exists() else None,
                  "download_status": "ok" if f.exists() else failures.get(e, "not downloaded")}
        entries.append(complete_entry_row(e, metals, meta[e], comp_metal, status))
        ents = complete_entity_rows(meta[e], clusters)
        entities += ents
        entries[-1]["reaction"] = ";".join(sorted({x for r in ents if r["reaction"] for x in r["reaction"].split(";")}))
        entries[-1]["enzyme"] = ";".join(sorted({x for r in ents if r["enzyme"] for x in r["enzyme"].split(";")}))
    ent_df = pd.DataFrame(entities)
    for n in CLUSTER_LEVELS:
        ent_df[f"cluster{n}"] = ent_df[f"cluster{n}"].astype("Int64")
    entry_df = pd.DataFrame(entries)
    for c in ("citation_pubmed", "citation_year", "deposited_atom_count", "polymer_entity_count",
              "nonpolymer_entity_count", "file_bytes"):
        entry_df[c] = entry_df[c].astype("Int64")
    for c in ("r_free", "r_work", "resolution"):
        entry_df[c] = entry_df[c].astype("float64")
    COMPLETE.mkdir(parents=True, exist_ok=True)
    entry_df.to_parquet(COMPLETE / "entries.parquet", index=False)
    ent_df.to_parquet(COMPLETE / "entities.parquet", index=False)
    prot = ent_df[ent_df.polymer_type == "Protein"]
    fetched = sorted({json.loads(line)["fetched_utc"][:10] for line in
                      (CCACHE / "entry_meta.jsonl").read_text().splitlines() if line.strip().endswith("}")})
    update_provenance("metadata", {
        "endpoint": GRAPHQL_URL, "batch_size": META_BATCH, "query": " ".join(COMPLETE_QUERY.split()),
        "fetch_dates_utc": fetched, "entries_with_metadata": int(entry_df.metadata_found.sum()),
        "entries_without_metadata": entry_df.entry_id[~entry_df.metadata_found].tolist(),
        "polymer_entities": len(ent_df), "cluster_files": cluster_info,
        "cluster_id": "0-based line index in clusters-by-entity-<N>.txt (null: entity not in the file, e.g. "
                      "released after the clustering, or not a protein)",
        "protein_entities_without_cluster30": int(prot.cluster30.isna().sum()),
        "reaction_labels": "EC lineage of the polymer entity matched against EC_CLASSES of bio/pdb_sweep.py"})
    print(f"{len(entry_df)} entries ({int(entry_df.downloaded.sum())} downloaded), {len(ent_df)} polymer entities "
          f"-> data/complete/entries.parquet, entities.parquet", flush=True)
    print(entry_df.method.value_counts().to_string(), flush=True)


def mode_complete(steps):
    """Complete d-block sweep; steps: search, meta, download, tables (default: all, in this order)."""
    t0 = time.time()
    steps = steps or list(COMPLETE_STEPS)
    if any(s not in COMPLETE_STEPS for s in steps):
        sys.exit(f"unknown step in {steps}; choose from {COMPLETE_STEPS}")
    COMPLETE.mkdir(parents=True, exist_ok=True)
    metals = complete_search() if "search" in steps else entry_metals()
    if not metals:
        sys.exit("no search results cached: run the search step first")
    if "meta" in steps:
        complete_meta(list(metals))
    if "download" in steps:
        complete_download(list(metals))
    if "tables" in steps:
        complete_tables(metals)
        update_provenance("pdb_sweep_complete", {
            "finished_utc": utcnow(), "steps": steps, "seconds_last_run": round(time.time() - t0),
            "endpoints": {"search": SEARCH_URL, "graphql": GRAPHQL_URL, "files": FILE_URL, "ccd": CCD_URL,
                          "clusters": CLUSTER_FILE_URL, "holdings_current": HOLDINGS_URL,
                          "holdings_status": STATUS_URL},
            "software": software_versions(), "code": code_version()})


if __name__ == "__main__":
    if len(sys.argv) < 2 or sys.argv[1] not in ("enzymes", "all", "download", "complete"):
        sys.exit(__doc__)
    if sys.argv[1] == "complete":
        mode_complete(sys.argv[2:])
    else:
        {"enzymes": mode_enzymes, "all": mode_all}.get(sys.argv[1], lambda: mode_download(sys.argv[2:]))()
