"""How often has a coordination motif been reported as a single-atom catalyst? OpenAlex + Europe PMC counts.

    python lit/novelty.py                 # motifs of the DFT batches (dft/batches.py)
    python lit/novelty.py Ni-N2S2 Mo-S4   # given motifs (metal-donors, site suffixes are ignored)

For each motif the name variants (NiN2S2, Ni-N2S2, NiS2N2, Ni-S2N2, NiN2S2 with "1" subscripts, ...) are
searched together with a single-atom context clause, once alone and once with "graphene" and with "hydrogen
evolution". OpenAlex `search` covers titles, abstracts and the full texts it holds; Europe PMC covers the
life-science-indexed part of the literature with open full texts. Zero hits means "no record under these
names", not "never made": the hits of the small counts are listed for checking by hand (and with Valency).
Answers are cached in data/cache/lit/ so reruns are free.
"""
import hashlib
import itertools
import json
import os
import re
import sys
import time

import pandas as pd
import requests

HERE = os.path.dirname(os.path.abspath(__file__))
CACHE = os.path.join(HERE, "..", "data", "cache", "lit")
OUT = os.path.join(HERE, "..", "data", "novelty.csv")
CONTEXT = '("single-atom" OR "single atom" OR "atomically dispersed" OR "single-site")'
EXTRA = {"all": "", "graphene": ' AND (graphene OR "N-doped carbon" OR "nitrogen-doped carbon")',
         "HER": ' AND ("hydrogen evolution" OR HER)'}
N_TITLES = 8


def variants(motif):
    """Name variants of 'Ni-N2S2' (site suffixes such as ':DV:cis5' are dropped)."""
    metal, formula = motif.split(":")[0].split("-")
    groups = re.findall(r"([A-Z][a-z]?)(\d*)", formula)
    names = set()
    for order in itertools.permutations(groups):
        for ones in (False, True):
            d = "".join(el + (n or ("1" if ones and len(groups) > 1 else "")) for el, n in order)
            names.update({metal + d, metal + "-" + d})
    return sorted(names)


def cached_get(url, params):
    key = hashlib.sha1((url + json.dumps(params, sort_keys=True)).encode()).hexdigest()
    path = os.path.join(CACHE, key + ".json")
    if os.path.exists(path):
        return json.load(open(path))
    for attempt in range(4):
        try:
            r = requests.get(url, params=params, timeout=60)
            if r.status_code == 200:
                os.makedirs(CACHE, exist_ok=True)
                json.dump(r.json(), open(path, "w"))
                return r.json()
        except requests.RequestException:
            pass
        time.sleep(2 ** attempt)
    return None


def openalex(query):
    j = cached_get("https://api.openalex.org/works",
                   {"search": query, "per-page": N_TITLES, "select": "doi,title,publication_year"})
    if j is None:
        return None, []
    hits = ["%s %s (%s)" % (w["publication_year"], (w["title"] or "").strip(), w["doi"] or "no doi")
            for w in j["results"]]
    return j["meta"]["count"], hits


def europepmc(query):
    j = cached_get("https://www.ebi.ac.uk/europepmc/webservices/rest/search",
                   {"query": query, "format": "json", "pageSize": N_TITLES})
    return None if j is None else j.get("hitCount")


def main(motifs):
    rows = []
    for motif in motifs:
        names = " OR ".join('"%s"' % v for v in variants(motif))
        row = {"motif": motif.split(":")[0], "names": names}
        for tag, extra in EXTRA.items():
            q = "(%s) AND %s%s" % (names, CONTEXT, extra)
            row["openalex_" + tag], hits = openalex(q)
            row["europepmc_" + tag] = europepmc(q)
            if tag == "all":
                row["top_hits"] = " | ".join(hits)
        rows.append(row)
        print("%-10s OpenAlex %5s (graphene %4s, HER %4s)  Europe PMC %5s" % (
            row["motif"], row["openalex_all"], row["openalex_graphene"], row["openalex_HER"], row["europepmc_all"]))
    df = pd.DataFrame(rows).drop_duplicates("motif")
    if os.path.exists(OUT):     # keep motifs from earlier runs
        old = pd.read_csv(OUT)
        df = pd.concat([old[~old.motif.isin(df.motif)], df])
    df.to_csv(OUT, index=False)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        sys.path.insert(0, os.path.join(HERE, "..", "dft"))
        from batches import GRAPHENE_BATCH1
        main(list(dict.fromkeys(m.split(":")[0] for m, _ in GRAPHENE_BATCH1)))
