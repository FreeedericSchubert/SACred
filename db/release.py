"""Assemble the Zenodo release: one folder with the database, structure files, results, code and metadata.

    python db/build.py && python db/release.py      # -> data/zenodo/sac-sites-v<VERSION>/ (rebuilt from scratch)

Layout of the folder:
    database/   sac_motifs.sqlite.xz, parquet.zip (every table), csv.zip (every table as .csv.gz)
    structures/ extended XYZ (plain text: ASE seeks within the file, which is very slow on gzip) with energies
                (eV) and forces (eV/A): every DFT relaxation step per
                pseudopotential family, the 3x3x1 single points, the relaxed GFN2-xTB structures
    results/    dG_H tables, stability, xTB scan and calibration, literature check, the enzyme site set
    code/       snapshot of the repository at the build commit (git archive)
    README.md, DATASHEET.md, SCHEMA.md, CHANGELOG.md, LICENSE.md, CITATION.cff, zenodo_metadata.json,
    MANIFEST.sha256
The cluster account and the local user name are replaced by '~' in the code snapshot. Authors, affiliations,
ORCIDs, contact and DOI are placeholders (<...>) that the authors fill in before the upload.
"""
import glob
import json
import os
import re
import shutil
import sqlite3
import subprocess
import sys
import tarfile
import tempfile
import zipfile
from datetime import date

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
REL = os.path.join(DATA, "release")
sys.path.insert(0, HERE)
import build                                     # noqa: E402  (VERSION, DB path)
import export                                    # noqa: E402  (structure writers, settings, manifest)

TITLE = "SACred: transition-metal sites of the Protein Data Bank as templates for single-atom catalysts"
RESULTS = ["dgh_graphene.csv", "dgh_mos2.csv", "stability_graphene.csv", "novelty.csv", "literature_check.json",
           "literature_check_deep.json",
           "xtb_scan.csv", "xtb_calibration.csv", "xtb_candidates.csv", "xtb_selected.csv", "enzyme_sites.csv",
           "supports_summary.md", "xtb_supports.csv", "xtb_reactions.csv", "xtb_orr_calibration.csv",
           "xtb_orr_calibration_fit.json", "xtb_offsets.csv", "xtb_library_comparison.csv", "dft_orr.csv", "dft_orr.md", "dft_orr_sites.csv",
           "xtb_stability.csv", "xtb_stability_descriptors.csv", "xtb_stability_descriptors.md",
           "geometry_mismatch.csv", "geometry_mismatch.md", "geometry_three_way.csv", "geometry_three_way.md",
           "literature_stability_lean.json", "literature_stability_deep.json", "literature_check_final.json", "dft_conditions.csv", "dft_conditions.md",
           "dft_points.csv", "dft_points.md", "dft_points_dgh.csv"]
SCRUB = [("~", "~"), ("<user>", "<user>"), ("~", "~"), ("<user>", "<user>"),
         ("<earlier_mos2_project>", "<earlier_mos2_project>")]
KEYWORDS = ["single-atom catalysts", "metalloenzymes", "Protein Data Bank", "coordination chemistry", "graphene",
            "MoS2", "hydrogen evolution", "density functional theory", "GFN2-xTB", "machine-learned force fields"]


def size(path):
    n = (sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(path) for f in fs)
         if os.path.isdir(path) else os.path.getsize(path))
    return "%.1f MB" % (n / 1e6) if n > 1e5 else "%.0f kB" % (n / 1e3)


def zip_dir(src, dst):
    with zipfile.ZipFile(dst, "w", zipfile.ZIP_STORED) as z:        # members are already compressed
        for f in sorted(os.listdir(src)):
            z.write(os.path.join(src, f), arcname=os.path.join(os.path.basename(dst)[:-4], f))


def code_snapshot(dst):
    """git archive of HEAD with identifiers scrubbed from every text file."""
    commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    with tempfile.TemporaryDirectory() as tmp:
        tar = os.path.join(tmp, "head.tar")
        subprocess.run(["git", "-C", ROOT, "archive", "--format=tar", "-o", tar, "HEAD"], check=True)
        src = os.path.join(tmp, "src")
        with tarfile.open(tar) as t:
            t.extractall(src, filter="data")
        changed = 0
        for dirpath, _, files in os.walk(src):
            for f in files:
                p = os.path.join(dirpath, f)
                try:
                    text = open(p, encoding="utf-8").read()
                except (UnicodeDecodeError, OSError):
                    continue                                       # binary (figures)
                new = text
                for a, b in SCRUB:
                    new = new.replace(a, b)
                if new != text:
                    open(p, "w", encoding="utf-8").write(new)
                    changed += 1
        name = "sac-pipeline-%s" % commit
        out = os.path.join(dst, name + ".tar.gz")
        with tarfile.open(out, "w:gz") as t:
            t.add(src, arcname=name)
    return out, commit, changed


def db_counts():
    con = sqlite3.connect(build.DB)
    counts = {t: con.execute('select count(*) from "%s"' % t).fetchone()[0]
              for t, in con.execute("select name from sqlite_master where type='table' order by name")}
    info = dict(con.execute("select key, value from dataset_info"))
    return counts, info


def readme(target, counts, info, files, struct_counts, commit):
    dft_done = counts.get("dft_results") and sqlite3.connect(build.DB).execute(
        "select count(*) from dft_results where dG_H is not null").fetchone()[0]
    lit = sqlite3.connect(build.DB).execute(
        "select motif, verdict, round(dG_H, 2) from literature_check order by abs(dG_H)").fetchall()
    as_of = json.loads(info.get("data_as_of", "{}"))
    con = sqlite3.connect(build.DB)
    n_shells = con.execute("select count(*) from motifs where set_name = 'complete' and kind = 'first_shell'").fetchone()[0]
    n_xtb_sites = con.execute("select count(distinct model_id) from xtb_scan").fetchone()[0]
    lines = [
        "# %s" % TITLE, "",
        "Version %s, pre-release (%s). `<AUTHOR 1>`, `<AUTHOR 2>` (`<AFFILIATION>`). DOI: `<DOI>`." % (
            build.VERSION, date.today().isoformat()), "",
        "Every transition-metal site in the Protein Data Bank with its first coordination shell, mapped onto doped "
        "graphene vacancies as a template for single-atom catalysts, with a GFN2-xTB screening, DFT hydrogen "
        "adsorption free energies (ΔG_H) for the most promising sites, a literature check, and the structure files "
        "of all calculations for training machine-learned force fields.", "",
        "**Status:** %s. DFT results as of %s." % (info.get("status", ""), as_of.get("dft_results")), "",
        "## In numbers", "",
        "| | |", "| --- | --- |",
        "| PDB entries with a transition-metal component (RCSB data version 2026_41) | %s |" % format(counts["processing_log"], ","),
        "| metal sites, every copy | %s |" % format(counts["sites"], ","),
        "| donor, weak and metal contacts | %s |" % format(counts["contacts"], ","),
        "| distinct first-shell motifs (metal and donor set, complete set) | %s |" % format(n_shells, ","),
        "| graphene sites screened with GFN2-xTB / H positions | %s / %s |" % (format(n_xtb_sites, ","), format(counts["xtb_scan"], ",")),
        "| sites with a three-donor face on MoS₂ S / in an S vacancy / on silica O (every copy) | %s |" % " / ".join(
            format(con.execute("select count(*) from support_sites where support = ? and face_size = 3", (s,)).fetchone()[0], ",")
            for s in ("MoS2", "MoS2vS", "SiO2")),
        "| SAC models with DFT runs / H positions finished | %d / %d |" % (counts["sac_models"], dft_done),
        "| DFT structures in the structure files (relaxation steps + single points) | %s |" % format(sum(struct_counts.values()), ","),
        "", "## Files", "", "| File | Contents | Size |", "| --- | --- | --- |"]
    lines += ["| `%s` | %s | %s |" % (f, desc, size(os.path.join(target, f))) for f, desc in files]
    lines += [
        "", "## Quick start", "",
        "```python",
        "import lzma, shutil, sqlite3, pandas as pd",
        "with lzma.open('database/sac_motifs.sqlite.xz') as fi, open('sac_motifs.sqlite', 'wb') as fo:",
        "    shutil.copyfileobj(fi, fo)",
        "con = sqlite3.connect('sac_motifs.sqlite')",
        "pd.read_sql('select metal, count(*) n from sites where is_representative group by metal order by n desc', con)",
        "",
        "from ase.io import read            # structure files: energies in eV, forces in eV/A",
        "frames = read('structures/dft_graphene_relax_gbrv.extxyz', index=':')",
        "frames[0].info, frames[0].get_potential_energy(), frames[0].get_forces().shape",
        "```", "",
        "Every table and column is described in `SCHEMA.md` and in the table `column_docs`; `DATASHEET.md` describes "
        "how the data were collected and processed, what they can and cannot be used for, and their limitations.", "",
        "## Methods in brief", "",
        "1. Sweep: every PDB entry with a chemical component containing a d-block metal (RCSB Search and Data APIs).",
        "2. First shells with gemmi: element-specific donor cutoffs, crystal symmetry, alternate conformers, weak "
        "contacts kept separately; checked against the PDB's own metal-link annotation (`validation_metalc`).",
        "3. Graphene mapping: four nearly square-planar donors into a double vacancy, otherwise the best three-donor "
        "face into a single vacancy; His N becomes pyridinic N, Cys S an S dopant, carboxylate or water O an O dopant.",
        "4. GFN2-xTB scan of the common donor environments with nine enzyme metals, calibrated against DFT.",
        "5. DFT (Quantum ESPRESSO 7.1, PBE-D3(BJ)): ΔG_H = ΔE_H + 0.24 eV; settings per dataset in `dft_methods`.",
        "6. Literature check of the sites near ΔG_H = 0 (abstract level).", "",
        "## Structure files for force-field training", "",
        "- The GBRV files (graphene, bulk metals) and the pslibrary file (MoS₂) have no common energy zero: fit them "
        "separately or with separate reference energies.",
        "- Energies and forces include the D3(BJ) correction; with smearing the energy is the free energy the forces "
        "belong to. Relaxation steps are correlated and near minima; the xTB structures are minima at the xTB level.",
        "- Each structure carries `run`, `motif`, `state` (`bare`, `H_M` H on the metal, `H_S` H on an S donor, `ref`, "
        "`bulk`), `step` and `final` in its info line.", "",
        "## Literature check of the sites near the optimum", "",
        "| Site | ΔG_H (eV) | Verdict |", "| --- | --- | --- |"]
    lines += ["| %s | %+.2f | %s |" % (m, g, v.replace("_", " ")) for m, v, g in lit]
    lines += [
        "", "Search by meaning and by name, abstracts of the closest papers, then their open-access full texts, citations "
        "and screening series (e.g. M-NxO4-x); paywalled full texts and supporting information were not read. \"Related "
        "only\" means the same elements with another donor count or support; \"reported computational only\" includes "
        "large DFT screenings that computed the site for other reactions. Details and DOIs: `results/literature_check.json`, "
        "`results/literature_check_deep.json` and the table `literature_papers`.", "",
        "## Licence and citation", "",
        "Database, tables and results: CC BY 4.0. Structure files (`structures/`): CC0 1.0. Source structures: wwPDB "
        "archive, CC0 1.0. Code snapshot: MIT. See `LICENSE.md`; citation metadata in `CITATION.cff`. Calculations on the "
        "Curta HPC cluster of Freie Universität Berlin.", "",
        "All values are drafts from a two-day hackathon and need expert checking. Code snapshot: commit %s." % commit]
    open(os.path.join(target, "README.md"), "w").write("\n".join(lines) + "\n")


def metadata(target, info):
    creators = [{"name": "<LAST NAME>, <FIRST NAME>", "affiliation": "<AFFILIATION>", "orcid": "<ORCID>"},
                {"name": "<LAST NAME>, <FIRST NAME>", "affiliation": "<AFFILIATION>", "orcid": "<ORCID>"}]
    desc = ("<p>Every transition-metal site in the Protein Data Bank with its first coordination shell, mapped onto "
            "doped graphene vacancies as templates for single-atom catalysts; GFN2-xTB screening, DFT hydrogen "
            "adsorption free energies for the most promising sites, a literature check, and the structure files of all "
            "calculations (energies and forces) for training machine-learned force fields.</p><p>%s</p>"
            "<p>Licences: database, tables and results CC BY 4.0; structure files (structures/) CC0 1.0; code snapshot (code/) "
            "MIT. See LICENSE.md, README.md and DATASHEET.md.</p>" % info.get("status", ""))
    json.dump({"metadata": {
        "title": "%s (v%s)" % (TITLE, build.VERSION), "upload_type": "dataset", "description": desc,
        "creators": creators, "keywords": KEYWORDS, "license": "cc-by-4.0", "version": build.VERSION,
        "access_right": "open", "language": "eng", "publication_date": "<YYYY-MM-DD>",
        "notes": "Pre-release from a two-day hackathon; values are drafts. Fill in the <...> fields before upload.",
        "related_identifiers": [{"identifier": "<code repository URL, once public>", "relation": "isSupplementedBy",
                                 "scheme": "url"}]}}, open(os.path.join(target, "zenodo_metadata.json"), "w"), indent=1)
    open(os.path.join(target, "CITATION.cff"), "w").write("\n".join([
        "cff-version: 1.2.0", 'message: "If you use this dataset, please cite it as below."',
        'title: "%s"' % TITLE, 'version: "%s"' % build.VERSION, "type: dataset",
        "authors:", '  - family-names: "<LAST NAME>"', '    given-names: "<FIRST NAME>"',
        '    affiliation: "<AFFILIATION>"', '    orcid: "https://orcid.org/<ORCID>"',
        '  - family-names: "<LAST NAME>"', '    given-names: "<FIRST NAME>"',
        '    affiliation: "<AFFILIATION>"', '    orcid: "https://orcid.org/<ORCID>"',
        "license:", "  - CC-BY-4.0", "  - CC0-1.0", "  - MIT", 'doi: "<DOI>"', 'date-released: "<YYYY-MM-DD>"', "keywords:"]
        + ["  - %s" % k for k in KEYWORDS]) + "\n")
    open(os.path.join(target, "LICENSE.md"), "w").write(
        "# Licences\n\n"
        "- Database, tables and results (`database/`, `results/`): Creative Commons Attribution 4.0 International "
        "(CC BY 4.0), https://creativecommons.org/licenses/by/4.0/\n"
        "- Structure files (`structures/`: DFT relaxation steps, single points, xTB structures): CC0 1.0 Universal "
        "(public domain dedication), https://creativecommons.org/publicdomain/zero/1.0/\n"
        "- Code snapshot in `code/`: MIT License (`LICENSE` in the snapshot).\n"
        "- Source structures from the wwPDB archive: CC0 1.0 (https://www.wwpdb.org/about/usage-policies).\n")
    open(os.path.join(target, "CHANGELOG.md"), "w").write(
        "# Changelog\n\n"
        "## %s (%s, hackathon release)\n\n"
        "- DFT: H adsorption on graphene sites from the xTB scan and library, metal formation energies and "
        "dissolution potentials against the metal-free vacancies (`dft_results`, `dft_stability`), O₂ reduction "
        "intermediates on Co–N₄ and Co–N₂OC checked against Cho et al. 2022 (`results/dft_orr.csv`); %s.\n"
        "- GFN2-xTB library: every N/S/O/C arrangement in graphene (20 single-vacancy, 76 double-vacancy) for the nine "
        "enzyme metals plus Cr, Ru and Pt; per-metal calibration with leave-one-out errors (`xtb_offsets`).\n"
        "- O₂, CO₂ and N₂ intermediates with GFN2-xTB (`xtb_reactions`), with the O₂ check against published DFT.\n"
        "- Stability and synthesis descriptors: xTB insertion and escape energies for every arrangement, NEB "
        "paths for the DFT sites, descriptor analysis (`results/xtb_stability*`, `results/xtb_neb.csv`).\n"
        "- Geometry of the first shell: enzyme against support model and free complex (`results/geometry_*`).\n"
        "- Literature: novelty of the sites near the optimum, and a stability/synthesis sweep with in-depth reads "
        "(`results/literature_*`).\n\n"
        "## 0.3 (22 Oct 2026, pre-release)\n\n"
        "- Complete PDB sweep (every entry with a transition-metal component) with validation against the PDB's "
        "metal-link annotation.\n"
        "- 58 SAC models with DFT runs (28 picked from enzymes, 30 from the xTB scan).\n"
        "- GFN2-xTB scan of 423 graphene sites with per-metal calibration (`xtb_scan`, `xtb_calibration`).\n"
        "- Stability references (bulk metals done, metal-free vacancies pending) and MoS₂ results.\n"
        "- Literature check of the sites near the optimum (`literature_check`, `literature_papers`).\n"
        "- Every site mapped onto MoS₂ (lattice S, S vacancy) and silica (siloxide O) (`support_sites`, "
        "`support_motifs`), with a first GFN2-xTB look at hydrogen binding on those supports (`xtb_supports`).\n"
        "- Structure files (energies, forces) of all DFT relaxation steps and xTB structures; code snapshot.\n\n"
        "## 0.2 (21 Oct 2026, internal)\n\n- Complete set, technical validation, datasheet.\n\n"
        "## 0.1 (20 Oct 2026, internal)\n\n- Enzyme and pdb30 site sets, first DFT batch.\n" % (
            build.VERSION, date.today().isoformat(), info.get("status", "")))


def neb_set(out):
    """Every image of the GFN2-xTB NEB paths (data/screen/xtb_neb/<motif>.xyz) with motif and image number."""
    from ase.io import read, write
    images = []
    for xyz in sorted(glob.glob(os.path.join(DATA, "screen", "xtb_neb", "*.xyz"))):
        for k, at in enumerate(read(xyz, index=":", format="extxyz")):
            at.info.update(dataset="xtb-neb", motif=os.path.basename(xyz)[:-4], image=k)
            images.append(at)
    write(out, images, format="extxyz")
    return len(images)


def neb_table(out):
    """One row per NEB site: barrier and reaction energy (eV), or why there is none."""
    import pandas as pd
    rows = [dict(json.load(open(f)), motif=os.path.basename(f)[:-5])
            for f in sorted(glob.glob(os.path.join(DATA, "screen", "xtb_neb", "*.json")))]
    if rows:
        pd.DataFrame(rows).to_csv(out, index=False, float_format="%.4f")


def main():
    target = os.path.join(DATA, "zenodo", "sac-sites-v%s" % build.VERSION)
    assert target.startswith(os.path.join(DATA, "zenodo") + os.sep)
    if os.path.exists(target):
        shutil.rmtree(target)                                    # generated by this script only
    for d in ("database", "structures", "results", "code"):
        os.makedirs(os.path.join(target, d))
    files = []

    with open(build.DB, "rb") as fi, open(os.path.join(target, "database", "sac_motifs.sqlite.xz"), "wb") as fo:
        subprocess.run(["xz", "-T0", "-6", "-c"], stdin=fi, stdout=fo, check=True)
    zip_dir(os.path.join(REL, "parquet"), os.path.join(target, "database", "parquet.zip"))
    zip_dir(os.path.join(REL, "csv"), os.path.join(target, "database", "csv.zip"))
    files += [("database/sac_motifs.sqlite.xz", "SQLite database, every table, with views (xz)"),
              ("database/parquet.zip", "every table as Parquet"), ("database/csv.zip", "every table as gzipped CSV")]

    struct_counts = {}
    with tempfile.TemporaryDirectory() as tmp:
        runs = export.RUNS
        jobs = [("dft_graphene_relax_gbrv", lambda o: export.dft_set(os.path.join(runs, "graphene", "*", "*"), "graphene-gbrv", o)[0],
                 "every ionic step of the graphene DFT relaxations (GBRV)"),
                ("dft_bulk_vcrelax_gbrv", lambda o: export.dft_set(os.path.join(runs, "bulk", "*"), "bulk-gbrv", o)[0],
                 "bulk-metal cell relaxations (GBRV, stability references)"),
                ("dft_mos2_relax_pslibrary", lambda o: export.dft_set(os.path.join(runs, "mos2", export.MOS2_GLOB), "mos2-pslibrary", o)[0],
                 "MoS₂ relaxations (pslibrary PAW)"),
                ("dft_graphene_scf_3x3x1_gbrv", export.scf_set, "3×3×1 k-point single points behind ΔG_H"),
                ("xtb_graphene_relaxed", export.xtb_set, "relaxed GFN2-xTB structures of the scan"),
                ("xtb_supports_relaxed", export.xtb_supports_set,
                 "relaxed GFN2-xTB structures of metal sites on MoS₂ (3×3) and on a silsesquioxane silica model"),
                ("xtb_reactions_relaxed", lambda o: export.xtb_dir_set(
                    "xtb_reactions", "xtb-reactions", ("kind", "site", "state", "converged", "steps", "M_X", "bound"), o),
                 "relaxed GFN2-xTB O₂, CO₂ and N₂ intermediates and gas molecules"),
                ("xtb_empty_sites_relaxed", lambda o: export.xtb_dir_set(
                    "xtb_sites", "xtb-sites", ("site", "order", "formula", "converged"), o),
                 "relaxed GFN2-xTB vacancies without the metal (stability references)"),
                ("xtb_adatoms_relaxed", lambda o: export.xtb_dir_set(
                    "xtb_adatoms", "xtb-adatoms", ("metal", "start", "converged", "height"), o),
                 "relaxed GFN2-xTB metal adatoms on clean graphene (escape energies)"),
                ("xtb_neb_paths", neb_set, "GFN2-xTB climbing-image NEB paths, site to adatom (every image)")]
        for name, write, desc in jobs:
            struct_counts[name] = write(os.path.join(target, "structures", name + ".extxyz"))
            files.append(("structures/%s.extxyz" % name, "%s: %s structures" % (desc, format(struct_counts[name], ","))))
    json.dump(export.SETTINGS, open(os.path.join(target, "structures", "settings.json"), "w"), indent=1)
    files.append(("structures/settings.json", "calculation settings per dataset"))

    for f in RESULTS:
        if os.path.exists(os.path.join(DATA, f)):
            shutil.copy2(os.path.join(DATA, f), os.path.join(target, "results", f))
    neb_table(os.path.join(target, "results", "xtb_neb.csv"))
    files.append(("results/", "dG_H, stability, xTB scan and calibration, literature check, enzyme site set"))
    code, commit, scrubbed = code_snapshot(os.path.join(target, "code"))
    files.append(("code/" + os.path.basename(code), "repository snapshot (identifiers in paths replaced by ~)"))

    for f in ("DATASHEET.md", "SCHEMA.md"):
        shutil.copy2(os.path.join(HERE, f), os.path.join(target, f))
    counts, info = db_counts()
    metadata(target, info)
    readme(target, counts, info, files, struct_counts, commit)
    n = export.manifest(target)
    total = sum(os.path.getsize(os.path.join(dp, f)) for dp, _, fs in os.walk(target) for f in fs)
    print("release %s: %d files, %.0f MB, structures %s, code %s (%d files scrubbed) -> %s" % (
        build.VERSION, n, total / 1e6, struct_counts, commit, scrubbed, target))


if __name__ == "__main__":
    main()
