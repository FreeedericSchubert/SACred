"""Copy the database, the result tables and the structure data for force-field training to one folder.

    python db/build.py && python db/export.py ~/Nextcloud/SAC_database

Layout of the target folder (see the README.md written there):
    database/   the release database (SQLite, Parquet, CSV) with DATASHEET.md and SCHEMA.md
    results/    dG_H, stability and literature tables
    xtb/        GFN2-xTB scan tables and every relaxed xTB structure (extxyz: energy, forces, charges)
    dft/        every ionic step of the Quantum ESPRESSO relaxations (extxyz: energy, forces, magnetic
                moments), one file per pseudopotential family, plus the 3x3x1 single points
Energies are in eV, forces in eV/A. Files of different pseudopotential families (GBRV vs pslibrary) have no
common energy zero and must not be pooled in one fit without separate reference energies.
"""
import datetime
import glob
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys

from ase.io import read, write

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(ROOT, "data")
RUNS = os.path.join(ROOT, "dft", "runs")
RESULTS = ["dgh_graphene.csv", "dgh_mos2.csv", "stability_graphene.csv", "novelty.csv", "enzyme_sites.csv"]
XTB_TABLES = ["xtb_scan.csv", "xtb_calibration.csv", "xtb_offsets.csv", "xtb_candidates.csv", "xtb_selected.csv",
              "xtb_library_comparison.csv", "xtb_supports.csv", "xtb_reactions.csv", "xtb_orr_calibration.csv",
              "xtb_stability.csv", "geometry_mismatch.csv"]
MOS2_GLOB = "3[0134]*"   # MoS2 runs of the dataset (runs 32* are not part of it)
MOS2_LABELS = {   # MoS2 runs have no meta.json: (motif, state)
    "30_h2_alone": ("H2", "ref"), "31_mos2_h": ("MoS2", "H_S"),
    "33_mos2_fe": ("Fe-S3 on MoS2, above Mo", "bare"),
    "33h_mos2_fe_hollow": ("Fe-S3 on MoS2, over the hollow", "bare"), "34_mos2_fe_h": ("Fe-S3 on MoS2, above Mo", "H_M")}
SETTINGS = {
    "graphene-gbrv": "QE 7.1, PBE-D3(BJ), GBRV USPP 40/200 Ry, spin-polarised, MV smearing 0.01 Ry, dipole correction, "
                     "6x6 graphene cell (c = 20 A), Gamma point, one far C fixed",
    "graphene-gbrv-scf331": "as graphene-gbrv, single point with 3x3x1 k points at the relaxed geometry",
    "bulk-gbrv": "QE 7.1, PBE-D3(BJ), GBRV USPP 40/200 Ry, spin-polarised, MV smearing 0.01 Ry, vc-relax, "
                 "16x16x16 (bcc, fcc) or 16x16x10 (hcp) k points",
    "mos2-pslibrary": "QE 7.1, PBE-D3(BJ), pslibrary 1.0.0 PAW 50/400 Ry (Fe: GBRV USPP), Gaussian smearing "
                      "0.0037 Ry, dipole correction, 6x6 MoS2 cell (c = 35 A), Gamma point, Mo fixed",
    "xtb": "GFN2-xTB (tblite 0.7), periodic, Gamma point, spin-restricted, Fermi smearing 1000 K, "
           "BFGS to 0.05 eV/A, one far C fixed; relaxed structures only",
    "xtb-supports": "GFN2-xTB (tblite 0.7), spin-restricted, Fermi smearing 1000 K, BFGS to 0.05 eV/A; MoS2: 3x3 cut "
                    "of the DFT sheet, periodic, Gamma point, Mo fixed; SiO2: Me7Si7O9(O)2(OH)M molecule; relaxed "
                    "structures only",
    "xtb-reactions": "settings of xtb; O2, CO2 and N2 intermediates (OH, O, OOH, COOH, CO, OCHO, N2, N2H, NH2, NH3) on "
                     "the metal, gas molecules in the support's own cell; relaxed structures only",
    "xtb-sites": "settings of xtb; the 96 donor arrangements of graphene without a metal (screen/xtb_stability.py)",
    "xtb-adatoms": "settings of xtb; clean 6x6 graphene and each metal on it, started on a hollow and on top of a C",
}


def run_outputs(run):
    """relax.out pieces of one run in time order: relax.out.1, relax.out.2, ..., relax.out."""
    parts = sorted(glob.glob(os.path.join(run, "relax.out.*")), key=lambda p: int(p.rsplit(".", 1)[1]))
    return parts + [os.path.join(run, "relax.out")] * os.path.exists(os.path.join(run, "relax.out"))


def frames(path):
    """All ionic steps of a pw.x output; a step cut off by the time limit is dropped (and dft_set drops images
    without energy and forces, e.g. the bare final coordinates of a vc-relax)."""
    try:
        return read(path, index=":", format="espresso-out")
    except Exception:
        n = len(re.findall(r"^!", open(path).read(), re.M))
        for k in range(n - 1, 0, -1):
            try:
                return read(path, index=":%d" % k, format="espresso-out")
            except Exception:
                continue
    return []


def dft_set(pattern, dataset, out):
    images, runs, skipped = [], 0, []
    for run in sorted(glob.glob(pattern)):
        meta = json.load(open(os.path.join(run, "meta.json"))) if os.path.exists(os.path.join(run, "meta.json")) else {}
        if not meta and os.path.basename(run) in MOS2_LABELS:
            meta = dict(zip(("motif", "state"), MOS2_LABELS[os.path.basename(run)]))
        pieces = run_outputs(run)
        if not pieces:
            continue
        finished = "End final coordinates" in open(pieces[-1]).read()
        steps = []
        for p in pieces:
            got = [a for a in frames(p) if a.calc is not None and {"energy", "forces"} <= set(a.calc.results)]
            if not got:
                skipped.append(os.path.relpath(p, RUNS))
            steps += got
        for k, at in enumerate(steps):
            at.info.update(run=os.path.relpath(run, RUNS), motif=meta.get("motif", ""), state=meta.get("state", ""),
                           step=k, final=bool(finished and k == len(steps) - 1), dataset=dataset)
        images += steps
        runs += 1
    write(out, images, format="extxyz")
    return len(images), runs, skipped


def scf_set(out):
    images = []
    for path in sorted(glob.glob(os.path.join(RUNS, "graphene", "*", "*", "scf.out"))):
        text = open(path).read()
        if "convergence has been achieved" not in text or "JOB DONE" not in text:
            continue
        run = os.path.dirname(path)
        meta = json.load(open(os.path.join(run, "meta.json"))) if os.path.exists(os.path.join(run, "meta.json")) else {}
        at = read(path, format="espresso-out", index=-1)
        at.info.update(run=os.path.relpath(run, RUNS), motif=meta.get("motif", ""), state=meta.get("state", ""),
                       dataset="graphene-gbrv-scf331")
        images.append(at)
    write(out, images, format="extxyz")
    return len(images)


def xtb_set(out):
    images = []
    for js in sorted(glob.glob(os.path.join(DATA, "screen", "xtb", "*", "*.json"))):
        res = json.load(open(js))
        xyz = js[:-5] + ".xyz"
        if not os.path.exists(xyz) or res.get("energy") is None:
            continue
        at = read(xyz, format="extxyz")
        at.info.update({k: v for k, v in res.items() if k in ("motif", "state", "converged", "steps", "max_M_X", "H_partner")
                        and v is not None}, dataset="xtb")
        images.append(at)
    write(out, images, format="extxyz")
    return len(images)


def xtb_supports_set(out):
    """Relaxed GFN2-xTB structures of screen/xtb_supports.py (metal sites on MoS2 and on the silica model)."""
    images = []
    for js in sorted(glob.glob(os.path.join(DATA, "screen", "xtb_supports", "*", "*.json"))):
        res = json.load(open(js))
        xyz = js[:-5] + ".xyz"
        if not os.path.exists(xyz) or res.get("energy") is None:
            continue
        at = read(xyz, format="extxyz")
        at.info.update({k: v for k, v in res.items() if k in ("support", "metal", "site", "state", "converged", "steps",
                                                              "H_partner") and v not in (None, "")},
                       dataset="xtb-supports")
        images.append(at)
    write(out, images, format="extxyz")
    return len(images)


def xtb_dir_set(subdir, dataset, keys, out):
    """Every relaxed structure under data/screen/<subdir> (one .json + .xyz each) with the given result keys."""
    images = []
    for js in sorted(glob.glob(os.path.join(DATA, "screen", subdir, "**", "*.json"), recursive=True)):
        res = json.load(open(js))
        xyz = js[:-5] + ".xyz"
        if not os.path.exists(xyz) or res.get("energy") is None:
            continue
        at = read(xyz, format="extxyz")
        at.info.update({k: v for k, v in res.items() if k in keys and v not in (None, "")}, dataset=dataset)
        images.append(at)
    write(out, images, format="extxyz")
    return len(images)


def copy(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    (shutil.copytree(src, dst, dirs_exist_ok=True) if os.path.isdir(src) else shutil.copy2(src, dst))


def readme(target, counts, commit):
    lines = [
        "# SAC database: enzyme metal sites as single-atom catalysts",
        "",
        "Draft data from a two-day hackathon (exported %s, code commit %s). All values still need expert checking."
        % (datetime.datetime.now().strftime("%Y-%m-%d %H:%M"), commit),
        "",
        "## Contents",
        "",
        "- `database/`: every transition-metal site of the PDB (RCSB data version 2026_41) with its first shell, weak "
        "contacts, graphene mapping, entries, entities, motif statistics, SAC models, DFT results and literature counts. "
        "SQLite (one file, with views), the same tables as Parquet and gzipped CSV; `DATASHEET.md` describes the data, "
        "`SCHEMA.md` every column.",
        "- `results/`: hydrogen adsorption free energies (graphene, MoS2), metal stability (formation energy, "
        "dissolution potential; filled as the runs finish), literature counts, the enzyme site set.",
        "- `xtb/`: the GFN2-xTB scan and library (%d relaxed structures in `xtb_structures.extxyz`), its calibration "
        "against DFT, and the tables of the support, reaction, stability and geometry screens; further relaxed "
        "structures: %s." % (counts["xtb"], ", ".join("`%s` %d" % kv for kv in counts.get("xtb_more", {}).items())),
        "- `dft/`: every ionic step of the DFT relaxations as extended XYZ with energy (eV), forces (eV/A) and, where "
        "printed, magnetic moments:",
    ]
    for name, (n, runs) in counts["dft"].items():
        lines.append("  - `%s`: %d structures from %d runs" % (name, n, runs))
    lines += [
        "  - `graphene_scf_3x3x1.extxyz`: %d single points with 3x3x1 k points (the energies behind dG_H)" % counts["scf"],
        "",
        "Each structure carries `run`, `motif` (e.g. `Fe-NS2:SV` = Fe with one N and two S donors in a single vacancy), "
        "`state` (`bare`, `H_M` = H on the metal, `H_S` = H on an S donor, `ref`, `bulk`), `step` and `final` "
        "(last step of a finished relaxation) in its info line.",
        "",
        "## Settings",
        "",
    ]
    lines += ["- `%s`: %s" % (k, v) for k, v in SETTINGS.items()]
    lines += [
        "",
        "## Notes for force-field training",
        "",
        "- The GBRV files (graphene, bulk) and the pslibrary file (MoS2) have no common energy zero: fit them separately "
        "or with separate reference energies per dataset. Their forces are both PBE-D3 but differ slightly with the "
        "pseudopotentials.",
        "- Energies and forces include the D3(BJ) dispersion correction; with smearing the energy is the free energy "
        "that the forces belong to.",
        "- Relaxation steps are strongly correlated and lie near minima; MD needs more distorted configurations "
        "(e.g. rattled structures or short MD snapshots labelled with DFT).",
        "- The xTB structures are relaxed minima at the xTB level (forces near zero), useful for pre-training or "
        "coverage of composition space, not as DFT labels.",
        "- Only the metal site and its support are in these cells: there is no protein second shell yet.",
        "",
        "## Sources and licence",
        "",
        "Structures from the RCSB PDB (CC0). Derived data: licence to be decided by the authors. See `MANIFEST.sha256` "
        "for file checksums.",
    ]
    open(os.path.join(target, "README.md"), "w").write("\n".join(lines) + "\n")


def manifest(target):
    out = []
    for dirpath, _, files in sorted(os.walk(target)):
        for f in sorted(files):
            p = os.path.join(dirpath, f)
            if f == "MANIFEST.sha256" or os.path.basename(dirpath).startswith("."):
                continue
            h = hashlib.sha256()
            with open(p, "rb") as fh:
                for block in iter(lambda: fh.read(1 << 20), b""):
                    h.update(block)
            out.append("%s  %s" % (h.hexdigest(), os.path.relpath(p, target)))
    open(os.path.join(target, "MANIFEST.sha256"), "w").write("\n".join(out) + "\n")
    return len(out)


def main(target):
    os.makedirs(target, exist_ok=True)
    rel = os.path.join(DATA, "release")
    for name in os.listdir(rel):
        copy(os.path.join(rel, name), os.path.join(target, "database", name))
    for name in ("DATASHEET.md", "SCHEMA.md"):
        copy(os.path.join(ROOT, "db", name), os.path.join(target, "database", name))
    for name in RESULTS:
        copy(os.path.join(DATA, name), os.path.join(target, "results", name))
    for name in XTB_TABLES:
        copy(os.path.join(DATA, name), os.path.join(target, "xtb", name))
    counts = {"xtb": xtb_set(os.path.join(target, "xtb", "xtb_structures.extxyz")), "dft": {},
              "xtb_more": {
                  "xtb_supports_structures.extxyz": xtb_supports_set(os.path.join(target, "xtb", "xtb_supports_structures.extxyz")),
                  "xtb_reactions_structures.extxyz": xtb_dir_set(
                      "xtb_reactions", "xtb-reactions", ("kind", "site", "state", "converged", "steps", "M_X", "bound"),
                      os.path.join(target, "xtb", "xtb_reactions_structures.extxyz")),
                  "xtb_empty_sites.extxyz": xtb_dir_set("xtb_sites", "xtb-sites", ("site", "order", "formula", "converged"),
                                                        os.path.join(target, "xtb", "xtb_empty_sites.extxyz")),
                  "xtb_adatoms.extxyz": xtb_dir_set("xtb_adatoms", "xtb-adatoms", ("metal", "start", "converged", "height"),
                                                    os.path.join(target, "xtb", "xtb_adatoms.extxyz"))}}
    os.makedirs(os.path.join(target, "dft"), exist_ok=True)
    for name, pattern, dataset in (
            ("graphene_relax_gbrv.extxyz", os.path.join(RUNS, "graphene", "*", "*"), "graphene-gbrv"),
            ("bulk_vcrelax_gbrv.extxyz", os.path.join(RUNS, "bulk", "*"), "bulk-gbrv"),
            ("mos2_relax_pslibrary.extxyz", os.path.join(RUNS, "mos2", MOS2_GLOB), "mos2-pslibrary")):
        n, runs, skipped = dft_set(pattern, dataset, os.path.join(target, "dft", name))
        counts["dft"][name] = (n, runs)
        print("%-30s %6d structures from %3d runs%s" % (name, n, runs, "; unreadable: %s" % skipped if skipped else ""))
    counts["scf"] = scf_set(os.path.join(target, "dft", "graphene_scf_3x3x1.extxyz"))
    json.dump(SETTINGS, open(os.path.join(target, "dft", "settings.json"), "w"), indent=1)
    commit = subprocess.run(["git", "-C", ROOT, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    readme(target, counts, commit)
    print("xtb structures %d, scf single points %d; %d files in MANIFEST.sha256 -> %s" % (
        counts["xtb"], counts["scf"], manifest(target), target))


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "data", "export"))
