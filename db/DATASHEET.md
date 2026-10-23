# Datasheet: transition-metal first shells of the PDB mapped to single-atom catalyst sites

Version 0.3, pre-release (21 October 2026). Structure after Gebru et al., "Datasheets for Datasets", Commun. ACM
64, 86 (2021). Fields in angle brackets (`<AUTHOR 1>`) are filled in by the authors before the release. Column-level documentation: `SCHEMA.md` (generated) and the table
`column_docs` inside the database.

## Motivation

- **Purpose.** Metalloenzymes reduce protons, O₂, CO₂ and N₂ at metal sites whose first coordination shells
  (e.g. Fe–N₄, Ni–S₂N₂, Mo–S₄) can serve as templates for single-atom catalysts (SACs) on graphene, MoS₂ or
  silica. Bio-inspired SAC design has so far started from a few hand-picked enzymes. This dataset lists every
  transition-metal site of the selected PDB entries with its first shell, a geometric mapping to graphene
  vacancy sites, and, for a subset of motifs, DFT hydrogen adsorption energies and literature counts.
- **Created by.** `<AUTHOR 1>`, `<AUTHOR 2>` (`<AFFILIATION>`), a two-person hackathon project. Built with the code in this repository (`bio/`, `dft/`,
  `lit/`, `db/`).
- **Funding.** `<FUNDING, or "none">`. Calculations on the Curta HPC cluster of Freie Universität Berlin.

## Composition

- **Instances.** One row per metal site (`sites`), with one row per donor atom (`donors`), weak contacts and
  metal–metal contacts, plus the PDB entries and polymer entities they come from (`entries`, `entities`).
  Aggregates: motif counts (`motifs`). Derived results: SAC models built from motifs (`sac_models`), DFT results
  (`dft_results`, settings in `dft_methods`), metal stability of the sites near ΔG_H = 0 (`dft_stability`:
  formation energy against the bulk metal and dissolution potential), literature counts (`literature`), the
  GFN2-xTB screening of 423 graphene sites with its DFT calibration (`xtb_scan`, `xtb_calibration`) and a
  literature check of the sites near ΔG_H = 0 (`literature_check`, `literature_papers`). The same first shells
  are also mapped onto MoS₂ (three lattice S, or a metal in an S vacancy) and onto silica (siloxide O)
  (`support_sites`, `support_motifs`, `bio/supports.py`), with a first GFN2-xTB look at hydrogen binding on the
  bare sites (`xtb_supports`).
- **Structure files** (release only, extended XYZ with energies in eV and forces in eV/Å): every ionic step of the
  DFT relaxations, one file per pseudopotential family (graphene and bulk metals: GBRV; MoS₂: pslibrary), the
  3×3×1 k-point single points behind ΔG_H, and every relaxed xTB structure. Intended uses include fine-tuning
  machine-learned force fields; the GBRV and pslibrary files have no common energy zero.
- **Site sets.** `site_sets` marks every site with the sets it belongs to.
  - `complete`: every PDB entry with a chemical component that contains a d-block metal (RCSB data version
    2026_41, released up to 2026-10-14): 53,940 entries, all downloaded and processed without failure;
    291,977 metal sites (first model, every copy of a site), 198,687 of them representative (one per copy group,
    lowest metal B-factor). 1,526,693 contacts: 1,246,430 donors, 114,622 weak contacts, 13,082 atoms of
    elements without a donor cutoff, 152,559 metal-metal contacts.
  - `enzyme` (1,811 sites): enzymes of H₂ evolution (EC 1.12), O₂ reduction (laccase, ascorbate and bilirubin
    oxidase, ferroxidase, cytochrome c, bd and bo₃ oxidases), CO₂ reduction or activation (CO dehydrogenases,
    formate dehydrogenases, acetyl-CoA synthase, formylmethanofuran dehydrogenase, methyl-coenzyme M reductase,
    carbonic anhydrase) and N₂ fixation (nitrogenase); ≤ 2.5 Å; one representative per 95 % sequence-identity
    cluster.
  - `pdb30` (7,142 sites): entries with a Fe, Ni, Co, Cu, Mo, W, Mn or V component; ≤ 2.0 Å; one representative
    per 30 % cluster.
  The `enzyme` and `pdb30` sites are all contained in `complete`, with identical mappings, and are always
  the representative copy there.
- **Metal class.** `sites.metal_class` separates the metals of enzymes (V, Mn, Fe, Co, Ni, Cu, Zn, Mo, W:
  233,648 sites in 51,157 entries) from the others (58,329 sites in 3,449 entries), which are almost always
  soaks, drugs, phasing agents or additives; osmium alone has 41,785 sites, mostly osmium hexammine in
  ribosome crystals. The flag is per element, not a check of the single site.
- **Labels.** Reaction and enzyme labels come from EC numbers of the host polymer entity (or the entry);
  graphene mappings are computed, not observed; DFT values are model results.
- **Missing information.** Hydrogen positions are not used. Only the first model of multi-model entries is
  read. Entries without a transition metal (apo forms, metal-free formate dehydrogenases) appear in `entries`
  but contribute no sites.
- **Redundancy.** Every copy of a site is a row; `copy_group` and `is_representative` select one per group.
  Sequence clusters at 30/50/95/100 % identity (RCSB) allow redundancy filtering across entries.
- **Personal data.** None. The dataset holds no author names; entries are referenced by PDB id.

## Collection process

- RCSB PDB Search API v2 and Data API (GraphQL), queried on 2026-10-20; structures (mmCIF) from
  files.rcsb.org; metal-containing chemical components from the wwPDB Chemical Component Dictionary;
  sequence clusters from the RCSB cluster files. Queries and thresholds: `data/sweep_summary.md`,
  `bio/pdb_sweep.py`.
- Literature counts: OpenAlex and Europe PMC REST APIs, queried 2026-10-20/21 with name variants of each motif
  and a single-atom context term (`lit/novelty.py`). Counts are indicators, not proof that a motif was or was
  not reported.

## Preprocessing

- First shells with gemmi: hydrogens removed; highest-occupancy alternative conformer; neighbour search with
  crystal symmetry; donor cutoffs N/O 2.6, S 2.8, Se 2.9, C 2.3, Cl 2.8 Å; contacts up to cutoff + 0.5 Å are kept
  as weak contacts; minority states (occupancy < 0.5 × metal occupancy) and alternative ligand copies are
  separated from donors. Rules added after inspecting structures are listed in `data/sweep_summary.md`.
- Graphene mapping: DV = most square-planar subset of four donors (score ≤ 55°, the minimum between the two
  modes of the score distribution); SV = three-donor face with the most labile donor left out.
- Validation (`bio/validate.py`, `data/validation/validation.md`):
  - Against the PDB's own metal-link annotation (`_struct_conn`, metalc) of the same files: 78.4 % of the sites
    are annotated; for those, 82.9 % have exactly the annotated inter-residue donors, precision 0.979 (our
    donors that are annotated) and recall 0.944 (annotated partners we count). Matched distances agree to
    0.005 Å for 99 % of the links (median difference 0.000 Å), which checks the symmetry handling. Most missed
    partners are annotated links longer than our cutoffs; they are kept as weak contacts.
  - Threshold sensitivity: moving every donor cutoff by ±0.1 Å changes the donor formula of 6.8–8.2 % of the
    sites; 8.3 % of the four-coordinate sites lie within 10° of the double-vacancy threshold.
  - 19 of 20 literature cases reproduce the expected first shell (`data/sanity_checks.csv`); the failing case
    (1WUI) is an oxidised, partly occupied state.
- Raw data: the PDB files are not redistributed; every site carries its PDB id, chain, residue and atom name.

## Uses

- Intended: choosing SAC motifs by biological precedent, statistics of metal coordination in enzymes,
  training data for coordination-geometry or catalyst models, starting structures for DFT.
- Not intended without further checks: claims about the physiological relevance of a single site
  (adventitious metals from crystallisation are included), or about catalytic performance (ΔG_H is a
  thermodynamic first filter).

## Distribution

- Zenodo record (DOI `<reserved on upload>`): the SQLite database (xz-compressed), every table as Parquet and as
  gzipped CSV (one zip each), the structure files (extended XYZ), the result tables, a snapshot of the
  code, this datasheet, `SCHEMA.md`, checksums. Rebuild: `python db/build.py && python db/release.py`.
- License: database, tables and results CC BY 4.0; structure files (`structures/`) CC0 1.0; code MIT. Source data
  are CC0 (wwPDB).

## Maintenance

- Versioned releases; every version is rebuilt from the pipeline (`python db/build.py`) and records the
  query date and code commit in `dataset_info`.
- Contact: `<CONTACT EMAIL>`

## Known limitations

- Labels follow EC annotations; enzymes without EC numbers are missed in the `enzyme` set.
- The `pdb30` set keeps one representative per cluster, which need not be the chain binding the metal.
- Physiological relevance of single sites is not annotated yet (planned: UniProt metal-binding annotations
  mapped through SIFTS).
- Metal sites in very recent entries may lack cluster ids.
- DFT: one spin state per structure, PBE without U, no solvent or potential dependence. Version 0.3 holds the
  runs finished at the time of the build (`dataset_info.status`, `data_as_of`); the rest follow in later versions.
- The literature check reads abstracts (PubMed, arXiv, chemRxiv via Valency); a motif reported only in a full text
  or supporting information can be missed. "Not found" means no prior report was found, not proven novelty.
- SAC models are lattice-embedded: the enzyme's donor atoms become heteroatoms at the vacancy edge (His N →
  pyridinic N, Cys S → S, carboxylate or water O → O), so no amino-acid fragments are grafted or capped;
  second-shell effects (hydrogen bonds, protein electrostatics) are not represented.
- Stability is thermodynamic only (metal formation energy against the bulk metal, dissolution potential at
  pH 0); no diffusion or clustering barriers.
