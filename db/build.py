"""Build the release database: PDB metal sites, their first shells, SAC mappings, DFT and literature results.

    python db/build.py          # -> data/release/sac_motifs.sqlite, data/release/{parquet,csv}/<table>.*, db/SCHEMA.md

Everything is rebuilt from the pipeline outputs, so the database never holds anything the scripts cannot
regenerate:
  data/complete/*.parquet + provenance.json   bio/ mode `complete`: every PDB entry with a d-block metal
  data/hits_*.csv, enzyme_sites.csv, pdb_sites_all.csv, motifs_*.csv, entities.csv   the curated subsets
  dft/batches.py, data/dgh_graphene.csv       SAC models and DFT results
  data/stability_graphene.csv                 metal formation energies and dissolution potentials
  data/novelty.csv                            literature counts
  data/xtb_scan.csv, xtb_calibration.csv, xtb_offsets.csv   GFN2-xTB screening of the graphene sites and its DFT
                                              calibration
  data/literature_check.json                  literature check of the sites near the optimum (Valency)
  data/complete/supports.parquet, motifs_supports.parquet   the same first shells mapped to MoS2 and silica
  data/xtb_supports.csv                       GFN2-xTB hydrogen binding on MoS2 and silica sites
  data/xtb_reactions.csv, xtb_orr_calibration.csv   GFN2-xTB O2, CO2 and N2 intermediates and the O calibration
The complete set is the master site table; the `enzyme` and `pdb30` sets are subsets of it (same site_id,
same mapping, always the representative copy) and are marked in `site_sets`. Column descriptions live in
COLUMNS below; they go into the table `column_docs` and into db/SCHEMA.md, so the file documents itself.
"""
import json
import os
import re
import sqlite3
import subprocess
import sys
from datetime import date, datetime

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
DATA = os.path.join(ROOT, "data")
COMPLETE = os.path.join(DATA, "complete")
OUT = os.path.join(DATA, "release")
DB = os.path.join(OUT, "sac_motifs.sqlite")
sys.path.insert(0, os.path.join(ROOT, "dft"))
from batches import GRAPHENE_BATCH1, GRAPHENE_BATCH2   # noqa: E402

VERSION = "0.4"
SETS = {"enzyme": ("enzyme_sites.csv", "hits_enzymes.csv", "motifs_enzyme.csv",
                   "HER/ORR/CO2RR/N2-fixation enzymes by EC class, <= 2.5 A, 95 % sequence-identity representatives"),
        "pdb30": ("pdb_sites_all.csv", "hits_all.csv", "motifs_pdb.csv",
                  "entries with a Fe/Ni/Co/Cu/Mo/W/Mn/V component, <= 2.0 A, 30 % sequence-identity representatives")}
COMPLETE_SET = "every PDB entry with a d-block metal component, any resolution and method, every copy of a site"
CONTACT_KINDS = {"donors": "donor", "weak": "weak", "other": "other", "metals": "metal"}
# metals of enzymes; the others (Os/Ir hexammines, Pt drugs, Hg/Au/Ta phasing agents, Cd additives, ...) are
# almost always soaks or additives. A coarse flag: physiological relevance of a single site is not checked.
BIOLOGICAL = {"V", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Mo", "W"}

SITE_DOCS = {
    "site_id": "<pdb_id>_<auth chain>_<metal residue>_<auth seq id + insertion code>_<metal atom>",
    "pdb_id": "PDB entry id", "model": "model number (only the first model is read)",
    "metal": "element of the metal",
    "metal_class": "biological (V, Mn, Fe, Co, Ni, Cu, Zn, Mo, W) or non-biological (other d-block metals: mostly "
    "soaks, drugs, phasing agents, additives); a coarse flag, not a check of the single site",
    "chain": "author chain id of the metal",
    "label_asym_id": "mmCIF label_asym_id (subchain) of the metal",
    "metal_res": "residue (chemical component) holding the metal, e.g. HEM, SF4, ZN",
    "seqid": "author residue number with insertion code", "auth_seq_id": "author residue number",
    "icode": "insertion code", "metal_atom": "atom name of the metal", "altloc": "alternative-location id read",
    "metal_entity": "entity id of the metal's residue", "x": "metal x (A, deposited frame)",
    "y": "metal y (A)", "z": "metal z (A)", "metal_occ": "occupancy of the metal",
    "metal_b": "B-factor of the metal (A^2)", "metal_b_rel": "metal B / median B of the entry",
    "copy_group": "sites that are copies of each other in one entry (same metal, residue, donor residues/atoms)",
    "is_representative": "1 for the copy with the lowest metal B in its copy group",
    "n_copies": "size of the copy group", "chains": "chains of the copy group",
    "cn": "coordination number: donors within the element cutoffs (N/O 2.6, S 2.8, Se 2.9, C 2.3, Cl 2.8 A)",
    "formula": "donor elements, ordered N, S, O, C, Se, Cl (e.g. N2S2)",
    "raw_motif": "<metal>-<formula>, the first-shell motif",
    "graphene_motif": "preferred graphene site: <metal>-<donors>(DV|SV[,isomer])[+ax:<axial donors>]",
    "preferred_site": "DV (double vacancy, 4 in-plane donors), SV (single vacancy, 3 donors), none (CN <= 2)",
    "dv_formula": "the most square-planar 4-donor subset", "dv_isomer": "cis/trans for 2+2 and 2+1+1 subsets",
    "dv_score": "RMS(cis - 90) + RMS(trans - 180) of that subset (deg); DV preferred if <= 55",
    "dv_tau4": "tau4 of that subset", "dv_classes": "donor classes of that subset",
    "axial_formula": "donors left axial in the DV mapping", "axial_classes": "their classes",
    "sv_formula": "the 3-donor face of the SV mapping", "sv_score": "RMS(angles - 90) of that face (deg)",
    "sv_face": "1 if the 3 donors form a face (all angles <= 140 deg)", "sv_classes": "donor classes of the face",
    "sv_removed": "donors left out of the face (the most labile first)",
    "tau4": "tau4 (Yang et al. 2007) for CN 4: 0 square planar, 1 tetrahedral",
    "tau5": "tau5 (Addison et al. 1984) for CN 5: 0 square pyramidal, 1 trigonal bipyramidal",
    "mean_dist": "mean metal-donor distance (A)", "donor_classes": "classes of all donors",
    "donors": "readable donor list: class chain:residue:seqid:atom(* = symmetry mate) distance",
    "n_metal_neighbors": "other metals at 2.0-3.0 A (cluster nuclearity)", "sym_donors": "donors from symmetry images",
    "resolution": "resolution of the entry (A)", "method": "experimental method of the entry",
    "host_entity": "polymer entity providing most protein donors (else the metal's chain)",
    "cluster30": "RCSB 30 % sequence cluster of the host entity (0-based line of clusters-by-entity-30.txt)",
    "cluster50": "50 % cluster", "cluster95": "95 % cluster", "cluster100": "100 % cluster",
    "reaction": "HER/ORR/CO2RR/NRR label (see label_source)", "enzyme": "enzyme label of the matched EC class",
    "query_ec": "matched EC class", "ec": "EC numbers of the host entity", "description": "host entity description",
    "organism": "host entity organism", "taxonomy_id": "NCBI taxonomy id of the host entity",
    "uniprot": "UniProt accession(s) of the host entity",
    "label_source": "entity: labels from the host entity; entry: union over the entry"}

COLUMNS = {
    "dataset_info": {"key": "property of the dataset", "value": "its value (JSON for structured values)"},
    "entries": {
        "pdb_id": "PDB entry id", "title": "entry title", "method": "experimental method (short)",
        "exptl_methods": "all experimental methods", "resolution": "best resolution (A)",
        "resolution_combined": "all resolutions reported", "r_free": "R-free (X-ray)", "r_work": "R-work (X-ray)",
        "deposit_date": "deposition date", "initial_release_date": "first release",
        "revision_date": "latest revision at download", "revision": "latest revision number at download",
        "citation_doi": "DOI of the primary citation", "citation_pubmed": "PubMed id of the primary citation",
        "citation_year": "year of the primary citation", "citation_journal": "journal abbreviation",
        "deposited_atom_count": "atoms in the deposited model", "polymer_entity_count": "polymer entities",
        "nonpolymer_entity_count": "non-polymer entities", "metals": "d-block metals with a component in the entry",
        "metal_comp_ids": "components that contain them", "file_bytes": "size of the mmCIF.gz used",
        "reaction": "reaction labels of the entry's entities", "enzyme": "enzyme labels",
        "in_enzyme_set": "1 if in the enzyme set", "in_pdb30_set": "1 if in the pdb30 set",
        "in_complete_set": "1 if in the complete set"},
    "entities": {
        "pdb_id": "PDB entry id", "entity_id": "polymer entity id", "rcsb_entity": "<pdb_id>_<entity_id>",
        "auth_asym_ids": "author chain ids", "asym_ids": "label (mmCIF) chain ids", "polymer_type": "Protein, DNA, ...",
        "description": "entity description (pdbx_description)", "ec": "EC numbers as deposited",
        "ec_lineage": "EC lineage ids (RCSB)", "reaction": "HER/ORR/CO2RR/NRR if the EC lineage matches a searched class",
        "query_ec": "matched EC class(es)", "enzyme": "enzyme label", "organism": "source organism(s)",
        "taxonomy_id": "NCBI taxonomy id(s)", "uniprot": "UniProt accession(s)",
        "cluster30": "RCSB 30 % sequence cluster (0-based line; null if not clustered)", "cluster50": "50 % cluster",
        "cluster95": "95 % cluster", "cluster100": "100 % cluster"},
    "sites": SITE_DOCS,
    "site_sets": {"site_id": "site", "set_name": "complete, enzyme or pdb30"},
    "contacts": {
        "site_id": "site", "kind": "donor (first shell), weak (not a donor, see reason), other (element without "
        "donor cutoff), metal (metal neighbour at 2.0-3.0 A)", "rank": "order by distance within its kind",
        "element": "element", "contact_class": "donor class, e.g. N_His, S_Cys, O_water, N_porphyrin, S_inorganic, "
        "S_dithiolene, C_CO/CN, C_carbide (empty for metals)", "auth_chain": "author chain id",
        "label_asym_id": "mmCIF subchain", "resname": "residue name", "seqid": "author residue number",
        "icode": "insertion code", "atom": "atom name", "distance": "distance to the metal (A)",
        "dx": "x(contact) - x(metal), symmetry applied (A)", "dy": "y difference (A)", "dz": "z difference (A)",
        "occupancy": "occupancy", "b_factor": "B-factor (A^2)",
        "symmetry_image": "0 = deposited model, k > 0 = gemmi cell image k (crystal symmetry, unapplied NCS)",
        "polymer": "1 if part of a polymer", "entity_id": "entity id",
        "reason": "why a weak/other contact is not a donor: beyond_cutoff, minor_occupancy, alt_disorder, "
        "alt_overlap, protein_carbon, his_second_ring_n, element_not_counted ('+'-joined)",
        "neighbor_site_id": "for metals: site_id of the neighbouring metal (asymmetric-unit copy)"},
    "motifs": {
        "set_name": "site set", "kind": "graphene or first_shell", "motif": "motif string", "metal": "metal",
        "n_sites": "sites", "n_entries": "entries", "n_clusters30": "30 % clusters",
        "n_sites_rep": "representative sites (complete set)", "n_entries_rep": "entries with a representative site",
        "n_clusters30_rep": "clusters, representative sites", "n_xray_entries": "X-ray entries",
        "n_em_entries": "electron-microscopy entries", "reactions": "reaction labels among the sites",
        "example_descriptions": "host-entity descriptions", "example_pdb": "best-resolution entries (resolution)",
        "best_resolution": "A", "median_tau4": "median tau4", "median_dv_tau4": "median tau4 of the DV subset",
        "median_dv_score": "deg", "median_sv_score": "deg", "frac_mononuclear": "share without a metal neighbour",
        "donor_classes_top": "most common donor-class combination"},
    "processing_log": {
        "pdb_id": "PDB entry id", "status": "ok / failed", "error": "exception if failed",
        "seconds": "extraction time", "n_sites": "sites found", "n_representative": "representative sites",
        "file_revision": "revision of the mmCIF used", "file_revision_date": "its date", "n_models": "models in file",
        "crystal_symmetry": "1 if searched with crystal symmetry", "n_atoms_model1": "atoms in model 1",
        "n_metal_atoms_file": "metal atoms in the file (all conformers)", "n_skipped_altloc": "not the selected conformer",
        "n_metal_atoms": "metal atoms considered", "n_skipped_zero_occupancy": "skipped: occupancy 0",
        "n_skipped_minor_alternative": "skipped: alternative position of another metal",
        "crystal_frame_standardized": "1 if a non-standard crystal frame was standardised for the search"},
    "sac_models": {
        "model_id": "<support>:<motif as built>", "support": "graphene (later MoS2, SiO2)",
        "motif": "as built by materials/graphene/sites.py, e.g. Ni-N2S2:DV:cis5",
        "sweep_motif": "the same site in the notation of sites.graphene_motif (without axial ligands)",
        "metal": "metal", "donors": "in-plane donor formula", "site_type": "DV or SV",
        "isomer": "cis5/cis6/trans or explicit d1..d4 order", "axial": "axial ligand below the sheet",
        "batch": "DFT batch (graphene3 = picked by the xTB scan)", "origin": "enzyme site the motif was taken from, or the xTB pick"},
    "dft_methods": {"method_id": "id", "description": "settings (JSON)"},
    "dft_results": {
        "model_id": "SAC model", "state": "bare, H_M (H on the metal), H_S (H on an S donor)",
        "method_id": "settings", "energy_eV": "total energy of the 3x3x1 SCF (eV)",
        "magnetization": "total magnetization (muB/cell)", "dE_H": "E(*H) - E(*) - E(H2)/2 (eV)",
        "dG_H": "dE_H + 0.24 eV", "H_partner": "atom the H ended up on", "H_dist": "its distance (A)",
        "max_M_X_bare": "largest metal-donor distance of the relaxed bare site (A)",
        "M_height_bare": "metal height above the sheet (A)"},
    "dft_stability": {
        "model_id": "SAC model", "method_id": "settings",
        "bare_ref": "bare configuration used (the lowest found, e.g. bare_Sdown)",
        "E_bare": "energy of the site with the metal (eV)",
        "E_site": "energy of the same doped vacancy without the metal (eV)",
        "E_bulk_atom": "energy per atom of the relaxed bulk metal, dense k mesh (eV)",
        "E_f": "E_bare - E_site - E_bulk_atom (eV); < 0: the single atom is more stable in the site than in "
               "its bulk metal",
        "U_diss": "dissolution potential U0(M/M^n+) - E_f/n (V vs SHE); > 0: stable under HER conditions"},
    "dft_mos2": {
        "run": "H on 6x6 MoS2: on a basal S (31, clean) or on an Fe bound to "
               "three S (34 clean, 36 charged against a plate)",
        "charge_e": "extra electrons per cell from the plate (QE gate; < 0: holes, acceptor-like), empty without",
        "host": "relaxed host: clean MoS2 from an earlier MoS2 project with the same settings, or computed here (Fe sites)",
        "E_eV": "final energy of the Gamma relaxation (eV; gate runs: without the QE 7.1 plate constant)",
        "magnetization": "total magnetization (muB/cell)", "dE_H": "E(host+H) - E(host) - E(H2)/2 (eV)",
        "dG_H": "dE_H + 0.24 eV"},
    "literature": {
        "motif": "metal-donors", "names": "name variants searched",
        "openalex_all": "OpenAlex works matching a name AND a single-atom term",
        "openalex_graphene": "... AND graphene / N-doped carbon", "openalex_HER": "... AND hydrogen evolution",
        "europepmc_all": "Europe PMC hits, same query", "europepmc_graphene": "", "europepmc_HER": "",
        "top_hits": "first OpenAlex hits: year title (DOI)"},
    "xtb_scan": {
        "model_id": "SAC model (graphene:<motif>)", "motif": "as built by materials/graphene/sites.py",
        "source": "scan (enzyme donor environments with the nine enzyme metals, and the DFT motifs) or library (every "
                  "N/S/O/C arrangement, screen/xtb_scan.py library)",
        "state": "H_M (H on the metal) or H_S (H on an S donor)", "method_id": "settings",
        "dG_H": "dE_H + 0.24 eV from GFN2-xTB (eV)", "dE_H": "E(*H) - E(*) - E(H2)/2 (eV)",
        "dG_H_calibrated": "dG_H plus one offset per metal (H on the metal) or one offset for H on S (xtb_offsets: "
                           "fitted on the sites that also have DFT values); empty for metals with fewer than two "
                           "DFT pairs",
        "H_partner": "atom the H ended up on", "H_dist": "its distance (A)",
        "max_M_X": "largest metal-donor distance with H (A)", "max_M_X_bare": "the same for the bare site (A)",
        "M_height": "metal height above the sheet (A)", "converged": "relaxation with H converged",
        "bare_converged": "bare relaxation converged", "steps": "BFGS steps", "seconds": "wall time (s)",
        "error": "error message if the relaxation failed"},
    "xtb_offsets": {
        "state": "H_M (H on the metal) or H_S (H on an S donor; one offset for every metal)",
        "metal": "metal (empty for H_S)", "n_pairs": "xTB/DFT pairs the offset is fitted on",
        "offset": "mean of dG_H(DFT) - dG_H(xTB) over those pairs (eV)",
        "loo_mae": "leave-one-out mean absolute error: each pair predicted from the offset of the others (eV)",
        "used": "1 if the offset is applied (at least two pairs)"},
    "xtb_calibration": {
        "motif_dft": "SAC model as in dft_results", "state_dft": "H position", "dG_H_dft": "DFT value (eV)",
        "dG_H_xtb": "GFN2-xTB value of the same structure (eV)", "dG_H_cal": "xTB after the global linear fit (eV)",
        "H_partner_dft": "H partner in DFT", "H_partner_xtb": "H partner in xTB"},
    "literature_check": {
        "motif": "metal-donors (donor order N, S, O, C)", "model_id": "SAC model with that first shell",
        "dG_H": "DFT dG_H of that model (eV)",
        "verdict": "reported_experimental, reported_computational_only, related_only (same elements, other donor "
                   "count or support) or not_found: same metal and exactly these donors on a carbon support",
        "confidence": "low, medium, high", "her_reported": "whether hydrogen evolution was studied for this motif",
        "summary": "what the reader found", "challenge_refuted": "whether the skeptic overturned the verdict",
        "challenge_notes": "what the skeptic searched and found",
        "method": "how the check was done (abstract level)"},
    "literature_papers": {
        "motif": "as in literature_check", "role": "closest (reader) or counter_evidence (skeptic)",
        "title": "paper title", "year": "publication year", "doi": "DOI (all registered at doi.org, checked)",
        "note": "coordination the abstract states, or why it does not match", "match": "exact or related"},
    "validation_metalc": {
        "site_id": "site", "annotated": "1 if the PDB annotates metal links (struct_conn metalc) for this metal atom",
        "n_ours": "our donors from other residues (comparable with the annotation)",
        "n_annotated": "annotated partners (non-metal, unique)", "n_both": "in both",
        "n_missed": "annotated partners we do not count as donors", "n_extra": "our donors without an annotated link",
        "n_intra": "donors inside the metal's own component (never in struct_conn)",
        "exact": "1 if our inter-residue donors equal the annotated partners",
        "missed_reasons": "why each missed partner is not a donor for us (weak contact reason, no cutoff, > 3.4 A)",
        "extra_classes": "donor classes of our unannotated donors"},
    "support_sites": {
        "site_id": "site", "support": "MoS2 (metal on three lattice S, 3.16 A apart), MoS2vS (metal in an S vacancy, "
                                      "on three Mo; only metal-cluster corners map) or SiO2 (siloxide O from Si-OH "
                                      "condensation, 2.6-3.2 A apart)",
        "face_size": "donors (MoS2vS: metal neighbours) the support supplies: 3, else 2, else 1",
        "face_classes": "donor classes of the face (MoS2vS: elements of the three metal neighbours)",
        "score": "RMS deviation of the face's X...X distances from the support's (A); 0 for one donor",
        "xx_mean": "mean X...X distance on the face (A)",
        "apical": "formula of the donors that stay on the metal above the surface (a thiol, a silane linker or "
                  "the solvent has to supply them)",
        "apical_classes": "their donor classes",
        "removed": "donors that would sit in the surface (MoS2vS: non-S donors in the place of the S ring)",
        "ring": "MoS2vS only: donors at the height of the S ring around the vacancy",
        "site_motif": "the bare single-atom site, e.g. Fe-S3(MoS2)", "motif": "site_motif + apical donors, e.g. "
                                                                         "Fe-S3(MoS2)+ax:S"},
    "support_motifs": {
        "set_name": "complete (all sites) or enzyme (sites of HER/ORR/CO2RR/NRR enzymes)", "support": "as in "
        "support_sites", "motif": "as in support_sites", "site_motif": "the bare single-atom site",
        "metal": "metal", "n_sites": "sites, every copy", "n_entries": "entries", "n_clusters30": "30 % sequence "
        "clusters (or entry + host entity)", "n_sites_rep": "representative sites", "n_entries_rep": "entries with "
        "a representative site", "n_clusters30_rep": "clusters with a representative site",
        "reactions": "reaction labels of the sites", "example_descriptions": "polymer descriptions, best "
        "resolution first", "example_pdb": "entries (resolution, A)", "best_resolution": "A",
        "n_sites_hq": "representative X-ray sites at <= 2.0 A (the geometry columns use only these)",
        "median_score": "median score (A)", "median_xx": "median X...X on the face (A)",
        "frac_fit": "share with score <= 0.25 A (MoS2, MoS2vS) or <= 0.10 A (SiO2)",
        "face_classes_top": "most common face", "removed_top": "most common removed donors",
        "graphene_motif_top": "most common graphene mapping of the same sites", "frac_mononuclear": "share "
        "without metal neighbours"},
    "xtb_supports": {
        "support": "MoS2 (3x3 cut of the DFT sheet, metal on three S over a Mo or a hollow) or SiO2 (metal in the "
                   "open corner of Me7Si7O9(OH)3: two siloxides and one silanol)",
        "metal": "metal (clean = MoS2 without metal, H on a basal S)", "state": "H_M (H on the metal), H_S (H on "
        "S, clean MoS2) or H_O (H on a siloxide O)", "site_bare": "the lower bare site, the reference",
        "site_H": "site of the lower structure with H", "dE_H": "E(*H) - E(*) - E(H2)/2 (eV)",
        "dG_H": "dE_H + 0.24 eV, GFN2-xTB (eV)", "dG_H_offset": "dG_H plus the per-metal offset of xtb_calibration "
        "(graphene, H on the metal); assumes xTB's error belongs to the metal, not the support (eV)",
        "H_partner": "atom the H ended up on", "H_dist": "its distance (A)",
        "M_X_bare": "the four shortest metal-S (MoS2) or metal-O (SiO2) distances, bare (A)",
        "M_X_H": "the same with H (A)", "M_height_bare": "MoS2: metal height above the top S layer (A)",
        "converged": "both relaxations converged", "seconds": "wall time of the structures (s)",
        "dG_H_dft": "DFT value of the same site where there is one (eV)"},
    "xtb_reactions": {
        "kind": "graphene, MoS2 or SiO2", "site": "graphene motif as built (materials/graphene/sites.py) or the metal",
        "metal": "metal", "calibration_only": "1 = a double-vacancy site of Cho et al. 2022 computed only for the O "
        "calibration",
        "dG_H": "dG_H with H on the metal, for selectivity against H2 evolution (eV): DFT where computed, else "
                "GFN2-xTB with its per-metal offset, else raw xTB",
        "dG_H_source": "DFT, xTB calibrated or xTB raw",
        "U_ORR_cal": "limiting potential of O2 -> H2O (4 e-) from the calibrated O intermediates (V vs RHE)",
        "U_OER_cal": "limiting potential of O2 evolution, calibrated (V vs RHE; equilibrium 1.23 V)",
        "U_H2O2_cal": "limiting potential of O2 -> H2O2 (2 e-), calibrated (V vs RHE; equilibrium 0.70 V)",
        "U_CO": "limiting potential of CO2 -> CO via *COOH and *CO (V vs RHE; equilibrium -0.10 V), uncalibrated",
        "dG_CO_release": "*CO -> CO(g) (eV); above ~0.5 eV CO poisons the site",
        "U_HCOOH": "limiting potential of CO2 -> HCOOH via *OCHO (V vs RHE; equilibrium -0.12 V), uncalibrated",
        "U_NRR": "limiting potential of N2 -> NH3 from the first (*N2 -> *N2H) and last (*NH2 -> *NH3) protonation "
                 "(V vs RHE), uncalibrated",
        "dG_N2_N2H": "first protonation (eV)", "dG_NH2_NH3": "last protonation (eV)",
        "U_ORR": "as U_ORR_cal from raw xTB", "U_OER": "raw xTB", "U_H2O2": "raw xTB",
        "dE_OH": "E(*OH) - E(*) - E(H2O) + E(H2)/2, xTB (eV; the scale of the calibration data)",
        "dE_O": "E(*O) - E(*) - E(H2O) + E(H2) (eV)", "dE_OOH": "E(*OOH) - E(*) - 2 E(H2O) + 3/2 E(H2) (eV)",
        "dG_OH": "free energy of *OH vs H2O and H2 (eV)", "dG_O": "*O (eV)", "dG_OOH": "*OOH (eV)",
        "dG_OH_cal": "*OH after the calibration (eV)", "dG_O_cal": "*O, calibrated (eV)",
        "dG_OOH_cal": "*OOH, calibrated (eV)",
        "dG_COOH": "CO2 + H+ + e- -> *COOH at 0 V vs RHE (eV)", "dG_CO": "CO(g) -> *CO (eV)",
        "dG_OCHO": "CO2 + H+ + e- -> *OCHO (eV)", "dG_N2": "N2 -> *N2, end-on (eV)",
        "dG_N2H": "N2 + H+ + e- -> *N2H (eV)", "dG_NH2": "1/2 N2 + 2 (H+ + e-) -> *NH2 (eV)",
        "dG_NH3": "1/2 N2 + 3 (H+ + e-) -> *NH3 (eV)"},
    "xtb_orr_calibration": {
        "site": "double-vacancy site (Cho et al., Front. Chem. 2022, Catalysis-Hub AraComputational2022)",
        "adsorbate": "OH, O or OOH", "dE_dft": "their PBE-D3 reaction energy vs H2O and H2 (eV)",
        "dE_xtb": "ours from GFN2-xTB, same definition (eV)"},
    "column_docs": {"table_name": "table", "column_name": "column", "description": "meaning, units"},
}

DFT_METHOD = {"code": "Quantum ESPRESSO 7.1", "functional": "PBE", "dispersion": "D3(BJ)",
              "pseudopotentials": "GBRV ultrasoft v1.5", "ecutwfc_Ry": 40, "ecutrho_Ry": 200,
              "smearing": "Marzari-Vanderbilt 0.01 Ry", "spin": "polarised, starting moment on the metal",
              "cell": "6x6 graphene (a = 2.467 A), 20 A height, dipole correction",
              "relaxation": "Gamma point, forces < 1e-3 Ry/bohr", "energies": "SCF with 3x3x1 k-points",
              "dG_H": "dE_H + 0.24 eV (Norskov et al. 2005)"}
ORDER = ["N", "S", "O", "C", "Se", "Cl"]     # donor order of bio/first_shell.py formulas


def read(name):
    path = os.path.join(DATA, name)
    if not os.path.exists(path):
        return None
    return pd.read_csv(path, dtype={"entry_id": str, "entity_id": str, "host_entity": str, "cluster30": str,
                                    "seqid": str}, keep_default_na=False, na_values=[""])


def site_id(r):
    return "%s_%s_%s_%s_%s" % (r.entry_id, r.chain, r.metal_res, r.seqid, r.metal_atom)


def contacts(sites):
    """One row per donor, weak contact, other atom and metal neighbour, from the shell_json records."""
    rows = []
    for sid, js in zip(sites.site_id, sites.shell_json):
        rec = json.loads(js)
        for key, kind in CONTACT_KINDS.items():
            for k, c in enumerate(rec.get(key) or [], 1):
                rows.append((sid, kind, k, c["el"], c.get("cls", ""), c["auth_chain"], c["label_asym"], c["resname"],
                             c["seqid"], c["icode"], c["atom"], c["dist"], c["dx"], c["dy"], c["dz"], c["occ"], c["b"],
                             c["sym"], int(bool(c["polymer"])), c["entity"], c.get("reason", ""), c.get("site_id")))
    return pd.DataFrame(rows, columns=list(COLUMNS["contacts"]))


def canonical(donors):
    """'C2S2' -> 'S2C2': the element order the sweep uses."""
    n = {el: int(k or 1) for el, k in re.findall(r"([A-Z][a-z]?)(\d*)", donors)}
    return "".join(el + (str(n[el]) if n[el] > 1 else "") for el in ORDER if el in n)


def sweep_isomer(iso):
    """Isomer of a built DV site in the sweep's terms: cis5/cis6 -> cis; an explicit d1..d4 order (e.g. SNSO)
    is cis when its repeated donor sits on a pentagon (d1, d3) or hexagon (d1, d2) pair, trans on d1, d4."""
    if iso in ("cis5", "cis6", "trans"):
        return "," + iso.rstrip("56")
    order = re.findall(r"[A-Z][a-z]?", iso or "")
    if len(order) == 4:
        pairs = [(i, j) for i in range(4) for j in range(i + 1, 4) if order[i] == order[j]]
        if len(pairs) == 1:
            return ",trans" if pairs[0] in ((0, 3), (1, 2)) else ",cis"
    return ""


def sac_models():
    """Every graphene model that has DFT runs: batches 1 and 2 (dft/batches.py) and batch 3, the xTB picks."""
    rows = []
    sel = read("xtb_selected.csv")
    xtb = [(m, "xTB scan: calibrated dG_H %+.2f eV (%s)" % (c, st)) for m, c, st in
           zip(sel.motif, sel.dG_H_cal, sel.state)] if sel is not None else []
    for batch, motifs in (("graphene1", GRAPHENE_BATCH1), ("graphene2", GRAPHENE_BATCH2), ("graphene3", xtb)):
        for motif, origin in motifs:
            head, *opts = motif.split(":")
            metal, donors = head.split("-")
            site = opts[0]
            iso = next((o for o in opts[1:] if not o.startswith("ax=")), "")
            ax = next((o[3:] for o in opts[1:] if o.startswith("ax=")), "")
            rows.append({"model_id": "graphene:" + motif, "support": "graphene", "motif": motif,
                         "sweep_motif": "%s-%s(%s%s)" % (metal, canonical(donors), site, sweep_isomer(iso)),
                         "metal": metal, "donors": donors, "site_type": site, "isomer": iso, "axial": ax,
                         "batch": batch, "origin": origin})
    return pd.DataFrame(rows)


XTB_METHOD = {"code": "GFN2-xTB via tblite 0.7 (ASE calculator)", "cell": "6x6 graphene, c = 20 A, periodic, Gamma point",
              "electrons": "spin-restricted, Fermi smearing 1000 K", "relaxation": "BFGS to 0.05 eV/A, at most 400 "
              "steps, one far carbon fixed", "dG_H": "dE_H + 0.24 eV, H2 relaxed in the same cell"}
XTB_SUPPORTS_METHOD = {"code": "GFN2-xTB via tblite 0.7 (ASE calculator)",
                       "MoS2": "3x3 cut of the relaxed 6x6 DFT sheet (a = 3.16 A), c = 19 A, periodic, Gamma point "
                               "(K folded in), Mo layer fixed; metal on three S over a Mo and over a hollow",
                       "SiO2": "Me7Si7O9(O)2(OH)M molecule, no periodic cell", "electrons": "spin-restricted, Fermi "
                       "smearing 1000 K", "relaxation": "BFGS to 0.05 eV/A, at most 400 steps",
                       "dG_H": "dE_H + 0.24 eV, H2 relaxed in the same cell (MoS2) or in vacuum (SiO2)"}


def xtb_tables():
    """GFN2-xTB scan with the per-metal calibrated value, and the xTB-DFT pairs it was calibrated on."""
    x, pairs = read("xtb_scan.csv"), read("xtb_calibration.csv")
    if x is None or pairs is None:
        return pd.DataFrame(columns=list(COLUMNS["xtb_scan"])), pd.DataFrame(columns=list(COLUMNS["xtb_calibration"]))
    o = read("xtb_offsets.csv")                # screen/calibrate.py: metals with enough DFT pairs
    o = o[o.used]
    offset = o[o.state == "H_M"].set_index("metal").offset
    offset_s = o[o.state == "H_S"].offset.iloc[0]
    metal = x.motif.str.split("-").str[0]
    x = x.assign(model_id="graphene:" + x.motif, method_id="gfn2-xtb-6x6",
                 dG_H_calibrated=x.dG_H + pd.Series([offset_s if s == "H_S" else offset.get(m) for m, s in
                                                     zip(metal, x.state)], index=x.index, dtype=float))
    return x.reindex(columns=list(COLUMNS["xtb_scan"])), pairs.reindex(columns=list(COLUMNS["xtb_calibration"]))


def support_tables():
    """bio/supports.py outputs -> support_sites (one row per site and support), support_motifs; and the xTB
    hydrogen binding on those supports."""
    path = os.path.join(COMPLETE, "supports.parquet")
    if not os.path.exists(path):
        return tuple(pd.DataFrame(columns=list(COLUMNS[n])) for n in ("support_sites", "support_motifs",
                                                                       "xtb_supports"))
    s = pd.read_parquet(path)
    parts = []
    for sup in ("MoS2", "MoS2vS", "SiO2"):
        cols = {c: c[len(sup) + 1:] for c in s.columns if c.startswith(sup + "_") and
                not (sup == "MoS2" and c.startswith("MoS2vS_"))}
        part = s[["site_id", *cols]].rename(columns=cols)
        part = part[part["motif"].notna()].assign(support=sup)
        parts.append(part.rename(columns={"k": "face_size", "site": "site_motif"}))
    sites = pd.concat(parts, ignore_index=True)
    sites["face_size"] = sites.face_size.astype("Int64")
    motifs = pd.read_parquet(os.path.join(COMPLETE, "motifs_supports.parquet")).rename(
        columns={"kind": "support", "site": "site_motif"})
    x = read("xtb_supports.csv")
    x = pd.DataFrame(columns=list(COLUMNS["xtb_supports"])) if x is None else x
    return (sites.reindex(columns=list(COLUMNS["support_sites"])),
            motifs.reindex(columns=list(COLUMNS["support_motifs"])),
            x.reindex(columns=list(COLUMNS["xtb_supports"])))


def literature_check(dft):
    """data/literature_check.json -> one row per motif and one per cited paper."""
    path = os.path.join(DATA, "literature_check.json")
    if not os.path.exists(path):
        return (pd.DataFrame(columns=list(COLUMNS["literature_check"])),
                pd.DataFrame(columns=list(COLUMNS["literature_papers"])))
    order = ["N", "S", "O", "C", "Se", "Cl"]

    def key(head):
        metal, donors = head.split("-")
        n = {el: int(k or 1) for el, k in re.findall(r"([A-Z][a-z]?)(\d*)", donors)}
        return metal + "-" + "".join(el + (str(n[el]) if n[el] > 1 else "") for el in order if el in n)

    best = (dft.dropna(subset=["dG_H"]).assign(key=lambda d: d.model_id.str[len("graphene:"):].str.split(":").str[0].map(key),
                                                absdg=lambda d: d.dG_H.abs())
            .sort_values("absdg").drop_duplicates("key").set_index("key"))
    method = ("Valency (PubMed, arXiv, chemRxiv), 21 Oct 2026: one search agent (titles), one reader (at most 8 abstracts) "
              "and one skeptic per motif; then a deep check of the closest papers (open-access full texts, citation graph, "
              "series notation such as M-NxO4-x) with a cross-motif critic; paywalled full texts and SI not read")
    rows, papers = [], []
    for r in json.load(open(path)):
        v, c = r["verdict"], r.get("challenge") or {}
        k = key(r["motif"])
        rows.append({"motif": k, "model_id": best.model_id.get(k), "dG_H": best.dG_H.get(k), "verdict": v["verdict"],
                     "confidence": v["confidence"], "her_reported": v["her_reported"], "summary": v["summary"],
                     "challenge_refuted": c.get("refuted"), "challenge_notes": c.get("notes"), "method": method})
        for p in v["closest"]:
            papers.append({"motif": k, "role": "closest", "title": p["title"], "year": p.get("year"), "doi": p.get("doi"),
                           "note": p.get("coordination"), "match": p.get("match")})
        for p in c.get("counter_evidence", []):
            papers.append({"motif": k, "role": "counter_evidence", "title": p["title"], "year": p.get("year"),
                           "doi": p.get("doi"), "note": p.get("why"), "match": "related"})
    return (pd.DataFrame(rows).reindex(columns=list(COLUMNS["literature_check"])),
            pd.DataFrame(papers).reindex(columns=list(COLUMNS["literature_papers"])))


def file_time(name):
    path = os.path.join(DATA, name)
    return datetime.fromtimestamp(os.path.getmtime(path)).strftime("%Y-%m-%d %H:%M") if os.path.exists(path) else None


def git_commit():
    try:
        sha = subprocess.check_output(["git", "-C", ROOT, "rev-parse", "HEAD"], text=True).strip()
        dirty = subprocess.call(["git", "-C", ROOT, "diff", "--quiet", "HEAD"]) != 0
        return sha + ("+uncommitted" if dirty else "")
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def validation_summary(v):
    """Agreement with the PDB's metal-link annotation and threshold sensitivity (bio/validate.py), as JSON."""
    a = v[v.annotated == True]   # noqa: E712 (nullable booleans)
    out = {"sites_annotated": int(len(a)), "sites": int(len(v)), "exact_agreement": float(a.exact.mean()),
           "precision": float(a.n_both.sum() / a.n_ours.sum()), "recall": float(a.n_both.sum() / a.n_annotated.sum())}
    cut = os.path.join(DATA, "validation", "cutoff_sensitivity.csv")
    if os.path.exists(cut):
        out["cutoff_sensitivity"] = pd.read_csv(cut).to_dict("records")
    return json.dumps(out) if len(v) else ""


def load_sites():
    """sites (complete set, without the JSON), contacts, site_sets."""
    s = pd.read_parquet(os.path.join(COMPLETE, "sites.parquet"))
    cont = contacts(s)
    s = s.rename(columns={"entry_id": "pdb_id"})
    s["metal_class"] = s.metal.map(lambda m: "biological" if m in BIOLOGICAL else "non-biological")
    sets = [pd.DataFrame({"site_id": s.site_id, "set_name": "complete"})]
    for name, (sites_file, _, _, _) in SETS.items():
        sub = read(sites_file)
        ids = pd.Series([site_id(r) for r in sub.itertuples()])
        missing = (~ids.isin(s.site_id)).sum()
        if missing:
            raise SystemExit("%d %s sites are not in the complete set: rebuild bio/ first" % (missing, name))
        sets.append(pd.DataFrame({"site_id": ids, "set_name": name}))
    return s[list(SITE_DOCS)], cont, pd.concat(sets).drop_duplicates()


def load_entries():
    e = pd.read_parquet(os.path.join(COMPLETE, "entries.parquet")).rename(columns={"entry_id": "pdb_id"})
    e["in_complete_set"] = 1
    for name, (_, hits_file, _, _) in SETS.items():
        h = read(hits_file).rename(columns={"entry_id": "pdb_id"})
        flag = "in_%s_set" % name
        e[flag] = e.pdb_id.isin(h.pdb_id).astype(int)
        extra = h[~h.pdb_id.isin(e.pdb_id)].drop_duplicates("pdb_id")        # entries without any metal
        if len(extra):
            extra = extra[["pdb_id", "title", "method", "resolution"]].assign(**{flag: 1, "in_complete_set": 0})
            e = pd.concat([e, extra], ignore_index=True)
    e = e.fillna({"in_enzyme_set": 0, "in_pdb30_set": 0})
    return e.reindex(columns=list(COLUMNS["entries"]))


def load_entities():
    c = pd.read_parquet(os.path.join(COMPLETE, "entities.parquet")).rename(columns={"entry_id": "pdb_id"})
    old = read("entities.csv").rename(columns={"entry_id": "pdb_id"})
    key = c.pdb_id + "_" + c.entity_id.astype(str)
    old = old[~(old.pdb_id + "_" + old.entity_id.astype(str)).isin(key)].copy()
    # the older table names an unclustered entity by its own id; here that is null, as in the complete set
    old["cluster30"] = pd.to_numeric(old.cluster30, errors="coerce").astype("Int64")
    return pd.concat([c, old], ignore_index=True).reindex(columns=list(COLUMNS["entities"]))


def load_motifs():
    parts = [pd.read_parquet(os.path.join(COMPLETE, "motifs.parquet")).assign(set_name="complete")]
    for name, (_, _, motifs_file, _) in SETS.items():
        m = read(motifs_file).rename(columns={"n_sequence_clusters": "n_clusters30",
                                              "example_enzymes": "example_descriptions"})
        parts.append(m.assign(set_name=name))
    return pd.concat(parts, ignore_index=True).reindex(columns=list(COLUMNS["motifs"]))


def main():
    os.makedirs(OUT, exist_ok=True)
    prov = json.load(open(os.path.join(COMPLETE, "provenance.json")))
    tables = {}
    tables["sites"], tables["contacts"], tables["site_sets"] = load_sites()
    tables["entries"] = load_entries()
    tables["entities"] = load_entities()
    tables["motifs"] = load_motifs()
    log = pd.read_parquet(os.path.join(COMPLETE, "entry_log.parquet")).rename(columns={"entry_id": "pdb_id"})
    tables["processing_log"] = log.reindex(columns=list(COLUMNS["processing_log"]))
    tables["sac_models"] = sac_models()
    tables["dft_methods"] = pd.DataFrame([{"method_id": "gbrv-pbe-d3-6x6", "description": json.dumps(DFT_METHOD)}])
    dgh = read("dgh_graphene.csv")
    if dgh is not None and len(dgh):
        dgh = dgh.assign(model_id="graphene:" + dgh.motif, method_id="gbrv-pbe-d3-6x6")
        dgh = dgh.rename(columns={"E_H": "energy_eV", "mag_H": "magnetization"})
        tables["dft_results"] = dgh.reindex(columns=list(COLUMNS["dft_results"]))
    else:
        tables["dft_results"] = pd.DataFrame(columns=list(COLUMNS["dft_results"]))
    stab = read("stability_graphene.csv")
    if stab is not None and len(stab):
        stab = stab.assign(model_id="graphene:" + stab.motif, method_id="gbrv-pbe-d3-6x6")
    tables["dft_stability"] = (stab if stab is not None else pd.DataFrame()).reindex(
        columns=list(COLUMNS["dft_stability"]))
    mos2 = read("dgh_mos2.csv")
    tables["dft_mos2"] = (mos2 if mos2 is not None else pd.DataFrame()).reindex(columns=list(COLUMNS["dft_mos2"]))
    val = os.path.join(DATA, "validation", "metalc_sites.parquet")
    tables["validation_metalc"] = (pd.read_parquet(val) if os.path.exists(val) else pd.DataFrame()).reindex(
        columns=list(COLUMNS["validation_metalc"]))
    lit = read("novelty.csv")
    tables["literature"] = (lit if lit is not None else pd.DataFrame()).reindex(columns=list(COLUMNS["literature"]))
    tables["xtb_scan"], tables["xtb_calibration"] = xtb_tables()
    tables["dft_methods"] = pd.concat([tables["dft_methods"], pd.DataFrame(
        [{"method_id": "gfn2-xtb-6x6", "description": json.dumps(XTB_METHOD)}])], ignore_index=True)
    tables["literature_check"], tables["literature_papers"] = literature_check(tables["dft_results"])
    tables["support_sites"], tables["support_motifs"], tables["xtb_supports"] = support_tables()
    for name in ("xtb_reactions", "xtb_orr_calibration", "xtb_offsets"):
        t = read(name + ".csv")
        t = pd.DataFrame() if t is None else t
        tables[name] = t.reindex(columns=[c for c in COLUMNS[name] if c in t.columns or t.empty])
    tables["dft_methods"] = pd.concat([tables["dft_methods"], pd.DataFrame(
        [{"method_id": "gfn2-xtb-supports", "description": json.dumps(XTB_SUPPORTS_METHOD)}])], ignore_index=True)
    ver = prov["search"]["rcsb_data_version_before"]
    tables["dataset_info"] = pd.DataFrame(list({
        "title": "Transition-metal first shells of the PDB mapped to single-atom catalyst sites (working title)",
        "version": VERSION, "build_date": date.today().isoformat(), "code_commit": git_commit(),
        "pdb_data_version": json.dumps(ver),
        "site_sets": json.dumps(dict({"complete": COMPLETE_SET}, **{k: v[3] for k, v in SETS.items()})),
        "counts": json.dumps({"entries": len(tables["entries"]), "sites": len(tables["sites"]),
                              "representative_sites": int(tables["sites"].is_representative.sum()),
                              "contacts": len(tables["contacts"])}),
        "source_data_license": "PDB archive data: CC0 1.0 (wwPDB)",
        "license": "CC BY 4.0 for the derived data (proposed; to be confirmed by the authors before release)",
        "data_as_of": json.dumps({"dft_results": file_time("dgh_graphene.csv"), "xtb_scan": file_time("xtb_scan.csv"),
                                  "literature_check": file_time("literature_check.json"),
                                  "support_sweep": file_time(os.path.join("complete", "supports.parquet")),
                                  "xtb_supports": file_time("xtb_supports.csv")}),
        "provenance_complete_set": json.dumps(prov),
        "technical_validation": validation_summary(tables["validation_metalc"]),
        "status": "hackathon release: DFT for %d of %d H positions finished, metal formation energies for %d "
                  "sites; values are drafts that need expert checking" % (
                      int(tables["dft_results"].dG_H.notna().sum()), len(tables["dft_results"]),
                      int(tables["dft_stability"].E_f.notna().sum()) if "dft_stability" in tables else 0),
    }.items()), columns=["key", "value"])
    tables["column_docs"] = pd.DataFrame([(t, c, d) for t, cols in COLUMNS.items() for c, d in cols.items()],
                                         columns=["table_name", "column_name", "description"])

    if os.path.exists(DB):
        os.remove(DB)
    con = sqlite3.connect(DB)
    for name, df in tables.items():
        df.to_sql(name, con, index=False, chunksize=50000)
    con.executescript("""
        CREATE UNIQUE INDEX ix_sites ON sites(site_id);
        CREATE INDEX ix_sites_pdb ON sites(pdb_id);
        CREATE INDEX ix_sites_motif ON sites(graphene_motif);
        CREATE INDEX ix_sites_raw ON sites(raw_motif);
        CREATE INDEX ix_contacts ON contacts(site_id, kind);
        CREATE INDEX ix_sets ON site_sets(site_id, set_name);
        CREATE INDEX ix_support_sites ON support_sites(site_id, support);
        CREATE UNIQUE INDEX ix_entries ON entries(pdb_id);
        CREATE INDEX ix_entities ON entities(pdb_id, entity_id);
        CREATE VIEW donors AS SELECT * FROM contacts WHERE kind = 'donor';
        CREATE VIEW weak_contacts AS SELECT * FROM contacts WHERE kind = 'weak';
        CREATE VIEW metal_neighbors AS SELECT * FROM contacts WHERE kind = 'metal';
        CREATE VIEW v_sites AS
            SELECT s.*, e.title, e.initial_release_date, e.citation_doi,
                   (SELECT group_concat(set_name) FROM site_sets ss WHERE ss.site_id = s.site_id) AS sets
            FROM sites s JOIN entries e USING (pdb_id);
        CREATE VIEW v_sac_overview AS
            SELECT m.model_id, m.motif, m.sweep_motif, m.origin,
                   (SELECT COUNT(*) FROM sites s JOIN site_sets ss USING (site_id)
                     WHERE ss.set_name = 'enzyme' AND s.graphene_motif LIKE m.sweep_motif || '%') AS enzyme_sites,
                   (SELECT COUNT(*) FROM sites s WHERE s.is_representative = 1
                     AND s.graphene_motif LIKE m.sweep_motif || '%') AS pdb_sites_representative,
                   (SELECT COUNT(DISTINCT s.pdb_id) FROM sites s WHERE s.is_representative = 1
                     AND s.graphene_motif LIKE m.sweep_motif || '%') AS pdb_entries,
                   (SELECT MIN(d.dG_H) FROM dft_results d WHERE d.model_id = m.model_id) AS dG_H_min,
                   (SELECT l.openalex_all FROM literature l WHERE l.motif = m.metal || '-' || m.donors) AS openalex_sac_hits
            FROM sac_models m;
    """)
    con.commit()
    counts = {n: con.execute("SELECT COUNT(*) FROM %s" % n).fetchone()[0] for n in tables}
    con.close()
    export(tables)
    write_schema(counts)
    print("wrote %s (%.0f MB)" % (DB, os.path.getsize(DB) / 1e6))
    for n, c in counts.items():
        print("  %-16s %9d rows" % (n, c))


def export(tables):
    """Every table as Parquet (typed) and gzipped CSV (universal) next to the database."""
    for fmt in ("parquet", "csv"):
        os.makedirs(os.path.join(OUT, fmt), exist_ok=True)
    for name, df in tables.items():
        df = df.copy()
        for col in df.columns[df.dtypes == object]:      # Parquet needs one type per column
            kinds = {type(v) for v in df[col].dropna()}
            if len(kinds) > 1:
                print("  %s.%s: mixed %s -> str" % (name, col, sorted(k.__name__ for k in kinds)))
                df[col] = df[col].map(lambda v: v if pd.isna(v) else str(v))
        df.to_parquet(os.path.join(OUT, "parquet", name + ".parquet"), index=False)
        df.to_csv(os.path.join(OUT, "csv", name + ".csv.gz"), index=False, compression="gzip")


def write_schema(counts):
    lines = ["# Database schema", "",
             "Generated by `db/build.py` from `COLUMNS`; the same text is in the table `column_docs`. Views: "
             "`donors`, `weak_contacts`, `metal_neighbors` (rows of `contacts` by kind), `v_sites` (sites with entry "
             "data and set names), `v_sac_overview` (per SAC model: matching sites, best dG_H, literature hits).", ""]
    for t, cols in COLUMNS.items():
        lines += ["## %s (%d rows)" % (t, counts.get(t, 0)), "", "| column | description |", "| --- | --- |"]
        lines += ["| %s | %s |" % (c, d) for c, d in cols.items()] + [""]
    with open(os.path.join(HERE, "SCHEMA.md"), "w") as f:
        f.write("\n".join(lines))


if __name__ == "__main__":
    main()
