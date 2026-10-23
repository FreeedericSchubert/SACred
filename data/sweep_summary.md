# PDB sweep: transition-metal first shells for graphene SAC motifs

Generated 2026-10-20 by `bio/pdb_sweep.py` (search, download) and `bio/first_shell.py` (extraction, motifs,
checks) from the RCSB PDB state of that day. All outputs are drafts that need expert review.

## Searches
* Enzyme search (RCSB Search API v2): EC lineage per class, resolution <= 2.5 A, polymer entities grouped at
  95 % sequence identity, best resolution as representative: **384 representative
  entities in 281 entries**.
* EC classes without any hit at <= 2.5 A: 1.10.3.11 (0 entities at any resolution), 1.17.1.10 (1 entities at any resolution), 1.17.1.11 (0 entities at any resolution), 1.17.2.3 (0 entities at any resolution), 1.19.6.1 (0 entities at any resolution).
* Whole-PDB sweep: entries with a bound component containing Fe/Ni/Co/Cu/Mo/W/Mn/V, resolution <= 2.0 A,
  protein entities grouped at 30 % identity: **2572 unique representative entries**
  (below the 6000 cap, all
  selected and downloaded). The chemical-service formula search (match_subset) returned only 130 Fe
  components and missed SF4/FES/F3S/ICS/CLF, so the metal component ids were taken from the wwPDB Chemical
  Component Dictionary; rcsb_entry_info.nonpolymer_bound_components is not searchable, the filter uses
  rcsb_nonpolymer_entity_container_identifiers.nonpolymer_comp_id.

| reaction | EC | enzyme | entities <=2.5 A | 95% representatives | entries |
|---|---|---|---|---|---|
| HER | 1.12 | hydrogenase | 412 | 53 | 29 |
| ORR | 1.10.3.2 | laccase | 189 | 38 | 38 |
| ORR | 1.10.3.3 | ascorbate oxidase | 3 | 1 | 1 |
| ORR | 1.3.3.5 | bilirubin oxidase | 35 | 2 | 2 |
| ORR | 1.16.3.1 | ferroxidase | 412 | 59 | 58 |
| ORR | 7.1.1.9 | cytochrome c oxidase | 335 | 34 | 18 |
| ORR | 1.9.3.1 | cytochrome c oxidase | 325 | 28 | 12 |
| ORR | 7.1.1.7 | cytochrome bd oxidase | 15 | 3 | 1 |
| ORR | 7.1.1.3 | bo3 quinol oxidase | 4 | 3 | 2 |
| ORR | 1.10.3.11 | alternative oxidase | 0 | 0 | 0 |
| CO2RR | 1.2.7.4 | Ni-Fe CO dehydrogenase | 99 | 9 | 9 |
| CO2RR | 1.2.5.3 | Mo-Cu CO dehydrogenase | 24 | 6 | 2 |
| CO2RR | 1.17.1.9 | formate dehydrogenase | 45 | 14 | 14 |
| CO2RR | 1.17.1.10 | formate dehydrogenase | 0 | 0 | 0 |
| CO2RR | 1.17.1.11 | formate dehydrogenase | 0 | 0 | 0 |
| CO2RR | 1.17.5.3 | formate dehydrogenase | 1 | 1 | 1 |
| CO2RR | 1.17.2.3 | formate dehydrogenase | 0 | 0 | 0 |
| CO2RR | 1.17.98.3 | formate dehydrogenase | 3 | 1 | 1 |
| CO2RR | 1.17.99.7 | formate dehydrogenase | 1 | 1 | 1 |
| CO2RR | 1.2.1.2 | formate dehydrogenase | 67 | 18 | 13 |
| CO2RR | 2.3.1.169 | acetyl-CoA synthase | 33 | 3 | 3 |
| CO2RR | 1.2.7.12 | formylmethanofuran dehydrogenase | 5 | 4 | 2 |
| CO2RR | 2.8.4.1 | methyl-coenzyme M reductase | 99 | 34 | 13 |
| CO2RR | 4.2.1.1 | carbonic anhydrase | 1566 | 60 | 60 |
| NRR | 1.18.6.1 | nitrogenase | 268 | 35 | 21 |
| NRR | 1.19.6.1 | nitrogenase | 0 | 0 | 0 |

| metal | CCD components | protein entities <=2.0 A | 30% clusters | representative entries | sites (metal) in sweep |
|---|---|---|---|---|---|
| Fe | 217 | 9236 | 1114 | 868 | 2866 |
| Ni | 41 | 2021 | 515 | 464 | 781 |
| Co | 75 | 1020 | 337 | 313 | 591 |
| Cu | 53 | 2398 | 341 | 267 | 644 |
| Mo | 79 | 310 | 58 | 42 | 154 |
| W | 44 | 150 | 67 | 57 | 352 |
| Mn | 26 | 3139 | 664 | 600 | 1284 |
| V | 66 | 149 | 63 | 61 | 145 |

## Sites
Enzyme set: 1811 unique transition-metal sites in 229 entries (253 DV, 1475 SV,
83 none).
Whole-PDB set: 7142 unique sites in 2567 entries.

| reaction | hit entries | entries with TM sites | unique TM sites |
|---|---|---|---|
| CO2RR | 109 | 84 | 562 |
| HER | 29 | 27 | 497 |
| NRR | 21 | 21 | 261 |
| ORR | 122 | 97 | 491 |

## Top 25 graphene motifs, enzyme set (sorted by n_entries; the 83 "none" sites, CN <= 2, not listed)
| motif | n_sites | n_entries | n_sequence_clusters | reactions | frac_mononuclear | example_pdb |
|---|---|---|---|---|---|---|
| Fe-S3(SV) | 1017 | 73 | 44 | CO2RR;HER;NRR;ORR | 0.02 | 4U9H(0.89) 5JSK(0.95) 9GZ4(0.96) 3U7Q(1.00) 7ADR(1.00) 9FPH(1.01) |
| Cu-N3(SV) | 88 | 48 | 13 | ORR | 1.0 | 6RI6(0.93) 3PXL(1.20) 7PUH(1.30) 2Q9O(1.30) 5Z1X(1.38) 6IQX(1.43) |
| Cu-N2S(SV) | 40 | 39 | 11 | ORR | 1.0 | 6RI6(0.93) 6EVG(1.10) 3PXL(1.20) 2Q9O(1.30) 7PUH(1.30) 5Z1X(1.38) |
| Fe-NS2(SV) | 40 | 34 | 15 | CO2RR;HER;NRR;ORR | 0.2 | 4U9H(0.89) 5JSK(0.95) 3U7Q(1.00) 9FPH(1.01) 6N59(1.02) 6FPO(1.05) |
| Zn-NS2(SV) | 31 | 25 | 15 | CO2RR;HER;ORR | 1.0 | 6LAC(1.30) 2W3Q(1.34) 3LAS(1.40) 9EAT(1.43) 5SWC(1.45) 3QY1(1.54) |
| Zn-N3(SV) | 20 | 19 | 6 | CO2RR | 1.0 | 3K34(0.90) 8Q3U(1.10) 9F2O(1.12) 4CNX(1.23) 7Q0D(1.24) 5HPJ(1.50) |
| Fe-N4(DV)+ax:N2 | 23 | 17 | 7 | CO2RR;ORR | 1.0 | 7COH(1.30) 7W3E(1.45) 5B1A(1.50) 1KQF(1.60) 1V54(1.80) 3S8F(1.80) |
| Cu-N2O(SV) | 19 | 16 | 7 | ORR | 1.0 | 6RI6(0.93) 7PUH(1.30) 5Z1X(1.38) 3T9W(1.50) 2QT6(1.50) 7PEN(1.60) |
| Cu-NS2(SV) | 34 | 16 | 5 | ORR | 0.03 | 7COH(1.30) 4W9Z(1.30) 7W3E(1.45) 5B1A(1.50) 2CUA(1.60) 1V54(1.80) |
| Ni-S3(SV) | 13 | 13 | 6 | CO2RR;HER | 0.0 | 4U9H(0.89) 5JSK(0.95) 9R52(1.65) 4OMF(1.71) 1E3D(1.80) 6QGR(1.84) |
| Fe-S2C2(DV,cis)+ax:C | 13 | 13 | 5 | HER | 0.0 | 4U9H(0.89) 5JSK(0.95) 6N59(1.02) 3AYX(1.18) 6SZK(1.20) 6SYO(1.25) |
| Fe-O3(SV) | 22 | 12 | 11 | HER;ORR | 0.91 | 6B8F(1.06) 5OK4(1.29) 3BVF(1.50) 4IWK(1.65) 1E3D(1.80) 1Z6O(1.91) |
| Zn-S3(SV) | 12 | 11 | 6 | CO2RR;HER;ORR | 0.75 | 7COH(1.30) 7W3E(1.45) 5B1A(1.50) 4OMF(1.71) 1V54(1.80) 1OAO(1.90) |
| Fe-S2C(SV) | 83 | 11 | 4 | CO2RR;NRR | 0.0 | 3U7Q(1.00) 7ADR(1.00) 4WES(1.08) 9SQ3(1.21) 8BOQ(1.55) 9HBC(1.82) |
| Fe-O4(DV)+ax:N | 12 | 10 | 5 | ORR | 1.0 | 5J93(1.10) 8PP3(1.55) 4IWK(1.65) 6L55(1.78) 7USN(1.79) 6LJG(1.80) |
| Fe-N4(DV)+ax:N | 16 | 9 | 6 | ORR | 1.0 | 1V54(1.80) 7OY2(2.06) 8UGH(2.10) 9RZV(2.23) 9WCY(2.26) 1M56(2.30) |
| Fe-N4(DV)+ax:NS | 10 | 8 | 6 | ORR | 1.0 | 1W2L(1.30) 7OY2(2.06) 8UGH(2.10) 9RZV(2.23) 9WCY(2.26) 8OVD(2.30) |
| Fe-N4(DV)+ax:NO | 8 | 8 | 3 | ORR | 1.0 | 7COH(1.30) 5B1A(1.50) 3S8F(1.80) 2GSM(2.00) 7N9Z(2.19) 3HB3(2.25) |
| Ni-N4(DV)+ax:SO | 9 | 8 | 1 | CO2RR | 1.0 | 5A0Y(1.10) 9QR3(1.34) 1E6Y(1.60) 7NKG(1.60) 5A8W(1.80) 5N1Q(1.90) |
| Fe-NO3(DV)+ax:O2 | 11 | 7 | 7 | HER;ORR | 1.0 | 5JSK(0.95) 8W1D(1.30) 4AM5(1.58) 6QGR(1.84) 5N5E(2.03) 5DA5(2.06) |
| Fe-NO2(SV) | 14 | 7 | 6 | ORR | 0.86 | 8PP3(1.55) 4AM5(1.58) 1NFV(1.95) 6LP5(1.98) 7Y74(1.98) 2FKZ(2.00) |
| Mo-NS2O(DV,cis)+ax:SO | 9 | 7 | 2 | NRR | 0.0 | 3U7Q(1.00) 4WES(1.08) 9SQ3(1.21) 1QGU(1.60) 9HBC(1.82) 9HBN(1.82) |
| Fe-N4(DV)+ax:S2 | 7 | 7 | 2 | ORR | 1.0 | 8SQO(1.55) 4AM5(1.58) 9BTS(1.85) 26RV(1.89) 3UOI(1.90) 1NFV(1.95) |
| Zn-N2O2(DV,cis)+ax:N | 7 | 6 | 6 | CO2RR | 1.0 | 5OGJ(1.06) 6RQQ(1.28) 4XFW(1.52) 4COQ(1.55) 1Y7W(1.86) 1KEQ(1.88) |
| Ni-S3O(DV)+ax:S | 7 | 6 | 2 | HER | 0.14 | 4UQL(1.22) 6SYO(1.25) 7ODH(1.34) 8DQV(1.52) 5XVD(1.57) 3MYR(2.10) |

## Top 15 graphene motifs, whole-PDB sweep (662 "none" sites not listed)
| motif | n_sites | n_entries | n_sequence_clusters | frac_mononuclear | example_pdb |
|---|---|---|---|---|---|
| Fe-S3(SV) | 1552 | 251 | 245 | 0.03 | 30OH(0.43) 30OR(0.43) 5D8V(0.48) 1YK4(0.69) 4U9H(0.89) 6LK1(0.90) |
| Mn-O4(DV)+ax:O2 | 361 | 219 | 216 | 0.96 | 2B97(0.75) 2JFR(0.83) 4PSS(0.85) 3M5Q(0.93) 2AU7(1.05) 3W42(1.06) |
| Fe-N4(DV)+ax:N2 | 239 | 117 | 106 | 1.0 | 3X34(0.76) 1X8P(0.85) 1J0P(0.91) 3QR7(0.94) 1M1Q(0.97) 9GTJ(0.99) |
| Mn-NO3(DV)+ax:O2 | 116 | 93 | 89 | 0.98 | 1MUW(0.86) 1NLS(0.94) 4ZA9(1.01) 4QRN(1.07) 31RQ(1.08) 7O91(1.10) |
| Cu-N3(SV) | 127 | 89 | 61 | 1.0 | 6RI6(0.93) 9KWT(0.96) 4A7U(0.98) 3F7L(0.99) 6Z5Y(1.01) 3ZIY(1.01) |
| Mn-O4(DV)+ax:O | 106 | 87 | 87 | 0.95 | 9MBI(1.09) 9T4C(1.21) 8ZXX(1.23) 3RL4(1.29) 6GIU(1.39) 4M5R(1.40) |
| Mn-O3(SV) | 105 | 86 | 86 | 0.9 | 1MUW(0.86) 3X0O(1.09) 9MBI(1.09) 9T4C(1.21) 3RL4(1.29) 8ONX(1.30) |
| Fe-N4(DV)+ax:NS | 98 | 85 | 78 | 1.0 | 4EIC(0.84) 1C75(0.97) 6K9J(0.98) 3ZIY(1.01) 4EIE(1.03) 4N0K(1.05) |
| Fe-N4(DV)+ax:NO | 65 | 60 | 56 | 1.0 | 5YCE(0.77) 4O6U(0.89) 3M5Q(0.93) 6RK0(0.99) 4W7K(1.05) 7CKA(1.06) |
| Fe-NO3(DV)+ax:NO | 68 | 59 | 58 | 1.0 | 8AXY(1.05) 1DS1(1.08) 9P9O(1.12) 3TM8(1.28) 2RDQ(1.31) 4Z6M(1.35) |
| Mn-NO3(DV)+ax:NO | 64 | 59 | 55 | 0.98 | 7B1S(0.99) 3P4H(1.10) 6F4Q(1.12) 6EY1(1.20) 7QOE(1.20) 3GVE(1.25) |
| Cu-N2S(SV) | 67 | 58 | 50 | 1.0 | 2OV0(0.75) 8K9N(0.86) 5AKR(0.87) 2XJH(0.92) 6RI6(0.93) 9KWT(0.96) |
| Cu-N2O(SV) | 71 | 57 | 49 | 1.0 | 6RI6(0.93) 5O2X(0.95) 8Q9X(1.05) 5MSZ(1.10) 6REK(1.10) 4E9X(1.14) |
| Ni-NO3(DV)+ax:NO | 59 | 54 | 52 | 1.0 | 5NFK(0.98) 6W1G(1.14) 9FFJ(1.27) 6GG1(1.30) 2GW3(1.40) 2XQO(1.40) |
| Fe-N2O2(DV,cis)+ax:O2 | 59 | 50 | 49 | 1.0 | 4XPX(1.03) 5LUN(1.08) 9MKJ(1.10) 5UQ6(1.18) 7ETK(1.22) 3TM8(1.28) |

## Sanity checks (19/20 passed)
| check | entry | expected | passed | observed | note |
|---|---|---|---|---|---|
| [NiFe]-H2ase DvMF Ni (Ni-R) | 4U9H | Ni: S4 = 4 Cys (2 bridging) | True | Ni EJ2 601/NI: S4 [S_Cys S_Cys S_Cys S_Cys] -> Ni-S3(SV) |  |
| [NiFe]-H2ase DvMF Fe (Ni-R) | 4U9H | Fe: C3S2 = 2 CN + CO + 2 bridging Cys | True | Fe EJ2 601/FE: S2C3 [C_CO/CN C_CO/CN C_CO/CN S_Cys S_Cys] -> Fe-S2C2(DV,cis)+ax:C |  |
| [NiFe]-H2ase DvMF Ni (Ni-A, alt) | 1WUI | Ni: S4 = 4 Cys (2 bridging) | False | Ni NFC 1004/NI: S4O4 [O_CSO O_CSO O_NFC O_NFC S_CSO S_CSO S_Cys S_Cys] -> Ni-S3O(DV)+ax:SO3 | Ni-A (oxidised, unready) model: Cys84/Cys546 partly oxidised to sulfenates (CSO) and bridging O species superposed at partial occupancy; the structure differs, S4 holds for the active Ni-R state (4U9H) |
| [FeFe]-H2ase CpI H-cluster Fe | 3C8Y | Fe_p/Fe_d: CO/CN carbon + S donors | True | Fe HCN 575/FE1: S3C3 [C_CO/CN C_CO/CN C_CO/CN S_Cys S_HCN S_HCN] -> Fe-S2C2(DV,cis)+ax:SC; Fe HCN 575/FE2: S2OC3 [C_CO/CN C_CO/CN C_CO/CN O_HCN S_HCN S_HCN] ->  |  |
| Nitrogenase MoFe Mo | 3U7Q | Mo: S3 N O2 (His, homocitrate) | True | Mo ICS 6496/MO1: NS3O2 [N_His O_HCA O_HCA S_inorganic S_inorganic S_inorganic] -> Mo-NS2O(DV,cis)+ax:SO |  |
| Nitrogenase FeMo-co belt Fe | 3U7Q | 6 Fe: S3 C (inorganic S, carbide) | True | Fe ICS 6496/FE2: S3C [C_carbide S_inorganic S_inorganic S_inorganic] -> Fe-S2C(SV); Fe ICS 7496/FE3: S3C [C_carbide S_inorganic S_inorganic S_inorganic] -> Fe-S |  |
| ACS distal Ni | 1RU3 | Ni: N2S2 (2 backbone N, 2 Cys), cis | True | Ni NI 735/NI: N2S2 [N_backbone N_backbone S_Cys S_Cys] -> Ni-N2S2(DV,cis) (weak: C_GLY:2.74) |  |
| ACS distal Ni (alt) | 1OAO | Ni: N2S2 (2 backbone N, 2 Cys), cis | True | Ni NI 1732/NI: N2S2 [N_backbone N_backbone S_Cys S_Cys] -> Ni-N2S2(DV,cis) (weak: C_GLY:2.78) |  |
| CcO heme a3 Fe | 1V54 | Fe: N4 porphyrin + N(His) axial | True | Fe HEA 516/FE: N5 [N_His N_porphyrin N_porphyrin N_porphyrin N_porphyrin] -> Fe-N4(DV)+ax:N |  |
| CcO CuB | 1V54 | CuB: N3 (3 His) | True | Cu CU 517/CU: N3 [N_His N_His N_His] -> Cu-N3(SV) |  |
| CcO heme a3 Fe (alt) | 5B1A | Fe: N4 porphyrin + N(His) axial | True | Fe HEA 602/FE: N5O [N_His N_porphyrin N_porphyrin N_porphyrin N_porphyrin O_PER] -> Fe-N4(DV)+ax:NO |  |
| CcO CuB (alt) | 5B1A | CuB: N3 (3 His) | True | Cu CU 603/CU: N3O [N_His N_His N_His O_PER] -> Cu-N3O(DV) |  |
| Laccase T1 Cu | 6RI6 | T1 Cu: N2S (2 His, Cys) (+S_Met) | True | Cu CU 504/CU: N2S [N_His N_His S_Cys] -> Cu-N2S(SV) |  |
| Laccase T3 Cu | 6RI6 | T3 Cu: N3 (3 His) | True | Cu CU 501/CU: N3O [N_His N_His N_His O_water] -> Cu-N3(SV) (weak: O_OXY:2.41(occ0.20) O_OXY:2.54(occ0.20) O_water:2.92); Cu CU 502/CU: N3 [N_His N_His N_His] -> |  |
| Ascorbate oxidase T1 Cu (Met) | 1AOZ | T1 Cu: N2S + S_Met | True | Cu CU 701/CU: N2S [N_His N_His S_Cys] -> Cu-N2S(SV) (weak: S_Met:2.90) |  |
| MCR Ni F430 | 1MRO | Ni: N4 (F430) + O axial | True | Ni F43 700/NI: N4SO [N_porphyrin N_porphyrin N_porphyrin N_porphyrin O_GLN S_COM] -> Ni-N4(DV)+ax:SO |  |
| MCR Ni F430 (alt) | 1HBN | Ni: N4 (F430) + O axial | True | Ni F43 1550/NI: N4SO [N_porphyrin N_porphyrin N_porphyrin N_porphyrin O_GLN S_COM] -> Ni-N4(DV)+ax:SO |  |
| Carbonic anhydrase II Zn | 3KS3 | Zn: N3 (His) + O_water | True | Zn ZN 262/ZN: N3O [N_His N_His N_His O_water] -> Zn-N3(SV) (weak: O_water:2.57(occ0.29)) |  |
| FDH-H Mo (oxidized) | 1FDO | Mo: S4 dithiolene + Se (SeCys) (+S/O) | True | Mo 6MO 803/MO: S4OSe [O_water S_dithiolene S_dithiolene S_dithiolene S_dithiolene Se_Sec] -> Mo-S3O(DV)+ax:SSe |  |
| FDH-H Mo (reduced, alt) | 1AA6 | Mo: S4 dithiolene + Se (SeCys) (+S/O) | True | Mo 4MO 803/MO: S4Se [S_dithiolene S_dithiolene S_dithiolene S_dithiolene Se_Sec] -> Mo-S4(DV)+ax:Se |  |

## Assumptions and thresholds
* Donor cutoffs (A): N 2.6, O 2.6, S 2.8, Se 2.9, C 2.3, Cl 2.8; metal-metal <= 3.0 A recorded
  as metal_neighbors; F/Br/I and other elements are not counted as donors. Contacts up to cutoff + 0.5 A
  are listed as weak_contacts (e.g. Cu-S(Met) 2.9-3.3 A of T1 Cu), not as donors.
* Hydrogens removed; per residue the altloc with the highest mean occupancy is kept (ties: first letter); ligands
  modelled as separate residues with different altlocs are resolved the same way where they overlap (< 1 A).
* Metals closer than 2.0 A to another metal are alternative positions: only the higher-occupancy one is
  kept. Contacts with occupancy < 0.5 x metal occupancy are minority states (weak_contacts).
* Protein carbon atoms are never donors; of the two His ring N only the closer one can bind; two donors less
  than 25.0 deg apart (L-M-L) are copies of one ligand (the original is kept); for a half-occupied
  metal off special positions, half-occupied atoms of symmetry mates are the other orientation of a ligand
  disordered across a symmetry element (both go to weak_contacts as "alt").
* Symmetry: neighbour search with all crystal symmetry operators (sym_donors counts donors from images);
  EM entries without a crystal cell are searched without symmetry.
* Donor classes from residue/atom names; S_inorganic = S without non-metal neighbour, S_dithiolene =
  S-C-C-S with two thiolate S, C_CO/CN = C bonded to one terminal N/O, C_carbide = C without non-metal
  neighbour in a multi-metal cluster, N_porphyrin = >= 4 ring N (bonded to two C) of one ligand around the
  metal (or of a known tetrapyrrole); everything else <element>_<residue>.
* Copies of a site in several chains (same metal atom, same donor residues/atoms; without protein donors also
  the same nearest protein residue) are counted once (n_copies).
* DV: most square-planar 4-subset, score = RMS(cis - 90) + RMS(trans - 180) deg; a tetrapyrrole always defines
  the plane; subsets within 5.0 deg of the best leave the more labile donors axial. DV is preferred when
  CN >= 4 and score <= 55.0 deg, set at the minimum between the two modes (square planar/octahedral
  ~12 deg, tetrahedral ~87 deg) of the score distribution: the KDE minimum of the 7231 CN >= 4 sites of both
  sets lies at 54.5 deg, and the split agrees with tau4(best 4-subset) <= 0.5 for 98.5% of them.
* SV (CN >= 3): 3-subset forming a face (all angles <= 140.0 deg) whose complement holds the most
  labile donors (water > other O > longest bond); score = RMS(angles - 90) deg. CN <= 2: none.
* Isomer for 2+2 and 2+1+1 DV compositions: trans if the like donors are > 135.0 deg apart.
* Reaction/enzyme labels: EC lineage of the polymer entity that provides most protein donors (label_source =
  entity), else the union over the entry (label_source = entry). n_sequence_clusters uses RCSB 30 % clusters.

## Limitations
* Representatives only: 95 % (enzymes) and 30 % (sweep) clusters, chosen by best resolution; the sweep
  representative of a cluster need not be the chain that binds the metal (entry-level component filter).
* Formate dehydrogenase classes 1.17.1.9/1.2.1.2 contain metal-free NAD-dependent enzymes; ferroxidase 1.16.3.1
  is dominated by ferritins; these entries add few or no sites.
* Zn and Cr sites in the sweep are incidental (Zn/Cr were not sweep metals) and therefore incomplete.
* No functional filter: the sweep contains adventitious ions from crystallisation (e.g. Mn/Co/Ni with water and
  carboxylate shells), polyoxometalates and metallocene ligands (each Cp ring counts as 5 C donors).
* Very recent entries may not be in the RCSB cluster files yet and then form their own cluster (e.g. 30OH/30OR,
  two refinements of one rubredoxin data set).
* Fe-S cluster iron dominates site counts (frac_mononuclear shows cluster sites); geometric DV/SV mapping does not
  consider oxidation state, spin or electronic structure.
