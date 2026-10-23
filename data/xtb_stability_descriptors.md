# What holds a single metal atom in carbon: descriptors from GFN2-xTB

E_ins = E(metal in the site) - E(empty site); more negative holds the metal more strongly. 1143 sites, 12 metals, 96 donor arrangements. Drafts: xTB, checked against the DFT formation energies.

## 1. xTB against DFT

| site | E_f DFT (eV) | E_ins xTB (eV) | E_f from xTB (eV) |
|---|---|---|---|
| Co-N3C:DV | -2.67 | -104.72 | -3.19 |
| Co-N2OC:DV:NONC | -1.57 | -103.81 | -2.27 |
| Co-S2C:SV | +1.51 | -99.33 | +2.21 |
| Co-S2O:SV | +2.58 | -98.44 | +3.10 |
| Cr-S2O2:DV:trans | +3.17 | -54.00 | - |
| Cu-N2C:SV | +0.47 | -108.00 | +0.62 |
| Cu-N3:SV | +0.68 | -108.09 | +0.53 |
| Cu-N2S:SV | +1.44 | -106.65 | +1.97 |
| Cu-NOC:SV | +1.58 | -107.50 | +1.13 |
| Cu-NS2:SV | +1.63 | -106.67 | +1.95 |
| Cu-NSC:SV | +1.72 | -106.45 | +2.18 |
| Cu-NSO:SV | +1.98 | -106.67 | +1.96 |
| Cu-N2O:SV | +2.29 | -107.41 | +1.22 |
| Cu-SO2:SV | +2.49 | -105.92 | +2.71 |
| Fe-C2S2:DV:cis5 | +0.01 | -87.13 | +1.37 |
| Fe-N2S2:DV:cis5 | +0.50 | -87.00 | +1.50 |
| Fe-N2SO:DV:NSNO | +1.14 | -87.16 | +1.33 |
| Fe-SO2C:DV:COOS | +1.22 | -87.22 | +1.28 |
| Fe-N2O2:DV:trans | +1.53 | -88.02 | +0.48 |
| Fe-S2C:SV | +1.56 | -85.52 | +2.98 |
| Fe-NS2O:DV:SNSO | +1.63 | -86.37 | +2.13 |
| Fe-S3:SV | +2.27 | -84.23 | +4.27 |
| Fe-N2S:SV | +2.29 | -86.28 | +2.22 |
| Fe-NSO2:DV:ONSO | +2.37 | -86.66 | +1.84 |
| Fe-NS2:SV | +2.52 | -85.31 | +3.18 |
| Fe-S2O2:DV:cis5 | +3.05 | -85.63 | +2.86 |
| Fe-NO3:DV | +3.11 | -87.07 | +1.42 |
| Fe-O4:DV | +4.56 | -86.01 | +2.49 |
| Fe-SO2:SV | +5.02 | -85.07 | +3.43 |
| Mo-NS2O:DV:SNSO | +1.15 | -57.19 | +5.87 |
| Mo-NS3:DV | +2.38 | -65.81 | -2.75 |
| Mo-O4:DV | +7.42 | -55.24 | +7.82 |
| Ni-N2S2:DV:cis5 | +0.25 | -132.41 | +0.41 |
| Ni-N3:SV | +0.52 | -132.61 | +0.21 |
| Ni-S3:SV | +1.85 | -129.21 | +3.60 |
| Ni-N2S:SV | +1.86 | -131.19 | +1.62 |
| Ni-NS2:SV | +2.15 | -130.48 | +2.34 |
| Ni-N2O:SV | +2.42 | -132.06 | +0.76 |
| Ni-S2O:SV | +2.76 | -129.95 | +2.86 |
| Ru-S2O2:DV:trans | +4.82 | -85.95 | - |
| V-O4:DV | +3.96 | -53.39 | - |
| W-O4:DV | +9.02 | -66.01 | - |
| Zn-N3O:DV | -0.95 | -18.50 | -1.26 |
| Zn-N3:SV | -0.16 | -17.20 | +0.04 |
| Zn-NO2C:DV:CNOO | -0.07 | -17.71 | -0.47 |
| Zn-NS2:SV | +1.27 | -15.47 | +1.77 |

Within each metal (sites ordered the same way?): Co 4 sites, Spearman 1.00, spread of E_f - E_ins 0.71 eV; Cu 9 sites, Spearman 0.67, spread of E_f - E_ins 0.50 eV; Fe 15 sites, Spearman 0.60, spread of E_f - E_ins 1.22 eV; Mo 3 sites, Spearman 0.50, spread of E_f - E_ins 4.93 eV; Ni 7 sites, Spearman 0.50, spread of E_f - E_ins 1.00 eV; Zn 4 sites, Spearman 0.80, spread of E_f - E_ins 0.43 eV.
With one constant per metal (in-sample) the mean error is 0.86 eV over 42 sites.

## 2. Within one metal

E_ins relative to the metal's mean, 1143 sites, one model for all metals. Cross-validated R2 (5 folds over the arrangements): 0.60 with donors, vacancy type and strain; 0.60 with donors and vacancy type alone.

| descriptor | weight (eV per standard deviation) |
|---|---|
| n_N | +0.15 |
| n_S | +1.33 |
| n_O | +1.01 |
| DV | -0.70 |
| M_height | +0.10 |
| max_M_X | -0.16 |

Per metal: Co (96 sites): CV R2 0.91, strongest Co-N2C2:DV:cis5, Co-N3C:DV, Co-N4:DV; Cr (96 sites): CV R2 0.89, strongest Cr-N4:DV, Cr-N2C2:DV:trans, Cr-N3C:DV; Cu (94 sites): CV R2 0.86, strongest Cu-N2C2:DV:cis5, Cu-N2C2:DV:trans, Cu-N3C:DV; Fe (96 sites): CV R2 0.88, strongest Fe-N2C2:DV:cis5, Fe-N3C:DV, Fe-N4:DV; Mn (96 sites): CV R2 0.92, strongest Mn-N3C:DV, Mn-N2C2:DV:trans, Mn-N2C2:DV:cis5; Mo (91 sites): CV R2 0.34, strongest Mo-S4:DV, Mo-S3C:DV, Mo-NS3:DV; Ni (96 sites): CV R2 0.90, strongest Ni-N2C2:DV:cis5, Ni-N2C2:DV:trans, Ni-N3C:DV; Pt (96 sites): CV R2 0.85, strongest Pt-N2C2:DV:cis5, Pt-N2C2:DV:trans, Pt-N3C:DV; Ru (96 sites): CV R2 0.89, strongest Ru-N2C2:DV:cis5, Ru-N2C2:DV:trans, Ru-N3C:DV; V (96 sites): CV R2 0.88, strongest V-N2C2:DV:trans, V-N2C2:DV:cis5, V-N3C:DV; W (96 sites): CV R2 0.91, strongest W-N2C2:DV:cis5, W-N2C2:DV:trans, W-NC3:DV; Zn (94 sites): CV R2 0.73, strongest Zn-N4:DV, Zn-N3C:DV, Zn-N2C2:DV:cis5.

## 3. Across metals

Mean DFT E_f per metal against the bulk cohesive energy: Co -0.04 eV (E_coh 4.39), Cr 3.17 eV (E_coh 4.10), Cu 1.58 eV (E_coh 3.49), Fe 2.19 eV (E_coh 4.28), Mo 3.65 eV (E_coh 6.82), Ni 1.69 eV (E_coh 4.44), Ru 4.82 eV (E_coh 6.74), V 3.96 eV (E_coh 5.31), W 9.02 eV (E_coh 8.90), Zn 0.02 eV (E_coh 1.35); Spearman 0.77 over 10 metals.

## 4. The enzyme prior

Spearman of E_ins (within the metal) with log(1 + PDB entries of that metal with that donor set): 0.15 (p = 0.000, 1143 sites; 445 of them with the donor set seen in enzymes). Negative means donor sets that enzymes use with a metal hold it more strongly.

## 5. Escape energy (lower bound of the barrier out of the site)

An hour of pyrolysis gets an atom over 3.5 eV at 800 C and 4.2 eV at 1000 C (1e13 /s).

| metal | sites | median (eV) | max (eV) | above 3.5 eV |
|---|---|---|---|---|
| Co | 96 | 3.81 | 6.72 | 54 |
| Cr | 96 | 4.29 | 7.69 | 61 |
| Cu | 94 | 2.93 | 6.86 | 35 |
| Fe | 96 | 3.15 | 5.95 | 39 |
| Mn | 96 | 4.08 | 7.14 | 60 |
| Mo | 91 | -0.03 | 6.44 | 9 |
| Ni | 96 | 3.28 | 6.68 | 41 |
| Pt | 96 | 3.13 | 7.24 | 42 |
| Ru | 96 | 3.57 | 6.97 | 49 |
| V | 96 | 4.11 | 7.01 | 61 |
| W | 96 | 3.66 | 7.50 | 52 |
| Zn | 94 | 1.72 | 4.53 | 9 |

## 6. Does a forced geometry cost stability?

199 sites have both an xTB E_ins and enzyme sites to compare with (screen/geometry_mismatch.py). Spearman with E_ins relative to the metal's mean (positive: the more the support changes the enzyme's geometry, the more weakly it holds the metal):

- angles changed (RMS, deg): 0.36 (p = 0.000)
- metal shift out of the donor plane (A): 0.40 (p = 0.000)
- bond-length change (mean |model - enzyme|, A): -0.05 (p = 0.470)
- donors lost in the model: 0.00 (p = 0.951)

Cross-validated R2 on these 199 sites: 0.39 with donors and vacancy type, 0.39 with the mismatch added.

DFT check on the 34 sites with E_f and enzyme sites: Spearman of E_f with the angle change 0.21 (p = 0.24); Zn-N3O:DV E_f -0.95 eV, angles 20 deg, Zn-N3:SV E_f -0.16 eV, angles 10 deg, Ni-N2S2:DV:cis5 E_f +0.25 eV, angles 6 deg, Cu-N2C:SV E_f +0.47 eV, angles 20 deg, Fe-N2S2:DV:cis5 E_f +0.50 eV, angles 3 deg, Ni-N3:SV E_f +0.52 eV, angles 11 deg, Cu-N3:SV E_f +0.68 eV, angles 21 deg, Fe-N2SO:DV:NSNO E_f +1.14 eV, angles 4 deg, Zn-NS2:SV E_f +1.27 eV, angles 67 deg, 3 donors lost, Cu-N2S:SV E_f +1.44 eV, angles 13 deg, Fe-N2O2:DV:trans E_f +1.53 eV, angles 7 deg, Cu-NS2:SV E_f +1.63 eV, angles 18 deg, Fe-NS2O:DV:SNSO E_f +1.63 eV, angles 3 deg, Ni-S3:SV E_f +1.85 eV, angles 8 deg, Ni-N2S:SV E_f +1.86 eV, angles 6 deg, Cu-NSO:SV E_f +1.98 eV, angles 32 deg, Ni-NS2:SV E_f +2.15 eV, angles 30 deg, Fe-S3:SV E_f +2.27 eV, angles 12 deg, Cu-N2O:SV E_f +2.29 eV, angles 33 deg, Fe-N2S:SV E_f +2.29 eV, angles 12 deg, Fe-NSO2:DV:ONSO E_f +2.37 eV, angles 16 deg, Mo-NS3:DV E_f +2.38 eV, angles 15 deg, Ni-N2O:SV E_f +2.42 eV, angles 20 deg, Cu-SO2:SV E_f +2.49 eV, angles 57 deg, Fe-NS2:SV E_f +2.52 eV, angles 10 deg, Co-S2O:SV E_f +2.58 eV, angles 14 deg, Ni-S2O:SV E_f +2.76 eV, angles 23 deg, Fe-S2O2:DV:cis5 E_f +3.05 eV, angles 14 deg, Fe-NO3:DV E_f +3.11 eV, angles 7 deg, V-O4:DV E_f +3.96 eV, angles 13 deg, Fe-O4:DV E_f +4.56 eV, angles 15 deg, Fe-SO2:SV E_f +5.02 eV, angles 25 deg, Mo-O4:DV E_f +7.42 eV, angles 15 deg, W-O4:DV E_f +9.02 eV, angles 14 deg.

