# How much the support changes the first shell

39343 enzyme sites (X-ray <= 2.0 A, mononuclear, nine enzyme metals) on 202 graphene models (52 DFT, 150 xTB). Model minus enzyme; medians with the 25-75 % range.

In 11 % of the enzyme sites the model lets go of at least one donor (longer than the donor cutoff + 0.3 A): DV 0 %, SV 22 %. The table below covers the sites whose model keeps all its donors.

| | sites | M-N (A) | M-S (A) | M-O (A) | angles, RMS (deg) | metal out of donor plane (A) | tau4 (DV) |
|---|---|---|---|---|---|---|---|
| DV, cofactor | 5698 | -0.14 (-0.16 to -0.12) | -0.01 (-0.19 to +0.03) | -0.09 (-0.22 to +0.04) | 3 (1-6) | -0.06 (-0.17 to -0.02) | -0.06 (-0.14 to -0.03) |
| DV, side chains | 14865 | -0.22 (-0.32 to -0.09) | -0.11 (-0.25 to +0.26) | -0.10 (-0.22 to +0.08) | 14 (7-22) | -0.02 (-0.17 to +0.66) | -0.07 (-0.21 to +0.46) |
| SV, cofactor | 128 | -0.31 (-0.37 to -0.17) | -0.19 (-0.27 to -0.10) | +0.22 (+0.05 to +0.29) | 16 (9-23) | +0.10 (-0.18 to +0.44) | - |
| SV, side chains | 14463 | -0.12 (-0.19 to +0.12) | +0.62 (+0.48 to +0.65) | +0.30 (+0.03 to +0.45) | 33 (13-47) | +0.81 (+0.26 to +1.29) | - |

Models that let go of a donor (by enzyme sites behind them): Zn-NS2:SV (DFT, 3078 sites, 3 lost); Zn-O3:SV (DFT, 1057 sites, 3 lost); Mo-S3O:DV (xTB, 35 sites, 1 lost); Mo-S2O2:DV:cis5 (xTB, 16 sites, 1 lost); Zn-NSC:SV (xTB, 3 sites, 1 lost).

Largest angular mismatch among the models that keep their donors (>= 10 enzyme sites):

- Zn-S2O:SV (xTB, 28 sites): 57 deg; metal 0.55 A out of the donor plane in enzymes, 2.32 A in the model
- Zn-S3:SV (xTB, 2544 sites): 56 deg; metal 0.74 A out of the donor plane in enzymes, 2.52 A in the model
- Zn-SO2:SV (xTB, 50 sites): 52 deg; metal 0.71 A out of the donor plane in enzymes, 2.24 A in the model
- Zn-NSO:SV (xTB, 266 sites): 48 deg; metal 0.73 A out of the donor plane in enzymes, 2.08 A in the model
- Zn-NO2:SV (xTB, 1823 sites): 44 deg; metal 0.72 A out of the donor plane in enzymes, 1.91 A in the model
- Zn-O4:DV (xTB, 958 sites): 41 deg; metal 0.21 A out of the donor plane in enzymes, 1.58 A in the model
- W-O3:SV (xTB, 40 sites): 40 deg; metal 0.61 A out of the donor plane in enzymes, 1.74 A in the model
- Zn-N2S:SV (xTB, 808 sites): 40 deg; metal 0.74 A out of the donor plane in enzymes, 1.74 A in the model
