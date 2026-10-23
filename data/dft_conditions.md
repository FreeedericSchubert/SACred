# Co-N2OC under H2-evolution conditions, and DFT escape energies

| site | quantity | state | what | value (eV) |
|---|---|---|---|---|
| Co-N2OC:DV:NONC | dG_H | Hdonor_C | H on the C donor | -0.02 |
| Co-N2OC:DV:NONC | dG_H | Hdonor_O | H on the O donor | +0.03 |
| Co-N2OC:DV:NONC | dG_H | CH_HM | H on Co, C donor hydrogenated | +0.61 |
| Co-N2OC:DV:NONC | dG_H | H_M | H on Co (reference) | +0.03 |
| Co-N2OC:DV:NONC | dG_cap | H2 | H capping of the emptied site (2 H) | -2.82 |
| Co-N2OC:DV:NONC | dG_dem | H2 | acid demetalation, pH 0 | -1.81 |
| Co-N2OC:DV:NONC | dG_cap | H4 | H capping of the emptied site (4 H) | -1.34 |
| Co-N2OC:DV:NONC | dG_dem | H4 | acid demetalation, pH 0, +0.0 V vs RHE | -0.32 |
| Co-N2OC:DV:NONC | dG_dem | H4 | acid demetalation, pH 0, -0.2 V vs RHE | -0.72 |
| Co-N4:DV | dG_cap | H2 | H capping of the emptied site (2 H) | running |
| Co-N4:DV | dG_dem | H2 | acid demetalation, pH 0 | running |
| Zn-NS2:SV | E_escape | ESC_end | escape energy, DFT (xTB: +0.34 eV, NEB barrier 0.38 eV) | running |
| Ni-NS2:SV | E_escape | ESC_end | escape energy, DFT (xTB: +0.72 eV, NEB barrier 0.96 eV) | running |
| Cu-N3:SV | E_escape | ESC_end | escape energy, DFT (xTB: +3.01 eV, NEB barrier 3.16 eV) | running |

dG_dem > 0: the metal stays in the site; dG_H on a donor below about 0 eV: the donor is hydrogenated under H2 evolution. E_escape: metal from its site to an adatom on clean graphene about 5 A away (metal kept on its hollow site laterally), the lower bound of the escape barrier.
