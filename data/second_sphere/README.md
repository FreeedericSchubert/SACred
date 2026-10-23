# Second coordination sphere: start structures for the force-field optimisation

`second_sphere_start.extxyz`: the three MD sites (Zn–NS₂, Ni–NS₂, Co–N₂OC; DFT-relaxed, 6 × 6 graphene cell),
each as `bare` and `H_M` (H on the metal), with the polar groups that surround the metal's first shell in a
source enzyme, as small neutral molecules. Both states of a site carry the same molecules, so
ΔG_H(second sphere) = E(H_M) − E(bare) − E(H₂)/2 + 0.24 eV is a clean difference.

| site | source enzyme (PDB site) | second-sphere molecules |
|---|---|---|
| Zn–NS₂ (SV) | β-carbonic anhydrase (6LAC_A_ZN_201) | 2 water, 2 ammonia (backbone N), formaldehyde (backbone O) |
| Ni–NS₂ (SV) | lactate racemase, Ni pincer (8EZH_B_NI_514) | methanol (Tyr OH), ammonia (Trp NE1), imidazole (His), ammonia (backbone N) |
| Co–N₂OC (DV) | organophosphate hydrolase (3HTW_A_CO_2) | acetic acid (Asp), ammonia (backbone N), formaldehyde (backbone O) |

Info line per frame: `motif`, `state`, `source_site`, `enzyme`, `second_sphere` (residue, atom: molecule),
`metal_index`, `donor_indices`, `fixed_index` (the far carbon fixed in DFT). Per-atom `tags`: 0 = the SAC site,
k = the k-th second-sphere molecule (e.g. fix tag 0 at first, or restrain one heavy atom per molecule to mimic the
protein scaffold, so the molecules do not just drift off into the vacuum).

How they were built (screen/second_sphere.py): the enzyme is laid onto the model by its metal and first-shell donors
(best donor assignment); every residue or water with an N, O or S atom 3.0–6.5 Å from the metal that is not a donor
becomes a neutral molecule on that atom (water; Asp/Glu → acetic acid; His → imidazole; Ser/Thr/Tyr → methanol;
Cys/Met → methanethiol; backbone and Asn/Gln O → formaldehyde; Lys, Arg, backbone and Asn/Gln N → ammonia), turned
to keep clear of the sheet; groups that would sit inside the sheet or clash are left out (at most six).
Starting geometries, not relaxed.

After the force-field optimisation, the DFT single points:
`python dft/make_inputs.py points <optimised.extxyz> ss1` (keep `motif` and `state` in the info line).

Caveats: the force field was trained on graphene SAC frames only (no free molecules or hydrogen bonds), so the
DFT single point (energy and residual force) is the check. For Co–N₂OC the enzyme source is weak: in the PDB, the
C donor of N₂OC environments is mostly a bound small molecule (Tris, cacodylate, formate, carbonate), not protein.
