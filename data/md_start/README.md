# Escape MD at 1000 K: start structures, script and DFT reference

`md_start_sites.extxyz`: three DFT-relaxed bare sites in the 6 × 6 graphene cell of the training data
(PBE-D3, Quantum ESPRESSO; c = 20 Å), chosen for clearly different stability. Per frame in the info line:
`motif`, `E_f_dft_eV` (metal in the site vs bulk metal; < 0 = no driving force to cluster), `escape_energy_xtb_eV`
(site → adatom on clean graphene), `neb_barrier_xtb_eV`, expected residence times `tau_<T>K_s` from
τ = 1 / (10¹³ s⁻¹ · exp(−E/kT)) (NEB barrier where there is one, else the escape energy), `metal_index`,
`donor_indices`, `fixed_index` (the far carbon that was fixed in DFT; fix it or remove the centre-of-mass drift).

| site | DFT E_f (eV) | barrier (eV) | τ at 1000 K | τ at 1500 K | τ at 2000 K | expectation |
|---|---|---|---|---|---|---|
| Zn–NS₂ (single vacancy) | +1.27 | 0.38 (NEB, xTB); 0.12 in DFT along that path | 8 ps | 2 ps | 1 ps | leaves almost at once |
| Ni–NS₂ (single vacancy) | +2.15 | 0.96 (NEB, xTB); up to 1.96 in DFT along that path | 7 ns (xTB) | 0.16 ns | 26 ps | xTB: leaves in long runs; DFT: should stay at 1000 K |
| Co–N₂OC (double vacancy) | −1.57 | 6.0 (NEB, xTB) | ~10⁸ years | ~3 weeks | ~30 s | stays |

Escape criterion used in `screen/md_escape.py`: fewer than two donors within their bond cutoff + 0.5 Å
(N, O 2.6 Å; S 2.8 Å; C 2.3 Å) for 1 ps. Drafts: barriers from GFN2-xTB.
Temperatures above ~3000 K leave the force field's training range (relaxation frames near minima) and approach
the melting of graphene itself; use 1000 K, and 1500–2000 K to speed up escapes for an Arrhenius check.

## Files

- `md_start_sites.extxyz`, `md_start_sites.csv`: the three sites (above).
- `md_escape_mace.py`: standalone MD with ASE + MACE (Langevin, escape detection, a frame every 50 fs, summary.csv):
  `python md_escape_mace.py --model <mace.model> --device cuda --T 1000 --ps 100 --replicas 3`.
- `escape_paths/`: the GFN2-xTB NEB paths of Zn–NS₂ and Ni–NS₂ (`xtb_neb_path_*.extxyz`, 9 images each) and DFT
  single points on them (`dft_single_points_on_xtb_paths.extxyz`: energy and forces, PBE-D3, Γ point, 17 of 18
  images; `escape_path_energies.csv`). These are the configurations the metal passes on its way out and are not in
  the training data: check the force field on them first (`python md_escape_mace.py --check --model ...`); if it
  reproduces the DFT profile (Zn: flat, ≤ 0.12 eV; Ni: rising to ~1.5–2 eV), the escapes in the MD can be trusted.
  The DFT forces on these xTB geometries are large (1–2.7 eV/Å), so the DFT numbers are energies along the xTB path,
  an upper bound of the DFT barrier, not a DFT NEB.
