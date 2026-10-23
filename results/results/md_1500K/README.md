# Local 1500 K MD analysis — provenance and gaps

Bundled actual extxyz coordinate frames: {
  "Zn_bonded_0_to_3p15ps.extxyz": {
    "frames": 64,
    "start_ps": 0.0,
    "end_ps": 3.15,
    "skipped_lines": 0
  },
  "Ni_3p1_to_20ps_PARTIAL.extxyz": {
    "frames": 339,
    "start_ps": 3.1,
    "end_ps": 20.0,
    "skipped_lines": 4600
  }
}

- `Zn_bonded_0_to_3p15ps.extxyz`: correct bonded-start Zn data.
- `Ni_3p1_to_20ps_PARTIAL.extxyz`: partial Ni only; do not imply early frames are included.
- **Co raw trajectory absent locally**: Co is represented in static visualization and log-observed 4/4 contacts only; request `pyrolysis_1500K/Co-N2OC_DV_NONC_1500K_r0.extxyz` from Modal.
- `metal_donor_distances.csv`, `analysis_summary.json`: computed from saved Zn and Ni atomic coordinates, with 3D periodic minimum-image distances. Original concatenated export was incomplete and included invalid Zn run; avoid raw concatenation.
- Images contain log-only observations for Co, not invented XYZ coordinates.

See `../../MODAL_DOWNLOAD_CHECKLIST.md` at the ZIP top level for download commands and other missing outputs.
