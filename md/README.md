# 1500 K single-atom catalyst MD — MACE v2

From the repository root (after Modal is configured and model exists in the Modal Volume):

```bash
modal run md/run_1500K_modal.py --mode smoke
modal run md/run_1500K_modal.py --mode production
modal run md/run_zn_bonded_modal.py --target-ps 5
```

Launchers in this bundle resolve `data/md_start/` and `data/second_sphere/` relative to the repo root. The original inputs are **already included in the team's repository** and are not duplicated in this upload overlay.

IMPORTANT: The original DFT-relaxed Zn starting configuration was already outside the donor-bond cutoffs; its original apparent escape at ~0.01 ps is invalid. The bonded-start Zn experiment uses first xTB NEB image instead. Its first loss of at least two contacts is at ~0.15 ps, all contacts lost by ~0.25 ps; a continuous 1 ps loss condition is confirmed at ~1.15 ps.

Ni: partial coordination loss around 3.5 ps; retained two metal–donor contacts through 20 ps. The bundled extxyz for Ni only contains complete frames from 3.1 ps onwards. Co: 4/4 donor contacts in available MD logs through 5 ps; its full extxyz is NOT available locally.

Results are exploratory: single trajectory per site, different run lengths and initial geometries, and MACE not validated for 1500 K transition-state configurations. Never present a single event time as a measured mean lifetime.

Script note: `scripts/analyze_pyrolysis_original.py`, `scripts/make_figures_original.py`, and `scripts/make_movie_*.py` are copies of local analysis/visualization scripts and contain some absolute `/mnt/data` paths. Update those paths for portable reruns. MD launchers are updated to resolve repository-relative directories.
