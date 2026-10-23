# MACE v2 interatomic potential — our hackathon contribution

Files:
- `prepare_data.py`: processes a DFT-labelled `dft.zip` into `train.xyz`, `valid.xyz`, `test.xyz` and `manifest.csv` using motif-disjoint splits.
- `audit_v2.py`: audits prepared splits and their atomic compositions.
- `train_v2_modal_fixed.py`: fine-tuning of MACE-MP-0 small in Modal.
- `resume_v2_modal.py`: resume from existing Modal checkpoint.
- `test_v2_modal.py`: independent held-out evaluation (overall and atom-role-specific RMSE).

**Critical provenance warning:** the local `dft.zip` in our all-files archive contains **6,614 frames and 52 motifs**, NOT the **9,413 frames** (8,388 train, 654 validation, 371 held-out test) used for the reported v2 model. The source dataset for v2 was a different/later DFT input archive (referred to as `dft_3.zip` in older scripts), which is NOT included here. Do not claim that this local `dft.zip` exactly reproduces the final trained model. Recover the original 9,413-frame input and its split/manifest before claiming full training reproducibility. Data `mlff_data_v2/` is expected in the repository root but should not be committed without considering size and provenance.

Example for data preparation (ONLY with verified full input archive):

```bash
python ml/mace/prepare_data.py /path/to/full_DFT_archive.zip --output mlff_data_v2
python ml/mace/audit_v2.py mlff_data_v2
modal run ml/mace/train_v2_modal_fixed.py
modal run ml/mace/test_v2_modal.py
```

The `Modal` code accesses `graphene-sac-mace-results` and uses `/results/full_v2_resume/checkpoints/GrapheneSAC_v2_run-42.model` for evaluation. Check license/redistribution terms before uploading the `.model` file.

Reported test metrics from prior completed evaluation (raw CSV/JSON output still to be downloaded from Modal): energy RMSE 4.9018 meV/atom; overall force RMSE 80.36311 meV/Å. Metal-specific force errors can be much larger, especially Zn 421 meV/Å.
