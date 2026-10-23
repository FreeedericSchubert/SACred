"""Train MACE small on full motif-disjoint train split (Modal A10G).

Usage after preparing mlff_data_v2/train.xyz and mlff_data_v2/valid.xyz:
    modal run train_v2_modal_fixed.py

Uses a NEW /results/full_v2 destination; does not overwrite pilot checkpoints.
Does not load the held-out test set into the training image.
Cloud GPU usage is billable.
"""
from pathlib import Path
import modal

app = modal.App("graphene-sac-mace-full-v2")
volume = modal.Volume.from_name("graphene-sac-mace-results", create_if_missing=True)
image = (
    modal.Image.from_registry("pytorch/pytorch:2.6.0-cuda12.4-cudnn9-runtime")
    .pip_install("mace-torch>=0.3.13,<0.4", "ase")
    .add_local_file("mlff_data_v2/train.xyz", "/data/train.xyz")
    .add_local_file("mlff_data_v2/valid.xyz", "/data/valid.xyz")
)

@app.function(image=image, gpu="A10G", cpu=4, memory=16000,
              timeout=14400, volumes={"/results": volume})
def train():
    import os
    import shutil
    import subprocess
    from pathlib import Path
    import torch

    for remote_file in ("/data/train.xyz", "/data/valid.xyz"):
        if not Path(remote_file).is_file():
            raise FileNotFoundError(f"Missing mounted training data: {remote_file}")
    assert torch.cuda.is_available(), "CUDA GPU unavailable"
    print("GPU:", torch.cuda.get_device_name(0), flush=True)
    os.environ["OMP_NUM_THREADS"] = "4"
    os.environ["MKL_NUM_THREADS"] = "4"
    work = Path("/tmp/mace_full_v2")
    work.mkdir(parents=True, exist_ok=True)

    cmd = [
        "mace_run_train",
        "--name=GrapheneSAC_v2",
        "--foundation_model=small",
        "--multiheads_finetuning=False",
        "--train_file=/data/train.xyz",
        "--valid_file=/data/valid.xyz",
        "--energy_key=energy",
        "--forces_key=forces",
        "--E0s=estimated",
        "--energy_weight=10.0",
        "--forces_weight=10.0",
        "--lr=0.001",
        "--weight_decay=0.0",
        "--scaling=rms_forces_scaling",
        "--batch_size=1",
        "--valid_batch_size=1",
        "--max_num_epochs=12",
        "--patience=5",
        "--ema",
        "--ema_decay=0.999",
        "--amsgrad",
        "--clip_grad=1.0",
        "--default_dtype=float32",
        "--device=cuda",
        "--seed=42",
    ]
    print("Training full v2; no held-out test data loaded", flush=True)
    try:
        subprocess.run(cmd, cwd=work, check=True)
    finally:
        dest = Path("/results/full_v2")
        dest.mkdir(parents=True, exist_ok=True)
        for dname in ("models", "logs", "checkpoints", "results"):
            src = work / dname
            if src.exists():
                shutil.copytree(src, dest / dname, dirs_exist_ok=True)
        volume.commit()
        print("Saved training artifacts in graphene-sac-mace-results/full_v2", flush=True)

@app.local_entrypoint()
def main():
    # Check local data only on the user's computer, not when the module
    # is imported again inside Modal's remote container.
    for f in ("mlff_data_v2/train.xyz", "mlff_data_v2/valid.xyz"):
        if not Path(f).is_file():
            raise FileNotFoundError(
                f"Missing local {f}. Run prepare_data.py first or change to its folder."
            )
    train.remote()
