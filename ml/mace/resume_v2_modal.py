"""Resume the interrupted GrapheneSAC_v2 training from saved Modal checkpoint.

Usage after preparing mlff_data_v2/train.xyz and mlff_data_v2/valid.xyz:
    modal run --detach resume_v2_modal.py

Reads old checkpoints from /results/full_v2/checkpoints and writes to /results/full_v2_resume. Original training files are not overwritten.
Does not load the held-out test set into the training image.
Cloud GPU usage is billable.
"""
from pathlib import Path
import modal

app = modal.App("graphene-sac-mace-v2-resume")
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
    source_checkpoints = Path("/results/full_v2/checkpoints")
    saved = list(source_checkpoints.glob("GrapheneSAC_v2*.pt"))
    if not saved:
        raise FileNotFoundError(
            "No GrapheneSAC_v2 .pt checkpoint found under "
            "/results/full_v2/checkpoints. First run: "
            "modal volume ls graphene-sac-mace-results full_v2/checkpoints/"
        )
    shutil.copytree(source_checkpoints, work / "checkpoints", dirs_exist_ok=True)
    print("Found saved checkpoints:", [p.name for p in saved], flush=True)
    print("Resuming with --restart_latest, maximum TOTAL epochs=12", flush=True)

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
        "--restart_latest",
    ]
    print("Resuming full v2; no held-out test data loaded", flush=True)
    try:
        subprocess.run(cmd, cwd=work, check=True)
    finally:
        dest = Path("/results/full_v2_resume")
        dest.mkdir(parents=True, exist_ok=True)
        for dname in ("models", "logs", "checkpoints", "results"):
            src = work / dname
            if src.exists():
                shutil.copytree(src, dest / dname, dirs_exist_ok=True)
        volume.commit()
        print("Saved training artifacts in graphene-sac-mace-results/full_v2_resume", flush=True)

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
