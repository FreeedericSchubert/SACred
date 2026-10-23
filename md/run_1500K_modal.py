"""Modal launcher for team's three DFT-relaxed SACs at 1500 K.

Run inside this bundle directory:
    modal run --detach run_1500K_modal.py --mode smoke
    modal run --detach run_1500K_modal.py --mode production
    modal run --detach run_1500K_modal.py --mode custom --zn-ps 10 --ni-ps 200 --co-ps 10
    modal run --detach run_1500K_modal.py --mode custom --replicas 3 --zn-ps 10 --ni-ps 20 --co-ps 10

Smoke trajectories can be extended automatically: same T, dt, replica and output path.
All outputs go to Volume graphene-sac-mace-results/pyrolysis_1500K.
"""
from pathlib import Path
import modal

app = modal.App("graphene-sac-1500k-escape")
volume = modal.Volume.from_name("graphene-sac-mace-results", create_if_missing=False)
image = (
    modal.Image.from_registry("pytorch/pytorch:2.6.0-cuda12.4-cudnn9-runtime")
    .pip_install("mace-torch==0.3.16", "ase", "numpy")
    .add_local_dir("md_start", "/app/md_start")
)

MODEL = "/results/full_v2_resume/checkpoints/GrapheneSAC_v2_run-42.model"
STARTS = "/app/md_start/md_start_sites.extxyz"
OUT = "/results/pyrolysis_1500K"


@app.function(image=image, gpu="A10G", cpu=4, memory=16000,
              timeout=14400, volumes={"/results": volume})
def run_site(site, target_ps, replica, temp, dt_fs, path_check):
    import importlib.util
    import os
    import sys
    from pathlib import Path
    from ase.io import read
    import torch
    assert torch.cuda.is_available(), "Expected NVIDIA GPU for MACE"
    if not Path(MODEL).is_file():
        raise FileNotFoundError("Trained model absent from Modal Volume: " + MODEL)
    sys.path.insert(0, "/app/md_start")
    import md_escape_mace

    os.environ["OMP_NUM_THREADS"] = "4"
    atoms_list = read(STARTS, index=":")
    matches = [a for a in atoms_list if str(a.info.get("motif")) == site]
    if len(matches) != 1:
        raise RuntimeError("Expected one starting structure for {}, found {}".format(site, len(matches)))
    print("GPU:", torch.cuda.get_device_name(0), "site:", site, flush=True)
    calc = md_escape_mace.mace(MODEL, "cuda")
    Path(OUT).mkdir(parents=True, exist_ok=True)
    if path_check and site == "Zn-NS2:SV" and replica == 0:
        md_escape_mace.check(calc, save_csv=str(Path(OUT) / "dft_path_check.csv"))
        volume.commit()
    result = md_escape_mace.run(
        matches[0], temp, target_ps, replica, calc, OUT,
        dt=dt_fs, checkpoint_ps=0.5, trajectory_ps=0.05,
        diagnostic_ps=0.01, on_checkpoint=volume.commit)
    volume.commit()
    print("FINISHED", result, flush=True)
    return result


@app.local_entrypoint()
def main(mode: str = "smoke", replicas: int = 1, zn_ps: float = 10.0,
         ni_ps: float = 20.0, co_ps: float = 5.0, temp: int = 1500,
         dt_fs: float = 0.25, check_dft_paths: bool = False):
    """One A10G per simultaneous site; a few dollars may accrue on long Ni runs."""
    if mode not in ("smoke", "production", "custom"):
        raise ValueError("Choose --mode smoke, production or custom")
    if mode == "smoke":
        durations = {"Zn-NS2:SV": 0.25,
                     "Ni-NS2:SV": 0.25,
                     "Co-N2OC:DV:NONC": 0.25}
    else:
        durations = {"Zn-NS2:SV": zn_ps,
                     "Ni-NS2:SV": ni_ps,
                     "Co-N2OC:DV:NONC": co_ps}
    if replicas < 1 or dt_fs <= 0 or any(d <= 0 for d in durations.values()):
        raise ValueError("Replicas, timestep and durations must be positive")
    jobs = [(site, ps, replica, temp, dt_fs, check_dft_paths)
            for site, ps in durations.items() for replica in range(replicas)]
    print("Submitting {} trajectory jobs at {} K ({} fs timestep)".format(
        len(jobs), temp, dt_fs))
    print("Durations:", durations)
    for result in run_site.starmap(jobs):
        print("RESULT:", result)
