"""Run 1500 K Zn-NS2 MACE MD from a BONDED Zn geometry (xTB NEB image 0).

Place alongside run_1500K_modal.py in mace_1500K_ready/ (with md_start/).
This uses the team's existing md_escape_mace.py, trained MACE v2 and donor criteria.
Outputs to a NEW Volume directory, pyrolysis_1500K_znbonded, to avoid
mistaking the already-dissociated DFT seed run for a new escape event.

Run: modal run --detach run_zn_bonded_modal.py --target-ps 5
Read saved data: modal volume ls graphene-sac-mace-results pyrolysis_1500K_znbonded
"""
from pathlib import Path
import modal

app = modal.App("graphene-sac-zn-bonded-1500k")
volume = modal.Volume.from_name("graphene-sac-mace-results", create_if_missing=False)
image = (
    modal.Image.from_registry("pytorch/pytorch:2.6.0-cuda12.4-cudnn9-runtime")
    .pip_install("mace-torch==0.3.16", "ase", "numpy")
    .add_local_dir("md_start", "/app/md_start")
)
MODEL = "/results/full_v2_resume/checkpoints/GrapheneSAC_v2_run-42.model"
OUT = "/results/pyrolysis_1500K_znbonded"


@app.function(image=image, gpu="A10G", cpu=4, memory=16000, timeout=14400,
              volumes={"/results": volume})
def run_zn(target_ps: float, replica: int):
    import sys
    import numpy as np
    from ase.io import read
    import torch

    sys.path.insert(0, "/app/md_start")
    import md_escape_mace

    if not torch.cuda.is_available():
        raise RuntimeError("CUDA GPU unavailable")
    if not Path(MODEL).is_file():
        raise FileNotFoundError(MODEL)

    # Preserve original metadata (motif, metal_index, donor_indices,
    # fixed_index) while replacing the already-dissociated Zn positions.
    seed = read("/app/md_start/md_start_sites.extxyz", index=0)
    bonded = read(
        "/app/md_start/escape_paths/xtb_neb_path_Zn-NS2_SV.extxyz", index=0
    )
    if seed.get_chemical_symbols() != bonded.get_chemical_symbols():
        raise ValueError("Zn seed and NEB image have different atom ordering")
    seed.set_cell(bonded.cell, scale_atoms=False)
    seed.set_positions(bonded.positions)
    seed.info["source"] = "xTB NEB image 0 (bonded Zn start; NOT the detached DFT-relaxed seed)"

    im = int(seed.info["metal_index"])
    donors = md_escape_mace.parse_atom_indices(seed.info["donor_indices"])
    distances = seed.get_distances(im, donors, mic=True)
    cutoffs = np.array([
        md_escape_mace.CUTOFF.get(seed[i].symbol, 2.6) + md_escape_mace.LOOSE
        for i in donors
    ])
    bound = int((distances < cutoffs).sum())
    print("Zn bonded start: distances (A):", distances.round(3).tolist(),
          "donor cutoffs:", cutoffs.tolist(), "bound donors:", bound, flush=True)
    if bound < 2:
        raise RuntimeError("Starting geometry already counts as escaped")

    calc = md_escape_mace.mace(MODEL, "cuda")
    Path(OUT).mkdir(parents=True, exist_ok=True)
    result = md_escape_mace.run(
        seed, 1500, target_ps, replica, calc, OUT,
        dt=0.25, checkpoint_ps=0.5, trajectory_ps=0.05,
        diagnostic_ps=0.01, on_checkpoint=volume.commit
    )
    volume.commit()
    print("FINISHED", result, flush=True)
    return result


@app.local_entrypoint()
def main(target_ps: float = 5.0, replica: int = 0):
    if target_ps <= 0 or replica < 0:
        raise ValueError("target_ps > 0 and replica >= 0 required")
    result = run_zn.remote(target_ps, replica)
    print("RESULT", result)
