"""Optional hackathon exploratory force-field optimisation of 6 second-sphere structures.

    modal run --detach optimise_second_sphere_modal.py

Output to graphene-sac-mace-results/second_sphere_v2/optimised.extxyz .
For a DFT single-point, team's instructions use:
    python dft/make_inputs.py points optimised.extxyz ss1

WARNING: isolated H-bonded organic molecules were not included in MACE v2 training;
         single-point DFT verification is needed for meaningful conclusions.
"""
from pathlib import Path
import modal

app = modal.App("graphene-sac-second-sphere-v2")
volume = modal.Volume.from_name("graphene-sac-mace-results", create_if_missing=False)
image = (
    modal.Image.from_registry("pytorch/pytorch:2.6.0-cuda12.4-cudnn9-runtime")
    .pip_install("mace-torch==0.3.16", "ase", "numpy")
    .add_local_file("second_sphere/second_sphere_start.extxyz", "/app/second_sphere_start.extxyz")
)
MODEL = "/results/full_v2_resume/checkpoints/GrapheneSAC_v2_run-42.model"
OUT = "/results/second_sphere_v2"


@app.function(image=image, gpu="A10G", cpu=4, memory=16000,
              timeout=14400, volumes={"/results": volume})
def optimise(fmax, steps):
    import json
    import numpy as np
    from pathlib import Path
    from ase.io import read, write
    from ase.constraints import FixAtoms
    from ase.optimize import FIRE
    from mace.calculators import MACECalculator

    out = Path(OUT)
    out.mkdir(parents=True, exist_ok=True)
    if not Path(MODEL).is_file():
        raise FileNotFoundError(MODEL)
    calc = MACECalculator(model_paths=MODEL, device="cuda", default_dtype="float32")
    results = []
    complete = []
    starts = read("/app/second_sphere_start.extxyz", index=":")
    for a0 in starts:
        atoms = a0.copy()
        motif = str(a0.info["motif"])
        state = str(a0.info["state"])
        # Freezing one heavy atom per molecule approximates enzyme scaffold anchoring.
        # Unlike full molecule fixation, the other atoms can still relax.
        tags = atoms.get_tags()
        anchors = [int(a0.info["fixed_index"])]
        for tag in sorted(set(tags)):
            if tag == 0:
                continue
            heavy = [i for i in range(len(atoms))
                     if tags[i] == tag and atoms[i].symbol != "H"]
            if heavy:
                anchors.append(heavy[0])
        atoms.set_constraint(FixAtoms(sorted(set(anchors))))
        atoms.calc = calc
        folder_name = (motif + "_" + state).replace(":", "_")
        folder = out / folder_name
        folder.mkdir(exist_ok=True)
        opt = FIRE(atoms, logfile=str(folder / "optimisation.log"))
        status = {"motif": motif, "state": state, "fixed_atoms": anchors}
        try:
            e_before = float(atoms.get_potential_energy())
            opt.run(fmax=fmax, steps=steps)
            e_after = float(atoms.get_potential_energy())
            max_free_force = float(np.linalg.norm(atoms.get_forces(), axis=1).max())
            status.update({"E_before_MACE_eV": e_before,
                           "E_after_MACE_eV": e_after,
                           "max_force_eV_A": max_free_force,
                           "steps": opt.nsteps,
                           "converged": max_free_force <= fmax})
            # Store a single structure per motif/state, preserving tags + info fields.
            saved = atoms.copy()
            saved.calc = None
            saved.set_constraint()  # no frozen DOFs in subsequent DFT calculations
            saved.info["mace_E_after_eV"] = e_after
            complete.append(saved)
            write(str(folder / "optimised.extxyz"), saved, format="extxyz")
        except Exception as exc:
            status["error"] = type(exc).__name__ + ": " + str(exc)
        finally:
            results.append(status)
            (folder / "status.json").write_text(json.dumps(status, indent=2))
            volume.commit()
            print("Second-sphere", status, flush=True)
    if complete:
        write(str(out / "optimised.extxyz"), complete, format="extxyz")
    (out / "summary.json").write_text(json.dumps(results, indent=2))
    volume.commit()
    return results


@app.local_entrypoint()
def main(fmax: float = 0.10, steps: int = 150):
    print(optimise.remote(fmax, steps))
