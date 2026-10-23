"""Independent, motif-held-out DFT test of the trained GrapheneSAC MACE v2 model.

Run from ~/Downloads/modal_sac_pilot/mace_full_modal_step:
    modal run test_v2_modal.py

Inputs:
    Local: mlff_data_v2/test.xyz (371 held-out structures)
    Modal Volume graphene-sac-mace-results:
        full_v2_resume/checkpoints/GrapheneSAC_v2_run-42.model

Outputs (persisted in the SAME Modal Volume, under independent location):
    full_v2_test/summary.json
    full_v2_test/by_structure.csv
    full_v2_test/by_element.csv
    full_v2_test/by_role.csv
    full_v2_test/by_metal.csv

No training; does not access training/validation data; does not overwrite model.
"""

from pathlib import Path
import modal

app = modal.App("graphene-sac-mace-v2-independent-test")
volume = modal.Volume.from_name("graphene-sac-mace-results")
image = (
    modal.Image.from_registry("pytorch/pytorch:2.6.0-cuda12.4-cudnn9-runtime")
    .pip_install("mace-torch==0.3.16", "ase", "numpy")
    .add_local_file("mlff_data_v2/test.xyz", "/data/test.xyz")
)


@app.function(
    image=image,
    gpu="A10G",
    cpu=4,
    memory=16000,
    timeout=3600,
    volumes={"/results": volume},
)
def evaluate():
    import csv
    import json
    import math
    import os
    import time
    from collections import defaultdict

    import numpy as np
    import torch
    from ase.data import chemical_symbols
    from ase.io import iread
    from mace.calculators import MACECalculator

    os.environ["OMP_NUM_THREADS"] = "4"
    torch.set_num_threads(4)
    if not torch.cuda.is_available():
        raise RuntimeError("No CUDA GPU detected")

    model_path = Path("/results/full_v2_resume/checkpoints/GrapheneSAC_v2_run-42.model")
    test_path = Path("/data/test.xyz")
    for path in (model_path, test_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required file not found: {path}")

    print(f"GPU: {torch.cuda.get_device_name(0)}", flush=True)
    print(f"MODEL: {model_path}", flush=True)
    print(f"TEST FILE: {test_path}", flush=True)

    calc = MACECalculator(
        model_paths=str(model_path),
        device="cuda",
        default_dtype="float32",
    )

    # Metrics accumulate squared errors; all forces are 3D Cartesian components.
    buckets = defaultdict(lambda: {"n_components": 0, "sq": 0.0, "sq_reference": 0.0, "abs": 0.0, "n_atoms": 0})
    energy_errors_per_atom = []
    frame_rows = []
    metal_groups = defaultdict(lambda: {"e": [], "fsq": 0.0, "fn": 0})
    metal_set = {24, 25, 26, 27, 28, 29, 30, 42, 44, 74, 78}
    start = time.monotonic()

    def add_forces(key, df, fref):
        df = np.asarray(df, dtype=np.float64).reshape(-1, 3)
        fref = np.asarray(fref, dtype=np.float64).reshape(-1, 3)
        if not len(df):
            return
        b = buckets[key]
        b["n_atoms"] += len(df)
        b["n_components"] += df.size
        b["sq"] += float(np.sum(df ** 2))
        b["abs"] += float(np.sum(np.abs(df)))
        b["sq_reference"] += float(np.sum(fref ** 2))

    for k, atoms in enumerate(iread(str(test_path), index=":"), start=1):
        # Read actual DFT labels BEFORE attaching MACE, otherwise ASE would
        # overwrite the reference SinglePointCalculator with our predictor.
        try:
            reference_energy = float(atoms.get_potential_energy())
            reference_forces = np.asarray(atoms.get_forces(), dtype=np.float64).copy()
        except Exception as exc:
            raise RuntimeError(f"Test frame {k} lacks readable DFT labels: {exc}") from exc

        if not np.isfinite(reference_energy) or not np.isfinite(reference_forces).all():
            raise ValueError(f"Non-finite DFT labels in frame {k}")
        reference_info = dict(atoms.info)
        atoms.calc = calc
        predicted_energy = float(atoms.get_potential_energy())
        predicted_forces = np.asarray(atoms.get_forces(), dtype=np.float64)
        if not np.isfinite(predicted_energy) or not np.isfinite(predicted_forces).all():
            raise ValueError(f"Non-finite model prediction in frame {k}")
        if reference_forces.shape != predicted_forces.shape:
            raise ValueError(f"Force shape mismatch in frame {k}")

        z = atoms.numbers.astype(int)
        n = len(z)
        eerr = (predicted_energy - reference_energy) / n
        energy_errors_per_atom.append(eerr)
        diff = predicted_forces - reference_forces
        add_forces("overall", diff, reference_forces)

        metal_indices = np.where(np.isin(z, list(metal_set)))[0]
        metal_types = sorted(set(z[metal_indices]))
        metal_label = "+".join(chemical_symbols[x] for x in metal_types) or "none"
        mg = metal_groups[metal_label]
        mg["e"].append(eerr)
        mg["fsq"] += float(np.sum(diff ** 2))
        mg["fn"] += diff.size
        add_forces("role:metal" , diff[metal_indices], reference_forces[metal_indices])

        if len(metal_indices):
            # ASE minimum-image distances account for periodic boundary conditions.
            dmin = np.minimum.reduce([
                atoms.get_distances(int(m), np.arange(n), mic=True)
                for m in metal_indices
            ])
        else:
            dmin = np.full(n, float("inf"))

        for element in np.unique(z):
            mask = z == element
            symbol = chemical_symbols[int(element)]
            add_forces(f"element:{symbol}", diff[mask], reference_forces[mask])

            # Geometric role groups are approximate; these do not assert
            # DFT bonding, coordination number or surface adsorption.
            if int(element) in (6, 7, 8, 16):
                radius = {6: 2.6, 7: 2.7, 8: 2.7, 16: 3.1}[int(element)]
                donor = mask & (dmin <= radius)
                add_forces(f"role:near_metal_{symbol}", diff[donor], reference_forces[donor])
                if int(element) == 6:
                    near_c = mask & (dmin <= 4.0)
                    far_c = mask & (dmin > 4.0)
                    add_forces("role:graphene_C_near_4A", diff[near_c], reference_forces[near_c])
                    add_forces("role:graphene_C_far_4A", diff[far_c], reference_forces[far_c])
            elif int(element) == 1:
                near_h = mask & (dmin <= 3.5)
                add_forces("role:H_within_3.5A_of_metal", diff[near_h], reference_forces[near_h])

        frame_rows.append({
            "frame": k,
            "n_atoms": n,
            "metal": metal_label,
            "reference_energy_eV": reference_energy,
            "predicted_energy_eV": predicted_energy,
            "energy_error_meV_per_atom": round(1000 * eerr, 6),
            "force_rmse_meV_per_A": round(1000 * math.sqrt(np.mean(diff ** 2)), 6),
            "metadata_state": str(reference_info.get("state", "")),
            "metadata_motif": str(reference_info.get("motif", "")),
        })
        if k == 1 or k % 50 == 0:
            print(f"Evaluated {k} frames ({time.monotonic()-start:.1f} s)", flush=True)

    def format_bucket(key):
        b = buckets[key]
        denom = b["n_components"]
        if not denom:
            return None
        rmse = 1000 * math.sqrt(b["sq"] / denom)
        reference_rms = 1000 * math.sqrt(b["sq_reference"] / denom)
        return {
            "group": key.split(":", 1)[-1],
            "atom_observations": b["n_atoms"],
            "force_rmse_meV_per_A": round(rmse, 5),
            "force_mae_meV_per_A": round(1000 * b["abs"] / denom, 5),
            "reference_force_rms_meV_per_A": round(reference_rms, 5),
            "relative_force_rmse_pct": round(100 * rmse / reference_rms, 3) if reference_rms > 0 else None,
        }

    if not frame_rows:
        raise ValueError("No test frames found")

    overall = format_bucket("overall")
    summary = {
        "model_volume_path": str(model_path),
        "test_file": str(test_path),
        "frames": len(frame_rows),
        "total_atom_observations": int(overall["atom_observations"]),
        "energy_rmse_meV_per_atom": round(1000 * math.sqrt(np.mean(np.square(energy_errors_per_atom))), 5),
        "energy_mae_meV_per_atom": round(1000 * float(np.mean(np.abs(energy_errors_per_atom))), 5),
        "force_rmse_meV_per_A": overall["force_rmse_meV_per_A"],
        "force_mae_meV_per_A": overall["force_mae_meV_per_A"],
        "relative_force_rmse_pct": overall["relative_force_rmse_pct"],
        "elapsed_seconds": round(time.monotonic() - start, 2),
        "role_classification": "Geometric distance cutoff from metal atom, with minimum-image PBC; approximate only",
    }

    element_rows = [r for key in sorted(buckets) if key.startswith("element:") if (r := format_bucket(key))]
    role_rows = [r for key in sorted(buckets) if key.startswith("role:") if (r := format_bucket(key))]
    metal_rows = []
    for key, group in sorted(metal_groups.items()):
        metal_rows.append({
            "metal": key,
            "frames": len(group["e"]),
            "energy_rmse_meV_per_atom": round(1000 * math.sqrt(np.mean(np.square(group["e"]))), 5),
            "force_rmse_meV_per_A": round(1000 * math.sqrt(group["fsq"] / group["fn"]), 5),
        })

    dest = Path("/results/full_v2_test")
    dest.mkdir(parents=True, exist_ok=True)

    def save_csv(filename, rows):
        if not rows:
            return
        with (dest / filename).open("w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)

    with (dest / "summary.json").open("w") as f:
        json.dump(summary, f, indent=2)
    save_csv("by_structure.csv", frame_rows)
    save_csv("by_element.csv", element_rows)
    save_csv("by_role.csv", role_rows)
    save_csv("by_metal.csv", metal_rows)
    volume.commit()

    print("\n=== INDEPENDENT HELD-OUT TEST ===", flush=True)
    print(json.dumps(summary, indent=2), flush=True)
    print("\n=== FORCE RMSE BY ELEMENT (meV/A) ===", flush=True)
    for r in element_rows:
        print(f"  {r['group']:>4}: {r['force_rmse_meV_per_A']:8.2f}  ({r['atom_observations']} atoms)", flush=True)
    print("\n=== FORCE RMSE BY ROLE (meV/A) ===", flush=True)
    for r in role_rows:
        print(f"  {r['group']:>28}: {r['force_rmse_meV_per_A']:8.2f}  ({r['atom_observations']} atoms)", flush=True)
    print("\n=== ERRORS BY METAL ===", flush=True)
    for r in metal_rows:
        print(f"  {r['metal']:>5}: {r['frames']:3d} frames, F_RMSE={r['force_rmse_meV_per_A']:.2f} meV/A", flush=True)
    print("\nSaved JSON/CSV results to graphene-sac-mace-results/full_v2_test", flush=True)


@app.local_entrypoint()
def main():
    test_file = Path("mlff_data_v2/test.xyz")
    if not test_file.is_file():
        raise FileNotFoundError(
            "Missing mlff_data_v2/test.xyz. Run this from mace_full_modal_step/ "
            "where mlff_data_v2/ is located."
        )
    evaluate.remote()
