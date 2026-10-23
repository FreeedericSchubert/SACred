"""Fe single-atom site on clean MoS2, bound to three S like an Fe of an Fe-S cluster: the step before doping
from underneath (charged plate).

    python3 mos2_fe.py       # on Curta, in ~/sac_dft: writes runs/mos2/<step>/relax.in, then bash submit.sh runs/mos2

Steps: 33 Fe above a Mo (top-Mo site), 33h Fe over the empty hexagon (hollow site, the other three-S face),
34 H on the Fe of 33; 35/36 the same as 33/34 with the sheet charged against a plate below (QE gate, the plate,
barrier and dipole layer of the gate setup in the earlier MoS2 project, meta.json "gate"), q_<q> = extra electrons per
6x6 cell (> 0: donor-like, < 0: acceptor-like). Their energies need the gate correction of qe_inputs.py
(dft/analyse.py applies it). Same settings and relaxed MoS2 (step 02) as mos2_h.py. The pslibrary PAW set has no Fe, so
Fe is the GBRV ultrasoft pseudopotential of the graphene runs (it needs 40/200 Ry, here 50/400); pseudo_mos2/
links both sets. The Fe sits over the Mo next to the S that carries H in run 31. Gentler SCF mixing than the writer's default, as for the
magnetic graphene sites.
"""
import math
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mos2_h                                    # noqa: E402  (puts qe_inputs on the path)
import qe_inputs                                 # noqa: E402

GBRV_FE = os.path.join(HERE, "pseudo_gbrv", "fe_pbe_v1.5.uspp.F.UPF")
PSEUDO_DIR = os.path.join(HERE, "pseudo_mos2")
D_FES, D_FEH = 2.2, 1.55     # starting Fe-S and Fe-H distances (A)
SPIN = {"Fe": 0.2}           # ~3 muB of the 16 valence electrons
GATE_Q = (-0.3, 0.3)         # extra electrons per cell, ~1e13 cm^-2
MIXING = ("  mixing_beta = 0.3\n  electron_maxstep = 200\n",
          "  mixing_beta = 0.1\n  mixing_ndim = 16\n  mixing_mode = 'local-TF'\n  scf_must_converge = .false.\n"
          "  electron_maxstep = 600\n")

qe_inputs.PSEUDO["Fe"] = (55.845, os.path.basename(GBRV_FE))
qe_inputs.ORDER.append("Fe")
qe_inputs.ZVAL["Fe"] = 16


def link_pseudos():
    os.makedirs(PSEUDO_DIR, exist_ok=True)
    paw = os.path.join(mos2_h.QE_MOS2, "pseudo")
    for src in [os.path.join(paw, f) for f in os.listdir(paw)] + [GBRV_FE]:
        dst = os.path.join(PSEUDO_DIR, os.path.basename(src))
        if not os.path.exists(dst):
            os.symlink(src, dst)


def lateral(a, b, cell):
    """In-plane vector b - a under the minimum image (cell: 2D lattice in the first two rows)."""
    (ax, ay), (bx, by) = cell[0][:2], cell[1][:2]
    det = ax * by - ay * bx
    dx, dy = b[1] - a[1], b[2] - a[2]
    f1, f2 = (dx * by - dy * bx) / det, (ax * dy - ay * dx) / det
    f1, f2 = f1 - round(f1), f2 - round(f2)
    return f1 * ax + f2 * bx, f1 * ay + f2 * by


def run(name, atoms, cell, gate=None):
    """gate: meta.json "gate" plus the charge (extra electrons per cell); the dipole layer then sits in the barrier."""
    d = os.path.join(mos2_h.OUT, name)
    path = os.path.join(d, "relax.in")
    if os.path.exists(path):
        print(name, "exists")
        return
    emaxpos = gate["emaxpos"] if gate else qe_inputs.sawtooth([a[3] for a in atoms], cell[2][2])
    qe_inputs.write_pw(path, "fe", "relax", cell, atoms, "gamma", emaxpos=emaxpos, fixed=("Mo",),
                       conv_thr="1.0d-7", max_hours=20, pseudo_dir=PSEUDO_DIR, magnetization=SPIN, gate=gate)
    text = open(path).read()
    assert MIXING[0] in text
    with open(path, "w") as f:
        f.write(text.replace(*MIXING))
    with open(os.path.join(d, "gamma_only"), "w") as f:
        f.write("energies from the Gamma relaxation\n")
    with open(os.path.join(d, "resources"), "w") as f:
        f.write(mos2_h.RESOURCES)
    print(name, "written")


def main():
    import json
    link_pseudos()
    meta = json.load(open(os.path.join(mos2_h.QE_MOS2, "meta.json")))
    cell = meta["flat"]["cell"]
    z_mo = meta["methyl"]["z_mo"]
    mos2 = mos2_h.relaxed("02_mos2_flatcell")
    sx, sy, _ = meta["methyl"]["site"]
    s = min((a for a in mos2 if a[0] == "S" and a[3] > z_mo), key=lambda a: (a[1] - sx) ** 2 + (a[2] - sy) ** 2)
    mo = min((a for a in mos2 if a[0] == "Mo"), key=lambda a: sum(v * v for v in lateral(s, a, cell)))
    vx, vy = lateral(mo, s, cell)                                   # Mo -> S, ~1.83 A
    z_fe = s[3] + math.sqrt(D_FES ** 2 - vx * vx - vy * vy)
    top = ("Fe", mo[1], mo[2], z_fe)
    hollow = ("Fe", mo[1] - vx, mo[2] - vy, z_fe)                   # across the Mo from that S
    run("33_mos2_fe", mos2 + [top], cell)
    run("33h_mos2_fe_hollow", mos2 + [hollow], cell)
    run("34_mos2_fe_h", mos2 + [top, ("H", top[1], top[2], z_fe + D_FEH)], cell)
    for q in GATE_Q:
        gate = dict(meta["gate"], charge=q)
        run(os.path.join("35_gate_mos2_fe", qe_inputs.gate_dir(q)), mos2 + [top], cell, gate)
        run(os.path.join("36_gate_mos2_fe_h", qe_inputs.gate_dir(q)), mos2 + [top, ("H", top[1], top[2], z_fe + D_FEH)],
            cell, gate)


if __name__ == "__main__":
    main()
