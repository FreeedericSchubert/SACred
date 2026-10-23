"""H on the basal S of clean MoS2, and the H2 reference in the same cell.

    python3 mos2_h.py        # on Curta, in ~/sac_dft: writes runs/mos2/<step>/relax.in, then bash submit.sh runs/mos2

Reuses an earlier MoS2 project on Curta unchanged (QE_MOS2, pure Python 3.6): its pw.x writer, the pslibrary PAW
set (50/400 Ry, PBE-D3(BJ), 6x6 MoS2, Gamma, Mo fixed) and its relaxed MoS2 (step 02). H starts on the top S at the
site stored in that project's meta.json. dE_H = E(H on MoS2) - E(MoS2) - E(H2)/2. Energies come from the Gamma
relaxations (marker file gamma_only: the 6x6 Gamma point holds the K valley).
"""
import json
import os
import sys

QE_MOS2 = os.path.expanduser("~/<earlier_mos2_project>")      # the earlier MoS2 project on Curta
sys.path.insert(0, QE_MOS2)
from qe_inputs import final_positions, sawtooth, write_pw   # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "runs", "mos2")
D_SH = 1.36             # starting S-H distance (A)
SPIN = {"H": 0.5}       # H radical: its electron goes into MoS2
RESOURCES = "NTASKS=24\nMEM=64G\nHOURS=20\n"     # read by submit.sh


def relaxed(step):
    return final_positions(os.path.join(QE_MOS2, step, "pw.out"), converged=True)


def run(name, atoms, cell, spin=SPIN):
    d = os.path.join(OUT, name)
    if os.path.exists(os.path.join(d, "relax.in")):
        print(name, "exists")
        return
    write_pw(os.path.join(d, "relax.in"), "h", "relax", cell, atoms, "gamma",
             emaxpos=sawtooth([a[3] for a in atoms], cell[2][2]), fixed=("Mo",), conv_thr="1.0d-7",
             max_hours=20, pseudo_dir=os.path.join(QE_MOS2, "pseudo"), magnetization=spin)
    with open(os.path.join(d, "gamma_only"), "w") as f:
        f.write("energies from the Gamma relaxation\n")
    with open(os.path.join(d, "resources"), "w") as f:
        f.write(RESOURCES)
    print(name, "written")


def main():
    meta = json.load(open(os.path.join(QE_MOS2, "meta.json")))
    cell = meta["flat"]["cell"]
    z_mo = meta["methyl"]["z_mo"]
    mos2 = relaxed("02_mos2_flatcell")
    sx, sy, _ = meta["methyl"]["site"]
    s = min((a for a in mos2 if a[0] == "S" and a[3] > z_mo), key=lambda a: (a[1] - sx) ** 2 + (a[2] - sy) ** 2)
    h = [("H", s[1], s[2], s[3] + D_SH)]
    centre = [(cell[0][i] + cell[1][i]) / 2 for i in (0, 1)]
    run("30_h2_alone", [("H", centre[0], centre[1], s[3] + 1.0), ("H", centre[0], centre[1], s[3] + 1.75)],
        cell, spin=None)
    run("31_mos2_h", mos2 + h, cell)


if __name__ == "__main__":
    main()
