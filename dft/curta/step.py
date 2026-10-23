"""Steps of one run directory on Curta (pure Python 3.6, no packages): relax -> restart -> SCF.

    python3 step.py relaxed    # exit 0 if relax.out finished with final coordinates
    python3 step.py restart    # unfinished relax.out: move it aside, put its last geometry into relax.in
    python3 step.py scf        # write scf.in (3x3x1 k) from relax.in and the final geometry
    python3 step.py done       # exit 0 if scf.out converged and finished (gamma_only runs: if relaxed)
"""
import os
import re
import sys

KPTS_SCF = "3 3 1 0 0 0"


def read(path):
    with open(path) as f:
        return f.read()


def nat(text):
    return int(re.search(r"nat\s*=\s*(\d+)", text).group(1))


def last_positions(out, n):
    """Last ATOMIC_POSITIONS block of a pw.x output (list of 'El x y z' lines) or None."""
    i = out.rfind("ATOMIC_POSITIONS")
    if i < 0:
        return None
    lines = out[i:].splitlines()[1:n + 1]
    return [" ".join(l.split()[:4]) for l in lines]


def replace_positions(inp, pos):
    """New positions into an input, keeping the if_pos flags of the old lines."""
    head, tail = inp.split("ATOMIC_POSITIONS angstrom\n")
    old = tail.splitlines()
    n = len(pos)
    new = []
    for p, o in zip(pos, old[:n]):
        flags = o.split()[4:]
        new.append("  " + p + ("  " + " ".join(flags) if flags else ""))
    return head + "ATOMIC_POSITIONS angstrom\n" + "\n".join(new) + "\n" + "\n".join(old[n:]) + "\n"


def relaxed():
    return os.path.exists("relax.out") and "End final coordinates" in read("relax.out")


def main(cmd):
    if cmd == "relaxed":
        return 0 if relaxed() else 1
    if cmd == "restart":
        if not os.path.exists("relax.out"):
            return 0
        k = 1
        while os.path.exists("relax.out.%d" % k):
            k += 1
        out = read("relax.out")
        os.rename("relax.out", "relax.out.%d" % k)
        inp = read("relax.in")
        pos = last_positions(out, nat(inp))
        if pos:
            with open("relax.in", "w") as f:
                f.write(replace_positions(inp, pos))
            print("restart from the last geometry of relax.out.%d" % k)
        return 0
    if cmd == "scf":
        if not relaxed():
            return 1
        out = read("relax.out")
        final = out[out.rfind("Begin final coordinates"):]
        inp = read("relax.in")
        scf = replace_positions(inp, last_positions(final, nat(inp)))
        scf = scf.replace("calculation = 'relax'", "calculation = 'scf'")
        scf = scf.replace("K_POINTS gamma", "K_POINTS automatic\n  " + KPTS_SCF)
        with open("scf.in", "w") as f:
            f.write(scf)
        return 0
    if cmd == "done":
        if os.path.exists("gamma_only"):
            single = "calculation = 'scf'" in read("relax.in") and os.path.exists("relax.out") \
                and "JOB DONE" in read("relax.out") and "convergence has been achieved" in read("relax.out")
            return 0 if relaxed() or single else 1     # a Gamma single point (NEB images) is done when it converged
        ok = os.path.exists("scf.out") and "JOB DONE" in read("scf.out") \
            and "convergence has been achieved" in read("scf.out")
        return 0 if ok else 1
    print(__doc__)
    return 2


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else ""))
