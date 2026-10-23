"""pw.x inputs for the single-atom-site screen (GBRV ultrasoft pseudopotentials, PBE-D3(BJ), spin-polarised).

Settings, and why:
- GBRV v1.5 (Garrity et al., Comput. Mater. Sci. 81, 446 (2014)): one family for every metal, tested at
  40/200 Ry, so all motifs share one cutoff and the screen stays cheap.
- PBE + D3(BJ), as in our earlier MoS2 runs; Marzari-Vanderbilt smearing 0.01 Ry (graphene is a semimetal).
- Spin-polarised with a starting moment on the metal; nosym so H on top can tilt off the axis.
- Dipole correction in the vacuum (H and axial ligands sit on one side only).
- Relax at Gamma (the 6x6 Gamma point contains graphene's K point), then SCF with 3x3x1 k at the relaxed
  geometry; energies are taken from the SCF.
"""
import numpy as np

PSEUDO_DIR = "~/sac_dft/pseudo_gbrv"
PSEUDO = {   # element: (mass, file, valence electrons)
    "H": (1.008, "h_pbe_v1.4.uspp.F.UPF", 1), "C": (12.011, "c_pbe_v1.2.uspp.F.UPF", 4),
    "N": (14.007, "n_pbe_v1.2.uspp.F.UPF", 5), "O": (15.999, "o_pbe_v1.2.uspp.F.UPF", 6),
    "S": (32.06, "s_pbe_v1.4.uspp.F.UPF", 6), "Se": (78.971, "se_pbe_v1.uspp.F.UPF", 6),
    "Si": (28.085, "si_pbe_v1.uspp.F.UPF", 4),
    "V": (50.942, "v_pbe_v1.4.uspp.F.UPF", 13), "Cr": (51.996, "cr_pbe_v1.5.uspp.F.UPF", 14),
    "Mn": (54.938, "mn_pbe_v1.5.uspp.F.UPF", 15), "Fe": (55.845, "fe_pbe_v1.5.uspp.F.UPF", 16),
    "Co": (58.933, "co_pbe_v1.2.uspp.F.UPF", 17), "Ni": (58.693, "ni_pbe_v1.4.uspp.F.UPF", 18),
    "Cu": (63.546, "cu_pbe_v1.2.uspp.F.UPF", 19), "Zn": (65.38, "zn_pbe_v1.uspp.F.UPF", 20),
    "Mo": (95.95, "mo_pbe_v1.uspp.F.UPF", 14), "W": (183.84, "w_pbe_v1.2.uspp.F.UPF", 14),
    "Ru": (101.07, "ru_pbe_v1.2.uspp.F.UPF", 16), "Pt": (195.084, "pt_pbe_v1.4.uspp.F.UPF", 16),
}
ECUT = (40, 200)          # Ry, GBRV recommendation
MOMENT = {"V": 3, "Cr": 4, "Mn": 4, "Fe": 3, "Co": 2, "Ni": 1, "Cu": 1, "Zn": 1, "Mo": 2, "W": 2, "Ru": 2, "Pt": 1}   # starting moment (muB);
# every metal gets one: odd-electron sites (e.g. Zn-N3 in a single vacancy) need nspin = 2 even for d10 metals
EMAXPOS, EOPREG = 0.85, 0.05     # sawtooth of the dipole correction: 17-18 A in a 20 A cell, in the vacuum
KPTS_SCF = (3, 3, 1)
MAX_SECONDS = 6 * 3600 - 1200    # pw.x stops cleanly 20 min before the 6 h job limit


def pw_input(atoms, calc="relax", fixed=(), prefix="pw", kpts=None, moments=None):
    """pw.x input text; fixed = atom indices kept in place. kpts = k mesh of a bulk cell (no dipole correction,
    symmetry on); moments = starting moments (muB) per element instead of MOMENT (metal-free sites)."""
    species = list(dict.fromkeys(atoms.get_chemical_symbols()))
    moments = MOMENT if moments is None else moments
    mag = {el: moments[el] / PSEUDO[el][2] for el in species if el in moments}
    spin = bool(mag)
    slab = kpts is None
    lines = ["&CONTROL", "  calculation = '%s'" % calc, "  prefix = '%s'" % prefix,
             "  pseudo_dir = '%s'" % PSEUDO_DIR, "  tprnfor = .true.", "  forc_conv_thr = 1.0d-3",
             "  etot_conv_thr = 1.0d-4", "  nstep = 200", "  max_seconds = %d" % MAX_SECONDS]
    if slab:
        lines += ["  tefield = .true.", "  dipfield = .true."]
    lines += ["/", "&SYSTEM", "  ibrav = 0", "  nat = %d" % len(atoms), "  ntyp = %d" % len(species),
              "  ecutwfc = %d" % ECUT[0], "  ecutrho = %d" % ECUT[1],
              "  occupations = 'smearing'", "  smearing = 'mv'", "  degauss = 0.01",
              "  vdw_corr = 'dft-d3'", "  dftd3_version = 4"]
    if slab:
        lines += ["  nosym = .true.", "  edir = 3", "  emaxpos = %.2f" % EMAXPOS, "  eopreg = %.2f" % EOPREG,
                  "  eamp = 0.0"]
    if spin:
        lines.append("  nspin = 2")
        lines += ["  starting_magnetization(%d) = %.3f" % (species.index(el) + 1, m) for el, m in mag.items()]
    # gentle mixing: magnetic Fe/Cu sites oscillated at beta 0.2 and pw.x then aborted the relaxation (batch 1);
    # scf_must_converge = .false. lets a relaxation continue past one unconverged step (energies come from the SCF)
    lines += ["/", "&ELECTRONS", "  conv_thr = 1.0d-7", "  mixing_beta = 0.1", "  mixing_ndim = 16",
              "  scf_must_converge = .false.", "  mixing_mode = 'local-TF'", "  electron_maxstep = 600", "/",
              "&IONS", "  ion_dynamics = 'bfgs'", "/"]
    if calc == "vc-relax":
        lines += ["&CELL", "  cell_dynamics = 'bfgs'", "  press_conv_thr = 0.5", "/"]
    lines += ["", "ATOMIC_SPECIES"]
    lines += ["  %-3s %8.3f  %s" % (el, PSEUDO[el][0], PSEUDO[el][1]) for el in species]
    lines += ["", "CELL_PARAMETERS angstrom"]
    lines += ["  %14.8f %14.8f %14.8f" % tuple(v) for v in atoms.cell]
    lines += ["", "ATOMIC_POSITIONS angstrom"]
    for i, a in enumerate(atoms):
        flags = "  0 0 0" if i in fixed else ""
        lines.append("  %-3s %14.8f %14.8f %14.8f%s" % (a.symbol, *a.position, flags))
    lines.append("")
    if kpts:
        lines += ["K_POINTS automatic", "  %d %d %d 0 0 0" % tuple(kpts)]
    elif calc == "scf":
        lines += ["K_POINTS automatic", "  %d %d %d 0 0 0" % KPTS_SCF]
    else:
        lines += ["K_POINTS gamma"]
    return "\n".join(lines) + "\n"


def anchor(atoms, imetal):
    """Index of the carbon farthest (in plane, minimum image) from the metal: fixed to stop the sheet drifting."""
    d = atoms.get_distances(imetal, range(len(atoms)), mic=True, vector=True)
    r = np.linalg.norm(d[:, :2], axis=1)
    r[[a.symbol != "C" for a in atoms]] = -1
    return int(np.argmax(r))
