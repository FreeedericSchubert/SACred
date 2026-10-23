"""The MoS2 part in one figure: the structures computed so far and how doping from below is modelled.

    python figures/mos2_views.py          # -> figures/mos2.png

a, c, d are side views of the 6x6 cells along the zigzag direction (x), all at the same height scale: H on a
basal S of clean MoS2 (31), the Fe-S3 site with H on Fe (34) and that site with the charged plate of the QE gate
(positions of the gate setup in the earlier MoS2 project, its meta.json "gate"). Each geometry is the last step of
its relaxation, or the starting geometry before the run starts. b: top views of the two Fe sites (33 above Mo, 33h
over the hollow). e: the potential energy of an electron along z from the plate and its counter charge in the sheet, for 0.3 e per cell
(barrier, field between plate and sheet, dipole layer; the atoms' own potential is left out).
"""
import os
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from ase.io import read
from matplotlib.colors import to_rgb
from matplotlib.patches import Circle, Rectangle

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RUNS = os.path.join(ROOT, "dft", "runs", "mos2")
DGH = os.path.join(ROOT, "data", "dgh_mos2.csv")
C_CELL = 35.0                       # A
GATE = {"zgate": 0.25716, "block_1": 0.20216, "block_2": 0.31431, "eopreg": 0.05}   # fractions of c
BARRIER = 0.2 * 13.605693           # eV (block_height 0.2 Ry)
Q = 0.3                             # e per 6x6 cell in the sketch
EPS0 = 1 / (4 * np.pi * 14.399645)  # e^2 / (eV A)
Z_MO = 15.07                        # Mo plane (A), the centre of the sheet's doping charge
Z_RANGE = (6.0, 21.2)
COLOURS = {"Mo": "#8a68b0", "S": "#e6c229", "Fe": "#c4682f", "H": "#f4c3cf", "C": "#7b5a3c", "N": "#7f93c4"}
RADII = {"Mo": 0.62, "S": 0.52, "Fe": 0.62, "H": 0.32, "C": 0.36, "N": 0.36}     # A, drawing only
BONDS = {("Mo", "S"): 2.7, ("Fe", "S"): 2.6, ("Fe", "H"): 1.9, ("H", "S"): 1.6, ("C", "C"): 1.75,
         ("C", "N"): 1.6, ("C", "H"): 1.2}
PLATE = "#b2306a"
INK, MUTED = "#222222", "#555555"


def structure(run):
    """(atoms, 'relaxed' | 'relaxation running' | 'starting geometry')."""
    out = os.path.join(RUNS, run, "relax.out")
    if os.path.exists(out) and "ATOMIC_POSITIONS" in open(out).read():
        state = "relaxed" if "End final coordinates" in open(out).read() else "relaxation running"
        return read(out, format="espresso-out", index=-1), state
    return read(os.path.join(RUNS, run, "relax.in"), format="espresso-in"), "starting geometry"


def shade(colour, f):
    """Colour mixed with white by f (0 = unchanged): atoms further back are lighter."""
    return tuple(c + (1 - c) * f for c in to_rgb(colour))


def bond_vectors(atoms):
    """(i, j, minimum-image vector i -> j) for every bonded pair, both directions."""
    sym = atoms.get_chemical_symbols()
    out = []
    for i in range(len(atoms)):
        d = atoms.get_distances(i, range(len(atoms)), mic=True, vector=True)
        for j in range(len(atoms)):
            cut = BONDS.get(tuple(sorted((sym[i], sym[j]))))
            if j != i and cut and np.linalg.norm(d[j]) < cut:
                out.append((i, j, d[j]))
    return out


def side_view(ax, atoms, labels=()):
    """Projection along x: y to the right, z up; atoms drawn back to front, periodic images at the y edges."""
    atoms = atoms.copy()
    atoms.wrap()
    ly = atoms.cell[1][1]
    pos, sym = atoms.positions, atoms.get_chemical_symbols()
    for i, j, d in bond_vectors(atoms):
        for shift in (-ly, 0, ly):
            ax.plot([pos[i, 1] + shift, pos[i, 1] + shift + d[1] / 2], [pos[i, 2], pos[i, 2] + d[2] / 2],
                    color="#8c8c8c", lw=1.5, solid_capstyle="butt", zorder=2)
    x = pos[:, 0]
    depth = (x - x.min()) / max(np.ptp(x), 1e-6)
    for k in np.argsort(x):
        for shift in (-ly, 0, ly):
            y = pos[k, 1] + shift
            if -1.5 < y < ly + 1.5:
                ax.add_patch(Circle((y, pos[k, 2]), RADII[sym[k]], facecolor=shade(COLOURS[sym[k]], 0.35 * (1 - depth[k])),
                                    edgecolor=shade("#000000", 0.45), lw=0.5, zorder=3 + depth[k]))
    for k, text, dy in labels:
        ax.text(pos[k, 1] + 0.95, pos[k, 2] + dy, text, fontsize=12, color=INK, va="center", zorder=10,
                bbox=dict(fc="white", ec="none", alpha=0.75, pad=0.6))
    ax.set_xlim(-0.6, ly + 0.6)
    ax.set_ylim(*Z_RANGE)
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.tick_params(labelsize=10, colors=MUTED)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color("#999999")
    ax.set_yticks([8, 12, 16, 20])
    ax.set_ylabel("z (Å)", fontsize=10.5, color=MUTED)


def gate_bands(ax, x0, x1, labels=True):
    """Barrier, dipole layer and charged plate of the QE gate, at the positions of the earlier gate runs."""
    b1, b2, zg = (GATE[k] * C_CELL for k in ("block_1", "block_2", "zgate"))
    dip = GATE["eopreg"] * C_CELL
    ax.add_patch(Rectangle((x0, b1), x1 - x0, b2 - b1, fc="#ececec", ec="none", zorder=0))
    ax.add_patch(Rectangle((x0, b1), x1 - x0, dip, fc="none", ec="#bdbdbd", hatch="////", lw=0, zorder=0.5))
    ax.plot([x0, x1], [zg, zg], color=PLATE, lw=3.2, solid_capstyle="butt", zorder=1)
    if labels:
        style = dict(fontsize=10.5, va="center", zorder=10, bbox=dict(fc="white", ec="none", alpha=0.8, pad=0.5))
        ax.text(x0 + 0.3, zg + 0.62, "plate, −%.1f e" % Q, color=PLATE, **style)
        ax.text(x0 + 0.3, b2 - 0.55, "barrier %.1f eV" % BARRIER, color=MUTED, **style)
        ax.text(x0 + 0.3, b1 + dip / 2, "dipole layer", color=MUTED, **style)
        ax.text(x0 + 0.3, 12.35, "sheet: +%.1f e (holes)" % Q, color=PLATE, **style)


def top_view(ax, atoms, title):
    """Top view of the Fe and its surroundings (minimum image around the Fe), drawn bottom to top."""
    fe = [a.index for a in atoms if a.symbol == "Fe"][0]
    d = atoms.get_distances(fe, range(len(atoms)), mic=True, vector=True)
    keep = [k for k in range(len(atoms)) if np.hypot(*d[k, :2]) < 5.6 and atoms[k].symbol != "H"]
    local = atoms[keep]
    local.positions = d[keep] + atoms.positions[fe]
    local.pbc = False
    p, sym = local.positions - local.positions[[local.get_chemical_symbols().index("Fe")]], local.get_chemical_symbols()
    clip = Circle((0, 0), 4.15, transform=ax.transData)
    for i in range(len(local)):
        for j in range(i + 1, len(local)):
            cut = BONDS.get(tuple(sorted((sym[i], sym[j]))))
            if cut and np.linalg.norm(p[i] - p[j]) < cut:
                ax.plot(p[[i, j], 0], p[[i, j], 1], color="#8c8c8c", lw=1.4, clip_path=clip,
                        zorder=2 if "Fe" not in (sym[i], sym[j]) else 8)
    for k in np.argsort(p[:, 2]):
        atom = Circle(p[k, :2], RADII[sym[k]], facecolor=COLOURS[sym[k]],
                      edgecolor=shade("#000000", 0.45), lw=0.5, zorder=3 + p[k, 2] / 10 + (5 if sym[k] == "Fe" else 0))
        ax.add_patch(atom)
        atom.set_clip_path(clip)
    ax.text(0, -0.05, "Fe", fontsize=10, ha="center", va="center", color="white", fontweight="bold", zorder=20)
    ax.set_xlim(-4.3, 4.3)
    ax.set_ylim(-4.3, 4.3)
    ax.set_aspect("equal")
    ax.set_axis_off()
    ax.add_patch(Circle((0, 0), 4.15, fc="none", ec="#bbbbbb", lw=1, zorder=30))
    ax.set_title(title, fontsize=11.5, color=INK, pad=4)


def potential(q):
    """Electron potential energy (eV) along z from the plate and the sheet's counter charge: zero in the vacuum
    above the sheet, a field between plate and sheet, the dipole layer restoring periodicity (barrier not added)."""
    b1, b2, zg = (GATE[k] * C_CELL for k in ("block_1", "block_2", "zgate"))
    top = b1 + GATE["eopreg"] * C_CELL
    z = np.linspace(*Z_RANGE, 2000)
    u_plate = q / (EPS0 * 18.96 * 16.42) * (Z_MO - zg)             # e * E * d; > 0: plate negative
    u = np.where(z >= Z_MO, 0.0, np.where(z >= zg, u_plate * (Z_MO - z) / (Z_MO - zg), u_plate))
    u = np.where(z < top, np.clip((z - b1) / (top - b1), 0, 1) * u_plate, u)
    return z, u, u_plate


def figure():
    fig = plt.figure(figsize=(13, 9.6))
    boxes = {"a": [0.05, 0.555, 0.27, 0.30],
             "c": [0.05, 0.10, 0.27, 0.30], "d": [0.37, 0.10, 0.27, 0.30], "e": [0.72, 0.10, 0.25, 0.30]}
    ax = {k: fig.add_axes(b) for k, b in boxes.items()}
    dgh = {}
    if os.path.exists(DGH):
        import pandas as pd
        dgh = pd.read_csv(DGH).set_index("run").dG_H.dropna().to_dict()

    a, sa = structure("31_mos2_h")
    side_view(ax["a"], a, [(len(a) - 1, "H", 0.0)])
    d, sd = structure("34_mos2_fe_h")
    fe, hd = ([x.index for x in d if x.symbol == el][0] for el in ("Fe", "H"))
    for k in ("c", "d"):
        side_view(ax[k], d, [(fe, "Fe", 0.0), (hd, "H", 0.0)])
    ly = d.cell[1][1]
    gate_bands(ax["d"], -0.6, ly + 0.6)

    # e: potential energy along z, both signs
    fz = ax["e"]
    b1, b2, zg = (GATE[k] * C_CELL for k in ("block_1", "block_2", "zgate"))
    fz.axhspan(b1, b2, color="#ececec", zorder=0)
    fz.axhspan(13.48, 16.64, color=shade(COLOURS["Mo"], 0.75), zorder=0)
    fz.axhline(zg, color=PLATE, lw=1.5, zorder=1)
    fz.axvline(0, color="#bbbbbb", lw=0.8, zorder=1)
    for q, ls, colour, lab, ha in ((Q, "-", PLATE, "plate −%.1f e\nacceptor-like" % Q, "left"),
                                   (-Q, "--", "#4a4a4a", "plate +%.1f e\ndonor-like" % Q, "right")):
        z, u, u_plate = potential(q)
        fz.plot(u, z, ls=ls, color=colour, lw=2, zorder=3)
        u_mid = np.interp(12.4, z, u)                      # label beside the curve, between plate and sheet
        fz.text(u_mid + (0.12 if q > 0 else -0.12), 12.4, lab, fontsize=10, color=colour, ha=ha, va="center")
    fz.text(0.6, 7.6, "barrier +%.1f eV" % BARRIER, fontsize=10, color=MUTED, va="center")
    fz.text(-1.75, zg + 0.5, "plate", fontsize=10, color=PLATE, va="center")
    fz.text(-1.75, 15.06, "MoS$_2$", fontsize=11, color=COLOURS["Mo"], va="center")
    fz.text(0.12, 19.0, "no field on the\nreactive side", fontsize=10, color=MUTED, va="center")
    fz.set_ylim(*Z_RANGE)
    fz.set_xlim(-1.8, 1.8)
    fz.set_yticks([8, 12, 16, 20])
    fz.set_ylabel("z (Å)", fontsize=10.5, color=MUTED)
    fz.set_xlabel("potential energy of an electron (eV)", fontsize=10.5, color=MUTED)
    fz.tick_params(labelsize=10, colors=MUTED)
    for side in ("top", "right"):
        fz.spines[side].set_visible(False)

    # b: two Fe sites from above
    for k, (run, title) in enumerate((("33_mos2_fe", "Fe above Mo (33)"), ("33h_mos2_fe_hollow", "Fe over hollow (33h)"))):
        top_view(fig.add_axes([0.40 + k * 0.29, 0.565, 0.25, 0.27]), structure(run)[0], title)

    dg = dgh.get("31_mos2_h")
    texts = {
        "a": ("Clean MoS$_2$", "H on a basal S: ΔG$_H$ = %s (%s)." % ("%+.2f eV" % dg if dg is not None else "pending", sa)),
        "b": ("Fe–S$_3$ site, from above", "Fe bound to three S, like an Fe of an Fe–S cluster. Both sites: %s; "
              "the lower one is the reference for ΔG$_H$." % structure("33_mos2_fe")[1].replace("starting geometry", "queued")),
        "c": ("Fe–S$_3$ site, H on Fe", "Clean MoS$_2$: %s, queued." % sd),
        "d": ("Charged plate below", "QE gate for doping from below: ±%.1f e per cell in the sheet (10¹³ cm⁻²), "
              "the opposite charge on a plate behind a barrier. Proposed." % Q),
        "e": ("What an electron feels", "Field only between plate and sheet; the dipole layer keeps the cell "
              "periodic; the barrier keeps electrons off the plate. Sketch for ±%.1f e, without the atoms' own "
              "potential." % Q),
    }
    anchors = {"a": boxes["a"], "b": [0.40, 0.555, 0.55, 0.30], "c": boxes["c"], "d": boxes["d"],
               "e": [boxes["e"][0] - 0.035, *boxes["e"][1:]]}
    for tag, (title, note) in texts.items():
        x, y, w, h = anchors[tag]
        fig.text(x - 0.02, y + h + 0.035, tag, fontsize=17, fontweight="bold", va="bottom", color=INK)
        fig.text(x + 0.005, y + h + 0.035, title, fontsize=13.5, va="bottom", color=INK)
        fig.text(x - 0.02, y - (0.075 if tag == "e" else 0.02), textwrap.fill(note, 50), fontsize=10.5,
                 color="#4a4a4a", va="top")
    handles = [plt.Line2D([], [], marker="o", ls="", ms=9, mfc=COLOURS[e], mec=shade("#000000", 0.45), label=e)
               for e in ("Mo", "S", "Fe", "H")]
    fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=11, bbox_to_anchor=(0.5, -0.025),
               handletextpad=0.2, columnspacing=1.4)
    out = os.path.join(HERE, "mos2.png")
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    figure()
