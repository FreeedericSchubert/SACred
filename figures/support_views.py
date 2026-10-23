"""MoS2 and silica examples: enzyme -> first shell -> the site on the support (GFN2-xTB) -> dG_H of every metal on
that support. Panels a and b are those of figures/enzyme_views.py; panel c renders the relaxed xTB structures of
screen/xtb_supports.py (bare and with H on the metal), panel d plots data/xtb_supports.csv.

    python figures/support_views.py              # figures/nife_mos2.png, figures/ca_silica.png
    python figures/support_views.py ca_silica    # one example

MoS2: the [NiFe] hydrogenase nickel. Its three-S face goes onto the S of the sheet; the fourth Cys has no
partner there, so its place is open for H. Silica: the cobalt form of the gamma-carbonic anhydrase of
Methanosarcina thermophila (1QRE). The three O it binds besides its His (bicarbonate and two waters) take the
places of the two siloxide O and the silanol O of the silsesquioxane model; the three His would have to come from
linkers, in the bare model their side is open for H.
"""
import os
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase.io import read as ase_read, write as ase_write
from matplotlib.lines import Line2D
from matplotlib.patches import ConnectionPatch, Ellipse, FancyArrowPatch

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)
import enzyme_views as ev   # noqa: E402

XTB = os.path.join(ROOT, "data", "screen", "xtb_supports")
TABLE = os.path.join(ROOT, "data", "xtb_supports.csv")
PATCH_RADIUS = 7.2      # A: MoS2 atoms within this lateral distance of the metal (its images sit 9.5 A away)
INK, MUTED, BAND = "#222222", "#5a5a5a", "#e3efe6"
DOT = "#2b5d8a"         # one hue for the dots; raw and offset differ by fill, not by colour

EXAMPLES = {
    "nife_mos2": dict(
        ev.EXAMPLES["nife"], support="MoS2", metal_el="Ni", site="top",
        shell_text="4 Cys, Ni–S 2.18–2.54 Å. On MoS$_2$ three of them (S···S 3.4 Å) map onto three S of the "
                   "sheet, 3.16 Å apart; the fourth Cys has no partner there, so its place is open for H.",
        motif_title="Ni on three S of MoS$_2$",
        motif_text="GFN2-xTB, relaxed in a 3×3 cut of the DFT sheet; Ni over a Mo (the lower of the two sites), "
                   "bare and with H on Ni. A thiol at the sheet could give back the fourth S.",
        chart_title="H on the metal, every metal on MoS$_2$"),
    "ca_silica": dict(
        pdb="1QRE", title="γ-Carbonic anhydrase, Co form", subtitle="PDB 1QRE, 1.46 Å · CO$_2$ hydration",
        chains="homotrimer, one colour per subunit; the site sits between two subunits", assembly=True,
        keep="resn HOH and (all within 2.7 of (chain A and resn CO and resi 214))",
        centre="(chain A and (resn CO+BCT or (resn HOH and (all within 2.7 of (resn CO and resi 214)))))",
        centre_residues="(chain A and resi 81+122) or (chain B and resi 117)",
        metal="chain A and resn CO and resi 214",
        donors="((chain A and resi 81 and name ND1) or (chain A and resi 122 and name NE2) or "
               "(chain B and resi 117 and name NE2) or (chain A and resn BCT and name O1) or "
               "(chain A and resn HOH and resi 322+323))",
        metal_bonds=[("chain A and resn CO and resi 214",
                      "((chain A and resi 81 and name ND1) or (chain A and resi 122 and name NE2) or "
                      "(chain B and resi 117 and name NE2) or (chain A and resn BCT and name O1) or "
                      "(chain A and resn HOH and resi 322+323))")],
        shell_title="First shell of Co: Co–N$_3$O$_3$",
        shell_text="3 His (one from the next subunit), bicarbonate and two waters, Co–O 2.39–2.57 Å. The three O "
                   "are the exchangeable side: on silica they become two Si–O–Co and one Si–OH···Co.",
        support="SiO2", metal_el="Co", site="T7",
        motif_title="Co on silica: silsesquioxane model",
        motif_text="GFN2-xTB, Me$_7$Si$_7$O$_9$(OH)$_3$ with Co in the open corner (two siloxides, one silanol). "
                   "The His side is left open for H here; imidazole linkers could supply it.",
        chart_title="H on the metal, every metal on silica"),
}


# ---------------------------------------------------------------- panels a and b (enzyme_views, plus assemblies)

def load_enzyme(cmd, ex):
    """enzyme_views.load_enzyme, plus: the biological assembly as chains A, B, C, ... and kept waters."""
    if not ex.get("assembly"):
        return ORIGINAL_LOAD(cmd, ex)
    cmd.set("assembly", "1")
    cmd.load(os.path.join(ev.PDB, ex["pdb"] + ".cif.gz"), "asm")
    cmd.set("assembly", "")
    cmd.split_states("asm", prefix="cp")
    copies = sorted(n for n in cmd.get_object_list() if n.startswith("cp"))
    for k, name in enumerate(copies):
        cmd.alter(name, "chain='%s'" % "ABCDEFGH"[k])
    cmd.create("enz", " or ".join(copies))
    for name in ["asm"] + copies:
        cmd.delete(name)
    cmd.remove("enz and hydro")
    cmd.remove("enz and not alt ''+A")
    cmd.alter("enz", "alt=''")
    cmd.remove("enz and solvent and not (%s)" % ex["keep"])
    for a, b in ex["metal_bonds"]:
        cmd.bond("enz and " + a, "enz and " + b)
    return "enz and ((%s) or ((%s) and sidechain))" % (ex["centre"], ex["centre_residues"])


ORIGINAL_LOAD = ev.load_enzyme
ev.load_enzyme = load_enzyme
ev.SHELL_LABEL = ("'O water' if resn == 'HOH' else ('O bicarbonate' if resn == 'BCT' else "
                  "'%s %s%s' % (elem.capitalize(), resn.capitalize(), resi))")


# ---------------------------------------------------------------- panel c: the xTB structures

def mos2_patch(atoms, im):
    """Periodic 3x3 cell -> one round patch around the metal (its images lie outside), metal up."""
    a1, a2 = atoms.cell[0], atoms.cell[1]
    centre = atoms.positions[im, :2]
    pos, sym = [], []
    for i in (-1, 0, 1):
        for j in (-1, 0, 1):
            for at in atoms:
                p = at.position + i * a1 + j * a2
                if np.linalg.norm(p[:2] - centre) <= PATCH_RADIUS:
                    pos.append(p)
                    sym.append(at.symbol)
    from ase import Atoms
    return Atoms(sym, positions=pos)


def upright(atoms, up):
    """Rotate so that `up` points along +y (PyMOL looks down -z with y up when the view is not reset)."""
    up = up / np.linalg.norm(up)
    atoms = atoms.copy()
    atoms.rotate(up, (0, 1, 0), center="COP")
    return atoms


def panel_site(xyz, support, path, tilt):
    atoms = ase_read(xyz)
    im = [a.index for a in atoms if a.symbol not in ("Mo", "S", "Si", "O", "C", "H")][0]
    if support == "MoS2":
        patch = mos2_patch(atoms, im)
        patch = upright(patch, np.array([0.0, 0.0, 1.0]))
    else:
        patch = atoms.copy()
        patch.pbc = False
        cage = [a.index for a in patch if a.symbol == "Si"]
        patch = upright(patch, patch.positions[im] - patch.positions[cage].mean(axis=0))
    patch.pbc = False
    pdb = os.path.join(ev.TMP, "support_site.pdb")
    ase_write(pdb, patch, format="proteindatabank")
    cmd = ev.pm()
    cmd.load(pdb, "site")
    metal = "site and not elem Mo+S+Si+O+C+H"
    cmd.unbond(metal, "site")
    cmd.bond(metal, "site and elem S+O within 2.6 of (%s)" % metal)
    cmd.bond(metal, "site and elem H within 1.8 of (%s)" % metal)
    cmd.unbond("site and elem Mo", "site and elem Mo")
    cmd.hide("everything")
    cmd.show("sticks", "site")
    cmd.show("spheres", "site")
    cmd.set("stick_radius", 0.14)
    cmd.set("sphere_scale", 0.26)
    cmd.set("sphere_scale", 0.42, metal)
    cmd.set("sphere_scale", 0.2, "site and elem H")
    ev.colour_sticks(cmd, "site")
    cmd.color("violetpurple", "site and elem Mo")
    cmd.color("yelloworange", "site and elem S")
    cmd.color("wheat", "site and elem Si")
    cmd.color("gray85", "site and elem H")
    cmd.label(metal, "elem.capitalize()")
    cmd.label("site and elem H and (all within 1.8 of (%s))" % metal, "'H'")
    cmd.set("label_position", (0, 0, 4))
    cmd.set("label_size", 22)
    view = list(cmd.get_view())
    view[:9] = [1, 0, 0, 0, 1, 0, 0, 0, 1]    # no rotation: y (the support normal after upright) points up
    cmd.set_view(view)
    cmd.turn("x", tilt)                      # look a little from above
    cmd.zoom("site", buffer=0.6, complete=1)
    cmd.clip("slab", 80)
    ev.render(cmd, path)


# ---------------------------------------------------------------- panel d: dG_H of every metal

def panel_chart(ax, support, highlight):
    t = pd.read_csv(TABLE)
    h = t[(t.support == support) & (t.state == "H_M")].copy()
    h["key"] = h.dG_H_offset.fillna(h.dG_H)
    h = h.sort_values("key", ascending=False).reset_index(drop=True)
    rows = list(h.itertuples())
    clean = t[(t.support == support) & (t.metal == "clean")]
    lo = min(-3.6 if support == "MoS2" else -1.6, h.dG_H.min() - 0.3)
    hi = max(2.3, h[["dG_H", "dG_H_offset"]].max().max() + 0.3)
    ax.axvspan(-0.2, 0.2, color=BAND, zorder=0, lw=0)
    ax.axvline(0, color="#9aa59d", lw=1, zorder=1)
    for y, r in enumerate(rows):
        if pd.notna(r.dG_H_offset):
            ax.plot([r.dG_H, r.dG_H_offset], [y, y], color="#b8c0c6", lw=1.5, zorder=2)
            ax.plot(r.dG_H_offset, y, "o", ms=9, color=DOT, zorder=3)
        ax.plot(r.dG_H, y, "o", ms=9, mfc="white", mec=DOT, mew=1.8, zorder=4)
    labels = [r.metal for r in rows]
    ticks = list(range(len(rows)))
    if len(clean):                         # the clean sheet as a reference row below the metals
        c = clean.iloc[0]
        ax.plot(c.dG_H, -1.4, "o", ms=9, mfc="white", mec=DOT, mew=1.8, zorder=4)
        ax.plot(c.dG_H_dft, -1.4, "D", ms=7.5, color=INK, zorder=4)
        ax.axhline(-0.7, color="#d0d0d0", lw=0.8)
        labels.insert(0, "S, clean sheet")
        ticks.insert(0, -1.4)
    ax.set_yticks(ticks)
    ax.set_yticklabels(labels, fontsize=11, color=INK)
    for tick, lab in zip(ax.get_yticklabels(), labels):
        if lab == highlight:
            tick.set_fontweight("bold")
    ax.set_xlim(lo, hi)
    ax.set_ylim(min(ticks) - 0.6, len(rows) - 0.4)
    ax.set_xlabel("ΔG$_H$ (eV)", fontsize=11, color=INK)
    ax.tick_params(axis="x", labelsize=10, colors=MUTED)
    ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#e6e6e6", lw=0.8, zorder=0)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#b0b0b0")
    handles = [Line2D([], [], ls="", marker="o", ms=8, mfc="white", mec=DOT, mew=1.8, label="xTB"),
               Line2D([], [], ls="", marker="o", ms=8, color=DOT, label="xTB + graphene offset of the metal")]
    if len(clean):
        handles.append(Line2D([], [], ls="", marker="D", ms=7, color=INK, label="DFT, 6×6 sheet"))
    ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=1, frameon=False,
              fontsize=9.5, handletextpad=0.4, labelcolor=INK)
    return h


# ---------------------------------------------------------------- the figure

def figure(name):
    ex = EXAMPLES[name]
    os.makedirs(ev.TMP, exist_ok=True)
    p = {k: os.path.join(ev.TMP, "%s_%s.png" % (name, k))
         for k in ("cartoon", "centre", "inset", "inset_centre", "shell", "site0", "site1")}
    ev.panel_enzyme(ex, [p["cartoon"], p["centre"], p["inset"], p["inset_centre"]])
    ev.panel_shell(ex, p["shell"])
    sup, m, site = ex["support"], ex["metal_el"], ex["site"]
    for k, state in enumerate(("bare", "H_M")):
        panel_site(os.path.join(XTB, sup, "%s_%s_%s.xyz" % (m, site, state)), sup, p["site%d" % k],
                   tilt=14 if sup == "MoS2" else 18)

    fig = plt.figure(figsize=(20.5, 6.6))
    w, h, y0 = 0.225, 0.62, 0.17
    xs = (0.0, 0.265, 0.53, 0.815)
    ax_a, ax_b = (fig.add_axes([x, y0, w, h]) for x in xs[:2])
    layer = plt.imread(p["centre"])
    ax_a.imshow(plt.imread(p["cartoon"]))
    ax_a.imshow(ev.halo(layer))
    ys, xs_ = np.nonzero(layer[..., 3] > 0.1)
    cx, cy, r = xs_.mean(), ys.mean(), 0.6 * max(np.ptp(xs_), np.ptp(ys)) + 25
    ax_a.add_patch(Ellipse((cx, cy), 2 * r, 2 * r, fill=False, ls="--", lw=1.2, color="#333333"))
    ax_in = fig.add_axes([0.14, 0.52, 0.105, 0.33])
    ax_in.imshow(plt.imread(p["inset"]))
    ax_in.imshow(ev.halo(plt.imread(p["inset_centre"])))
    ax_in.set_xticks([]), ax_in.set_yticks([])
    for sp in ax_in.spines.values():
        sp.set_color("#333333"), sp.set_linewidth(1.2), sp.set_linestyle("--")
    for corner, dy in (((0, 1), -r), ((0, 0), r)):
        fig.add_artist(ConnectionPatch(xyA=(cx, cy + dy), coordsA=ax_a.transData, xyB=corner,
                                       coordsB=ax_in.transAxes, color="#333333", lw=0.8, ls="--"))
    ax_b.imshow(plt.imread(p["shell"]))
    row = load_row(sup, m)
    bonds = ", ".join("%.2f" % d for d in sorted(json_list(row.M_X_bare))[:3])
    el = "S" if sup == "MoS2" else "O"
    views = ["bare: %s–%s %s Å%s" % (m, el, bonds, ", %s %.1f Å above the S plane" % (m, row.M_height_bare)
                                     if sup == "MoS2" else ""),
             "with H on %s; ΔG$_H$ = %+.2f eV (xTB), %+.2f eV with offset" % (m, row.dG_H, row.dG_H_offset)]
    hc = h / 2
    for k, label in enumerate(views):
        top = y0 + h - k * hc
        fig.text(xs[2], top, label, fontsize=10.5, va="top", color=INK)
        ax = fig.add_axes([xs[2], top - hc, w, hc - 0.05])
        ax.imshow(ev.autocrop(plt.imread(p["site%d" % k])))
        ax.set_axis_off()
    ax_d = fig.add_axes([xs[3] + 0.035, y0 + 0.09, 0.15, h - 0.24])
    panel_chart(ax_d, sup, m)
    texts = [("a", ex["title"], ex["subtitle"] + "\n" + ex["chains"] + "; active centre as sticks"),
             ("b", ex["shell_title"], ex["shell_text"]),
             ("c", ex["motif_title"], ex["motif_text"]),
             ("d", ex["chart_title"], "Band: ±0.2 eV. The offset is the per-metal xTB error fitted on graphene "
                                      "against DFT (none for Mn, Mo, W). xTB only; spin-restricted.")]
    for x, (tag, title, note) in zip(xs, texts):
        fig.text(x, 0.86, tag, fontsize=17, fontweight="bold", va="bottom")
        fig.text(x + 0.016, 0.86, title, fontsize=13.5, va="bottom")
        wrapped = "\n".join(textwrap.fill(line, 50) for line in note.split("\n"))
        fig.text(x, 0.15, wrapped, fontsize=10.5, color="#4a4a4a", va="top")
    for x0, x1 in ((0.232, 0.258), (0.497, 0.523), (0.762, 0.788)):
        fig.add_artist(FancyArrowPatch((x0, 0.48), (x1, 0.48), transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=22, lw=2, color="#888888"))
    for a in (ax_a, ax_b):
        a.set_axis_off()
    out = os.path.join(HERE, "%s.png" % name)
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    print("wrote", out)


def json_list(s):
    return [float(x) for x in str(s).strip("[]").split(",") if x.strip()]


def load_row(support, metal):
    t = pd.read_csv(TABLE)
    r = t[(t.support == support) & (t.metal == metal) & (t.state == "H_M")]
    return r.iloc[0] if len(r) else None


if __name__ == "__main__":
    for n in sys.argv[1:] or EXAMPLES:
        figure(n)
