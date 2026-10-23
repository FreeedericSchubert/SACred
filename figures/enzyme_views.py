"""Enzyme -> first shell -> graphene site, as one figure per example (PyMOL renders, matplotlib layout).

    python figures/enzyme_views.py                 # every example in EXAMPLES -> figures/<name>.png
    python figures/enzyme_views.py nife            # one example

Panel a: the whole enzyme as cartoon (one colour per chain), the active centre as sticks with the metals as
spheres. The centre is buried, so it is rendered as its own layer on top of the cartoon (with a white halo)
and repeated as a zoomed inset: sticks of the centre inside the surrounding cartoon. Panel b: the first
coordination shell of the catalytic metal only (metal + donor atoms, as sticks), the motif the PDB sweep
extracts. Panel c: the graphene site(s) the motif maps to (data/enzyme_sites.csv), each in a given state (bare,
H on the metal, H on S): the relaxed DFT structure when dft/runs has it, else the starting geometry from
materials/graphene/sites.py; the label says which, and gives dG_H from data/dgh_graphene.csv. Metal-donor
distances that opened beyond 2.5 A during the relaxation are drawn dashed with their length.
"""
import json
import os
import sys
import textwrap

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pymol
import pandas as pd
from ase.io import read as ase_read, write as ase_write
from matplotlib.patches import ConnectionPatch, Ellipse, FancyArrowPatch
from scipy.ndimage import maximum_filter

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, os.path.join(ROOT, "materials", "graphene"))
from sites import build, with_h   # noqa: E402

PDB = os.path.join(ROOT, "data", "pdb")
RUNS = os.path.join(ROOT, "dft", "runs", "graphene")
DGH = os.path.join(ROOT, "data", "dgh_graphene.csv")
OPENED = 2.5            # A: a metal-donor distance beyond this is drawn dashed in panel c
SHELL_LABEL = "'%s %s%s' % (elem.capitalize(), resn.capitalize(), resi)"   # e.g. S Cys453
TMP = os.path.join(HERE, "tmp")
SIZE = (1400, 1200)     # px per rendered panel
CHAIN_COLOURS = ["lightblue", "wheat", "palegreen", "lightpink", "lightorange", "palecyan"]
CARBON = "gray55"       # carbon of sticks
HALO_PX = 7             # white rim around the active-centre layer in panel a

EXAMPLES = {
    "nife": dict(
        pdb="4U9H", title="[NiFe] hydrogenase", subtitle="PDB 4U9H, 0.89 Å · H$_2$ evolution",
        chains="large subunit blue, small subunit wheat",
        centre="resn EJ2",                                       # Ni, Fe, CO, 2 CN
        centre_residues="chain L and resi 81+84+546+549",        # the four Cys
        metal="chain L and resn EJ2 and name NI",
        donors="chain L and resi 81+84+546+549 and name SG",
        metal_bonds=[("chain L and resn EJ2 and name NI", "chain L and resi 81+84+546+549 and name SG"),
                     ("chain L and resn EJ2 and name FE", "chain L and resi 84+549 and name SG")],
        shell_title="First shell of Ni: Ni–S$_4$",
        shell_text="4 Cys, Ni–S 2.18–2.54 Å; τ$_4$ = 0.54, between square planar and tetrahedral. "
                   "Cys84 and Cys549 also bind the Fe.",
        motifs=[("Ni-S4:DV", "H_M", "double vacancy: Ni–S$_4$, all four S kept, flattened"),
                ("Ni-S3:SV", "H_M", "single vacancy: Ni–S$_3$, Cys549 left out, its place open for H")],
        motif_title="Graphene analogues",
        motif_text="Both are computed, with H on Ni (grey). The sweep's planarity test puts this site "
                   "(τ$_4$ = 0.54) on the single-vacancy side.",
    ),
    "laccase": dict(
        pdb="6RI6", title="Laccase", subtitle="PDB 6RI6, 0.93 Å · O$_2$ reduction",
        chains="Steccherinum murashkinskyi laccase 2, one chain",
        centre="resn CU",                                        # T1, T2 and the T3 pair
        centre_residues="chain A and resi 452+402+112+67+110+454+65+400+458+397+453",
        metal="chain A and resn CU and resi 504",               # T1 Cu
        donors="chain A and ((resi 458+397 and name ND1) or (resi 453 and name SG))",
        metal_bonds=[("chain A and resn CU and resi 504",
                      "chain A and ((resi 458+397 and name ND1) or (resi 453 and name SG))"),
                     ("chain A and resn CU and resi 501", "chain A and resi 452+402+112 and name NE2"),
                     ("chain A and resn CU and resi 502",
                      "chain A and ((resi 67 and name ND1) or (resi 110+454 and name NE2))"),
                     ("chain A and resn CU and resi 503", "chain A and resi 65+400 and name NE2")],
        shell_title="First shell of the T1 Cu: Cu–N$_2$S",
        shell_text="2 His (Nδ) and Cys; Cu–N 1.99 and 2.04 Å, Cu–S 2.19 Å. Three-coordinate: this fungal "
                   "laccase has no axial Met.",
        motifs=[("Cu-N2S:SV", "bare", "graphene single vacancy"),
                ("Cu-N2S:SV", "H_M", "with H on Cu")],
        motif_title="Graphene analogue: Cu–N$_2$S",
        sac_turn=-76,
        motif_text="PBE-D3, one spin state. When H binds, the Cu–S bond opens (dashed) and the S flips below "
                   "the sheet. The bare site prefers the S below as well, so ΔG$_H$ ends at +0.36 eV.",
    ),
}


def pm():
    pymol.finish_launching(["pymol", "-qc"])
    from pymol import cmd
    cmd.reinitialize()
    cmd.bg_color("white")
    cmd.set("fetch_path", TMP)                 # ligand definitions PyMOL downloads go to the scratch folder
    for k, v in {"ray_opaque_background": 1, "antialias": 2, "ray_trace_mode": 1, "ray_trace_gain": 0.06,
                 "ray_shadows": 0, "specular": 0.15, "ambient": 0.45, "direct": 0.6, "depth_cue": 0,
                 "orthoscopic": 1, "cartoon_fancy_helices": 1, "valence": 0, "label_font_id": 7,
                 "label_color": "black", "label_size": 24, "float_labels": 1}.items():
        cmd.set(k, v)
    return cmd


def clean(cmd, obj):
    cmd.remove(obj + " and hydro")
    cmd.remove(obj + " and not alt ''+A")      # keep one conformer
    cmd.alter(obj, "alt=''")
    cmd.remove(obj + " and solvent")


def colour_sticks(cmd, sel):
    cmd.color(CARBON, "(%s) and elem C" % sel)
    cmd.util.cnc(sel)                          # N, O, S, metals by element


def render(cmd, path, opaque=True, size=SIZE):
    cmd.set("ray_opaque_background", 1 if opaque else 0)
    cmd.png(path, width=size[0], height=size[1], dpi=300, ray=1)


def load_enzyme(cmd, ex):
    cmd.load(os.path.join(PDB, ex["pdb"] + ".cif.gz"), "enz")
    clean(cmd, "enz")
    for a, b in ex["metal_bonds"]:
        cmd.bond("enz and " + a, "enz and " + b)
    return "enz and ((%s) or ((%s) and sidechain))" % (ex["centre"], ex["centre_residues"])


def panel_enzyme(ex, paths):
    """paths: cartoon, active-centre layer (transparent background), inset cartoon, inset centre layer."""
    cmd = pm()
    site = load_enzyme(cmd, ex)
    cmd.hide("everything")
    for c, colour in zip(cmd.get_chains("enz and polymer.protein"), CHAIN_COLOURS):
        cmd.color(colour, "enz and polymer.protein and chain %s" % c)
    colour_sticks(cmd, site)
    cmd.orient("enz and polymer.protein")
    cmd.zoom("enz and polymer.protein", buffer=-2, complete=1)
    cmd.show("cartoon", "enz and polymer.protein")
    cmd.set("cartoon_transparency", 0.15)
    render(cmd, paths[0])
    cmd.hide("cartoon")
    cmd.show("sticks", site)
    cmd.show("spheres", "(%s) and metals" % site)
    cmd.set("stick_radius", 0.55)
    cmd.set("sphere_scale", 0.75)
    render(cmd, paths[1], opaque=False)
    # inset: the centre in its pocket; cartoon and centre again as two layers, the Cys now from CA on
    near = "enz and polymer.protein and byres (all within 12 of (%s))" % site
    full = "enz and ((%s) or ((%s) and (sidechain or name CA)))" % (ex["centre"], ex["centre_residues"])
    colour_sticks(cmd, full)
    cmd.hide("everything")
    cmd.show("cartoon", near)
    cmd.set("cartoon_transparency", 0.2)
    cmd.zoom(full, buffer=3.0, complete=1)
    cmd.clip("slab", 34, full)
    render(cmd, paths[2])
    cmd.hide("cartoon")
    cmd.show("sticks", full)
    cmd.show("spheres", "(%s) and metals" % full)
    cmd.set("stick_radius", 0.28)
    cmd.set("sphere_scale", 0.4)
    render(cmd, paths[3], opaque=False)


def panel_shell(ex, path):
    cmd = pm()
    load_enzyme(cmd, ex)
    cmd.create("shell", "enz and ((%s) or (%s))" % (ex["metal"], ex["donors"]))
    cmd.delete("enz")
    cmd.unbond("shell", "shell")
    cmd.bond("shell and metals", "shell and not metals")
    cmd.hide("everything")
    cmd.show("sticks", "shell")
    cmd.set("stick_radius", 0.17)
    colour_sticks(cmd, "shell")
    cmd.label("shell and metals", "elem.capitalize()")
    cmd.label("shell and not metals", SHELL_LABEL)
    cmd.set("label_position", (0, 0, 4))
    cmd.orient("shell")
    cmd.zoom("shell", buffer=2.4, complete=1)
    render(cmd, path)


def site_structure(motif, state):
    """(atoms, metal index, donor indices, relaxed?) for a graphene site in a state (bare, H_M, H_S)."""
    run = os.path.join(RUNS, motif.replace(":", "_").replace("=", ""), state)
    out = os.path.join(run, "relax.out")
    if os.path.exists(out) and "End final coordinates" in open(out).read():
        meta = json.load(open(os.path.join(run, "meta.json")))
        return ase_read(out, format="espresso-out", index=-1), meta["imetal"], meta["idonors"], True
    atoms, info = build(motif)
    if state != "bare":
        atoms = with_h(atoms, info, state[2:])
    return atoms, info["imetal"], info["idonors"], False


def panel_graphene(motif, state, path, turn=-62):
    atoms, im, idon, relaxed = site_structure(motif, state)
    centre = atoms.positions[im]
    d = np.linalg.norm((atoms.positions - centre)[:, :2], axis=1)
    keep = [i for i in range(len(atoms)) if d[i] < 6.2]
    for _ in range(3):                                     # drop edge carbons with a single neighbour left
        sub = atoms[keep]
        nn = (sub.get_all_distances(mic=False) < 1.6).sum(axis=1) - 1   # no periodic images: a patch has edges
        keep = [k for k, n, a in zip(keep, nn, sub) if n >= 2 or a.symbol != "C"]
    patch = atoms[keep]
    patch.pbc = False                                      # PyMOL must not see a cell
    pdb = os.path.join(TMP, "graphene_site.pdb")
    ase_write(pdb, patch, format="proteindatabank")
    cmd = pm()
    cmd.load(pdb, "sac")
    cmd.unbond("sac and (metals or elem H)", "sac")
    cmd.bond("sac and metals", "sac and not elem C+H within %.1f of (sac and metals)" % OPENED)
    cmd.bond("sac and metals", "sac and elem C within 2.2 of (sac and metals)")
    cmd.bond("sac and elem H", "sac and not elem H within 1.7 of (sac and elem H)")
    cmd.hide("everything")
    cmd.show("sticks", "sac and not elem H")
    cmd.set("stick_radius", 0.17)
    cmd.show("spheres", "sac and elem H")
    cmd.set("sphere_scale", 0.3, "sac and elem H")
    colour_sticks(cmd, "sac")
    cmd.color("gray85", "sac and elem H")
    cmd.label("sac and metals", "elem.capitalize()")
    cmd.label("sac and elem H", "'H'")
    cmd.set("label_position", (0, 0, 4))
    opened = []                                            # bonds that opened: dashed; lengths go into the label
    for k in idon:
        dist = atoms.get_distance(im, k, mic=True)
        if k in keep and dist > OPENED:
            cmd.distance("open%d" % k, "sac and index %d" % (keep.index(im) + 1), "sac and index %d" % (keep.index(k) + 1))
            cmd.hide("labels", "open%d" % k)
            opened.append("%s–%s %.2f Å" % (atoms[im].symbol, atoms[k].symbol, dist))
    cmd.set("dash_color", "black")
    cmd.set("dash_width", 3.0)
    cmd.set("dash_radius", 0.05)
    cmd.set("label_size", 22)
    cmd.orient("sac and not elem H")
    cmd.turn("x", turn)                      # oblique view: the metal sits above the sheet
    cmd.zoom("sac", buffer=0.5, complete=1)
    cmd.clip("slab", 80)                      # nothing cut by the near plane; margins are cropped later
    render(cmd, path)
    return relaxed, opened


def dg_h(motif, state):
    """(dG_H, provisional?) from data/dgh_graphene.csv; provisional while another bare configuration is pending."""
    if state.startswith("bare") or not os.path.exists(DGH):
        return None, False
    t = pd.read_csv(DGH, keep_default_na=False, na_values=[""])
    r = t[(t.motif == motif) & (t.state == state)].dropna(subset=["dG_H"])
    if not len(r):
        return None, False
    pending = r.get("bare_states_pending", pd.Series([""])).iloc[0]
    return float(r.dG_H.iloc[0]), bool(pd.notna(pending) and str(pending).strip())


def autocrop(img, margin=12):
    """Trim the white margin of a render (the graphene views are wide and flat)."""
    ys, xs = np.nonzero((img[..., :3] < 0.97).any(axis=-1))
    return img[max(ys.min() - margin, 0):ys.max() + margin, max(xs.min() - margin, 0):xs.max() + margin]


def halo(layer):
    """Active-centre layer with a white rim, so it reads on top of the cartoon."""
    alpha = layer[..., 3]
    rim = maximum_filter(alpha, size=2 * HALO_PX + 1)
    out = np.ones_like(layer)
    out[..., 3] = rim
    a = alpha[..., None]
    out[..., :3] = layer[..., :3] * a + (1 - a)
    return out


def figure(name):
    ex = EXAMPLES[name]
    os.makedirs(TMP, exist_ok=True)
    p = {k: os.path.join(TMP, "%s_%s.png" % (name, k))
         for k in ("cartoon", "centre", "inset", "inset_centre", "shell", "sac0", "sac1")}
    panel_enzyme(ex, [p["cartoon"], p["centre"], p["inset"], p["inset_centre"]])
    panel_shell(ex, p["shell"])
    views = [panel_graphene(motif, state, p["sac%d" % k], ex.get("sac_turn", -62))
             for k, (motif, state, _) in enumerate(ex["motifs"])]

    fig = plt.figure(figsize=(16, 6.4))
    w, h, y0, xs = 0.29, 0.62, 0.17, (0.0, 0.355, 0.71)       # panel boxes match the 1400 x 1200 renders
    ax_a, ax_b, ax_c = (fig.add_axes([x, y0, w, h]) for x in xs)
    ax_c.set_axis_off()
    hc = h / len(ex["motifs"])                                 # graphene views stacked in panel c
    layer = plt.imread(p["centre"])
    ax_a.imshow(plt.imread(p["cartoon"]))
    ax_a.imshow(halo(layer))
    ys, xs_ = np.nonzero(layer[..., 3] > 0.1)                  # where the centre sits in panel a
    cx, cy, r = xs_.mean(), ys.mean(), 0.6 * max(np.ptp(xs_), np.ptp(ys)) + 25
    ax_a.add_patch(Ellipse((cx, cy), 2 * r, 2 * r, fill=False, ls="--", lw=1.2, color="#333333"))
    ax_in = fig.add_axes([0.185, 0.52, 0.135, 0.33])
    ax_in.imshow(plt.imread(p["inset"]))
    ax_in.imshow(halo(plt.imread(p["inset_centre"])))
    ax_in.set_xticks([]), ax_in.set_yticks([])
    for sp in ax_in.spines.values():
        sp.set_color("#333333"), sp.set_linewidth(1.2), sp.set_linestyle("--")
    for corner, dy in (((0, 1), -r), ((0, 0), r)):
        fig.add_artist(ConnectionPatch(xyA=(cx, cy + dy), coordsA=ax_a.transData, xyB=corner,
                                       coordsB=ax_in.transAxes, color="#333333", lw=0.8, ls="--"))
    ax_b.imshow(plt.imread(p["shell"]))
    for k, (motif, state, label) in enumerate(ex["motifs"]):
        top = y0 + h - k * hc                                  # label line, then the view below it
        dg, provisional = dg_h(motif, state)
        relaxed, opened = views[k]
        label += (", relaxed (DFT)" if relaxed else ", starting geometry (DFT running)") + \
            ("; ΔG$_H$ = %+.2f eV%s" % (dg, " (provisional)" if provisional else "") if dg is not None else "") + \
            ("; %s (dashed)" % ", ".join(opened) if opened else "")
        fig.text(xs[2], top, label, fontsize=10.5, va="top", color="#222222")
        ax = fig.add_axes([xs[2], top - hc, w, hc - 0.05])
        ax.imshow(autocrop(plt.imread(p["sac%d" % k])))
        ax.set_axis_off()
    texts = [("a", ex["title"], ex["subtitle"] + "\n" + ex["chains"] + "; active centre as sticks"),
             ("b", ex["shell_title"], ex["shell_text"]),
             ("c", ex["motif_title"], ex["motif_text"])]
    for ax, x, (tag, title, note) in zip((ax_a, ax_b, ax_c), xs, texts):
        ax.set_axis_off()
        fig.text(x, 0.86, tag, fontsize=17, fontweight="bold", va="bottom")
        fig.text(x + 0.018, 0.86, title, fontsize=13.5, va="bottom")
        wrapped = "\n".join(textwrap.fill(line, 58) for line in note.split("\n"))
        fig.text(x, 0.15, wrapped, fontsize=10.5, color="#4a4a4a", va="top")
    for x0, x1 in ((0.300, 0.345), (0.655, 0.700)):              # reduction steps a -> b -> c
        fig.add_artist(FancyArrowPatch((x0, 0.48), (x1, 0.48), transform=fig.transFigure,
                                       arrowstyle="-|>", mutation_scale=22, lw=2, color="#888888"))
    out = os.path.join(HERE, "%s.png" % name)
    fig.savefig(out, dpi=200, bbox_inches="tight", facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    for n in sys.argv[1:] or EXAMPLES:
        figure(n)
