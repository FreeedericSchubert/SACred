"""Activity-stability figure of the one-page summary: dG_H against E_f for every DFT site with both values."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "..", "data")
d = pd.read_csv(os.path.join(DATA, "dgh_graphene.csv")).dropna(subset=["dG_H"])
d = d.assign(a=d.dG_H.abs()).sort_values("a").drop_duplicates("motif")[["motif", "dG_H"]]
s = pd.read_csv(os.path.join(DATA, "stability_graphene.csv")).dropna(subset=["E_f"])[["motif", "E_f", "U_diss"]]
m = d.merge(s, on="motif")
good = (m.dG_H.abs() <= 0.2) & (m.E_f < 0) & (m.U_diss > 0)
sub = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
label = lambda x: x.split(":")[0].replace("-", "–", 1).translate(sub)     # noqa: E731
YMAX = 5.0
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.5})
fig, ax = plt.subplots(figsize=(3.6, 2.45), dpi=300)
ax.add_patch(plt.Rectangle((-0.2, -3), 0.4, 3, facecolor="#e6e9eb", edgecolor="#9aa5ab", lw=0.6, zorder=0))
ax.axhline(0, color="#9aa5ab", lw=0.6, ls="--", zorder=1)
ax.axvline(0, color="#9aa5ab", lw=0.6, ls="--", zorder=1)
y = m.E_f.clip(upper=YMAX)
ax.scatter(m.dG_H[~good], y[~good], s=10, color="#9aa5ab", zorder=2, lw=0)
ax.scatter(m.dG_H[good], y[good], s=22, color="#2c6bd3", zorder=3, lw=0)
for r in m[good].itertuples():
    ax.annotate(label(r.motif), (r.dG_H, r.E_f), xytext=(0.32, r.E_f), va="center", fontsize=7.5, fontweight="bold",
                color="#1c2427", arrowprops=dict(arrowstyle="-", color="#9aa5ab", lw=0.5, shrinkA=0, shrinkB=3))
ax.text(-0.22, -2.85, "target", ha="right", va="bottom", fontsize=6.5, color="#59666b")
n_off = int((m.E_f > YMAX).sum())
if n_off:
    ax.text(1.55, YMAX - 0.15, "%d sites above 5 eV at the top edge" % n_off, ha="right", va="top", fontsize=6, color="#59666b")
ax.set_xlim(-1.2, 1.6)
ax.set_ylim(-3, YMAX + 0.2)
ax.set_xlabel(r"$\Delta G_\mathrm{H}$ (eV), ideal 0")
ax.set_ylabel(r"$E_\mathrm{f}$ vs bulk metal (eV)")
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
ax.tick_params(length=2)
fig.tight_layout(pad=0.3)
fig.savefig(os.path.join(HERE, "figure.png"))
print("%d sites, %d near-optimal (|dG_H| <= 0.2), %d pass; E_f > 0: %d" % (
    len(m), int((m.dG_H.abs() <= 0.2).sum()), int(good.sum()), int((m.E_f > 0).sum())))
