#!/usr/bin/env python3
"""
make_thesis_figures.py  (v3 - elegant, grey-free, no label collisions)
----------------------------------------------------------------------
Nine figures, each with a DISTINCT jewel/earth palette (no grey as a data colour,
no palette reused across figures), warm ivory backgrounds, and every label given
its own airspace so nothing overlaps. All numbers hard-coded from the real
Phase 1-5 results.

Run:  python3 make_thesis_figures.py   ->   ./figures/*.png  (300 dpi)
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import FancyBboxPatch, Patch

os.makedirs("figures", exist_ok=True)
mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 11.5,
    "axes.titlesize": 14.5, "axes.titleweight": "bold",
    "axes.labelsize": 11.5, "axes.labelweight": "medium",
    "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False,
    "axes.grid": True, "grid.alpha": 0.22, "grid.linewidth": 0.6, "axes.axisbelow": True,
})

# Distinct palette per figure. No grey as a data colour anywhere.
PAL = {
 "f1": dict(a="#2D6A4F", b="#D98324", bg="#FBF8F3", ink="#1B2D2A", tie="#8A9B8F"),         # forest vs amber
 "f2": dict(seq=["#264653", "#2A9D8F", "#E9C46A"], bg="#FBF7EF", ink="#22333B"),            # teal->sand funnel
 "f3": dict(a="#5A189A", b="#9D4EDD", c="#C77DFF", bg="#FAF6FD", ink="#2A0A45",
            ok="#3A86FF", hard="#C1121F"),                                                  # violet ladder
 "f4": dict(cmap="viridis", hot="#E5383B", cool="#0077B6", bg="#F5FAF8", ink="#1B3A4B"),    # jewel scatter
 "f5": dict(b1="#457B9D", b2="#E07A5F", b3="#2A9D8F", bg="#FBF7F2", ink="#1D3557",
            sig="#6A4C93", nsig="#B9AE96"),                                                 # slate/terracotta/teal
 "f6": dict(conf="#2A9D8F", weak="#E9C46A", inact="#E07A5F", unt="#457B9D",
            bar="#5A189A", bg="#F4FAF8", ink="#12403B"),                                     # teal/gold/coral/blue
 "f7": dict(fam="#2A9D8F", miss="#C1121F", tot="#EADFC8", bg="#FBF8F2", ink="#22403B"),      # cream track
 "f8": dict(bgb="#7C9C99", hit="#6A4C93", null="#E5383B", bg="#F7F5FC", ink="#2C1A47"),      # sage/plum/red
 "f9": dict(within="#1D7874", cross="#EE964B", yes="#2A9D8F", no="#D0393B",
            bg="#F4F9F8", ink="#0B2E2B"),                                                    # pine/apricot
}

def card(ax, x, y, text, color, ha="left", va="center", fs=9.5, pad=0.42, alpha=0.97):
    ax.annotate(text, (x, y), ha=ha, va=va, fontsize=fs, color=color, zorder=20,
                fontweight="bold", linespacing=1.32,
                bbox=dict(boxstyle=f"round,pad={pad}", fc="white", ec=color, lw=1.3, alpha=alpha))

def leader(ax, x0, y0, x1, y1, color):
    ax.annotate("", xy=(x1, y1), xytext=(x0, y0), zorder=19,
                arrowprops=dict(arrowstyle="-", color=color, lw=1.8,
                                connectionstyle="arc3,rad=0.18"))
    ax.scatter([x1], [y1], s=44, color=color, zorder=21, edgecolor="white", lw=1.2)

def tint(fig, ax, color):
    fig.patch.set_facecolor("white"); ax.set_facecolor(color)

def save(fig, name):
    p = f"figures/{name}.png"; fig.savefig(p, facecolor="white"); plt.close(fig); print("wrote", p)

# ---------------------------------------------------------------- FIG 1 (forest/amber)
def fig1():
    P = PAL["f1"]; fig, ax = plt.subplots(figsize=(9, 5.5)); tint(fig, ax, P["bg"])
    groups = ["n = 33\nstrict", "n = 33\nloose", "n = 23\nstrict", "n = 23\nloose"]
    holo = [87.9, 87.9, 82.6, 82.6]; apo = [75.8, 78.8, 69.6, 73.9]
    x = np.arange(len(groups)); w = 0.36
    b1 = ax.bar(x-w/2, holo, w, label="Holo (bound)", color=P["a"], edgecolor="white", lw=2, zorder=3)
    b2 = ax.bar(x+w/2, apo, w, label="Apo (unbound)", color=P["b"], edgecolor="white", lw=2, zorder=3)
    for b in list(b1)+list(b2):
        ax.text(b.get_x()+b.get_width()/2, b.get_height()+1.4, f"{b.get_height():.1f}",
                ha="center", va="bottom", fontsize=10, fontweight="bold", color=P["ink"])
    for i in range(len(groups)):
        ax.plot([i-w/2, i-w/2, i+w/2, i+w/2],
                [holo[i]+6, holo[i]+8.5, holo[i]+8.5, apo[i]+6],
                color=P["tie"], lw=1.1, zorder=2)
        card(ax, i, holo[i]+10.2, f"\u2212{holo[i]-apo[i]:.1f} pts", P["ink"], ha="center", fs=9, pad=0.30)
    ax.set_xticks(x); ax.set_xticklabels(groups, fontsize=10.5)
    ax.set_ylabel("Top-3 detection rate  (%)"); ax.set_ylim(0, 112); ax.set_yticks(range(0, 101, 20))
    ax.set_title("Pocket detection collapses in the unbound state", color=P["ink"], pad=34)
    ax.text(0.5, 1.045, "holo beats apo in every panel  \u2022  zero apo-only discordant targets",
            transform=ax.transAxes, ha="center", fontsize=10, style="italic", color=P["b"])
    ax.legend(loc="lower center", ncol=2, frameon=True, fontsize=10.5, bbox_to_anchor=(0.5, -0.22))
    ax.grid(axis="x", alpha=0)
    save(fig, "fig1_detection_asymmetry")

# ---------------------------------------------------------------- FIG 2 (teal->sand funnel)
def fig2():
    P = PAL["f2"]; fig, ax = plt.subplots(figsize=(10, 5.4)); tint(fig, ax, P["bg"])
    stages = ["Raw geometric\nmatches", "Pocket-supported\n(P2Rank cavity)", "High-confidence\n(tiered)"]
    vals = [1107, 681, 58]; cols = P["seq"]; maxv = 1107; ys = [2, 1, 0]
    for i in range(len(vals)-1):
        l_hi=(maxv-vals[i])/2; r_hi=l_hi+vals[i]; l_lo=(maxv-vals[i+1])/2; r_lo=l_lo+vals[i+1]
        ax.add_patch(plt.Polygon([[l_hi, ys[i]-0.33],[r_hi, ys[i]-0.33],
                                  [r_lo, ys[i+1]+0.33],[l_lo, ys[i+1]+0.33]],
                                 closed=True, color=cols[i+1], alpha=0.22, zorder=1))
    for yi, v, c, s in zip(ys, vals, cols, stages):
        left=(maxv-v)/2
        ax.barh(yi, v, left=left, color=c, edgecolor="white", lw=3, height=0.66, zorder=3)
        ax.text(maxv/2, yi+0.44, s, ha="center", va="bottom", color=c,           # NAME above bar
                fontweight="bold", fontsize=11, zorder=4, linespacing=1.1)
        ax.text(maxv/2, yi, f"{v:,}", ha="center", va="center", color="white",   # COUNT inside bar
                fontweight="bold", fontsize=17, zorder=5)
        pct="100% of raw" if v==maxv else f"{v/maxv*100:.1f}% of raw"
        card(ax, maxv+55, yi, pct, c, ha="left", fs=10.5, pad=0.4)
    card(ax, maxv/2, 1.5, "\u00d7 0.62 pass", cols[1], ha="center", fs=8.5, pad=0.26)
    card(ax, maxv/2, 0.5, "\u00d7 0.085 pass", cols[2], ha="center", fs=8.5, pad=0.26)
    ax.set_yticks([]); ax.set_xlim(-30, maxv+360); ax.set_ylim(-0.7, 2.95); ax.grid(False)
    ax.set_xlabel("Number of hits")
    ax.set_title("Off-target funnel \u2014 a ~20-fold concentration to 58 candidates", color=P["ink"], pad=14)
    save(fig, "fig2_funnel")

# ---------------------------------------------------------------- FIG 3 (violet ladder)
def fig3():
    P = PAL["f3"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.0), gridspec_kw={"width_ratios": [1, 1.35]})
    tint(fig, ax1, P["bg"]); tint(fig, ax2, P["bg"])
    cats = ["Rank 1", "Top 2\u20135", "Outside\ntop 5"]; vals = [16, 6, 1]
    cols = [P["a"], P["b"], P["c"]]
    bars = ax1.bar(cats, vals, color=cols, edgecolor="white", lw=2, zorder=3, width=0.62)
    for b, v in zip(bars, vals):
        ax1.text(b.get_x()+b.get_width()/2, v+0.3, str(v), ha="center", fontweight="bold",
                 fontsize=12, color=P["ink"])
    ax1.set_ylabel("Targets (of 23)"); ax1.set_ylim(0, 19)
    ax1.set_title("Self-recovery: all 23 recover", color=P["ink"], pad=26)
    ax1.text(0.5, 1.04, "median top-hit RMSD 0.48 \u00c5  \u2022  18/23 sub-\u00c5ngstr\u00f6m",
             transform=ax1.transAxes, ha="center", fontsize=9.5, style="italic", color=P["b"])
    ax1.grid(axis="x", alpha=0)
    tg = ["MEK7", "GR", "MEK1", "USP15", "PTGR1", "USP7", "LTA4H"]
    dca = [15.71, 8.46, 4.43, 4.59, 2.12, 25.08, 6.10]; fdrank = [1, 1, 2, 1, 1, 1, 1]
    yy = np.arange(len(tg))[::-1]
    ax2.hlines(yy, 0, dca, color=P["c"], lw=3, zorder=2)
    ax2.scatter(dca, yy, s=170, color=P["a"], edgecolor="white", lw=2, zorder=3)
    ax2.axvline(4.0, ls=(0, (4, 3)), color=P["hard"], lw=1.5, zorder=1)
    card(ax2, 4.0, 6.8, "4 \u00c5 detection\ncutoff", P["hard"], ha="center", fs=8.5, pad=0.3)
    for yi, d, r in zip(yy, dca, fdrank):
        ax2.text(d+0.7, yi, f"DCA {d:.1f} \u00c5", va="center", fontsize=9, color=P["a"], fontweight="bold")
        ax2.scatter(-1.8, yi, s=340, marker="o", color=P["ok"], zorder=3, edgecolor="white", lw=1.5)
        ax2.text(-1.8, yi, f"rank {r}", ha="center", va="center", fontsize=7.5,
                 color="white", fontweight="bold", zorder=4)
    ax2.set_yticks(yy); ax2.set_yticklabels(tg, fontsize=11, fontweight="bold")
    ax2.set_xlim(-3.4, 30); ax2.set_xlabel("P2Rank DCA (\u00c5) \u2014 higher = worse detection")
    ax2.grid(axis="y", alpha=0)
    ax2.set_title("The seven sites P2Rank misses are all\nrecovered by FoldDisco (blue circle = self-rank)", color=P["ink"], pad=12)
    save(fig, "fig3_self_recovery")

# ---------------------------------------------------------------- FIG 4 (viridis, callouts parked)
def fig4():
    P = PAL["f4"]; fig, ax = plt.subplots(figsize=(9.5, 6.8)); tint(fig, ax, P["bg"])
    tg = ["DHFR","UCK2","GR","PDE5A","MEK1","SYK","PDE4A","AKR1B1","USB1","RORC","USP15",
          "USP28","ALK","PTGR1","PXR","CA2","USP7","MEK7","FKBP1A","PIN1","AR","LTA4H","NUDT1"]
    ov = [1.00,0.65,0.72,0.94,0.44,0.92,1.00,0.90,1.00,0.50,0.83,0.93,1.00,0.88,0.92,1.00,
          1.00,1.00,0.79,1.00,1.00,1.00,0.82]
    bb = [0.99,2.08,2.52,1.24,2.20,0.56,0.33,1.35,0.73,0.94,1.14,1.50,0.50,0.58,0.97,0.32,
          2.63,1.82,0.87,0.53,0.70,0.15,1.02]
    sc = ax.scatter(ov, bb, s=135, c=bb, cmap=P["cmap"], edgecolor="white", lw=1.5, zorder=4)
    special = {"MEK7", "MEK1"}
    # targets clustered at overlap=1.00 need staggered offsets to avoid overlap
    offs = {"DHFR":(6,4),"PDE4A":(6,-9),"USB1":(-30,4),"ALK":(6,3),"CA2":(-24,-9),
            "USP7":(8,2),"PIN1":(6,4),"AR":(6,-9),"LTA4H":(8,2),"AKR1B1":(6,3),
            "SYK":(6,3),"PDE5A":(6,3),"GR":(7,2),"RORC":(6,3),"USP15":(6,3),
            "USP28":(6,3),"PXR":(6,-9),"UCK2":(6,3),"FKBP1A":(-40,2),"NUDT1":(6,3)}
    for x, y, t in zip(ov, bb, tg):
        if t in special: continue
        dx, dy = offs.get(t, (6, 3))
        ax.annotate(t, (x, y), fontsize=8.3, color=P["ink"], xytext=(dx, dy),
                    textcoords="offset points", zorder=5)
    ax.scatter([1.00], [1.82], s=440, facecolors="none", edgecolors=P["hot"], lw=3, zorder=5)
    ax.scatter([0.44], [2.20], s=440, facecolors="none", edgecolors=P["cool"], lw=3, zorder=5)
    card(ax, 0.60, 0.60,
         "MEK7 \u2014 overlap 1.00,\nyet cavity undetectable\n(DCA 15.7 \u00c5)\ngeometry \u2260 detectability",
         P["hot"], ha="left", va="center", fs=9)
    leader(ax, 0.71, 0.86, 1.00, 1.82, P["hot"])
    card(ax, 0.585, 3.05,
         "MEK1 \u2014 all three axes\nagree (true collapse):\noverlap 0.44, bb 2.20 \u00c5",
         P["cool"], ha="left", va="center", fs=9)
    leader(ax, 0.55, 2.85, 0.44, 2.22, P["cool"])
    ax.set_xlabel("Apo spatial overlap  (binding-residue persistence)")
    ax.set_ylabel("Backbone RMSD  (\u00c5, global motion)")
    ax.set_xlim(1.06, 0.38); ax.set_ylim(0.0, 3.5)
    ax.set_title("Three axes of conformational change disagree \u2014 informatively", color=P["ink"], pad=14)
    cb = plt.colorbar(sc, pad=0.02); cb.set_label("Backbone RMSD (\u00c5)")
    save(fig, "fig4_conformational_axes")

# ---------------------------------------------------------------- FIG 5 (slate/terracotta/teal + plum)
def fig5():
    P = PAL["f5"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.0), gridspec_kw={"width_ratios": [1.35, 1]})
    tint(fig, ax1, P["bg"]); tint(fig, ax2, P["bg"])
    tg = ["NUDT1","DHFR","USB1","USP7","SYK","ALK","PDE4A","PTGR1","AKR1B1","USP15","CA2"]
    jac = [0.000,0.011,0.029,0.035,0.050,0.070,0.090,0.120,0.150,0.250,0.318]
    cols = [P["b1"] if j < 0.1 else P["b2"] if j < 0.2 else P["b3"] for j in jac]
    ax1.bar(range(len(tg)), jac, color=cols, edgecolor="white", lw=1.6, zorder=3, width=0.7)
    ax1.axhline(0.066, ls=(0, (5, 3)), color=P["ink"], lw=1.6, zorder=2)
    card(ax1, len(tg)-1.2, 0.115, "median 0.066", P["ink"], ha="right", fs=9.5, pad=0.32)
    ax1.set_xticks(range(len(tg))); ax1.set_xticklabels(tg, rotation=42, ha="right", fontsize=9)
    ax1.set_ylabel("Holo\u2013apo Jaccard overlap"); ax1.set_ylim(0, 0.36)
    ax1.set_title("Off-target lists diverge between conformations", color=P["ink"], pad=24)
    ax1.text(0.5, 1.03, "reliable subset, n = 11  \u2022  colour = divergence band",
             transform=ax1.transAxes, ha="center", fontsize=9.5, style="italic", color=P["b1"])
    ax1.grid(axis="x", alpha=0)
    leg = [Patch(color=P["b1"], label="Jaccard < 0.10"), Patch(color=P["b2"], label="0.10\u20130.20"),
           Patch(color=P["b3"], label="> 0.20")]
    ax1.legend(handles=leg, loc="upper left", frameon=True, fontsize=8.5)
    axes = ["spatial\noverlap", "backbone\nRMSD", "pocket\nRMSD"]; rho = [0.094, -0.482, -0.251]; pv = [0.782, 0.139, 0.457]
    bcol = [P["sig"] if p < 0.05 else P["nsig"] for p in pv]; yy = np.arange(3)[::-1]
    ax2.barh(yy, rho, color=bcol, edgecolor="white", lw=1.6, zorder=3, height=0.6)
    ax2.axvline(0, color=P["ink"], lw=1.2)
    for yi, r, p in zip(yy, rho, pv):
        ax2.text(r+(0.04 if r >= 0 else -0.04), yi, f"\u03c1={r:+.2f}\np={p:.2f}", va="center",
                 ha="left" if r >= 0 else "right", fontsize=9, fontweight="bold", color=P["ink"])
    ax2.set_yticks(yy); ax2.set_yticklabels(axes, fontsize=10)
    ax2.set_xlim(-0.85, 0.55); ax2.set_xlabel("Spearman \u03c1  (Jaccard vs axis)")
    ax2.grid(axis="y", alpha=0)
    ax2.set_title("Divergence is not explained by\npocket movement (all correlations n.s.)", color=P["ink"], pad=12)
    ax2.text(0.5, -0.22, "faded bars = not significant", transform=ax2.transAxes,
             ha="center", fontsize=9, style="italic", color=P["nsig"])
    save(fig, "fig5_dual_query_divergence")

# ---------------------------------------------------------------- FIG 6 (donut label fixed)
def fig6():
    P = PAL["f6"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.4), gridspec_kw={"width_ratios": [1, 1.25]})
    tint(fig, ax1, P["bg"]); tint(fig, ax2, P["bg"])
    sizes = [7, 4, 12, 35]
    labels = ["Confirmed \u22641 \u00b5M", "Weak 1\u201310 \u00b5M", "Tested-inactive", "Untested"]
    cols = [P["conf"], P["weak"], P["inact"], P["unt"]]
    wedges, _ = ax1.pie(sizes, colors=cols, startangle=90, counterclock=False, radius=1.0,
                        wedgeprops=dict(width=0.38, edgecolor="white", lw=2.5))
    # center label: ONLY "58 / hits", nothing else near it
    ax1.text(0, 0, "58\nhits", ha="center", va="center", fontsize=18, fontweight="bold", color=P["ink"])
    # legend well BELOW the donut so it never overlaps the ring or the centre
    ax1.legend(wedges, [f"{l}  ({s})" for l, s in zip(labels, sizes)],
               loc="upper center", bbox_to_anchor=(0.5, -0.08), ncol=2, frameon=True,
               fontsize=9.5, columnspacing=1.2, handlelength=1.2)
    ax1.set_title("ChEMBL binding status", color=P["ink"], pad=18)
    ax1.text(0.5, 1.0, "7 sub-\u00b5M of 23 measured = 30%", transform=ax1.transAxes,
             ha="center", fontsize=9.5, style="italic", color=P["conf"])
    names = ["MEK2  (selumetinib)", "MAP3K2  (crizotinib)", "AURKA  (crizotinib)",
             "YES1  (crizotinib)", "AKR1B10  (sulindac)", "FRK  (crizotinib)"]
    pot = [52, 72, 89.2, 100, 350, 1000]; yy = np.arange(len(names))[::-1]
    ax2.hlines(yy, 20, pot, color="#BEEAD9", lw=3, zorder=2)
    ax2.scatter(pot, yy, s=155, color=P["bar"], edgecolor="white", lw=2, zorder=3)
    ax2.set_xscale("log")
    for yi, p in zip(yy, pot):
        ax2.text(p*1.18, yi, f"{p:g} nM", va="center", fontsize=9.5, fontweight="bold", color=P["bar"])
    ax2.axvline(1000, ls=(0, (4, 3)), color=P["inact"], lw=1.5)
    card(ax2, 1000, 5.7, "1 \u00b5M", P["inact"], ha="center", fs=8.5, pad=0.3)
    ax2.set_yticks(yy); ax2.set_yticklabels(names, fontsize=9.5)
    ax2.set_xlabel("Measured potency  (nM, log scale)"); ax2.set_xlim(30, 2600); ax2.grid(axis="y", alpha=0)
    ax2.set_title("Seven confirmed sub-micromolar off-targets", color=P["ink"], pad=18)
    ax2.text(1.0, 1.0, "incl. documented crizotinib targets AURKA + YES1",
             transform=ax2.transAxes, ha="right", fontsize=9, style="italic", color=P["bar"])
    save(fig, "fig6_chembl_binding")

# ---------------------------------------------------------------- FIG 7 (cream track, no grey)
def fig7():
    P = PAL["f7"]; fig, ax = plt.subplots(figsize=(10, 5.6)); tint(fig, ax, P["bg"])
    drugs = ["Pentoxifylline", "Selumetinib", "Crizotinib", "Sulindac", "Tacrolimus",
             "Imatinib", "Ibrutinib", "Raloxifene", "Vorinostat"]
    found = [3, 1, 2, 1, 1, 0, 0, 0, 0]; total = [3, 1, 3, 3, 3, 4, 4, 2, 3]; fam = [True]*5+[False]*4
    yy = np.arange(len(drugs))[::-1]
    for yi, f, t, isfam in zip(yy, found, total, fam):
        ax.barh(yi, t, color=P["tot"], edgecolor="white", lw=1.5, height=0.66, zorder=2)
        ax.barh(yi, f, color=P["fam"] if isfam else P["miss"], edgecolor="white", lw=1.5, height=0.66, zorder=3)
        ax.text(t+0.09, yi, f"{f}/{t}", va="center", fontsize=10.5, fontweight="bold", color=P["ink"])
    ax.axhline(3.5, color=P["ink"], lw=1, ls=(0, (3, 3)), alpha=0.5)
    card(ax, 4.5, 6.5, "family off-targets\n\u2192 recovered", P["fam"], ha="center", fs=9, pad=0.32)
    card(ax, 4.5, 1.5, "unrelated folds\n\u2192 missed", P["miss"], ha="center", fs=9, pad=0.32)
    ax.set_yticks(yy); ax.set_yticklabels(drugs, fontsize=10.5)
    ax.set_xlabel("Documented off-targets recovered @ high confidence"); ax.set_xlim(0, 5.7); ax.grid(axis="y", alpha=0)
    ax.set_title("High-confidence recall splits by fold-relatedness  (8/26 = 0.31)", color=P["ink"], pad=26)
    leg = [Patch(color=P["fam"], label="Family off-target (found)"),
           Patch(color=P["miss"], label="Unrelated-fold (missed)"),
           Patch(color=P["tot"], label="Total documented")]
    ax.legend(handles=leg, loc="lower right", frameon=True, fontsize=9.5)
    ax.text(0.5, 1.045, "the fold-unrelated boundary is method-independent (BLAST and Foldseek also recover 0)",
            transform=ax.transAxes, ha="center", fontsize=9, style="italic", color=P["miss"])
    save(fig, "fig7_recall_split")

# ---------------------------------------------------------------- FIG 8 (sage/plum/red)
def fig8():
    P = PAL["f8"]; fig, ax = plt.subplots(figsize=(9, 5.6)); tint(fig, ax, P["bg"])
    cats = ["Background\n(pocket-bearing\nlibrary)", "All high\nhits", "Family-\nremoved"]
    rate = [5.3, 27.1, 0.0]; cols = [P["bgb"], P["hit"], P["null"]]
    bars = ax.bar(cats, rate, color=cols, edgecolor="white", lw=2.2, width=0.6, zorder=3)
    for b, r in zip(bars, rate):
        ax.text(b.get_x()+b.get_width()/2, r+0.7, f"{r:.1f}%", ha="center", fontweight="bold",
                fontsize=12, color=P["ink"])
    ax.annotate("", xy=(1, 27.1), xytext=(0, 5.3), zorder=4,
                arrowprops=dict(arrowstyle="-", color=P["hit"], lw=1.6, connectionstyle="arc3,rad=-0.3"))
    card(ax, 0.5, 20.5, "5.1\u00d7 enrichment\nFisher p = 2.6\u00d710\u207b\u2076", P["hit"], ha="center", fs=10)
    card(ax, 2, 6.6, "OR \u2248 1,  p = 1\nunderpowered (n = 8)", P["null"], ha="center", fs=8.5, pad=0.32)
    ax.set_ylabel("Known drug targets  (% of set)"); ax.set_ylim(0, 34); ax.grid(axis="x", alpha=0)
    ax.set_title("Hits are enriched 5.1\u00d7 for known drug targets", color=P["ink"], pad=16)
    ax.text(0.5, 1.0, "paralogue-driven \u2014 the family-removed set is null",
            transform=ax.transAxes, ha="center", fontsize=9.5, style="italic", color=P["null"])
    save(fig, "fig8_enrichment_null")

# ---------------------------------------------------------------- FIG 9 (pine/apricot)
def fig9():
    P = PAL["f9"]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5.0), gridspec_kw={"width_ratios": [0.85, 1.5]})
    tint(fig, ax1, P["bg"]); tint(fig, ax2, P["bg"])
    w, _ = ax1.pie([50, 8], colors=[P["within"], P["cross"]], startangle=90, counterclock=False,
                   wedgeprops=dict(width=0.40, edgecolor="white", lw=2.5))
    ax1.text(0, 0, "58\nhits", ha="center", va="center", fontsize=15, fontweight="bold", color=P["ink"])
    ax1.legend(w, ["Within-family (paralogue)  50", "Cross-family  8"],
               loc="upper center", bbox_to_anchor=(0.5, -0.04), frameon=True, fontsize=9.5)
    ax1.set_title("58 high-confidence hits\n86% within-family", color=P["ink"], pad=14)
    probes = ["ChEMBL\nbinding", "Recall", "Enrichment\nvs null", "Pathway\nenrichment"]
    txt = [["signal", "none"], ["recovered", "0 recall"], ["5.1\u00d7", "null"], ["family term", "too small"]]
    signal = [[1, 0], [1, 0], [1, 0], [1, 0]]
    for i in range(4):
        for j in range(2):
            c = P["yes"] if signal[i][j] else P["no"]
            ax2.add_patch(FancyBboxPatch((j+0.06, i+0.06), 0.88, 0.88,
                          boxstyle="round,pad=0.02,rounding_size=0.08", fc=c, ec="white", lw=2.5,
                          alpha=0.92, zorder=2))
            ax2.text(j+0.5, i+0.5, txt[i][j], ha="center", va="center", color="white",
                     fontweight="bold", fontsize=10.5, zorder=3)
    ax2.set_xlim(0, 2); ax2.set_ylim(0, 4)
    ax2.set_xticks([0.5, 1.5]); ax2.set_xticklabels(["Within-family\n(expected)", "Cross-family\n(novel)"],
                                                    fontsize=10.5, fontweight="bold")
    ax2.set_yticks([i+0.5 for i in range(4)]); ax2.set_yticklabels(probes, fontsize=10)
    ax2.invert_yaxis(); ax2.grid(False)
    for s in ax2.spines.values(): s.set_visible(False)
    ax2.set_title("All four probes agree: signal is within-family,\nno novel cross-family biology",
                  color=P["ink"], pad=12)
    save(fig, "fig9_family_split_verdict")

FIGS = [fig1, fig2, fig3, fig4, fig5, fig6, fig7, fig8, fig9]
if __name__ == "__main__":
    for f in FIGS:
        try: f()
        except Exception as e: print("ERROR in", f.__name__, ":", repr(e))
    print("\nAll figures in ./figures/")
