#!/usr/bin/env python3
"""Generate the methods workflow figure (fig0_workflow.png) for the thesis.
The figure shows the full pipeline: the dual holo/apo query, the distinct roles of
FoldDisco (structural search) and P2Rank (druggability filter), the funnel of candidate
counts (1,107 -> 681 -> 58), and the four biological validation probes.
Run from a folder containing a writable figures/ subdirectory:
    python3 make_workflow_figure.py
"""
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
mpl.rcParams.update({"font.family":"DejaVu Sans","savefig.dpi":300,"savefig.bbox":"tight"})

INK="#1B2D2A"; BG="#FBF8F3"; MUTED="#6B7C78"
STAGE="#2D6A4F"; SEARCH="#5A189A"; FILTER="#0077B6"; VALID="#2A9D8F"; DATA="#B5651D"

fig, ax = plt.subplots(figsize=(12, 8.4))
fig.patch.set_facecolor("white"); ax.set_facecolor(BG)
ax.set_xlim(0,100); ax.set_ylim(0,100); ax.axis("off")

def box(x,y,w,h,title,color,sub=None,textcol="white",fs=11,subfs=8.6):
    b=FancyBboxPatch((x-w/2,y-h/2),w,h,boxstyle="round,pad=0.3,rounding_size=1.4",
                     fc=color,ec="white",lw=2,zorder=3)
    ax.add_patch(b)
    if sub:
        ax.text(x,y+h*0.16,title,ha="center",va="center",color=textcol,fontsize=fs,
                fontweight="bold",zorder=4)
        ax.text(x,y-h*0.24,sub,ha="center",va="center",color=textcol,fontsize=subfs,
                zorder=4,linespacing=1.25,alpha=0.95)
    else:
        ax.text(x,y,title,ha="center",va="center",color=textcol,fontsize=fs,
                fontweight="bold",zorder=4,linespacing=1.25)

def arrow(x0,y0,x1,y1,color=INK,lw=2.4,style="-|>"):
    ax.add_patch(FancyArrowPatch((x0,y0),(x1,y1),arrowstyle=style,mutation_scale=20,
                 color=color,lw=lw,zorder=2,shrinkA=3,shrinkB=3))

def tag(x,y,text,color):
    ax.text(x,y,text,ha="center",va="center",fontsize=8.2,color=color,style="italic",
            fontweight="bold",zorder=5)

# ---- Title ----
ax.text(50,97,"Proteome-wide structural off-target pipeline",ha="center",
        fontsize=14.5,fontweight="bold",color=INK)
ax.text(50,92.5,"FoldDisco searches for pocket geometry;  P2Rank filters for druggable cavities  (distinct, complementary roles)",
        ha="center",fontsize=9.3,color=MUTED,style="italic")

# ================= LEFT COLUMN: one-time library prep =================
tag(16,86,"BUILD ONCE",DATA)
box(16,80,26,9,"AlphaFold human proteome","#8A5A2B",sub="23,586 models → 20,550 first-fragment\n(99.31% reviewed coverage)",fs=10)
arrow(16,75.5,16,70)
box(16,64,26,9,"P2Rank pocket library",FILTER,sub="cavity prediction per structure\n16,325 proteins with ≥1 pocket",fs=10)

# ================= CENTRE: per-drug query =================
tag(50,86,"PER DRUG TARGET",STAGE)
box(50,80,30,9.5,"Query pocket = ligand-contact residues",STAGE,
    sub="defined from the solved complex\n(not a predicted pocket)",fs=10)
# dual conformation split
arrow(50,75,40,69); arrow(50,75,60,69)
box(38,64,17,8,"HOLO query","#3A7D5C",sub="bound state",fs=9.6,subfs=8)
box(62,64,17,8,"APO query","#52946F",sub="unbound state",fs=9.6,subfs=8)
tag(50,58.5,"two conformations searched separately",MUTED)

# ================= FoldDisco search =================
arrow(38,60,46,52); arrow(62,60,54,52)
box(50,46,34,10,"FoldDisco motif search",SEARCH,
    sub='"where else does this residue geometry recur?"\nranked hits across all 20,550 models',fs=11,subfs=8.4)
# library feeds the search (single clean elbow from bottom of library column)
arrow(16,59.5,16,49,color=MUTED,lw=1.8,style="-")
arrow(16,49,33,47,color=MUTED,lw=1.8,style="->")
ax.text(21,50.5,"searched against",ha="center",fontsize=7.6,color=MUTED,style="italic")

# ================= P2Rank filter (distinct role) =================
arrow(50,41,50,35.5)
box(50,30,34,9,"P2Rank druggability filter",FILTER,
    sub="keep hits whose matched residues fall in a real cavity",fs=10.5,subfs=8.2)

# ================= confidence tier / funnel =================
arrow(50,25.5,50,20.5)
box(50,15.5,34,9,"Confidence tier",STAGE,
    sub="1,107 raw → 681 pocket-supported → 58 high-confidence",fs=10.5,subfs=8.4)

# ================= RIGHT: validation probes =================
tag(88,60,"VALIDATE",VALID)
for i,(yy,txt) in enumerate([(52,"ChEMBL binding"),(45,"Documented\npolypharmacology"),
                              (38,"Enrichment vs null"),(31,"Functional pathways")]):
    box(88,yy,20,5.6,txt,VALID,fs=8.8)
# arrow from confidence tier to validation stack
arrow(67,15.5,88,27.5,color=VALID,lw=2)
ax.text(80,20,"58 candidates",ha="center",fontsize=8,color=VALID,style="italic",fontweight="bold")

# ---- outcome banner ----
box(50,5.5,66,6.5,"Outcome: specific detector of active-site-family off-targets","#1B2D2A",
    sub="high precision within family · fold-unrelated off-targets lie beyond the pocket signal",
    fs=10.5,subfs=8.3)

plt.savefig("figures/fig0_workflow.png",dpi=300,facecolor="white",bbox_inches="tight")
print("saved figures/fig0_workflow.png")
