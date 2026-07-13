#!/usr/bin/env python3
"""
dual_query_divergence.py  (roadmap Part 3.5)
--------------------------------------------
Does querying the bound (holo) vs unbound (apo) pocket change the off-target list,
and is that change predicted by how much the pocket physically moved (Phase 3 axes)?

Reads offtarget_hits_annotated.tsv (all hits, both conformations) and, if present,
self_recovery_summary.tsv (self-rank shift). Conformational axes are embedded
(sources: apo spatial overlap from the manifest-build log; backbone_rmsd and
pocket_rmsd from dataset_final_nonredundant.csv).

Outputs (all four honest views of one question):
  1. dual_query_divergence.tsv  - FULL 23-row table, low-n rows flagged (view b)
  2. correlation on the reliable subset (both arms >= MIN_HITS non-self)  (view a)
  3. sensitivity strip across cutoffs (>=3 / >=5 / >=8)
  4. leverage check (re-run dropping the highest-movement targets)

p-values are permutation-based (numpy only; no distributional assumption).
"""
import csv, os, sys
import numpy as np

ANNOT   = "offtarget_hits_annotated.tsv"
SELF    = "self_recovery_summary.tsv"
OUT     = "dual_query_divergence.tsv"
MIN_HITS = 5           # a target's Jaccard is "reliable" iff BOTH arms have >= this
N_PERM  = 20000
RNG     = np.random.default_rng(0)

# --- Phase 3 conformational axes (embedded; see docstring for sources) --------
SPATIAL_OVERLAP = {"6VCJ":1.00,"1UEI":0.65,"1NHZ":0.72,"6L6E":0.94,"4U7Z":0.44,
"1XBB":0.92,"3TVX":1.00,"3RX3":0.90,"5V1M":1.00,"6Q6O":0.50,"6GH9":0.83,"8HJE":0.93,
"2XP2":1.00,"2Y05":0.88,"7AXA":0.92,"1AVN":1.00,"7XPY":1.00,"6YG2":1.00,"4ODR":0.79,
"3TC5":1.00,"1GS4":1.00,"4R7L":1.00,"9XRI":0.82}
BACKBONE_RMSD = {"6VCJ":0.99,"1UEI":2.08,"1NHZ":2.52,"6L6E":1.24,"4U7Z":2.20,"1XBB":0.56,
"3TVX":0.33,"3RX3":1.35,"5V1M":0.73,"6Q6O":0.94,"6GH9":1.14,"8HJE":1.50,"2XP2":0.50,
"2Y05":0.58,"7AXA":0.97,"1AVN":0.32,"7XPY":2.63,"6YG2":1.82,"4ODR":0.87,"3TC5":0.53,
"1GS4":0.70,"4R7L":0.15,"9XRI":1.02}
POCKET_RMSD = {"6VCJ":1.03,"1UEI":2.44,"1NHZ":2.68,"6L6E":1.85,"4U7Z":3.37,"1XBB":0.87,
"3TVX":2.10,"3RX3":4.85,"5V1M":1.12,"6Q6O":0.57,"6GH9":0.93,"8HJE":2.32,"2XP2":0.63,
"2Y05":1.34,"7AXA":1.12,"1AVN":0.33,"7XPY":0.80,"6YG2":2.01,"4ODR":0.87,"3TC5":1.57,
"1GS4":0.68,"4R7L":0.77,"9XRI":0.93}
# fallback self-ranks (used only if self_recovery_summary.tsv is absent)
SELF_FALLBACK = {"6VCJ":(2,2),"1UEI":(1,2),"1NHZ":(1,1),"6L6E":(1,1),"4U7Z":(2,164),
"1XBB":(1,1),"3TVX":(3,4),"3RX3":(1,1),"5V1M":(1,1),"6Q6O":(8,1),"6GH9":(1,2),"8HJE":(2,1),
"2XP2":(2,2),"2Y05":(1,1),"7AXA":(1,1),"1AVN":(1,1),"7XPY":(1,1),"6YG2":(1,1),"4ODR":(2,3),
"3TC5":(1,1),"1GS4":(1,1),"4R7L":(1,1),"9XRI":(1,1)}
# -----------------------------------------------------------------------------

def rankdata(a):
    a = np.asarray(a, float); order = a.argsort(); ranks = np.empty_like(order, float)
    ranks[order] = np.arange(len(a))
    # average ties
    _, inv, cnt = np.unique(a, return_inverse=True, return_counts=True)
    sums = np.zeros(len(cnt)); np.add.at(sums, inv, ranks)
    return (sums/cnt)[inv] + 1.0

def spearman(x, y):
    x, y = np.asarray(x,float), np.asarray(y,float)
    if len(x) < 3 or np.all(x==x[0]) or np.all(y==y[0]): return np.nan
    rx, ry = rankdata(x), rankdata(y)
    rx-=rx.mean(); ry-=ry.mean()
    d = np.sqrt((rx*rx).sum()*(ry*ry).sum())
    return float((rx*ry).sum()/d) if d>0 else np.nan

def perm_p(x, y, rho, n=N_PERM):
    if np.isnan(rho): return np.nan
    x=np.asarray(x,float); y=np.asarray(y,float); c=0
    for _ in range(n):
        if abs(spearman(x, RNG.permutation(y))) >= abs(rho)-1e-12: c+=1
    return (c+1)/(n+1)

def load_self():
    d={}
    if os.path.exists(SELF):
        with open(SELF) as f:
            for r in csv.DictReader(f, delimiter="\t"):
                try: d[r["target"]]=(_num(r["holo_rank"]), _num(r["apo_rank"]))
                except: pass
        if d: return d
    return SELF_FALLBACK

def _num(x):
    try: return int(x)
    except: 
        try: return float(x)
        except: return None

# ---- build per-target holo/apo hit dicts (uniprot -> idf), self excluded -----
holo, apo = {}, {}
with open(ANNOT) as f:
    for r in csv.DictReader(f, delimiter="\t"):
        if r["is_self"].strip() in ("1","True","true"): continue
        t=r["query_target"]; conf=r["query_conf"]; u=r["hit_uniprot"]
        try: idf=float(r["idf"])
        except: idf=0.0
        d = holo if conf=="holo" else apo
        d.setdefault(t, {})[u]=idf

targets = sorted(set(holo)|set(apo))
selfr = load_self()

rows=[]
for t in targets:
    h=holo.get(t,{}); a=apo.get(t,{})
    sh=set(h); sa=set(a); shared=sh&sa; union=sh|sa
    jac = len(shared)/len(union) if union else np.nan
    # overlap coefficient: robust to list-size asymmetry (shared / smaller list)
    ocoef = len(shared)/min(len(sh),len(sa)) if (sh and sa) else np.nan
    # ranking stability on shared hits (spearman of idf among shared)
    if len(shared)>=3:
        sl=sorted(shared)
        rho_sh = spearman([h[u] for u in sl],[a[u] for u in sl])
    else:
        rho_sh = np.nan
    hr = selfr.get(t,(None,None))
    shift = (hr[1]-hr[0]) if (hr[0] is not None and hr[1] is not None) else None
    ratio = (len(a)/len(h)) if len(h)>0 else np.nan
    reliable = (len(h)>=MIN_HITS and len(a)>=MIN_HITS)
    rows.append(dict(target=t, n_holo=len(h), n_apo=len(a), n_shared=len(shared),
        jaccard=jac, overlap_coef=ocoef, reliable=reliable, spearman_shared=rho_sh,
        self_rank_holo=hr[0], self_rank_apo=hr[1], self_rank_shift=shift,
        hit_ratio=ratio, spatial_overlap=SPATIAL_OVERLAP.get(t),
        backbone_rmsd=BACKBONE_RMSD.get(t), pocket_rmsd=POCKET_RMSD.get(t)))

# ---- write full table (view b) ----------------------------------------------
cols=["target","n_holo","n_apo","n_shared","jaccard","overlap_coef","reliable","spearman_shared",
      "self_rank_holo","self_rank_apo","self_rank_shift","hit_ratio",
      "spatial_overlap","backbone_rmsd","pocket_rmsd"]
def fmt(v):
    if v is None: return ""
    if isinstance(v,float): return "" if np.isnan(v) else f"{v:.3f}"
    return str(v)
with open(OUT,"w") as f:
    f.write("\t".join(cols)+"\n")
    for r in rows: f.write("\t".join(fmt(r[c]) for c in cols)+"\n")

# ---- descriptive ------------------------------------------------------------
def subset(min_hits):
    return [r for r in rows if r["n_holo"]>=min_hits and r["n_apo"]>=min_hits]

rel=subset(MIN_HITS)
jvals=[r["jaccard"] for r in rel]
print(f"=== Part 3.5 dual-query divergence ===")
print(f"Full table: {OUT}  ({len(rows)} targets; {len(rel)} reliable at both arms >= {MIN_HITS})\n")
print(f"Jaccard over reliable subset (n={len(rel)}): median {np.median(jvals):.3f}, "
      f"range {min(jvals):.3f}-{max(jvals):.3f}")
ovals=[r["overlap_coef"] for r in rel]
print(f"Overlap coefficient (size-robust) over reliable subset: median {np.median(ovals):.3f}, range {min(ovals):.3f}-{max(ovals):.3f}")
print("Most conformation-sensitive (lowest Jaccard):")
for r in sorted(rel,key=lambda x:x["jaccard"])[:5]:
    print(f"   {r['target']:5s} J={r['jaccard']:.3f}  holo={r['n_holo']} apo={r['n_apo']} "
          f"shared={r['n_shared']}  overlap={r['spatial_overlap']} bbRMSD={r['backbone_rmsd']}")

# ---- correlation: does movement predict divergence? (view a) ----------------
def corr(sub, axis):
    xs=np.array([r[axis] for r in sub]); ys=np.array([r["jaccard"] for r in sub])
    rho=spearman(xs,ys); return rho, perm_p(xs,ys,rho), len(sub)

print(f"\n=== Correlation: Jaccard vs conformational axes (reliable subset, n={len(rel)}) ===")
print("(expect: +ve vs spatial_overlap [stable pocket->similar list]; -ve vs RMSDs)")
for axis in ["spatial_overlap","backbone_rmsd","pocket_rmsd"]:
    rho,p,n=corr(rel,axis); print(f"   Jaccard ~ {axis:16s}: rho={rho:+.3f}  perm_p={p:.4f}  (n={n})")

# ---- sensitivity strip across cutoffs ---------------------------------------
print(f"\n=== Sensitivity: rho(Jaccard ~ backbone_rmsd) at different hit cutoffs ===")
for mh in (3,5,8):
    sub=subset(mh); rho,p,n=corr(sub,"backbone_rmsd")
    print(f"   cutoff >= {mh}: rho={rho:+.3f}  perm_p={p:.4f}  (n={n})")

# ---- leverage check: drop the highest-movement targets ----------------------
print(f"\n=== Leverage check (reliable subset, drop 3 highest backbone_rmsd) ===")
movers=sorted(rel,key=lambda x:-x["backbone_rmsd"])[:3]
mv={r["target"] for r in movers}
print("   dropped movers:", ", ".join(f"{r['target']}({r['backbone_rmsd']})" for r in movers))
keep=[r for r in rel if r["target"] not in mv]
for axis in ["spatial_overlap","backbone_rmsd","pocket_rmsd"]:
    rho,p,n=corr(keep,axis); print(f"   without movers: Jaccard ~ {axis:16s}: rho={rho:+.3f} perm_p={p:.4f} (n={n})")
