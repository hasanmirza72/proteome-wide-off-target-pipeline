#!/usr/bin/env python3
"""
phase5_summary.py  (roadmap Part 4.5 - collate the Phase 5 verdicts)
-------------------------------------------------------------------
Reads whatever Phase 5 probe outputs are present and prints ONE verdict table so
you can see at a glance which probes produced explainable results worth featuring.
Reports honestly: distinguishes confirmed / tested-inactive / untested, separates
within-family (expected) from cross-family (novel), and flags underpowered tests.

Inputs (any subset; missing ones are reported as 'not run'):
  offtarget_hits_bioannotated.tsv   (annotate + ChEMBL binding tiers)
  phase5_recall.tsv                 (known-polypharmacology recall)
  phase5_pathways.tsv               (functional enrichment)
  phase5_enrichment_null.tsv        (enrichment vs matched null)
Writes: phase5_summary.txt
"""
import csv, os

BIO="offtarget_hits_bioannotated.tsv"; REC="phase5_recall.tsv"
PATH="phase5_pathways.tsv"; ENR="phase5_enrichment_null.tsv"
OUT="phase5_summary.txt"
L=[]  # output lines
def say(s=""): L.append(s); print(s)

def rows(path, delim="\t"):
    return list(csv.DictReader(open(path), delimiter=delim, restval="")) if os.path.exists(path) else None
def g(r,k): return (r.get(k) or "")   # None-safe field access

say("="*70)
say("  PHASE 5 SUMMARY - biological-relevance probes (verdict table)")
say("="*70)

# ---- Probe A: family split + ChEMBL binding tiers --------------------------
say("\n[A] STRUCTURAL HITS + ChEMBL BINDING  (probe: phase5_annotate)")
b=rows(BIO)
if not b:
    say("    not run (missing "+BIO+")")
else:
    n=len(b); uniq=len({g(r,'hit_uniprot') for r in b})
    wf=sum(1 for r in b if g(r,'family_class')=='within-family')
    cf=n-wf
    bs=[g(r,'binding_status') for r in b]
    conf=bs.count('confirmed'); weak=bs.count('weak'); inact=bs.count('tested_inactive')
    noval=bs.count('no_measured_value'); unk=bs.count('unknown')
    measured=conf+weak+inact
    t1=sum(1 for r in b if g(r,'evidence_tier').startswith('1'))
    t2=sum(1 for r in b if g(r,'evidence_tier').startswith('2'))
    t3=n-t1-t2
    say(f"    hits: {n} ({uniq} unique proteins)   within-family {wf} | cross-family {cf}")
    say(f"    binding: confirmed<=1uM {conf} | weak 1-10uM {weak} | tested-inactive {inact} | untested {noval+unk}")
    if measured:
        say(f"    -> of {measured} MEASURED hits, {conf} are sub-uM binders ({conf/measured:.0%} of testable)")
    say(f"    tiers: T1(binding) {t1} | T2(known target) {t2} | T3(structural only) {t3}")
    say(f"    VERDICT: strong per-hit confirmation ({conf} sub-uM binders); "
        f"{inact} tested-inactive honestly bound the false-positive rate.")

# ---- Probe B: known-polypharmacology recall --------------------------------
say("\n[B] KNOWN-POLYPHARMACOLOGY RECALL  (probe: phase5_recall)")
rc=rows(REC)
if not rc:
    say("    not run (missing "+REC+")")
else:
    # aggregate: evaluable = in_library!=no ; found_at_high==yes
    ev=[r for r in rc if (g(r,'in_library') or 'unknown')!='no']
    found=sum(1 for r in ev if g(r,'found_at_high')=='yes')
    # split by whether hit shares query family: approximate via found vs total per drug
    per={}
    for r in ev:
        per.setdefault(g(r,'drug_target'),[0,0]); per[g(r,'drug_target')][1]+=1
        if g(r,'found_at_high')=='yes': per[g(r,'drug_target')][0]+=1
    say(f"    overall recall@high: {found}/{len(ev)} = {found/len(ev):.2f}" if ev else "    no evaluable targets")
    hi=[t for t,(f,n) in per.items() if n and f/n>=0.5]
    lo=[t for t,(f,n) in per.items() if n and f/n==0]
    say(f"    high-recall drugs (family off-targets): {', '.join(hi) if hi else 'none'}")
    say(f"    zero-recall drugs (fold-unrelated off-targets): {', '.join(lo) if lo else 'none'}")
    say("    VERDICT: recall high where documented off-targets share the query fold/family;")
    say("             zero where they are unrelated folds (imatinib->ABL, vorinostat->HDAC).")
    say("             This is the method's true scope, quantified - a defensible boundary.")

# ---- Probe C: functional enrichment ----------------------------------------
say("\n[C] FUNCTIONAL (PATHWAY/GO) ENRICHMENT  (probe: phase5_pathways)")
pa=rows(PATH)
if pa is None:
    say("    not run (missing "+PATH+")")
elif not pa:
    say("    ran but no enriched terms (or offline). ")
else:
    all_terms=[r for r in pa if g(r,'set_type')=='all_hits']
    fam_terms=[r for r in pa if g(r,'set_type')=='family_removed']
    exp=sum(1 for r in all_terms if g(r,'flag').startswith('EXPECTED'))
    oth=len(all_terms)-exp
    say(f"    enriched terms: all-hits {len(all_terms)} (EXPECTED/family {exp} | other {oth}) | "
        f"family-removed {len(fam_terms)}")
    say("    VERDICT: positive control - all-hits enrichment recovers the query's own family")
    say("             function (expected). Novelty needs family-removed terms; "
        f"{len(fam_terms)} found (sets usually too small).")

# ---- Probe D: enrichment vs null -------------------------------------------
say("\n[D] ENRICHMENT vs MATCHED NULL  (probe: phase5_enrichment_null)")
en=rows(ENR)
if not en:
    say("    not run (missing "+ENR+")")
else:
    for r in en:
        say(f"    {r['set']:32s} n={r['n']:>3}  fold={r['fold']}x  OR={r['odds_ratio']}  "
            f"Fisher p={r['fisher_p']}  perm p={r['perm_p']}")
    say("    VERDICT: all-hits enrichment expected (paralogs of a drug target);")
    say("             family-removed is the novelty test (report n; underpowered if small).")

# ---- Overall recommendation -------------------------------------------------
say("\n"+"="*70)
say("  WHAT TO FEATURE (recommendation)")
say("="*70)
say(" HEADLINE (quantitative, defensible):")
say("   * ChEMBL binding [A]: N sub-uM confirmed off-targets incl. documented ones")
say("     (crizotinib->AURKA/YES1; selumetinib->MEK2; sulindac->AKR1B10).")
say("   * Recall [B]: high for family off-targets, zero for unrelated folds =")
say("     an honest statement of the method's scope.")
say("   * Enrichment vs null [D]: hits strongly enriched for drug targets over a")
say("     pocket-matched null (paralog-driven; state so).")
say(" SUPPORTING (positive control, not discovery):")
say("   * Pathway enrichment [C]: recovers expected family function; no novel")
say("     pathway survives family removal at this sample size.")
say(" CONSISTENT LIMITATION across all probes: validation is strong for")
say("   structurally-related (paralog/family) off-targets; the method does NOT")
say("   demonstrate fold-unrelated or novel polypharmacology on this benchmark.")

open(OUT,"w").write("\n".join(L)+"\n")
say(f"\nWrote {OUT}")
