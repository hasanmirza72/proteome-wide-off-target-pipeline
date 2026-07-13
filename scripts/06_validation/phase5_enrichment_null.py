#!/usr/bin/env python3
"""
phase5_enrichment_null.py  (roadmap Part 4.3 - enrichment vs a matched null)
---------------------------------------------------------------------------
Honest headline metric: are the high-confidence structural hits enriched for
KNOWN DRUG TARGETS relative to a matched random draw from the SEARCHABLE,
POCKET-BEARING library (the 16,325 proteins in library_pockets.tsv)?

Property tested: is_known_drug_target = the protein has a ChEMBL drug MECHANISM.
It is drug-INDEPENDENT and defined identically for hits and background, so ChEMBL
coverage gaps affect both arms equally (the property is robust to incompleteness).

WHY THE BACKGROUND IS THE 16,325 POCKET-BEARING SET (not all 20,550): the hits are
already filtered to pocket-supported proteins, so drawing the null from all proteins
would make enrichment circular with that filter. Sampling from pocket-bearing
proteins controls for "has a cavity" and isolates "is a real drug target".

Runs the test TWICE: all high hits, and FAMILY-REMOVED (cross-family only). All-hits
enrichment is expected (paralogs of a drug target are drug targets); the
family-removed test is the novelty test and is reported with its (small) n.

One-time cost: annotates N_BG random library proteins via ChEMBL (cached/resumable).
Stats: Fisher's exact (odds ratio + p) primary; permutation p as a distribution-free
check. Reads offtarget_hits_bioannotated.tsv + library_pockets.tsv.
Writes phase5_enrichment_null.tsv.
"""
import csv, os, json, time, sys, hashlib
import numpy as np

BIO   = "offtarget_hits_bioannotated.tsv"
LIB   = "library_pockets.tsv"
OUT   = "phase5_enrichment_null.tsv"
CACHE = "chembl_cache"; os.makedirs(CACHE, exist_ok=True)  # shared with annotate
BASE  = "https://www.ebi.ac.uk/chembl/api/data"
N_BG  = 1000          # random background proteins to annotate (one-time, cached)
N_PERM= 20000
RNG   = np.random.default_rng(0)

def _get(url):
    key=os.path.join(CACHE, hashlib.md5(url.encode()).hexdigest()+".json")
    if os.path.exists(key): return json.load(open(key))
    try:
        import requests
        for a in range(4):
            r=requests.get(url, headers={"Accept":"application/json"}, timeout=30)
            if r.status_code==200:
                d=r.json(); json.dump(d, open(key,"w")); return d
            if r.status_code in (429,500,502,503): time.sleep(2*(a+1)); continue
            return {"_error":r.status_code}
        return {"_error":"retries"}
    except Exception as e:
        return {"_offline":str(e)}

def chembl_target(uniprot):
    d=_get(f"{BASE}/target?target_components__accession={uniprot}&format=json&limit=1")
    if isinstance(d,dict) and d.get("targets"):
        return d["targets"][0].get("target_chembl_id"), None
    return None, (d.get("_offline") or d.get("_error") if isinstance(d,dict) else "err")

def is_known_drug_target(uniprot):
    """1 if the protein has a ChEMBL drug mechanism; 0 if not; None if offline."""
    tid,err=chembl_target(uniprot)
    if err and tid is None: return None
    if not tid: return 0
    d=_get(f"{BASE}/mechanism?target_chembl_id={tid}&format=json&limit=1")
    if isinstance(d,dict) and "mechanisms" in d:
        return 1 if d["mechanisms"] else 0
    if isinstance(d,dict) and (d.get("_offline")): return None
    return 0

def fisher(a,b,c,d):
    """2x2 [[a,b],[c,d]] -> (odds_ratio, p_two_sided), self-contained (no scipy).
    p via hypergeometric tail sum (Fisher's exact); robust to zero cells."""
    from math import lgamma, exp
    def logC(n,k):
        if k<0 or k>n: return float("-inf")
        return lgamma(n+1)-lgamma(k+1)-lgamma(n-k+1)
    n=a+b+c+d; r1=a+b; c1=a+c
    if 0 in (a,b,c,d):
        orr=((a+0.5)*(d+0.5))/((b+0.5)*(c+0.5))   # Haldane-Anscombe correction
    else:
        orr=(a*d)/(b*c)
    if n==0 or r1 in (0,n) or c1 in (0,n):
        return orr, float("nan")
    def hyp(k): return exp(logC(c1,k)+logC(n-c1,r1-k)-logC(n,r1))
    p_obs=hyp(a); tol=1e-9
    lo=max(0,r1-(n-c1)); hi=min(r1,c1)
    p=sum(hyp(k) for k in range(lo,hi+1) if hyp(k)<=p_obs+tol)
    return orr, min(1.0,p)

def main():
    if not os.path.exists(BIO): sys.exit(f"missing {BIO}")
    if not os.path.exists(LIB): sys.exit(f"missing {LIB} (the 16,325 pocket-bearing background)")

    # --- hits: unique accessions + known-drug-target flag from the annotated file ---
    hit_known={}; cross=set()
    for r in csv.DictReader(open(BIO), delimiter="\t"):
        acc=r["hit_uniprot"].strip()
        hit_known[acc]= 1 if r["is_known_drug_target"]=="yes" else 0
        if r["family_class"]=="cross-family": cross.add(acc)
    all_hits=sorted(hit_known); cross_hits=sorted(cross)

    # --- background pool: random pocket-bearing library proteins, annotated once ---
    lib=[]
    with open(LIB) as f:
        rd=csv.reader(f, delimiter="\t"); next(rd,None)
        for row in rd:
            if row and row[0].strip(): lib.append(row[0].strip())
    lib=sorted(set(lib))
    hitset=set(all_hits)
    pool=[a for a in lib if a not in hitset]            # background excludes the hits
    RNG.shuffle(pool)
    sample=pool[:N_BG]
    print(f"Annotating {len(sample)} random background proteins (cached; one-time)...")
    bg=[]; offline=False
    for i,acc in enumerate(sample):
        v=is_known_drug_target(acc)
        if v is None: offline=True; break
        bg.append(v)
        if (i+1)%100==0: print(f"   {i+1}/{len(sample)}")
    if offline:
        print("\n[!] ChEMBL not reachable - run on a networked node. Cache makes reruns free.")
        sys.exit(0)

    bg=np.array(bg); bg_rate=bg.mean(); bg_yes=int(bg.sum()); bg_no=len(bg)-bg_yes

    def test(hits, label):
        vals=np.array([hit_known[a] for a in hits])
        hy=int(vals.sum()); hn=len(vals)-hy; hrate=vals.mean() if len(vals) else float("nan")
        orr,p = fisher(hy,hn,bg_yes,bg_no)
        fold = (hrate/bg_rate) if bg_rate>0 else float("nan")
        # permutation: draw len(hits) from background pool, count known targets
        obs=hy; k=len(hits); c=0
        for _ in range(N_PERM):
            draw=bg[RNG.choice(len(bg), size=k, replace=False)]
            if draw.sum() >= obs: c+=1
        permp=(c+1)/(N_PERM+1)
        print(f"\n[{label}]  n={len(hits)}  known-target={hy} ({hrate:.2%})   "
              f"background={bg_rate:.2%}")
        print(f"   fold enrichment = {fold:.2f}x   odds ratio = {orr:.2f}   "
              f"Fisher p = {p:.4g}   perm p = {permp:.4g}")
        return dict(label=label,n=len(hits),known=hy,rate=hrate,bg_rate=bg_rate,
                    fold=fold,odds=orr,fisher_p=p,perm_p=permp)

    print(f"\nBackground known-drug-target rate: {bg_yes}/{len(bg)} = {bg_rate:.2%}")
    r_all=test(all_hits,"ALL_HITS")
    r_cross=test(cross_hits,"FAMILY_REMOVED (cross-family)") if len(cross_hits)>=2 else None

    with open(OUT,"w") as f:
        f.write("set\tn\tknown_targets\thit_rate\tbg_rate\tfold\todds_ratio\tfisher_p\tperm_p\n")
        for r in [r_all,r_cross]:
            if r: f.write("\t".join(str(x) for x in [r['label'],r['n'],r['known'],
                f"{r['rate']:.3f}",f"{r['bg_rate']:.3f}",f"{r['fold']:.2f}",
                f"{r['odds']:.2f}",f"{r['fisher_p']:.4g}",f"{r['perm_p']:.4g}"])+"\n")
    print(f"\nWrote {OUT}")
    print("\nHONEST READING:")
    print(" - ALL_HITS enrichment is EXPECTED: the hits are dominated by paralogs of the query")
    print("   (itself a drug target), so they inherit drug-target status. It confirms the hits")
    print("   are biologically non-random, but it is not novel-target discovery.")
    print(" - FAMILY_REMOVED is the novelty test. With only", len(cross_hits),
          "cross-family hits it is")
    print("   underpowered; report the number but do not over-read a null result.")
    print(" - Property defined identically for hits and background => robust to ChEMBL gaps.")

if __name__ == "__main__":
    main()
