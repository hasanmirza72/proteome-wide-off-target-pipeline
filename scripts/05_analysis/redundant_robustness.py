#!/usr/bin/env python3
"""
redundant_robustness.py  (roadmap Part 2.6 robustness sub-study)
----------------------------------------------------------------
Question: for a protein co-crystallised with several different drugs, how much does
the off-target hit list depend on WHICH drug's pocket defined the query?

Reads the non-redundant representatives' holo hit lists from
offtarget_hits_annotated.tsv and the extra redundant entries from
hits_redundant/<target>_holo.tsv. For each protein it computes all-pairs Jaccard
and overlap coefficient across its drugs, reports the mean within-protein agreement,
and prints each entry's ligand atom count and hit count so the ligand-size confound
is visible. Compare the numbers here to the cross-conformation divergence of Part
3.5 (median Jaccard ~0.07).

CASE STUDY, not a proteome claim: 3 deeply-sampled proteins (CA2, PDE5A, PXR).
"""
import csv, os, re
from itertools import combinations
import numpy as np

ANNOT   = "offtarget_hits_annotated.tsv"
REDDIR  = "hits_redundant"
OUT     = "redundant_robustness.tsv"

# protein -> list of (entry, ligand, n_lig_atoms, source, candidate_source)
# source: 'annot' = read from offtarget_hits_annotated.tsv ; 'redund' = hits_redundant/
GROUPS = {
    "CA2 (P00918)": [
        ("1AVN","HSM", 8,"annot","drug-matched"),
        ("1OQ5","CEL",26,"redund","drug-matched"),
        ("2AW1","COX",22,"redund","drug-matched"),
        ("2POU","I7A",16,"redund","drug-matched"),
        ("1CIL","ETS",19,"redund","drug-matched"),
    ],
    "PDE5A (O76074)": [
        ("6L6E","E6L",34,"annot","drug-matched"),
        ("1XOZ","CIA",29,"redund","drug-matched"),
        ("1UHO","VDN",34,"redund","any-ligand-pocket"),
    ],
    "PXR/NR1I2 (O75469)": [
        ("7AXA","CL6",50,"annot","drug-matched"),
        ("7AXK","EST",20,"redund","drug-matched"),
        ("6HJ2","P06",35,"redund","drug-matched"),
    ],
}
PROT_UNIPROT = {"CA2 (P00918)":"P00918","PDE5A (O76074)":"O76074","PXR/NR1I2 (O75469)":"O75469"}

_AF = re.compile(r"AF-([A-Za-z0-9]+)-F\d+-model")
def tid_to_uniprot(tid):
    b = os.path.basename(tid.strip()); m = _AF.search(b)
    if m: return m.group(1)
    for e in (".pdb.gz",".cif.gz",".pdb",".cif"):
        if b.endswith(e): return b[:-len(e)]
    return b

def hits_from_annot(entry, self_uni):
    s = set()
    if not os.path.exists(ANNOT): return s
    with open(ANNOT) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["query_target"]==entry and r["query_conf"]=="holo" and r["is_self"].strip() not in ("1","True","true"):
                s.add(r["hit_uniprot"].strip())
    return s

def hits_from_redund(entry, self_uni):
    s = set(); path = os.path.join(REDDIR, f"{entry}_holo.tsv")
    if not os.path.exists(path): return None            # not run yet
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    if len(lines) < 2: return s
    delim = "\t" if "\t" in lines[0] else None
    header = lines[0].split(delim) if delim else lines[0].split()
    ti = header.index("tid") if "tid" in header else 0
    for ln in lines[1:]:
        p = ln.split(delim) if delim else ln.split()
        if ti < len(p):
            acc = tid_to_uniprot(p[ti])
            if acc != self_uni:                          # exclude on-target self
                s.add(acc)
    return s

def jaccard(a,b): return len(a&b)/len(a|b) if (a|b) else float("nan")
def ocoef(a,b):   return len(a&b)/min(len(a),len(b)) if (a and b) else float("nan")

def main():
    out_rows = ["\t".join(["protein","drug_a","drug_b","n_a","n_b","shared","jaccard","overlap_coef"])]
    print("=== Part 2.6 redundant robustness: same protein, different drugs ===")
    print("(compare to Part 3.5 cross-conformation median Jaccard ~0.07)\n")
    missing = []
    all_means = []
    for prot, entries in GROUPS.items():
        uni = PROT_UNIPROT[prot]
        sets = {}
        print(f"--- {prot} ---")
        for (entry,lig,natom,src,csrc) in entries:
            s = hits_from_annot(entry,uni) if src=="annot" else hits_from_redund(entry,uni)
            if s is None:
                missing.append(entry); print(f"   {entry} ({lig}): [hits_redundant/{entry}_holo.tsv NOT FOUND - run run_redundant.sh]"); continue
            sets[entry]=s
            flag = "" if csrc=="drug-matched" else f"  [{csrc}]"
            print(f"   {entry:5s} {lig:4s} lig_atoms={natom:<3d} hits={len(s):<4d}{flag}")
        # all-pairs within this protein
        js=[]; 
        for a,b in combinations(sets,2):
            J=jaccard(sets[a],sets[b]); O=ocoef(sets[a],sets[b])
            js.append(J)
            out_rows.append("\t".join([prot,a,b,str(len(sets[a])),str(len(sets[b])),
                              str(len(sets[a]&sets[b])),f"{J:.3f}",f"{O:.3f}"]))
        valid=[j for j in js if j==j]  # drop nan
        if valid:
            m=float(np.mean(valid)); all_means.append((prot,m,len(valid)))
            print(f"   -> mean pairwise Jaccard across {len(sets)} drugs "
                  f"({len(valid)} valid pairs): {m:.3f}\n")
        else:
            print("   -> insufficient hits for a within-protein comparison "
                  "(too few off-target hits per drug)\n")

    with open(OUT,"w") as f: f.write("\n".join(out_rows)+"\n")
    print("=== summary ===")
    for prot,m,npair in all_means:
        print(f"   {prot:22s} mean within-protein Jaccard = {m:.3f}  ({npair} pairs)")
    print(f"\nWrote {OUT}")
    if missing:
        print(f"\n[!] {len(missing)} entries not yet queried: {', '.join(missing)}")
        print("    Run build_redundant_manifest.py -> run_redundant.sh first.")
    else:
        print("\nNote the ligand-atom column: larger ligands define larger motifs, so some")
        print("within-protein divergence is ligand-size-driven, not pocket-instability-driven.")

if __name__ == "__main__":
    main()
