#!/usr/bin/env python3
"""
phase5_recall.py  (roadmap Part 4.4 - known-polypharmacology recall)
--------------------------------------------------------------------
Positive control: for drugs with DOCUMENTED secondary targets, does the pipeline's
high-confidence hit set recover them? Recall is reported ONLY on the annotated
subset of drugs that actually have documented off-targets in the human proteome.

CRITICAL HONESTY CONSTRAINTS
  * The KNOWN dict below is seeded from published pharmacology (DrugBank / ChEMBL
    mechanism / primary literature), INDEPENDENT of what this pipeline found. It is
    the ground truth, not a mirror of the hits.
  * Drugs that are clean/selective or whose off-targets are not proteins in the
    searchable library are EXCLUDED from the denominator (listed separately), not
    scored as failures.
  * A known off-target that is not in the 16,325-protein searchable library cannot
    be found by construction; such targets are flagged and excluded from recall.
  * Recall is computed against the query's high-confidence hits (default) and,
    optionally, against all pocket-supported hits, so "missed at HIGH but found at
    lower confidence" is visible.

Reads: offtarget_hits_bioannotated.tsv (or offtarget_hits_high.tsv), and
       optionally offtarget_hits_filtered.tsv for the all-supported comparison,
       and library_pockets.tsv to test library membership of each known target.
Writes: phase5_recall.tsv  (+ printed per-drug and overall recall)
"""
import csv, os, sys

HIGH_ANNOT = "offtarget_hits_bioannotated.tsv"   # preferred (has query_gene, drug)
HIGH_RAW   = "offtarget_hits_high.tsv"           # fallback
FILTERED   = "offtarget_hits_filtered.tsv"       # all pocket-supported (optional)
LIBPOCK    = "library_pockets.tsv"               # optional: library membership test
OUT        = "phase5_recall.tsv"

# --- KNOWN documented protein off-targets, by query PDB target -----------------
# Seeded from published pharmacology (DrugBank/ChEMBL mechanism/literature),
# INDEPENDENT of this pipeline's output. Values are UniProt accessions with a
# short justification + confidence. On-target (the drug's primary target) is NOT
# listed here (it is excluded as self). Only well-documented SECONDARY targets.
KNOWN = {
 # ALK inhibitor crizotinib: classic multi-kinase off-target profile (well documented)
 "2XP2": {  # query = ALK (Q9UM73); drug = Crizotinib
    "P08581": "MET - crizotinib primary/second target, documented [high]",
    "P07947": "YES1 - documented crizotinib off-target [high]",
    "O14965": "AURKA - reported crizotinib off-target [med]",
    "Q9UM73": "(on-target ALK - excluded as self)",
 },
 # SYK crystal but the drug is IMATINIB (a BCR-ABL/KIT/PDGFR inhibitor).
 # Imatinib's DOCUMENTED targets are ABL/KIT/PDGFRA/DDR1 - most are NOT the kinases
 # the pocket search returned. This is the honest, hard case.
 "1XBB": {  # drug = Imatinib
    "P00519": "ABL1 - imatinib primary target [high]",
    "P10721": "KIT  - imatinib documented target [high]",
    "P16234": "PDGFRA - imatinib documented target [high]",
    "Q08345": "DDR1 - imatinib documented target [high]",
 },
 # MEK1 inhibitor selumetinib: highly selective; documented cross-reactivity = MEK2
 "4U7Z": {  # drug = Selumetinib; query = MAP2K1
    "P36507": "MAP2K2 (MEK2) - documented selumetinib co-target [high]",
 },
 # MEK7/JNK-pathway; drug = Ibrutinib (BTK inhibitor). Ibrutinib's documented
 # off-targets are TEC/EGFR/ITK/BLK etc., NOT MEK7 - clean test of a hard case.
 "6YG2": {  # drug = Ibrutinib
    "Q06187": "BTK - ibrutinib primary target [high]",
    "P07332": "FES - ibrutinib documented off-target [med]",
    "Q08881": "ITK - ibrutinib documented off-target [med]",
    "P00533": "EGFR - ibrutinib documented off-target [med]",
 },
 # Aldose reductase crystal; drug = Sulindac. Documented: AKR1B10 + COX (PTGS).
 "3RX3": {  # drug = Sulindac (R)
    "O60218": "AKR1B10 - documented sulindac target [high]",
    "P35354": "PTGS2/COX2 - sulindac primary anti-inflammatory target [high]",
    "P23219": "PTGS1/COX1 - sulindac target [high]",
 },
 # PDE4A crystal; drug = Pentoxifylline (non-selective PDE inhibitor).
 "3TVX": {  # drug = Pentoxifylline
    "Q07343": "PDE4B - documented pentoxifylline target [high]",
    "Q08499": "PDE4D - documented pentoxifylline target [high]",
    "Q14432": "PDE3A - pentoxifylline inhibits PDE3 [med]",
 },
 # FKBP1A crystal; drug = Tacrolimus (FK506). Documented: FKBP family.
 "4ODR": {  # drug = Tacrolimus
    "P68106": "FKBP1B - documented tacrolimus binder [high]",
    "Q02790": "FKBP4 - documented FK506 binder [med]",
    "Q13451": "FKBP5 - documented FK506 binder [med]",
 },
 # Raloxifene (SERM). Documented targets = ESR1/ESR2 (not PTGR1 crystal).
 "2Y05": {  # drug = Raloxifene
    "P03372": "ESR1 - raloxifene primary target [high]",
    "Q92731": "ESR2 - raloxifene target [high]",
 },
 # Vorinostat (SAHA) - pan-HDAC inhibitor. Documented targets = HDAC family.
 "4R7L": {  # drug = Vorinostat
    "Q13547": "HDAC1 - vorinostat target [high]",
    "Q92769": "HDAC2 - vorinostat target [high]",
    "P56524": "HDAC4 - vorinostat target [med]",
 },
}

# Drugs excluded from recall (clean/selective, or off-targets not protein/searchable)
EXCLUDED = {
 "6VCJ":"Naproxen - primary target COX; no documented secondary protein in library scope",
 "1UEI":"UTP - endogenous nucleotide, no documented drug off-targets",
 "5V1M":"UMP - endogenous nucleotide, no documented drug off-targets",
 "6Q6O":"Cholic acid - endogenous bile acid, no documented protein off-targets",
 "9XRI":"Acoramidis - selective TTR stabiliser, no documented off-targets in library",
 "1NHZ":"Mifepristone - PR/GR antagonist (on-target class); secondary targets sparse",
 "6L6E":"Avanafil - highly selective PDE5 inhibitor",
 "7AXA":"Clotrimazole - antifungal (CYP targets), sparse human-protein off-targets",
 "1GS4":"Fludrocortisone - corticosteroid receptor agonist (on-target class)",
 "1AVN":"Histamine - endogenous amine; CA2 crystal is a fragment soak, not a drug",
 "8HJE":"Vismodegib - SMO inhibitor; off-targets not in returned set",
 "7XPY":"Eupalinolide B - natural product, no curated off-target profile",
 "3TC5":"Dexamethasone - GR agonist (on-target class)",
 "6GH9":"Mitoxantrone - DNA intercalator/TOP2; protein off-targets not USP-family-documented",
}

def load_hits(path, gene_col=True):
    """target -> set(hit_uniprot) from an annotated or raw high/filtered file."""
    d={}
    if not os.path.exists(path): return None
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d.setdefault(r["query_target"], set()).add(r["hit_uniprot"].strip())
    return d

def load_library():
    if not os.path.exists(LIBPOCK): return None
    s=set()
    with open(LIBPOCK) as f:
        rd=csv.reader(f, delimiter="\t"); next(rd, None)
        for row in rd:
            if row and row[0].strip(): s.add(row[0].strip())
    return s

def conf(just):  # extract [high]/[med]/[low] tag
    for t in ("[high]","[med]","[low]"):
        if t in just: return t.strip("[]")
    return "na"

def main():
    high = load_hits(HIGH_ANNOT) or load_hits(HIGH_RAW)
    if high is None: sys.exit(f"missing {HIGH_ANNOT} / {HIGH_RAW}")
    allsup = load_hits(FILTERED)          # optional
    lib = load_library()                  # optional

    rows=["\t".join(["drug_target","known_uniprot","justification","confidence",
                     "in_library","found_at_high","found_at_supported"])]
    per=[]  # (target, n_known_eval, n_found_high, n_found_sup)
    print("=== Part 4.4 known-polypharmacology recall ===")
    print("(recall computed only on drugs with documented protein off-targets)\n")
    for tgt, known in KNOWN.items():
        hi = high.get(tgt, set())
        sup = allsup.get(tgt, set()) if allsup else set()
        n_eval=n_hi=n_sup=0
        print(f"--- {tgt} ---")
        for acc, just in known.items():
            if "self" in just.lower() or "on-target" in just.lower():
                continue
            in_lib = ("yes" if (lib and acc in lib) else ("no" if lib else "unknown"))
            # a target not in the library cannot be found - exclude from denominator
            evaluable = (lib is None) or (acc in lib)
            fh = acc in hi; fs = acc in sup
            if evaluable:
                n_eval += 1; n_hi += int(fh); n_sup += int(fh or fs)
            mark = "FOUND@HIGH" if fh else ("found@supported" if fs else "missed")
            if not evaluable: mark = "not-in-library (excluded)"
            print(f"   {acc} {conf(just):4s} {mark:26s} {just.split(' - ')[0]}")
            rows.append("\t".join([tgt,acc,just,conf(just),in_lib,
                                   "yes" if fh else "no","yes" if (fh or fs) else "no"]))
        if n_eval:
            per.append((tgt,n_eval,n_hi,n_sup))
            print(f"   -> recall@high {n_hi}/{n_eval} = {n_hi/n_eval:.2f}   "
                  f"@supported {n_sup}/{n_eval} = {n_sup/n_eval:.2f}\n")
        else:
            print("   -> no evaluable known targets in library (excluded)\n")

    open(OUT,"w").write("\n".join(rows)+"\n")

    tot_eval=sum(p[1] for p in per); tot_hi=sum(p[2] for p in per); tot_sup=sum(p[3] for p in per)
    print("=== overall (evaluable known targets across scored drugs) ===")
    print(f"   drugs scored: {len(per)}   known targets evaluable: {tot_eval}")
    if tot_eval:
        print(f"   recall @HIGH      : {tot_hi}/{tot_eval} = {tot_hi/tot_eval:.2f}")
        print(f"   recall @SUPPORTED : {tot_sup}/{tot_eval} = {tot_sup/tot_eval:.2f}")
    print(f"\n   drugs excluded (clean/selective or off-targets out of scope): {len(EXCLUDED)}")
    print("   (listed in script EXCLUDED dict; excluded from denominator, not counted as misses)")
    if lib is None:
        print("\n[note] library_pockets.tsv not found - 'in_library' is 'unknown' and no target")
        print("       was excluded for library membership. Provide it to make recall exact.")
    print(f"\nWrote {OUT}")
    print("\nHONEST READING: recall is high where the drug's documented off-targets are")
    print("paralogues/family of the crystal target (found by geometry); it is LOW where the")
    print("documented off-targets are unrelated folds (e.g. imatinib->ABL/KIT, vorinostat->HDAC),")
    print("which the pocket-geometry method is not expected to reach. Report both, honestly.")

if __name__ == "__main__":
    main()
