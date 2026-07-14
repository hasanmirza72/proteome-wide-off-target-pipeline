#!/usr/bin/env python3
"""
sensitivity_analysis.py — threshold sensitivity for the high-confidence tier.

Re-applies the high-confidence filter to the existing annotated hits at a range of
threshold values, one axis at a time, and reports how the high-confidence count
changes. This shows whether the 58 high-confidence candidates are robust to the exact
threshold choices, which is the sensitivity analysis requested in review.

It reads offtarget_hits_annotated.tsv (the full annotated hit list, before tiering) and
re-derives the high-confidence count. It does NOT re-run the pipeline; it re-filters the
existing per-hit values, so the baseline count reproduces the reported number exactly.

Run from your project root (where offtarget_hits_annotated.tsv lives):
    python3 sensitivity_analysis.py
"""
import csv

INFILE = "offtarget_hits_annotated.tsv"

# Baseline high-confidence thresholds (must match offtarget_overlap_filter.py exactly)
BASE = dict(nodes=4, overlap=3, frac=0.30, rmsd=2.0, idf=3.0)

def to_float(x, default):
    try: return float(x)
    except (TypeError, ValueError): return default

def load_hits(path):
    hits = []
    with open(path) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r.get("is_self", "0") in ("1", "True", "true"):
                continue  # exclude self-hits, as the pipeline does
            hits.append(dict(
                nodes=to_float(r.get("node_count"), 0),
                overlap=to_float(r.get("overlap_count"), 0),
                frac=to_float(r.get("overlap_frac"), 0),
                rmsd=to_float(r.get("min_rmsd"), 99.0),
                idf=to_float(r.get("idf"), 0),
            ))
    return hits

def count_high(hits, th):
    n = 0
    for h in hits:
        if (h["nodes"] >= th["nodes"] and h["overlap"] >= th["overlap"]
                and h["frac"] >= th["frac"] and h["rmsd"] <= th["rmsd"]
                and h["idf"] >= th["idf"]):
            n += 1
    return n

def main():
    hits = load_hits(INFILE)
    base_n = count_high(hits, BASE)
    print("=" * 60)
    print(f"BASELINE high-confidence count: {base_n}")
    print(f"(thresholds: nodes>={BASE['nodes']}, overlap>={BASE['overlap']}, "
          f"frac>={BASE['frac']}, rmsd<={BASE['rmsd']}, idf>={BASE['idf']})")
    print("=" * 60)

    # vary each threshold one at a time, holding the others at baseline
    sweeps = {
        "idf":     [2.0, 2.5, 3.0, 3.5, 4.0, 5.0],
        "nodes":   [3, 4, 5, 6],
        "overlap": [2, 3, 4],
        "frac":    [0.20, 0.30, 0.40, 0.50],
        "rmsd":    [1.5, 2.0, 2.5, 3.0],
    }
    for axis, values in sweeps.items():
        print(f"\n--- varying {axis} (others at baseline) ---")
        for v in values:
            th = dict(BASE); th[axis] = v
            n = count_high(hits, th)
            delta = n - base_n
            marker = "  <- baseline" if v == BASE[axis] else ""
            print(f"   {axis} = {v:<5} -> {n:3d} high-confidence  "
                  f"({'+' if delta>=0 else ''}{delta} vs baseline){marker}")

    print("\n" + "=" * 60)
    print("Interpretation: if the count changes only modestly across nearby")
    print("threshold values, the high-confidence set is robust to the exact")
    print("choices. Report the baseline count and the range in the thesis.")

if __name__ == "__main__":
    main()
