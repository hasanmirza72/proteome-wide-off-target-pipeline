#!/usr/bin/env python3
"""
check_manifest.py -- pre-flight validation of queries_manifest.tsv.
Run BEFORE launching the gate. Reports, per target: field count, holo/apo
residue counts, chains used, and any structural flags that matter for a
single-chain AlphaFold index.
"""
import csv
import re
import sys

MANIFEST = sys.argv[1] if len(sys.argv) > 1 else "queries_manifest.tsv"
TINY = 5   # motifs smaller than this are low-specificity

def chains(res_field):
    cs = []
    for tok in res_field.split(","):
        m = re.match(r"^([A-Za-z]+)", tok.strip())
        if m and m.group(1) not in cs:
            cs.append(m.group(1))
    return cs

def n_res(res_field):
    return len([t for t in res_field.split(",") if t.strip()])

rows, bad_fields = [], []
with open(MANIFEST) as f:
    header = f.readline()
    for i, line in enumerate(f, start=2):
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) != 6:
            bad_fields.append((i, len(parts), line[:60]))
            continue
        target, uni, hpdb, hres, apdb, ares = parts
        rows.append(dict(target=target, uni=uni, hpdb=hpdb, hres=hres, apdb=apdb, ares=ares))

print(f"{'target':<7}{'uniprot':<9}{'holo_n':>7}{'holo_ch':>9}"
      f"{'apo_n':>7}{'apo_ch':>9}  flags")
print("-" * 70)

n_multi = n_switch = n_empty_apo = n_tiny = 0
for r in rows:
    hc, ac = chains(r["hres"]), chains(r["ares"])
    hn, an = n_res(r["hres"]), n_res(r["ares"])
    flags = []
    if len(hc) > 1:
        flags.append("MULTI_CHAIN_HOLO"); n_multi += 1
    if len(ac) > 1:
        flags.append("MULTI_CHAIN_APO")
    if an == 0:
        flags.append("EMPTY_APO"); n_empty_apo += 1
    elif set(hc) != set(ac):
        flags.append("chain-switch"); n_switch += 1
    if hn < TINY:
        flags.append("tiny-motif"); n_tiny += 1
    print(f"{r['target']:<7}{r['uni']:<9}{hn:>7}{','.join(hc):>9}"
          f"{an:>7}{','.join(ac) or '-':>9}  {' '.join(flags)}")

print("-" * 70)
print(f"targets: {len(rows)} | bad-field rows: {len(bad_fields)} | "
      f"multi-chain holo: {n_multi} | chain-switch: {n_switch} | "
      f"empty-apo: {n_empty_apo} | tiny-motif: {n_tiny}")
for (ln, nf, txt) in bad_fields:
    print(f"  !! line {ln}: {nf} fields (need 6): {txt}...")
print("\nNote: chain-switch is expected (holo/apo are different PDBs). "
      "MULTI_CHAIN_* motifs cannot fully match a single-chain AF model.")
