#!/usr/bin/env python3
"""
Extract FoldDisco RAW (unfiltered) ranked hits per query from the hits/ folder.
CORRECTED for the real format: accession is inside the tid path as AF-<ACC>-F1-model_v6.pdb,
and rows are already sorted by IDF (best first) = rank order.

Produces one file: folddisco_raw_hits.tsv
  columns: query_pdb  conf  rank  hit_uniprot  idf

Run from your project folder (where hits/ lives):
    python3 extract_raw_ranks.py
Then upload folddisco_raw_hits.tsv
"""
import os, glob, csv, re

HITS_DIR = "hits"
OUT = "folddisco_raw_hits.tsv"

# accession sits between 'AF-' and '-F' :  AF-P00918-F1-model_v6.pdb -> P00918
ACC_RE = re.compile(r'AF-([A-Z0-9]+)-F\d')

rows = []
files = sorted(glob.glob(os.path.join(HITS_DIR, "*.tsv")))
if not files:
    print(f"ERROR: no .tsv files in {HITS_DIR}/ — run from the project folder."); raise SystemExit(1)

for path in files:
    base = os.path.basename(path).replace(".tsv", "")     # 1AVN_apo
    pdb, conf = base.rsplit("_", 1) if "_" in base else (base, "na")
    with open(path) as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)                              # skip header row
        # locate the tid column (the path) and idf column
        tid_col = next((i for i,c in enumerate(header) if c.lower() in ("tid","target","path")), 0)
        idf_col = next((i for i,c in enumerate(header) if c.lower() == "idf"), 1)
        rank = 0
        for parts in reader:
            if not parts or len(parts) <= tid_col: continue
            m = ACC_RE.search(parts[tid_col])
            if not m:                                      # skip anything unparseable
                continue
            acc = m.group(1)
            rank += 1
            idf = parts[idf_col] if idf_col < len(parts) else ""
            rows.append((pdb, conf, rank, acc, idf))

with open(OUT, "w") as out:
    out.write("query_pdb\tconf\trank\thit_uniprot\tidf\n")
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")

print(f"Wrote {OUT}: {len(rows)} rows from {len(files)} files in {HITS_DIR}/")

# sanity: show a couple of parsed rows
print("\nsample parsed rows:")
for r in rows[:3]:
    print("  ", r)

# THE KEY QUESTION: does imatinib (1XBB) RAW search find ABL1/KIT/PDGFRA/DDR1?
print("\n=== Does FoldDisco RAW search find imatinib's documented targets? ===")
for known, name in [("P00519","ABL1"),("P10721","KIT"),("P16234","PDGFRA"),("Q08345","DDR1")]:
    ranks = [rk for (p,c,rk,a,idf) in rows if p=="1XBB" and a==known]
    if ranks:
        print(f"  {name} ({known}): PRESENT at rank(s) {sorted(ranks)}")
    else:
        print(f"  {name} ({known}): absent from raw hits")

# also report how deep 1XBB's list goes, so we know top-k context
depth = max((rk for (p,c,rk,a,idf) in rows if p=="1XBB"), default=0)
print(f"\n  (1XBB raw hit list depth: {depth} hits per conformation)")
