#!/usr/bin/env python3
"""
Definitive target-count audit. Answers three questions:
  1. How many non-redundant targets does EACH pipeline file actually contain?
  2. Is Q5SLE7 present or absent in each?
  3. Which is the single consistent N (should be 22 or 23) to use in the thesis?

Run from inside your project folder (where all these files live):
    python3 audit_target_count.py

It reads whatever files exist and skips the rest, so it is safe to run as-is.
"""
import csv, os, re, gzip

SUSPECT = "Q5SLE7"  # the target we are checking

# (filename, how-to-extract-accessions). Each extractor returns a set of accessions.
def col_from_csv(path, candidates):
    """Read a CSV/TSV and pull accessions from the first matching column name."""
    if not os.path.exists(path): return None
    delim = '\t' if path.endswith(('.tsv','.txt')) else ','
    accs = set()
    with open(path, newline='') as f:
        r = csv.DictReader(f, delimiter=delim)
        if not r.fieldnames: return set()
        col = next((c for c in candidates if c in r.fieldnames), None)
        if col is None:
            # fall back: first column
            col = r.fieldnames[0]
        for row in r:
            v = (row.get(col) or '').strip()
            if v: accs.add(v.split('_')[0].split('|')[-1] if '|' in v else v.split('_')[0])
    return accs

def first_field_lines(path):
    """Pull the first whitespace/tab field from each non-header line."""
    if not os.path.exists(path): return None
    accs = set()
    with open(path) as f:
        for i,line in enumerate(f):
            if i==0 and not re.match(r'^[OPQ][0-9]', line):  # skip header
                continue
            tok = re.split(r'[\t, ]', line.strip())[0]
            tok = tok.split('_')[0]
            if tok: accs.add(tok)
    return accs

ACC_RE = re.compile(r'\b[OPQ][0-9][A-Z0-9]{3}[0-9]\b')  # UniProt accession pattern
def grep_accessions(path):
    """Extract anything that looks like a UniProt accession from any text file."""
    if not os.path.exists(path): return None
    accs = set()
    opener = gzip.open if path.endswith('.gz') else open
    try:
        with opener(path, 'rt', errors='ignore') as f:
            for line in f:
                accs.update(ACC_RE.findall(line))
    except Exception:
        return None
    return accs

checks = [
    ("dataset_final_nonredundant.csv", lambda: col_from_csv("dataset_final_nonredundant.csv", ["uniprot","accession","query"])),
    ("apo_holo_pairs_nonredundant.csv", lambda: col_from_csv("apo_holo_pairs_nonredundant.csv", ["uniprot","accession","query"])),
    ("verified_alphafold_targets.csv", lambda: col_from_csv("verified_alphafold_targets.csv", ["uniprot","accession","query"])),
    ("queries_manifest.tsv", lambda: grep_accessions("queries_manifest.tsv")),
    ("self_recovery_summary.tsv", lambda: grep_accessions("self_recovery_summary.tsv")),
    ("phase5_recall.tsv", lambda: grep_accessions("phase5_recall.tsv")),
    ("offtarget_hits_annotated.tsv", lambda: grep_accessions("offtarget_hits_annotated.tsv")),
    ("benchmark_targets.fasta", lambda: grep_accessions("benchmark_targets.fasta")),
    ("library_20550.fasta (is target IN library?)", lambda: grep_accessions("library_20550.fasta")),
    ("blast_hits.tsv (queries)", lambda: {ln.split('\t')[0] for ln in open("blast_hits.tsv")} if os.path.exists("blast_hits.tsv") else None),
    ("foldseek_hits.tsv (queries)", lambda: {ln.split('\t')[0] for ln in open("foldseek_hits.tsv")} if os.path.exists("foldseek_hits.tsv") else None),
]

print(f"{'FILE':<48} {'N_targets':>10} {SUSPECT:>10}")
print("="*72)
libset = None
results = {}
for name, fn in checks:
    try:
        accs = fn()
    except Exception as e:
        accs = None
    if accs is None:
        print(f"{name:<48} {'(missing)':>10}")
        continue
    results[name] = accs
    has = "PRESENT" if SUSPECT in accs else "absent"
    # for the query files, count only plausible targets (their own query sets are small);
    # for big files, this is just the accession universe, shown for the suspect check
    label = len(accs)
    print(f"{name:<48} {label:>10} {has:>10}")
    if name.startswith("library_20550"):
        libset = accs

print("="*72)
# focused verdict on the query sets (these define the real N)
for key in ["blast_hits.tsv (queries)","foldseek_hits.tsv (queries)",
            "self_recovery_summary.tsv","offtarget_hits_annotated.tsv"]:
    if key in results:
        s = results[key]
        print(f"{key:<48} distinct queries/targets = {len(s)}")

if libset is not None:
    print()
    print(f"Is {SUSPECT} in the 20,550 library FASTA? ",
          "YES" if SUSPECT in libset else "NO  <-- explains why it never self-recovers")

print("\nVERDICT GUIDANCE:")
print("  - If the query/self-recovery/offtarget files all show 22 and Q5SLE7 is 'absent',")
print("    then N=22 is the true count and Q5SLE7 was never searched by any method.")
print("  - If any FoldDisco file DOES contain Q5SLE7, tell that to your assistant;")
print("    it would mean the library FASTA is missing it but FoldDisco saw it (a fixable")
print("    library-build inconsistency), which is a different situation.")
