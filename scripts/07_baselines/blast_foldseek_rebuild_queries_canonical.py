#!/usr/bin/env python3
"""
Rebuild BOTH baseline query inputs from the CANONICAL non-redundant file, so the
baseline searches exactly the same 23 targets as FoldDisco (removes the Q5SLE7
contaminant, restores P62942/FKBP1A).

Builds:
  benchmark_targets.fasta        (for BLAST)   - first-chain sequence per target
  foldseek_queries/<ACC>.pdb     (for Foldseek)- first-chain structure per target

Source of truth:
  dataset_final_nonredundant.csv  (columns include 'uniprot' and 'holo_pdb')
  dataset_monomers_pdb/           (folder of holo .pdb files)

Usage:
  python3 rebuild_queries_canonical.py
Then re-run BLAST and Foldseek on the rebuilt inputs.
"""
import csv, os

CANON = 'dataset_final_nonredundant.csv'   # <-- the CORRECT source (not apo_holo_pairs_*)
PDBDIR = 'dataset_monomers_pdb'
OUT_FASTA = 'benchmark_targets.fasta'
OUT_QDIR  = 'foldseek_queries'

AA = {'CYS':'C','ASP':'D','SER':'S','GLN':'Q','LYS':'K','ILE':'I','PRO':'P','THR':'T',
      'PHE':'F','ASN':'N','GLY':'G','HIS':'H','LEU':'L','ARG':'R','TRP':'W','ALA':'A',
      'VAL':'V','GLU':'E','TYR':'Y','MET':'M'}

def first_chain_seq(path):
    seq, first, seen = [], None, set()
    with open(path) as f:
        for line in f:
            if not line.startswith('ATOM'): continue
            if line[12:16].strip() != 'CA': continue
            ch = line[21]
            if first is None: first = ch
            elif ch != first: break
            alt = line[16].strip()
            if alt not in ('','A'): continue
            key = (ch, line[22:26].strip(), line[26].strip())
            if key in seen: continue
            seen.add(key)
            seq.append(AA.get(line[17:20].strip(),'X'))
    return "".join(seq)

def first_chain_struct(inp, outp):
    first = None
    with open(inp) as fi, open(outp,'w') as fo:
        for line in fi:
            if line.startswith(('ATOM','HETATM','TER','ANISOU')):
                ch = line[21]
                if first is None: first = ch
                if ch != first: continue
                fo.write(line)
            elif line.startswith(('HEADER','CRYST','MODEL','ENDMDL')):
                fo.write(line)
        fo.write("END\n")

def find_pdb(pdb):
    for cand in (f"{pdb}.pdb", f"{pdb.lower()}.pdb", f"{pdb.upper()}.pdb"):
        p = os.path.join(PDBDIR, cand)
        if os.path.exists(p): return p
    return None

os.makedirs(OUT_QDIR, exist_ok=True)
# clear any old (possibly contaminated) query files
for old in os.listdir(OUT_QDIR):
    if old.endswith('.pdb'): os.remove(os.path.join(OUT_QDIR, old))

targets, missing = [], []
with open(CANON, newline='') as f:
    r = csv.DictReader(f)
    assert 'uniprot' in r.fieldnames and 'holo_pdb' in r.fieldnames, \
        f"Expected 'uniprot' and 'holo_pdb' columns; found {r.fieldnames}"
    for row in r:
        acc = row['uniprot'].strip()
        pdb = row['holo_pdb'].strip()
        src = find_pdb(pdb)
        if not src:
            missing.append((acc,pdb)); continue
        targets.append((acc, src))

with open(OUT_FASTA,'w') as fa:
    for acc, src in targets:
        seq = first_chain_seq(src)
        if seq:
            fa.write(f">{acc}\n{seq}\n")
        first_chain_struct(src, os.path.join(OUT_QDIR, f"{acc}.pdb"))

accs = sorted(a for a,_ in targets)
print(f"Rebuilt from {CANON}")
print(f"  targets written : {len(targets)}")
print(f"  FASTA           : {OUT_FASTA}")
print(f"  query structures: {OUT_QDIR}/")
print(f"  Q5SLE7 present? : {'YES (BAD)' if 'Q5SLE7' in accs else 'no (good)'}")
print(f"  P62942 present? : {'yes (good)' if 'P62942' in accs else 'NO (BAD)'}")
if missing:
    print("  MISSING holo PDBs:")
    for a,p in missing: print(f"    {a} -> {p}")
print("\nNext:")
print("  bash scripts/run_blast.sh        # rebuild db + search on corrected 23")
print("  bash scripts/run_foldseek.sh     # point at foldseek_queries/ (corrected 23)")
print("  # then re-verify: cut -f1 blast_hits.tsv | sort -u | wc -l   -> 23, incl P62942 not Q5SLE7")
