#!/usr/bin/env python3
"""
Audit a FASTA you already built, to catch the multi-chain / repeat problem
BEFORE running BLAST. Flags any sequence that looks like a concatenated repeat
or is suspiciously long.

Usage:
    python3 audit_fasta.py library_20550.fasta
    python3 audit_fasta.py benchmark_targets.fasta
"""
import sys

def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None:
                yield name, "".join(seq)
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name is not None:
        yield name, "".join(seq)

def looks_repeated(seq):
    """Heuristic: does the sequence contain an internal exact repeat of its own head?"""
    n = len(seq)
    if n < 60:
        return False, None
    head = seq[:40]                      # first 40 residues
    second = seq.find(head, 1)           # does that 40-mer recur later?
    if second != -1:
        # estimate copy number from the spacing
        period = second
        copies = round(n / period) if period else 1
        return True, copies
    return False, None

def main(path):
    total = 0
    flagged = 0
    lengths = []
    print(f"Auditing {path}\n" + "="*60)
    for name, seq in read_fasta(path):
        total += 1
        lengths.append(len(seq))
        rep, copies = looks_repeated(seq)
        xfrac = seq.count("X") / max(1, len(seq))
        problems = []
        if rep:
            problems.append(f"REPEAT (~{copies} copies)")
        if xfrac > 0.05:
            problems.append(f"{xfrac*100:.0f}% X (unknown residues)")
        if len(seq) > 1500:
            problems.append(f"very long ({len(seq)} aa)")
        if problems:
            flagged += 1
            print(f"  ⚠ {name}: len={len(seq)}  ->  {', '.join(problems)}")
    lengths.sort()
    med = lengths[len(lengths)//2] if lengths else 0
    print("="*60)
    print(f"sequences: {total}   flagged: {flagged}")
    print(f"length  min={min(lengths) if lengths else 0}  median={med}  max={max(lengths) if lengths else 0}")
    if flagged == 0:
        print("✓ No repeat/length anomalies detected.")
    else:
        print("✗ Some sequences look multi-chain or repeated — re-extract with first-chain-only.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 audit_fasta.py <file.fasta>"); sys.exit(1)
    main(sys.argv[1])
