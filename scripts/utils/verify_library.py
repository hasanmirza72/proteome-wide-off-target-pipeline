#!/usr/bin/env python3
"""
Decide whether the 739 'very long' flags in library_20550.fasta are:
  (a) harmless genuinely-long single-chain AlphaFold models  -> library is FINE, or
  (b) real multi-chain repeats that need re-extraction.

The distinguishing test: a repeat contains its own head sequence again internally.
A genuinely long protein does not.

Usage:
    python3 verify_library.py library_20550.fasta
"""
import sys

def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        line = line.rstrip()
        if line.startswith(">"):
            if name is not None: yield name, "".join(seq)
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name is not None: yield name, "".join(seq)

def is_real_repeat(seq):
    """A concatenated-chain repeat contains an exact copy of its first 40 residues later on."""
    if len(seq) < 80: return False
    head = seq[:40]
    return seq.find(head, 1) != -1

def main(path):
    total = long_count = repeat_count = 0
    repeat_examples = []
    long_examples = []
    for name, seq in read_fasta(path):
        total += 1
        n = len(seq)
        if n > 1500:
            long_count += 1
            if is_real_repeat(seq):
                repeat_count += 1
                if len(repeat_examples) < 10:
                    repeat_examples.append((name, n))
            else:
                if len(long_examples) < 5:
                    long_examples.append((name, n))
    print(f"Total sequences         : {total}")
    print(f"Long (>1500 aa)         : {long_count}")
    print(f"  -> REAL repeats       : {repeat_count}   <-- these are the only real problem")
    print(f"  -> genuinely long OK  : {long_count - repeat_count}")
    print()
    if repeat_count == 0:
        print("✓✓ VERDICT: No real repeats. The 'very long' flags are genuine long")
        print("   single-chain AlphaFold models. library_20550.fasta is CORRECT as-is.")
        print("   You do NOT need to re-extract it. Proceed to makeblastdb.")
        if long_examples:
            print("\n   (examples of harmless long proteins:)")
            for nm,n in long_examples: print(f"     {nm}: {n} aa")
    else:
        print(f"✗ VERDICT: {repeat_count} sequences are real multi-chain repeats.")
        print("   Re-extract ONLY these with the first-chain-only function.")
        print("\n   repeat examples:")
        for nm,n in repeat_examples: print(f"     {nm}: {n} aa")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 verify_library.py library_20550.fasta"); sys.exit(1)
    main(sys.argv[1])
