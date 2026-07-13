#!/usr/bin/env python3
"""
build_library_pockets.py  (roadmap Part 1.5)
--------------------------------------------
Flattens the per-protein P2Rank *_predictions.csv files (run over the AlphaFold
proteome) into ONE UniProt-keyed table used by the Part 3.3 overlap filter.

Output: library_pockets.tsv, one row per predicted pocket:
    uniprot  pocket_rank  probability  score  center_x  center_y  center_z  residues

`residues` is a comma-joined, FoldDisco-style list (chain+number, e.g. A123,A125)
so it compares directly against FoldDisco's matching_residues field.

The UniProt accession is extracted from each prediction filename, which for the
AFDB looks like  AF-<ACCESSION>-F1-model_v4.pdb.gz_predictions.csv.
"""

import csv
import glob
import os
import re
import sys

# --- Configuration -----------------------------------------------------------
P2RANK_DIR = "./alphafold_predictions"   # dir with *_predictions.csv
OUT_TSV    = "library_pockets.tsv"
# -----------------------------------------------------------------------------

# Official UniProt accession pattern (6 or 10 chars).
ACC_RE = re.compile(r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}")


def uniprot_from_name(fname):
    m = ACC_RE.search(fname.upper())
    return m.group(0) if m else None


def norm_residues(residue_ids_field):
    """P2Rank 'A_123 A_124' -> 'A123,A124' (skip malformed tokens)."""
    out = []
    for tok in residue_ids_field.split():
        tok = tok.strip()
        if "_" in tok:
            chain, _, num = tok.partition("_")
            if chain and num:
                out.append(f"{chain}{num}")
    return ",".join(out)


def get(row, *names, default=""):
    """Fetch a column tolerant of leading/trailing spaces in P2Rank headers."""
    for n in names:
        if n in row and row[n] is not None:
            return row[n].strip()
    return default


def main():
    files = sorted(glob.glob(os.path.join(P2RANK_DIR, "*_predictions.csv")))
    if not files:
        sys.exit(f"No *_predictions.csv under {P2RANK_DIR}")
    print(f"Flattening {len(files)} prediction files...")

    out = ["\t".join(["uniprot", "pocket_rank", "probability", "score",
                      "center_x", "center_y", "center_z", "residues"])]
    n_pockets = 0
    n_no_acc = 0
    n_files_with_pockets = 0

    for path in files:
        acc = uniprot_from_name(os.path.basename(path))
        if not acc:
            n_no_acc += 1
            continue
        with open(path, newline="") as f:
            reader = csv.DictReader(f)
            # strip whitespace from header names once
            if reader.fieldnames:
                reader.fieldnames = [c.strip() for c in reader.fieldnames]
            had = False
            for row in reader:
                residues = norm_residues(get(row, "residue_ids"))
                if not residues:
                    continue
                had = True
                n_pockets += 1
                out.append("\t".join([
                    acc,
                    get(row, "rank"),
                    get(row, "probability"),
                    get(row, "score"),
                    get(row, "center_x"),
                    get(row, "center_y"),
                    get(row, "center_z"),
                    residues,
                ]))
            if had:
                n_files_with_pockets += 1

    with open(OUT_TSV, "w") as f:
        f.write("\n".join(out) + "\n")

    print(f"Done. {n_pockets} pockets from {n_files_with_pockets} proteins -> {OUT_TSV}")
    if n_no_acc:
        print(f"  ({n_no_acc} files had no UniProt accession in their name — check naming)")


if __name__ == "__main__":
    main()
