import os
import re

# --- Configuration ---
PDB_DIR = "alphafold_structures_pdb_only"
SWISSPROT_FILE = "swissprot_human.txt"

# Regex to extract just the UniProt ID from F1 files
af_regex = re.compile(r"AF-([A-Z0-9]+)-F1-model")

def main():
    # 1. Load the official Swiss-Prot IDs
    sp_ids = set()
    with open(SWISSPROT_FILE, "r") as f:
        for line in f:
            if line.strip():
                sp_ids.add(line.strip())

    # 2. Extract the AlphaFold F1 IDs from your directory
    af_ids = set()
    for fname in os.listdir(PDB_DIR):
        match = af_regex.search(fname)
        if match:
            af_ids.add(match.group(1))

    # 3. Calculate the intersection (Proteins that are in BOTH sets)
    reviewed_af_models = af_ids.intersection(sp_ids)
    
    # 4. Math
    total_sp = len(sp_ids)
    total_af = len(af_ids)
    overlap = len(reviewed_af_models)
    coverage = (overlap / total_sp) * 100

    print("=== SWISS-PROT COVERAGE ANALYSIS ===")
    print(f"Official Swiss-Prot Human Proteins: {total_sp:,}")
    print(f"Your AlphaFold F1 Models:           {total_af:,}")
    print(f"AlphaFold Models that are Reviewed: {overlap:,}")
    print(f"\nFinal Coverage Percentage:          {coverage:.2f}%")

if __name__ == "__main__":
    main()
