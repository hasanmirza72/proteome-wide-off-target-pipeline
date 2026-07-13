import os
import re
import csv

# --- Configuration ---
PDB_DIR = "alphafold_structures_pdb_only"
LIB_FILE = "library_pockets.tsv"

# Regex to explicitly extract ONLY F1 models
af_regex = re.compile(r"AF-([A-Z0-9]+)-F1-model")

def main():
    print("Auditing F1 structures against library_pockets.tsv...\n")
    
    # 1. Get all F1 Accessions
    f1_ids = set()
    if not os.path.exists(PDB_DIR):
        print(f"[!] Directory not found: {PDB_DIR}")
        return
        
    for fname in os.listdir(PDB_DIR):
        match = af_regex.search(fname)
        if match:
            f1_ids.add(match.group(1))
            
    # 2. Get Accessions in the Flattened Library
    lib_ids = set()
    if os.path.exists(LIB_FILE):
        with open(LIB_FILE, "r") as f:
            reader = csv.reader(f, delimiter="\t")
            next(reader, None)  # skip header
            for row in reader:
                if row and row[0].strip():
                    lib_ids.add(row[0].strip())
    else:
        print(f"[!] File not found: {LIB_FILE}")
        return
                    
    # 3. Calculate the Difference
    missing_ids = f1_ids - lib_ids
    
    print("=== POCKET FILTER AUDIT ===")
    print(f"Total F1 models (Starting set):        {len(f1_ids):,}")
    print(f"Proteins with >= 1 pocket (Library):   {len(lib_ids):,}")
    print(f"Proteins with ZERO pockets (Filtered): {len(missing_ids):,}\n")
    
    # 4. Save the proof
    out_file = "zero_pocket_proteins.txt"
    with open(out_file, "w") as f:
        for uid in sorted(missing_ids):
            f.write(f"{uid}\n")
            
    print(f"[ok] Saved list of {len(missing_ids)} zero-pocket proteins to '{out_file}'")
    print("     (Biologically, these are likely intrinsically disordered or lack cavities.)")

if __name__ == "__main__":
    main()
