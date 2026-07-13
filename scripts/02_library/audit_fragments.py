import os
import re
from collections import defaultdict

# --- Configuration ---
PDB_DIR = "alphafold_structures_pdb_only"

# Regex to extract the Accession ID and the Fragment Number
# e.g., AF-Q8WZ42-F2-model_v6.pdb -> Group 1: Q8WZ42, Group 2: 2
af_regex = re.compile(r"AF-([A-Z0-9]+)-F(\d+)-model")

def main():
    if not os.path.exists(PDB_DIR):
        print(f"[!] Directory not found: {PDB_DIR}")
        return

    print("Scanning AlphaFold directory for multi-fragment proteins...\n")
    
    # Dictionary to hold the fragments for each accession
    # Format: {'Q8WZ42': [1, 2, 3, 4, ...]}
    accession_fragments = defaultdict(list)
    
    total_files = 0
    
    for fname in os.listdir(PDB_DIR):
        if not fname.endswith(".pdb"):
            continue
            
        total_files += 1
        match = af_regex.search(fname)
        if match:
            acc = match.group(1)
            frag_num = int(match.group(2))
            accession_fragments[acc].append(frag_num)

    # --- Analyze the Results ---
    multi_frag_proteins = {}
    total_f1 = 0
    total_extra_frags = 0

    for acc, frags in accession_fragments.items():
        if 1 in frags:
            total_f1 += 1
        
        # If a protein has more than 1 fragment (i.e., it has F2, F3...)
        if len(frags) > 1:
            multi_frag_proteins[acc] = sorted(frags)
            total_extra_frags += (len(frags) - 1)

    # --- Print the Summary ---
    print("=== FRAGMENT AUDIT SUMMARY ===")
    print(f"Total raw .pdb files scanned:       {total_files}")
    print(f"Total unique proteins (F1 found):   {total_f1}")
    print(f"Total 'extra' fragments (F2+):      {total_extra_frags}")
    
    math_check = total_f1 + total_extra_frags
    print(f"\nMath Check: {total_f1} + {total_extra_frags} = {math_check} (Should match {total_files})")
    
    print(f"\nNumber of giant proteins split into multiple files: {len(multi_frag_proteins)}")

    # --- Save the proof ---
    out_file = "multi_fragment_proteins.txt"
    with open(out_file, "w") as f:
        f.write("Accession\tFragments\n")
        # Sort by the number of fragments (largest proteins first!)
        for acc, frags in sorted(multi_frag_proteins.items(), key=lambda item: len(item[1]), reverse=True):
            f.write(f"{acc}\t{frags}\n")

    print(f"\n[ok] Saved the detailed list of giant proteins to '{out_file}'")

if __name__ == "__main__":
    main()
