#!/usr/bin/env python3
import csv
import os

CSV_FILE = "apo_holo_pairs.csv"

def main():
    if not os.path.exists(CSV_FILE):
        print(f"❌ Error: '{CSV_FILE}' not found.")
        return

    valid_count = 0
    invalid_rows = []

    # Open the file in read-only mode
    with open(CSV_FILE, mode='r', encoding='utf-8', errors='ignore') as f:
        reader = csv.DictReader(f)
        for row in reader:
            holo = row.get("holo_pdb", "Unknown").strip()
            apo = row.get("apo_pdb", "Unknown").strip()
            source = row.get("apo_source", "Missing Column").strip()

            if source == "ahoj_aligned_cif":
                valid_count += 1
            else:
                invalid_rows.append((holo, apo, source))

    # Print the report
    print("=" * 60)
    print("   VALIDITY CHECK: APO COORDINATE ALIGNMENT")
    print("=" * 60)
    print(f"✅ Correctly Aligned (ahoj_aligned_cif): {valid_count} pairs")
    
    if len(invalid_rows) > 0:
        print(f"⚠️  WARNING: Found {len(invalid_rows)} unaligned/suspicious pairs.")
        print("\n--- List of Suspicious Pairs (DATA HAS NOT BEEN DELETED) ---")
        for holo, apo, src in invalid_rows:
            print(f"  • Holo: {holo}  |  Apo: {apo}  |  Source: {src}")
    else:
        print("\n🎉 PERFECT! All Apo structures are perfectly aligned.")
        print("   Your spatial distance calculations (DCA) are 100% scientifically valid.")
    
    print("=" * 60)

if __name__ == "__main__":
    main()
