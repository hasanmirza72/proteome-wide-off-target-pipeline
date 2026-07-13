#!/usr/bin/env python3
import csv

CSV_FILE = "p2rank_evaluation.csv"

def is_top3(rank_str):
    """Checks if the rank is 1, 2, or 3."""
    try:
        rank = int(rank_str)
        return 1 <= rank <= 3
    except (ValueError, TypeError):
        return False

def main():
    both_success = []
    holo_only_success = []
    apo_only_success = []
    failed_both = []

    try:
        with open(CSV_FILE, newline="") as f:
            reader = csv.DictReader(f)
            
            for row in reader:
                holo_pdb = row["holo"].upper()
                apo_pdb = row["apo"].upper()
                holo_rank = row["holo_hit_rank"]
                apo_rank = row["apo_hit_rank"]
                
                # Check Top-3 status
                holo_hit = is_top3(holo_rank)
                apo_hit = is_top3(apo_rank)
                
                # Format a clean string for printing
                record = f"  • {holo_pdb} (Apo: {apo_pdb}) | Holo Rank: {holo_rank}, Apo Rank: {apo_rank}"
                
                if holo_hit and apo_hit:
                    both_success.append(record)
                elif holo_hit and not apo_hit:
                    holo_only_success.append(record)
                elif apo_hit and not holo_hit:
                    apo_only_success.append(record)
                else:
                    failed_both.append(record)
                    
    except FileNotFoundError:
        print(f"❌ Error: Could not find '{CSV_FILE}'.")
        return

    # --- PRINT THE RESULTS ---
    print("=" * 80)
    print("   P2RANK TOP-3 SUCCESS ANALYSIS")
    print("=" * 80)
    
    print(f"\n✅ SUCCEEDED IN BOTH ({len(both_success)} structures)")
    print("   (The pocket is obvious whether the drug is there or not)")
    for r in both_success: print(r)
        
    print(f"\n⚠️ SUCCEEDED IN HOLO ONLY ({len(holo_only_success)} structures)")
    print("   (The pocket likely collapsed or closed up in the empty Apo state)")
    for r in holo_only_success: print(r)
        
    print(f"\n🤯 SUCCEEDED IN APO ONLY ({len(apo_only_success)} structures)")
    print("   (Anomalies: The algorithm preferred the empty shape over the drug-bound shape)")
    for r in apo_only_success: print(r)
        
    print(f"\n❌ FAILED IN BOTH ({len(failed_both)} structures)")
    print("   (These are exceptionally difficult, cryptic, or shallow pockets)")
    for r in failed_both: print(r)
    
    print("\n" + "=" * 80)

if __name__ == "__main__":
    main()
