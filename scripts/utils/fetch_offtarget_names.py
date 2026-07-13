import urllib.request
import json
import csv
import os

# --- Configuration ---
HITS_FILE = "offtarget_hits_high.tsv"

def main():
    if not os.path.exists(HITS_FILE):
        print(f"[!] Could not find {HITS_FILE}. Make sure you are running this from the main directory.")
        return

    # 1. Extract the unique UniProt IDs from the 3rd column
    uniprot_ids = set()
    with open(HITS_FILE, "r") as f:
        # Assuming it's tab-separated based on the .tsv extension
        reader = csv.reader(f, delimiter="\t")
        next(reader, None) # Skip the header
        for row in reader:
            if len(row) >= 3:
                # The hit_uniprot is usually the 3rd column (index 2)
                uid = row[2].strip()
                if uid and uid != "is_self": # basic cleanup
                    uniprot_ids.add(uid)

    print(f"Found {len(uniprot_ids)} unique off-target UniProt IDs.\nFetching names from UniProt API...\n")
    print(f"{'UniProt':<10} {'Gene':<15} {'Protein Name'}")
    print("-" * 80)

    # 2. Fetch the names
    for uid in sorted(uniprot_ids):
        url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
        try:
            req = urllib.request.urlopen(url)
            data = json.loads(req.read())
            
            # Try to get the Recommended Name, fallback to Submission Name
            try:
                prot_name = data['proteinDescription']['recommendedName']['fullName']['value']
            except KeyError:
                try:
                    prot_name = data['proteinDescription']['submissionNames'][0]['fullName']['value']
                except (KeyError, IndexError):
                    prot_name = "Unknown"
                    
            # Try to get the Gene Name
            try:
                gene_name = data['genes'][0]['geneName']['value']
            except (KeyError, IndexError):
                gene_name = "N/A"
                
            print(f"{uid:<10} {gene_name:<15} {prot_name}")
            
        except Exception as e:
            print(f"{uid:<10} Error fetching from UniProt (Might be an obsolete/merged ID)")

if __name__ == "__main__":
    main()
