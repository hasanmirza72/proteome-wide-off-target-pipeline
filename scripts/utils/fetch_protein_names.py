import csv
import json
import urllib.request
import sys

manifest = "queries_manifest.tsv"

print(f"{'Target':<8} {'UniProt':<10} {'Gene':<10} {'Protein Name'}")
print("-" * 80)

with open(manifest) as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        pdb = row.get('target_id', '').strip()
        uid = row.get('self_uniprot', '').strip()
        
        if not uid or pdb.startswith("#"):
            continue
            
        url = f"https://rest.uniprot.org/uniprotkb/{uid}.json"
        
        try:
            req = urllib.request.urlopen(url)
            data = json.loads(req.read())
            
            # Extract Protein Name
            try:
                prot_name = data['proteinDescription']['recommendedName']['fullName']['value']
            except KeyError:
                try:
                    prot_name = data['proteinDescription']['submissionNames'][0]['fullName']['value']
                except (KeyError, IndexError):
                    prot_name = "Unknown"
                    
            # Extract Gene Name
            try:
                gene_name = data['genes'][0]['geneName']['value']
            except (KeyError, IndexError):
                gene_name = "N/A"
                
            print(f"{pdb:<8} {uid:<10} {gene_name:<10} {prot_name}")
            
        except Exception as e:
            print(f"{pdb:<8} {uid:<10} Error fetching from UniProt")

