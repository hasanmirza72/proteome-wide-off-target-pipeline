import os
import glob
import sys

# Standard 3-letter to 1-letter amino acid code mapping
# Corrected sequence extractor: FIRST CHAIN ONLY, with insertion-code handling.

AA_MAPPING = {
    'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
    'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
}

def extract_sequence_from_pdb(pdb_path):
    """
    Extract the 1-letter sequence from the FIRST chain only.
    - reads CA atoms (one per residue)
    - stops at the first chain break (so multi-chain crystals give one copy)
    - de-duplicates on (chain, resseq, icode) so altloc/insertion codes
      never double-count a residue
    - ignores HETATM (only ATOM records), so ligands/waters are excluded
    """
    sequence = []
    first_chain = None
    seen = set()
    with open(pdb_path, 'r') as f:
        for line in f:
            if not line.startswith('ATOM'):
                continue
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue
            chain_id = line[21]              # column 22 = chain identifier
            resseq   = line[22:26].strip()   # residue sequence number
            icode    = line[26].strip()      # insertion code
            altloc   = line[16].strip()      # alternate location indicator

            # lock onto the first chain we encounter; stop when it changes
            if first_chain is None:
                first_chain = chain_id
            elif chain_id != first_chain:
                break

            # skip alternate conformations beyond the first (keep '' or 'A')
            if altloc not in ('', 'A'):
                continue

            key = (chain_id, resseq, icode)
            if key in seen:
                continue
            seen.add(key)

            res_name = line[17:20].strip()
            sequence.append(AA_MAPPING.get(res_name, 'X'))
    return "".join(sequence)


# quick self-test if run directly on one file:
if __name__ == "__main__":
    import sys
    if len(sys.argv) == 2:
        seq = extract_sequence_from_pdb(sys.argv[1])
        print(f"length: {len(seq)}")
        print(seq)
def build_fasta_from_folder(pdb_folder, output_fasta_path):
    """Iterates through your folder of PDBs and writes a consolidated FASTA file."""
    pdb_files = glob.glob(os.path.join(pdb_folder, '*.pdb'))
    
    if not pdb_files:
        print(f"Error: No .pdb files found in {pdb_folder}")
        return
    
    print(f"Processing {len(pdb_files)} PDB files from {pdb_folder}...")
    
    with open(output_fasta_path, 'w') as fasta_out:
        for count, pdb_path in enumerate(pdb_files, 1):
            # Clean extraction: A0A024R1R8.pdb -> A0A024R1R8
            uniprot_id = os.path.splitext(os.path.basename(pdb_path))[0]
            sequence = extract_sequence_from_pdb(pdb_path)
            
            if sequence:
                fasta_out.write(f">{uniprot_id}\n{sequence}\n")
            
            if count % 5000 == 0:
                print(f"  Processed {count} files...")
                
    print(f"Successfully created: {output_fasta_path}")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python pdb_to_fasta.py <input_pdb_folder> <output_fasta_name.fasta>")
        sys.exit(1)
        
    input_folder = sys.argv[1]
    output_fasta = sys.argv[2]
    build_fasta_from_folder(input_folder, output_fasta)
