#!/usr/bin/env python3
"""
Strip ligands/waters from holo PDBs to produce 'computationally apo' inputs for P2Rank.

Upgrades over a plain 'reject everything HETATM' filter:
  1. First model only  - NMR/multi-model entries otherwise give P2Rank several
     overlaid copies of the protein.
  2. Single altloc      - keeps one conformer per atom (altloc ' ' or 'A'); leaving
     all altlocs in place puts duplicate atoms on the surface.
  3. Keep modified polymer residues (MSE, phospho-serine, etc.) - these are part of
     the chain. Deleting them (they are HETATM records) punches gaps into the
     backbone right where the pocket often is. We keep them so the surface stays
     continuous, while still removing the drug, buffers, ions and water.

Note: P2Rank predicts from the protein surface only, so stripping does not change
predictions relative to the untouched holo - it just makes inputs uniform. Keep your
ORIGINAL holo files: the ligand in them is your ground-truth pocket for scoring.
"""

import os
from Bio.PDB import PDBParser, PDBIO, Select

INPUT_DIR = "dataset_monomers_pdb"
OUTPUT_DIR = "dataset_monomers_cleaned"

# Modified amino acids / nucleotides that are genuinely part of the polymer chain.
KEEP_MODIFIED = {
    "MSE", "SEP", "TPO", "PTR", "CSO", "CSD", "KCX", "LLP", "MLY", "M3L",
    "HYP", "PCA", "CAS", "CME", "OCS", "CSX", "CGU", "SNN", "MLZ", "ALY",
    "FME", "AYA", "PHD", "SAC", "CSS", "CSW", "NEP", "HIC", "TYS",
}


class CleanProtein(Select):
    def accept_model(self, model):
        return model.id == 0                      # first model only

    def accept_residue(self, residue):
        hetflag = residue.id[0]
        if hetflag == " ":                        # standard residue
            return 1
        if hetflag == "W":                        # water
            return 0
        return 1 if residue.resname.strip() in KEEP_MODIFIED else 0  # keep polymer mods

    def accept_atom(self, atom):
        alt = atom.get_altloc()
        return alt in (" ", "A")                  # single conformer


def main():
    if not os.path.exists(INPUT_DIR):
        print(f"Error: '{INPUT_DIR}' not found."); return
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    parser = PDBParser(QUIET=True)
    io = PDBIO()
    pdb_files = [f for f in os.listdir(INPUT_DIR) if f.endswith(".pdb")]

    print("=" * 80)
    print(f"Stripping ligands from {len(pdb_files)} holo structures")
    print("=" * 80)

    ok = 0
    for name in pdb_files:
        inp = os.path.join(INPUT_DIR, name)
        out = os.path.join(OUTPUT_DIR, name.replace(".pdb", "_clean.pdb"))
        try:
            structure = parser.get_structure("p", inp)
            io.set_structure(structure)
            io.save(out, CleanProtein())
            ok += 1
            print(f"  cleaned {name} -> {os.path.basename(out)}")
        except Exception as e:
            print(f"  [!] failed {name}: {e}")

    print("-" * 80)
    print(f"Done: {ok}/{len(pdb_files)} cleaned -> ./{OUTPUT_DIR}/")


if __name__ == "__main__":
    main()
