#!/usr/bin/env python3
"""
Check whether a drug site is inter-subunit (quaternary).

Load the FULL biological assembly (NOT a monomer extraction). For each copy of the
ligand, report which protein chains contribute residues within the contact cutoff.
If a single ligand copy contacts >1 chain, the pocket is an interface pocket that a
monomeric AlphaFold model cannot reproduce -> legitimate special case / exclusion.

Usage:
    python check_quaternary_site.py 2Y05_assembly1.pdb RAL
"""
import sys
from Bio.PDB import MMCIFParser, PDBParser, NeighborSearch, is_aa

CUTOFF = 4.5  # A, ligand-heavy-atom to protein-heavy-atom


def load(path):
    p = MMCIFParser(QUIET=True) if path.lower().endswith((".cif", ".mmcif")) \
        else PDBParser(QUIET=True)
    return p.get_structure("s", path)


def main(path, resname):
    resname = resname.strip().upper()
    model = next(iter(load(path)))

    protein_atoms = [a for a in model.get_atoms()
                     if is_aa(a.get_parent(), standard=False)
                     and (a.element or "").strip() != "H"]
    ns = NeighborSearch(protein_atoms)

    lig_copies = [(ch.id, res) for ch in model for res in ch
                  if res.resname.strip().upper() == resname]
    if not lig_copies:
        print(f"No copies of ligand {resname} found.")
        return

    print(f"Ligand {resname}: {len(lig_copies)} copy/copies. Cutoff {CUTOFF} A.\n")
    for lig_chain, res in lig_copies:
        contacts = {}
        for atom in res:
            if (atom.element or "").strip() == "H":
                continue
            for nb in ns.search(atom.coord, CUTOFF):
                pr = nb.get_parent()
                contacts.setdefault(pr.get_parent().id, set()).add((pr.id[1], pr.resname))
        chains = sorted(contacts)
        verdict = "INTER-SUBUNIT (quaternary)" if len(chains) > 1 else "single-chain"
        print(f"  {resname} on chain {lig_chain} -> contacts chains {chains}  [{verdict}]")
        for ch in chains:
            residues = sorted(contacts[ch])
            shown = [f"{n}{rn}" for n, rn in residues][:12]
            print(f"      chain {ch}: {len(residues)} residues {shown}"
                  + (" ..." if len(residues) > 12 else ""))
    print()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python check_quaternary_site.py <full_assembly_file> <LIG_RESNAME>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
