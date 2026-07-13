#!/usr/bin/env python3
"""
build_redundant_manifest.py  (roadmap Part 2.6 robustness sub-study)
-------------------------------------------------------------------
Builds holo pocket motifs for the EXTRA redundant co-crystal entries (same protein,
different drug) so their off-target hit lists can be compared against the
non-redundant representative. Only holo motifs are needed for the robustness study.

Reuses the exact, tested motif logic from build_queries_manifest.py (byres 4.0 A,
single best-bound ligand copy). Writes redundant_manifest.tsv in the same 6-column
format run_offtarget.sh reads (apo columns left empty -> apo query skipped).

The 8 entries (holo PDB in dataset_monomers_pdb/, verified on disk):
  CA2   (P00918): 1OQ5/CEL, 2AW1/COX, 2POU/I7A, 1CIL/ETS   (+ 1AVN/HSM non-redundant)
  PDE5A (O76074): 1XOZ/CIA, 1UHO/VDN                        (+ 6L6E/E6L)
  PXR   (O75469): 7AXK/EST, 6HJ2/P06                        (+ 7AXA/CL6)
"""
import os, sys
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, is_aa

HOLO_DIR = "dataset_monomers_pdb"
OUT_TSV  = "redundant_manifest.tsv"
CUTOFF   = 4.0
MIN_LIG_CONTACTS = 1

# (target_id, uniprot, ligand, candidate_source)  -- extra redundant entries only
ENTRIES = [
    ("1OQ5", "P00918", "CEL", "drug-matched"),
    ("2AW1", "P00918", "COX", "drug-matched"),
    ("2POU", "P00918", "I7A", "drug-matched"),
    ("1CIL", "P00918", "ETS", "drug-matched"),
    ("1XOZ", "O76074", "CIA", "drug-matched"),
    ("1UHO", "O76074", "VDN", "any-ligand-pocket"),
    ("7AXK", "O75469", "EST", "drug-matched"),
    ("6HJ2", "O75469", "P06", "drug-matched"),
]


def load_first_model(path):
    low = path.lower(); stem = low[:-3] if low.endswith(".gz") else low
    parser = MMCIFParser(QUIET=True) if stem.endswith((".cif", ".mmcif")) else PDBParser(QUIET=True)
    opener = open  # holo files here are plain .pdb
    with opener(path, "rt") as f:
        return next(iter(parser.get_structure("s", f)))


def heavy_coords(residue):
    cs = [a.coord for a in residue if (a.element or "").strip().upper() != "H"]
    return np.array(cs, dtype=float) if cs else np.empty((0, 3))


def res_label(chain_id, res):
    _, seq, icode = res.id; icode = (icode or "").strip()
    return f"{chain_id}{seq}{icode}", (chain_id, int(seq), icode)


def protein_residues(model):
    out = []
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=False):
                continue
            hc = heavy_coords(res)
            if hc.shape[0] == 0:
                continue
            label, sortkey = res_label(chain.id, res)
            ca = res["CA"].coord if "CA" in res else None
            out.append((label, sortkey, ca, hc))
    return out


def choose_ligand_copy(model, resname):
    resname = resname.strip().upper()
    prot = protein_residues(model)
    prot_atoms = np.vstack([hc for (_, _, _, hc) in prot]) if prot else np.empty((0, 3))
    copies = []
    for chain in model:
        for res in chain:
            if res.resname.strip().upper() != resname:
                continue
            hc = heavy_coords(res)
            if hc.shape[0] == 0:
                continue
            if prot_atoms.shape[0]:
                d = np.linalg.norm(hc[:, None, :] - prot_atoms[None, :, :], axis=2)
                n = int((d.min(axis=1) <= CUTOFF).sum())
            else:
                n = 0
            copies.append((n, hc))
    if not copies:
        return None, 0, 0
    copies.sort(key=lambda t: t[0], reverse=True)
    return copies[0][1], len(copies), copies[0][0]


def pocket_by_contact(model, lig):
    hits = []
    for (label, sortkey, ca, hc) in protein_residues(model):
        if np.linalg.norm(hc[:, None, :] - lig[None, :, :], axis=2).min() <= CUTOFF:
            hits.append((label, sortkey))
    hits.sort(key=lambda t: t[1])
    return [lab for (lab, _) in hits]


def main():
    out = ["\t".join(["target_id", "self_uniprot", "holo_pdb",
                      "holo_residues", "apo_pdb", "apo_residues"])]
    n_ok = 0
    print(f"Building redundant holo motifs (cutoff {CUTOFF} A)\n")
    for target, uni, lig, src in ENTRIES:
        hp = os.path.join(HOLO_DIR, f"{target}.pdb")
        if not os.path.exists(hp):
            print(f"  [SKIP] {target}: holo PDB not found ({hp})"); continue
        try:
            m = load_first_model(hp)
            coords, ncopy, ncon = choose_ligand_copy(m, lig)
            if coords is None:
                print(f"  [SKIP] {target}: ligand {lig} not found"); continue
            if ncopy > 1:
                print(f"  [note] {target}: {ncopy} copies of {lig}; kept one with {ncon} contacts")
            res = pocket_by_contact(m, coords)
            if not res:
                print(f"  [SKIP] {target}: empty pocket"); continue
            out.append("\t".join([target, uni, hp, ",".join(res), "", ""]))
            n_ok += 1
            flag = "" if src == "drug-matched" else f"  [{src}]"
            print(f"  [OK]  {target} ({uni}, {lig}): {len(res)} holo res{flag}")
        except Exception as e:
            print(f"  [ERROR] {target}: {e}")
    with open(OUT_TSV, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"\nDone. {n_ok} entries -> {OUT_TSV}")


if __name__ == "__main__":
    main()
