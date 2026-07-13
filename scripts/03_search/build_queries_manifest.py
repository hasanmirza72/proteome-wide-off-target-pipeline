#!/usr/bin/env python3
"""
build_queries_manifest.py
-------------------------
Builds queries_manifest.tsv for the FoldDisco self-recovery / off-target search.

Output columns (tab-separated), consumed by the query drivers:
    target_id  self_uniprot  holo_pdb  holo_residues  apo_pdb  apo_residues

Residue strings are comma-joined FoldDisco labels: chain + auth_seq_id (+ icode),
e.g. "A7,A15,A32". Holo residues use holo numbering; apo residues use apo
numbering (each list must match its own -p file when queried).

Apo residue-set modes (roadmap Part 2.5):
  IDENTITY_FIXED_APO = False (default)  per-state contact set (matches Part 2.2)
  IDENTITY_FIXED_APO = True             same physical residues as holo, located in
                                        apo by nearest-Ca; isolates conformation

Handles gzip-compressed inputs (.pdb.gz / .cif.gz).
The apo file MUST be the AHoJ-ALIGNED apo (in the holo frame); the script warns
if the drug sits nowhere near the apo protein, and reports a *spatial* overlap
(fraction of holo pocket Ca that have an apo Ca within SPATIAL_TOL) so you can
tell a genuinely collapsed pocket from a mere residue-numbering difference.
"""

import csv
import gzip
import os
import sys
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, is_aa

# --- Configuration -----------------------------------------------------------
DATA_CSV  = "dataset_final_nonredundant.csv"
HOLO_DIR  = "dataset_monomers_pdb"
OUT_TSV   = "queries_manifest.tsv"

CUTOFF    = 4.0     # A; pocket = protein residues within CUTOFF of the drug
IDENTITY_FIXED_APO = False
MAP_TOL   = 3.0     # A; max Ca-Ca distance to accept a holo->apo residue match
SPATIAL_TOL = 2.0   # A; Ca-Ca distance that counts a holo pocket residue as
                    # "still present" in the apo (used for the overlap report)
MIN_LIG_CONTACTS = 1
UNALIGNED_WARN = 8.0  # A; nearest apo-atom-to-drug distance above which we warn
# -----------------------------------------------------------------------------


def _open_text(path):
    try:
        # First, try to read it as a gzip file, even if it doesn't end in .gz
        with gzip.open(path, "rt") as f:
            f.read(1) # Try reading 1 byte to see if it crashes
        return gzip.open(path, "rt")
    except OSError:
        # If it crashes, it's not a gzip file. Open it normally.
        return open(path, "rt")

def load_first_model(path):
    low = path.lower()
    stem = low[:-3] if low.endswith(".gz") else low
    parser = MMCIFParser(QUIET=True) if stem.endswith((".cif", ".mmcif")) \
        else PDBParser(QUIET=True)
    with _open_text(path) as f:
        structure = parser.get_structure("s", f)
    return next(iter(structure))          # first model (NMR-safe)


def heavy_coords(residue):
    cs = [a.coord for a in residue if (a.element or "").strip().upper() != "H"]
    return np.array(cs, dtype=float) if cs else np.empty((0, 3))


def res_label(chain_id, res):
    _, seq, icode = res.id
    icode = (icode or "").strip()
    return f"{chain_id}{seq}{icode}", (chain_id, int(seq), icode)


def protein_residues(model):
    """(label, sortkey, ca_coord|None, heavy_coords) for amino acids."""
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


def apo_ca_table(model):
    """(labels, Nx3 Ca array) for apo amino acids that have a Ca."""
    labels, cas = [], []
    for (lab, sk, ca, hc) in protein_residues(model):
        if ca is not None:
            labels.append(lab)
            cas.append(ca)
    cas = np.array(cas, dtype=float) if cas else np.empty((0, 3))
    return labels, cas


def choose_ligand_copy(model, resname):
    """Heavy-atom coords of the best-bound copy of `resname`."""
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
                n_contact = int((d.min(axis=1) <= CUTOFF).sum())
            else:
                n_contact = 0
            copies.append((n_contact, hc))
    if not copies:
        return None, 0, 0
    copies.sort(key=lambda t: t[0], reverse=True)
    return copies[0][1], len(copies), copies[0][0]


def pocket_by_contact(model, lig_coords):
    """(label, sortkey, ca) for residues within CUTOFF of lig_coords."""
    hits = []
    for (label, sortkey, ca, hc) in protein_residues(model):
        d = np.linalg.norm(hc[:, None, :] - lig_coords[None, :, :], axis=2).min()
        if d <= CUTOFF:
            hits.append((label, sortkey, ca))
    hits.sort(key=lambda t: t[1])
    return hits


def spatial_overlap(holo_pocket, apo_cas, tol):
    """Fraction of holo pocket residues (with Ca) whose nearest apo Ca <= tol.
    Meaningful regardless of numbering: a collapsed/unaligned apo -> low value."""
    n_total = n_kept = 0
    for (label, sortkey, ca) in holo_pocket:
        if ca is None or apo_cas.shape[0] == 0:
            continue
        n_total += 1
        if np.linalg.norm(apo_cas - ca[None, :], axis=1).min() <= tol:
            n_kept += 1
    return (n_kept / n_total) if n_total else 0.0


def map_holo_to_apo(holo_pocket, apo_labels, apo_cas):
    """Same physical residues as holo, located in apo by nearest-Ca (<= MAP_TOL)."""
    mapped, unmapped, used = [], 0, set()
    for (label, sortkey, ca) in holo_pocket:
        if ca is None or apo_cas.shape[0] == 0:
            unmapped += 1
            continue
        d = np.linalg.norm(apo_cas - ca[None, :], axis=1)
        j = int(d.argmin())
        if d[j] <= MAP_TOL and apo_labels[j] not in used:
            used.add(apo_labels[j])
            mapped.append(apo_labels[j])
        else:
            unmapped += 1
    return mapped, unmapped


def main():
    if not os.path.exists(DATA_CSV):
        sys.exit(f"Input CSV not found: {DATA_CSV}")
    with open(DATA_CSV) as f:
        rows = list(csv.DictReader(f))

    out = ["\t".join(["target_id", "self_uniprot", "holo_pdb",
                      "holo_residues", "apo_pdb", "apo_residues"])]
    n_ok = 0
    print(f"Mode: {'IDENTITY-FIXED apo' if IDENTITY_FIXED_APO else 'per-state contact apo'}"
          f" | cutoff {CUTOFF} A\n")

    for row in rows:
        target   = row["holo"].strip()
        uniprot  = row.get("uniprot", "").strip()
        ligand   = row["ligand"].strip()
        apo_file = row["apo_file"].strip()

        holo_path = os.path.join(HOLO_DIR, f"{target}.pdb")
        if not os.path.exists(holo_path):
            for alt in (f"{target.lower()}.pdb", f"{target}.pdb.gz",
                        f"{target.lower()}.pdb.gz", f"{target}.cif", f"{target}.cif.gz"):
                cand = os.path.join(HOLO_DIR, alt)
                if os.path.exists(cand):
                    holo_path = cand
                    break
        if not os.path.exists(holo_path) or not os.path.exists(apo_file):
            print(f"  [SKIP] {target}: missing holo/apo file")
            continue

        try:
            holo_model = load_first_model(holo_path)
            apo_model  = load_first_model(apo_file)

            lig_coords, n_copies, n_contacts = choose_ligand_copy(holo_model, ligand)
            if lig_coords is None:
                print(f"  [SKIP] {target}: ligand {ligand} not found in holo")
                continue
            if n_copies > 1:
                print(f"  [note] {target}: {n_copies} copies of {ligand}; "
                      f"kept the one with {n_contacts} protein contacts")
            if n_contacts < MIN_LIG_CONTACTS:
                print(f"  [WARN] {target}: chosen ligand copy has {n_contacts} contacts")

            holo_pocket = pocket_by_contact(holo_model, lig_coords)
            holo_labels = [lab for (lab, _, _) in holo_pocket]
            if not holo_labels:
                print(f"  [SKIP] {target}: empty holo pocket (frame/ligand issue?)")
                continue

            apo_labels_all, apo_cas_all = apo_ca_table(apo_model)

            # apo-alignment sanity: nearest apo atom to the ghost drug
            apo_all = protein_residues(apo_model)
            apo_atoms = np.vstack([hc for (_, _, _, hc) in apo_all]) if apo_all \
                else np.empty((0, 3))
            min_apo_lig = (np.linalg.norm(apo_atoms[:, None, :] - lig_coords[None, :, :],
                                          axis=2).min() if apo_atoms.shape[0] else float("inf"))
            aligned = min_apo_lig <= UNALIGNED_WARN
            if not aligned:
                print(f"  [WARN] {target}: nearest apo atom is {min_apo_lig:.1f} A from the "
                      f"drug — apo is NOT aligned into the holo frame (fix apo_file path)")

            # spatial overlap: how much of the holo pocket is structurally present in apo
            overlap = spatial_overlap(holo_pocket, apo_cas_all, SPATIAL_TOL) if aligned else 0.0

            if IDENTITY_FIXED_APO:
                apo_labels, n_unmapped = map_holo_to_apo(holo_pocket, apo_labels_all, apo_cas_all)
                if n_unmapped:
                    print(f"  [note] {target}: {n_unmapped}/{len(holo_labels)} holo residues "
                          f"unmatched in apo within {MAP_TOL} A (moved/collapsed)")
            else:
                apo_labels = [lab for (lab, _, _) in pocket_by_contact(apo_model, lig_coords)]

            if aligned and overlap < 0.5:
                print(f"  [note] {target}: apo spatial overlap {overlap:.2f} — pocket looks "
                      f"collapsed/displaced in the apo model (a real conformational change)")

            holo_res = ",".join(holo_labels)
            apo_res  = ",".join(apo_labels)
            if not apo_res:
                print(f"  [note] {target}: apo residue list empty — apo query will be skipped")

            out.append("\t".join([target, uniprot, holo_path, holo_res, apo_file, apo_res]))
            n_ok += 1
            ov = f"{overlap:.2f}" if aligned else "n/a"
            print(f"  [OK]  {target}: holo {len(holo_labels)} res, apo {len(apo_labels)} res, "
                  f"apo spatial overlap {ov}")

        except Exception as e:
            print(f"  [ERROR] {target}: {e}")

    with open(OUT_TSV, "w") as f:
        f.write("\n".join(out) + "\n")
    print(f"\nDone. {n_ok} targets written to {OUT_TSV}")


if __name__ == "__main__":
    main()
