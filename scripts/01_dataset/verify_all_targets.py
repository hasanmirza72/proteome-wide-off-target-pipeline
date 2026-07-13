#!/usr/bin/env python3
"""
Batch verification: for each holo/drug target, label which P2Rank pocket on the
protein's AlphaFold model is the true (drug-defined) binding site, via
sequence-guided superposition + DCA, and score detection the Phase-1 way.

Scope: this is the GROUND-TRUTH LABELLER. It tells you the true AF pocket id per
target; you then cross-reference that against FoldDisco's self-recovery output.
It does not run FoldDisco itself.
"""

import csv
import gzip
import os
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, Superimposer, is_aa
from Bio.Align import PairwiseAligner, substitution_matrices

# --- CONFIGURATION ---
DATASET_CSV = "dataset_final_nonredundant.csv"
HOLO_DIR = "dataset_monomers_pdb"
AF_STRUCT_DIR = "alphafold_structures"
AF_PRED_DIR = "alphafold_predictions"
OUT_CSV = "verified_alphafold_targets.csv"

TOP_N = 3          # Phase-1 Top-N ranking window
STRICT = 4.0       # strict DCA threshold (A)
LOOSE = 5.5        # relaxed DCA threshold (A)
MIN_CHAIN_CA = 10  # minimum Calpha count to treat a chain as a protein chain

_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def load_structure(path):
    opener = gzip.open if path.endswith(".gz") else open
    base = path[:-3] if path.endswith(".gz") else path
    parser = MMCIFParser(QUIET=True) if base.endswith((".cif", ".mmcif")) \
        else PDBParser(QUIET=True)
    with opener(path, "rt") as fh:
        return parser.get_structure("s", fh)


def resolve_holo_path(holo_id):
    """Try common extensions / casings for the holo file."""
    for stem in (holo_id, holo_id.lower(), holo_id.upper()):
        for ext in (".pdb", ".pdb.gz", ".cif", ".cif.gz", ".ent"):
            p = os.path.join(HOLO_DIR, stem + ext)
            if os.path.exists(p):
                return p
    return None


def chain_ca_seq(chain):
    cas, seq = [], []
    for res in chain:
        if not is_aa(res, standard=False) or "CA" not in res:
            continue
        cas.append(res["CA"])
        seq.append(_3TO1.get(res.resname.strip().upper(), "X"))
    return cas, "".join(seq)


def ligand_heavy_atoms(model, resname):
    """All heavy atoms of every copy of `resname` in the first model."""
    coords = []
    rn = resname.strip().upper()
    for chain in model:
        for res in chain:
            if res.resname.strip().upper() == rn:
                for atom in res:
                    el = (atom.element or "").strip()
                    if el and el != "H":
                        coords.append(atom.coord)
    return np.array(coords)


def choose_binding_chain(model, lig_coords):
    """Protein chain (>= MIN_CHAIN_CA Calpha) whose backbone is nearest the ligand."""
    best_id, best_d = None, np.inf
    for chain in model:
        ca_xyz = np.array([r["CA"].coord for r in chain
                           if is_aa(r, standard=False) and "CA" in r])
        if len(ca_xyz) < MIN_CHAIN_CA:
            continue
        d = np.min(np.linalg.norm(ca_xyz[:, None, :] - lig_coords[None, :, :], axis=2))
        if d < best_d:
            best_d, best_id = d, chain.id
    return best_id


def get_af_chain(model):
    """AlphaFold monomers are chain A; fall back to the first protein chain."""
    if "A" in model:
        return model["A"]
    for chain in model:
        if any(is_aa(r, standard=False) and "CA" in r for r in chain):
            return chain
    return None


def pair_by_alignment(seqH, casH, seqM, casM):
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -11
    aligner.extend_gap_score = -1
    aln = aligner.align(seqH, seqM)[0]
    fixed, moving = [], []
    for (h0, h1), (m0, m1) in zip(*aln.aligned):
        for k in range(h1 - h0):
            fixed.append(casH[h0 + k])
            moving.append(casM[m0 + k])
    return fixed, moving


def parse_p2rank(csv_path):
    pockets = []
    with open(csv_path) as fh:
        header = [h.strip() for h in fh.readline().split(",")]
        idx = {h: i for i, h in enumerate(header)}
        for line in fh:
            parts = line.split(",")
            if len(parts) <= idx["center_z"]:
                continue
            pockets.append({
                "rank": int(parts[idx["rank"]].strip()),
                "center": np.array([
                    float(parts[idx["center_x"]]),
                    float(parts[idx["center_y"]]),
                    float(parts[idx["center_z"]]),
                ]),
            })
    return pockets


def process_target(row):
    holo_id = row["holo"].strip()
    ligand = row["ligand"].strip()
    uniprot = row.get("uniprot", "").strip()
    if not uniprot:
        return {"status": "Error: No UniProt ID in CSV"}

    holo_path = resolve_holo_path(holo_id)
    if holo_path is None:
        return {"status": "Error: Holo file not found"}

    af_base = f"AF-{uniprot}-F1-model_v6.pdb.gz"
    af_path = os.path.join(AF_STRUCT_DIR, af_base)
    p2rank_csv = os.path.join(AF_PRED_DIR, f"{af_base}_predictions.csv")
    if not os.path.exists(af_path) or not os.path.exists(p2rank_csv):
        return {"status": "Error: AlphaFold files missing"}

    try:
        holo_model = next(iter(load_structure(holo_path)))
        af_model = next(iter(load_structure(af_path)))

        # Ligand first, independently of chain layout.
        lig_coords = ligand_heavy_atoms(holo_model, ligand)
        if lig_coords.size == 0:
            return {"status": f"Error: Ligand {ligand} not found / no heavy atoms"}

        # Protein chain = the one actually binding the drug.
        holo_chain_id = choose_binding_chain(holo_model, lig_coords)
        if holo_chain_id is None:
            return {"status": "Error: No protein chain found in holo"}
        af_chain = get_af_chain(af_model)
        if af_chain is None:
            return {"status": "Error: No protein chain in AF model"}

        casH, seqH = chain_ca_seq(holo_model[holo_chain_id])
        casM, seqM = chain_ca_seq(af_chain)

        fixed, moving = pair_by_alignment(seqH, casH, seqM, casM)
        if len(fixed) < MIN_CHAIN_CA:
            return {"status": "Error: Too few matched residues"}

        sup = Superimposer()
        sup.set_atoms(fixed, moving)        # maps AF (moving) -> holo (fixed)
        rot, tran = sup.rotran

        pockets = parse_p2rank(p2rank_csv)
        if not pockets:
            return {"status": "No pockets predicted in AF model",
                    "af_rmsd": round(sup.rms, 3)}

        rank_dca = []
        for p in pockets:
            c = np.dot(p["center"], rot) + tran     # AF centre -> holo frame
            dca = float(np.min(np.linalg.norm(lig_coords - c, axis=1)))
            rank_dca.append((p["rank"], round(dca, 3)))

        # LABEL: geometrically closest pocket overall (any rank).
        true_rank, true_dca = min(rank_dca, key=lambda rd: rd[1])

        # DETECTION (Phase-1 Top-N): any pocket ranked <= TOP_N within threshold.
        top_dcas = [d for (r, d) in rank_dca if r <= TOP_N]
        det_strict = any(d <= STRICT for d in top_dcas)
        det_loose = any(d <= LOOSE for d in top_dcas)

        # Ranking failure: true site is close but ranked outside Top-N.
        ranking_failure = (true_dca <= STRICT) and (true_rank > TOP_N)

        return {
            "status": "Success",
            "holo_chain_used": holo_chain_id,
            "n_matched_residues": len(fixed),
            "af_rmsd": round(sup.rms, 3),
            "n_pockets": len(pockets),
            "true_af_pocket_rank": true_rank,
            "true_af_pocket_dca": true_dca,
            "detected_top3_strict": det_strict,
            "detected_top3_loose": det_loose,
            "ranking_failure_strict": ranking_failure,
        }

    except Exception as e:
        return {"status": f"Error: {type(e).__name__}: {e}"}


def main():
    print("Starting AlphaFold verification...")
    with open(DATASET_CSV) as f:
        rows = list(csv.DictReader(f))

    results, n_success, n_detect = [], 0, 0
    for i, row in enumerate(rows):
        print(f"[{i+1}/{len(rows)}] {row['holo']} (drug {row['ligand']})...")
        out = process_target(row)
        results.append({**row, **out})
        if out["status"] == "Success":
            n_success += 1
            n_detect += int(out["detected_top3_strict"])
            flag = "DETECTED" if out["detected_top3_strict"] else (
                "RANKING-FAIL" if out["ranking_failure_strict"] else "MISS")
            print(f"   OK  chain {out['holo_chain_used']} | true pocket "
                  f"rank {out['true_af_pocket_rank']} DCA {out['true_af_pocket_dca']}A | {flag}")
        else:
            print(f"   {out['status']}")

    # Union of keys so error rows and success rows can share one CSV.
    fieldnames = []
    for r in results:
        for k in r:
            if k not in fieldnames:
                fieldnames.append(k)
    with open(OUT_CSV, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        for r in results:
            w.writerow({k: r.get(k, "") for k in fieldnames})

    print("\n" + "=" * 55)
    print(f"Mapped to AF: {n_success}/{len(rows)} | "
          f"Detected@Top3 (strict): {n_detect}/{n_success}")
    print(f"Saved: {OUT_CSV}")
    print("=" * 55)


if __name__ == "__main__":
    main()
