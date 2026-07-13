#!/usr/bin/env python3
"""
Clean the apo-holo pair list before final evaluation.

Applies three filters that a thesis examiner would expect, and writes both a
cleaned set and a non-redundant (one-per-UniProt) subset so you can report the pair.

Filters
  1. Ligand-size: drop pairs whose "drug" is really an ion/fragment (< MIN_LIG_ATOMS
     heavy atoms, e.g. OH, CO3, GLY). Their DCA ground truth is meaningless.
  2. Structural quality: drop pairs where the apo does not really contain the pocket
     (tm_score < TM_MIN or mapped_binding_residues_percent < MAPPED_MIN).
  3. Redundancy: for the non-redundant subset, keep ONE representative holo per
     UniProt (best resolution, then highest mapped%, then highest tm_score).

n_lig_atoms is taken from p2rank_evaluation.csv if present (identical to what the
evaluator scored), else counted directly from the holo PDB, else left unknown.
"""

import csv
import os

PAIRS_CSV = "apo_holo_pairs.csv"
EVAL_CSV = "p2rank_evaluation.csv"      # optional, for n_lig_atoms
HOLO_PDB_DIR = "dataset_monomers_pdb"   # fallback source for ligand atom count

OUT_CLEAN = "apo_holo_pairs_clean.csv"
OUT_NR = "apo_holo_pairs_nonredundant.csv"

MIN_LIG_ATOMS = 6
TM_MIN = 0.50
MAPPED_MIN = 80.0
BLACKLIST_LIG = {"OH", "CO3", "GLY", "ACT", "EDO", "GOL", "SO4", "PO4",
                 "NO3", "FLC", "DMS", "IPA", "MPD", "CL", "NA", "ZN", "MG", "CA", "K"}


def fnum(v, default=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def count_ligand_atoms_from_pdb(holo_pdb, lig):
    for cand in (f"{holo_pdb}.pdb", f"{holo_pdb.lower()}.pdb"):
        path = os.path.join(HOLO_PDB_DIR, cand)
        if os.path.exists(path):
            lig = (lig or "").strip().upper()
            n = 0
            with open(path, encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if line.startswith(("HETATM", "ATOM")) and line[17:20].strip().upper() == lig:
                        n += 1
            return n
    return None


def load_lig_atom_counts():
    counts = {}
    if os.path.exists(EVAL_CSV):
        with open(EVAL_CSV, encoding="utf-8", errors="ignore", newline="") as f:
            for r in csv.DictReader(f):
                v = fnum(r.get("n_lig_atoms"))
                if v is not None:
                    counts[r.get("holo", "").upper()] = int(v)
    return counts


def main():
    if not os.path.exists(PAIRS_CSV):
        print(f"Error: '{PAIRS_CSV}' not found."); return

    with open(PAIRS_CSV, encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f)
        fieldnames = reader.fieldnames
        rows = list(reader)

    eval_counts = load_lig_atom_counts()

    kept, removed = [], []
    for r in rows:
        holo = r["holo_pdb"].strip()
        lig = r["holo_ligand"].strip()
        tm = fnum(r.get("tm_score"))
        mapped = fnum(r.get("mapped_binding_residues_percent"))

        n_lig = eval_counts.get(holo.upper())
        if n_lig is None:
            n_lig = count_ligand_atoms_from_pdb(holo, lig)
        r["n_lig_atoms"] = n_lig if n_lig is not None else ""

        reasons = []
        if lig.upper() in BLACKLIST_LIG:
            reasons.append(f"non-drug ligand '{lig}'")
        if n_lig is not None and n_lig < MIN_LIG_ATOMS:
            reasons.append(f"ligand too small ({n_lig} atoms)")
        if tm is not None and tm < TM_MIN:
            reasons.append(f"low TM-score ({tm})")
        if mapped is not None and mapped < MAPPED_MIN:
            reasons.append(f"low mapped% ({mapped})")

        if reasons:
            removed.append((r, reasons))
        else:
            kept.append(r)

    # non-redundant: one representative per UniProt
    def better(a, b):
        ra, rb = fnum(a.get("resolution"), 1e9), fnum(b.get("resolution"), 1e9)
        if ra != rb:
            return a if ra < rb else b
        ma, mb = fnum(a.get("mapped_binding_residues_percent"), 0), fnum(b.get("mapped_binding_residues_percent"), 0)
        if ma != mb:
            return a if ma > mb else b
        ta, tb = fnum(a.get("tm_score"), 0), fnum(b.get("tm_score"), 0)
        return a if ta >= tb else b

    rep = {}
    for r in kept:
        key = (r.get("uniprot") or "").strip() or f"__{r['holo_pdb']}"
        rep[key] = r if key not in rep else better(rep[key], r)
    nonredundant = list(rep.values())

    out_fields = fieldnames + (["n_lig_atoms"] if "n_lig_atoms" not in fieldnames else [])
    for path, data in ((OUT_CLEAN, kept), (OUT_NR, nonredundant)):
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=out_fields)
            w.writeheader()
            w.writerows(data)

    # report
    print("=" * 70)
    print("   PAIR CLEANING REPORT")
    print("=" * 70)
    print(f"input pairs:        {len(rows)}")
    print(f"removed:            {len(removed)}")
    for r, reasons in removed:
        print(f"   - {r['holo_pdb']} ({r['holo_ligand']}) -> apo {r['apo_pdb']}: {'; '.join(reasons)}")
    print(f"kept (clean):       {len(kept)}   -> {OUT_CLEAN}")

    uni = {}
    for r in kept:
        uni.setdefault((r.get("uniprot") or "").strip(), []).append(r["holo_pdb"])
    dupes = {u: v for u, v in uni.items() if u and len(v) > 1}
    print(f"unique UniProts in clean set: {len([u for u in uni if u])}")
    if dupes:
        print("redundant proteins (collapsed in non-redundant subset):")
        for u, v in sorted(dupes.items(), key=lambda kv: -len(kv[1])):
            print(f"   - {u}: {len(v)}x  ({', '.join(v)})")
    print(f"non-redundant pairs: {len(nonredundant)}   -> {OUT_NR}")

    # role audit
    holo_ids = {r["holo_pdb"].upper() for r in kept}
    apo_ids = {r["apo_pdb"].upper() for r in kept}
    both = sorted(holo_ids & apo_ids)
    print("role audit:", f"[!] BOTH roles: {', '.join(both)}" if both
          else "clean (no PDB is both holo and apo)")

    print("\nNext: point the evaluator's PAIRS_CSV at each file and rerun:")
    print(f"   PAIRS_CSV = '{OUT_CLEAN}'   (main result)")
    print(f"   PAIRS_CSV = '{OUT_NR}'      (non-redundant sensitivity check)")


if __name__ == "__main__":
    main()
