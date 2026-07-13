#!/usr/bin/env python3
"""
Evaluate P2Rank predictions clean (v4 - PyMOL Engine)
- Uses headless PyMOL to perfectly parse both .pdb and .cif files.
- Accurately calculates exact biological binding residues.
- Reports strict (4.0A) and loose (5.5A) detection rates.
"""

import csv
import glob
import math
import os
import sys

# Try to import PyMOL
try:
    import pymol
    from pymol import cmd
    # Start PyMOL in headless (no GUI) mode, quietly
    pymol.finish_launching(['pymol', '-cq'])
    _PYMOL_AVAILABLE = True
except ImportError:
    print("❌ ERROR: PyMOL python module not found.")
    print("Please ensure you are in a conda environment where 'pymol' is installed (e.g., conda install -c conda-forge pymol).")
    sys.exit(1)

# ---------------- paths ----------------
PAIRS_CSV = "apo_holo_pairs_clean.csv"
HOLO_PDB_DIR = "dataset_monomers_pdb"          
HOLO_PRED_DIR = "p2rank_out_holo"              
APO_PRED_DIR = "p2rank_out_apo"                

DCA_STRICT = 4.0
DCA_LOOSE = 5.5
TOP_N = 3
CONTACT_DIST = 4.0          
OUT_CSV = "p2rank_evaluation.csv"


# ================= ground truth: drug atoms =================
def read_ligand_atoms(holo_pdb_path, lig_resname):
    """Extracts True ligand coordinates using plain text to avoid PyMOL selection bugs with weird ligand names."""
    lig = (lig_resname or "").strip().upper()
    coords = []
    with open(holo_pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM")) and line[17:20].strip().upper() == lig:
                try:
                    coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    pass
    return coords


# ================= PyMOL Residue Engine =================
def get_true_binding_residues_pymol(struct_path, struct_name, ligand_atoms):
    """
    Uses PyMOL to load the structure, create pseudoatoms representing the drug, 
    and asks PyMOL to find all protein residues within CONTACT_DIST.
    """
    if not ligand_atoms or not os.path.exists(struct_path):
        return set()

    cmd.delete("all") # Clear PyMOL workspace
    
    # 1. Load the protein (handles PDB or CIF perfectly)
    cmd.load(struct_path, "protein")
    
    # 2. Build the "ghost drug" inside PyMOL using the coordinates
    for i, (x, y, z) in enumerate(ligand_atoms):
        cmd.pseudoatom("ghost_drug", pos=[x, y, z])
        
    # 3. Ask PyMOL to select all polymer residues near the ghost drug
    cmd.select("binding_site", f"polymer and (protein within {CONTACT_DIST} of ghost_drug)")
    
    # 4. Extract the Chain and Residue Number for the P2Rank comparison
    binding_residues = set()
    
    # Iterate over the selection to get chain and resi
    space = {'binding_residues': binding_residues}
    cmd.iterate("binding_site and name CA", "binding_residues.add(f'{chain}_{resi}')", space=space)
    
    return space['binding_residues']


# ================= P2Rank predictions =================
def read_pockets(pred_csv):
    pockets = []
    with open(pred_csv, encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.reader(f)
        header = [h.strip() for h in next(reader)]
        col = {n: i for i, n in enumerate(header)}
        if not all(k in col for k in ("rank", "center_x", "center_y", "center_z")):
            return []
        ridx = col.get("residue_ids")
        for row in reader:
            if not row or not row[col["rank"]].strip():
                continue
            try:
                pk = {
                    "rank": int(float(row[col["rank"]])),
                    "center": (float(row[col["center_x"]]), float(row[col["center_y"]]), float(row[col["center_z"]])),
                    "residues": (row[ridx].split() if ridx is not None and ridx < len(row) else []),
                }
            except (ValueError, IndexError):
                continue
            pockets.append(pk)
    return sorted(pockets, key=lambda p: p["rank"])

def find_pred_file(pred_dir, pdb_id):
    pid = pdb_id.lower()
    cands = [p for p in glob.glob(os.path.join(pred_dir, "*_predictions.csv")) if pid in os.path.basename(p).lower()]
    return cands[0] if cands else None


# ================= geometry / metrics =================
def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)

def centroid(atoms):
    n = len(atoms)
    return (sum(a[0] for a in atoms)/n, sum(a[1] for a in atoms)/n, sum(a[2] for a in atoms)/n)

def score_structure(pockets, ligand_atoms, true_res):
    r = {"hit_strict": 0, "hit_loose": 0, "best_dca": None, "best_dcc": None,
         "closest_rank": None, "res_recall": None, "res_jaccard": None, "n_true_res": len(true_res)}
    if not pockets or not ligand_atoms:
        return r
    cen = centroid(ligand_atoms)
    best_dca = best_dcc = math.inf
    closest = None
    for p in pockets:
        dca = min(dist(p["center"], a) for a in ligand_atoms)
        dcc = dist(p["center"], cen)
        best_dca = min(best_dca, dca)
        best_dcc = min(best_dcc, dcc)
        if closest is None or dca < closest["_dca"]:
            closest = dict(p, _dca=dca)
        if dca <= DCA_STRICT and r["hit_strict"] == 0:
            r["hit_strict"] = p["rank"]
        if dca <= DCA_LOOSE and r["hit_loose"] == 0:
            r["hit_loose"] = p["rank"]
    r["best_dca"] = round(best_dca, 2)
    r["best_dcc"] = round(best_dcc, 2)
    r["closest_rank"] = closest["rank"]
    if true_res and closest["residues"]:
        pred = set(closest["residues"])
        inter = len(true_res & pred)
        r["res_recall"] = round(inter / len(true_res), 3)
        union = len(true_res | pred)
        r["res_jaccard"] = round(inter / union, 3) if union else 0.0
    return r

def mcnemar_exact(b, c):
    n = b + c
    if n == 0: return 1.0
    k = min(b, c)
    tail = sum(math.comb(n, i) for i in range(k + 1)) * (0.5 ** n)
    return min(1.0, 2 * tail)


# ================= main =================
def main():
    if not os.path.exists(PAIRS_CSV):
        print(f"Error: '{PAIRS_CSV}' not found."); return

    with open(PAIRS_CSV, encoding="utf-8", errors="ignore", newline="") as f:
        pairs = list(csv.DictReader(f))

    rows = []
    for p in pairs:
        holo, apo, lig = p["holo_pdb"], p["apo_pdb"], p["holo_ligand"]
        apo_file = p.get("apo_file", "")

        holo_file = os.path.join(HOLO_PDB_DIR, f"{holo}.pdb")
        if not os.path.exists(holo_file):
            holo_file = os.path.join(HOLO_PDB_DIR, f"{holo.lower()}.pdb")
        
        ligand_atoms = read_ligand_atoms(holo_file, lig) if os.path.exists(holo_file) else []

        # USE PYMOL ENGINE FOR RESIDUES!
        holo_true = get_true_binding_residues_pymol(holo_file, "holo", ligand_atoms)
        
        # Determine the correct apo path (it might be in a different folder depending on how you segregated it)
        apo_full_path = apo_file if os.path.exists(apo_file) else os.path.join("dataset_apo_pdb", os.path.basename(apo_file))
        apo_true = get_true_binding_residues_pymol(apo_full_path, "apo", ligand_atoms)

        row = {"holo": holo, "apo": apo, "ligand": lig, "n_lig_atoms": len(ligand_atoms)}
        for tag, pred_dir, match_id, true_res in (
                ("holo", HOLO_PRED_DIR, holo, holo_true),
                ("apo", APO_PRED_DIR, apo, apo_true)):
            pf = find_pred_file(pred_dir, match_id)
            if not pf:
                for suffix in ("hit_strict", "hit_loose", "best_dca", "best_dcc",
                               "closest_rank", "res_recall", "res_jaccard", "n_true_res"):
                    row[f"{tag}_{suffix}"] = "NA"
                row[f"{tag}_hit_strict"] = "no_pred_file"
                continue
            m = score_structure(read_pockets(pf), ligand_atoms, true_res)
            for kk, vv in m.items():
                row[f"{tag}_{kk}"] = "NA" if vv is None else vv
        rows.append(row)

    if rows:
        with open(OUT_CSV, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader()
            w.writerows(rows)

    def success_flags(tag, field):
        out = {}
        for r in rows:
            v = r.get(f"{tag}_{field}")
            out[r["holo"]] = (1 if 1 <= v <= TOP_N else 0) if isinstance(v, int) else None
        return out

    print("=" * 74)
    print(f"P2RANK EVALUATION   (Top-{TOP_N};  strict DCA<={DCA_STRICT}A  loose DCA<={DCA_LOOSE}A)")
    print("=" * 74)

    for thr_field, label in (("hit_strict", f"strict (<= {DCA_STRICT} A)"),
                             ("hit_loose", f"loose  (<= {DCA_LOOSE} A)")):
        h = success_flags("holo", thr_field)
        a = success_flags("apo", thr_field)
        paired = [k for k in h if h[k] is not None and a[k] is not None]
        n = len(paired)
        hs = sum(h[k] for k in paired)
        as_ = sum(a[k] for k in paired)
        print(f"\n--- {label} --- (paired, scorable: {n})")
        if n:
            print(f"  holo Top-{TOP_N}: {hs}/{n} = {100*hs/n:.1f}%")
            print(f"  apo  Top-{TOP_N}: {as_}/{n} = {100*as_/n:.1f}%")
            b = sum(1 for k in paired if h[k] == 1 and a[k] == 0)
            c = sum(1 for k in paired if h[k] == 0 and a[k] == 1)
            print(f"  discordant: holo-only={b}, apo-only={c}  "
                  f"-> McNemar exact p = {mcnemar_exact(b, c):.3f}")

    def mean_recall(tag):
        vals = [r[f"{tag}_res_recall"] for r in rows
                if isinstance(r.get(f"{tag}_res_recall"), (int, float))
                and isinstance(r.get(f"{tag}_hit_strict"), int) and r[f"{tag}_hit_strict"] != 0]
        return (sum(vals)/len(vals), len(vals)) if vals else (None, 0)

    print("\n--- residue overlap of closest pocket (strict hits only) ---")
    for tag in ("holo", "apo"):
        mr, k = mean_recall(tag)
        print(f"  {tag}: no scorable residue data" if mr is None
              else f"  {tag}: mean binding-residue recall = {mr:.2f}  (over {k} hits)")

    cmd.quit() # Safely close the headless PyMOL instance

if __name__ == "__main__":
    main()
