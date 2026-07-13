#!/usr/bin/env python3
import csv, glob, math, os, sys

try:
    import pymol
    from pymol import cmd
    pymol.finish_launching(['pymol', '-cq'])
except ImportError:
    print("❌ ERROR: PyMOL python module not found."); sys.exit(1)

PAIRS_CSV = "apo_holo_pairs_nonredundant.csv" # CHANGE THIS TO NONREDUNDANT FOR YOUR SECOND RUN
HOLO_PDB_DIR = "dataset_monomers_pdb"          
HOLO_PRED_DIR = "p2rank_out_holo"              
APO_PRED_DIR = "p2rank_out_apo"                
DCA_STRICT, DCA_LOOSE, TOP_N, CONTACT_DIST = 4.0, 5.5, 3, 4.0          
OUT_CSV = "p2rank_evaluation.csv"

def read_ligand_atoms(holo_pdb_path, lig_resname):
    lig = (lig_resname or "").strip().upper()
    coords = []
    with open(holo_pdb_path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM")) and line[17:20].strip().upper() == lig:
                try: coords.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError: pass
    return coords

def get_true_binding_residues_pymol(struct_path, struct_name, ligand_atoms):
    if not ligand_atoms or not os.path.exists(struct_path): return set()
    cmd.delete("all")
    cmd.load(struct_path, "protein")
    for i, (x, y, z) in enumerate(ligand_atoms): cmd.pseudoatom("ghost_drug", pos=[x, y, z])
    # THE CRITICAL FIX: Added 'byres'
    cmd.select("binding_site", f"byres (polymer within {CONTACT_DIST} of ghost_drug)")
    space = {'binding_residues': set()}
    cmd.iterate("binding_site and name CA", "binding_residues.add(f'{chain}_{resi}')", space=space)
    return space['binding_residues']

def read_pockets(pred_csv):
    pockets = []
    with open(pred_csv, encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.reader(f)
        header = [h.strip() for h in next(reader)]
        col = {n: i for i, n in enumerate(header)}
        if not all(k in col for k in ("rank", "center_x", "center_y", "center_z")): return []
        ridx = col.get("residue_ids")
        for row in reader:
            if not row or not row[col["rank"]].strip(): continue
            try:
                pockets.append({
                    "rank": int(float(row[col["rank"]])),
                    "center": (float(row[col["center_x"]]), float(row[col["center_y"]]), float(row[col["center_z"]])),
                    "residues": (row[ridx].split() if ridx is not None and ridx < len(row) else []),
                })
            except (ValueError, IndexError): continue
    return sorted(pockets, key=lambda p: p["rank"])

def find_pred_file(pred_dir, pdb_id):
    pid = pdb_id.lower()
    cands = [p for p in glob.glob(os.path.join(pred_dir, "*_predictions.csv")) if pid in os.path.basename(p).lower()]
    return cands[0] if cands else None

def dist(a, b): return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)
def centroid(atoms): n = len(atoms); return (sum(a[0] for a in atoms)/n, sum(a[1] for a in atoms)/n, sum(a[2] for a in atoms)/n)

def score_structure(pockets, ligand_atoms, true_res):
    r = {"hit_strict": 0, "hit_loose": 0, "best_dca": None, "closest_rank": None, "res_recall": None, "n_true_res": len(true_res)}
    if not pockets or not ligand_atoms: return r
    best_dca = math.inf
    closest = None
    for p in pockets:
        dca = min(dist(p["center"], a) for a in ligand_atoms)
        best_dca = min(best_dca, dca)
        if closest is None or dca < closest["_dca"]: closest = dict(p, _dca=dca)
        if dca <= DCA_STRICT and r["hit_strict"] == 0: r["hit_strict"] = p["rank"]
        if dca <= DCA_LOOSE and r["hit_loose"] == 0: r["hit_loose"] = p["rank"]
    r["best_dca"] = round(best_dca, 2)
    r["closest_rank"] = closest["rank"]
    if true_res and closest["residues"]:
        inter = len(true_res & set(closest["residues"]))
        r["res_recall"] = round(inter / len(true_res), 3)
    return r

def main():
    if not os.path.exists(PAIRS_CSV): print(f"Error: '{PAIRS_CSV}' not found."); return
    with open(PAIRS_CSV, encoding="utf-8", errors="ignore", newline="") as f: pairs = list(csv.DictReader(f))
    rows = []
    for p in pairs:
        holo, apo, lig = p["holo_pdb"], p["apo_pdb"], p["holo_ligand"]
        apo_file = p.get("apo_file", "")
        holo_file = os.path.join(HOLO_PDB_DIR, f"{holo}.pdb")
        if not os.path.exists(holo_file): holo_file = os.path.join(HOLO_PDB_DIR, f"{holo.lower()}.pdb")
        ligand_atoms = read_ligand_atoms(holo_file, lig) if os.path.exists(holo_file) else []
        holo_true = get_true_binding_residues_pymol(holo_file, "holo", ligand_atoms)
        apo_full = apo_file if os.path.exists(apo_file) else os.path.join("dataset_apo_pdb", os.path.basename(apo_file))
        apo_true = get_true_binding_residues_pymol(apo_full, "apo", ligand_atoms)
        row = {"holo": holo, "apo": apo, "ligand": lig, "n_lig_atoms": len(ligand_atoms)}
        for tag, pred_dir, match_id, true_res in (("holo", HOLO_PRED_DIR, holo, holo_true), ("apo", APO_PRED_DIR, apo, apo_true)):
            pf = find_pred_file(pred_dir, match_id)
            if not pf:
                for suffix in ("hit_strict", "hit_loose", "best_dca", "closest_rank", "res_recall", "n_true_res"): row[f"{tag}_{suffix}"] = "NA"
                continue
            m = score_structure(read_pockets(pf), ligand_atoms, true_res)
            for kk, vv in m.items(): row[f"{tag}_{kk}"] = "NA" if vv is None else vv
        rows.append(row)

    print("=" * 74)
    print(f"P2RANK EVALUATION (v5 - PyMOL byres fix applied)")
    print("=" * 74)
    
    n_total = len(rows)
    for tag in ["holo", "apo"]:
        strict_hits = sum(1 for r in rows if isinstance(r.get(f"{tag}_hit_strict"), int) and 1 <= r[f"{tag}_hit_strict"] <= TOP_N)
        found_but_misranked = sum(1 for r in rows if isinstance(r.get(f"{tag}_hit_strict"), int) and r[f"{tag}_hit_strict"] > TOP_N and r.get(f"{tag}_best_dca", 99) <= DCA_STRICT)
        print(f"\n--- {tag.upper()} ---")
        print(f"  Top-{TOP_N} Success: {strict_hits}/{n_total} = {100*strict_hits/n_total:.1f}%")
        print(f"  Detected but Mis-Ranked (>Top 3): {found_but_misranked}")
        
        recalls = [r[f"{tag}_res_recall"] for r in rows if isinstance(r.get(f"{tag}_res_recall"), (int, float)) and isinstance(r.get(f"{tag}_hit_strict"), int) and 1 <= r[f"{tag}_hit_strict"] <= TOP_N]
        if recalls: print(f"  Mean Binding-Residue Recall: {sum(recalls)/len(recalls):.2f} (n={len(recalls)})")

    cmd.quit()

if __name__ == "__main__": main()
