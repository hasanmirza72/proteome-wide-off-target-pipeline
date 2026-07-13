#!/bin/bash
# =============================================================================
# run_redundant.sh  (roadmap Part 2.6 robustness sub-study)
# Runs the SAME off-target search on the extra redundant co-crystal entries
# (holo pocket per drug) against the EXISTING index -- no re-indexing.
# Same filters as run_offtarget.sh so the lists are directly comparable.
# =============================================================================
#SBATCH -A berzelius-2026-12
#SBATCH -n 32
#SBATCH -t 01:00:00
#SBATCH -J fd_redundant
#SBATCH -o logs/redundant_%j.out
#SBATCH -e logs/redundant_%j.err

set -euo pipefail

FOLDDISCO=/proj/berzelius-2021-29/users/x_miali/new_startegy/tools/folddisco/target/release/folddisco
INDEX=./folddisco_af_index
MANIFEST=redundant_manifest.tsv
OUTDIR=hits_redundant
THREADS=32
DIST=0.5; ANGLE=5.0; TOP=100
CONNECTED_NODE=3; MAX_RMSD=2.5; MIN_PLDDT=50
COLS=tid,idf,node_count,edge_count,min_rmsd,total_match_count,plddt,matching_residues,query_residues

mkdir -p "$OUTDIR" logs

while IFS=$'\t' read -r target self_uniprot holo_pdb holo_res apo_pdb apo_res; do
    [[ -z "${target:-}" || "$target" == \#* ]] && continue
    [[ -z "$holo_pdb" || -z "$holo_res" ]] && { echo "[skip] $target"; continue; }
    echo "=== $target (self=$self_uniprot) ==="
    "$FOLDDISCO" query -p "$holo_pdb" -q "$holo_res" -i "$INDEX" -t "$THREADS" \
        -d "$DIST" -a "$ANGLE" --per-structure --top "$TOP" \
        --connected-node "$CONNECTED_NODE" --rmsd "$MAX_RMSD" --plddt "$MIN_PLDDT" \
        --header --format-output "$COLS" -o "$OUTDIR/${target}_holo.tsv"
    echo "  -> $OUTDIR/${target}_holo.tsv"
done < <(tail -n +2 "$MANIFEST")

echo "Done. Next: python3 redundant_robustness.py"
