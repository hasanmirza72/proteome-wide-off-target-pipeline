#!/bin/bash
# =============================================================================
# run_offtarget.sh
# -----------------------------------------------------------------------------
# The real off-target search (roadmap Part 3.2). Same manifest, same query form
# as the self-recovery gate, but tuned for REPORTING rather than locating:
#   - --top 100            keep only the strongest 100 hits per query
#   - quality filters      require a coherent, well-superposing motif
# The query protein's OWN AF model (self) is NOT removed here; that is done in
# offtarget_overlap_filter.py (Part 3.4), so the raw hits stay auditable.
#
# RUN THIS ONLY AFTER the self-recovery gate passes for most targets.
# =============================================================================

#SBATCH -A berzelius-2026-12        # <-- verify your project account
#SBATCH -n 32                    # <-- adjust to your cluster policy
#SBATCH -t 04:00:00
#SBATCH -J fd_offtarget
#SBATCH -o logs/offtarget_%j.out
#SBATCH -e logs/offtarget_%j.err

set -euo pipefail

# ---- Config -----------------------------------------------------------------
FOLDDISCO=/proj/berzelius-2021-29/users/x_miali/new_startegy/tools/folddisco/target/release/folddisco
INDEX=./folddisco_af_index
MANIFEST=queries_manifest.tsv
OUTDIR=hits_offtarget
THREADS=32

# Same extended-search tolerance you used for the gate (keep holo == apo).
DIST=0.5
ANGLE=5.0

TOP=100

# Quality filters for reported hits. Tune after inspecting the gate output.
#   --connected-node : matched residues must form a connected motif of >= N nodes
#   --rmsd           : drop matches whose superposition RMSD exceeds this (A)
#   --plddt          : ignore matches whose region is below this AF confidence
CONNECTED_NODE=3
MAX_RMSD=2.5
MIN_PLDDT=50

COLS=tid,idf,node_count,edge_count,min_rmsd,total_match_count,plddt,matching_residues,query_residues
# -----------------------------------------------------------------------------

mkdir -p "$OUTDIR" logs

run_query () {
    local pdb="$1" res="$2" out="$3" label="$4"
    if [[ -z "$pdb" || -z "$res" ]]; then
        echo "    [skip] $label (no structure/residues)"
        return
    fi
    if [[ ! -f "$pdb" ]]; then
        echo "    [WARN] $label: structure not found: $pdb" >&2
        return
    fi
    "$FOLDDISCO" query \
        -p "$pdb" \
        -q "$res" \
        -i "$INDEX" \
        -t "$THREADS" \
        -d "$DIST" -a "$ANGLE" \
        --per-structure \
        --top "$TOP" \
        --connected-node "$CONNECTED_NODE" \
        --rmsd "$MAX_RMSD" \
        --plddt "$MIN_PLDDT" \
        --header \
        --format-output "$COLS" \
        -o "$out"
    echo "    [ok]   $label -> $out"
}

while IFS=$'\t' read -r target self_uniprot holo_pdb holo_res apo_pdb apo_res; do
    [[ -z "${target:-}" || "$target" == \#* ]] && continue
    echo "=== $target  (self = $self_uniprot) ==="
    run_query "$holo_pdb" "$holo_res" "$OUTDIR/${target}_holo.tsv" "holo"
    run_query "$apo_pdb"  "$apo_res"  "$OUTDIR/${target}_apo.tsv"  "apo"
done < <(tail -n +2 "$MANIFEST")

echo
echo "Off-target queries done. Next:"
echo "  1) python3 build_library_pockets.py       # once, builds library_pockets.tsv"
echo "  2) python3 offtarget_overlap_filter.py    # annotates + filters these hits"
