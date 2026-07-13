#!/bin/bash
# =============================================================================
# run_self_recovery.sh
# -----------------------------------------------------------------------------
# FoldDisco self-recovery GATE (roadmap Part 3.1).
# For each target: query its holo pocket (and its apo pocket) against the
# whole-proteome AlphaFold index, and save the ranked per-structure hit list.
# The question this answers: does a target's OWN AlphaFold model come back
# near the top of its own pocket query? Parse the outputs afterwards with
# summarize_self_recovery.py.
#
# NOTE: this only RUNS the queries. Rank/pass-fail is computed by the parser,
# so you can re-tune the pass criterion without re-querying.
# =============================================================================

# ---- SLURM (optional) -------------------------------------------------------
# These queries are fast (seconds each, ~46 total), so you can equally run this
# inside an interactive `salloc`/`interactive` session. If you do submit it,
# VERIFY the account and node request for your cluster (Berzelius shown).
#SBATCH -A berzelius-2026-12        # <-- verify your project account
#SBATCH -n 32                    # <-- Berzelius allocates CPUs per GPU fraction; adjust to your policy
#SBATCH -t 02:00:00
#SBATCH -J fd_selfrec
#SBATCH -o logs/selfrec_%j.out
#SBATCH -e logs/selfrec_%j.err

set -euo pipefail

# ---- Config: edit these -----------------------------------------------------
FOLDDISCO=/proj/berzelius-2021-29/users/x_miali/new_startegy/tools/folddisco/target/release/folddisco
INDEX=./folddisco_af_index          # your v6 index (no extension needed)
MANIFEST=queries_manifest.tsv
OUTDIR=hits
THREADS=32

# Extended-search tolerances. KEEP THESE IDENTICAL FOR HOLO AND APO so that
# conformation is the only variable in the dual-query comparison.
DIST=0.5      # -d  distance threshold (Angstrom)   [folddisco default 0.5]
ANGLE=5.0     # -a  angle threshold (degrees)        [folddisco default 5.0]

# For the GATE we keep a generous window so the self can be LOCATED even when it
# ranks poorly (e.g. collapsed 6YG2). For the later off-target search you'd use 100.
TOP=1000

# per-structure columns to emit (all valid per-structure keys)
COLS=tid,idf,node_count,edge_count,min_rmsd,total_match_count,matching_residues,query_residues
# -----------------------------------------------------------------------------

mkdir -p "$OUTDIR" logs

run_query () {
    # $1 pdb  $2 residues  $3 output file  $4 label
    local pdb="$1" res="$2" out="$3" label="$4"
    if [[ -z "$pdb" || -z "$res" ]]; then
        echo "    [skip] $label  (no structure/residues in manifest)"
        return
    fi
    if [[ ! -f "$pdb" ]]; then
        echo "    [WARN] $label: structure file not found: $pdb" >&2
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
        --header \
        --format-output "$COLS" \
        -o "$out"
    echo "    [ok]   $label -> $out"
}

# Read the manifest (skip the header row), one tab-delimited record at a time.
while IFS=$'\t' read -r target self_uniprot holo_pdb holo_res apo_pdb apo_res; do
    [[ -z "${target:-}" || "$target" == \#* ]] && continue
    echo "=== $target  (self = $self_uniprot) ==="
    run_query "$holo_pdb" "$holo_res" "$OUTDIR/${target}_holo.tsv" "holo"
    run_query "$apo_pdb"  "$apo_res"  "$OUTDIR/${target}_apo.tsv"  "apo"
done < <(tail -n +2 "$MANIFEST")

echo
echo "All queries finished. Now summarize:"
echo "    python3 summarize_self_recovery.py"
