#!/bin/bash
# Quick spot-checks to confirm the true target count and Q5SLE7's status.
# Run each line from inside your project folder.
#SBATCH -A berzelius-2026-12
#SBATCH -n 32
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --job-name=protein_check

echo "==================================================================="
echo "1. How many non-redundant targets in the canonical dataset file?"
echo "==================================================================="
wc -l dataset_final_nonredundant.csv
echo "   (subtract 1 for header = target count)"
echo "   Is Q5SLE7 in it?"
grep -c "Q5SLE7" dataset_final_nonredundant.csv

echo
echo "==================================================================="
echo "2. How many targets did FoldDisco self-recovery actually record?"
echo "==================================================================="
wc -l self_recovery_summary.tsv
echo "   (subtract 1 for header)"
echo "   Is Q5SLE7 in it?"
grep -c "Q5SLE7" self_recovery_summary.tsv

echo
echo "==================================================================="
echo "3. Is Q5SLE7 in the searchable library FASTA at all?"
echo "==================================================================="
grep -c "Q5SLE7" library_20550.fasta
echo "   (0 = not in library = cannot self-recover in ANY method)"

echo
echo "==================================================================="
echo "4. Distinct targets in each FoldDisco output (the real N):"
echo "==================================================================="
echo -n "  self_recovery_summary.tsv : "
tail -n +2 self_recovery_summary.tsv | grep -oE '[OPQ][0-9][A-Z0-9]{3}[0-9]' | sort -u | wc -l
echo -n "  phase5_recall.tsv         : "
tail -n +2 phase5_recall.tsv | grep -oE '[OPQ][0-9][A-Z0-9]{3}[0-9]' | sort -u | wc -l
echo -n "  offtarget_hits_annotated  : "
tail -n +2 offtarget_hits_annotated.tsv | cut -f1 | sed 's/_.*//' | sort -u | wc -l

echo
echo "==================================================================="
echo "5. Where did Q5SLE7 get dropped? Trace it back through the pipeline:"
echo "==================================================================="
for f in dataset_final.csv dataset_final_nonredundant.csv dataset_excluded.csv \
         verified_alphafold_targets.csv queries_manifest.tsv \
         apo_holo_pairs_nonredundant.csv library_20550.fasta \
         self_recovery_summary.tsv; do
  if [ -f "$f" ]; then
    c=$(grep -c "Q5SLE7" "$f" 2>/dev/null || echo 0)
    echo "  $f : $c"
  fi
done
echo "  ^ The last file where count>0 is where Q5SLE7 survived;"
echo "    the first file where it's 0 is where it was dropped."
