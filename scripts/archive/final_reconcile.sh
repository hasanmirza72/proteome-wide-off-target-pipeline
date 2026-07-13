#!/bin/bash
# DEFINITIVE check: does the baseline query set == the canonical 23 non-redundant targets?
# Run from your project folder.

echo "=== A. What accessions are in benchmark_targets.fasta (baseline queries)? ==="
grep '>' benchmark_targets.fasta | sed 's/>//' | sort > /tmp/baseline_q.txt
echo "count: $(wc -l < /tmp/baseline_q.txt)"

echo
echo "=== B. What are the canonical 23 non-redundant targets? ==="
# Try to find the uniprot column automatically; adjust if needed.
head -1 dataset_final_nonredundant.csv
echo "--- extracting accessions (UniProt pattern) ---"
tail -n +2 dataset_final_nonredundant.csv | grep -oE '[OPQ][0-9][A-Z0-9]{3}[0-9]' | sort -u > /tmp/canon_q.txt
echo "count: $(wc -l < /tmp/canon_q.txt)"

echo
echo "=== C. THE RECONCILIATION ==="
echo "--- In baseline but NOT in canonical 23 (should be empty; Q5SLE7 = contamination): ---"
comm -23 /tmp/baseline_q.txt /tmp/canon_q.txt

echo "--- In canonical 23 but MISSING from baseline (should be empty!): ---"
comm -13 /tmp/baseline_q.txt /tmp/canon_q.txt

echo
echo "=== D. Same check against the FoldDisco self-recovery target set ==="
grep -oE '[OPQ][0-9][A-Z0-9]{3}[0-9]' self_recovery_summary.tsv | sort -u > /tmp/folddisco_q.txt
echo "FoldDisco self-recovery targets: $(wc -l < /tmp/folddisco_q.txt)"
echo "--- In FoldDisco but missing from baseline (baseline should search these!): ---"
comm -13 /tmp/baseline_q.txt /tmp/folddisco_q.txt
echo "--- In baseline but not in FoldDisco (contamination): ---"
comm -23 /tmp/baseline_q.txt /tmp/folddisco_q.txt

echo
echo "=== VERDICT ==="
echo "If section C/D show Q5SLE7 as 'extra' AND some real accession as 'missing',"
echo "then benchmark_targets.fasta has a swap: rebuild it from the CANONICAL file,"
echo "re-run BLAST + Foldseek on the corrected 23, and THEN score."
echo "If both 'missing' lists are empty, the baseline set is correct and we score now."
