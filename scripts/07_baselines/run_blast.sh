#!/bin/bash
# BLAST sequence baseline for the off-target comparison (Section 5.5 / 9.4).
# Library and queries already audited: 20,550 single-chain models, 23 benchmark targets.
#SBATCH -A berzelius-2026-12
#SBATCH -n 32
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --job-name=blast_baseline

set -euo pipefail

# 1. Build the indexed BLAST database (-parse_seqids for clean accession handling)
makeblastdb -in library_20550.fasta -dbtype prot -parse_seqids -out library_blastdb

# 2. Proteome-wide protein BLAST.
#    - outfmt 6 columns: query, subject, %identity, query-coverage, e-value, bitscore
#    - max_target_seqs 500 (not 100): retrieve generously, then YOU keep the top 100
#      by e-value when scoring, so the well-known max_target_seqs ordering caveat
#      cannot silently drop a true hit.
#    - num_threads: use available cores
blastp -query benchmark_targets.fasta \
       -db library_blastdb \
       -outfmt "6 qseqid sseqid pident qcovs evalue bitscore" \
       -max_target_seqs 500 \
       -num_threads 8 \
       -out blast_hits.tsv

echo "Done. Wrote blast_hits.tsv"
echo "Columns: qseqid  sseqid  pident  qcovs  evalue  bitscore"
echo "Next: sort by evalue within each qseqid and keep top 100 when scoring."
