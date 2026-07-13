#!/bin/bash
# Foldseek structural baseline for the off-target comparison (Section 5.5 / 9.4).
# Library: 20,550 single-chain AlphaFold models (af_by_accession/)
# Queries: 23 benchmark holo monomer structures (dataset_monomers_pdb/)
#SBATCH -A berzelius-2026-12
#SBATCH -n 32
#SBATCH -t 01:00:00
#SBATCH --cpus-per-task=8
#SBATCH --nodes=1
#SBATCH --job-name=flodseek_baseline
set -euo pipefail

echo "Building the Foldseek database..."
foldseek createdb af_by_accession/ library_foldseekdb

echo "Running the Foldseek search..."
# Key fairness settings vs the naive defaults:
#   -e 10           : loosen the E-value cutoff (default 0.001 would silently drop
#                     weak-but-real hits). BLAST default is 10; match it, then trim
#                     by rank yourself when scoring, so both methods report equally wide.
#   --max-seqs 500  : retrieve generously; keep top 100 by e-value at scoring time.
#   -a / extra cols : report fraction-identical and coverage so short spurious
#                     structural matches can be filtered the same way as BLAST.
foldseek easy-search foldseek_queries/ library_foldseekdb \
         foldseek_hits.tsv tmp \
         --format-output "query,target,evalue,bits,alntmscore,fident,qcov,tcov" \
         -e 10 \
         --max-seqs 500 \
         --threads 8

echo "Done. Wrote foldseek_hits.tsv"
echo "Columns: query target evalue bits alntmscore fident qcov tcov"
echo
echo "SANITY CHECK — confirm target IDs are clean accessions:"
echo "  cut -f2 foldseek_hits.tsv | head"
echo "  (should look like A0A024R1R8, not A0A024R1R8_A or A0A024R1R8.pdb)"
