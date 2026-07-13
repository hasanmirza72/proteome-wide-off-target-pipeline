# Files to add from your Berzelius machine

This package already contains 9 data files. To complete it, copy the remaining files
below from your project root into this directory. All are confirmed present on your
machine.

## Core files (needed — referenced by thesis tables/figures)

| File | Referenced by |
|------|---------------|
| `dataset_final_nonredundant.csv` | Table 1, Table 4 (the 23 targets) |
| `dataset_excluded.csv` | Table 2 (curation exclusions) |
| `p2rank_evaluation.csv` | Table 3, Table 4 (detection) |
| `self_recovery_summary.tsv` | Table 15, Table 16 (self-recovery) |
| `dual_query_divergence.tsv` | Table 5, Figure 4 (divergence) |
| `redundant_robustness.tsv` | Chapter 5 (robustness) |
| `phase5_enrichment_null.tsv` | Table 11 (enrichment) |
| `phase5_pathways.tsv` | Table 12 (pathways) |

## Recommended extras (strengthen the package)

| File | Why include |
|------|-------------|
| `dataset_final.csv` | the fuller dataset (33 scorable) before non-redundancy collapse |
| `verified_alphafold_targets.csv` | ground-truth pocket labels underpinning detection scoring |
| `library_pockets.tsv` | the P2Rank pocket library (16,325 proteins) — the filter's input |
| `phase5_summary.txt` | human-readable summary of all four validation probes |

## Optional (small reference lists)

| File | Content |
|------|---------|
| `apo_holo_pairs_nonredundant.csv` | the apo/holo pairing step |
| `zero_pocket_proteins.txt` | the 4,225 proteins with no predicted pocket |
| `multi_fragment_proteins.txt` | the giant multi-fragment proteins |

## One command to copy everything (run from your project root)

```bash
cp dataset_final_nonredundant.csv dataset_excluded.csv p2rank_evaluation.csv \
   self_recovery_summary.tsv dual_query_divergence.tsv redundant_robustness.tsv \
   phase5_enrichment_null.tsv phase5_pathways.tsv \
   dataset_final.csv verified_alphafold_targets.csv library_pockets.tsv \
   phase5_summary.txt apo_holo_pairs_nonredundant.csv \
   zero_pocket_proteins.txt multi_fragment_proteins.txt \
   supplementary_data/
```

## Do NOT include (too large / binary / not data)

These are excluded on purpose — they are gigabyte-scale structures, binary search
indexes, API caches, or raw per-query output, not datasets a reader needs:

`af_by_accession/`, `alphafold_structures*/`, `alphafold_predictions/`,
`alphafold_pockets_clean/`, `dataset_*_pdb/`, `p2rank_out_*/`,
`library_blastdb.*`, `library_foldseekdb*`, `folddisco_af_index*`,
`chembl_cache/`, `enrichr_cache/`, `hits/`, `hits_offtarget/`, `hits_redundant/`,
`*.tar.gz`, `logs/`, `tmp/`.

Once the files are copied, delete this file — the package is then complete.
