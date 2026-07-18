# Supplementary Data

Datasets underlying the MSc thesis *Proteome-Wide Structural Prediction of Drug
Off-Targets*. These files contain the full data behind the tables and figures in the
thesis, so the results can be inspected and verified independently.

Author: Mirza Muhammad Hasan Ali · Stockholm University / SciLifeLab

## Contents

The package is organised by pipeline stage. Each file is a plain tab-separated (`.tsv`)
or comma-separated (`.csv`) table that opens in Excel, R, Python, or any text editor.

### 1. Benchmark dataset

| File | Description | Thesis reference |
|------|-------------|------------------|
| `dataset_final_nonredundant.csv` | The 23 non-redundant benchmark targets (holo/apo pairs, drug, family) | Table 4 |
| `dataset_excluded.csv` | The pairs removed during curation, with reasons | Table 2 |
| `queries_manifest.tsv` | The ligand-contact-residue query motifs searched by FoldDisco | Methods §9.3 |
| `dataset_final.csv` | The fuller dataset (33 scorable pairs) before non-redundancy collapse | Chapter 3 (Table 2)  |
| `verified_alphafold_targets.csv` | AlphaFold-model pocket detection per target | Chapter 3 (suppl. to Tables 3–4) |
| `library_pockets.tsv` | The P2Rank pocket library (16,325 proteins) | Chapter 2 (§2.4, Table 1) |

### 2. Detection and self-recovery

| File | Description | Thesis reference |
|------|-------------|------------------|
| `p2rank_evaluation.csv` | Per-target P2Rank pocket detection (holo and apo) | Table 3, Table 4 |
| `self_recovery_summary.tsv` | FoldDisco self-recovery ranks (positive control) | Table 9, Table 10 |

### 3. Off-target predictions

| File | Description | Thesis reference |
|------|-------------|------------------|
| `folddisco_raw_hits.tsv` | All raw FoldDisco geometric matches (unfiltered) | Section 5.2 |
| `offtarget_hits_annotated.tsv` | Raw hits with gene/family annotation (1,107 hits) | Table 8 (funnel) |
| `offtarget_hits_filtered.tsv` | Pocket-supported hits (681) | Table 8 (funnel) |
| `offtarget_hits_high.tsv` | High-confidence hits (58 candidates) | Table 7, Table 8 |
| `offtarget_hits_bioannotated.tsv` | High-confidence hits with ChEMBL/SIDER/openFDA annotation | Table 13, Table 14 |

### 4. Conformational analysis

| File | Description | Thesis reference |
|------|-------------|------------------|
| `dual_query_divergence.tsv` | Holo/apo hit-list divergence (Jaccard) vs conformational axes | Table 5, Figure 3-4 |
| `redundant_robustness.tsv` | Robustness to which drug's pocket defines the query | Chapter 4.5 (Table 6) |

### 5. Biological validation

| File | Description | Thesis reference |
|------|-------------|------------------|
| `phase5_recall.tsv` | Documented off-target recall | Table 15 |
| `phase5_enrichment_null.tsv` | Enrichment against a pocket-matched null | Table 16 |
| `phase5_pathways.tsv` | GO/KEGG/Reactome functional enrichment | Table 17 |

### 6. Baseline comparison

| File | Description | Thesis reference |
|------|-------------|------------------|
| `blast_hits.tsv` | BLAST+ sequence-similarity baseline results | Tables 11-12 |
| `foldseek_hits.tsv` | Foldseek structure-similarity baseline results | Tables 11-12 |
| `self_recovery_by_method.tsv` | Self-recovery of all 23 targets per method (median self-rank, count at rank 1, top 5, not returned) | Table 9 |
| `self_recovery_hard_targets.tsv` | Self-rank of the seven detection-hard targets, per method | Table 10 |
| `baseline_comparison.tsv` | Documented off-target recall by method at rank cutoffs 10, 30, 100, with median rank of recovered targets | Table 11 | 
| `baseline_split.tsv`| That recall@100 split into three tiers by the off-target's structural relationship to the query pocket | Table 12 |

## File formats

- All tables are UTF-8 encoded, tab-separated (`.tsv`) or comma-separated (`.csv`).
- The first row of every file is a header naming the columns.
- Missing values are empty cells (not `NA` or `0`), except where a column explicitly
  encodes a status such as `no_measured_value`.

See `DATA_DICTIONARY.md` for the meaning of the columns in the main tables.

## How these were produced

Every file is the output of a script in the accompanying code repository, run on the
Berzelius compute cluster. The pipeline that generates them is documented in the
repository's `docs/PIPELINE.md`. Live web-service data (ChEMBL, UniProt, etc.) was
accessed in July 2026.
