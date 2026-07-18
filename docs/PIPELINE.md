# Pipeline: full ordered walkthrough

Run the stages in order. Each command reads the output of earlier stages. External tools
(FoldDisco, P2Rank, BLAST+, Foldseek) must be installed and on `PATH`.

## Stage 1 — Dataset construction
```
python scripts/01_dataset/sider_pdb_miner.py      # mine drug-bound human structures (RCSB + SIDER/openFDA/PubChem)
python scripts/01_dataset/ahoj_apo_miner.py       # retrieve matching apo structures (AHoJ)
python scripts/01_dataset/clean_pairs.py          # curate apo/holo pairs
python scripts/01_dataset/finalize_dataset.py     # final curation + exclusions
python scripts/01_dataset/verify_all_targets.py   # ground-truth pocket labeller (true AF pocket per target)
python scripts/01_dataset/check_quaternary_site.py # flag interface (inter-subunit) pockets
```
Produces: `dataset_final_nonredundant.csv` (23 non-redundant targets), `dataset_final.csv`.

## Stage 2 — Library and pockets
```
python scripts/02_library/strip_ligands.py           # make computationally-apo inputs for P2Rank
python scripts/02_library/library_20550_fasta.py     # sequences for the 20,550 models
python scripts/02_library/build_library_pockets.py    # P2Rank pocket table
python scripts/02_library/analyze_top3_pockets.py     # P2Rank detection evaluation
```
Produces: `library_20550.fasta`, `library_pockets.tsv`, `p2rank_evaluation.csv`.

## Stage 3 — FoldDisco search
```
python scripts/03_search/build_queries_manifest.py    # ligand-contact-residue queries
python scripts/03_search/build_redundant_manifest.py  # redundant-set manifest
bash   scripts/03_search/run_self_recovery.sh         # self-recovery (positive control)
bash   scripts/03_search/run_offtarget.sh             # off-target search
bash   scripts/03_search/run_redundant.sh             # redundant-set search
python scripts/03_search/summarize_self_recovery.py
python scripts/03_search/evaluate_predictions_nonredundant.py         # detection (N=23)
python scripts/03_search/evaluate_predictions_residue_call_nonred.py  # residue-level metric
```
Produces: `self_recovery_summary.tsv`, `hits/`, `hits_offtarget/`.

## Stage 4 — Off-target filtering
```
python scripts/04_filter/offtarget_overlap_filter.py  # pocket support + confidence tiers
```
Produces: `offtarget_hits_annotated.tsv` → `offtarget_hits_filtered.tsv` → `offtarget_hits_high.tsv` (58).

## Stage 5 — Conformational analysis
```
python scripts/05_analysis/dual_query_divergence.py   # holo/apo divergence (Spearman, permutation)
python scripts/05_analysis/redundant_robustness.py    # redundant-set robustness
```

## Stage 6 — Biological validation (Phase 5)
```
python scripts/06_validation/phase5_annotate.py         # ChEMBL/SIDER/openFDA annotation
python scripts/06_validation/phase5_recall.py           # documented off-target recall
python scripts/06_validation/phase5_enrichment_null.py  # enrichment vs matched null
python scripts/06_validation/phase5_pathways.py         # GO/KEGG/Reactome via Enrichr
python scripts/06_validation/phase5_summary.py          # consolidated summary
```
Produces: `offtarget_hits_bioannotated.tsv`, `phase5_recall.tsv`, `phase5_enrichment_null.tsv`,
`phase5_pathways.tsv`, `phase5_summary.txt`.

## Stage 7 — Baselines
```
python scripts/07_baselines/blast_foldseek_rebuild_queries_canonical.py
bash   scripts/07_baselines/run_blast.sh
bash   scripts/07_baselines/run_foldseek.sh
python scripts/07_baselines/score_baselines.py
```
Produces: `blast_hits.tsv`, `foldseek_hits.tsv`.

## Stage 8 — Figures
```
python scripts/figures/make_thesis_figures.py
python scripts/figures/make_workflow_figure.py
```

## Utilities (`scripts/utils/`)
Audit and verification helpers used during development: FASTA/fragment/target-count audits,
manifest checks, library verification, name lookups, dependency scan. Not required to run the
pipeline, but useful for checking data integrity.

## Archive (`scripts/archive/`)
Earlier iterations kept for provenance: `evaluate_predictions.py` and `_clean.py` variants were
run on pre-final datasets (raw and clean pairs); the thesis reports the non-redundant versions.
