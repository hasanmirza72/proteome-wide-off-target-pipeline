# Scripts

Organised by pipeline stage. Run `organize_repo.sh` from your project root to populate
these folders from your working `scripts/` directory. See `../docs/PIPELINE.md` for the
full ordered walkthrough.

- `01_dataset/` — build and curate the drug-bound benchmark
- `02_library/` — AlphaFold library and P2Rank pockets
- `03_search/` — FoldDisco index, queries, self-recovery, detection
- `04_filter/` — pocket-support and confidence filtering
- `05_analysis/` — conformational / dual-query analysis
- `06_validation/` — ChEMBL, recall, enrichment, pathways
- `07_baselines/` — BLAST and Foldseek
- `figures/` — figure generation
- `utils/` — audit and verification helpers
- `archive/` — superseded earlier iterations (kept for provenance)
