# Figures

Scripts that generate the thesis figures. All values are computed from the pipeline's
result files (or hard-coded from the verified Phase 1-5 outputs); the scripts only draw.

- `make_workflow_figure.py` — the methods workflow overview (Figure 1 / fig0_workflow.png)
- `make_thesis_figures.py`  — the nine result figures (detection, self-recovery,
  dual-query divergence, funnel, ChEMBL binding, recall, enrichment, family verdict, etc.)

Each script writes PNGs into a `figures/` subdirectory. Run from a folder where that
subdirectory is writable:
```bash
python3 make_workflow_figure.py
python3 make_thesis_figures.py
```
