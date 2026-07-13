# How to build and publish this repository

This scaffold contains the structure, documentation, and cleaned scripts for the
pipeline. Follow these steps on the machine where your project lives (the cluster or
wherever your `scripts/` folder and result files are).

## What's in this scaffold

- `README.md` — the repo's front page (what the project does, how to run it)
- `docs/PIPELINE.md` — the full step-by-step pipeline walkthrough
- `requirements.txt` — Python dependencies with exact versions
- `.gitignore` — tells Git to skip large data (structures, indexes, caches)
- `LICENSE` — MIT license
- `organize_repo.sh` — copies your 48 scripts into the clean stage folders
- `corrected_scripts/` — 15 cleaned scripts (portable paths, no internal references,
  DrugBank comment fixed). `organize_repo.sh` uses these automatically.
- `results/`, `scripts/` — folders with short README notes

Note: this scaffold does NOT contain your ~33 other scripts — those are on your
machine. `organize_repo.sh` pulls them in.

## Steps

### 1. Extract the scaffold into your project root
```bash
cd ~/model_building_strategy        # wherever your scripts/ folder is
mkdir -p repo
tar xzf offtarget_repo_scaffold.tar.gz -C repo
```

### 2. Run the organizer
```bash
bash repo/organize_repo.sh
```
This copies every script into its stage folder (01_dataset, 02_library, ...), using the
corrected version wherever one exists. Your original `scripts/` folder is left untouched.

### 3. Add your small result files (recommended)
```bash
cp dataset_final_nonredundant.csv offtarget_hits_high.tsv repo/results/
cp phase5_recall.tsv phase5_enrichment_null.tsv phase5_pathways.tsv repo/results/
cp self_recovery_summary.tsv blast_hits.tsv foldseek_hits.tsv repo/results/
```
Small tables are fine to commit; large data is blocked by `.gitignore`.

### 4. Add the figure scripts if needed
```bash
cp make_thesis_figures.py make_workflow_figure.py repo/scripts/figures/
```

### 5. Initialize Git and commit
```bash
cd repo
git init
git add .
git status          # check: no gigabyte files should be listed
git commit -m "Initial commit: proteome-wide off-target pipeline"
```

### 6. Create the GitHub repo and push
On github.com, click **New repository**, give it a name, and do NOT let GitHub add a
README (you already have one). Then:
```bash
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO.git
git push -u origin main
```

## Safety check before pushing

Run `git status` (step 5) and confirm you do NOT see any of these:
`af_by_accession/`, `alphafold_structures/`, `library_*db*`, `chembl_cache/`,
`*.tar.gz`. If you do, the `.gitignore` isn't being applied — stop and check you're in
the `repo/` folder. Committing gigabytes of structures is the one mistake to avoid.
