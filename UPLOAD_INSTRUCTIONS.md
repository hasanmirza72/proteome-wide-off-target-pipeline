# Uploading the updated repo to GitHub (clean rebuild + force push)

Your clean-slate approach is correct. Here is the corrected full sequence, with the
fixes noted. Run these in your Berzelius terminal.

## 1. Go to your project root and delete the old repo folder

```bash
cd ~/model_building_strategy
rm -rf repo
```

## 2. Extract the new scaffold

```bash
mkdir -p repo
tar xzf offtarget_repo_scaffold.tar.gz -C repo
```

## 3. Run the organizer FROM THE PROJECT ROOT (not from inside repo/)

The organizer reads your original `scripts/` folder (at the project root) and copies
into `repo/scripts/`. So stay in the project root when you run it:

```bash
bash repo/organize_repo.sh
```

If you see many `[skip] ... not found` lines, you are in the wrong folder — make sure
you are in `~/model_building_strategy`, where your original `scripts/` folder lives.

## 4. Copy result files, supplementary data, and figures

The scaffold already includes the 10 figures and figure scripts, so those are done.
You still need to add the small result tables, the supplementary data, and the figures
(if the organizer did not find them). Run all of this from the project root:

```bash
# result tables
cp dataset_final_nonredundant.csv offtarget_hits_high.tsv \
   phase5_recall.tsv phase5_enrichment_null.tsv phase5_pathways.tsv \
   self_recovery_summary.tsv blast_hits.tsv foldseek_hits.tsv \
   repo/results/ 2>/dev/null

# supplementary data (complete the set — the scaffold already has 9 of these)
cp dataset_final_nonredundant.csv dataset_excluded.csv p2rank_evaluation.csv \
   self_recovery_summary.tsv dual_query_divergence.tsv redundant_robustness.tsv \
   phase5_enrichment_null.tsv phase5_pathways.tsv \
   dataset_final.csv verified_alphafold_targets.csv library_pockets.tsv \
   phase5_summary.txt \
   repo/supplementary_data/ 2>/dev/null

# figures (only needed if you regenerated them; scaffold already has them)
cp repo/scripts/figures/figures/*.png repo/scripts/figures/figures/ 2>/dev/null
```

## 5. (Optional) delete the staging folders so the repo is tidy

The `corrected_scripts/` folder and `organize_repo.sh` have done their job once step 3
has run. You can remove them before committing:

```bash
rm -rf repo/corrected_scripts repo/organize_repo.sh
```

## 6. Initialize Git and CHECK STATUS — the critical safety step

```bash
cd repo
git init
git add .
git status
```

STOP AND LOOK at the list. You should see scripts, figures, supplementary_data, README,
docs. You should NOT see any of these:
`alphafold_structures/`, `af_by_accession/`, `chembl_cache/`, `library_blastdb.*`,
`folddisco_af_index*`, `*.tar.gz`. If any giant folder appears, the `.gitignore` is not
being applied — stop and check you are inside `repo/` (where `.gitignore` lives).

Quick sanity check on how many files are staged (should be a few dozen, not thousands):

```bash
git status --short | wc -l
```

## 7. Commit

```bash
git commit -m "Rebuild: pipeline, figures, and supplementary data"
```

## 8. Force-push to overwrite the old version on GitHub

```bash
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git
git push -f -u origin main
```

Notes:
- If `git remote add origin` says "remote origin already exists", set the URL instead:
  `git remote set-url origin https://github.com/YOUR_USERNAME/YOUR_REPO_NAME.git`
- When prompted for a password, paste your Personal Access Token (PAT), not your
  account password.
- `-f` (force) is what lets the new version overwrite the old one on GitHub.

You are done when you see `Writing objects: 100%` and `main -> main`.

## Alternative to force-push (safer, if you prefer)

Force-push discards the old history. If you would rather keep it, you can instead pull
and merge first — but for a solo thesis repo where you just want the new version to
win, the force-push above is simplest and fine.
