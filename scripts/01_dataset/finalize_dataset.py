#!/usr/bin/env python3
"""
FINALIZE: one script that produces the definitive dataset + evaluation.

It starts from apo_holo_pairs.csv and, in documented stages:
  A. curation filters  - non-drug/small ligands, weak apo (low TM / low mapped%)
  B. scorability       - a pair is kept only if BOTH holo & apo P2Rank predictions
                         exist, the drug has atoms, and (if PyMOL is available) the
                         holo ligand actually contacts the protein (real pocket)
  C. metrics           - DCA strict(4.0)/loose(5.5), best DCA/DCC, closest rank,
                         detected-but-mis-ranked, and byres binding-residue recall
  D. statistics        - Top-3 rates with Wilson 95% CI, exact McNemar (paired),
                         discordant-pair listing, on BOTH the full set and a
                         non-redundant one-per-UniProt subset

Outputs:
  dataset_final.csv              <- THE dataset to use downstream (scorable pairs + metrics)
  dataset_final_nonredundant.csv <- one representative per UniProt
  dataset_excluded.csv           <- every dropped pair, with stage + reason
  final_report.txt               <- the numbers/text for the thesis (also printed)

Recall needs PyMOL; if it is not importable, all detection stats still run and the
recall columns become NA (they are not required for the headline result).
"""

import csv, glob, math, os

# ---------------- config ----------------
RAW_PAIRS = "apo_holo_pairs.csv"
HOLO_PDB_DIR = "dataset_monomers_pdb"
APO_DIR = "dataset_apo_pdb"
HOLO_PRED_DIR = "p2rank_out_holo"
APO_PRED_DIR = "p2rank_out_apo"

DCA_STRICT, DCA_LOOSE, TOP_N, CONTACT_DIST = 4.0, 5.5, 3, 4.0
MIN_LIG_ATOMS, TM_MIN, MAPPED_MIN = 6, 0.50, 80.0
BLACKLIST_LIG = {"OH", "CO3", "GLY", "ACT", "EDO", "GOL", "SO4", "PO4", "NO3",
                 "FLC", "DMS", "IPA", "MPD", "CL", "NA", "ZN", "MG", "CA", "K"}

FINAL_CSV = "dataset_final.csv"
FINAL_NR_CSV = "dataset_final_nonredundant.csv"
EXCLUDED_CSV = "dataset_excluded.csv"
REPORT_TXT = "final_report.txt"

try:
    import pymol
    from pymol import cmd
    pymol.finish_launching(['pymol', '-cq'])
    _PYMOL = True
except Exception:
    _PYMOL = False


# ---------------- helpers ----------------
def fnum(v, d=None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return d


def holo_path(holo):
    for c in (f"{holo}.pdb", f"{holo.lower()}.pdb"):
        p = os.path.join(HOLO_PDB_DIR, c)
        if os.path.exists(p):
            return p
    return None


def read_ligand_atoms(path, lig):
    lig = (lig or "").strip().upper()
    out = []
    if not path:
        return out
    with open(path, encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith(("HETATM", "ATOM")) and line[17:20].strip().upper() == lig:
                try:
                    out.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                except ValueError:
                    pass
    return out


def true_binding_residues(path, ligand_atoms):
    if not (_PYMOL and ligand_atoms and path and os.path.exists(path)):
        return None  # None = uncomputable (PyMOL off), set() = computed-but-empty
    cmd.delete("all")
    cmd.load(path, "prot")
    for (x, y, z) in ligand_atoms:
        cmd.pseudoatom("ghost", pos=[x, y, z])
    cmd.select("bs", f"byres (polymer within {CONTACT_DIST} of ghost)")
    space = {"res": set()}
    cmd.iterate("bs and name CA", "res.add(f'{chain}_{resi}')", space=space)
    return space["res"]


def read_pockets(pred_csv):
    pockets = []
    with open(pred_csv, encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.reader(f)
        header = [h.strip() for h in next(reader)]
        col = {n: i for i, n in enumerate(header)}
        if not all(k in col for k in ("rank", "center_x", "center_y", "center_z")):
            return []
        ridx = col.get("residue_ids")
        for row in reader:
            if not row or not row[col["rank"]].strip():
                continue
            try:
                pockets.append({
                    "rank": int(float(row[col["rank"]])),
                    "center": (float(row[col["center_x"]]), float(row[col["center_y"]]),
                               float(row[col["center_z"]])),
                    "residues": (row[ridx].split() if ridx is not None and ridx < len(row) else []),
                })
            except (ValueError, IndexError):
                continue
    return sorted(pockets, key=lambda p: p["rank"])


def find_pred_file(pred_dir, pdb_id):
    pid = pdb_id.lower()
    cands = [p for p in glob.glob(os.path.join(pred_dir, "*_predictions.csv"))
             if pid in os.path.basename(p).lower()]
    return cands[0] if cands else None


def dist(a, b):
    return math.sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2 + (a[2]-b[2])**2)


def centroid(atoms):
    n = len(atoms)
    return (sum(a[0] for a in atoms)/n, sum(a[1] for a in atoms)/n, sum(a[2] for a in atoms)/n)


def score(pockets, ligand_atoms, true_res):
    r = {"hit_strict": 0, "hit_loose": 0, "best_dca": None, "best_dcc": None,
         "closest_rank": None, "res_recall": None,
         "n_true_res": (len(true_res) if isinstance(true_res, set) else None)}
    if not pockets or not ligand_atoms:
        return r
    cen = centroid(ligand_atoms)
    bd = bc = math.inf
    closest = None
    for p in pockets:
        dca = min(dist(p["center"], a) for a in ligand_atoms)
        dcc = dist(p["center"], cen)
        bd, bc = min(bd, dca), min(bc, dcc)
        if closest is None or dca < closest["_dca"]:
            closest = dict(p, _dca=dca)
        if dca <= DCA_STRICT and r["hit_strict"] == 0:
            r["hit_strict"] = p["rank"]
        if dca <= DCA_LOOSE and r["hit_loose"] == 0:
            r["hit_loose"] = p["rank"]
    r["best_dca"], r["best_dcc"] = round(bd, 2), round(bc, 2)
    r["closest_rank"] = closest["rank"]
    if isinstance(true_res, set) and true_res and closest["residues"]:
        r["res_recall"] = round(len(true_res & set(closest["residues"])) / len(true_res), 3)
    return r


def mcnemar_exact(b, c):
    n = b + c
    if n == 0:
        return 1.0
    k = min(b, c)
    return min(1.0, 2 * sum(math.comb(n, i) for i in range(k + 1)) * 0.5 ** n)


def wilson(k, n, z=1.96):
    if n == 0:
        return (0.0, 0.0)
    p = k / n
    d = 1 + z*z/n
    c = (p + z*z/(2*n)) / d
    h = z/d * math.sqrt(p*(1-p)/n + z*z/(4*n*n))
    return (round(100*(c-h), 1), round(100*(c+h), 1))


# ---------------- statistics over a set of scored rows ----------------
def compute_stats(rows):
    n = len(rows)
    out = {"n": n}
    for thr in ("strict", "loose"):
        field = f"hit_{thr}"
        h = [1 if 1 <= r[f"holo_{field}"] <= TOP_N else 0 for r in rows]
        a = [1 if 1 <= r[f"apo_{field}"] <= TOP_N else 0 for r in rows]
        hs, as_ = sum(h), sum(a)
        b = sum(1 for i in range(n) if h[i] == 1 and a[i] == 0)
        c = sum(1 for i in range(n) if h[i] == 0 and a[i] == 1)
        out[thr] = {
            "holo": hs, "apo": as_,
            "holo_ci": wilson(hs, n), "apo_ci": wilson(as_, n),
            "holo_only": b, "apo_only": c, "p": mcnemar_exact(b, c),
            "holo_only_pairs": [rows[i]["holo"] for i in range(n) if h[i] == 1 and a[i] == 0],
            "apo_only_pairs": [rows[i]["holo"] for i in range(n) if h[i] == 0 and a[i] == 1],
        }
    for tag in ("holo", "apo"):
        misrank = [r["holo"] for r in rows
                   if r[f"{tag}_hit_strict"] > TOP_N and fnum(r[f"{tag}_best_dca"], 99) <= DCA_STRICT]
        rec = [r[f"{tag}_res_recall"] for r in rows
               if isinstance(r.get(f"{tag}_res_recall"), (int, float))
               and 1 <= r[f"{tag}_hit_strict"] <= TOP_N]
        out[f"{tag}_misrank"] = misrank
        out[f"{tag}_recall"] = (round(sum(rec)/len(rec), 2), len(rec)) if rec else (None, 0)
    return out


def fmt_block(title, s):
    L = [f"### {title}   (paired n = {s['n']})"]
    for thr, lbl in (("strict", f"strict DCA<={DCA_STRICT}A"), ("loose", f"loose DCA<={DCA_LOOSE}A")):
        d = s[thr]
        L.append(f"  {lbl}:")
        L.append(f"    holo Top-{TOP_N}: {d['holo']}/{s['n']} = {100*d['holo']/s['n']:.1f}%  95%CI{d['holo_ci']}")
        L.append(f"    apo  Top-{TOP_N}: {d['apo']}/{s['n']} = {100*d['apo']/s['n']:.1f}%  95%CI{d['apo_ci']}")
        L.append(f"    discordant holo-only={d['holo_only']} {d['holo_only_pairs']}, "
                 f"apo-only={d['apo_only']} {d['apo_only_pairs']}  McNemar p={d['p']:.3f}")
    for tag in ("holo", "apo"):
        mr, k = s[f"{tag}_recall"]
        rline = "recall NA (PyMOL off)" if mr is None else f"mean binding-residue recall={mr} (n={k})"
        L.append(f"  {tag}: {rline}; detected-but-mis-ranked: {s[f'{tag}_misrank'] or 'none'}")
    return "\n".join(L)


# ---------------- main ----------------
def main():
    if not os.path.exists(RAW_PAIRS):
        print(f"Error: '{RAW_PAIRS}' not found."); return
    with open(RAW_PAIRS, encoding="utf-8", errors="ignore", newline="") as f:
        reader = csv.DictReader(f)
        base_fields = reader.fieldnames
        raw = list(reader)

    excluded = []       # (holo, apo, stage, reason)
    scored = []         # final scorable rows (dict with base cols + metrics)

    for p in raw:
        holo, apo, lig = p["holo_pdb"].strip(), p["apo_pdb"].strip(), p["holo_ligand"].strip()
        hp = holo_path(holo)
        ligand_atoms = read_ligand_atoms(hp, lig)
        n_lig = len(ligand_atoms)

        # ---- Stage A: curation ----
        tm = fnum(p.get("tm_score"))
        mapped = fnum(p.get("mapped_binding_residues_percent"))
        reasons = []
        if lig.upper() in BLACKLIST_LIG:
            reasons.append(f"A:non-drug ligand '{lig}'")
        if hp and n_lig < MIN_LIG_ATOMS:
            reasons.append(f"A:ligand too small ({n_lig} atoms)")
        if tm is not None and tm < TM_MIN:
            reasons.append(f"A:low TM ({tm})")
        if mapped is not None and mapped < MAPPED_MIN:
            reasons.append(f"A:low mapped% ({mapped})")
        if reasons:
            excluded.append((holo, apo, "curation", "; ".join(reasons))); continue

        # ---- Stage B: scorability ----
        if hp is None:
            excluded.append((holo, apo, "scorability", "B:holo PDB file missing")); continue
        if n_lig == 0:
            excluded.append((holo, apo, "scorability", "B:no ligand atoms (name mismatch)")); continue
        holo_pf = find_pred_file(HOLO_PRED_DIR, holo)
        apo_pf = find_pred_file(APO_PRED_DIR, apo)
        if not holo_pf:
            excluded.append((holo, apo, "scorability", "B:no holo prediction file")); continue
        if not apo_pf:
            excluded.append((holo, apo, "scorability", "B:no apo prediction file")); continue

        apo_file = p.get("apo_file", "")
        apo_full = apo_file if os.path.exists(apo_file) else os.path.join(APO_DIR, os.path.basename(apo_file))

        holo_true = true_binding_residues(hp, ligand_atoms)
        # real-pocket check (only enforce when PyMOL could actually compute it)
        if isinstance(holo_true, set) and len(holo_true) == 0:
            excluded.append((holo, apo, "scorability", "B:ligand contacts 0 residues (surface/crystal artifact)"))
            continue
        apo_true = true_binding_residues(apo_full, ligand_atoms)

        hm = score(read_pockets(holo_pf), ligand_atoms, holo_true)
        am = score(read_pockets(apo_pf), ligand_atoms, apo_true)

        rec = dict(p)
        rec.update({"holo": holo, "apo": apo, "ligand": lig, "n_lig_atoms": n_lig})
        for tag, m in (("holo", hm), ("apo", am)):
            for k, v in m.items():
                rec[f"{tag}_{k}"] = "" if v is None else v
        scored.append(rec)

    # ---- non-redundant subset (best resolution per UniProt) ----
    def better(a, b):
        ra, rb = fnum(a.get("resolution"), 1e9), fnum(b.get("resolution"), 1e9)
        if ra != rb:
            return a if ra < rb else b
        return a if fnum(a.get("mapped_binding_residues_percent"), 0) >= \
            fnum(b.get("mapped_binding_residues_percent"), 0) else b
    rep = {}
    for r in scored:
        key = (r.get("uniprot") or "").strip() or f"__{r['holo']}"
        rep[key] = r if key not in rep else better(rep[key], r)
    nonredundant = list(rep.values())

    # ---- write outputs ----
    metric_cols = ["holo", "apo", "ligand", "n_lig_atoms",
                   "holo_hit_strict", "holo_hit_loose", "holo_best_dca", "holo_best_dcc",
                   "holo_closest_rank", "holo_res_recall", "holo_n_true_res",
                   "apo_hit_strict", "apo_hit_loose", "apo_best_dca", "apo_best_dcc",
                   "apo_closest_rank", "apo_res_recall", "apo_n_true_res"]
    extra = [c for c in base_fields if c not in metric_cols]
    fields = metric_cols + extra
    for path, data in ((FINAL_CSV, scored), (FINAL_NR_CSV, nonredundant)):
        with open(path, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
            w.writeheader()
            w.writerows(data)
    with open(EXCLUDED_CSV, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["holo_pdb", "apo_pdb", "stage", "reason"])
        w.writerows(excluded)

    # ---- report ----
    holo_ids = {r["holo"].upper() for r in scored}
    apo_ids = {r["apo"].upper() for r in scored}
    both = sorted(holo_ids & apo_ids)
    uniq = len({(r.get("uniprot") or r["holo"]).strip() for r in scored})

    lines = []
    lines.append("=" * 78)
    lines.append("   FINAL DATASET & EVALUATION REPORT")
    lines.append("=" * 78)
    lines.append(f"PyMOL available (residue recall): {_PYMOL}")
    lines.append(f"input pairs: {len(raw)}   excluded: {len(excluded)}   "
                 f"FINAL scorable: {len(scored)}   unique proteins: {uniq}")
    lines.append("")
    lines.append("Exclusions (also in dataset_excluded.csv):")
    if excluded:
        for h, a, stage, why in excluded:
            lines.append(f"   - {h} -> {a}  [{stage}]  {why}")
    else:
        lines.append("   none")
    lines.append("")
    lines.append("Role audit: " + ("[!] BOTH roles: " + ", ".join(both) if both
                                    else "clean (no PDB is both holo and apo)"))
    lines.append("")
    lines.append(fmt_block("FULL FINAL SET", compute_stats(scored)))
    lines.append("")
    lines.append(fmt_block("NON-REDUNDANT (one per UniProt)", compute_stats(nonredundant)))
    lines.append("")
    lines.append("Recall definition: fraction of that structure's OWN binding residues "
                 "(any protein atom <=4A of the drug, byres) recovered by the closest "
                 "predicted pocket. Holo residues from the holo frame; apo residues from "
                 "the AHoJ-aligned apo frame using the transferred holo ligand.")
    lines.append("")
    lines.append(f"Files: {FINAL_CSV} (use downstream), {FINAL_NR_CSV}, {EXCLUDED_CSV}")
    report = "\n".join(lines)
    print(report)
    with open(REPORT_TXT, "w", encoding="utf-8") as f:
        f.write(report + "\n")

    if _PYMOL:
        cmd.quit()


if __name__ == "__main__":
    main()
