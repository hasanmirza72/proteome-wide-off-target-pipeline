#!/usr/bin/env python3
"""
offtarget_overlap_filter.py  (roadmap Part 3.3 + 3.4)
-----------------------------------------------------
Joins FoldDisco off-target hits to the P2Rank library_pockets table and asks,
for each hit: do the residues FoldDisco matched actually sit inside a
P2Rank-predicted cavity on that target protein?

  - pocket_supported : matched residues overlap a P2Rank pocket  (higher confidence)
  - non_pocket       : match lands on buried / non-pocket surface (flag, keep, lower)
  - is_self          : the query protein's own AF model (the on-target; excluded
                       from the reported off-target list, per Part 3.4)

Inputs : queries_manifest.tsv, library_pockets.tsv, hits_offtarget/<target>_{holo,apo}.tsv
Outputs: offtarget_hits_annotated.tsv (everything, auditable)
         offtarget_hits_filtered.tsv  (self removed, pocket_supported only)
"""

import csv
import os
import re
import sys

# --- Configuration -----------------------------------------------------------
MANIFEST   = "queries_manifest.tsv"
LIBRARY    = "library_pockets.tsv"
HITDIR     = "hits_offtarget"
OUT_ALL    = "offtarget_hits_annotated.tsv"
OUT_FILT   = "offtarget_hits_filtered.tsv"
OUT_HIGH   = "offtarget_hits_high.tsv"     # self removed, HIGH confidence only (Phase 4 input)
MIN_OVERLAP_RES = 1     # >= this many matched residues inside a pocket => supported

# === Confidence tiering (added for Phase 4) ==================================
# A hit is only as trustworthy as its geometry. `pocket_supported` alone admits
# 2-residue, single-overlap, low-IDF junk, so we layer a high/medium/low tier on
# top of it. HIGH is the set you take into the enrichment test.
#   HIGH   : substantial motif AND real pocket overlap AND tight geometry AND rare
#   MEDIUM : pocket_supported but not all HIGH bars met
#   LOW    : pocket_supported but weak (small motif / 1-residue overlap)
#   none   : is_self or non_pocket (blank confidence)
HIGH_MIN_NODES   = 4      # matched motif size (folddisco node_count)
HIGH_MIN_OVERLAP = 3      # matched residues inside the P2Rank pocket
HIGH_MIN_FRAC    = 0.30   # fraction of matched residues inside the pocket
HIGH_MAX_RMSD    = 2.0    # superposition RMSD ceiling (A)
HIGH_MIN_IDF     = 3.0    # rarity floor (absolute IDF)
LOW_MAX_NODES    = 3      # <= this with only 1 overlap => LOW
# -----------------------------------------------------------------------------


RES_RE = re.compile(r"^([A-Za-z]+)(-?\d+)")

# === PATCH (abspath index) ===================================================
# With `--id abspath`, a hit's `tid` is the full structure path, not a bare
# UniProt accession. Turn it back into the accession for the library join.
_AF_RE = re.compile(r"AF-([A-Za-z0-9]+)-F\d+-model")
def tid_to_uniprot(tid):
    base = os.path.basename(tid.strip())
    m = _AF_RE.search(base)
    if m:
        return m.group(1)                      # AF-<ACC>-F1-model_v6.pdb -> <ACC>
    for ext in (".pdb.gz", ".cif.gz", ".pdb", ".cif"):
        if base.endswith(ext):
            return base[: -len(ext)]           # <ACC>.pdb -> <ACC>
    return base
# =============================================================================


def res_key(tok):
    # === PATCH: real folddisco output tags residues as "A71:0.31" -> drop :rmsd
    m = RES_RE.match(tok.split(":")[0].strip())
    return (m.group(1), int(m.group(2))) if m else None


def res_set(field):
    # === PATCH: split on BOTH ',' and ';' (multi-match rows use ';' between
    # sub-matches); '_' placeholder tokens simply don't match and are dropped.
    return {k for k in (res_key(t) for t in re.split(r"[;,]", field)) if k}


def load_library(path):
    """uniprot -> list of dicts(rank, prob, residues:set)."""
    lib = {}
    if not os.path.exists(path):
        sys.exit(f"library not found: {path} (run build_library_pockets.py first)")
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            acc = row["uniprot"].strip()
            lib.setdefault(acc, []).append({
                "rank": row.get("pocket_rank", "").strip(),
                "prob": row.get("probability", "").strip(),
                "residues": res_set(row.get("residues", "")),
            })
    return lib


def load_manifest(path):
    self_of = {}
    with open(path) as f:
        r = csv.DictReader(f, delimiter="\t")
        for row in r:
            self_of[row["target_id"].strip()] = row.get("self_uniprot", "").strip()
    return self_of


def read_hits(path):
    """Yield dicts keyed by folddisco column names (assumes --header)."""
    if not os.path.exists(path):
        return
    with open(path) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    if len(lines) < 2:
        return
    delim = "\t" if "\t" in lines[0] else None
    header = lines[0].split(delim) if delim else lines[0].split()
    for ln in lines[1:]:
        parts = ln.split(delim) if delim else ln.split()
        yield {header[i]: (parts[i] if i < len(parts) else "") for i in range(len(header))}


def best_overlap(matched, pockets):
    """Return (best_count, best_frac, best_rank, best_prob) over a protein's pockets."""
    best = (0, 0.0, "", "")
    denom = len(matched) or 1
    for p in pockets:
        c = len(matched & p["residues"])
        if c > best[0]:
            best = (c, c / denom, p["rank"], p["prob"])
    return best


def _to_float(x, default):
    try:
        return float(x)
    except (TypeError, ValueError):
        return default


def confidence_tier(status, node_count, overlap_cnt, overlap_frac, min_rmsd, idf):
    """high / medium / low for pocket_supported hits; '' otherwise."""
    if status != "pocket_supported":
        return ""
    nodes = _to_float(node_count, 0)
    rmsd  = _to_float(min_rmsd, 99.0)     # missing rmsd => treat as bad
    idfv  = _to_float(idf, 0)
    if (nodes >= HIGH_MIN_NODES and overlap_cnt >= HIGH_MIN_OVERLAP
            and overlap_frac >= HIGH_MIN_FRAC and rmsd <= HIGH_MAX_RMSD
            and idfv >= HIGH_MIN_IDF):
        return "high"
    if nodes <= LOW_MAX_NODES and overlap_cnt <= 1:
        return "low"
    return "medium"


def main():
    self_of = load_manifest(MANIFEST)
    lib = load_library(LIBRARY)

    cols = ["query_target", "query_conf", "hit_uniprot", "is_self",
            "idf", "node_count", "min_rmsd", "plddt", "n_matched",
            "best_pocket_rank", "best_pocket_prob", "overlap_count",
            "overlap_frac", "status", "confidence"]
    all_rows = ["\t".join(cols)]
    filt_rows = ["\t".join(cols)]
    high_rows = ["\t".join(cols)]

    per_query = {}   # (target,conf) -> [n_total, n_self, n_supported, n_high]

    for target, self_acc in self_of.items():
        for conf in ("holo", "apo"):
            hitfile = os.path.join(HITDIR, f"{target}_{conf}.tsv")
            key = (target, conf)
            per_query.setdefault(key, [0, 0, 0, 0])
            for h in read_hits(hitfile):
                tid = h.get("tid", "").strip()
                if not tid:
                    continue
                acc = tid_to_uniprot(tid)          # === PATCH (abspath index) ===
                per_query[key][0] += 1
                is_self = (acc == self_acc)         # compare accessions, not paths
                matched = res_set(h.get("matching_residues", ""))
                cnt, frac, prank, pprob = best_overlap(matched, lib.get(acc, []))
                supported = cnt >= MIN_OVERLAP_RES
                if is_self:
                    status = "is_self"
                    per_query[key][1] += 1
                elif supported:
                    status = "pocket_supported"
                    per_query[key][2] += 1
                else:
                    status = "non_pocket"

                tier = confidence_tier(status, h.get("node_count", ""),
                                       cnt, frac, h.get("min_rmsd", ""),
                                       h.get("idf", ""))
                if tier == "high":
                    per_query[key][3] += 1

                row = "\t".join(str(x) for x in [
                    target, conf, acc, int(is_self),   # === PATCH: accession, not path
                    h.get("idf", ""), h.get("node_count", ""),
                    h.get("min_rmsd", ""), h.get("plddt", ""), len(matched),
                    prank, pprob, cnt, f"{frac:.2f}", status, tier,
                ])
                all_rows.append(row)
                if (not is_self) and supported:
                    filt_rows.append(row)
                if (not is_self) and tier == "high":
                    high_rows.append(row)

    with open(OUT_ALL, "w") as f:
        f.write("\n".join(all_rows) + "\n")
    with open(OUT_FILT, "w") as f:
        f.write("\n".join(filt_rows) + "\n")
    with open(OUT_HIGH, "w") as f:
        f.write("\n".join(high_rows) + "\n")

    print(f"Wrote {OUT_ALL} ({len(all_rows)-1} rows), "
          f"{OUT_FILT} ({len(filt_rows)-1} pocket-supported), "
          f"{OUT_HIGH} ({len(high_rows)-1} HIGH-confidence)\n")
    print(f"{'query':<14}{'total':>7}{'self':>6}{'pocket-ok':>11}{'high':>6}")
    for (target, conf), (tot, slf, sup, hi) in sorted(per_query.items()):
        if tot:
            print(f"{target+'/'+conf:<14}{tot:>7}{slf:>6}{sup:>11}{hi:>6}")


if __name__ == "__main__":
    main()
