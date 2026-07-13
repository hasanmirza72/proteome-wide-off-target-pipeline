#!/usr/bin/env python3
"""
summarize_self_recovery.py
--------------------------
Reads the per-target FoldDisco hit files produced by run_self_recovery.sh and,
for each target, finds the rank at which the target's OWN AlphaFold model
(its self_uniprot) appears in its own pocket query. Emits a pass/fail table
for holo and apo.

Pass logic (three tiers, so you keep the nuance from the progress report):
  PASS      self found within the top PASS_RANK
  MARGINAL  self found, but ranked below PASS_RANK (still recovered, just not top)
  FAIL      self not present in the hit list at all (or no/empty output)

Config via env vars (optional):
  PASS_RANK   integer, default 5
  MANIFEST    default queries_manifest.tsv
  HITDIR      default hits
  OUT         default self_recovery_summary.tsv

Assumes the hit files were written with --header (they are, by the driver).
Rank is 1-based over the data rows, in whatever order folddisco sorted them.
"""

import os
import re
import sys

# === PATCH (abspath index) ===================================================
# The index was built with `--id abspath`, so FoldDisco's `tid` is now the full
# path of the structure file (e.g. /proj/.../AF-P00374-F1-model_v6.pdb), NOT a
# bare UniProt accession. This helper turns any tid back into the accession so
# the `tid == self_uniprot` comparison still works.
_AF_RE = re.compile(r"AF-([A-Za-z0-9]+)-F\d+-model")
def tid_to_uniprot(tid):
    base = os.path.basename(tid.strip())
    m = _AF_RE.search(base)
    if m:
        return m.group(1)                      # AF-<ACC>-F1-model_v6.pdb -> <ACC>
    for ext in (".pdb.gz", ".cif.gz", ".pdb", ".cif"):
        if base.endswith(ext):
            return base[: -len(ext)]           # <ACC>.pdb -> <ACC>
    return base                                # already a bare accession
# =============================================================================

MANIFEST  = os.environ.get("MANIFEST", "queries_manifest.tsv")
HITDIR    = os.environ.get("HITDIR", "hits")
OUT       = os.environ.get("OUT", "self_recovery_summary.tsv")
PASS_RANK = int(os.environ.get("PASS_RANK", "5"))


def read_manifest(path):
    rows, header = [], None
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line or line.lstrip().startswith("#"):
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def _split(line, delim):
    return line.split(delim) if delim else line.split()


def self_rank(hitfile, self_id):
    """Return dict: status, rank, idf, node_count, min_rmsd for the self row."""
    blank = {"status": "", "rank": "", "idf": "", "node_count": "", "min_rmsd": ""}
    if not os.path.exists(hitfile):
        return {**blank, "status": "NO_OUTPUT"}
    with open(hitfile) as f:
        lines = [ln.rstrip("\n") for ln in f if ln.strip()]
    if len(lines) < 2:                       # header only or empty
        return {**blank, "status": "EMPTY"}

    delim = "\t" if "\t" in lines[0] else None
    header = _split(lines[0], delim)
    idx = {name: i for i, name in enumerate(header)}
    # our --format-output leads with tid, so column 0 is a safe fallback
    ti = idx.get("tid", 0)

    for rank, ln in enumerate(lines[1:], start=1):
        p = _split(ln, delim)
        if ti >= len(p):
            continue
        tid = p[ti]
        # === PATCH (abspath index): compare accession parsed from the path ===
        if tid_to_uniprot(tid) == self_id:
            def get(col):
                j = idx.get(col)
                return p[j] if (j is not None and j < len(p)) else ""
            return {
                "status": "FOUND",
                "rank": rank,
                "idf": get("idf"),
                "node_count": get("node_count"),
                "min_rmsd": get("min_rmsd"),
            }
    return {**blank, "status": "NOT_IN_TOP"}


def verdict(info):
    if info["status"] != "FOUND":
        return "FAIL"
    return "PASS" if int(info["rank"]) <= PASS_RANK else "MARGINAL"


def main():
    if not os.path.exists(MANIFEST):
        sys.exit(f"Manifest not found: {MANIFEST}")
    rows = read_manifest(MANIFEST)

    cols = ["target", "self_uniprot",
            "holo_rank", "holo_nodes", "holo_rmsd", "holo_idf", "holo_verdict",
            "apo_rank",  "apo_nodes",  "apo_rmsd",  "apo_idf",  "apo_verdict"]
    out_lines = ["\t".join(cols)]

    n = len(rows)
    n_holo_pass = n_apo_pass = 0
    for r in rows:
        target  = r.get("target_id") or r.get("target") or ""
        self_id = r.get("self_uniprot", "")

        h = self_rank(os.path.join(HITDIR, f"{target}_holo.tsv"), self_id)
        a = self_rank(os.path.join(HITDIR, f"{target}_apo.tsv"),  self_id)
        hv, av = verdict(h), verdict(a)
        n_holo_pass += (hv == "PASS")
        n_apo_pass  += (av == "PASS")

        # show status word instead of a blank rank when the self was not found
        h_rank = h["rank"] if h["status"] == "FOUND" else h["status"]
        a_rank = a["rank"] if a["status"] == "FOUND" else a["status"]

        out_lines.append("\t".join(str(x) for x in [
            target, self_id,
            h_rank, h["node_count"], h["min_rmsd"], h["idf"], hv,
            a_rank, a["node_count"], a["min_rmsd"], a["idf"], av,
        ]))

    with open(OUT, "w") as f:
        f.write("\n".join(out_lines) + "\n")

    print(f"Wrote {OUT}  (PASS threshold: self within top {PASS_RANK})")
    print(f"HOLO self-recovery PASS: {n_holo_pass}/{n}")
    print(f"APO  self-recovery PASS: {n_apo_pass}/{n}")
    print()
    # pretty-print to the terminal
    widths = [max(len(str(row.split(chr(9))[i])) for row in out_lines)
              for i in range(len(cols))]
    for row in out_lines:
        cells = row.split("\t")
        print("  ".join(c.ljust(widths[i]) for i, c in enumerate(cells)))


if __name__ == "__main__":
    main()
