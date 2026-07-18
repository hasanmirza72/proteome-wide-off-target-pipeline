#!/usr/bin/env python3
"""
score_baselines.py — compute the three-method baseline comparison.

Reproduces the documented off-target recall comparison between FoldDisco, BLAST, and
Foldseek that is reported in the thesis, at common rank cutoffs, and the split of that
recall by the off-target's structural relationship to the query pocket.

Inputs (from the project root):
    folddisco_raw_hits.tsv          FoldDisco raw ranked hits (from extract_raw_ranks.py)
    blast_hits_corrected.tsv        BLAST+ hits    (from run_blast.sh, canonical queries)
    foldseek_hits_corrected.tsv     Foldseek hits  (from run_foldseek.sh, canonical queries)

Outputs (all also printed):
    self_recovery_by_method.tsv     self-recovery of all 23 targets, per method
    self_recovery_hard_targets.tsv  self-rank of the seven detection-hard targets
    baseline_comparison.tsv         documented off-target recall at three cutoffs
    baseline_split.tsv              that recall split by structural relationship

The documented off-target ground truth is the same set used by phase5_recall.py, and is
independent of anything the pipeline retrieved. The self-hit (the query protein's own
model) is excluded from every method's list before ranking.

Run:  python3 score_baselines.py
"""
import re, csv, statistics
from collections import defaultdict

FD    = "folddisco_raw_hits.tsv"
BLAST = "blast_hits.tsv"
FSEEK = "foldseek_hits.tsv"
CUTOFFS = (10, 30, 100)
CONF = "holo"     # query conformation used for FoldDisco; baselines are holo-derived

# query PDB -> query UniProt accession (all 23 non-redundant benchmark targets)
QUERY = {"6VCJ":"P00374","1UEI":"Q9BZX2","1NHZ":"P04150","6L6E":"O76074","4U7Z":"Q02750",
         "1XBB":"P43405","3TVX":"P27815","3RX3":"P15121","5V1M":"Q9BQ65","6Q6O":"P51449",
         "6GH9":"Q9Y4E8","8HJE":"Q96RU2","2XP2":"Q9UM73","2Y05":"Q14914","7AXA":"O75469",
         "1AVN":"P00918","7XPY":"Q93009","6YG2":"O14733","4ODR":"P62942","3TC5":"Q13526",
         "1GS4":"P10275","4R7L":"P09960","9XRI":"P36639"}

# the seven targets whose pocket P2Rank fails to detect (Chapter 3), for Table 10
HARD = [("6YG2","MEK7"), ("1NHZ","GR (NR3C1)"), ("4U7Z","MEK1"), ("6GH9","USP15"),
        ("2Y05","PTGR1"), ("7XPY","USP7"), ("4R7L","LTA4H")]

# Documented off-targets per query drug (same ground truth as phase5_recall.py),
# each tagged with its structural relationship to the query pocket.
#   pocket_similar : shares local binding-site geometry / same enzyme family
#   fold_unrelated : occupies a fold unrelated to the query pocket
KNOWN = {
 # tier: pocket_similar     = shares local binding-site geometry with the query
 #       fold_pocket_diff   = same broad family/fold, but the pocket itself differs
 #       fold_unrelated     = occupies a fold unrelated to the query pocket
 "2XP2": {"P08581":"pocket_similar","P07947":"pocket_similar","O14965":"pocket_similar"},
 "1XBB": {"P00519":"fold_pocket_diff","P10721":"fold_pocket_diff",
          "P16234":"fold_pocket_diff","Q08345":"fold_pocket_diff"},
 "4U7Z": {"P36507":"pocket_similar"},
 "6YG2": {"Q06187":"fold_unrelated","P07332":"fold_unrelated",
          "Q08881":"fold_unrelated","P00533":"fold_unrelated"},
 "3RX3": {"O60218":"pocket_similar","P35354":"fold_unrelated","P23219":"fold_unrelated"},
 "3TVX": {"Q07343":"pocket_similar","Q08499":"pocket_similar","Q14432":"pocket_similar"},
 "4ODR": {"P68106":"pocket_similar","Q02790":"pocket_similar","Q13451":"pocket_similar"},
 "2Y05": {"P03372":"fold_unrelated","Q92731":"fold_unrelated"},
 "4R7L": {"Q13547":"fold_unrelated","Q92769":"fold_unrelated","P56524":"fold_unrelated"},
}

TIER_LABEL = {"pocket_similar":   "Pocket-similar",
              "fold_pocket_diff": "Fold-related, pocket-dissimilar",
              "fold_unrelated":   "Fold-unrelated"}


def load_folddisco(path, conf):
    """query PDB -> ranked list of hit accessions (file order = rank order)."""
    d = defaultdict(list)
    with open(path) as f:
        next(f, None)
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5 or p[1] != conf:
                continue
            d[p[0]].append(p[3])
    return d


def load_tabular(path, parse_sp=False):
    """query accession -> ranked list of hit accessions (file order = rank order).

    BLAST reports one line per high-scoring segment pair, so a subject protein can
    appear several times for one query. Only the first (best-scoring) occurrence of
    each subject is kept, so that a rank counts distinct proteins rather than
    alignment segments. Foldseek reports one line per target and is unaffected.
    """
    d = defaultdict(list)
    seen = defaultdict(set)
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) < 2:
                continue
            target = p[1]
            if parse_sp:
                m = re.search(r"\|([A-Z0-9]+)\|", target)
                target = m.group(1) if m else target
            if target in seen[p[0]]:
                continue
            seen[p[0]].add(target)
            d[p[0]].append(target)
    return d


def rank_of(hits, target, self_acc):
    """1-based rank of `target` after removing the self-hit; None if absent."""
    r = 0
    for acc in hits:
        if acc == self_acc:
            continue
        r += 1
        if acc == target:
            return r
    return None


def collect(hits_by_query, key_is_pdb):
    """-> {(pdb, target): rank or None} across every documented off-target."""
    out = {}
    for pdb, targets in KNOWN.items():
        self_acc = QUERY[pdb]
        key = pdb if key_is_pdb else self_acc
        hits = hits_by_query.get(key, [])
        for t in targets:
            out[(pdb, t)] = rank_of(hits, t, self_acc)
    return out


def self_rank(hits, self_acc):
    """1-based rank of the query protein's own model in its own hit list."""
    for i, acc in enumerate(hits, start=1):
        if acc == self_acc:
            return i
    return None


def self_recovery(fd_hits, bl_hits, fs_hits):
    """-> {method: {pdb: self-rank or None}} over all 23 benchmark targets."""
    out = {"FoldDisco": {}, "BLAST": {}, "Foldseek": {}}
    for pdb, acc in QUERY.items():
        out["FoldDisco"][pdb] = self_rank(fd_hits.get(pdb, []), acc)
        out["BLAST"][pdb]     = self_rank(bl_hits.get(acc, []), acc)
        out["Foldseek"][pdb]  = self_rank(fs_hits.get(acc, []), acc)
    return out


def main():
    methods = {
        "FoldDisco": collect(load_folddisco(FD, CONF), key_is_pdb=True),
        "BLAST":     collect(load_tabular(BLAST, parse_sp=True), key_is_pdb=False),
        "Foldseek":  collect(load_tabular(FSEEK), key_is_pdb=False),
    }
    n_total = sum(len(v) for v in KNOWN.values())

    # ---- self-recovery (Table 9) ----
    fd_hits = load_folddisco(FD, CONF)
    bl_hits = load_tabular(BLAST, parse_sp=True)
    fs_hits = load_tabular(FSEEK)
    sr = self_recovery(fd_hits, bl_hits, fs_hits)

    print("=" * 78)
    print(f"SELF-RECOVERY BY SEARCH METHOD  (n = {len(QUERY)}; query conformation = {CONF})")
    print("=" * 78)
    print(f"{'Method':<12}{'Median self-rank':>18}{'Rank 1 (n)':>12}{'Top 5 (n)':>11}{'Not returned':>14}")
    sr_rows = []
    for name in ("FoldDisco", "Foldseek", "BLAST"):
        ranks = sr[name]
        got = [r for r in ranks.values() if r is not None]
        med = statistics.median(got) if got else float("nan")
        r1 = sum(1 for r in got if r == 1)
        t5 = sum(1 for r in got if r <= 5)
        miss = sum(1 for r in ranks.values() if r is None)
        print(f"{name:<12}{med:>18g}{r1:>12}{t5:>11}{miss:>14}")
        sr_rows.append([name, f"{med:g}", r1, t5, miss])
    with open("self_recovery_by_method.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Method", "Median self-rank", "Rank 1 (n)", "Top 5 (n)", "Not returned (n)"])
        w.writerows(sr_rows)

    # ---- self-rank of the detection-hard targets (Table 10) ----
    print()
    print("=" * 78)
    print("SELF-RANK OF THE SEVEN DETECTION-HARD TARGETS, BY METHOD")
    print("=" * 78)
    print(f"{'Target':<14}{'FoldDisco':>11}{'Foldseek':>10}{'BLAST':>8}")
    hard_rows = []
    for pdb, label in HARD:
        row = [label] + [sr[m][pdb] for m in ("FoldDisco", "Foldseek", "BLAST")]
        fmt = lambda v: str(v) if v is not None else "absent"
        print(f"{row[0]:<14}{fmt(row[1]):>11}{fmt(row[2]):>10}{fmt(row[3]):>8}")
        hard_rows.append([row[0]] + [fmt(v) for v in row[1:]])
    with open("self_recovery_hard_targets.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Target", "FoldDisco", "Foldseek", "BLAST"])
        w.writerows(hard_rows)
    print()

    # ---- recall table ----
    print("=" * 78)
    print(f"DOCUMENTED OFF-TARGET RECALL  (n = {n_total}; FoldDisco query = {CONF})")
    print("=" * 78)
    header = ["Method"] + [f"Recall @{c}" for c in CUTOFFS] + ["Median rank"]
    print(f"{header[0]:<12}" + "".join(f"{h:>13}" for h in header[1:]))
    rows = []
    for name, ranks in methods.items():
        got = [r for r in ranks.values() if r is not None]
        counts = [sum(1 for r in got if r <= c) for c in CUTOFFS]
        within = [r for r in got if r <= max(CUTOFFS)]
        med = statistics.median(within) if within else float("nan")
        # %g keeps a genuine half-value (10.5) instead of rounding it away
        print(f"{name:<12}" + "".join(f"{c:>13}" for c in counts) + f"{med:>13g}")
        rows.append([name] + counts + [f"{med:g}"])
    with open("baseline_comparison.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t"); w.writerow(header); w.writerows(rows)

    # ---- split by structural relationship (three tiers) ----
    print()
    print("=" * 78)
    print(f"RECALL @{max(CUTOFFS)} SPLIT BY STRUCTURAL RELATIONSHIP TO THE QUERY POCKET")
    print("=" * 78)
    print(f"{'Relationship':<34}{'n':>3}" + "".join(f"{m:>12}" for m in methods))
    srows = []
    for tier in ("pocket_similar", "fold_pocket_diff", "fold_unrelated"):
        members = [(p, t_) for p, d in KNOWN.items() for t_, lab in d.items() if lab == tier]
        row = [TIER_LABEL[tier], len(members)]
        for name, ranks in methods.items():
            row.append(sum(1 for k in members
                           if ranks[k] is not None and ranks[k] <= max(CUTOFFS)))
        print(f"{row[0]:<34}{row[1]:>3}" + "".join(f"{v:>12}" for v in row[2:]))
        srows.append(row)
    with open("baseline_split.tsv", "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["Relationship", "n"] + list(methods))
        w.writerows(srows)

    # ---- median rank within the pocket-similar tier ----
    print()
    print("Median rank within the pocket-similar tier (recovered targets only):")
    t1 = [(p, t_) for p, d in KNOWN.items() for t_, lab in d.items() if lab == "pocket_similar"]
    for name, ranks in methods.items():
        got = [ranks[k] for k in t1 if ranks[k] is not None and ranks[k] <= max(CUTOFFS)]
        med = statistics.median(got) if got else float("nan")
        print(f"   {name:<12} median {med:>5g}   (recovered {len(got)}/{len(t1)})")

    # ---- targets recovered by one method alone ----
    print()
    print("Documented off-targets recovered by exactly one method (within the cutoff):")
    found = False
    for k in sorted({k for d in (methods.values()) for k in d}):
        inside = {m: (r is not None and r <= max(CUTOFFS)) for m, r in
                  ((m, methods[m][k]) for m in methods)}
        if sum(inside.values()) == 1:
            only = [m for m, v in inside.items() if v][0]
            print(f"   {k[0]}/{k[1]}: only {only} (rank {methods[only][k]})")
            found = True
    if not found:
        print("   none")

    print()
    print("Wrote self_recovery_by_method.tsv, self_recovery_hard_targets.tsv,")
    print("      baseline_comparison.tsv and baseline_split.tsv")


if __name__ == "__main__":
    main()
