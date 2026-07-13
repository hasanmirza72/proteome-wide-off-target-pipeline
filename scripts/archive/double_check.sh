#!/bin/bash
# Run this AFTER re-running BLAST and Foldseek on the rebuilt queries.
# It checks the QUERY column (col 1) specifically, and confirms self-hits.

echo "=== 1. Are the QUERIES correct? (P62942 in, Q5SLE7 out) ==="
echo "--- BLAST queries (col 1): ---"
echo -n "  P62942 as a query? "; cut -f1 blast_hits.tsv | grep -qx "P62942" && echo "YES (good)" || echo "NO (still broken)"
echo -n "  Q5SLE7 as a query?  "; cut -f1 blast_hits.tsv | grep -qx "Q5SLE7" && echo "YES (STILL CONTAMINATED - did BLAST re-run?)" || echo "no (good)"
echo "--- Foldseek queries (col 1): ---"
echo -n "  P62942 as a query? "; cut -f1 foldseek_hits.tsv | grep -qx "P62942" && echo "YES (good)" || echo "NO (still broken)"
echo -n "  Q5SLE7 as a query?  "; cut -f1 foldseek_hits.tsv | grep -qx "Q5SLE7" && echo "YES (STILL CONTAMINATED - did Foldseek re-run?)" || echo "no (good)"

echo
echo "=== 2. Distinct query count (must be 23 in each) ==="
echo -n "  BLAST:    "; cut -f1 blast_hits.tsv | sort -u | wc -l
echo -n "  Foldseek: "; cut -f1 foldseek_hits.tsv | sort -u | wc -l

echo
echo "=== 3. Does P62942 now SELF-recover? (it's in the library, so it should) ==="
echo -n "  BLAST P62942 self-hit: "
awk -F'\t' '$1=="P62942"{split($2,a,"|"); if(a[2]=="P62942") print "YES rank-ok ("$5")"}' blast_hits.tsv | head -1
echo -n "  Foldseek P62942 self-hit: "
awk -F'\t' '$1=="P62942" && $2=="P62942"{print "YES ("$3")"}' foldseek_hits.tsv | head -1

echo
echo "=== 4. All 23 queries have a self-hit now? (expect 23/23 each) ==="
echo -n "  BLAST self-hits:    "
awk -F'\t' '{split($2,a,"|"); if($1==a[2]) print $1}' blast_hits.tsv | sort -u | wc -l
echo -n "  Foldseek self-hits: "
awk -F'\t' '$1==$2{print $1}' foldseek_hits.tsv | sort -u | wc -l

echo
echo "If 1-4 all look good (P62942 in with self-hit, Q5SLE7 gone, 23 queries,"
echo "23/23 self-hits), the baseline data is finally clean. Upload both hit files."
