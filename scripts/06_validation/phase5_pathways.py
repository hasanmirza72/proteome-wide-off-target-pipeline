#!/usr/bin/env python3
"""
phase5_pathways.py  (roadmap Part 4.x - functional enrichment bridge)
--------------------------------------------------------------------
Pipes each drug's high-confidence off-target GENES into functional enrichment
(GO Biological Process / KEGG / Reactome) via Enrichr, to test the
structural-motif -> pathway -> phenotype bridge.

DESIGN FOR HONESTY (these are the whole point):
  * Enrichment is run PER DRUG (off-targets of different drugs are not pooled),
    because a pathway that "explains a side effect" must belong to one drug.
  * For every drug we run TWO sets: all high hits, and FAMILY-REMOVED (cross-family
    only). The paralog trap means all-hits enrichment is usually just the query's
    own family function (expected, not novel); only family-removed enrichment can
    support a NOVEL-pathway claim - and those sets are usually too small (flagged).
  * Terms matching the query's own family are labelled EXPECTED (positive control),
    not novel discovery.
  * Background: uses a custom background of the searchable library IF
    background_genes.txt (gene symbols, one per line) is provided; otherwise falls
    back to Enrichr's genome background, and SAYS SO - genome background inflates
    significance, so without the custom background treat results as descriptive /
    hypothesis-generating only.

Network: Enrichr / speedrichr API with on-disk caching. Offline -> writes nothing
but the plan and exits cleanly.

Reads : offtarget_hits_bioannotated.tsv (preferred) or offtarget_hits_high.tsv
Writes: phase5_pathways.tsv  (+ printed per-drug top terms)
"""
import csv, os, json, time, sys

BIO   = "offtarget_hits_bioannotated.tsv"
HIGH  = "offtarget_hits_high.tsv"
BGF   = "background_genes.txt"          # optional custom background (gene symbols)
OUT   = "phase5_pathways.tsv"
CACHE = "enrichr_cache"; os.makedirs(CACHE, exist_ok=True)
LIBRARIES = ["GO_Biological_Process_2023","KEGG_2021_Human","Reactome_2022"]
MIN_GENES = 4          # do not enrich sets smaller than this
TOP_TERMS = 5
PADJ = 0.05

# gene fallback map (only used if reading the raw high file without hit_gene col)
GENE = {"D6RBQ6":"USP17L17","D6RCP7":"USP17L19","O14733":"MAP2K7","O14965":"AURKA","O43283":"MAP3K13",
"O60218":"AKR1B10","O60658":"PDE8A","O94768":"STK17B","O94782":"USP1","O94966":"USP19","O95711":"LY86",
"O96017":"CHEK2","P07947":"YES1","P20718":"GZMH","P21802":"FGFR2","P31749":"AKT1","P36507":"MAP2K2",
"P42685":"FRK","P45974":"USP5","P52895":"AKR1C2","P68106":"FKBP1B","Q02338":"BDH1","Q02750":"MAP2K1",
"Q02779":"MAP3K10","Q07343":"PDE4B","Q08493":"PDE4C","Q08499":"PDE4D","Q13163":"MAP2K5","Q13370":"PDE3B",
"Q13882":"PTK6","Q14432":"PDE3A","Q14694":"USP10","Q16512":"PKN1","Q5T2L2":"AKR1C8","Q5TGU0":"TSPO2",
"Q5VVH2":"FKBP1C","Q86XF0":"DHFR2","Q8IY84":"NIM1K","Q8NFD2":"ANKK1","Q8NGB9":"OR4F6","Q96KB5":"PBK",
"Q96T53":"MBOAT4","Q9BQ65":"USB1","Q9H3Y6":"SRMS","Q9H7Z6":"KAT8","Q9HA47":"UCK1","Q9UPT9":"USP22","Q9Y2U5":"MAP3K2"}

# query family -> keywords that mark an enriched term as EXPECTED (not novel)
EXPECTED_KW = {
 "kinase":["kinase","phosphoryl","atp binding","mapk","erk"],
 "phosphodiesterase":["camp","cgmp","phosphodiesterase","nucleotide","cyclic"],
 "DUB":["ubiquitin","deubiquitin","proteasom"],
 "reductase":["reductase","aldo-keto","aldo keto","oxidoreduct","folate","nadp"],
 "isomerase-FKBP":["isomerase","prolyl","immunosupp","fk506","calcineurin"],
}
def qfamily(g):
    g=(g or "").upper()
    if g.startswith("PDE") or g=="USB1": return "phosphodiesterase"
    if g.startswith("USP"): return "DUB"
    if g.startswith("AKR") or g.startswith("DHFR"): return "reductase"
    if g.startswith("FKBP"): return "isomerase-FKBP"
    return "kinase"  # default for this benchmark's remaining queries

def _cache(k):
    return os.path.join(CACHE, k.replace("/","_").replace("?","_").replace("&","_").replace("=","_")+".json")

def enrichr(genes, library, background=None):
    """Return list of (term, adj_p, overlap_genes) or None if offline/failed."""
    ck=_cache(f"{library}__{background is not None}__{'_'.join(sorted(genes))[:120]}")
    if os.path.exists(ck): return json.load(open(ck))
    try:
        import requests
        if background:  # speedrichr custom background
            base="https://maayanlab.cloud/speedrichr/api"
            r=requests.post(f"{base}/addList",
                files={'list':(None,"\n".join(genes)),'description':(None,'hits')},timeout=60)
            uid=r.json()["userListId"]
            rb=requests.post(f"{base}/addbackground",data={'background':"\n".join(background)},timeout=120)
            bid=rb.json()["backgroundid"]
            re_=requests.post(f"{base}/backgroundenrich",
                data={'userListId':uid,'backgroundid':bid,'backgroundType':library},timeout=120)
            data=re_.json().get(library,[])
        else:           # classic Enrichr, genome background
            base="https://maayanlab.cloud/Enrichr"
            r=requests.post(f"{base}/addList",
                files={'list':(None,"\n".join(genes)),'description':(None,'hits')},timeout=60)
            uid=r.json()["userListId"]; time.sleep(0.4)
            re_=requests.get(f"{base}/enrich?userListId={uid}&backgroundType={library}",timeout=60)
            data=re_.json().get(library,[])
        # Enrichr row: [rank, term, pval, zscore, combined, genes, adj_pval, ...]
        out=[]
        for row in data:
            term=row[1]; padj=row[6]; ov=row[5]
            out.append([term, float(padj), ov])
        out.sort(key=lambda x:x[1])
        json.dump(out, open(ck,"w")); return out
    except Exception as e:
        return None

def load_hits():
    """drug_key -> dict(query_gene, drug, genes=set, cross_genes=set)"""
    groups={}
    if os.path.exists(BIO):
        rd=csv.DictReader(open(BIO),delimiter="\t")
        for r in rd:
            k=r["query_target"]; g=groups.setdefault(k,{"query_gene":r["query_gene"],
                "drug":r["drug_name"],"genes":set(),"cross":set()})
            hg=r["hit_gene"].strip(); g["genes"].add(hg)
            if r["family_class"]=="cross-family": g["cross"].add(hg)
    elif os.path.exists(HIGH):
        rd=csv.DictReader(open(HIGH),delimiter="\t")
        for r in rd:
            k=r["query_target"]; hg=GENE.get(r["hit_uniprot"].strip(),r["hit_uniprot"].strip())
            g=groups.setdefault(k,{"query_gene":k,"drug":k,"genes":set(),"cross":set()})
            g["genes"].add(hg)
    else:
        sys.exit(f"missing {BIO} / {HIGH}")
    return groups

def expected(term, qfam):
    kws=EXPECTED_KW.get(qfam,[]); t=term.lower()
    return any(k in t for k in kws)

def main():
    groups=load_hits()
    background=None
    if os.path.exists(BGF):
        background=[l.strip() for l in open(BGF) if l.strip()]
        print(f"[bg] custom background: {len(background)} genes from {BGF}")
    else:
        print("[bg] NO custom background file -> using Enrichr genome background.")
        print("     Results are DESCRIPTIVE only; genome background inflates significance.")
        print("     For a rigorous test, provide background_genes.txt (library gene symbols).\n")

    rows=["\t".join(["drug_target","drug","set_type","n_genes","library","term","adj_p","overlap","flag"])]
    offline=False
    for k,g in sorted(groups.items()):
        qfam=qfamily(g["query_gene"])
        for set_type, genes in (("all_hits",sorted(g["genes"])),
                                ("family_removed",sorted(g["cross"]))):
            if len(genes)<MIN_GENES:
                if set_type=="all_hits":
                    print(f"--- {k} ({g['drug']}, query {g['query_gene']}): {len(genes)} genes < {MIN_GENES}, skipped")
                elif g["cross"]:
                    print(f"      family_removed set too small ({len(genes)}) - cannot test novelty")
                continue
            print(f"--- {k} ({g['drug']}, query {g['query_gene']}) [{set_type}] {len(genes)} genes: {','.join(genes)}")
            any_hit=False
            for lib in LIBRARIES:
                res=enrichr(genes, lib, background)
                if res is None: offline=True; continue
                sig=[t for t in res if t[1]<PADJ][:TOP_TERMS]
                for term,padj,ov in sig:
                    fl="EXPECTED(family)" if expected(term,qfam) else "other"
                    ovs=";".join(ov) if isinstance(ov,list) else str(ov)
                    rows.append("\t".join([k,g["drug"],set_type,str(len(genes)),lib,
                                           term,f"{padj:.2e}",ovs,fl]))
                    print(f"      [{lib.split('_')[0]}] {term[:55]:55s} padj={padj:.1e} {fl}")
                    any_hit=True
            if not any_hit and not offline:
                print("      (no terms below padj<0.05)")
    open(OUT,"w").write("\n".join(rows)+"\n")
    if offline:
        print("\n[!] Enrichr not reachable - run on a networked node (cache makes reruns free).")
    print(f"\nWrote {OUT}")
    print("\nHONEST READING:")
    print(" - 'EXPECTED(family)' terms are a positive control: the hits ARE that family, so")
    print("   enriching for its function is trivially true, not a discovered off-target pathway.")
    print(" - A NOVEL bridge requires an enriched term in a FAMILY_REMOVED set. If those sets")
    print("   are too small to test (common here), state that the data cannot support a novel")
    print("   pathway claim rather than citing the all-hits (paralog-driven) enrichment.")

if __name__ == "__main__":
    main()
