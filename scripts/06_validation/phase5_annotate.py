#!/usr/bin/env python3
"""
phase5_annotate.py  (roadmap Part 4.1 — shared annotation layer for Phase 5)
---------------------------------------------------------------------------
Builds the enriched table every downstream Phase 5 probe reads. For each
high-confidence off-target hit it attaches:
  - query drug identity (+ SIDER side-effect count) and query gene
  - hit gene name and within/cross-family label
  - ChEMBL: is the hit a known drug target? does the query drug have a measured
    binding value (Ki/Kd/IC50/EC50) against it? best potency (nM)?
  - evidence tier (roadmap 4.5): 1 = structural + ChEMBL binding;
    2 = structural + known drug target; 3 = structural only

Network: hits the ChEMBL web API with on-disk caching (reruns are free/reproducible).
If offline, ChEMBL fields are marked 'no_network' and everything else still works.
ABSENCE OF A CHEMBL RECORD IS NOT EVIDENCE AGAINST A HIT (coverage gap) — do not
read a missing binding value as "does not bind".

Output: offtarget_hits_bioannotated.tsv
"""
import csv, os, json, time, hashlib, sys

HIGH   = "offtarget_hits_high.tsv"
OUT    = "offtarget_hits_bioannotated.tsv"
CACHE  = "chembl_cache"
BASE   = "https://www.ebi.ac.uk/chembl/api/data"
ACT_TYPES = {"Ki","Kd","IC50","EC50"}
# Potency thresholds (nM). Many panels report non-binders at a fixed ceiling
# (commonly 10,000 nM = 10 uM), so a value alone does NOT mean binding.
BIND_NM   = 1000.0     # <= this = genuine binding (confirmed)
CEILING_NM = 10000.0   # >= this = assay ceiling => "tested-inactive", NOT binding
os.makedirs(CACHE, exist_ok=True)

# --- query drug + gene per target (from dataset_monomers.json + verified names) ---
QUERY = {  # target: (gene, uniprot, drug_name, se_count, source)
 "6VCJ":("DHFR","P00374","Naproxen",0,"openFDA"),
 "1UEI":("UCK2","Q9BZX2","Uridine triphosphate",0,"openFDA"),
 "1NHZ":("NR3C1","P04150","Mifepristone",0,"openFDA"),
 "6L6E":("PDE5A","O76074","Avanafil",108,"SIDER"),
 "4U7Z":("MAP2K1","Q02750","Selumetinib",0,"openFDA"),
 "1XBB":("SYK","P43405","Imatinib",323,"SIDER"),
 "3TVX":("PDE4A","P27815","Pentoxifylline",94,"SIDER"),
 "3RX3":("AKR1B1","P15121","Sulindac",0,"openFDA"),
 "5V1M":("USB1","Q9BQ65","Uridine monophosphate",0,"openFDA"),
 "6Q6O":("RORC","P51449","Cholic acid",0,"openFDA"),
 "6GH9":("USP15","Q9Y4E8","Mitoxantrone",158,"SIDER"),
 "8HJE":("USP28","Q96RU2","Vismodegib",50,"SIDER"),
 "2XP2":("ALK","Q9UM73","Crizotinib",0,"openFDA"),
 "2Y05":("PTGR1","Q14914","Raloxifene",88,"SIDER"),
 "7AXA":("NR1I2","O75469","Clotrimazole",26,"SIDER"),
 "1AVN":("CA2","P00918","Histamine",96,"SIDER"),
 "7XPY":("USP7","Q93009","Eupalinolide B",0,"openFDA"),
 "6YG2":("MAP2K7","O14733","Ibrutinib",0,"openFDA"),
 "4ODR":("FKBP1A","P62942","Tacrolimus",0,"openFDA"),
 "3TC5":("PIN1","Q13526","Dexamethasone",0,"openFDA"),
 "1GS4":("AR","P10275","Fludrocortisone",76,"SIDER"),
 "4R7L":("LTA4H","P09960","Vorinostat",62,"SIDER"),
 "9XRI":("NUDT1","P36639","Acoramidis",0,"openFDA"),
}

# --- hit accession -> gene (UniProt-verified) ---
GENE = {
"D6RBQ6":"USP17L17","D6RCP7":"USP17L19","O14733":"MAP2K7","O14965":"AURKA","O43283":"MAP3K13",
"O60218":"AKR1B10","O60658":"PDE8A","O94768":"STK17B","O94782":"USP1","O94966":"USP19","O95711":"LY86",
"O96017":"CHEK2","P07947":"YES1","P20718":"GZMH","P21802":"FGFR2","P31749":"AKT1","P36507":"MAP2K2",
"P42685":"FRK","P45974":"USP5","P52895":"AKR1C2","P68106":"FKBP1B","Q02338":"BDH1","Q02750":"MAP2K1",
"Q02779":"MAP3K10","Q07343":"PDE4B","Q08493":"PDE4C","Q08499":"PDE4D","Q13163":"MAP2K5","Q13370":"PDE3B",
"Q13882":"PTK6","Q14432":"PDE3A","Q14694":"USP10","Q16512":"PKN1","Q5T2L2":"AKR1C8","Q5TGU0":"TSPO2",
"Q5VVH2":"FKBP1C","Q86XF0":"DHFR2","Q8IY84":"NIM1K","Q8NFD2":"ANKK1","Q8NGB9":"OR4F6","Q96KB5":"PBK",
"Q96T53":"MBOAT4","Q9BQ65":"USB1","Q9H3Y6":"SRMS","Q9H7Z6":"KAT8","Q9HA47":"UCK1","Q9UPT9":"USP22",
"Q9Y2U5":"MAP3K2",
}

KINASE_GENES = {"SYK","ALK","MAP2K1","MAP2K2","MAP2K5","MAP2K7","MAP3K2","MAP3K10","MAP3K13","AURKA",
"STK17B","CHEK2","YES1","FGFR2","AKT1","FRK","PTK6","PKN1","NIM1K","ANKK1","PBK","SRMS","UCK1","UCK2"}
def family(gene):
    g=(gene or "").upper()
    if g in KINASE_GENES: return "kinase"
    if g.startswith(("USP","DUB")): return "DUB"
    if g.startswith("PDE") or g=="USB1": return "phosphodiesterase"
    if g.startswith("AKR"): return "reductase"
    if g.startswith("DHFR"): return "reductase-DHFR"
    if g.startswith("FKBP"): return "isomerase-FKBP"
    if g.startswith("ROR") or g.startswith("NR"): return "nuclear-receptor"
    return "other:"+g

def fam_class(qgene, hgene):
    qf, hf = family(qgene), family(hgene)
    # kinase superfamily counts as within (shared ATP site); DHFR~DHFR, AKR~AKR etc.
    if qf==hf: return "within-family"
    if "kinase" in qf and "kinase" in hf: return "within-family"
    if qf.startswith("reductase") and hf.startswith("reductase"): return "within-family"
    return "cross-family"

# ---------------- ChEMBL web API (cached) --------------------------------------
def _get(url):
    key = os.path.join(CACHE, hashlib.md5(url.encode()).hexdigest()+".json")
    if os.path.exists(key):
        return json.load(open(key))
    try:
        import requests
        for attempt in range(4):
            r = requests.get(url, headers={"Accept":"application/json"}, timeout=30)
            if r.status_code==200:
                d=r.json(); json.dump(d, open(key,"w")); return d
            if r.status_code in (429,500,502,503): time.sleep(2*(attempt+1)); continue
            return {"_error":r.status_code}
        return {"_error":"retries"}
    except Exception as e:
        return {"_offline":str(e)}

def chembl_target(uniprot, cache_t={}):
    if uniprot in cache_t: return cache_t[uniprot]
    d=_get(f"{BASE}/target?target_components__accession={uniprot}&format=json&limit=1")
    tid=None
    if isinstance(d,dict) and d.get("targets"):
        tid=d["targets"][0].get("target_chembl_id")
    cache_t[uniprot]=(tid, d.get("_offline") or d.get("_error"))
    return cache_t[uniprot]

def chembl_molecule(name, cache_m={}):
    if name in cache_m: return cache_m[name]
    d=_get(f"{BASE}/molecule?pref_name__iexact={name.replace(' ','%20')}&format=json&limit=1")
    mid=None
    if isinstance(d,dict) and d.get("molecules"):
        mid=d["molecules"][0].get("molecule_chembl_id")
    cache_m[name]=(mid, d.get("_offline") or d.get("_error"))
    return cache_m[name]

def is_known_drug_target(tid):
    if not tid: return "unknown"
    d=_get(f"{BASE}/mechanism?target_chembl_id={tid}&format=json&limit=1")
    if isinstance(d,dict) and "mechanisms" in d:
        return "yes" if d["mechanisms"] else "no"
    return "unknown"

def best_activity(mid, tid):
    """Return (binding_status, best_nM, type). binding_status is one of:
       confirmed (<= BIND_NM), tested_inactive (>= CEILING_NM), weak (between),
       no_measured_value, unknown."""
    if not (mid and tid): return ("unknown","","")
    d=_get(f"{BASE}/activity?molecule_chembl_id={mid}&target_chembl_id={tid}&format=json&limit=200")
    if not (isinstance(d,dict) and "activities" in d): return ("unknown","","")
    best=None; btype=""
    for a in d["activities"]:
        t=a.get("standard_type"); v=a.get("standard_value"); u=a.get("standard_units")
        if t in ACT_TYPES and v not in (None,"") and u=="nM":
            try:
                fv=float(v)
                if best is None or fv<best: best=fv; btype=t
            except: pass
    if best is None: return ("no_measured_value","","")
    if best <= BIND_NM:      status="confirmed"
    elif best >= CEILING_NM: status="tested_inactive"   # assay ceiling, not a binder
    else:                    status="weak"              # measurable but 1-10 uM
    return (status, f"{best:.1f}", btype)

# ------------------------------------------------------------------------------
def tier(binds, known):
    # only GENUINE binding is Tier 1; ceiling/weak values are not confirmation
    if binds=="confirmed": return "1_binding_confirmed"
    if known=="yes":       return "2_known_drug_target"
    return "3_structural_only"

def main():
    if not os.path.exists(HIGH): sys.exit(f"missing {HIGH}")
    rows=list(csv.DictReader(open(HIGH), delimiter="\t"))
    cols=["query_target","query_gene","query_uniprot","drug_name","drug_se_count","drug_source",
          "hit_uniprot","hit_gene","family_class","idf","node_count","overlap_frac","confidence",
          "chembl_target_id","is_known_drug_target","binding_status","best_activity_nM",
          "activity_type","evidence_tier"]
    out=["\t".join(cols)]
    net_ok=True
    for r in rows:
        qt=r["query_target"]; q=QUERY.get(qt,("?","?","?",0,"?"))
        hacc=r["hit_uniprot"]; hgene=GENE.get(hacc,hacc)
        fc=fam_class(q[0],hgene)
        tid,terr=chembl_target(hacc)
        mid,merr=chembl_molecule(q[2])
        if terr or merr: net_ok=False
        known=is_known_drug_target(tid) if tid else "unknown"
        binds,act,atype=best_activity(mid,tid)
        out.append("\t".join(str(x) for x in [
            qt,q[0],q[1],q[2],q[3],q[4],hacc,hgene,fc,
            r.get("idf",""),r.get("node_count",""),r.get("overlap_frac",""),r.get("confidence",""),
            tid or "", known, binds, act, atype, tier(binds,known)]))
    open(OUT,"w").write("\n".join(out)+"\n")

    # summary
    nrows=len(rows)
    wf=sum(1 for l in out[1:] if l.split("\t")[8]=="within-family")
    cf=nrows-wf
    print(f"Wrote {OUT}  ({nrows} hits)")
    print(f"Family split: within={wf}  cross={cf}")
    if not net_ok:
        print("\n[!] ChEMBL not reachable (offline or blocked) — ChEMBL columns are 'unknown'/'no_network'.")
        print("    Re-run on a networked node; cached responses make reruns free.")
    else:
        b=[l.split('\t')[15] for l in out[1:]]        # binding_status column
        conf=b.count("confirmed"); inact=b.count("tested_inactive"); weak=b.count("weak")
        noval=b.count("no_measured_value"); unk=b.count("unknown")
        t1=sum(1 for l in out[1:] if l.split('\t')[-1].startswith('1'))
        t2=sum(1 for l in out[1:] if l.split('\t')[-1].startswith('2'))
        print(f"\nBinding (potency-thresholded, BIND<= {BIND_NM:.0f}nM, ceiling>= {CEILING_NM:.0f}nM):")
        print(f"   confirmed(<=1uM) = {conf}")
        print(f"   weak(1-10uM)     = {weak}")
        print(f"   tested-inactive(>=10uM, assay ceiling) = {inact}   <- NOT binding")
        print(f"   no measured value / untested            = {noval+unk}   <- coverage gap, NOT evidence against")
        print(f"Evidence tiers: Tier1(confirmed binding)={t1}  Tier2(known target)={t2}  "
              f"Tier3(structural only)={nrows-t1-t2}")
        print("\nReminder: 'tested-inactive' means measured >= ceiling (a negative); "
              "'untested' means ChEMBL has no value (not a negative).")

if __name__ == "__main__":
    main()
