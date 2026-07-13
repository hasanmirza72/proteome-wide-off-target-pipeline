#!/usr/bin/env python3
"""
SIDER + openFDA PDB miner  (corrected build)

Fixes vs. the three previous versions:
  * Script 2 crash  -> dict.get(k, default) does NOT return the default when the key
                       exists with value None. rcsb_entity_source_organism is
                       sometimes null, so [{}][0] became None[0]. Now handled safely.
  * Scripts 1 & 2   -> the ligand->CID bridge read a field that does not exist in the
                       RCSB chemcomp document, so every CID came back None. Because the
                       openFDA fallback is gated behind `if pubchem_cid:`, BOTH validators
                       were silently skipped -> 0 hits everywhere. Fixed with a schema
                       tolerant PubChem extractor (works no matter what the exact key is).
  * Script 3 error  -> nonpolymer_comp.chem_comp.pdbx_chem_comp_external_references is not
                       a valid GraphQL field. We stop asking GraphQL for CIDs and use the
                       REST chemcomp endpoint (parallelised) instead.
"""

import urllib.request
import urllib.error
import urllib.parse
import json
import time
import os
import concurrent.futures
import pandas as pd

# ---------------- CONFIGURATION ----------------
VERBOSE = True            # per-ligand logging so you can watch it work
DEBUG_FIRST_LIGAND = True # print the raw chemcomp JSON keys of the very first ligand once

TARGET_CLASSES = {
    "Kinase": "kinase",
    "Nuclear Receptor": "nuclear receptor",
    "Phosphodiesterase": "phosphodiesterase",
    "Protease": "protease",
    "Carbonic Anhydrase": "carbonic anhydrase",
    "Reductase": "reductase",
    "Transferase": "transferase",
    "Isomerase": "isomerase",
    "Hydrolase": "hydrolase",
}

TARGET_PER_CLASS = 10
DIR_MONOMERS = "dataset_monomers_pdb"
DIR_DIMERS = "dataset_dimers_pdb"

DRUG_NAMES_FILE = "drug_names.tsv"
SE_FILE = "meddra_all_se.tsv.gz"
DRUG_NAMES_URL = "http://sideeffects.embl.de/media/download/drug_names.tsv"
SE_URL = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz"

IGNORE_LIST = {
    "HOH", "SO4", "DMS", "EDO", "GOL", "CL", "NA", "MG", "ZN", "PO4",
    "IOD", "BR", "CA", "K", "F", "NI", "CU", "CO", "MN", "CD", "HG",
    "MO", "W", "AU", "PT", "YB",
    "PEG", "FMT", "ACY", "UNX", "UNL", "TRS", "MES", "HEPES", "EPE",
    "MPD", "DTT", "BME", "CIT", "FLC", "TLA", "BOG", "MPO", "MOP",
    "PIP", "POP", "CAC", "PGE", "PG4", "1PE", "NHE",
    "ATP", "ADP", "AMP", "GTP", "GDP", "GMP", "NAD", "NAP", "NDP", "NADP",
    "FAD", "FMN", "HEM", "HEC", "SAH", "SAM", "PLP", "COA", "ACO",
    "MYR", "PLM", "OLA", "STE", "OLC", "OLB", "CLR",
    "GLC", "FRU", "SUC", "TRE", "MAN", "GAL", "NAG", "NDG",
    "BGC", "BMA", "SIA", "UDP", "GSH", "GSN", "BHE", "ACT", "IMD",
}

# ---------------- ERROR LOGGING & CACHES ----------------
seen_errors = set()
def log_error_once(msg):
    if msg not in seen_errors:
        print(f"  [!] ERROR/WARNING: {msg}")
        seen_errors.add(msg)

cid_cache = {}     # comp_id -> int CID or None
title_cache = {}   # cid -> generic name or None
_debug_state = {"dumped": False}


def http_json(url, timeout=8, data=None):
    req = urllib.request.Request(
        url,
        data=data,
        headers={
            "User-Agent": "Bioinformatics Pipeline",
            "Content-Type": "application/json",
        },
    )
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return json.loads(r.read().decode())


# ---------------- 1. PubChem CID bridge (the key fix) ----------------
def _extract_pubchem_cid(obj):
    """
    Schema-tolerant scan of an RCSB chemcomp JSON document. It looks for ANY nested
    dict whose 'name' field literally says 'PubChem' and grabs the sibling id value.
    This works whether the schema calls it:
        rcsb_chem_comp_related / resource_name / resource_accession_code   (current)
        pdbx_chem_comp_external_references / reference_resource_name / resource_id
    or anything similar, so you never have to guess the exact key again.
    """
    name_keys = ("reference_resource_name", "resource_name", "provenance_source",
                 "database_name", "reference_database_name")
    id_keys = ("resource_id", "resource_accession_code", "reference_resource_id",
               "accession_code", "identifier", "reference_id")

    def walk(node):
        if isinstance(node, dict):
            label = None
            for k in name_keys:
                v = node.get(k)
                if isinstance(v, str):
                    label = v
                    break
            if label and label.strip().lower() == "pubchem":
                for k in id_keys:
                    v = node.get(k)
                    if v is not None:
                        digits = "".join(ch for ch in str(v) if ch.isdigit())
                        if digits:
                            return int(digits)
            for v in node.values():
                got = walk(v)
                if got is not None:
                    return got
        elif isinstance(node, list):
            for v in node:
                got = walk(v)
                if got is not None:
                    return got
        return None

    return walk(obj)


def get_pubchem_cid_from_pdb(comp_id):
    """Resolve a 3-letter PDB ligand code -> PubChem CID via the RCSB chemcomp REST doc."""
    if comp_id in cid_cache:
        return cid_cache[comp_id]
    try:
        data = http_json(f"https://data.rcsb.org/rest/v1/core/chemcomp/{comp_id}", timeout=8)

        if DEBUG_FIRST_LIGAND and not _debug_state["dumped"]:
            _debug_state["dumped"] = True
            print(f"  [debug] chemcomp top-level keys for {comp_id}: {sorted(data.keys())}")
            cid_preview = _extract_pubchem_cid(data)
            print(f"  [debug] extracted PubChem CID for {comp_id}: {cid_preview}")

        cid = _extract_pubchem_cid(data)
        cid_cache[comp_id] = cid
        return cid
    except Exception:
        cid_cache[comp_id] = None
        return None


def get_generic_name_from_pubchem(cid):
    if cid in title_cache:
        return title_cache[cid]
    try:
        url = (f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}"
               f"/property/Title/JSON")
        data = http_json(url, timeout=8)
        title = data["PropertyTable"]["Properties"][0]["Title"]
        title_cache[cid] = title
        time.sleep(0.15)  # PubChem politeness (max ~5 req/s)
        return title
    except Exception:
        title_cache[cid] = None
        return None


def check_openfda(drug_name, retries=1):
    if not drug_name or drug_name == "Unknown":
        return False
    clean = urllib.parse.quote(drug_name.lower().split()[0])
    url = (f"https://api.fda.gov/drug/label.json"
           f'?search=openfda.generic_name:"{clean}"&limit=1')
    try:
        data = http_json(url, timeout=8)
        results = data.get("results", [])
        return bool(results and "adverse_reactions" in results[0])
    except urllib.error.HTTPError as e:
        if e.code == 429 and retries > 0:
            print("      [!] openFDA rate limit — sleeping 5s...")
            time.sleep(5)
            return check_openfda(drug_name, retries=0)
        if e.code != 404:
            log_error_once(f"openFDA HTTP {e.code} for {drug_name}")
    except Exception:
        pass
    return False


# ---------------- 2. SIDER SETUP ----------------
def download_file(url, dest):
    if os.path.exists(dest):
        return True
    print(f"Downloading {dest}...")
    try:
        urllib.request.urlretrieve(url, dest)
        return True
    except Exception as e:
        log_error_once(f"Failed to download {dest}: {e}")
        return False


print("=" * 90)
print("   SIDER + openFDA PDB MINER (corrected)")
print("=" * 90)

sider_available = download_file(DRUG_NAMES_URL, DRUG_NAMES_FILE) and download_file(SE_URL, SE_FILE)
stitch_to_name = {}
se_counts_by_stitch = {}

if sider_available:
    print("Loading SIDER flat files into memory...")
    try:
        sider_names = pd.read_csv(DRUG_NAMES_FILE, sep="\t", header=None,
                                  names=["stitch_id", "drug_name"], on_bad_lines="skip")
        stitch_to_name = dict(zip(sider_names.stitch_id, sider_names.drug_name))

        se_cols = ["stitch_id_flat", "stitch_id_stereo", "umls_label",
                   "meddra_type", "umls_meddra", "side_effect_name"]
        sider_se = pd.read_csv(SE_FILE, sep="\t", compression="gzip", header=None,
                               names=se_cols, on_bad_lines="skip")
        sider_se_pt = sider_se[sider_se["meddra_type"] == "PT"]
        se_counts_by_stitch = (sider_se_pt.groupby("stitch_id_flat")["side_effect_name"]
                               .nunique().to_dict())
        print(f"Loaded {len(stitch_to_name)} SIDER drug names, "
              f"{len(se_counts_by_stitch)} with side-effect data.\n")
    except Exception as e:
        log_error_once(f"SIDER parsing failed: {e}")
        sider_available = False

if not sider_available:
    print("SIDER unavailable — relying on PubChem + openFDA only.\n")


# ---------------- 3. PDB SEARCH ----------------
def get_human_experimental_pdbs(class_keyword):
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entity_source_organism.scientific_name",
                    "operator": "exact_match", "value": "Homo sapiens"}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "exptl.method", "operator": "in",
                    "value": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY", "SOLUTION NMR"]}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_polymer_entity.pdbx_description",
                    "operator": "contains_words", "value": class_keyword}},
                {"type": "terminal", "service": "text", "parameters": {
                    "attribute": "rcsb_entry_info.nonpolymer_entity_count",
                    "operator": "greater", "value": 0}},
            ],
        },
        "request_options": {"return_all_hits": True},
        "return_type": "entry",
    }
    try:
        body = http_json("https://search.rcsb.org/rcsbsearch/v2/query", timeout=20,
                         data=json.dumps(query).encode())
        return [hit["identifier"] for hit in body.get("result_set", [])]
    except Exception as e:
        log_error_once(f"Search API failed for {class_keyword}: {e}")
        return []


# ---------------- 4. GRAPHQL METADATA (entry-level only, which is valid) ----------------
def fetch_metadata_batch(pdb_ids):
    query = """
    query($ids: [String!]!) {
      entries(entry_ids: $ids) {
        rcsb_id
        assemblies { pdbx_struct_assembly { oligomeric_count } }
        polymer_entities {
          rcsb_entity_source_organism { scientific_name }
          rcsb_polymer_entity { pdbx_description }
        }
        nonpolymer_entities {
          pdbx_entity_nonpoly { comp_id name }
        }
      }
    }
    """
    try:
        payload = json.dumps({"query": query, "variables": {"ids": pdb_ids}}).encode()
        data = http_json("https://data.rcsb.org/graphql", timeout=20, data=payload)
        if "errors" in data:
            log_error_once(f"GraphQL Error: {data['errors'][0]['message']}")
            return []
        return data.get("data", {}).get("entries", []) or []
    except Exception as e:
        log_error_once(f"Batch GraphQL Error: {e}")
        return []


# ---------------- helpers to survive null fields ----------------
def is_human(polys):
    for p in polys:
        for org in (p.get("rcsb_entity_source_organism") or []):
            if "homo sapiens" in str((org or {}).get("scientific_name", "")).lower():
                return True
    return False


def assembly_type_of(entry):
    for a in (entry.get("assemblies") or []):
        pa = a.get("pdbx_struct_assembly") or {}
        try:
            count = int(pa.get("oligomeric_count"))
        except (TypeError, ValueError):
            continue
        if count == 1:
            return "MONOMER"
        if count == 2:
            return "DIMER"
    return None


# ---------------- 5. MAIN ----------------
if __name__ == "__main__":
    os.makedirs(DIR_MONOMERS, exist_ok=True)
    os.makedirs(DIR_DIMERS, exist_ok=True)

    master_monomers, master_dimers = {}, {}
    seen_monomer_drugs, seen_dimer_drugs = set(), set()

    for p_class, keyword in TARGET_CLASSES.items():
        print(f"\nMining Class: {p_class}")
        master_monomers[p_class], master_dimers[p_class] = [], []

        hits = get_human_experimental_pdbs(keyword)
        if not hits:
            print(f"  No raw structures found for {p_class}.")
            continue
        print(f"  Discovered {len(hits)} raw human structures. Batch-evaluating...")

        batch_size = 50
        for i in range(0, len(hits), batch_size):
            if (len(master_monomers[p_class]) >= TARGET_PER_CLASS
                    and len(master_dimers[p_class]) >= TARGET_PER_CLASS):
                break

            entries = fetch_metadata_batch(hits[i:i + batch_size])
            time.sleep(0.2)

            # ---- parallel CID pre-fetch for every unknown ligand in this batch ----
            unknown = set()
            for entry in entries:
                if not entry:
                    continue
                for nonpoly in (entry.get("nonpolymer_entities") or []):
                    cid_code = (nonpoly.get("pdbx_entity_nonpoly") or {}).get("comp_id", "")
                    if cid_code and cid_code not in IGNORE_LIST and cid_code not in cid_cache:
                        unknown.add(cid_code)
            if unknown:
                with concurrent.futures.ThreadPoolExecutor(max_workers=8) as ex:
                    list(ex.map(get_pubchem_cid_from_pdb, unknown))

            # ---- evaluate each structure ----
            for entry in entries:
                if not entry:
                    continue
                pdb_id = entry.get("rcsb_id")
                polys = entry.get("polymer_entities") or []
                if not is_human(polys):
                    continue

                atype = assembly_type_of(entry)
                if not atype:
                    continue
                if atype == "MONOMER" and len(master_monomers[p_class]) >= TARGET_PER_CLASS:
                    continue
                if atype == "DIMER" and len(master_dimers[p_class]) >= TARGET_PER_CLASS:
                    continue

                protein_desc = ((polys[0].get("rcsb_polymer_entity") or {})
                                .get("pdbx_description", "Unknown")) if polys else "Unknown"

                matched_drug, matched_source, matched_se = None, "Unknown", 0

                for nonpoly in (entry.get("nonpolymer_entities") or []):
                    nonpoly_info = nonpoly.get("pdbx_entity_nonpoly") or {}
                    comp_id = nonpoly_info.get("comp_id", "")
                    raw_name = nonpoly_info.get("name", "Unknown")
                    if not comp_id or comp_id in IGNORE_LIST:
                        continue

                    cid = get_pubchem_cid_from_pdb(comp_id)  # O(1) from cache now
                    if not cid:
                        continue

                    # Pass 1: SIDER
                    if sider_available:
                        flat = f"CID1{cid:08d}"
                        if flat in se_counts_by_stitch:
                            generic = stitch_to_name.get(flat) or \
                                get_generic_name_from_pubchem(cid) or raw_name
                            matched_drug = comp_id
                            matched_source = f"SIDER ({str(generic).title()})"
                            matched_se = se_counts_by_stitch.get(flat, 0)
                            break

                    # Pass 2: openFDA
                    generic = get_generic_name_from_pubchem(cid)
                    if generic and check_openfda(generic):
                        matched_drug = comp_id
                        matched_source = f"openFDA ({generic})"
                        break

                    if VERBOSE:
                        print(f"    - {pdb_id}:{comp_id} (CID {cid}) not clinically validated")

                if not matched_drug:
                    continue

                if atype == "MONOMER" and matched_drug in seen_monomer_drugs:
                    continue
                if atype == "DIMER" and matched_drug in seen_dimer_drugs:
                    continue

                obj = {"pdb_id": pdb_id, "protein": protein_desc,
                       "drug_id": matched_drug, "source": matched_source,
                       "se_count": matched_se}

                if atype == "MONOMER":
                    master_monomers[p_class].append(obj)
                    seen_monomer_drugs.add(matched_drug)
                    dl_path = os.path.join(DIR_MONOMERS, f"{pdb_id}.pdb")
                    print(f"  [MONOMER] {pdb_id} | Drug: {matched_drug} | {matched_source}")
                else:
                    master_dimers[p_class].append(obj)
                    seen_dimer_drugs.add(matched_drug)
                    dl_path = os.path.join(DIR_DIMERS, f"{pdb_id}.pdb")
                    print(f"  [DIMER]   {pdb_id} | Drug: {matched_drug} | {matched_source}")

                try:
                    if not os.path.exists(dl_path):
                        urllib.request.urlretrieve(
                            f"https://files.rcsb.org/download/{pdb_id}.pdb", dl_path)
                except Exception:
                    pass

        print("-" * 80)
        print(f"  Summary {p_class}: "
              f"{len(master_monomers[p_class])} Monomers | {len(master_dimers[p_class])} Dimers")
        print("-" * 80)

    print("\n" + "=" * 90)
    print("PIPELINE COMPLETE")
    print("=" * 90)

    with open("dataset_monomers.json", "w") as f:
        json.dump(master_monomers, f, indent=4)
    with open("dataset_dimers.json", "w") as f:
        json.dump(master_dimers, f, indent=4)

    total_m = sum(len(v) for v in master_monomers.values())
    total_d = sum(len(v) for v in master_dimers.values())
    print(f"Monomers: dataset_monomers.json ({total_m}) -> ./{DIR_MONOMERS}/")
    print(f"Dimers:   dataset_dimers.json ({total_d}) -> ./{DIR_DIMERS}/")
