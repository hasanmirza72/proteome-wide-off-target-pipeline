#!/usr/bin/env python3
"""
AHoJ apo-holo pairing — FINAL.

For each holo (your PDB + drug), it:
  1. queries AHoJ (bare PDB id -> auto-detect pockets, the mode that works),
  2. keeps only the apo candidates that are unbound at the DRUG's pocket
     (matched via target_poi = chain_LIGAND_resnum),
  3. deduplicates candidates to one row per apo PDB (best row kept),
  4. ranks them: experimental > full binding site present > most similar to the
     holo (pocket_rmsd) > best resolution,
  5. assigns ONE UNIQUE apo per holo (no apo is reused across holos),
  6. downloads AHoJ's apo structure ALREADY ALIGNED to the holo (ideal for P2Rank,
     since predicted-pocket coordinates are then directly comparable),
  7. writes dataset_apo_pairs.json + apo_holo_pairs.csv.

Ranking weights completeness above similarity on purpose: if the apo is missing
binding-site residues, P2Rank can't fairly find the pocket and the comparison breaks.
Flip PREFER_SIMILARITY_OVER_RESOLUTION if you want geometry to outrank resolution.
"""

import json
import os
import sys
import csv
import time
import urllib.request
import urllib.error

# ---------------- CONFIG ----------------
INPUT_FILE = "dataset_monomers.json"
DIR_APO = "dataset_apo_pdb"
OUT_JSON = "dataset_apo_pairs.json"
OUT_CSV = "apo_holo_pairs.csv"
SITE_ROOT = "https://apoholo.cz"
API_BASE = SITE_ROOT + "/api"

POLL_S = 6
MAX_WAIT_S = 300
PREFER_SIMILARITY_OVER_RESOLUTION = False


# ---------------- HTTP ----------------
def http(url, method="GET", payload=None, timeout=30):
    data = json.dumps(payload).encode() if payload is not None else None
    headers = {"Accept": "application/json"}
    if data is not None:
        headers["Content-Type"] = "application/json"
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            return r.getcode(), r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()
    except Exception:
        return None, b""


def http_json(url, method="GET", payload=None, timeout=30):
    code, body = http(url, method, payload, timeout)
    try:
        return code, json.loads(body.decode(errors="replace"))
    except Exception:
        return code, None


# ---------------- AHoJ ----------------
def submit(pdb_id):
    _, j = http_json(f"{API_BASE}/job", method="POST",
                     payload={"job_name": f"m_{pdb_id}", "queries": pdb_id, "options": {}})
    return j.get("job_id") if isinstance(j, dict) else None


def wait_done(job_id):
    deadline = time.time() + MAX_WAIT_S
    while time.time() < deadline:
        _, job = http_json(f"{API_BASE}/job/{job_id}")
        if isinstance(job, dict):
            q = (job.get("queries") or [None])[0]
            if job.get("done") or (isinstance(q, dict) and q.get("status") == "done"):
                return q
        time.sleep(POLL_S)
    return None


# ---------------- candidate handling ----------------
def num(v, default):
    try:
        return float(v)
    except (TypeError, ValueError):
        return default


def poi_ligand(e):
    """target_poi like 'A_8E8_514' -> '8E8'."""
    parts = (e.get("target_poi") or "").split("_")
    return parts[1].upper() if len(parts) >= 2 else ""


def rank_key(e):
    sim = num(e.get("pocket_rmsd"), 1e9)
    res = num(e.get("resolution"), 1e9)
    tail = (sim, res) if PREFER_SIMILARITY_OVER_RESOLUTION else (res, sim)
    return (
        0 if not e.get("is_alphafold") else 1,           # experimental first
        -num(e.get("mapped_binding_residues_percent"), 0),  # full pocket present
        num(e.get("mapped_binding_residues_unobserved"), 1e9),  # no missing residues
        *tail,
    )


def candidates_for_drug(found_apo, drug_id, holo_pdb):
    """Return apo rows unbound at the drug's pocket, deduped to best row per apo PDB."""
    drug = (drug_id or "").upper()
    rows = [e for e in (found_apo or []) if isinstance(e, dict)]

    lig_rows = [e for e in rows if e.get("target_poi_type") == "lig"]
    matched = [e for e in lig_rows if poi_ligand(e) == drug]
    if matched:
        source = "drug-matched"
        pool = matched
    elif lig_rows:
        source = f"any-ligand-pocket (drug '{drug}' not among detected)"
        pool = lig_rows
    else:
        source = "any-pocket"
        pool = rows

    best = {}
    for e in pool:
        pid = (e.get("pdb_id") or "").lower()[:4]
        if not pid or pid == holo_pdb.lower()[:4]:
            continue
        if pid not in best or rank_key(e) < rank_key(best[pid]):
            best[pid] = e
    ranked = sorted(best.values(), key=rank_key)
    return ranked, source


# ---------------- download AHoJ aligned structure ----------------
def download_apo(entry, holo_pdb):
    apo = (entry.get("pdb_id") or "").upper()[:4]
    url_rel = entry.get("structure_file_url")
    if url_rel:
        url = f"{SITE_ROOT}/{url_rel.lstrip('/')}"
        code, body = http(url)
        if code == 200 and body:
            path = os.path.join(DIR_APO, f"{apo}_aligned_to_{holo_pdb.upper()}.cif")
            with open(path, "wb") as f:
                f.write(body)
            return path, "ahoj_aligned_cif"
    # fallback: raw structure from RCSB (NOT aligned to the holo)
    for ext, kw in ((".cif", "cif"), (".pdb", "pdb")):
        code, body = http(f"https://files.rcsb.org/download/{apo}{ext}")
        if code == 200 and body:
            path = os.path.join(DIR_APO, f"{apo}{ext}")
            with open(path, "wb") as f:
                f.write(body)
            return path, f"rcsb_{kw}_unaligned"
    return None, "download_failed"


# ---------------- MAIN ----------------
def main():
    if not os.path.exists(INPUT_FILE):
        print(f"Error: '{INPUT_FILE}' not found."); sys.exit(1)
    os.makedirs(DIR_APO, exist_ok=True)
    with open(INPUT_FILE) as f:
        master = json.load(f)

    holos = []  # {pdb, ligand, class, protein}
    for p_class, structs in master.items():
        for s in structs:
            holos.append({"pdb": s["pdb_id"], "ligand": s["drug_id"],
                          "class": p_class, "protein": s.get("protein", "Unknown")})
    print("=" * 90)
    print(f"   AHoJ FINAL PAIRING — {len(holos)} holo monomers")
    print("=" * 90)

    # Phase 1: submit everything
    print("\nSTEP 1: submitting jobs...")
    for h in holos:
        h["job_id"] = submit(h["pdb"])
        print(f"  {h['pdb']} ({h['ligand']}) -> {h['job_id']}")
        time.sleep(0.4)

    # Phase 2: collect + filter + rank candidates per holo
    print("\nSTEP 2: waiting for results and ranking candidates...")
    for h in holos:
        h["candidates"], h["cand_source"] = [], "no-job"
        if not h["job_id"]:
            continue
        q = wait_done(h["job_id"])
        if not isinstance(q, dict):
            print(f"  {h['pdb']}: timed out"); continue
        ranked, source = candidates_for_drug(q.get("found_apo"), h["ligand"], h["pdb"])
        h["candidates"], h["cand_source"] = ranked, source
        print(f"  {h['pdb']}: {len(ranked)} unique apo ({source})")

    # Phase 3: assign ONE UNIQUE apo per holo.
    # Order by fewest options first so constrained holos aren't starved.
    print("\nSTEP 3: assigning unique apo per holo...")
    order = sorted(holos, key=lambda h: (len(h["candidates"]) or 10**6,
                                         rank_key(h["candidates"][0]) if h["candidates"] else ()))
    used_apo = set()
    pairs, unmatched = [], []
    for h in order:
        pick = None
        for e in h["candidates"]:
            pid = (e.get("pdb_id") or "").lower()[:4]
            if pid not in used_apo:
                pick = e
                used_apo.add(pid)
                break
        if pick:
            pairs.append((h, pick))
        else:
            unmatched.append(h)

    # Phase 4: download aligned apo + write outputs
    print("\nSTEP 4: downloading aligned apo structures...")
    apo_dataset = {c: [] for c in master.keys()}
    csv_rows = []
    for h, e in pairs:
        apo = (e.get("pdb_id") or "").upper()[:4]
        path, kind = download_apo(e, h["pdb"])
        print(f"  {h['pdb']} ({h['ligand']}) -> apo {apo}  [{kind}]")
        rec = {
            "target_class": h["class"],
            "protein_desc": h["protein"],
            "holo_pdb": h["pdb"].upper(),
            "holo_ligand": h["ligand"],
            "apo_pdb": apo,
            "apo_file": path,
            "apo_source": kind,
            "candidate_source": h["cand_source"],
            "resolution": num(e.get("resolution"), None),
            "pocket_rmsd": num(e.get("pocket_rmsd"), None),
            "backbone_rmsd": num(e.get("rmsd"), None),
            "tm_score": num(e.get("tm_score"), None),
            "mapped_binding_residues_percent": num(e.get("mapped_binding_residues_percent"), None),
            "mapped_binding_residues_unobserved": num(e.get("mapped_binding_residues_unobserved"), None),
            "apo_pocket": e.get("pocket"),
            "uniprot": (e.get("uniprot_ids") or [None])[0],
            "n_candidates": len(h["candidates"]),
        }
        apo_dataset[h["class"]].append(rec)
        csv_rows.append(rec)

    with open(OUT_JSON, "w") as f:
        json.dump(apo_dataset, f, indent=4)
    if csv_rows:
        with open(OUT_CSV, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=list(csv_rows[0].keys()))
            w.writeheader()
            w.writerows(csv_rows)

    print("\n" + "=" * 90)
    print(f"DONE — {len(pairs)} unique apo-holo pairs; {len(unmatched)} holo without a free apo")
    if unmatched:
        print("  unmatched holo:", ", ".join(h["pdb"] for h in unmatched))
    print(f"  pairs -> {OUT_JSON} and {OUT_CSV}")
    print(f"  aligned apo structures -> ./{DIR_APO}/")
    print("=" * 90)


if __name__ == "__main__":
    main()
