# Data Dictionary

Column definitions for the main supplementary data tables.

## offtarget_hits_high.tsv — the 58 high-confidence candidates

| Column | Meaning |
|--------|---------|
| `query_target` | PDB ID of the query drug-target complex |
| `query_conf` | query conformation used (holo or apo) |
| `hit_uniprot` | UniProt accession of the candidate off-target |
| `is_self` | whether the hit is the query protein itself (self-recovery control) |
| `idf` | FoldDisco inverse-document-frequency score (motif rarity weight) |
| `node_count` | number of matched residues in the motif |
| `min_rmsd` | geometric RMSD of the matched motif (Å) |
| `plddt` | AlphaFold per-residue confidence at the matched site |
| `n_matched` | matched residue count after filtering |
| `best_pocket_rank` | rank of the P2Rank pocket containing the match |
| `best_pocket_prob` | P2Rank druggability probability of that pocket |
| `overlap_count` | number of query residues falling inside the pocket |
| `overlap_frac` | fraction of query residues inside the pocket |
| `status` | pocket-support status |
| `confidence` | final confidence tier (high) |

## offtarget_hits_bioannotated.tsv — high-confidence hits with biological annotation

| Column | Meaning |
|--------|---------|
| `query_target` | PDB ID of the query complex |
| `query_gene`, `query_uniprot` | gene name and accession of the query protein |
| `drug_name` | the drug bound in the query structure |
| `drug_se_count` | number of recorded side effects (from SIDER) |
| `drug_source` | source of the drug annotation (SIDER or openFDA) |
| `hit_uniprot`, `hit_gene` | accession and gene of the candidate off-target |
| `family_class` | within-family or cross-family relative to the query |
| `idf`, `node_count`, `overlap_frac`, `confidence` | as in the high-confidence table |
| `chembl_target_id` | ChEMBL target identifier, if the hit is in ChEMBL |
| `is_known_drug_target` | whether ChEMBL lists the hit as a known drug target |
| `binding_status` | confirmed, weak, tested-inactive, no-measured-value, or unknown |
| `best_activity_nM` | best measured potency in nanomolar, if any |
| `activity_type` | assay type of the best measurement (Ki, Kd, IC50, EC50) |
| `evidence_tier` | 1 = confirmed binding, 2 = known target, 3 = structural only |

## phase5_recall.tsv — documented off-target recall

| Column | Meaning |
|--------|---------|
| `drug_target` | PDB ID of the query complex |
| `known_uniprot` | accession of a documented off-target for that drug |
| `justification` | source/evidence for the documented off-target |
| `confidence` | confidence in the documented annotation |
| `in_library` | whether the documented target is present in the search library |
| `found_at_high` | whether it was recovered among the high-confidence hits |
| `found_at_supported` | whether it was recovered among pocket-supported hits |

## folddisco_raw_hits.tsv — all raw FoldDisco matches

| Column | Meaning |
|--------|---------|
| `query_pdb` | PDB ID of the query complex |
| `conf` | query conformation (holo or apo) |
| `rank` | rank of the hit in the FoldDisco output (by IDF, best first) |
| `hit_uniprot` | UniProt accession of the matched protein |
| `idf` | FoldDisco IDF score for the match |

## blast_hits.tsv — BLAST+ baseline (standard outfmt 6)

Columns are the standard BLAST tabular format: query accession, subject accession,
percent identity, alignment length, E-value, and bit score.

## foldseek_hits.tsv — Foldseek baseline

Columns follow the Foldseek tabular output: query accession, target accession, E-value,
alignment bit score, sequence identity, and the three structural similarity measures
(query/target coverage and TM-score).
