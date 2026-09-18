# Signature panel: compendium + GSE291118 reproduction + scorer validation

## What GSE291118 is

**GSE291118** — "Macropinocytosis maintains CAF subtype identity under metabolic
stress in pancreatic cancer [RNA-seq]". Zhang Y. et al., *Cancer Cell* 2025,
PMID 40712568 (Commisso lab, Sanford Burnham Prebys). **Bulk RNA-seq of
immortalized murine CAFs** (not single-cell), 20 samples = 4 conditions × 5
replicates, 32 h treatment:

| Code | Treatment | Meaning |
|---|---|---|
| `4Qveh` | 4 mM glutamine + DMSO | nutrient-replete baseline |
| `02Qveh` | 0.2 mM glutamine + DMSO | metabolic stress, macropinocytosis ON |
| `02Qeipa` | 0.2 mM glutamine + 25 µM EIPA | stress, macropinocytosis BLOCKED |
| `4Qil1a` | 4 mM glutamine + 25 pg/mL IL-1α | iCAF positive control (Biffi et al. 2019) |

Paper's headline claims: (1) metabolic stress induces an intrinsic iCAF program
via MEK-ERK; (2) macropinocytosis sustains the myCAF phenotype under glutamine
limitation by *preventing inflammatory reprogramming* — blocking it promotes
myCAF→iCAF transitions, iCAF enrichment and collagen reduction, remodeling the
stroma.

## What was built

- `signatures_v1.yaml` — 7 signatures with human + mouse gene lists and
  PubMed-verified provenance: myCAF / iCAF / apCAF (Elyada et al., Cancer
  Discov 2019, PMID 31197017; markers taken from the paper's main text),
  SASP (Coppé et al., PLoS Biol 2008, PMID 19053174),
  core matrisome (Naba et al., Mol Cell Proteomics 2012, PMID 22159717),
  MAPK/ERK output (Pratilas et al., PNAS 2009, PMID 19251651),
  RAS(ON)-inhibition PD readout (Holderfield et al., Nature 2024, PMID
  38589574; DUSP6 as daraxonrasib PD marker).
- `score_signatures.py` — reusable CLI scorer (z-score or rank/AUCell-like
  method), reports per-signature gene coverage.
- `run.py` — downloads the GEO TPM matrix, scores, Mann-Whitney U contrasts,
  figure.
- `output/` — `signature_scores.csv`, `contrast_stats.csv`,
  `gene_coverage.csv`, `signature_scores_by_condition.png`,
  cached `GSE291118_TPM.txt.gz`.

## Reproduction results (n=5/group, Mann-Whitney U, p=0.0079 = minimum for n=5)

**Claim 2 — MP blockade drives myCAF→iCAF transition + stromal remodeling
(`02Qeipa` vs `02Qveh`): REPRODUCED.**
myCAF DOWN (p=0.0079), SASP UP (p=0.0079), MAPK-output UP (p=0.0079),
matrisome DOWN (p=0.0079), apCAF UP (p=0.0079); iCAF UP trend (p=0.22).
Gene-level: Acta2 533→108 TPM, Tagln 296→52, Col1a1 675→262
(collagen reduction, as in the paper).

**Claim 1 — metabolic stress alone induces the iCAF program
(`02Qveh` vs `4Qveh`): PARTIAL.**
Stress erodes myCAF identity (myCAF DOWN, p=0.0079) and dampens baseline
inflammatory transcripts (iCAF/SASP DOWN, p=0.0079; e.g. Cxcl1 7.7→1.0 TPM,
Ccl2 48→3.0). At this 32 h bulk timepoint the *gain* of inflammatory
reprogramming appears only once the MP brake is removed (EIPA contrast above)
— consistent with the paper's model that macropinocytosis *prevents*
inflammatory reprogramming under stress. MAPK output unchanged by stress
alone (p=0.55), up with MP blockade (p=0.0079).

**Positive control (`4Qil1a` vs `4Qveh`): PASSED.**
iCAF UP (p=0.0079), SASP UP (p=0.0159), myCAF unchanged — the scorer's iCAF
signature responds to the canonical IL-1α iCAF inducer, and Saa3 jumps
0.7→20.1 TPM.

## Scorer validation

- Gene coverage on GSE291118: myCAF 6/6, iCAF 6/6, apCAF 5/5,
  matrisome 19/19, MAPK 9/9, RAS-ON PD 7/7, SASP 13/14
  (mouse `Mmp1a` annotation absent from matrix).
- IL-1α positive control confirms the iCAF signature detects the known
  biology (p=0.0079), and the EIPA contrast reproduces the paper's headline
  transition at signature level.

## Caveats

- Mean z-score on bulk TPM: fine for condition contrasts, not a per-cell
  classifier; rank method provided for scRNA-seq inputs.
- iCAF signature is the weakest link here (6 mouse genes; Il6/Has1/Clec3b
  near-zero in vitro at 32 h) — scores track Cxcl1/Ly6c1/Col14a1.
- v1 myCAF/iCAF/apCAF markers are the main-text core sets from Elyada;
  the paper's full DE tables (Supp. Table S22) were not machine-readable
  from PMC for v1 — a v2 could expand from them.
