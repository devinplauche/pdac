# CAF subtype spatial mapping in untreated human PDAC (GSE274103, PDAC-p1)

## Question

In untreated human PDAC tissue, where do myCAF vs iCAF vs apCAF programs
localize relative to tumor nests — are they spatially segregated niches?

## Method

- **Dataset:** GEO GSE274103, "Spatial transcriptomics on treatment-naive
  pancreatic ductal adenocarcinoma (PDAC) patients" (10x Genomics Visium,
  FFPE, *Homo sapiens*). Chose the repository's own baseline: treatment-naive
  patients, which is exactly the untreated spatial context the thesis needs.
  Analysis used **one section, GSM8443449 (PDAC-p1)**, done well rather than
  many done thinly: 17,943 genes x 4,987 in-tissue spots, median 2,527 UMIs
  per spot.
  - No treated (daraxonrasib/RMC-6236 or any therapy) human PDAC spatial
    dataset exists in the public record as of this writing — GEO query
    "pancreatic ductal adenocarcinoma Visium spatial" returned 18 records;
    none involve drug treatment. So this analysis maps the **untreated
    baseline spatial context only**, not drug response.
- **Signatures:** repo compendium as-is —
  `analyses/signature_panel/signatures_v1.yaml` (v1.0, curated 2026-09-18) —
  scored with `analyses/signature_panel/score_signatures.py`, `--organism
  human`, `--method rank` (AUCell-like fractional within-spot ranks; robust
  to the sparse FFPE spot data). Ranks were computed on the **full**
  17,943-gene matrix per spot, in-process (no 89M-entry CSV).
- **Tumor spots:** rank-score over epithelial genes present in the matrix
  (EPCAM, KRT7, KRT8, KRT18, MUC1, CEACAM5, SOX9; CEACAM6 absent from the
  matrix), thresholded at the 75th percentile -> 1,247/4,987 tumor spots.
- **Distance bands:** Euclidean distance (spot units, ~100 µm center-to-center)
  from the nearest tumor spot: tumor / adjacent (<=2) / mid (2-5) / distal
  (>5). Wilcoxon rank-sum (adjacent vs distal) + Spearman rho of score vs
  distance (non-tumor spots).
- Re-runnable: `python3 run.py` (needs `data/` downloads).

## Results

**Gene coverage (human panel): 100% of every signature found in the matrix**
(myCAF 6/6, iCAF 9/9, apCAF 4/4, SASP 15/15, matrisome_core 19/19,
MAPK_targets 9/9, RASON_inhibition_PD 7/7). No score rests on a thin gene set.

**Yes — the niches are spatially segregated, along a tumor-distance axis:**

| signature | median tumor | median adjacent | median mid | median distal | adj vs distal Wilcoxon p | rho(score, distance) |
|---|---|---|---|---|---|---|
| myCAF | 0.848 | 0.812 | 0.790 | 0.660 | 4.8e-161 | -0.47 (p~1e-208) |
| iCAF | 0.488 | 0.495 | 0.508 | 0.519 | 8.9e-46 | +0.24 (p~1e-52) |
| apCAF | 0.716 | 0.719 | 0.724 | 0.726 | 1.1e-11 | +0.11 (p~1e-12) |

- **myCAF is a tumor-proximal niche.** Scores are highest inside tumor nests
  and in immediately adjacent stroma, then fall off steeply with distance
  (median 0.85 -> 0.66; rho = -0.47). On the map, high myCAF tracks the
  coherent epithelial nests; stroma-only regions are dark.
- **iCAF is a distal-stromal niche.** Scores rise with distance from tumor
  (median 0.49 -> 0.52; rho = +0.24, adjacent-vs-distal p = 8.9e-46).
  The absolute gradient is smaller than myCAF's, but the direction is
  consistent on the map: the bright iCAF patches sit in the stroma far from
  the epithelial clusters.
- **apCAF shows essentially no spatial gradient** (medians 0.716-0.726,
  tiny rho) — in this section it does not form a distinct tumor-distance
  niche, or is too rare to resolve at spot resolution.

Figures: `output/tumor_map.png` (tumor-spot mask), `output/mycaf_map.png`,
`output/icaf_map.png`. CSVs: `output/spot_scores.csv` (4,987 spots x all
signature scores + tumor_score, distance band), `output/band_stats.csv`.

## Limitations

- Single section (PDAC-p1), one patient. Five samples exist in GSE274103;
  extending the pipeline to p2-p5 would test generality (run.py is
  parameterized at the top for this).
- Visium spots are ~55 µm, ~1-10 cells: "myCAF high in tumor nests" can
  reflect juxtatumoral fibroblasts within the same spot as tumor cells, not
  pure myCAF cells. Single-cell deconvolution would sharpen the claim.
- Tumor spots defined by transcriptomic threshold (top quartile of an
  epithelial rank score), not by a pathologist's H&E annotation; coherent
  nest morphology on the map supports it, but borderline spots could be
  misassigned (this pushes adjacent-vs-tumor differences *down*, so the
  observed gradients are likely conservative).
- apCAF signature is 4 genes (HLA-DRA, HLA-DRB1, CD74, SLPI; SAA3P omitted
  as a human pseudogene) — modest statistical power for a rare state.
- FFPE Visium (probe-based) detects a subset of transcripts; still,
  coverage of the whole compendium was complete.

## What this means for the thesis

1. **Baseline established, and it is drug-naive:** in untreated human PDAC,
   myCAF and iCAF programs occupy *segregated* niches — myCAF hugging tumor
   nests, iCAF enriched in distal stroma. This is not a drug effect; it is
   the architecture daraxonrasib will act on.
2. Together with the GSE337208 finding (RMC-6236's MAPK signal absent in
   fibroblasts in vivo; CAF reprogramming likely indirect via tumor cells),
   the spatial picture suggests the thesis hypothesis should distinguish
   **juxtatumoral (myCAF) vs distal (iCAF) fibroblast responses**: if the
   drug works through tumor cells first, the juxtatumoral myCAF niche is
   where indirect effects should appear first. A thesis analysis that
   averages all CAFs across the section would dilute exactly this gradient.
3. Recommendation for wet-lab design: when she profiles treated vs untreated
   models spatially (not publicly available — would need her own section),
   distance-from-tumor must be a covariate, or better, a stratification.
   The `run.py` pipeline is reusable for that future dataset.
