# CAF subtype signatures in bulk TCGA-PAAD: clinical correlates

## Question
In untreated human PDAC tumors, do bulk myCAF / iCAF / apCAF / SASP signature
scores associate with overall survival, pathologic stage, or immune
infiltration? A clinical-correlative angle for a thesis built on CAF-subtype
biology.

## Method
- **Cohort:** TCGA-PAAD via the GDC API (UCSC Xena downloads were blocked —
  S3 `AccessDenied` on every hub tried). 183 harmonized STAR-Counts
  RNA-seq files (178 primary tumor, 4 solid-tissue normal, 1 metastatic);
  clinical (vital status, days-to-death / days-to-last-follow-up,
  AJCC pathologic stage) from the GDC `/cases` endpoint, 185 cases.
  Analyses below use the 176 primary-tumor samples with OS data.
- **Scoring:** repo panel AS-IS (`analyses/signature_panel/signatures_v1.yaml`,
  `--organism human`, `--method zscore`), run with the panel's own
  `score_signatures.py` (log2(x+1), per-gene z across samples, mean z of
  member genes). Gene coverage was **100% for every signature**
  (myCAF 6/6, iCAF 9/9, apCAF 4/4, SASP 15/15; see `signature_coverage.csv`).
- **Immune/stromal readouts:** deliberately crude, naive mean-z scores —
  immune = PTPRC, CD3E, CD8A, CD4, MS4A1, CD68, NKG7, GZMB (8/8 present);
  stromal = COL1A1, COL3A1, DCN, PDPN, FAP, VIM (6/6 present). No
  deconvolution, no purity correction — labeled "naive" throughout.
- **Statistics:** Kaplan–Meier with median split + hand-rolled log-rank
  (Mantel–Haenszel, χ²(1); lifelines not installed); Kruskal–Wallis across
  collapsed stage groups (I/II/III/IV); Spearman correlations vs
  immune/stromal-naive scores. Effect sizes reported alongside p-values.
- **Reproducibility:** `run.py` rebuilds everything from `data/` (raw GDC
  files + clinical JSON live there; nothing in `output/` exceeds ~70 KB).

## Results
- **Survival — null.** Median-split log-rank for overall survival:
  iCAF χ²=0.26, **p=0.61** (median OS high 386d vs low 394d);
  myCAF χ²=0.16, **p=0.69** (high 468d vs low 308d);
  apCAF p=0.17; SASP p=0.65. n=88/88 per arm. KM curves
  (`fig_survival_km.png`) overlap. No CAF-subtype score stratifies OS in
  untreated bulk TCGA-PAAD.
- **Stage — null, and underpowered.** Kruskal–Wallis across stage groups:
  iCAF p=0.25, myCAF p=0.30 (`fig_stage_boxplot.png`). The cohort is 146/173
  stage II (84%), with n=3 stage III and n=4 stage IV — there is essentially
  no contrast to detect. This is a cohort-composition limit, not evidence
  of no biology.
- **Immune infiltration — strong positive correlations.** Spearman (n=176):
  iCAF vs immune-naive **r=0.69, p=7.2e-26**; apCAF vs immune-naive
  **r=0.76, p=9.2e-34**; SASP vs immune-naive r=0.57, p=1.5e-16;
  myCAF vs immune-naive r=0.49, p=3.4e-12. Against the stromal-naive score:
  myCAF r=0.93 (near-circular — myCAF markers are core stromal genes, so
  this mostly validates the score tracks stromal content), iCAF r=0.59,
  apCAF r=0.56, SASP r=0.65 (all p<1e-15). Heatmap in `fig_immune_corr.png`.
- **Cross-signature structure:** iCAF–SASP Spearman r=0.82 (shared IL6/CXCL
  chemokine genes — expected, not independent); myCAF correlates only
  modestly with iCAF (r=0.48) and apCAF (r=0.50), consistent with distinct
  fibroblast programs rather than one generic "stroma" axis.

## Limitations
1. **Bulk scoring cannot resolve cell-type proportions.** A high iCAF score
   can mean more iCAFs, more inflammatory signaling per fibroblast, or
   inflammatory tumor/epithelial cells expressing the same chemokines. All
   "CAF" claims here are really "CAF-signature-high bulk tumor" claims.
2. **TCGA-PAAD is untreated baseline, not daraxonrasib.** These patients
   received no RAS(ON) inhibition; the results say nothing about how
   daraxonrasib reshapes CAF states — they only describe the prognostic
   landscape the drug would act on.
3. **Survival is confounded and the split is crude.** Median split is
   arbitrary; no multivariable adjustment (grade, resection margin, adjuvant
   therapy, molecular subtype) was attempted. TCGA-PAAD OS is short and
   the cohort small (n=176, ~100 events) — modest effects could be missed.
4. **Immune/stromal scores are naive.** Mean-z of a handful of markers;
   no ESTIMATE, no deconvolution, no tumor-purity adjustment. The myCAF–
   stromal r=0.93 is partly definitional.
5. **Histologic heterogeneity.** TCGA-PAAD includes non-ductal cases
   (neuroendocrine, mucinous, etc.); analysis used all primary tumors.
   Stage analysis is dominated by resected stage-II disease.
6. **Data provenance.** GDC STAR-Counts (GENCODE v36) were used instead of
   the planned Xena HTSeq counts (Xena downloads denied); results are
   z-score based and should be robust to this, but the exact matrix is
   GDC-harmonized, not Xena-toil.

## What this means for the thesis
- The **null survival result is usable, not a failure**: bulk CAF-subtype
  scores do not stratify OS in untreated PDAC. That is the honest baseline
  against which any daraxonrasib-induced CAF reprogramming must be judged —
  if the drug shifts iCAF/myCAF balance, the clinical relevance of that
  shift cannot be assumed from untreated prognostic data, because none is
  detectable here.
- The **iCAF/apCAF–immune correlations are the positive finding**:
  inflammatory and antigen-presenting CAF programs track with immune
  infiltration signatures in bulk tumors. This supports framing CAF states
  as part of an inflamed-tumor ecosystem (relevant if daraxonrasib acts
  partly via microenvironmental reprogramming, per the GSE337208
  CAF-mining result), but bulk data cannot say whether CAFs recruit immune
  cells or vice versa.
- **Recommended next step:** test the same panel in a single-cell or
  deconvolved setting (or against published PDAC scRNA-seq CAF annotations)
  to check the bulk scores actually track CAF-subtype abundance before
  using them as drug-response readouts; and repeat the survival analysis
  with molecular-subtype (basal/classical) adjustment.

## Outputs
- `output/scores_clinical.csv` — per-sample scores + clinical merge
- `output/stats_summary.csv` — all test statistics, p-values, effect sizes
- `output/signature_coverage.csv` — gene coverage per signature
- `output/fig_survival_km.png`, `output/fig_stage_boxplot.png`,
  `output/fig_immune_corr.png`
- `output/run.log` — full run transcript
- `data/` — raw GDC downloads (~1 GB incl. tar), clinical JSON, file maps
  (**gitignore this directory**)
