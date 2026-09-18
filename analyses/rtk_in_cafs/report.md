# Do PDAC CAFs express the daraxonrasib-induced RTKs? (baseline check)

**Question.** Daraxonrasib (RMC-6236) treatment upregulates RTKs (EGFR, MET,
ERBB2/HER2, ERBB3) on PDAC tumor cells (RevMed AACR 2026, PR009). Step zero
for the resistant-CAF hypothesis: do cancer-associated fibroblasts express
these receptors at all in untreated human PDAC? If CAFs don't express them,
an RTK-upregulation mechanism in stroma is a non-starter.

**Answer: yes, but weakly.** All four RTKs are detectable in fibroblasts
across 7–9 independent human PDAC scRNA-seq datasets, at levels consistently
below malignant cells (malignant/fibroblast ratios ≈ 2–23×).

| gene | Fibroblasts (mean) | Malignant (mean) | ratio M/F |
|------|-------------------|-----------------|-----------|
| EGFR  | 0.08 | 0.19 | 2.3× |
| ERBB2 | 0.14 | 0.36 | 2.6× |
| ERBB3 | 0.05 | 1.05 | 23× |
| MET   | 0.01 | 0.28 | 20× |

Values: mean log-normalized expression per major-lineage cell type,
averaged across datasets (TISCH2). ERBB3 and MET show the largest
malignant-vs-fibroblast gaps — both are barely detectable in fibroblasts,
so a treatment-induced fold-change there would be easy to miss without
sensitive assays. EGFR and ERBB2 have the most robust baseline fibroblast
signal (detectable in all 7–9 datasets), making them the better candidates
for testing the CAF-upregulation hypothesis in vitro.

**Implication for the project.** The RTK-upregulation mechanism is
mechanistically plausible in stroma, not just tumor cells. If chronic
daraxonrasib exposure raises EGFR/ERBB2 on CAFs the way it does on tumor
cells, the conditioned-media experiment (resistant-CAF CM → drug-treated
PDAC cells) is the right next test. Prioritize EGFR/ERBB2 readouts; for
MET/ERBB3 expect low baseline and plan for high-sensitivity detection.

## Methods

- Source: TISCH2 (http://tisch.comp-genomics.org), PAAD collection,
  `Celltype_curated` (major-lineage) annotation.
- Queried the site's own `search-gene` API for EGFR, MET, ERBB2, ERBB3 with
  `plottype=heatmap`; parsed the returned per-dataset × per-cell-type
  expression matrix (7–9 datasets per gene; some genes absent from some
  datasets).
- The CELLxGENE Census S3 endpoint was unreachable from this sandbox, so
  TISCH2's precomputed values were used instead of raw counts.
- Reproduce: `python analyses/rtk_in_cafs/run.py` (writes `output/`).

## Outputs

- `output/rtk_by_celltype_all.csv` — full gene × dataset × cell-type table
- `output/rtk_fibroblast_vs_malignant.csv` — fibroblast vs malignant summary
- `output/rtk_fibroblast_vs_malignant.png` — bar plot
- `output/meta.json` — query provenance

## Caveats

- TISCH2 values are precomputed means; no access to per-cell distributions
  or differential-expression statistics.
- One dataset (PAAD_GSE162708) shows inverted fibroblast/malignant ordering
  for EGFR/ERBB2 — likely an annotation or cell-count artifact; the mean
  across datasets is robust to it.
- Baseline expression ≠ inducibility. This check only clears the
  prerequisite (receptors present); the upregulation itself must be tested
  in her daraxonrasib-treated CAF time course.
