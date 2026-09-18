# Mining GSE337208 for CAF subtype shifts under RMC-6236

## Question
Does the RAS(ON) inhibitor daraxonrasib (RMC-6236) shift CAF subtype identity
(myCAF / iCAF / apCAF) in vivo? GSE337208 (Mahadevan et al., *Nat Commun*
2026, PMID 42649211) is scRNA-seq of orthotopic murine KPC PDAC treated with
RMC-6236 or MRTX1133 ± anti-CTLA4/PD1 for 10–14 days. The paper's own analysis
is T-cell-centric; the fibroblast compartment was never scored — this analysis
fills that gap with the repo's validated signature panel
(`../signature_panel/signatures_v1.yaml`, compendium v1.0, mouse lists).

## Method
- Per-sample CellRanger count matrices from GEO (9 samples, 3 KPC models).
  No alignment; processed matrices only.
- QC per cell: ≥200 genes, ≥500 UMIs, ≤15% mitochondrial.
- Fibroblast identification: no author annotations exist, so k-means (k=6,
  deterministic seed) on a 16-gene panel (fibroblast: Col1a1/Col1a2/Col3a1/
  Dcn/Pdpn/Lum/Sparc; epithelial: Epcam/Krt19/Krt8/Krt18; immune: Ptprc/Cd3e/
  Cd68/Cd79a/Itgam), log-normalized. The fibroblast cluster was the one
  maximizing (fibroblast program − other programs); doublet cleanup dropped
  cells whose epithelial/immune program exceeded the fibroblast program.
- Scoring: repo `score_signatures.py` **used as-is** —
  (a) pseudobulk fibroblast counts per sample, z-score method across samples;
  (b) per-cell rank (AUCell-like) scores for the primary KPC2 RMC6236-vs-Veh
  pair, Mann–Whitney U across cells.
- Contrasts are strictly **within KPC model** (treatment is confounded with
  model across KPC1/KPC2/KPC3).

## Results

### Fibroblast yield per sample (the first surprise)
| Sample | Model / treatment | Fibroblasts | % of QC cells |
|---|---|---|---|
| GSM9850794 | KPC1 Veh | 1857 | 15.9% |
| GSM9850795 | KPC1 MRTX1133 | 195 | 2.6% |
| GSM9850796 | KPC3 anti-CTLA4 | 447 | 9.3% |
| GSM9850797 | KPC3 MRTX1133 | 596 | 12.4% |
| GSM9850798 | KPC3 MRTX1133+anti-CTLA4 | 383 | 5.6% |
| GSM9850799 | KPC2 RMC6236 | 53 | 0.9% |
| GSM9850800 | KPC2 Veh | 42 | 0.7% |
| GSM9850801 | KPC2 RMC6236+anti-CTLA4 | **0** | 0% |
| GSM9850802 | KPC2 RMC6236+anti-PD1 | **0** | 0% |

KRAS inhibition collapses the captured fibroblast compartment: 15.9%→2.6%
(KPC1, MRTX1133), and the RMC6236+checkpoint combos yielded **zero**
fibroblasts (Col1a1 detected in 0.02% of cells — the program is absent, not
just downregulated). Whether this is stromal depletion by effective therapy
or dissociation bias cannot be separated at n=1.

### Primary contrast: KPC2 RMC6236 vs Veh (53 vs 42 fibroblasts)
Cell-level rank scores, Mann–Whitney U (cell-level; pseudoreplicated —
1 sample/arm, see Limitations):

| Signature | median Veh | median RMC6236 | Δ | MWU p |
|---|---|---|---|---|
| myCAF | 0.626 | 0.548 | −0.078 | 0.012 |
| iCAF | 0.475 | 0.533 | +0.058 | 0.0005 |
| apCAF | 0.690 | 0.867 | +0.177 | 0.0002 |
| SASP | 0.553 | 0.550 | −0.004 | 0.84 |
| matrisome (core) | 0.849 | 0.854 | +0.005 | 0.73 |
| MAPK output | 0.627 | 0.627 | +0.000 | 0.42 |
| RAS(ON) inhibition PD | 0.622 | 0.625 | +0.003 | 0.94 |

RMC-6236 shifts fibroblasts **away from myCAF and toward iCAF/apCAF**
(iCAF and apCAF survive Bonferroni over 7 signatures; myCAF at p=0.012 does
not). SASP and the structural matrisome program do not move — the shift is
subtype identity, not senescence or bulk ECM loss. **MAPK-output and the
RMC-6236 PD readout (Dusp6/Spry/Etv) do not budge in fibroblasts** — the
drug's transcriptional PD signal is not detectable in CAFs at day 14.
Pseudobulk z-scores agree in direction (see `output/fig2_pseudobulk_heatmap.png`).

### Secondary contrasts (brief)
- **KPC1 MRTX1133 vs Veh:** all signatures fall in pseudobulk (myCAF,
  MAPK, PD, matrisome), but the Veh arm's fibroblast pool is 10× larger
  (1857 vs 195 cells) — compositional confounding makes this descriptive only.
- **KPC3 MRTX1133 vs anti-CTLA4:** modest iCAF/apCAF elevation under
  MRTX1133; the **MRTX1133+anti-CTLA4 combo suppresses apCAF (−1.25),
  MAPK-output (−0.91) and PD (−0.86)** in the remaining fibroblasts.
- KPC3 has no clean vehicle arm, so these are hypothesis-generating.

## Figures
- `output/fig1_kpc2_cell_scores.png` — per-cell signature distributions,
  KPC2 Veh vs RMC6236.
- `output/fig2_pseudobulk_heatmap.png` — fibroblast pseudobulk z-scores,
  all 7 samples (compare within model only).

## Limitations (read before citing)
1. **n=1 sample per arm.** Cell-level p-values are pseudoreplicated; there is
   no sample-level statistical inference anywhere here. Treat effect
   directions as observations, not findings.
2. **No author cell-type annotations.** Fibroblasts were identified by
   unsupervised clustering on a marker panel — clean clusters (fibroblast
   program mean 2.5–3.3 vs ≤0.4 for other programs), but a heuristic.
3. **Small fibroblast n in KPC2** (42–53 cells) — the primary contrast's
   home turf. KPC1/KPC3 have hundreds.
4. **Dropout:** iCAF genes are sparsely detected (Il6 26%, Cxcl1 17–25%,
   Ly6c1 2–19%, Has1 13–26% of fibroblasts); Mmp1a is absent from the
   annotation entirely; Csf2/Cxcl3 near-absent. The iCAF score rides on few
   genes per cell — the signature panel's own known weak link.
5. **10–14 day treatment** = short-term response, not acquired resistance.
6. **Zero fibroblasts in combo arms** — stromal depletion vs dissociation
   artifact is unresolvable here.
7. Different CellRanger references per sample batch; pseudobulk aligned on
   gene-symbol union.

## What this means for her daraxonrasib–CAF thesis
1. In the only public in-vivo dataset with CAFs under daraxonrasib, the drug
   pushes fibroblasts from a myofibroblastic toward inflammatory and
   antigen-presenting states — the same myCAF→iCAF direction reported for
   KRAS inhibition in appendiceal models (GSE336609), now in PDAC itself.
2. The drug's own MAPK pharmacodynamic readout does not register in
   fibroblasts, so CAF reprogramming here is likely indirect (tumor-cell
   signaling or microenvironmental cues), not direct MAPK suppression in CAFs.
3. Effective KRAS-inhibitor combinations nearly erase the captured
   fibroblast compartment — any stromal readout in her own combination
   experiments needs dissociation-bias controls, or she may mistake
   capture loss for biology.
4. These are single-sample observations, not proof: her own time-course
   CAF RNA-seq with replicates is exactly what would turn this hint into
   a result, and the validated scorer in this repo is ready to score it.

## Reproducibility
`python3 run.py` re-runs everything (expects `data/*.mtx.gz` from GEO;
`dl.sh` documents the download). `make_figures.py` regenerates the figures.
