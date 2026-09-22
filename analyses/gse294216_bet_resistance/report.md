# Transcriptional program of BET-inhibitor-vulnerable MRTX1133-resistant PDAC (GSE294216)

## Question
What is the transcriptional program distinguishing the MRTX1133-resistant state in human PDAC (Principe/Becker et al. (Munshi), *Mol Cancer Ther* 2026), and does it intersect (a) histone-acetylation machinery — readers (BRD2/3/4), writers (EP300/CREBBP), erasers (HDACs) — and (b) the repo's resistance themes: RTK upregulation, EMT, MYC?

## Method
- **Data:** GEO processed `featureCounts` count matrices — no alignment, processed only. Primary contrast: human **PANC1 (parental, n=3) vs PANC1K (MRTX1133-resistant, n=3)**, grown 48 h in 3D floating collagen. Matrices are already normalized (fractional counts); used as supplied.
- **DEG:** per-gene Welch two-sample t-test on log2(x+1) counts, BH-FDR. Calls: |log2FC| ≥ 1, p < 0.05 (n=3/group makes Mann–Whitney useless — two-sided min p = 0.1 — so it is documented but not used; FDR<0.10 numbers also reported for strictness).
- **Pathway summary:** hypergeometric overlap of DEG-up / DEG-down sets against inline hand-built Hallmark-style sets (EMT, MYC targets, RTK, MAPK, BCL2/apoptosis, FOSL1/AP-1 program, HAT writers, BET readers, HDAC erasers), over the detected-gene background (~26,475 genes).
- **Targeted tests:** BRD2/3/4, EP300, CREBBP, HDAC/SIRT family, FOSL1 axis.
- **BET-reversal check:** within-batch t-tests of PANC1K + JQ1 (n=3) and PANC1K + OTX-015 (n=3) vs untreated PANC1K (n=3); sign-concordance of treatment log2FC against the resistance DEG program (binomial vs 50%). Plus an EP300-inhibition (SGC-CBP30) spot-check of paper-relevant genes.
- Re-runnable: `python3 run.py`. Figures: `output/volcano_resistant_vs_parental.png`, `output/heatmap_acetylation_machinery.png`.

## Results

### 1. The resistance program: RTK rewiring + MAPK-feedback collapse, not EMT/MYC
- **993 genes up / 1,119 genes down** in resistant vs parental (p<0.05, |log2FC|≥1; stricter FDR<0.10: 262 up / 381 down).
- **RTK set strongly enriched among UP genes** (p = 1.4e-6, fold-enrichment 7.7): ERBB4, **FGFR2, FGFR3, FGFR4**, FLT1, FLT4, INSR, KIT, NTRK2. Meanwhile **EGFR itself is down** (−0.73, p=0.013) and AXL/KDR/PDGFRB/RET are down — this is a *swap* of RTK dependencies (FGFR axis on, EGFR-family off), not a generic RTK amplification. Matches the repo's RTK-upregulation resistance theme, with a specific FGFR flavor.
- **MAPK feedback genes collapse** (enriched among DOWN, p = 3.2e-5): SPRY2 (−4.35), SPRY4, DUSP5, DUSP6, ETV4, ETV5, FOS. Top down hits overall: FGF5 (−8.2), WISP1 (−8.1), FN1 (−5.1), SPARC (−3.0). Consistent with resistant cells no longer running MAPK feedback loops.
- **No classic EMT:** EMT set is enriched among *down* genes (p = 3.1e-4: FN1, MMP2, SPARC, TGFBI, PDGFRB, ZEB2); VIM −0.63, CDH1 −0.54 both down; only COL3A1/COL5A1/MMP9 up. **MYC program absent:** MYC mRNA down (−0.58, p=0.016); MYC-target set has zero DEG overlap.
- **BCL2-family:** BCL2 mRNA up +3.9 but noisy (p=0.067); BCL2L1 and MCL1 flat-to-down.

### 2. Acetylation machinery: rebalancing, not induction of EP300/BRD
- Core nodes are **essentially unchanged at mRNA**: EP300 −0.08 (ns), CREBBP −0.19 (ns), BRD2 −0.06 (ns), BRD3 −0.18 (ns), BRD4 −0.25 (nominal p=0.045, fdr=0.35 — weak).
- Instead, a **bidirectional eraser/writer rebalance**: UP — HDAC11 (+0.98, p=0.009), HDAC5 (+0.89, p=0.0017), SIRT7 (+0.72), SIRT2 (+0.65), reader CECR2 (+1.03, p=9.5e-4); DOWN — HDAC9 (−2.25, p=1.9e-4), HDAC7 (−1.09, p=0.0024), KAT2B (−1.33, p=0.016), KAT6B (−1.02, p=0.0074), HDAC8 (−0.79, p=0.0063).
- Reading: the "global shift toward histone acetylation" reported in the paper is **not driven by transcriptional upregulation of EP300 or BET readers** — it must be a net-activity phenomenon (post-translational / compensatory, e.g. the HDAC up-shifts read as compensation). The BET vulnerability is functional: resistant cells depend on reader activity without overexpressing the readers.

### 3. FOSL1 axis — a cautionary note
- FOSL1 mRNA is **down** in resistant vs parental (−1.86, p=0.0035), as are FOS, JUN, DUSP1/6. This does *not* contradict the paper: their claim is a *dependency* (siFOSL1 / EP300i resensitize resistant cells), and dependency ≠ steady-state mRNA direction. Notably, both BET inhibitors further suppress FOSL1 in resistant cells (JQ1 −0.36, p=0.042; OTX-015 −0.58, p=0.013) and EP300 inhibition nudges it down (−0.20, ns), placing FOSL1 inside the BET-dependent resistance program.

### 4. BET inhibitors invert the resistance-up program (both drugs agree)
- JQ1 reverses **844/993** resistance-up DEGs (binomial p ≈ 1e-118); OTX-015 reverses **862/993** (p ≈ 1e-133). The top resistance-up gene COL3A1 (+6.6) is crushed by both (JQ1 −3.08, p=1e-5; OTX −2.87, p=3e-5); FGFR3/FGFR4 are significantly suppressed by both BET inhibitors (OTX: FGFR3 −0.58 p=0.0013, FGFR4 −0.78 p=0.0011). Resistance-down genes are *not* restored (only ~14% move back up) — consistent with BET inhibition being broadly transcriptionally repressive rather than a full state reversion. Caveat: sign-concordance was computed across separate experimental batches.

### 5. H358 RMC-6236 comparison — future work
`analyses/resistance_signatures/output/` does not exist yet (only `data/` with raw series matrices), so the H358 RMC-6236-resistance signature comparison is **future work**, not run here. The GSE324125 series matrix for H358 is already in the repo and is the source to build it from.

## Limitations
- n=3/group: Welch t is the honest choice; nominal-p DEG sets are used for overlaps (FDR<0.10 calls are much smaller: 262/381). Treat exact gene ranks as soft, enrichment patterns as firm.
- Bulk RNA-seq (3D floating collagen monoculture): no cell-state deconvolution; tumor-cell-intrinsic program only.
- Matrices are pre-normalized featureCounts; normalization details undocumented beyond the file format.
- Gene sets are hand-defined inline approximations of Hallmark sets, not the MSigDB collections.
- Reversal analysis compares across experimental batches (within-batch treatment contrasts, cross-batch sign concordance).

## What this means for the thesis
- **MRTX1133, not daraxonrasib**; **human PDAC cell line**, bulk, in vitro. The GSE337208 finding — RMC-6236's MAPK signal not registering in fibroblasts in vivo, implying indirect CAF reprogramming — is untouched and complementary: this result is the *tumor-cell-intrinsic* side of KRAS-inhibitor resistance.
- The resistant state is an **RTK-rewiring program (FGFR2/3/4 up, EGFR down) with MAPK-feedback collapse** — convergent with the repo's RTK-upregulation resistance theme, and notably **not** an EMT or MYC program. Testable follow-ups: (a) FGFR inhibition as a resistance-breaking combo alongside BET inhibition; (b) measure HDAC/HAT *activity* and H3K27ac rather than machinery mRNA, since the mRNA shifts look compensatory; (c) in the co-culture experiment already recommended (tumor–CAF under drug), stain/score whether the FGFR-axis up program appears in tumor cells while CAF reprogramming stays indirect.
- BET-inhibitor vulnerability is validated at the program level: ~85% of the resistance-up transcriptome inverts under JQ1 or OTX-015. If daraxonrasib resistance shares this RTK-rewiring axis, BET combinations remain rational — but that cross-drug claim needs the H358 signature comparison (future work) and, ideally, a daraxonrasib-resistant PDAC model.
