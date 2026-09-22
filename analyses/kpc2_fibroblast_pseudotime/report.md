# Analysis B — Pseudotime of the KPC2 fibroblast myCAF → iCAF transition

## Question

Can the 95 KPC2 fibroblasts (53 RMC6236, 42 vehicle) be ordered along a
myCAF → iCAF axis, and which genes change along it?

## Method

- Reused the sibling cell-level signature scores
  (`../gse337208_caf_mining/output/kpc2_cell_signature_scores.csv`, mean
  fractional ranks). Reloaded both KPC2 matrices with the **identical**
  QC/normalization (n_genes ≥ 200, UMI ≥ 500, pct_mt ≤ 15; CP10k, log1p),
  re-identified fibroblasts with the sibling k-means rule, and mapped the
  95 scored cells back to filtered-matrix column indices (all 95 matched).
- **Pseudotime = iCAF score − myCAF score** per cell; cells binned into 4
  quartile windows (W1 = myCAF-like → W4 = iCAF-like). Sensitivity check:
  PC1 of PCA on the signature-gene expression matrix, correlated against
  the score-difference ordering (r = 0.57 — moderate agreement, so the
  assumed axis is a simplification, not a discovered trajectory).
- For every gene detected in ≥ 10% of the 95 fibroblasts (12,007 genes):
  Spearman ρ of log-normalized expression vs pseudotime rank, with nominal
  p-values and Benjamini–Hochberg FDR for reference. Hypothesis-generating
  only.

Run: `python3 run.py`. Outputs in `output/`.

## Results

Window means validate the ordering: myCAF 0.77 → 0.49, iCAF 0.47 → 0.61
from W1 to W4 (`output/fig1_pseudotime_scores.png`).

159 genes pass FDR < 0.05. Full table: `output/gene_trends.csv`;
window means for top genes: `output/window_means_top_genes.csv`;
trajectories: `output/fig2_top_gene_trends.png`.

### Top genes rising toward the iCAF end (ρ > 0)

| gene | ρ | FDR | notes |
|---|---|---|---|
| Clec3b | 0.53 | 5.5e-05 | iCAF signature member (expected) |
| Add3 | 0.50 | 2.6e-04 | actin capping, cytoskeletal remodeling |
| Dcn | 0.50 | 2.6e-04 | decorin; pan-fibroblast but iCAF-skewed |
| Nkain4 | 0.47 | 8.5e-04 | Na+/K+ ATPase interactor |
| Gbp2 | 0.45 | 2.1e-03 | **interferon-inducible GTPase** |
| Cd34 | 0.45 | 2.2e-03 | progenitor/stromal marker |
| H2-K1 | 0.45 | 2.3e-03 | **MHC class I** |
| Fth1 | 0.44 | 2.4e-03 | ferritin heavy chain |
| C3 | 0.43 | 4.3e-03 | **complement C3** |
| C1ra | 0.42 | 7.1e-03 | **complement C1r** |
| C1s1 | 0.41 | 9.2e-03 | **complement C1s** |
| Aldh1a1 | 0.41 | 8.7e-03 | retinoic-acid synthesis |
| Psmb10 | 0.41 | 8.8e-03 | immunoproteasome subunit |
| Fras1 | 0.41 | 7.9e-03 | ECM |
| Clec2d | 0.41 | 9.4e-03 | C-type lectin |

Notable pattern: the iCAF end is marked by **complement components
(C3, C1ra, C1s1)** and **interferon-associated genes (Gbp2, H2-K1,
Psmb10)** — an interferon/complement-flavored inflammatory program,
plus Ly6c1 (ρ = 0.34, FDR 0.05) and Has1 (ρ = 0.27).

Strikingly, **the canonical iCAF cytokines do NOT track the axis**:
Il6 ρ = −0.07 (p = 0.53), Cxcl1 ρ = −0.02 (p = 0.82). The iCAF-score rise
is carried by Clec3b/Ly6c1/Has1/complement, not by Il6/Cxcl1 — either
dropout of low-abundance cytokines or a genuinely Il6-independent
inflammatory flavor.

### Top genes falling toward the myCAF end (ρ < 0)

| gene | ρ | FDR | notes |
|---|---|---|---|
| Tagln | −0.63 | 9.7e-08 | myCAF signature member (expected) |
| Ncam1 | −0.61 | 2.3e-07 | neural cell adhesion molecule |
| Ndufa4l2 | −0.57 | 7.9e-06 | hypoxia-associated mitochondrial gene |
| Col12a1 | −0.56 | 8.7e-06 | myCAF signature member (expected) |
| Thbs2 | −0.55 | 1.3e-05 | myCAF signature member (expected) |
| Ptprn | −0.52 | 1.1e-04 | receptor phosphatase |
| Acta2 | −0.51 | 1.6e-04 | myCAF signature member (expected) |
| Tnc | −0.51 | 2.1e-04 | tenascin-C, ECM |
| Tpm4 | −0.49 | 4.6e-04 | tropomyosin, contractile apparatus |
| Slc16a3 | −0.49 | 5.5e-04 | lactate transporter (MCT4) |
| Tmem45a | −0.48 | 7.4e-04 | transmembrane protein |
| Col5a1 | −0.48 | 7.4e-04 | fibrillar collagen |
| Mmp14 | −0.47 | 8.4e-04 | membrane MMP |
| Tln2 | −0.47 | 8.4e-04 | talin-2, focal adhesions |
| Mmp2 | −0.47 | 8.5e-04 | gelatinase |

The myCAF end loses the contractile apparatus (Tagln, Acta2, Tpm4, Tln2),
focal-adhesion/ECM-remodeling machinery (Col12a1, Col5a1, Thbs2, Tnc,
Mmp2, Mmp14), plus Ncam1 and the hypoxia-linked Ndufa4l2.

## Limitations

- **95 cells, no real trajectory inference.** The myCAF → iCAF axis is
  *assumed* from the score difference, not learned from branching
  structure; the PC1 sensitivity check agrees only moderately (r = 0.57).
  There is no evidence these cells form a single continuous trajectory —
  they may be two overlapping populations.
- **Pseudotime is heavily confounded with treatment.** Window composition:
  W1 = 19 Veh / 5 RMC6236; W4 = 20 RMC6236 / 4 Veh. The "trends" therefore
  blend CAF state with drug response; a gene "rising along pseudotime"
  may simply be drug-induced. This analysis cannot separate the two.
- Il6/Cxcl1 dropout: absence of a trend for these cytokines may be
  technical (low-abundance transcripts) rather than biological.
- Nominal p-values from 12,007 Spearman tests on 95 non-independent cells;
  FDR values are reference only.

## What this means for the thesis

Even with all caveats, the gene lists give concrete, testable structure
to the iCAF shift beyond the signature scores:

1. **The drug-associated iCAF state has a complement/interferon flavor**
   (C3, C1ra, C1s1, Gbp2, H2-K1, Psmb10) rather than being purely
   Il6/Cxcl1-driven — at the single-cell level Il6 and Cxcl1 don't even
   track the axis. For the thesis, this reframes the iCAF program to test
   in co-culture: read out **complement and interferon-response genes**,
   not just Il6/Cxcl1, when asking whether tumor cells or drug drive the
   shift.
2. **The myCAF loss is a coordinated disassembly** of the contractile +
   ECM-remodeling program (Tagln/Acta2/Tpm4/Tln2/Mmp2/Mmp14/Col12a1) —
   consistent with a TGFβ-maintenance signal being withdrawn rather than
   an inflammatory signal being added (ties back to Analysis A's negative
   ligand result).
3. **Candidate validation targets**, in priority order: Ncam1 and
   Ndufa4l2 (strongest non-signature myCAF-end markers — is Ncam1 loss
   causal or correlative?); Gbp2/Irf-axis (is there an interferon input
   upstream of the complement program? Irf7 trends weakly, ρ = 0.29);
   C3 (does fibroblast-derived C3 feed back on tumor cells?).

Concrete experiment this motivates: tumor–CAF co-culture under RMC-6236
with an IL-1R blockade arm and a TGFβR-inhibition arm, reading out the
top-15 gene panels from each direction by qPCR — the panels here
(`output/gene_trends.csv`) are a ready-made readout list. But note the
confound: because pseudotime tracks treatment, the *first* validation
should be whether these genes move in **vehicle** co-culture at all
(i.e., are they drug effects or state effects?).
