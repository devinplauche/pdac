# Analysis A — Ligand–receptor search for the indirect-reprogramming hypothesis

## Question

The main GSE337208 mining found that RMC-6236's MAPK pharmacodynamic readout
(Dusp6/Spry/Etv) does **not** move in fibroblasts, suggesting CAF reprogramming
is indirect (via tumor cells / the microenvironment). Which ligands change in
the **non-fibroblast** compartments (tumor/epithelial, immune) under RMC-6236 vs
vehicle — KPC2 primary, MRTX1133 arms secondary — in a way that could drive
the fibroblast iCAF shift? Focus on known iCAF-inducing axes:
IL-1 (Il1a/Il1b), IL-6 family (Il6, Lif, Osm), Cxcl1, Csf2/Csf3, Tnf.

## Method

- Same loading/QC/normalization conventions as
  `../gse337208_caf_mining/run.py`: CellRanger matrices from the already
  downloaded `data/` (no re-download), QC (n_genes ≥ 200, UMI ≥ 500,
  pct_mt ≤ 15), CP10k + log1p, deterministic k-means (seed 0) on the
  16-gene fibroblast/epithelial/immune marker panel.
- Compartment labels from the sibling `cluster_annotation.csv`:
  fibroblast = `fibroblast_cluster`; epithelial = non-fibroblast clusters
  with epi_mean ≥ 1.0; immune = non-fibroblast, epi_mean < 1.0,
  imm_mean ≥ 0.7. (Cell-level compartment masks use the same per-cell
  program comparison the sibling used for doublet cleanup; the fibroblast
  counts were cross-checked against the sibling annotation.)
- Per-sample pseudobulk counts per compartment for KPC2 Veh (GSM9850800) vs
  KPC2 RMC6236 (GSM9850799); secondary: KPC1 Veh (GSM9850794) vs KPC1
  MRTX1133 (GSM9850795). For each of the 9 candidate ligand genes:
  log2 fold change of log2(CP10k + 1) pseudobulk expression, treated vs
  vehicle. **n = 1 sample/arm → descriptive only, no p-values.**
- Cognate receptor expression checked in the KPC2 fibroblast pseudobulk
  (`../gse337208_caf_mining/output/fibroblast_pseudobulk_counts.tsv`),
  as CP10k.
- Ligand → receptor map: Il1a/Il1b → Il1r1 (+ Il1r2 decoy);
  Il6 → Il6ra + Il6st (gp130); Lif → Lifr + Il6st; Osm → Osmr + Il6st;
  Cxcl1 → Cxcr2; Csf2 → Csf2ra + Csf2rb; Csf3 → Csf3r;
  Tnf → Tnfrsf1a + Tnfrsf1b.

Run: `python3 run.py`. Outputs in `output/`.

## Results

### Primary contrast: KPC2 RMC6236 vs Veh — NEGATIVE

**None of the 9 candidate ligands changes in either non-fibroblast
compartment.** All |log2FC| ≤ 0.16 (threshold for "up/down" was |0.5|;
everything classified "flat") — see `output/ligand_receptor_degs.csv` and
`output/fig1_ligand_log2fc.png`:

| compartment | largest |log2FC| |
|---|---|
| epithelial | Il1b 0.08, Lif −0.09 |
| immune | Tnf 0.16, Il1a 0.16, Osm −0.10 |

Additional observations:
- **Il6 is essentially absent from tumor and immune compartments**
  (CP10k 0.00–0.03 in both arms) — in this dataset Il6 is a fibroblast
  (iCAF) product, not a paracrine signal from tumor/immune cells.
- Il1b is moderately expressed by immune cells (7.3 vs 7.05 CP10k) but
  unchanged by drug.
- The epithelial compartment in the RMC6236 arm is small (261 cells vs
  ~560 in vehicle), limiting sensitivity there.

### Secondary contrast: KPC1 MRTX1133 vs Veh — positive control pattern

Here the immune compartment shows a **coordinated dampening** of
inflammatory ligands under MRTX1133: Il1b −2.99, Cxcl1 −0.92, Osm −0.78,
Tnf −0.60, Il1a −0.60 (log2FC), with epithelial Il1b also down (−0.73).
**Caveat:** this is a different model (KPC1/343P, not KPC2) and a different
drug (MRTX1133, KRAS-G12D-selective, not the pan-RAS RMC-6236) — it shows
that RAS inhibition *can* blunt the inflammatory milieu in one model, but
it is not evidence about RMC-6236 in KPC2.

### Receptor landscape in KPC2 fibroblasts (`receptor_expression_fibroblasts.csv`, `fig2_receptor_expression.png`)

| receptor | CP10k (Veh / RMC6236) | verdict |
|---|---|---|
| Il1r1 | 6.1 / 12.4 | **clearly expressed**, rises with drug |
| Il1r2 (decoy) | 0.2 / 0.6 | low |
| Osmr | 1.4 / 1.5 | present |
| Il6st (gp130) | 0.9 / 1.4 | present |
| Tnfrsf1a / Tnfrsf1b | 1.6–1.8 / 0.5–1.0 | present |
| Il6ra | 0.15 / 0.27 | **effectively absent** |
| Lifr | 0.18 / 0.19 | low |
| Cxcr2 | 0.02 / 0.36 | **effectively absent** |
| Csf2ra / Csf2rb / Csf3r | 0.06–0.65 | low |

So: fibroblasts are competent to receive IL-1 (Il1r1 high), OSM (Osmr +
gp130), and TNF signals, but the classic IL-6 axis is incomplete
(gp130 present, IL-6Rα absent — only trans-signaling via soluble IL-6R
would work), and Cxcl1→Cxcr2 signaling into fibroblasts is essentially
impossible at these expression levels.

## Limitations

- **n = 1 sample/arm.** All fold changes are descriptive; no statistical
  test is possible and single-sample idiosyncrasies cannot be ruled out.
- **Single post-treatment timepoint.** A transient ligand pulse (hours)
  would be invisible in this snapshot.
- **scRNA dropout.** Low-abundance cytokines (Il6-family, Csf2/3) may be
  expressed below detection; "absent" means "not detected," not "not there."
- Compartments are coarse k-means clusters; rare ligand-producing subsets
  (e.g. a specific myeloid population) are diluted in the immune pseudobulk.
- The KPC2 RMC6236 epithelial compartment is small (261 cells).

## What this means for the thesis

This is a **negative result, and negatives count**: in the KPC2 model, the
RMC-6236-driven fibroblast iCAF shift is **not** accompanied by increased
production of any canonical iCAF-inducing ligand by tumor or immune cells
at this timepoint. Together with the absent MAPK PD readout in fibroblasts,
it pushes the indirect-reprogramming hypothesis toward mechanisms other
than "more IL-1/IL-6/TNF from neighbors":

1. **Loss of a myCAF-maintaining signal** (e.g. TGFβ withdrawal) rather
   than gain of an iCAF-inducing one — consistent with the myCAF score
   drop without a MAPK move.
2. **A signal from an unmeasured source** (endothelial cells, nerves,
   or post-transcriptional/protein-level regulation invisible to scRNA-seq).
3. **Fibroblast-autonomous reprogramming** downstream of drug-altered
   tumor–stroma contact or mechanics.

For the wet-lab plan: the IL-1 axis remains the most testable candidate
mechanistically (Il1r1 is the best-expressed ligand receptor on these
fibroblasts, and it rises under drug) — but this dataset gives **no
evidence that RMC-6236 increases IL-1 ligand supply**. A co-culture
experiment should therefore test IL-1R blockade as a *permissive* rather
than *drug-driven* hypothesis, and should include a TGFβ-axis arm given
the myCAF-loss pattern. The KPC1 MRTX1133 immune-dampening pattern is a
reminder that the inflammatory milieu *is* drug-responsive in at least one
model — worth a direct KPC2 time-course (bulk or flow) rather than a
single scRNA snapshot.
