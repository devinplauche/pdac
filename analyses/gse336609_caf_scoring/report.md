# Cross-tissue validation of KRAS-inhibitor CAF reprogramming — GSE336609

## Question
Does KRAS inhibition drive the same myCAF→iCAF stromal shift in an
independent tissue/model as RMC-6236 did in PDAC KPC fibroblasts? Our
GSE337208 mining found, at single-cell resolution in KPC2 fibroblasts:
myCAF down (MWU p=0.012), iCAF up (p=0.0005), apCAF up (p=0.0002).
GSE336609 (Shen et al., *J Hematol Oncol* 2026) tests KRAS inhibitors in
human appendiceal adenocarcinoma PDX — a different tissue, different disease —
with human (tumor) and mouse (host stroma) reads extracted separately.

## Method
- **Data actually deposited (verified):** raw data / scRNA-seq were **not**
  deposited — the GEO record states "raw data are missing due to patient
  privacy concerns." The only files are four per-sample **bulk** count
  matrices (treated vs vehicle contrast files):
  AAP01 (KRAS G12D) ± MRTX1133 (3+3 mouse, 3+3 human) and AAP16 (KRAS G12V)
  ± RMC-6236 (2+3 mouse, 2+3 human). There is therefore **no fibroblast-level
  data** — this analysis scores the bulk stromal (mouse-read) signal, exactly
  the fallback the task specified.
- **Gene mapping:** matrices are Ensembl-ID rows. The ~150 signature/marker
  symbols were resolved to matrix Ensembl IDs via mygene.info (mouse) and the
  HGNC REST API (human); non-target rows keep Ensembl IDs (unique, so the
  rank distribution used by the scorer is unchanged). Coverage: **100% of
  signature genes found** in both matrices (myCAF 6/6, iCAF 6/6, apCAF 5/5,
  SASP 14/14, matrisome 19/19, MAPK 9/9, RASON 7/7).
- **Scoring:** repo `score_signatures.py` + `signatures_v1.yaml` **as-is**,
  `--organism mouse` (stroma) / `human` (tumor PD check), rank method
  (AUCell-like, per-sample) as primary, zscore as secondary.
- **Contrasts:** treated vs vehicle **per model** (AAP01: MRTX1133 n=3v3;
  AAP16: RMC-6236 n=2v3). Mann–Whitney U on per-sample rank scores
  (**sample-level**, not cell-level — no cells exist here; minimum achievable
  two-sided p is 0.10 for 3v3 and 0.20 for 2v3, so p<0.05 is unreachable by
  construction) plus member-gene coherence: median log2FC of signature genes
  and Wilcoxon signed-rank on member-gene log2FCs (gene-level, descriptive).
- **Compartment check:** fibroblast / myofibroblast / pericyte / immune /
  endothelial / host-epithelial marker fractions of mouse reads, to ground
  whether the mouse signal is fibroblast-dominated.
- Re-runnable: `python3 run.py` (downloads GEO files if missing; gene maps
  cached in `data/`).

## Results

### Compartment composition (mouse reads)
Among identifiable compartments the mouse signal is fibroblast-dominated:
fibroblast markers 1.3–3.6% of mouse reads vs immune 0.13–1.02%,
endothelial ≤0.06%, host-epithelial ≤0.12% (absolute %s are small because most
bulk reads are housekeeping/ribosomal, as expected). So "stromal signal"
reads as primarily fibroblast — with the infiltration caveat below.

### AAP01 stroma — MRTX1133 (G12D-selective), 3v3
| Signature | Δ score | MWU p | median member log2FC | signed-rank p |
|---|---|---|---|---|
| myCAF | −0.024 | 0.10 (= min) | −1.82 | **0.031** |
| iCAF | +0.017 | 0.40 | −0.02 | 1.00 |
| apCAF | +0.010 | 0.10 | +0.42 | 0.0625 |
| SASP | +0.038 | 0.70 | +0.22 | 0.43 |
| matrisome | −0.005 | 0.10 | −0.52 | **0.0039** |
| stromal MAPK | −0.017 | 0.70 | −0.33 | 0.039 |

- **myCAF collapses coherently:** all 6 member genes down — Acta2 −2.92,
  Tagln −2.66, Col12a1 −2.42, Igfbp3 −1.21, Thbs2 −1.14, Thy1 −0.80
  (signed-rank p=0.031). Pan-fibroblast markers are flat (Col1a1 −0.37,
  Dcn +0.21, Pdpn +0.31, Vim −0.30): the contractile **program** is down,
  fibroblasts are not depleted.
- **iCAF does NOT induce:** Il6 (+1.10) and Cxcl1 (+1.74) rise, but Clec3b,
  Has1, Ly6c1 are flat and Col14a1 falls (−1.43) — no coherent program.
- **apCAF "up" is Saa3/Slpi, not MHC-II:** Saa3 +2.91, Slpi +1.23, but
  H2-Aa/H2-Ab1/Cd74 barely move (+0.04 to +0.42). Acute-phase signal, not
  antigen presentation.
- Immune infiltration essentially unchanged (Ptprc +0.25).
- Tumor PD (human reads) is clean: MAPK_targets Δ −0.055, member −0.94,
  p=0.0039; RASON_inhibition_PD Δ −0.066, member −0.94, p=0.016 — all 7 genes
  down (FOSL1 −2.78, SPRY4 −1.41, CCND1 −1.08, ETV4 −0.94, DUSP6 −0.49).
  On-target KRAS inhibition confirmed in tumor.

### AAP16 stroma — RMC-6236 (pan-KRAS), 2v3
| Signature | Δ score | MWU p | median member log2FC | signed-rank p |
|---|---|---|---|---|
| myCAF | +0.040 | 0.20 (= min) | +0.97 | 0.44 |
| iCAF | +0.065 | 0.20 | +2.19 | 0.094 |
| apCAF | +0.001 | 1.00 | −0.67 | 0.81 |
| SASP | +0.104 | 0.20 | +1.15 | **0.017** |
| matrisome | +0.016 | 0.20 | +0.01 | 0.52 |
| stromal MAPK | −0.036 | 0.20 | −1.16 | **0.0039** |

- **myCAF does NOT suppress** — members split: Tagln +2.79, Col12a1 +2.28,
  Thy1 +1.96 up while Acta2 is flat (−0.01) and Igfbp3/Thbs2 fall. No coherent
  change (opposite of AAP01 and of PDAC).
- **Inflammatory program induces:** iCAF 5/6 genes up (Il6 +3.53, Cxcl1
  +3.23, Ly6c1 +2.56, Clec3b +1.81; Col14a1 −1.36), SASP coherent up
  (p=0.017). **But immune infiltration rises in parallel** (Ptprc +1.28,
  Itgam +1.65, Cd3e +0.72), and Ly6c1 is also a monocyte marker — in bulk the
  inflammatory signal cannot be assigned to fibroblasts vs myeloid cells.
  (The 9–12x cytokine induction exceeds the ~2.4x Ptprc rise, so a genuine
  inflammatory induction exists; its cellular source is ambiguous.)
- **apCAF does not replicate:** MHC-II genes go *down* (H2-Ab1 −1.03,
  H2-Aa −0.67, Cd74 −0.73) while Saa3 (+4.64)/Slpi (+3.64) spike — again
  acute-phase, not antigen presentation.
- **Stromal MAPK is directly suppressed** (member −1.16, p=0.0039) — expected:
  RMC-6236 is pan-KRAS and stroma is WT KRAS, so the drug hits stromal KRAS
  directly. (In AAP01 the stromal MAPK dip is marginal, p=0.039, consistent
  with an indirect effect of a G12D-selective drug.)
- Tumor PD (human reads) is **uninterpretable**: human libraries are
  210–454K reads / 6.5–7.8K genes (vs 1.3–15M / 17–26K in AAP01) — the tumor
  was essentially ablated ("dramatically reduced tumor cellularity" per the
  paper) and most MAPK genes are near-zero counts. This is a sparsity
  artifact, not a negative PD result.

### Figures
- `output/fig1_stromal_caf_scores.png` — per-sample rank scores, all 7
  signatures, both models (dots = samples, tick = median).
- `output/fig2_member_gene_log2fc.png` — per-gene log2FCs for myCAF / iCAF /
  apCAF members, both models.

## Limitations (read before citing)
1. **Bulk, not single-cell.** The deposited data cannot resolve fibroblasts
   from other host cells. The GSE337208 result was fibroblast-resolved; this
   is a stromal-bulk approximation — a weaker test by design.
2. **Tiny n.** 3v3 and 2v3: per-model MWU cannot reach p<0.05. All
   "significance" leans on gene-level coherence (signed-rank), where genes
   are not independent observations.
3. **Drug confound.** MRTX1133 (G12D-selective) vs RMC-6236 (pan-KRAS) differ
   in whether they can directly engage stromal WT KRAS — we *observe* direct
   stromal MAPK suppression only with RMC-6236. Model/drug/tissue are fully
   confounded (one drug per model).
4. **Infiltration confound (AAP16).** The iCAF-like signal co-occurs with
   increased Ptprc/Itgam; Ly6c1 marks monocytes too. Fibroblast reprogramming
   vs myeloid influx cannot be separated in bulk.
5. **apCAF non-replication is about signature composition:** Saa3/Slpi spike
   in both models while MHC-II genes don't — the "apCAF up" in PDAC was
   MHC-II-driven and is not reproduced here.
6. **Appendiceal ≠ PDAC.** Even a positive result would be cross-tissue
   corroboration, not PDAC validation. AAP16's tumor compartment is too
   ablated to assess tumor PD from these matrices.

## What this means for the thesis
- **No clean cross-tissue replication of the full myCAF→iCAF→apCAF program.**
  Do not cite GSE336609 as independent validation of the three-way shift.
- **The myCAF-suppression leg DOES replicate** under MRTX1133 in AAP01:
  contractile genes (Acta2/Tagln) collapse 4–7x with fibroblast markers
  intact — and it happens without direct stromal MAPK engagement
  (G12D-selective drug, WT-KRAS stroma), implying a **tumor-mediated**
  mechanism. This is the strongest cross-model corroboration and the
  cleanest conceptual match to the thesis: KRAS inhibition in tumor cells
  relaxes the myofibroblast program in neighboring stroma.
- **The iCAF-induction leg is model-dependent and confounded:** present as a
  bulk inflammatory program in AAP16 (RMC-6236) but inseparable from myeloid
  infiltration; absent as a coherent program in AAP01. The thesis's
  single-cell fibroblast-resolved data remains the authoritative evidence
  for iCAF induction — this dataset can't strengthen or weaken it.
- **Do not claim apCAF induction** from this dataset (Saa3/Slpi ≠ MHC-II).
- **Design lesson for the thesis wet-lab work:** the G12D-selective arm
  (AAP01) is the informative control — it isolates tumor-mediated stromal
  effects from direct stromal KRAS inhibition. Any follow-up co-culture
  experiment should include a G12D-selective inhibitor arm alongside
  RMC-6236 to separate these two mechanisms, exactly as this dataset
  accidentally does.
- The paper's own scRNA-seq claim ("shift from normal to inflammatory CAFs")
  is neither confirmed nor refuted here — those data were not deposited.

## Provenance
- GEO GSE336609, public 2026-06-25; series record confirms "raw data are
  missing due to patient privacy concerns."
- Panel: `../signature_panel/signatures_v1.yaml` (v1.0), scorer
  `../signature_panel/score_signatures.py` used as-is (rank method primary).
- Symbol resolution: mygene.info (mouse, 83/83) + HGNC REST (human, 56/56);
  cached in `data/symbol_to_ensembl_{mouse,human}.tsv`.
- Outputs: `output/` — `contrast_stats.csv`, `member_gene_log2fc.csv`,
  `mouse_scores_per_sample.csv`, `human_scores_per_sample.csv`,
  `mouse_gene_coverage.csv`, `human_gene_coverage.csv`,
  `mouse_compartment_fractions.csv`, `qc_table.csv`, 2 PNG figures.
  Large inputs/intermediates live in `data/` (gitignore).
