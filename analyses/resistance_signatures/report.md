# Public-data KRAS-inhibitor resistance expression signatures

**Run:** 2026-09-22 · `run.py` (re-runnable; reads `data/`, writes `output/`)
· Repo: `~/workspace/pdac`, branch main @ 93d77f57 — nothing committed, nothing pushed.

## Question

Can we derive bulk-RNA expression signatures of *acquired* KRAS-inhibitor
resistance from public GEO data — a PDAC signature from GSE269985 (MRTX1133,
G12D inhibitor, closest to daraxonrasib biology) and an RMC-6236 signature
from GSE324125 (H358 NSCLC, the only deposited RMC-6236-resistant line) — and
does the drug-specific signal transfer across tissue?

## Method

**Data acquisition (processed-first strategy).** Series matrices from GEO FTP
carried metadata only (no data rows) for all three accessions. Supplementary
files were the usable processed source:

- **GSE269985** — `GSE269985_RAW.tar` (13.4 MB): 30 per-sample raw-count CSVs.
  Used directly (no SRA). Included: **human arm** (GPL24676: Panc 02.03 and
  PANC-1, 3 parental + 3 resistant replicates each, n=12) and **mouse arm**
  (GPL24247: KPC lines 6694C2, 6499C4, 6419C5, 3+3 replicates each, n=18).
  Both are MRTX1133 dose-escalation resistant lines — the G12D-inhibitor arm
  the task prioritized. Excluded from this analysis: the clinical adagrasib/
  sotorasib G12C patient data (not deposited as processed matrices).
- **GSE269313** — companion in vivo KPC arm: **infeasible here with reason**.
  Deposited processed data are per-sample snRNAseq `.rds` objects (~90+ MB
  and growing) and this sandbox has no R or pyreadr to parse them. The
  download was stopped after confirming the format.
- **GSE324125** — `GSE324125_STAR_Gene_Counts_*.csv.gz` + `*_TPM_*.csv.gz`
  (~1.1 MB): 6 samples, **all RMC-6236-resistant H358**; 3× DOX− (control) vs
  3× DOX+ (tet-MUC1 shRNA silencing). **No parental H358 was deposited**, so a
  direct parental-vs-resistant DEG is impossible from deposited data. MUC1
  knockdown verified in counts (218–234 → 43–60).

**DEG:** Welch t-test on log2(CPM+1) (GSE269985) or log2(TPM+1) (GSE324125),
Benjamini–Hochberg FDR. Not DESeq2/edgeR — deliberately simple per the task.
DEGs were computed pooled across cell lines (6v6 human, 9v9 mouse), with
per-line log2FC concordance saved separately.

**Signatures:** |log2FC| ≥ 1 and FDR < 0.05, split UP (higher in resistant) /
DOWN. Scorer `score_signatures.py` used **unchanged** (rank method);
`signatures_v1.yaml` **not edited** — proposed additions are in
`output/signatures_v2_proposed_additions.yaml`.

## Results

### 1. PDAC KRASi-resistance signature (GSE269985 human arm, n=12)

- **13 UP + 1 DOWN** genes pass cutoffs (pooled 6v6 across two lines is
  conservative — line effects inflate variance).
- UP: SLITRK6, ZNF608, TLE4, MME, SEMA6B, PDE4B, C11orf86, PLCXD3, DSP, GPCPD1,
  NEO1, TNFRSF1B, FHL1. DOWN: ALPG.
- All 14 are direction-concordant in **both** PANC-1 and Panc 02.03
  (`perline_log2FC_GSE269985_human.csv`). Strongest: MME (+3.1–3.4),
  FHL1 (+3.7), TLE4 (+1.8).
- KRAS/MYC/MET/CDK6 mRNA are **not** significantly changed
  (KRAS log2FC +0.09, padj 1.0) — resistance amplifications in the source
  paper are genomic; bulk RNA doesn't capture them here.

### 2. KPC KRASi-resistance signature (GSE269985 mouse arm, n=18)

- **102 UP + 319 DOWN** genes.
- Top UP: Ifitm3 (+3.7), Sema3a (+3.6), Brsk1 (+2.5), Col7a1 (+2.2), Aim2
  (+2.1), Wnt10a (+1.9), Akt3 (+1.4), Nr4a1 (+1.2). Top DOWN: S100g (−4.1),
  Sptssb (−3.8), plus a block of detox/xenobiotic genes (Ugt2b34/35, Sult1c2,
  Sult1d1, Gstt1, Gsta1, Cyp2c65, Cbr3).
- **CDK6 up +2.0 (padj 0.069)** — just misses FDR but matches the source
  paper's CDK6-amplification theme.
- MAPK feedback genes (Dusp6, Spry2, Egr1) are flat — consistent with MAPK
  pathway reactivation rather than MAPK-off.

### 3. H358 RMC-6236 signature (GSE324125, n=6)

- No parental H358 exists in the deposit, so this is a **MUC1-C-supported
  program within RMC-R cells**, not a resistance-vs-parental signature:
  **93 genes DOWN on MUC1 silencing** (supported by MUC1-C) and **186 genes
  de-repressed** (up on silencing).
- MUC1-supported top genes: MAP2, NPNT, PYGM, HAPLN1, CHI3L1, SLITRK5,
  ALDOC, ADORA1, INSIG1, CA12, F2RL2. De-repressed: MAP2K3, ANKRD1, G0S2,
  SEMA7A, PLK3, LAMA3, ATP13A2.

### 4. Transfer tests — the drug-specific signal does NOT transfer across tissue

| Test | Result |
|---|---|
| T1: H358 MUC1-supported signature → GSE269985 human resistant vs parental | p = 0.82 — **no transfer** |
| T1: H358 MUC1-supported signature → GSE269985 mouse resistant vs parental | p = 0.54 — **no transfer** |
| T2: PDAC-UP signature → H358 RMC-R, DOX− vs DOX+ | means 0.647 vs 0.646, p = 1.0 — **no response to MUC1 silencing** |
| T2: KPC-UP signature → H358 RMC-R, DOX− vs DOX+ | p = 0.70 — **no transfer** |
| Positive controls (each signature on its own data) | human UP p = 0.002; mouse KPC-UP p = 0.0004; H358 own-data direction correct (0.86 vs 0.72, p = 0.10, n=3+3) |

- Gene-level overlap is likewise ~nil: H358 MUC1-supported ∩ KPC-UP human
  orthologs = {MAP2} only; ∩ human PDAC-UP = ∅.
- **Interpretation:** the RMC-6236 resistance program in H358 is lung-lineage
  and MUC1-C-centric; the MRTX1133 PDAC resistance program is a different,
  broader rewiring (xenobiotic-metabolism shutdown, Wnt/axon-guidance genes,
  CDK6). Cross-tissue transfer of a 6-sample-derived signature was always a
  long shot; the honest result is negative.

### 5. Comparison with the Aronchik et al. Nature Medicine 2026 resistance map

(GSE335219 catalog notes: KRAS amp 36%, MAPK 25%, RTK 9%, PI3K 9%, MYC amp,
RTK upregulation, partial EMT, no secondary KRAS mutations — clinical +
HPAC-xenograft daraxonrasib data.)

- **Concordant:** CDK6 up in mouse arm (+2.0, padj 0.069 — cell-cycle
  escape); Akt3 up (+1.4) touches PI3K-axis rewiring; flat MAPK-feedback
  genes fit MAPK reactivation.
- **Discordant/absent:** no KRAS/MYC/MET mRNA upregulation (amps are
  genomic, invisible in bulk RNA); no EMT signature (no VIM/ZEB1/SNAI signal
  in the significant sets); RTK upregulation not evident at RNA level.
- Bulk RNA-seq of resistant *lines* captures downstream rewiring, not the
  genomic amplifications that dominate the clinical map — a real biological
  limit of this modality, not just an analysis artifact.

## Limitations (be honest)

1. **Simple stats, small n.** Welch t-test on log2(CPM/TPM+1), no
   variance shrinkage (no DESeq2/edgeR/limma). n=3 per line-condition;
   GSE324125 is 3v3. Treat all signatures as hypothesis-generating.
2. **Human PDAC signature is tiny** (14 genes) because pooling two distinct
   lines is conservative. It is concordant across lines, but thin.
3. **No true H358 parental-vs-RMC-R comparison exists** in deposited data —
   the "RMC-6236 signature" is really a MUC1-C-dependency program inside
   resistant cells. Any claim about RMC-6236 *resistance* per se from
   GSE324125 overreaches the data.
4. **Cross-species mapping is naive** (mouse symbols upper-cased for human
   orthologs). Works for most genes but will mis-map a few.
5. **Bulk vs single-cell:** these are bulk line profiles; no CAF/stromal
   signal. The in vivo KPC arm (GSE269313, snRNAseq — the one with stroma)
   was infeasible here (RDS objects, no R).
6. **Transfer negatives are underpowered** (n=3+3 in H358), but the
   near-identical T2 means (0.647 vs 0.646) and near-zero gene overlap make
   "no transfer" the right read, not just "couldn't detect."

## What this means for the thesis

- **Usable now:** the KPC mouse-arm signature (102 UP / 319 DOWN, n=18, 3
  lines, good concordance) is the strongest public-data KRASi-resistance
  expression signature in our hands — propose it for v2 of the signature
  compendium (`output/signatures_v2_proposed_additions.yaml`), and score it
  on her CAF/RMC-6236 scRNA-seq to ask whether tumor cells (not CAFs) carry
  it. The GSE337208 finding (MAPK signal absent in fibroblasts in vivo)
  pairs naturally: test the resistance signature in tumor-cell vs CAF
  compartments separately.
- **Negative result worth keeping:** the only deposited RMC-6236-resistant
  line (H358, lung) yields a MUC1-C program that does **not** transfer to
  PDAC models. Don't use H358-derived resistance genes as priors for her
  PDAC work — and cite this negative when justifying PDAC-specific models.
- **Don't over-claim:** bulk line RNA-seq can't see the genomic amps that
  dominate the clinical daraxonrasib resistance map (Aronchik). Pair any
  expression signature with the genomic picture, and treat the 14-gene human
  signature as a seed, not a panel.
- **Still missing:** an in vivo PDAC resistance expression dataset with
  stroma (GSE269313 needs R; GSE336609 has RMC-6236 PDX + TME scRNA-seq but
  is short-term treatment, not acquired resistance). That gap is the honest
  next target, not more squeezing of these three series.

## Files

- `run.py` — full re-runnable pipeline
- `data/` — downloaded GEO files (13 MB GSE269985 tar + series matrices +
  GSE324125 processed CSVs). **Gitignore this.**
- `output/` —
  DEG tables: `deg_GSE269985_human_resistant_vs_parental.csv`,
  `deg_GSE269985_mouse_resistant_vs_parental.csv`,
  `deg_GSE324125_DOXplus_vs_DOXminus.csv`
  signatures: `sig_PDAC_KRASi_resistance_{UP,DOWN}.csv`,
  `sig_KPC_KRASi_resistance_{UP,DOWN}.csv` (with human-ortholog guesses),
  `sig_H358_RMCR_MUC1_{supported,derepressed}.csv`,
  `signatures_v2_proposed_additions.yaml`
  transfer: `scores_T1_*.csv`, `scores_T2_*.csv`, `transfer_T*_groupmeans*.csv`,
  `transfer_mannwhitney.csv`
  per-line concordance: `perline_log2FC_GSE269985_{human,mouse}.csv`
  expression matrices: `expr_GSE269985_{human,mouse}_CPM.csv`,
  `expr_GSE324125_TPM.csv`
  figures: `volcano_GSE269985_{human,mouse}.png`,
  `volcano_GSE324125.png`, `heatmap_GSE269985_{human,mouse}.png`,
  `heatmap_GSE324125.png`
