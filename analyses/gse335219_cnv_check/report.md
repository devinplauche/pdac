# Verifying Aronchik et al. Nature Medicine 2026 resistance claims from GSE335219

## Question
Can the resistance-genomics map from Aronchik et al. (*Nat Med* 2026;32(8):2865–2877;
mutant KRAS amp 36%, MAPK reactivation 25%, RTK 9%, PI3K 9%, MYC amp,
RTK upregulation, partial EMT, no secondary KRAS mutations) be verified from the
deposited primary data? The deposited data for GSE335219 are HPAC xenograft WES
MAFs only: parental P18 (GSM9808344, DMSO) vs resistant P28 (GSM9808345, selected
at 0.1 µM daraxonrasib). **The paper's stated frequencies come from ctDNA of 44
phase 1/2 patients that was never deposited — this is an n=1 pair and cannot
reproduce 36%/25%/9%/9% frequencies.** The verifiable questions are the
mechanistic ones: secondary KRAS mutations? MAPK/PI3K pathway mutations?

## Method
- Downloaded `GSE335219_RAW.tar` from GEO (68,229,120 bytes; MD5-size matched
  Content-Length; note: plain curl over the GEO FTP/HTTPS endpoint kept dying
  near 92% with `curl: (18) transfer closed`; `wget -c` over HTTPS completed it —
  see `run.sh`).
- Tar contains exactly two files, both MAF 2.4: `..._Parental-P18-..._tnscope.vcf.maf.gz`
  (372,090 rows) and `..._0-1uM-...-P28-..._tnscope.vcf.maf.gz` (378,624 rows).
  Sentieon TNscope, tumor-only with panel-of-normals (Matched_Norm = "NORMAL"
  placeholder; no real matched normal).
- **No copy-number segment data deposited.** Per task instructions, no copy-number
  calls were fabricated from MAFs; the analysis is mutation-level only.
- `run.py` filters to PASS filter, protein-coding consequence, gnomAD AF < 1% or
  blank (2,158 parental / 2,123 resistant variants), then keys variants on
  (chr, start, end, ref, alt) and splits shared / parental-private / resistant-private.
- Pathway panels (curated): MAPK {KRAS/NRAS/HRAS/BRAF/RAF1/ARAF/MAP2K1/MAP2K2/
  MAPK1/MAPK3/NF1/SOS1/SHOC2/PTPN11/RIT1/MRAS}, PI3K {PIK3CA/B/D/G/R1/R2/PTEN/
  AKT1/2/3/MTOR/TSC1/TSC2/RHEB}, RTK {EGFR/ERBB2/3/4/MET/FGFR1–4/PDGFRA/B/KIT/
  RET/ALK/ROS1/NTRK1–3/IGF1R/KDR/FLT1/FLT4/TEK/AXL/DDR2}, MYC {MYC/MYCN/MYCL}.

## Results
- Overlap: 1,600 shared / 558 parental-private / 523 resistant-private
  high-confidence coding variants (see `output/variant_comparison_summary.csv`,
  `output/private_mutations.csv`, `output/vaf_shared_vs_private.png`).
- **KRAS: no secondary KRAS mutations.** Both samples carry exactly one KRAS
  variant — shared KRAS p.G12D (chr12:25245350 C>T, c.35G>A). Zero
  resistant-private KRAS variants (see `output/kras_mutations.csv`).
  → Consistent with the paper's "no secondary KRAS mutations."
- **KRAS VAF rises parental→resistant: 0.49 (153/312 reads) → 0.82 (749/913).**
  This is *suggestive* of mutant-allele amplification in the resistant line but is
  NOT a copy-number call — noted as an observation consistent with, not a
  verification of, the paper's mutant-KRAS-amplification finding.
- **MAPK reactivation: resistant-private MAP2K1 p.K57E** (chr15:66435115 A>G,
  VAF 0.40, depth 146; activating MEK1 allele). Absent from the parental
  callset entirely (checked all FILTER states, not just PASS).
  → Mutation-level corroboration of the paper's MAPK-reactivation mechanism.
- **PI3K activation: resistant-private PIK3CA p.E542K** (chr3:179218294 G>A,
  VAF 0.37, depth 152; canonical hotspot). Absent from parental.
  → Mutation-level corroboration of the paper's PI3K mechanism.
- **RTK / MYC: zero resistant-private protein-coding variants** in the RTK or MYC
  panels at the point-mutation level — consistent with the paper framing those
  as copy-number/expression-level mechanisms (untestable here).
- Also resistant-private: GNAS p.R232H (VAF 0.35, depth 270; not a classic
  hotspot — treat cautiously); ARID1A and KMT2C low-VAF (<0.04) subclonal hits.
- Parental-private pathway hit: SHOC2 p.N550K at VAF 0.075 (subclonal loss).
- The bulk of private calls are dominated by MUC-family / HLA / KIR genes and
  pile at VAF ~0.05 (see figure) — the private sets are noisy; only the
  above high-VAF, hotspot, well-covered calls are interpreted.

## Limitations
- **Copy-number verification infeasible from deposited data.** GEO holds only
  MAFs; there are no segment files. Frequencies (36%/25%/9%/9%) come from
  undeposited 44-patient ctDNA and cannot be reproduced by an n=1 pair.
- Tumor-only WES callset (panel-of-normals, no matched normal): private variant
  sets contain mapping artifacts (MUC/HLA/KIR pileup); ~558 vs 523 private calls
  are mostly noise. Only hotspot, high-VAF calls were interpreted.
- Private = absent from parental *callset*; a low-frequency parental subclone
  below detection cannot be excluded.
- **WES re-analysis assessed and not pursued.** Raw reads exist in SRA
  (BioProject PRJNA1477593; SRX33858264 HPAC-DMSO = SRR39113092, 66.4M spots /
  19.9 Gb; SRX33858265 HPAC-daraxonrasib = SRR39113091, 57.6M spots / 17.3 Gb;
  ~5.9 GB each, ~11 GB total — under the 15 GB bar). Judged not worthwhile
  within the timebox: full re-alignment (no aligner/CNVkit installed, 2 vCPU),
  no matched normal, and human-xenograft mouse-contamination risk — for an n=1
  pair that still could not reproduce the clinical frequencies. Recorded here
  instead of attempted.

## What this means for the thesis
The deposited preclinical data corroborate three of Aronchik's mechanistic
claims at the mutation level in the HPAC resistant model: (1) no secondary KRAS
mutations; (2) acquired MAPK-pathway reactivation (MAP2K1 K57E); (3) acquired
PI3K-pathway activation (PIK3CA E542K). The headline KRAS-amplification
mechanism cannot be tested from the deposited MAFs (though the KRAS G12D VAF
0.49→0.82 shift is directionally consistent); RTK/MYC claims are
copy-number/expression mechanisms invisible to this data. For the thesis, these
two acquired point mutations are usable as concrete resistance-alteration
examples alongside the paper's ctDNA map — cite as n=1 preclinical corroboration,
not frequency evidence.
