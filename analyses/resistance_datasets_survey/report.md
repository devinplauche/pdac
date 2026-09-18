# RMC-6236 / daraxonrasib resistance datasets survey

**Date:** 2026-09-18
**Question:** Is there a public omics dataset of daraxonrasib-resistant PDAC that can serve as the tumor-cell comparator for the resistance chapter — or must she generate resistant Capan-1/KP-4 lines herself?
**Sources searched:** GEO (NCBI E-utilities, live), ArrayExpress/EBI Search (live), ENA portal API (live). Search script: `run.py` → `output/raw_hits.csv`. Curated table: `output/catalog.csv`.

## Bottom line

**No public omics of daraxonrasib-resistant PDAC models exists.** The flagship resistance paper (Aronchik et al., *Nature Medicine*, 11 Aug 2026; GEO GSE335219) deposited only **2 WES samples** — parental HPAC (human PDAC, KRAS G12D/WT) ± daraxonrasib, as MAF files. The 44-patient paired pre-/end-of-treatment ctDNA targeted-seq (800+ genes) and the human + murine preclinical resistant-model omics described in the paper are **not in the public deposit**. ArrayExpress and ENA return zero usable RMC-6236/daraxonrasib studies.

That is itself the answer to the comparator question: generating her own resistant lines is justified, because the field's central resistance dataset is not publicly reusable. The Nat Med paper is the citation backbone, not a data source.

## What's actually public (ranked by usefulness)

### True RMC-6236/daraxonrasib data
| Accession | What it is | Catch |
|---|---|---|
| **GSE337208** | scRNA-seq of murine KPC PDAC treated with **RMC-6236** or MRTX1133 ± anti-CTLA4/PD1 (9 samples, Jul 2026) | Short-term treatment (d10–14), not acquired resistance — but it's **in vivo PDAC with stromal/CAF populations under RMC-6236**. Best dataset for "what does RMC-6236 do to CAFs in a tumor". |
| **GSE336609** | Appendiceal adenocarcinoma organoids + PDX treated with RMC-6236 (G12V) / MRTX1133 (G12D); TME scRNA-seq with human/mouse reads separated | Wrong primary (appendiceal), raw data withheld (privacy) — but shows a **normal-to-iCAF shift after KRAS inhibition**, directly relevant to her CAF thesis. |
| **GSE324125** | H358 NSCLC **RMC-6236-resistant** (RMC-R) cells, bulk RNA-seq ± MUC1 shRNA (PMID 42415848) | Only true RMC-6236-resistant line with omics found — wrong tissue (lung). Usable as a RAS(ON)-resistance expression signature. |
| **GSE335219** | HPAC PDAC ± daraxonrasib WES (MAF); BioProject PRJNA1477593 | Parental only; the clinical resistance data behind the paper isn't deposited. |

### Closest PDAC substitutes (KRAS-inhibitor resistance, public, raw in SRA)
| Accession | What it is |
|---|---|
| **GSE269985 / GSE269313** (PMID 38975874) | **Best PDAC comparator.** Multi-model KRASi resistance: PDAC lines + organoids under MRTX1133, KPC mouse with emergent resistance, clinical KRAS G12C PDAC (adagrasib/sotorasib). KRAS/MYC/MET/EGFR/CDK6 amplifications, EMT, PI3K-AKT-mTOR — mirrors the daraxonrasib resistance themes. |
| **GSE298318** (PMID 42147356) | Acquired resistance to KRAS G12D inhibition in PDAC. |
| **GSE294216** | MRTX1133 resistance in PDAC via histone acetylation; BET-inhibitor vulnerability. |
| **GSE298925 / GSE298803** (PMID 40789946) | K-Ras-inhibitor-resistant PDAC, metabolic vulnerabilities (ZBTB11). |
| **GSE273116** | Pan-KRAS inhibitor MCB-294 (drug-class proxy for acute RMC-6236 response). |

Full details per dataset (publication, model, treatment, data type, resistant model yes/no, raw availability) are in `output/catalog.csv`.

## Recommended comparator strategy

1. **Now, zero bench cost:** use **GSE269985** (+ GSE298318) as the public PDAC resistance comparator — it's the closest biology with raw data, and its resistance themes (RTK amps, MYC, EMT) match the Nat Med daraxonrasib findings.
2. **For the CAF angle:** mine **GSE337208** (RMC-6236-treated KPC scRNA-seq) for fibroblast/iCAF shifts, and cite **GSE336609**'s normal→iCAF shift as precedent.
3. **Generate the resistant lines anyway:** no public daraxonrasib-resistant PDAC omics exists, so Capan-1/KP-4-derived resistant lines fill a genuine gap — and the Nat Med paper itself generated human + murine resistant models to corroborate the clinic, so her plan mirrors the field.
4. **Open check:** the Nat Med paper's data-availability statement may point to controlled access (dbGaP/EGA) for the ctDNA — worth one look before assuming it's fully closed.

## Reproducibility
`run.py` re-runs the live searches (GEO esearch/esummary over 8 query buckets, ENA portal, EBI ArrayExpress) and rewrites `output/raw_hits.csv`. `output/catalog.csv` is the human-curated layer on top — verify before reuse, especially publication links.

## Verification addendum (2026-09-18)

Follow-up verification pass resolving the open items from the original survey. No commits made; catalog.csv updated in place.

### 1. Papers behind the newest deposits — all three VERIFIED
- **GSE337208 →** Mahadevan KK et al., "Oncogenic Kras targeting with MRTX1133 or Daraxonrasib specifically synergize with anti-CTLA4 to promote anti-tumor immunity in pancreatic cancer," *Nature Communications* 2026;17(1):9105 (PMID 42649211, DOI 10.1038/s41467-026-75960-3). The paper's data-availability statement explicitly cites GSE337208. Match to the survey's description: confirmed — the survey's "Mahadevan et al. 2026 (KRAS + anti-CTLA4 synergy in PDAC)" attribution was correct.
- **GSE336609 →** Chowdhury S, Ito I, et al. (corresponding John Paul Shen), "KRAS inhibition is an effective therapy for appendiceal adenocarcinoma," *Journal of Hematology & Oncology* 2026;19(1):70 (PMID 42393724, PMCID PMC13536544, DOI 10.1186/s13045-026-01817-3; preprint bioRxiv 2026.04.07.717107). The GEO series title is verbatim the paper title. Match: confirmed — AA organoids + PDX ± MRTX1133/RMC-6236, human/mouse reads separated, normal-fibroblast → iCAF shift reported.
- **GSE294216 →** Principe DR, Becker JH, et al. (senior Hidayatullah G. Munshi), "Resistance to the KRASG12D Inhibitor MRTX1133 is Associated with Increased Sensitivity to BET Inhibition," *Molecular Cancer Therapeutics* 2026;25(6):1001 (PMID 41527355; DOI not retrieved this pass). Match: confirmed — acquired MRTX1133 resistance in human PDAC with histone-acetylation shift and BET-inhibitor re-sensitization.

### 2. Aronchik et al. Nature Medicine — data-availability statement
Paper: Aronchik I, Kar S, et al., *Nat Med* 2026;32(8):2865–2877 (PMID 42581230, PMCID PMC13472880 — open access; DOI 10.1038/s41591-026-04537-w). Data Availability statement (verbatim): "All data supporting the findings of this study are available within the article and its Supplementary Information. DNA sequencing data from the parental and resistant HPAC xenograft models have been deposited in the NCBI Gene Expression Omnibus under accession number GSE335219. Clinical and translational data were generated as part of the phase 1/2 study NCT05379985 and available as stated in 10.1056/NEJMoa2505783. Source data are provided with this paper."
- **No controlled access anywhere in the statement** — zero mentions of dbGaP, EGA, or any controlled-access repository. The only deposit is public GEO (GSE335219).
- The 44-patient paired ctDNA targeted-seq is **not deposited**; availability is deferred to the phase 1/2 trial paper (Wolpin et al., *NEJM* 2026;394:1790–1802, PMID 42090791), whose data-sharing statement could not be verified here (subscription-only, not in PMC) — treat as **UNVERIFIED**, likely sponsor-mediated request access typical of industry trials.

### 3. GSE335219 GEO record — claim VERIFIED, with refinements
Confirmed via E-utilities + GEO FTP listing: n_samples = 2 (GSM9808344 "HPAC DMSO", GSM9808345 "HPAC Daraxonrasib"), MAF files in `GSE335219_RAW.tar`, BioProject PRJNA1477593. Refinements vs the original survey text:
- The samples are **xenograft-model** data ("parental and resistant HPAC xenograft models"), not simply cell line ± drug; filenames indicate parental passage P18 vs resistant passage P28 selected in 0.1 µM daraxonrasib.
- **Raw reads also exist in SRA** (SRX33858264 HPAC DMSO; SRX33858265 HPAC Daraxonrasib) — the survey had noted only the MAFs.
- No human+murine preclinical omics in this deposit — consistent with the paper's statement naming only the HPAC xenograft deposit.

**Bottom line of this addendum:** every open item resolved except the NEJM trial paper's data-sharing statement (UNVERIFIED). The core verdict stands unchanged — no public daraxonrasib-resistant PDAC omics, nothing under controlled access, so generating her own resistant lines remains justified.
