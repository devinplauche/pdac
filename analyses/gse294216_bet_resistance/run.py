#!/usr/bin/env python3
"""GSE294216 (Principe/Becker et al. (Munshi), Mol Cancer Ther 2026):
transcriptional program of BET-inhibitor-vulnerable MRTX1133-resistant PDAC.

Data: GEO processed featureCounts matrices (no alignment; processed only).
- GSE294216_PANC1_PANC1K_countDatafeatureCounts.txt.gz:
    human PANC1 (parental, n=3) vs PANC1K (MRTX1133-resistant, n=3)
- GSE294216_PANC1k_PANC1KJQ1_countDatafeatureCounts.txt.gz:
    PANC1K (n=3) vs PANC1K + JQ1 (BET inhibitor, n=3)
- GSE294216_PANC1k_PANC1KOTX_countDatafeatureCounts.txt.gz:
    PANC1K (n=3) vs PANC1K + OTX-015 (BET inhibitor, n=3)

Pipeline:
  1. Load matrices (featureCounts-normalized counts; fractional values).
  2. DEG PANC1K vs PANC1 on log2(x+1): Welch two-sample t-test per gene,
     BH-FDR. (Mann-Whitney is useless at n=3/group: two-sided min p = 0.1,
     so it is documented but not used for calls.)
  3. Targeted tests: BRD readers, HAT writers, HDAC/SIRT erasers,
     FOSL1/AP-1 program (paper's pro-survival axis).
  4. Gene-set overlap: hypergeometric test of DEG-up/DEG-down sets vs
     inline Hallmark-style sets (EMT, MYC targets, RTK signaling,
     MAPK targets, apoptosis/BCL2) over the detected-gene background.
  5. Reversal check: do BET inhibitors (JQ1, OTX-015) invert the
     resistance DEG program in PANC1K? t-test PANC1K+treatment vs PANC1K;
     concordance test: fraction of resistance DEGs whose treatment log2FC
     has opposite sign (binomial test vs 0.5).

Re-runnable: python3 run.py  (expects data/*.txt.gz)
Outputs: output/*.csv, output/*.png
"""
import gzip
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind, hypergeom, binomtest

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)

F_MAIN = os.path.join(DATA, "GSE294216_PANC1_PANC1K_countDatafeatureCounts.txt.gz")
F_JQ1 = os.path.join(DATA, "GSE294216_PANC1k_PANC1KJQ1_countDatafeatureCounts.txt.gz")
F_OTX = os.path.join(DATA, "GSE294216_PANC1k_PANC1KOTX_countDatafeatureCounts.txt.gz")
F_SGC = os.path.join(DATA, "GSE294216_PANC1k_PANC1KSGC_countDatafeatureCounts.txt.gz")


# ---------- figures ----------


def load_matrix(path):
    """Return counts DataFrame (genes x samples) and column labels."""
    with gzip.open(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t")
    df = df.rename(columns={"Gene": "gene"})
    df = df.set_index("gene")
    df.columns = [c.strip('"') for c in df.columns]
    return df


def assign_groups_main(df):
    groups = []
    for c in df.columns:
        if "PANC1K." in c or "Panc1K." in c:
            groups.append("resistant")
        elif "Panc1." in c:
            groups.append("parental")
        else:
            groups.append("unknown")
    return np.array(groups)


def assign_groups_bet(df):
    groups = []
    for c in df.columns:
        lc = c.lower()
        if "jq1" in lc or "otx" in lc or "otx-015" in lc:
            groups.append("treated")
        else:
            groups.append("resistant_base")
    return np.array(groups)


def deg_table(df, groups, g1, g2, label):
    """Welch t-test on log2(x+1) normalized counts; BH-FDR."""
    m = np.log2(df.values + 1.0)
    a = m[:, groups == g1]
    b = m[:, groups == g2]
    lfc = a.mean(axis=1) - b.mean(axis=1)  # log2(g1/g2)
    t, p = ttest_ind(a, b, axis=1, equal_var=False)
    p = np.where(np.isnan(p), 1.0, p)
    res = pd.DataFrame({"gene": df.index, "log2FC": lfc, "p": p})
    # Benjamini-Hochberg
    order = np.argsort(p)
    n = len(p)
    fdr = np.empty(n)
    fdr[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    res["fdr"] = np.clip(fdr, 0, 1)
    res = res.sort_values("p").reset_index(drop=True)
    res.to_csv(os.path.join(OUT, f"{label}_deg.csv"), index=False)
    return res


# ---------------- inline gene sets ----------------
EMT = ["SNAI1", "SNAI2", "TWIST1", "TWIST2", "ZEB1", "ZEB2", "VIM", "CDH2",
       "FN1", "MMP2", "MMP9", "MMP14", "COL1A1", "COL3A1", "COL5A1", "SPARC",
       "TGFBI", "SERPINE1", "LOX", "LOXL2", "FAP", "PDGFRB", "ACTA2"]
MYC_TARGETS = ["MYC", "CDK4", "CCND2", "CCNE1", "ODC1", "NCL", "NPM1", "CAD",
               "DHFR", "HSPD1", "HSPE1", "MRPL12", "RPLP0", "RPS7", "EEF1B2",
               "EIF4A1", "EIF4E", "LDHA", "GAPDH", "TPI1", "PGK1", "ENO1",
               "PKM", "SORD", "SNRPA1", "HNRNPA1", "SRSF7", "PSMA2", "PSMB1",
               "PA2G4", "TFAM", "ATP5A1", "FBL", "DKC1", "NOP16", "NOP56",
               "BYSL", "EBNA1BP2", "MRPL3", "RFC4", "MCM4", "MCM5", "MCM6",
               "PCNA", "MCM7", "CDC20", "MAD2L1", "BUB1", "AURKB", "PLK1"]
RTK = ["EGFR", "ERBB2", "ERBB3", "ERBB4", "MET", "PDGFRA", "PDGFRB", "FGFR1",
       "FGFR2", "FGFR3", "FGFR4", "IGF1R", "INSR", "RET", "AXL", "MERTK",
       "TYRO3", "KIT", "FLT1", "FLT4", "KDR", "TEK", "NTRK1", "NTRK2",
       "NTRK3", "ROS1", "ALK", "EPHA2", "EPHB2", "DDR1", "DDR2"]
MAPK = ["KRAS", "HRAS", "NRAS", "BRAF", "ARAF", "RAF1", "MAP2K1", "MAP2K2",
        "MAPK1", "MAPK3", "MAPK14", "DUSP4", "DUSP5", "DUSP6", "SPRY2",
        "SPRY4", "ETV4", "ETV5", "CCND1", "FOS", "EGR1", "PHLDA1", "IER2"]
APOPTOSIS_BCL2 = ["BCL2", "BCL2L1", "MCL1", "BCL2A1", "BCL2L2", "BAX", "BAK1",
                  "BID", "BBC3", "PMAIP1", "BIM", "BAD", "BMF", "HRK",
                  "CASP3", "CASP8", "CASP9", "XIAP", "BIRC5", "BIRC2",
                  "BIRC3", "DIABLO", "APAF1", "CYCS"]
FOSL1_PROGRAM = ["FOSL1", "FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3",
                 "ATF4", "EGR1", "EGR2", "DUSP1", "DUSP6", "CCND1",
                 "BCL2L1", "MCL1", "IL6", "CXCL8", "PTGS2"]
WRITERS_HAT = ["EP300", "CREBBP", "KAT2A", "KAT2B", "KAT5", "KAT6A", "KAT6B",
               "KAT7", "KAT8", "HAT1", "GTF3C4", "TAF1", "TAF1L", "CLOCK",
               "ELP3", "NCOA1", "NCOA2", "NCOA3", "ATF2"]
READERS_BRD = ["BRD2", "BRD3", "BRD4", "BRDT", "BRD7", "BRD9", "BRPF1",
               "BRPF3", "BAZ2A", "BAZ2B", "CECR2", "ATAD2", "ATAD2B"]
ERASERS_HD = ["HDAC1", "HDAC2", "HDAC3", "HDAC4", "HDAC5", "HDAC6", "HDAC7",
              "HDAC8", "HDAC9", "HDAC10", "HDAC11", "SIRT1", "SIRT2",
              "SIRT3", "SIRT6", "SIRT7"]

SETS = {"EMT": EMT, "MYC_TARGETS": MYC_TARGETS, "RTK": RTK, "MAPK": MAPK,
        "APOPTOSIS_BCL2": APOPTOSIS_BCL2, "FOSL1_PROGRAM": FOSL1_PROGRAM,
        "HAT_WRITERS": WRITERS_HAT, "BET_READERS": READERS_BRD,
        "HDAC_ERASERS": ERASERS_HD}


def hypergeom_overlap(deg_genes, geneset, background):
    gs = set(geneset) & set(background)
    q = len(set(deg_genes) & gs)          # drawn successes
    m = len(gs)                            # successes in population
    n = len(background) - m                # failures
    k = len(set(deg_genes))                # draws
    if k == 0 or m == 0:
        return q, m, np.nan, np.nan
    p = hypergeom.sf(q - 1, m + n, m, k)   # P(X >= q)
    fe = (q / k) / (m / (m + n))
    return q, m, p, fe


def main():
    # ---------- 1. resistance DEG (human PANC1 vs PANC1K) ----------
    df = load_matrix(F_MAIN)
    g = assign_groups_main(df)
    print("main matrix:", df.shape, dict(zip(*np.unique(g, return_counts=True))))
    res = deg_table(df, g, "resistant", "parental", "resistant_vs_parental")

    # DEG calls: |log2FC| >= 1 and p < 0.05 (n=3/group: FDR is conservative;
    # report both cutoffs and use p<0.05 for set overlaps)
    up = res[(res.log2FC >= 1) & (res.p < 0.05)]
    dn = res[(res.log2FC <= -1) & (res.p < 0.05)]
    up_fdr = res[(res.log2FC >= 1) & (res.fdr < 0.10)]
    dn_fdr = res[(res.log2FC <= -1) & (res.fdr < 0.10)]
    background = res.gene.tolist()
    print(f"DEG up (p<0.05,|FC|>=2): {len(up)}; down: {len(dn)}")
    print(f"DEG up (FDR<0.10): {len(up_fdr)}; down: {len(dn_fdr)}")

    # ---------- 2. gene-set overlap ----------
    rows = []
    for name, genes in SETS.items():
        q_u, m, p_u, fe_u = hypergeom_overlap(up.gene, genes, background)
        q_d, _, p_d, fe_d = hypergeom_overlap(dn.gene, genes, background)
        hit_u = sorted(set(up.gene) & set(genes))
        hit_d = sorted(set(dn.gene) & set(genes))
        rows.append({"set": name, "set_size_detected": m,
                     "up_overlap": q_u, "up_p": p_u, "up_fold_enrich": fe_u,
                     "up_hits": ";".join(hit_u),
                     "dn_overlap": q_d, "dn_p": p_d, "dn_fold_enrich": fe_d,
                     "dn_hits": ";".join(hit_d)})
    ovl = pd.DataFrame(rows)
    ovl.to_csv(os.path.join(OUT, "geneset_overlap.csv"), index=False)

    # ---------- 3. acetylation machinery panel ----------
    mach = []
    for fam, genes in [("writer", WRITERS_HAT), ("reader", READERS_BRD),
                       ("eraser", ERASERS_HD)]:
        for gene in genes:
            hit = res[res.gene == gene]
            if len(hit):
                r = hit.iloc[0]
                mach.append({"gene": gene, "family": fam,
                             "log2FC_res_vs_par": r.log2FC,
                             "p": r.p, "fdr": r.fdr})
    mach_df = pd.DataFrame(mach).sort_values("log2FC_res_vs_par", ascending=False)
    mach_df.to_csv(os.path.join(OUT, "acetylation_machinery.csv"), index=False)
    print(mach_df[["gene", "family", "log2FC_res_vs_par", "p"]].to_string(index=False))

    # spotlight: paper's key genes
    for gene in ["FOSL1", "EP300", "CREBBP", "BRD2", "BRD3", "BRD4",
                 "MYC", "CDH1", "VIM", "ZEB1", "EGFR", "MET", "BCL2L1",
                 "MCL1", "BCL2"]:
        hit = res[res.gene == gene]
        if len(hit):
            r = hit.iloc[0]
            print(f"{gene:8s} log2FC={r.log2FC:+.2f} p={r.p:.2g} fdr={r.fdr:.2g}")
        else:
            print(f"{gene:8s} NOT DETECTED")

    # ---------- 4. BET-inhibitor reversal ----------
    for f, tag in [(F_JQ1, "jq1"), (F_OTX, "otx")]:
        d2 = load_matrix(f)
        g2 = assign_groups_bet(d2)
        print(tag, "matrix:", d2.shape,
              dict(zip(*np.unique(g2, return_counts=True))))
        t = deg_table(d2, g2, "treated", "resistant_base", f"panc1k_{tag}_vs_base")
        lut = dict(zip(t.gene, t.log2FC))
        for dset, dname in [(up, "up"), (dn, "down")]:
            vals = [(g_, lut.get(g_, np.nan)) for g_ in dset.gene]
            vals = [(a, b) for a, b in vals if np.isfinite(b)]
            if not vals:
                continue
            genes_, tl = zip(*vals)
            opposite = [1 if (b < 0 if dname == "up" else b > 0) else 0
                        for b in tl]
            bt = binomtest(sum(opposite), len(opposite), 0.5,
                           alternative="greater")
            print(f"{tag}: resistance-{dname} DEGs reversed by treatment: "
                  f"{sum(opposite)}/{len(opposite)} "
                  f"p_binom={bt.pvalue:.2g}")
            pd.DataFrame({"gene": genes_, f"{tag}_log2FC": tl,
                          "reversed": opposite}).to_csv(
                os.path.join(OUT, f"reversal_{tag}_{dname}.csv"), index=False)

    # ---------- 5. treatment spot-check on paper-relevant genes ----------
    SPOT = ["FOSL1", "FOS", "JUN", "JUNB", "EP300", "BRD4",
            "FGFR2", "FGFR3", "FGFR4", "COL3A1", "MYC", "SPRY2",
            "DUSP6", "EGFR"]
    spot_rows = []
    for f, tag, tr in [(F_SGC, "sgc-cbp30", ".SGC."),
                       (F_JQ1, "jq1", ".JQ1."),
                       (F_OTX, "otx-015", ".OTX.")]:
        d2 = load_matrix(f)
        g2 = np.array(["treated" if tr in c else "base" for c in d2.columns])
        m2 = np.log2(d2.values + 1.0)
        for gene in SPOT:
            if gene in d2.index:
                a = m2[d2.index.get_loc(gene)]
                lfc = a[g2 == "treated"].mean() - a[g2 == "base"].mean()
                p = ttest_ind(a[g2 == "treated"], a[g2 == "base"],
                              equal_var=False).pvalue
                spot_rows.append({"treatment": tag, "gene": gene,
                                  "log2FC_vs_resistant_base": lfc, "p": p})
    pd.DataFrame(spot_rows).to_csv(
        os.path.join(OUT, "treatment_spotcheck.csv"), index=False)
    volcano(res, up, dn)
    machinery_heatmap(df, g, mach_df)


def volcano(res, up, dn):
    fig, ax = plt.subplots(figsize=(7, 5))
    x = res.log2FC.values
    y = -np.log10(res.p.clip(lower=1e-300))
    ax.scatter(x, y, s=6, c="#9db8d2", alpha=0.7, edgecolors="none", label="ns")
    m_u = res.gene.isin(up.gene)
    m_d = res.gene.isin(dn.gene)
    ax.scatter(x[m_u], y[m_u], s=10, c="#d62728", alpha=0.8,
               edgecolors="none", label=f"up in resistant (n={m_u.sum()})")
    ax.scatter(x[m_d], y[m_d], s=10, c="#1f77b4", alpha=0.8,
               edgecolors="none", label=f"down in resistant (n={m_d.sum()})")
    for gene in ["FOSL1", "EP300", "BRD4", "VIM", "CDH1", "MYC", "ZEB1"]:
        hit = res[res.gene == gene]
        if len(hit):
            r = hit.iloc[0]
            ax.annotate(gene, (r.log2FC, -np.log10(max(r.p, 1e-300))),
                        fontsize=8, textcoords="offset points", xytext=(4, 4))
    ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c="gray")
    ax.axvline(1, ls="--", lw=0.8, c="gray")
    ax.axvline(-1, ls="--", lw=0.8, c="gray")
    ax.set_xlabel("log2FC PANC1K (resistant) / PANC1 (parental)")
    ax.set_ylabel("-log10 p (Welch t, log2 counts)")
    ax.set_title("GSE294216: MRTX1133-resistant vs parental PANC1 (n=3/group)")
    ax.legend(fontsize=8, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "volcano_resistant_vs_parental.png"), dpi=150)
    plt.close(fig)


def machinery_heatmap(df, groups, mach_df):
    genes = [g for g in mach_df.gene
             if abs(mach_df[mach_df.gene == g].iloc[0].log2FC_res_vs_par) >= 0.5]
    if not genes:
        genes = mach_df.gene.tolist()
    sub = np.log2(df.loc[genes].values + 1.0)
    z = (sub - sub.mean(axis=1, keepdims=True)) / (sub.std(axis=1, keepdims=True) + 1e-9)
    order = np.argsort(groups)
    labels = ["parental" if groups[i] == "parental" else "resistant"
              for i in order]
    fam = dict(zip(mach_df.gene, mach_df.family))
    fig, ax = plt.subplots(figsize=(6, max(4, 0.32 * len(genes))))
    im = ax.imshow(z[:, order], aspect="auto", cmap="RdBu_r",
                   vmin=-2.5, vmax=2.5)
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels([f"{gn} ({fam.get(gn, '?')})" for gn in genes],
                       fontsize=8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8, rotation=45, ha="right")
    ax.set_title("Histone-acetylation machinery: PANC1 vs PANC1K\n"
                 "(z-scored log2 counts; |log2FC|>=0.5)")
    fig.colorbar(im, ax=ax, label="z-score")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "heatmap_acetylation_machinery.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    main()
