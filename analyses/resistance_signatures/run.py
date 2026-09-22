#!/usr/bin/env python3
"""Derive public-data KRAS-inhibitor resistance expression signatures.

Datasets:
  GSE269985 (PMID 38975874) — human PDAC cell lines (Panc 02.03, PANC-1) and
      murine KPC lines (6694C2, 6499C4, 6419C5), parental vs MRTX1133-resistant
      (G12D inhibitor arm, closest to daraxonrasib biology). Processed
      per-sample count matrices from GEO RAW tar (no SRA needed).
  GSE269313 (PMID 38975874) — companion in vivo KPC arm. INFEASIBLE here:
      deposited processed data are per-sample snRNAseq RDS objects and the
      sandbox has no R/pyreadr to parse them (see report.md).
  GSE324125 (PMID 42415848) — human NSCLC H358 RMC-6236-resistant (RMC-R),
      tet-MUC1 shRNA DOX- (control) vs DOX+ (MUC1 silencing). NOTE: the
      deposited series has NO parental H358 — only the 6 RMC-R samples — so
      no direct parental-vs-resistant DEG is possible from deposited data.
      Derived signature is the MUC1-C-supported program within RMC-R cells.

Method: Welch t-test on log2(CPM+1) or log2(TPM+1), Benjamini-Hochberg FDR.
NOT DESeq2/edgeR (task requested a simple approach). Signatures = genes with
|log2FC| >= 1 and FDR < 0.05, split into UP (higher in resistant) / DOWN.
Transfer tests use the repo scorer (score_signatures.py) unchanged; proposed
v2 additions are written to output/ (signatures_v1.yaml is NOT edited).

Usage: python3 run.py   (re-runnable; reads data/, writes output/)
"""
import gzip
import io
import os
import subprocess
import sys
import tarfile

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "output")
SCORER = os.path.normpath(os.path.join(HERE, "..", "signature_panel",
                                        "score_signatures.py"))
os.makedirs(OUT, exist_ok=True)

FC_CUT = 1.0
FDR_CUT = 0.05


# ---------------- helpers ----------------

def collapse_duplicates(df):
    """Keep the row with max mean for duplicated gene symbols
    (e.g. Excel-mangled '2-Mar'/'1-Mar')."""
    if df.index.is_unique:
        return df
    means = df.mean(axis=1)
    keep = means.groupby(level=0).idxmax()
    return df.loc[keep.values]


def bh_fdr(p):
    p = np.asarray(p, dtype=float)
    n = len(p)
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / np.maximum(np.arange(1, n + 1), 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.clip(q, 0, 1)
    return out


def to_cpm(counts):
    lib = counts.sum(axis=0).replace(0, np.nan)
    return counts.div(lib, axis=1) * 1e6


def deg_welch(logexpr, labels, g1, g2, label1="group1", label2="group2"):
    """Welch t-test, genes x samples. log2FC = mean(g2) - mean(g1)."""
    cols1 = [c for c, l in zip(logexpr.columns, labels) if l == g1]
    cols2 = [c for c, l in zip(logexpr.columns, labels) if l == g2]
    a = logexpr[cols1].to_numpy(float)
    b = logexpr[cols2].to_numpy(float)
    with np.errstate(all="ignore"):
        _, p = stats.ttest_ind(b, a, axis=1, equal_var=False, nan_policy="omit")
    p = np.where(np.isnan(p), 1.0, p)
    res = pd.DataFrame({
        "gene": logexpr.index,
        f"mean_{label1}": np.nanmean(a, axis=1),
        f"mean_{label2}": np.nanmean(b, axis=1),
        "log2FC": np.nanmean(b, axis=1) - np.nanmean(a, axis=1),
        "p_value": p,
        "padj": bh_fdr(p),
        f"n_{label1}": len(cols1),
        f"n_{label2}": len(cols2),
    }).set_index("gene").sort_values("padj")
    return res


def make_signature(deg, direction="up"):
    if direction == "up":
        sel = deg[(deg["log2FC"] >= FC_CUT) & (deg["padj"] < FDR_CUT)]
    else:
        sel = deg[(deg["log2FC"] <= -FC_CUT) & (deg["padj"] < FDR_CUT)]
    return sel.index.tolist(), sel


def per_line_concordance(logexpr, meta, line_col="line", cond_col="condition"):
    """Per-cell-line resistant-vs-parental log2FC for the pooled signature."""
    out = {}
    for line, sub in meta.groupby(line_col):
        cols = sub.index
        p = [c for c in cols if sub.loc[c, cond_col] == "parental"]
        r = [c for c in cols if sub.loc[c, cond_col] == "resistant"]
        if p and r:
            out[line] = logexpr[r].mean(axis=1) - logexpr[p].mean(axis=1)
    return pd.DataFrame(out)


def zscore_rows(df):
    return df.sub(df.mean(axis=1), axis=0).div(
        df.std(axis=1, ddof=1).replace(0, np.nan), axis=0)


def volcano(deg, title, path):
    x = deg["log2FC"].to_numpy()
    y = -np.log10(np.clip(deg["padj"].to_numpy(), 1e-300, 1))
    sig = (np.abs(x) >= FC_CUT) & (deg["padj"].to_numpy() < FDR_CUT)
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(x[~sig], y[~sig], s=6, c="#999999", alpha=0.5, rasterized=True)
    ax.scatter(x[sig], y[sig], s=10, c="#c0392b", alpha=0.8, rasterized=True)
    ax.axvline(FC_CUT, ls="--", c="k", lw=0.8)
    ax.axvline(-FC_CUT, ls="--", c="k", lw=0.8)
    ax.axhline(-np.log10(FDR_CUT), ls="--", c="k", lw=0.8)
    for _, r in deg[sig].head(12).iterrows():
        ax.annotate(r.name, (r["log2FC"], -np.log10(max(r["padj"], 1e-300))),
                    fontsize=6, alpha=0.8)
    ax.set_xlabel("log2FC (resistant / parental)")
    ax.set_ylabel("-log10 FDR")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def heatmap_top(logexpr, deg, labels, path, title, top_n=40):
    sig = deg[(deg["padj"] < FDR_CUT) & (np.abs(deg["log2FC"]) >= FC_CUT)]
    top = sig.reindex(sig["log2FC"].abs().sort_values(ascending=False).index).head(top_n)
    genes = [g for g in top.index if g in logexpr.index]
    if not genes:
        print("  heatmap: no significant genes, skipping", path)
        return
    z = zscore_rows(logexpr.loc[genes])
    # order columns by group
    order = np.argsort([0 if l in ("parental", "DOX-") else 1 for l in labels])
    z = z.iloc[:, order]
    fig, ax = plt.subplots(figsize=(max(8, 0.4 * z.shape[1]),
                                    max(6, 0.22 * len(genes))))
    im = ax.imshow(z.to_numpy(), aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2,
                   interpolation="nearest")
    ax.set_yticks(range(len(genes)))
    ax.set_yticklabels(genes, fontsize=6)
    ax.set_xticks(range(z.shape[1]))
    ax.set_xticklabels([f"{c}\n{l}" for c, l in
                        zip(z.columns, [labels[i] for i in order])],
                       fontsize=6, rotation=45, ha="right")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, label="per-gene z-score")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def mouse_to_human(genes):
    """Naive ortholog guess (Kras -> KRAS). Documented limitation."""
    return [g.upper() for g in genes]


def run_scorer(input_csv, yaml_path, organism, out_csv, method="rank"):
    r = subprocess.run(
        [sys.executable, SCORER, "--input", input_csv, "--signatures", yaml_path,
         "--organism", organism, "--out", out_csv, "--method", method,
         "--no-log"],
        capture_output=True, text=True)
    print(r.stdout.strip().splitlines()[-1] if r.stdout.strip() else r.stderr)
    r.check_returncode()
    return pd.read_csv(out_csv, index_col=0)


def group_means(scores, labels):
    df = scores.copy()
    df["group"] = labels
    return df.groupby("group").mean(numeric_only=True)


# ---------------- dataset: GSE324125 ----------------

def load_gse324125():
    counts = pd.read_csv(
        os.path.join(DATA, "GSE324125_STAR_Gene_Counts_H358_RMC-R_DOX-_vs_DOX+.csv.gz"),
        index_col=0)
    tpm = pd.read_csv(
        os.path.join(DATA, "GSE324125_STAR_Gene_TPM_H358_RMC-R_DOX-_vs_DOX+.csv.gz"),
        index_col=0)
    counts = collapse_duplicates(counts)
    tpm = collapse_duplicates(tpm)
    common = counts.index.intersection(tpm.index)
    tpm = tpm.loc[common]
    labels = ["DOX+" if "Silencing" in c else "DOX-" for c in tpm.columns]
    return tpm, labels


# ---------------- dataset: GSE269985 ----------------

HUMAN_SAMPLES = ["GSM8332065", "GSM8332066", "GSM8332067", "GSM8332068",
                 "GSM8332069", "GSM8332070", "GSM8332071", "GSM8332072",
                 "GSM8332073", "GSM8332074", "GSM8332075", "GSM8332076"]
MOUSE_SAMPLES = ["GSM8332077", "GSM8332078", "GSM8332079", "GSM8332080",
                 "GSM8332081", "GSM8332082", "GSM8332083", "GSM8332084",
                 "GSM8332085", "GSM8332086", "GSM8332087", "GSM8332088",
                 "GSM8332089", "GSM8332090", "GSM8332091", "GSM8332092",
                 "GSM8332093", "GSM8332094"]


def load_gse269985_arm(tar_path, accessions):
    """Build genes x samples raw-count matrix from per-sample member files."""
    cols = {}
    with tarfile.open(tar_path, "r:") as tf:
        members = {m.name: m for m in tf.getmembers() if m.isfile()}
        for acc in accessions:
            hit = [n for n in members if n.startswith(acc + "_")]
            assert len(hit) == 1, (acc, hit)
            raw = gzip.decompress(tf.extractfile(members[hit[0]]).read())
            df = pd.read_csv(io.BytesIO(raw), index_col=0)  # gene | ensembl_id | sample
            scol = df.columns[-1]
            s = df[scol]
            s.name = acc
            cols[acc] = s
    mat = pd.DataFrame(cols)
    mat.index.name = "gene"
    return collapse_duplicates(mat)


def meta_gse269985_human():
    rows = []
    for i, acc in enumerate(HUMAN_SAMPLES):
        line = "Panc02.03" if i < 6 else "PANC-1"
        cond = "parental" if (i % 6) < 3 else "resistant"
        rows.append({"sample": acc, "line": line, "condition": cond})
    return pd.DataFrame(rows).set_index("sample")


def meta_gse269985_mouse():
    rows = []
    lines = ["6694C2", "6499C4", "6419C5"]
    for i, acc in enumerate(MOUSE_SAMPLES):
        line = lines[i // 6]
        cond = "parental" if (i % 6) < 3 else "resistant"
        rows.append({"sample": acc, "line": line, "condition": cond})
    return pd.DataFrame(rows).set_index("sample")


# ---------------- main ----------------

def main():
    summary = {}

    # ---- GSE324125: H358 RMC-R, MUC1 silencing ----
    print("== GSE324125 (H358 RMC-R, DOX- vs DOX+ MUC1 kd) ==")
    tpm324, lab324 = load_gse324125()
    log324 = np.log2(tpm324 + 1)
    deg324 = deg_welch(log324, lab324, "DOX-", "DOX+", "DOXminus", "DOXplus_MUC1kd")
    deg324.to_csv(os.path.join(OUT, "deg_GSE324125_DOXplus_vs_DOXminus.csv"))
    muc1_dep, _ = make_signature(deg324, "down")   # lower on kd = MUC1-supported
    muc1_derep, _ = make_signature(deg324, "up")   # higher on kd = de-repressed
    summary["GSE324125_n_genes"] = len(deg324)
    summary["GSE324125_MUC1_supported"] = len(muc1_dep)
    summary["GSE324125_MUC1_derepressed"] = len(muc1_derep)
    print(f"  genes={len(deg324)} MUC1-supported={len(muc1_dep)} "
          f"de-repressed={len(muc1_derep)}")
    pd.DataFrame({"gene": muc1_dep}).to_csv(
        os.path.join(OUT, "sig_H358_RMCR_MUC1_supported.csv"), index=False)
    pd.DataFrame({"gene": muc1_derep}).to_csv(
        os.path.join(OUT, "sig_H358_RMCR_MUC1_derepressed.csv"), index=False)
    volcano(deg324, "GSE324125: H358 RMC-R, MUC1 silencing (DOX+) vs control (DOX-)",
            os.path.join(OUT, "volcano_GSE324125.png"))
    heatmap_top(log324, deg324, lab324,
                os.path.join(OUT, "heatmap_GSE324125.png"),
                "GSE324125 top DEGs (MUC1 kd vs control)")
    tpm324.to_csv(os.path.join(OUT, "expr_GSE324125_TPM.csv"))

    # ---- GSE269985 human arm ----
    print("== GSE269985 human arm (Panc02.03 + PANC-1, MRTX1133) ==")
    tar985 = os.path.join(DATA, "GSE269985_RAW.tar")
    h_counts = load_gse269985_arm(tar985, HUMAN_SAMPLES)
    h_meta = meta_gse269985_human()
    h_cpm = to_cpm(h_counts)
    # light filter: mean CPM > 0.5
    keep = h_cpm.mean(axis=1) > 0.5
    h_cpm, h_counts = h_cpm[keep], h_counts[keep]
    h_log = np.log2(h_cpm + 1)
    h_labels = [h_meta.loc[c, "condition"] for c in h_log.columns]
    deg_h = deg_welch(h_log, h_labels, "parental", "resistant",
                      "parental", "resistant")
    deg_h.to_csv(os.path.join(OUT, "deg_GSE269985_human_resistant_vs_parental.csv"))
    h_up, _ = make_signature(deg_h, "up")
    h_dn, _ = make_signature(deg_h, "down")
    h_conc = per_line_concordance(h_log, h_meta)
    h_conc.to_csv(os.path.join(OUT, "perline_log2FC_GSE269985_human.csv"))
    summary.update({"GSE269985_human_genes": len(deg_h),
                    "GSE269985_human_up": len(h_up),
                    "GSE269985_human_down": len(h_dn)})
    print(f"  genes={len(deg_h)} up={len(h_up)} down={len(h_dn)}")
    pd.DataFrame({"gene": h_up}).to_csv(
        os.path.join(OUT, "sig_PDAC_KRASi_resistance_UP.csv"), index=False)
    pd.DataFrame({"gene": h_dn}).to_csv(
        os.path.join(OUT, "sig_PDAC_KRASi_resistance_DOWN.csv"), index=False)
    volcano(deg_h, "GSE269985 human: MRTX1133-resistant vs parental (Panc02.03 + PANC-1)",
            os.path.join(OUT, "volcano_GSE269985_human.png"))
    heatmap_top(h_log, deg_h, h_labels,
                os.path.join(OUT, "heatmap_GSE269985_human.png"),
                "GSE269985 human top DEGs (resistant vs parental)")
    h_cpm.to_csv(os.path.join(OUT, "expr_GSE269985_human_CPM.csv"))

    # ---- GSE269985 mouse arm ----
    print("== GSE269985 mouse arm (6694C2, 6499C4, 6419C5, MRTX1133) ==")
    m_counts = load_gse269985_arm(tar985, MOUSE_SAMPLES)
    m_meta = meta_gse269985_mouse()
    m_cpm = to_cpm(m_counts)
    keep = m_cpm.mean(axis=1) > 0.5
    m_cpm = m_cpm[keep]
    m_log = np.log2(m_cpm + 1)
    m_labels = [m_meta.loc[c, "condition"] for c in m_log.columns]
    deg_m = deg_welch(m_log, m_labels, "parental", "resistant",
                      "parental", "resistant")
    deg_m.to_csv(os.path.join(OUT, "deg_GSE269985_mouse_resistant_vs_parental.csv"))
    m_up, _ = make_signature(deg_m, "up")
    m_dn, _ = make_signature(deg_m, "down")
    m_conc = per_line_concordance(m_log, m_meta)
    m_conc.to_csv(os.path.join(OUT, "perline_log2FC_GSE269985_mouse.csv"))
    summary.update({"GSE269985_mouse_genes": len(deg_m),
                    "GSE269985_mouse_up": len(m_up),
                    "GSE269985_mouse_down": len(m_dn)})
    print(f"  genes={len(deg_m)} up={len(m_up)} down={len(m_dn)}")
    pd.DataFrame({"gene": m_up, "human_ortholog_guess": mouse_to_human(m_up)}).to_csv(
        os.path.join(OUT, "sig_KPC_KRASi_resistance_UP.csv"), index=False)
    pd.DataFrame({"gene": m_dn, "human_ortholog_guess": mouse_to_human(m_dn)}).to_csv(
        os.path.join(OUT, "sig_KPC_KRASi_resistance_DOWN.csv"), index=False)
    volcano(deg_m, "GSE269985 mouse: MRTX1133-resistant vs parental (3 KPC lines)",
            os.path.join(OUT, "volcano_GSE269985_mouse.png"))
    heatmap_top(m_log, deg_m, m_labels,
                os.path.join(OUT, "heatmap_GSE269985_mouse.png"),
                "GSE269985 mouse top DEGs (resistant vs parental)")
    m_cpm.to_csv(os.path.join(OUT, "expr_GSE269985_mouse_CPM.csv"))

    # ---- proposed v2 additions yaml (NOT editing signatures_v1.yaml) ----
    yaml_path = os.path.join(OUT, "signatures_v2_proposed_additions.yaml")
    with open(yaml_path, "w") as f:
        f.write("# PROPOSED additions to signatures_v1.yaml — NOT yet curated in.\n")
        f.write("# Derived 2026-09-22 by analyses/resistance_signatures/run.py\n")
        f.write("# Method: Welch t-test on log2(CPM/TPM+1), |log2FC|>=1, FDR<0.05.\n")
        f.write("# Mouse->human orthologs are naive upper() guesses (see report).\n")
        f.write("version: \"1.0+proposed\"\n\nsignatures:\n\n")
        def block(name, label, desc, human_genes, mouse_genes, direction, src):
            f.write(f"  {name}:\n")
            f.write(f"    label: \"{label}\"\n")
            f.write(f"    description: >-\n      {desc}\n")
            f.write(f"    source:\n      first_author: \"{src['author']}\"\n")
            f.write(f"      year: {src['year']}\n      journal: \"{src['journal']}\"\n")
            f.write(f"      pmid: {src['pmid']}\n")
            f.write(f"      geo: \"{src['geo']}\"\n")
            f.write(f"      method: \"Welch t-test on log2(CPM/TPM+1); |log2FC|>=1, FDR<0.05\"\n")
            f.write(f"    human: [{', '.join(human_genes)}]\n")
            f.write(f"    mouse: [{', '.join(mouse_genes)}]\n")
            f.write(f"    direction: \"{direction}\"\n\n")
        block("PDAC_KRASi_resistance_UP", "PDAC KRASi-resistance program (up)",
              "Genes consistently upregulated in MRTX1133-resistant vs parental "
              "PDAC lines (Panc 02.03 + PANC-1, GSE269985, human arm).",
              h_up, [g.capitalize() for g in h_up],
              "increases with acquired KRASi resistance",
              {"author": "Kumar et al.", "year": 2024,
               "journal": "Cancer Cell (PMID 38975874)", "pmid": 38975874,
               "geo": "GSE269985"})
        block("PDAC_KRASi_resistance_DOWN", "PDAC KRASi-resistance program (down)",
              "Genes consistently downregulated in MRTX1133-resistant vs parental "
              "PDAC lines (Panc 02.03 + PANC-1, GSE269985, human arm).",
              h_dn, [g.capitalize() for g in h_dn],
              "decreases with acquired KRASi resistance",
              {"author": "Kumar et al.", "year": 2024,
               "journal": "Cancer Cell (PMID 38975874)", "pmid": 38975874,
               "geo": "GSE269985"})
        block("KPC_KRASi_resistance_UP", "KPC KRASi-resistance program (up)",
              "Genes upregulated in MRTX1133-resistant vs parental murine KPC "
              "lines (6694C2, 6499C4, 6419C5, GSE269985, mouse arm). Human list "
              "is a naive upper-case ortholog guess.",
              mouse_to_human(m_up), m_up,
              "increases with acquired KRASi resistance",
              {"author": "Kumar et al.", "year": 2024,
               "journal": "Cancer Cell (PMID 38975874)", "pmid": 38975874,
               "geo": "GSE269985"})
        block("H358_RMCR_MUC1_supported", "MUC1-C-supported program in RMC-6236-resistant H358",
              "Genes downregulated when MUC1-C is silenced (DOX+) in "
              "RMC-6236-resistant H358 (GSE324125) — programs the resistant "
              "line maintains via MUC1-C. NOT a parental-vs-resistant signature "
              "(no parental H358 was deposited).",
              muc1_dep, [g.capitalize() for g in muc1_dep],
              "supported by MUC1-C in RMC-6236-resistant cells",
              {"author": "MUC1-C/RAS(ON) study", "year": 2026,
               "journal": "(PMID 42415848)", "pmid": 42415848,
               "geo": "GSE324125"})
    print("wrote", yaml_path)

    # ---- transfer tests with the repo scorer ----
    print("== transfer tests (rank-based scoring) ==")
    # T1: H358 MUC1-supported signature on PDAC human + mouse matrices
    s1h = run_scorer(os.path.join(OUT, "expr_GSE269985_human_CPM.csv"),
                     yaml_path, "human",
                     os.path.join(OUT, "scores_T1_H358sig_on_GSE269985human.csv"))
    s1m = run_scorer(os.path.join(OUT, "expr_GSE269985_mouse_CPM.csv"),
                     yaml_path, "mouse",
                     os.path.join(OUT, "scores_T1_H358sig_on_GSE269985mouse.csv"))
    t1h = group_means(s1h, h_labels)
    t1m = group_means(s1m, m_labels)
    t1h.to_csv(os.path.join(OUT, "transfer_T1_groupmeans_human.csv"))
    t1m.to_csv(os.path.join(OUT, "transfer_T1_groupmeans_mouse.csv"))
    print("T1 human:\n", t1h.round(3).to_string())
    print("T1 mouse:\n", t1m.round(3).to_string())

    # T2: PDAC + KPC resistance signatures on H358 RMC-R (DOX- vs DOX+)
    s2 = run_scorer(os.path.join(OUT, "expr_GSE324125_TPM.csv"),
                    yaml_path, "human",
                    os.path.join(OUT, "scores_T2_PDACsig_on_H358.csv"))
    t2 = group_means(s2, lab324)
    t2.to_csv(os.path.join(OUT, "transfer_T2_groupmeans.csv"))
    print("T2 (all samples are RMC-R; DOX+ = MUC1 kd):\n", t2.round(3).to_string())

    # Mann-Whitney on transfer scores
    mw = []
    for name, scores, labels, ga, gb in [
            ("T1_H358_MUC1_supported_on_GSE269985_human", s1h, h_labels, "resistant", "parental"),
            ("T1_H358_MUC1_supported_on_GSE269985_mouse", s1m, m_labels, "resistant", "parental"),
            ("T2_PDAC_UP_on_H358_DOXminus_vs_DOXplus", s2, lab324, "DOX-", "DOX+"),
            ("T2_PDAC_DOWN_on_H358_DOXminus_vs_DOXplus", s2, lab324, "DOX-", "DOX+")]:
        for sig in scores.columns:
            a = scores.loc[np.array(labels) == ga, sig]
            b = scores.loc[np.array(labels) == gb, sig]
            try:
                _, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            except Exception:
                p = np.nan
            mw.append({"test": name, "signature": sig,
                       f"mean_{ga}": a.mean(), f"mean_{gb}": b.mean(), "p": p})
    pd.DataFrame(mw).to_csv(os.path.join(OUT, "transfer_mannwhitney.csv"), index=False)
    print(pd.DataFrame(mw).round(4).to_string())

    pd.DataFrame([summary]).T.to_csv(os.path.join(OUT, "run_summary.csv"),
                                     header=["value"])
    print("SUMMARY:", summary)
    print("done.")


if __name__ == "__main__":
    main()
