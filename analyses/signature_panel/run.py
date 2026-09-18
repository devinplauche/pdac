#!/usr/bin/env python3
"""Reproduce the headline analysis of GSE291118 with the signature compendium.

GSE291118 (Zhang et al., Cancer Cell 2025, PMID 40712568):
"Macropinocytosis maintains CAF subtype identity under metabolic stress in
pancreatic cancer." Bulk RNA-seq (TPM) of immortalized murine CAFs, n=5 per
condition, 32 h treatment:

    02Qveh : 0.2 mM glutamine + DMSO        (metabolic stress, macropinocytosis ON)
    02Qeipa: 0.2 mM glutamine + 25 uM EIPA  (metabolic stress, macropinocytosis BLOCKED)
    4Qveh  : 4 mM glutamine + DMSO          (nutrient replete, baseline)
    4Qil1a : 4 mM glutamine + 25 pg/mL IL-1a (positive control: IL-1 drives iCAF)

Paper's headline claims to reproduce at signature level:
  1. Metabolic stress induces an intrinsic iCAF program via MEK-ERK
     -> 02Qveh vs 4Qveh: iCAF score UP, MAPK_targets UP, SASP UP.
  2. Blocking macropinocytosis promotes myCAF-to-iCAF transitions and reduces collagen
     -> 02Qeipa vs 02Qveh: myCAF DOWN, iCAF UP, matrisome_core DOWN.
  3. IL-1a drives iCAF (Biffi et al. 2019) -> 4Qil1a vs 4Qveh: iCAF UP (positive control).

Statistics: Mann-Whitney U per contrast (n=5/group); direction = median difference.
"""
import os
import urllib.request

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from score_signatures import load_signatures, score_zscore

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)

TPM_URL = ("https://ftp.ncbi.nlm.nih.gov/geo/series/GSE291nnn/GSE291118/"
           "suppl/GSE291118_TPM.txt.gz")
TPM_LOCAL = os.path.join(OUT, "GSE291118_TPM.txt.gz")

CONDS = {
    "02Qveh": "0.2 mM glutamine + DMSO (stress, MP on)",
    "02Qeipa": "0.2 mM glutamine + EIPA (stress, MP blocked)",
    "4Qveh": "4 mM glutamine + DMSO (baseline)",
    "4Qil1a": "4 mM glutamine + IL-1a (iCAF positive control)",
}
# (test_group, ref_group, expected direction per paper)
CONTRASTS = [
    ("02Qveh", "4Qveh", "iCAF UP, MAPK UP, SASP UP (stress induces iCAF program)"),
    ("02Qeipa", "02Qveh", "myCAF DOWN, iCAF UP, matrisome DOWN (MP block -> transition)"),
    ("4Qil1a", "4Qveh", "iCAF UP (IL-1a positive control)"),
]


def cond_of(col):
    core = col.rstrip("0123456789")  # drop replicate number: 02Qeipa01 -> 02Qeipa
    return {"02Qveh": "02Qveh", "02Qeipa": "02Qeipa",
            "4Qveh": "4Qveh", "4Qil1a": "4Qil1a"}.get(core)


def main():
    if not os.path.exists(TPM_LOCAL):
        print("downloading GSE291118 TPM matrix...")
        urllib.request.urlretrieve(TPM_URL, TPM_LOCAL)
    tpm = pd.read_csv(TPM_LOCAL, sep="\t")
    expr = tpm.set_index("gene_name")
    keep = [c for c in expr.columns if cond_of(c) in CONDS]
    expr = expr[keep]
    expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)
    groups = {c: cond_of(c) for c in expr.columns}
    assert sorted(set(groups.values())) == sorted(CONDS), set(groups.values())
    print(f"matrix: {expr.shape[0]} genes x {expr.shape[1]} samples")

    logexpr = np.log2(expr + 1)
    sigs, version = load_signatures(os.path.join(HERE, "signatures_v1.yaml"), "mouse")
    scores, coverage = {}, {}
    for name, s in sigs.items():
        sc, n = score_zscore(logexpr, s["genes"])
        scores[name] = sc
        coverage[name] = {"n_genes": len(s["genes"]), "n_found": n,
                          "missing": sorted(set(s["genes"]) - set(logexpr.index))}
    S = pd.DataFrame(scores)
    S["condition"] = S.index.map(groups)
    S.index.name = "sample"
    S.to_csv(os.path.join(OUT, "signature_scores.csv"))
    pd.DataFrame(coverage).T.to_csv(os.path.join(OUT, "gene_coverage.csv"))
    print("gene coverage:")
    for k, v in coverage.items():
        print(f"  {k}: {v['n_found']}/{v['n_genes']} found"
              + (f" (missing: {v['missing']})" if v["missing"] else ""))

    # ---- contrasts ----
    rows = []
    for test, ref, note in CONTRASTS:
        a = S[S["condition"] == test]
        b = S[S["condition"] == ref]
        for sig in scores:
            stat, p = mannwhitneyu(a[sig], b[sig], alternative="two-sided")
            rows.append({"contrast": f"{test} vs {ref}", "signature": sig,
                         "median_test": a[sig].median(), "median_ref": b[sig].median(),
                         "delta_median": a[sig].median() - b[sig].median(),
                         "direction": "UP" if a[sig].median() > b[sig].median() else "DOWN",
                         "p_mwu": p, "n": f"{len(a)}/{len(b)}",
                         "paper_expectation": note})
    R = pd.DataFrame(rows)
    R.to_csv(os.path.join(OUT, "contrast_stats.csv"), index=False)
    print(R[["contrast", "signature", "direction", "p_mwu"]].to_string(index=False,
          float_format="%.4f"))

    # ---- figure: signature scores per condition ----
    order = ["4Qveh", "02Qveh", "02Qeipa", "4Qil1a"]
    labels = ["4 mM Gln\n+ DMSO", "0.2 mM Gln\n+ DMSO", "0.2 mM Gln\n+ EIPA",
              "4 mM Gln\n+ IL-1a"]
    sigs_plot = ["myCAF", "iCAF", "apCAF", "MAPK_targets", "SASP", "matrisome_core",
                 "RASON_inhibition_PD"]
    titles = {"myCAF": "myCAF (Elyada)", "iCAF": "iCAF (Elyada)",
              "apCAF": "apCAF (Elyada)", "MAPK_targets": "MAPK/ERK output (Pratilas)",
              "SASP": "SASP (Coppe)", "matrisome_core": "Core matrisome (Naba)",
              "RASON_inhibition_PD": "RAS(ON) inhibition PD (Holderfield)"}
    fig, axes = plt.subplots(3, 3, figsize=(13, 10), sharex=True)
    axes = axes.flat
    for ax, sig in zip(axes, sigs_plot):
        data = [S[S["condition"] == c][sig].values for c in order]
        bp = ax.boxplot(data, labels=labels, patch_artist=True)
        for patch in bp["boxes"]:
            patch.set_facecolor("#aec7e8")
        ax.set_title(titles[sig], fontsize=10)
        ax.tick_params(axis="x", labelsize=8)
    for ax in axes[len(sigs_plot):]:
        ax.axis("off")
    fig.suptitle("GSE291118 signature scores per condition (mean z-score, n=5/group)",
                 fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "signature_scores_by_condition.png"), dpi=150)
    print("wrote output/: signature_scores.csv, gene_coverage.csv, "
          "contrast_stats.csv, signature_scores_by_condition.png")


if __name__ == "__main__":
    main()
