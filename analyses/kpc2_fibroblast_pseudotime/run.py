#!/usr/bin/env python3
"""ANALYSIS B: pseudotime of the KPC2 fibroblast myCAF->iCAF transition.

Question: can the 95 KPC2 fibroblasts (53 RMC6236, 42 vehicle) be ordered
along a myCAF->iCAF axis, and which genes change along it?

Method (deliberately simple, honestly labeled):
  1. Load the sibling cell-level signature scores
     (kpc2_cell_signature_scores.csv, fractional-rank means). The sibling
     run.py stored column indices into the QC-filtered matrix; reload both
     KPC2 matrices with the IDENTICAL QC/normalization (n_genes>=200,
     totals>=500, pct_mt<=15; CP10k, log1p) and map those indices back to
     filtered columns. Fibroblast cells are re-identified with the sibling
     k-means rule as a consistency check.
  2. Pseudotime = iCAF_score - myCAF_score (per cell). Cells binned into 4
     quartile windows (early=myCAF-like -> late=iCAF-like). A PCA-on-
     signature-genes ordering is computed as a sensitivity check and its
     correlation with the score-difference ordering is reported.
  3. For every gene detected in >=10% of the 95 fibroblasts: Spearman rho
     of log-normalized expression vs the continuous pseudotime ordering.
     Top genes in both directions (rising toward iCAF = iCAF drivers;
     falling = myCAF markers) are tabulated with nominal p-values and
     FDR (Benjamini-Hochberg) for reference -- with 95 cells these are
     hypothesis-generating only.
  4. Report window means for the top genes + treatment composition per
     window (the "transition" partly tracks treatment).

Caveat: 95 cells, no real trajectory inference; a myCAF->iCAF axis is
ASSUMED from the score difference, not learned from branching structure.

Re-runnable: python3 run.py   (writes output/ + figures)
"""
import gzip
import os

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
SIBLING = os.path.abspath(os.path.join(HERE, "..", "gse337208_caf_mining"))
DATA = os.path.join(SIBLING, "data")
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)

SAMPLES = {
    "GSM9850799": ("Hy-RMC6236", "KPC2", "RMC6236"),
    "GSM9850800": ("Hy-Veh", "KPC2", "Veh"),
}

FIBRO_MARKERS = ["Col1a1", "Col1a2", "Col3a1", "Dcn", "Pdpn", "Lum", "Sparc"]
EPI_MARKERS = ["Epcam", "Krt19", "Krt8", "Krt18"]
IMM_MARKERS = ["Ptprc", "Cd3e", "Cd68", "Cd79a", "Itgam"]
PANEL = {"fib": FIBRO_MARKERS, "epi": EPI_MARKERS, "imm": IMM_MARKERS}

N_WINDOWS = 4
MIN_DET_FRAC = 0.10


def kmeans(X, k, seed=0, iters=100):
    rng = np.random.default_rng(seed)
    C = X[rng.choice(X.shape[0], k, replace=False)]
    for _ in range(iters):
        d = ((X[:, None, :] - C[None, :, :]) ** 2).sum(-1)
        lab = d.argmin(1)
        newC = np.array([X[lab == i].mean(0) if (lab == i).any() else C[i]
                         for i in range(k)])
        if np.allclose(newC, C):
            break
        C = newC
    return lab


def load_sample(gsm):
    pfx, model, trt = SAMPLES[gsm]
    base = os.path.join(DATA, f"{gsm}_{pfx}")
    with gzip.open(base + "_features.tsv.gz", "rt") as f:
        feats = pd.read_csv(f, sep="\t", header=None,
                            names=["ensembl", "symbol", "type"])
    genes = feats["symbol"].values
    with gzip.open(base + "_barcodes.tsv.gz", "rt") as f:
        barcodes = [l.strip() for l in f]
    mat = mmread(base + "_matrix.mtx.gz").tocsr()
    return genes, mat, trt


def qc_normalize(genes, mat):
    totals = np.asarray(mat.sum(axis=0)).ravel()
    n_genes = np.asarray((mat > 0).sum(axis=0)).ravel()
    is_mt = np.array([g.startswith("mt-") for g in genes])
    mt_counts = np.asarray(mat[is_mt, :].sum(axis=0)).ravel()
    pct_mt = np.where(totals > 0, 100 * mt_counts / totals, 0)
    keep = (n_genes >= 200) & (totals >= 500) & (pct_mt <= 15)
    mat = mat[:, keep]
    totals = totals[keep]
    logn = (mat @ sparse.diags(1e4 / totals)).tocsr()
    logn.data = np.log1p(logn.data)
    return logn


def fibroblast_mask(genes, logn):
    gidx = {g: i for i, g in enumerate(genes)}
    panel_genes = [g for grp in PANEL.values() for g in grp if g in gidx]
    Xp = np.asarray(logn[[gidx[g] for g in panel_genes], :].toarray()).T
    Xz = (Xp - Xp.mean(0)) / np.maximum(Xp.std(0), 1e-8)
    lab = kmeans(Xz, 6, seed=0)
    prog_of = {}
    for grp, gl in PANEL.items():
        cols = [panel_genes.index(g) for g in gl if g in panel_genes]
        prog_of[grp] = Xp[:, cols].mean(1) if cols else np.zeros(Xp.shape[0])
    best, best_margin = -1, -1e9
    for i in range(6):
        m = lab == i
        if not m.any():
            continue
        margin = (prog_of["fib"][m].mean() -
                  max(prog_of["epi"][m].mean(), prog_of["imm"][m].mean()))
        if margin > best_margin:
            best_margin, best = margin, i
    is_fib = lab == best
    is_fib = is_fib & ~(prog_of["epi"] > prog_of["fib"])
    is_fib = is_fib & ~(prog_of["imm"] > prog_of["fib"])
    return is_fib


def main():
    scores = pd.read_csv(
        os.path.join(SIBLING, "output", "kpc2_cell_signature_scores.csv"))
    print(f"score table: {len(scores)} cells", flush=True)

    # assemble log-normalized expression for the scored cells
    cell_frames = []
    for gsm in SAMPLES:
        print(f"--- {gsm}", flush=True)
        genes, mat, trt = load_sample(gsm)
        logn = qc_normalize(genes, mat)
        is_fib = fibroblast_mask(genes, logn)
        fib_cols = np.where(is_fib)[0]
        sub = scores[scores["sample"] == gsm]
        # sibling stored filtered-matrix column indices in "cell"
        cols = sub["cell"].to_numpy().astype(int)
        assert set(cols) <= set(fib_cols), \
            f"cell indices mismatch for {gsm}"
        assert len(cols) == len(set(cols)), f"duplicate cells in {gsm}"
        X = logn.tocsc()[:, cols].toarray().T  # cells x genes
        df = pd.DataFrame(X, columns=genes)
        # collapse duplicated gene symbols (CellRanger dup symbols)
        df = df.T.groupby(level=0).max().T
        df["sample"] = gsm
        df["treatment"] = trt
        df["iCAF"] = sub["iCAF"].to_numpy()
        df["myCAF"] = sub["myCAF"].to_numpy()
        cell_frames.append(df.reset_index(drop=True))
        print(f"  {len(cols)} fibroblast cells matched", flush=True)
    cells = pd.concat(cell_frames, ignore_index=True)
    n = len(cells)
    print(f"total fibroblasts: {n}", flush=True)

    # pseudotime = iCAF - myCAF
    pt = cells["iCAF"] - cells["myCAF"]
    order = np.argsort(pt.to_numpy(), kind="stable")
    cells["pseudotime"] = pt
    cells["pseudotime_rank"] = np.argsort(np.argsort(pt.to_numpy()))
    cells["window"] = pd.qcut(cells["pseudotime"], N_WINDOWS,
                              labels=[f"W{i+1}" for i in range(N_WINDOWS)])
    # sensitivity: PCA on signature genes (the gene-level columns in scores)
    sig_gene_cols = [c for c in scores.columns
                     if c not in ("sample", "treatment", "cell", "myCAF",
                                  "iCAF", "apCAF", "SASP", "matrisome_core",
                                  "MAPK_targets", "RASON_inhibition_PD")]
    S = cells[sig_gene_cols].to_numpy()
    Sz = (S - S.mean(0)) / np.maximum(S.std(0), 1e-8)
    U, svals, _ = np.linalg.svd(Sz - Sz.mean(0), full_matrices=False)
    pc1 = U[:, 0] * svals[0]
    if np.corrcoef(pc1, pt)[0, 1] < 0:
        pc1 = -pc1
    cells["pc1_siggenes"] = pc1
    r_pca = float(np.corrcoef(pc1, pt)[0, 1])
    print(f"PC1-on-signature-genes vs (iCAF-myCAF): r = {r_pca:.3f}", flush=True)
    win_comp = (cells.groupby(["window", "treatment"]).size()
                    .unstack(fill_value=0))
    print(win_comp)

    # gene trends: Spearman rho vs pseudotime rank
    gene_cols = [c for c in cells.columns if c not in
                 ("sample", "treatment", "iCAF", "myCAF", "pseudotime",
                  "pseudotime_rank", "window", "pc1_siggenes")]
    det = (cells[gene_cols].to_numpy() > 0).mean(0)
    keep_g = np.array(gene_cols)[det >= MIN_DET_FRAC]
    print(f"testing {len(keep_g)} genes (det >= {MIN_DET_FRAC})", flush=True)
    prank = cells["pseudotime_rank"].to_numpy().astype(float)
    rows = []
    X = cells[keep_g].to_numpy()
    for j, g in enumerate(keep_g):
        rho, p = spearmanr(X[:, j], prank)
        rows.append({"gene": g, "spearman_rho": float(rho),
                     "p_nominal": float(p),
                     "detection_frac": round(float((X[:, j] > 0).mean()), 3)})
    trends = pd.DataFrame(rows).sort_values("spearman_rho")
    m = len(trends)
    trends["p_fdr_bh"] = (trends["p_nominal"].rank() * 1.0)
    trends = trends.sort_values("p_nominal")
    trends["p_fdr_bh"] = (trends["p_nominal"] * m /
                           np.arange(1, m + 1)).clip(upper=1.0)
    trends["p_fdr_bh"] = trends["p_fdr_bh"][::-1].cummin()[::-1]
    trends = trends.sort_values("spearman_rho")
    trends.to_csv(os.path.join(OUT, "gene_trends.csv"), index=False)

    cells[["sample", "treatment", "pseudotime", "pseudotime_rank",
           "window", "iCAF", "myCAF", "pc1_siggenes"]].to_csv(
        os.path.join(OUT, "cell_pseudotime.csv"), index=False)

    top_up = trends.tail(15).iloc[::-1]    # rising toward iCAF
    top_dn = trends.head(15)               # falling (myCAF markers)
    print("\nTOP RISING along myCAF->iCAF pseudotime:")
    print(top_up[["gene", "spearman_rho", "p_nominal", "p_fdr_bh",
                  "detection_frac"]].to_string(index=False))
    print("\nTOP FALLING along myCAF->iCAF pseudotime:")
    print(top_dn[["gene", "spearman_rho", "p_nominal", "p_fdr_bh",
                  "detection_frac"]].to_string(index=False))

    # window means for top genes
    wm = cells.groupby("window")[[*top_up["gene"], *top_dn["gene"]]].mean()
    wm.to_csv(os.path.join(OUT, "window_means_top_genes.csv"))

    # figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    col = {"Veh": "#2980b9", "RMC6236": "#c0392b"}
    ax = axes[0]
    for trt in ("Veh", "RMC6236"):
        s = cells[cells["treatment"] == trt]
        ax.scatter(s["myCAF"], s["iCAF"], c=col[trt], label=trt, s=40,
                   alpha=0.8, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("myCAF score (mean fractional rank)")
    ax.set_ylabel("iCAF score (mean fractional rank)")
    ax.set_title("KPC2 fibroblasts: myCAF vs iCAF\n"
                 "pseudotime = iCAF − myCAF (95 cells)")
    ax.legend()
    ax = axes[1]
    win_means = cells.groupby("window", observed=False)[["myCAF", "iCAF"]].mean()
    win_means.plot(kind="bar", ax=ax, color=["#27ae60", "#e67e22"],
                   edgecolor="black")
    ax.set_title("Signature scores across pseudotime windows")
    ax.set_xlabel("pseudotime window (W1=myCAF-like → W4=iCAF-like)")
    ax.set_ylabel("mean score")
    ax.legend(title="")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_pseudotime_scores.png"), dpi=150)
    plt.close(fig)

    plot_genes = [*top_up["gene"].head(4), *top_dn["gene"].head(4)]
    fig, axes = plt.subplots(2, 4, figsize=(14, 7), sharex=True)
    wlab = [f"W{i+1}" for i in range(N_WINDOWS)]
    for ax, g in zip(axes.ravel(), plot_genes):
        means = cells.groupby("window", observed=False)[g].mean().reindex(wlab)
        ses = cells.groupby("window", observed=False)[g].sem().reindex(wlab)
        ax.errorbar(range(N_WINDOWS), means, yerr=ses, fmt="-o",
                    color="#2c3e50", ecolor="#95a5a6", capsize=3)
        rho = trends.set_index("gene").loc[g, "spearman_rho"]
        ax.set_title(f"{g}  (ρ={rho:.2f})", fontsize=10)
        ax.set_xticks(range(N_WINDOWS))
        ax.set_xticklabels(wlab)
        ax.set_ylabel("mean log1p(CP10k)")
    fig.suptitle("Top genes changing along myCAF→iCAF pseudotime "
                 "(mean ± SEM per window)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_top_gene_trends.png"), dpi=150)
    plt.close(fig)

    print("pca_vs_score_diff_r:", round(r_pca, 3))
    print("DONE ->", OUT)


if __name__ == "__main__":
    main()
