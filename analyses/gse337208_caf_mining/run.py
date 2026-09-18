#!/usr/bin/env python3
"""Mine GSE337208 (murine KPC PDAC scRNA-seq, RMC-6236 / MRTX1133 +/- checkpoint
blockade) for CAF subtype shifts under RMC-6236.

Pipeline:
  1. Load per-sample CellRanger count matrices (genes x cells, sparse).
  2. Basic QC (n_genes, total counts, % mitochondrial).
  3. Log-normalize (CP10k, log1p); identify fibroblasts by pan-fibroblast
     marker score (Col1a1/Col1a2/Col3a1/Dcn/Pdpn/Lum), validated against
     epithelial (Epcam/Krt19) and immune (Ptprc) markers.
  4. Pseudobulk fibroblast profiles per sample -> repo scorer
     (score_signatures.py, zscore method, mouse lists) -> descriptive
     within-model contrasts (n=1 sample/arm: no sample-level p-values).
  5. Cell-level rank scores for the primary KPC2 RMC6236-vs-Veh pair
     (memory-efficient fractional ranks), Mann-Whitney U across cells
     (labeled as cell-level / pseudoreplicated), figures.

Re-runnable: python3 run.py  (expects data/*.tsv.gz + *.mtx.gz from GEO)
"""
import gzip
import os
import subprocess
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from scipy.stats import mannwhitneyu

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)

SCORER = os.path.abspath(os.path.join(HERE, "..", "signature_panel",
                                       "score_signatures.py"))
SIGYAML = os.path.abspath(os.path.join(HERE, "..", "signature_panel",
                                        "signatures_v1.yaml"))

# gsm -> (file prefix, KPC model, treatment label)
SAMPLES = {
    "GSM9850794": ("343P", "KPC1", "Veh"),
    "GSM9850795": ("MRTX1133343P", "KPC1", "MRTX1133"),
    "GSM9850796": ("KPC-T-CTLA4", "KPC3", "anti-CTLA4"),
    "GSM9850797": ("KPC-T-MRTX", "KPC3", "MRTX1133"),
    "GSM9850798": ("KPC-T-1133-CTLA4", "KPC3", "MRTX1133+anti-CTLA4"),
    "GSM9850799": ("Hy-RMC6236", "KPC2", "RMC6236"),
    "GSM9850800": ("Hy-Veh", "KPC2", "Veh"),
    "GSM9850801": ("Hy-RMC6236-aCTLA4", "KPC2", "RMC6236+anti-CTLA4"),
    "GSM9850802": ("Hy-RMC6236-aPD1", "KPC2", "RMC6236+anti-PD1"),
}

FIBRO_MARKERS = ["Col1a1", "Col1a2", "Col3a1", "Dcn", "Pdpn", "Lum",
                 "Sparc"]
EPI_MARKERS = ["Epcam", "Krt19", "Krt8", "Krt18"]
IMM_MARKERS = ["Ptprc", "Cd3e", "Cd68", "Cd79a", "Itgam"]
PANEL = {"fib": FIBRO_MARKERS, "epi": EPI_MARKERS, "imm": IMM_MARKERS}


def kmeans(X, k, seed=0, iters=100):
    """Deterministic k-means (numpy only). X: n x p."""
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
    mat = mmread(base + "_matrix.mtx.gz").tocsr()  # genes x cells
    assert mat.shape[0] == len(genes) and mat.shape[1] == len(barcodes)
    return genes, barcodes, mat, model, trt


def main():
    qc_rows, pb_counts, pb_meta = [], {}, []
    cluster_table = []
    cell_store = {}  # gsm -> dict for the KPC2 pair (cell-level analysis)

    for gsm, (pfx, model, trt) in SAMPLES.items():
        print(f"--- {gsm} ({model} {trt})", flush=True)
        genes, barcodes, mat, model, trt = load_sample(gsm)
        n_cells = mat.shape[1]
        totals = np.asarray(mat.sum(axis=0)).ravel()
        n_genes = np.asarray((mat > 0).sum(axis=0)).ravel()
        is_mt = np.array([g.startswith("mt-") for g in genes])
        mt_counts = np.asarray(mat[is_mt, :].sum(axis=0)).ravel()
        pct_mt = np.where(totals > 0, 100 * mt_counts / totals, 0)

        keep = (n_genes >= 200) & (totals >= 500) & (pct_mt <= 15)
        mat = mat[:, keep]
        totals, n_genes, pct_mt = totals[keep], n_genes[keep], pct_mt[keep]
        print(f"  cells: {n_cells} -> {keep.sum()} after QC", flush=True)

        # log-normalize
        scale = 1e4 / totals
        logn = (mat @ sparse.diags(scale)).tocsr()
        logn.data = np.log1p(logn.data)

        gidx = {g: i for i, g in enumerate(genes)}

        # --- fibroblast identification: k-means on marker panel ---
        panel_genes = [g for grp in PANEL.values() for g in grp if g in gidx]
        Xp = np.asarray(logn[[gidx[g] for g in panel_genes], :].toarray()).T
        Xz = (Xp - Xp.mean(0)) / np.maximum(Xp.std(0), 1e-8)
        lab = kmeans(Xz, 6, seed=0)
        prog_of = {}
        for grp, gl in PANEL.items():
            cols = [panel_genes.index(g) for g in gl if g in panel_genes]
            prog_of[grp] = Xp[:, cols].mean(1) if cols else np.zeros(Xp.shape[0])
        cluster_rows = []
        best, best_margin = -1, -1e9
        for i in range(6):
            m = lab == i
            fm = float(prog_of["fib"][m].mean()) if m.any() else 0.0
            em = float(prog_of["epi"][m].mean()) if m.any() else 0.0
            im = float(prog_of["imm"][m].mean()) if m.any() else 0.0
            margin = fm - max(em, im)
            cluster_rows.append({"sample": gsm, "cluster": i,
                                 "n_cells": int(m.sum()),
                                 "fib_mean": round(fm, 3),
                                 "epi_mean": round(em, 3),
                                 "imm_mean": round(im, 3)})
            if margin > best_margin:
                best_margin, best = margin, i
        bm = cluster_rows[best]
        has_fib = bm["fib_mean"] > 1.0 and best_margin > 1.0
        if not has_fib:
            # No fibroblast population captured in this sample (e.g. the
            # RMC6236+checkpoint combo arms). Record and skip downstream.
            print(f"  WARNING: no clean fibroblast cluster in {gsm}; "
                  f"best cluster {best}: {bm}", flush=True)
            for r in cluster_rows:
                r["fibroblast_cluster"] = False
                r["note"] = "no_fibroblast_cluster"
            cluster_table.extend(cluster_rows)
            qc_rows.append({"sample": gsm, "model": model, "treatment": trt,
                            "cells_raw": int(n_cells),
                            "cells_qc": int(keep.sum()),
                            "fibroblasts": 0,
                            "fibroblast_pct": 0.0,
                            "median_umi_qc": round(float(np.median(totals)), 1),
                            "median_genes_qc": round(float(np.median(n_genes)), 1)})
            pb_meta.append((gsm, model, trt, 0))
            continue
        for r in cluster_rows:
            r["fibroblast_cluster"] = (r["cluster"] == best)
            r["note"] = ""
        # doublet cleanup: drop fibroblast-cluster cells dominated by
        # epithelial or immune programs
        is_fib = lab == best
        is_fib = is_fib & ~(prog_of["epi"] > prog_of["fib"])
        is_fib = is_fib & ~(prog_of["imm"] > prog_of["fib"])
        n_fib = int(is_fib.sum())
        print(f"  fibroblast cluster {best}: {bm['n_cells']} cells "
              f"(fib={bm['fib_mean']}, epi={bm['epi_mean']}, "
              f"imm={bm['imm_mean']}) -> {n_fib} after doublet cleanup",
              flush=True)
        cluster_table.extend(cluster_rows)

        qc_rows.append({"sample": gsm, "model": model, "treatment": trt,
                        "cells_raw": int(n_cells),
                        "cells_qc": int(keep.sum()),
                        "fibroblasts": n_fib,
                        "fibroblast_pct": round(100 * n_fib / keep.sum(), 2),
                        "median_umi_qc": round(float(np.median(totals)), 1),
                        "median_genes_qc": round(float(np.median(n_genes)), 1)})
        # pseudobulk fibroblast counts
        pb_counts[gsm] = np.asarray(mat[:, is_fib].sum(axis=1)).ravel()
        pb_meta.append((gsm, model, trt, n_fib))
        if model == "KPC2" and trt in ("Veh", "RMC6236"):
            cell_store[gsm] = {"genes": genes, "logn": logn, "is_fib": is_fib,
                               "trt": trt}

    qc = pd.DataFrame(qc_rows)
    qc.to_csv(os.path.join(OUT, "qc_table.csv"), index=False)
    print(qc.to_string(index=False))
    pd.DataFrame(cluster_table).to_csv(
        os.path.join(OUT, "cluster_annotation.csv"), index=False)

    # ---- pseudobulk matrix -> repo scorer (zscore across samples) ----
    # NOTE: samples were processed with different CellRanger references
    # (feature files differ in length), so align on gene-symbol union.
    pb_series = {}
    for gsm, model, trt, nfib in pb_meta:
        if nfib == 0:
            print(f"  skipping pseudobulk for {gsm} (no fibroblasts)",
                  flush=True)
            continue
        pfx = SAMPLES[gsm][0]
        with gzip.open(os.path.join(
                DATA, f"{gsm}_{pfx}_features.tsv.gz"), "rt") as f:
            feats = pd.read_csv(f, sep="\t", header=None,
                                names=["ensembl", "symbol", "type"])
        s = pd.Series(pb_counts[gsm], index=feats["symbol"].values)
        s = s.groupby(level=0).sum()  # collapse duplicated symbols
        pb_series[gsm] = s
    pb = pd.concat(pb_series, axis=1).fillna(0).astype(int)
    pb.index.name = "gene"
    pb_path = os.path.join(OUT, "fibroblast_pseudobulk_counts.tsv")
    pb.reset_index().rename(columns={"index": "gene"}).to_csv(
        pb_path, sep="\t", index=False)
    scores_path = os.path.join(OUT, "pseudobulk_signature_scores.csv")
    cov_path = os.path.join(OUT, "gene_coverage.csv")
    r = subprocess.run(
        [sys.executable, SCORER, "--input", pb_path, "--signatures", SIGYAML,
         "--organism", "mouse", "--method", "zscore",
         "--out", scores_path, "--coverage-out", cov_path],
        capture_output=True, text=True)
    print(r.stdout)
    if r.returncode != 0:
        print(r.stderr, file=sys.stderr)
        sys.exit(1)
    scores = pd.read_csv(scores_path, index_col="sample")
    meta = pd.DataFrame(pb_meta, columns=["sample", "model", "treatment",
                                          "n_fibroblasts"]).set_index("sample")
    (scores.join(meta)).to_csv(os.path.join(OUT,
                                            "pseudobulk_scores_with_meta.csv"))

    # detection rates of signature genes in fibroblasts (dropout report)
    import yaml
    with open(SIGYAML) as f:
        comp = yaml.safe_load(f)["signatures"]
    det_rows = []
    union_genes = set(pb.index)
    for gsm, model, trt, nfib in pb_meta:
        for signame, s in comp.items():
            for g in s["mouse"]:
                det_rows.append({"sample": gsm, "signature": signame,
                                 "gene": g,
                                 "detected": g in union_genes})
    det = pd.DataFrame(det_rows)
    det.to_csv(os.path.join(OUT, "signature_gene_presence.csv"), index=False)
    # per-gene detection fraction across fibroblast cells (KPC2 pair)
    print("coverage table written")

    # ---- cell-level rank scores for KPC2 Veh vs RMC6236 ----
    # memory-efficient fractional ranks: only signature genes needed
    sig_genes = sorted({g for s in comp.values() for g in s["mouse"]})
    cell_scores = []
    for gsm, d in cell_store.items():
        genes, logn, is_fib, trt = (d["genes"], d["logn"], d["is_fib"],
                                    d["trt"])
        csc = logn.tocsc()
        n_genes_total = csc.shape[0]
        gidx = {g: i for i, g in enumerate(genes)}
        sig_idx = {g: gidx[g] for g in sig_genes if g in gidx}
        fib_cols = np.where(is_fib)[0]
        print(f"  ranking {len(fib_cols)} fibroblasts x "
              f"{len(sig_idx)} sig genes ({gsm})", flush=True)
        # per cell: fractional rank of each signature gene among all genes
        for j in fib_cols:
            col = csc.getcol(j)
            vals = col.data
            nnz = len(vals)
            n_zero = n_genes_total - nnz
            srt = np.sort(vals)
            row = {"sample": gsm, "treatment": trt, "cell": j}
            for g, gi in sig_idx.items():
                x = col[gi, 0]
                if x == 0:
                    # all zeros tie: mean rank of the zero block
                    rank = (n_zero + 1) / 2 / n_genes_total
                else:
                    lt = np.searchsorted(srt, x, side="left")
                    le = np.searchsorted(srt, x, side="right")
                    rank = (n_zero + (lt + le) / 2 + 1) / n_genes_total
                row[g] = rank
            cell_scores.append(row)
    cells = pd.DataFrame(cell_scores)
    # signature score = mean fractional rank of member genes
    for signame, s in comp.items():
        cols = [g for g in s["mouse"] if g in cells.columns]
        cells[signame] = cells[cols].mean(axis=1) if cols else np.nan
        print(f"{signame}: {len(cols)}/{len(s['mouse'])} genes in matrix")
    cells.to_csv(os.path.join(OUT, "kpc2_cell_signature_scores.csv"),
                 index=False)

    # MWU tests Veh vs RMC6236 (cell-level; pseudoreplicated - labeled as such)
    stat_rows = []
    for signame in comp:
        a = cells.loc[cells.treatment == "Veh", signame].dropna()
        b = cells.loc[cells.treatment == "RMC6236", signame].dropna()
        u, p = mannwhitneyu(a, b, alternative="two-sided")
        stat_rows.append({"signature": signame,
                          "n_cells_Veh": len(a), "n_cells_RMC6236": len(b),
                          "median_Veh": round(float(a.median()), 4),
                          "median_RMC6236": round(float(b.median()), 4),
                          "median_diff_RMC_minus_Veh":
                              round(float(b.median() - a.median()), 4),
                          "MWU_p_cell_level": float(p)})
    stats = pd.DataFrame(stat_rows)
    stats.to_csv(os.path.join(OUT, "kpc2_contrast_stats.csv"), index=False)
    print(stats.to_string(index=False))

    # detection fraction per signature gene in KPC2 fibroblasts
    detf = []
    for gsm, d in cell_store.items():
        genes, logn, is_fib = d["genes"], d["logn"], d["is_fib"]
        gidx = {g: i for i, g in enumerate(genes)}
        csc = logn.tocsc()[:, np.where(is_fib)[0]]
        for g in sig_genes:
            if g in gidx:
                col = csc.getrow(gidx[g])
                frac = (col.nnz / csc.shape[1])
                detf.append({"sample": gsm, "gene": g,
                             "fibroblast_detection_frac": round(frac, 3)})
    pd.DataFrame(detf).to_csv(os.path.join(OUT,
                                           "signature_gene_detection.csv"),
                              index=False)
    print("ALL DONE")


if __name__ == "__main__":
    main()
