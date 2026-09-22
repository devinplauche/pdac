#!/usr/bin/env python3
"""Map CAF subtype signatures on one untreated human PDAC Visium section.

Dataset: GEO GSE274103 ("Spatial transcriptomics on treatment-naive pancreatic
ductal adenocarcinoma (PDAC) patients"), sample GSM8443449 (PDAC-p1),
10x Visium on FFPE, processed matrix + spatial files downloaded from GEO.

Scoring: repo signature compendium (analyses/signature_panel/signatures_v1.yaml)
scored with analyses/signature_panel/score_signatures.py, --organism human,
--method rank (AUCell-like; good for sparse spot data).

Pipeline:
  1. Load 10x h5 (genes x spots), keep in-tissue spots.
  2. Score all signatures in-process with score_signatures.score_rank on the
     FULL gene matrix (fractional ranks are computed within each spot).
  3. Define tumor spots from an epithelial/tumor rank-score
     (EPCAM/KRT7/KRT8/KRT18/MUC1/CEACAM5/CEACAM6/SOX9, present genes only).
  4. For every spot, compute distance (spot units) to the nearest tumor spot.
  5. Compare myCAF / iCAF / apCAF scores in distance bands
     (tumor / adjacent<=2 / mid 2-5 / distal>5) with Wilcoxon rank-sum and
     Spearman correlation of score vs distance.
  6. Write spot scores + annotations CSV and 3 spatial PNGs.

Re-runnable: python3 run.py  (expects data/ files; writes output/).
"""
import os
import sys

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.sparse import csc_matrix
from scipy.stats import mannwhitneyu, spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)

H5 = os.path.join(DATA, "GSM8443449_PDAC-p1_filtered_feature_bc_matrix.h5")
POS = os.path.join(DATA, "C30_spatial", "tissue_positions.csv")
SIG_YAML = os.path.join(HERE, "..", "signature_panel", "signatures_v1.yaml")
SCORE_PY = os.path.join(HERE, "..", "signature_panel", "score_signatures.py")

sys.path.insert(0, os.path.dirname(SCORE_PY))
import score_signatures as ss  # noqa: E402

TUMOR_GENES = ["EPCAM", "KRT7", "KRT8", "KRT18", "MUC1",
               "CEACAM5", "CEACAM6", "SOX9"]
NIL = object()


def main():
    # ---- 1. load matrix ---------------------------------------------------
    with h5py.File(H5) as f:
        m = f["matrix"]
        names = m["features"]["name"][:].astype(str)
        barcodes = m["barcodes"][:].astype(str)
        shape = tuple(m["shape"][:])
        X = csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]),
                       shape=shape)  # genes x spots (CSC)
    pos = pd.read_csv(POS, index_col=0)
    keep = pos.index[pos["in_tissue"] == 1]
    keep = [b for b in barcodes if b in set(keep)]
    cols = np.array([np.flatnonzero(barcodes == b)[0] for b in keep])
    X = X[:, cols]
    pos = pos.loc[keep]
    genes = list(names)
    n_genes, n_spots = X.shape
    print(f"matrix: {n_genes} genes x {n_spots} tissue spots; "
          f"median UMIs/spot={float(np.median(np.asarray(X.sum(axis=0)).ravel())):.0f}")

    # ---- 2. rank scores (full-matrix ranks, in-process) -------------------
    dense = X.toarray().astype(np.float32)  # ~89M entries ~ 357 MB
    expr = pd.DataFrame(dense, index=genes, columns=keep)
    del dense
    sigs, version = ss.load_signatures(SIG_YAML, "human")
    print(f"compendium v{version}")
    scores, coverage = {}, {}
    for name, s in sigs.items():
        sc, n = ss.score_rank(expr, s["genes"])
        scores[name] = sc
        coverage[name] = (s["label"], len(s["genes"]), n)
        missing = sorted(set(s["genes"]) - set(genes))
        print(f"  {name:18s} {n:2d}/{len(s['genes'])} genes"
              + (f"  missing: {','.join(missing)}" if missing else ""))
    scores = pd.DataFrame(scores)

    # tumor/epithelial rank score from genes present in the matrix
    tumor_genes = [g for g in TUMOR_GENES if g in genes]
    tscore, tn = ss.score_rank(expr, tumor_genes)
    print(f"tumor score from {tn} genes: {tumor_genes}")

    # ---- 3. tumor spots ---------------------------------------------------
    # threshold: rank score above the 75th percentile (check map visually)
    thr = tscore.quantile(0.75)
    is_tumor = (tscore >= thr)
    print(f"tumor-spot threshold (75th pct) = {thr:.3f}; "
          f"n tumor spots = {is_tumor.sum()} / {n_spots}")

    # ---- 4. distance from tumor -------------------------------------------
    rc = pos[["array_row", "array_col"]].to_numpy(dtype=float)
    tumor_rc = rc[is_tumor.to_numpy()]
    # brute-force nearest-tumor distance in spot units
    dmin = np.full(n_spots, np.inf)
    CH = 500
    for i in range(0, n_spots, CH):
        d = np.sqrt(((rc[i:i + CH, None, :] - tumor_rc[None, :, :]) ** 2).sum(-1))
        dmin[i:i + CH] = d.min(axis=1)
    band = np.where(is_tumor.to_numpy(), "tumor",
             np.where(dmin <= 2, "adjacent", np.where(dmin <= 5, "mid", "distal")))
    print("band counts:", pd.Series(band).value_counts().to_dict())

    # ---- 5. stats ---------------------------------------------------------
    df = pd.DataFrame({
        "barcode": keep, "tumor_score": tscore.values,
        "is_tumor": is_tumor.to_numpy(), "dist_to_tumor": dmin, "band": band,
        **{c: scores[c].values for c in scores.columns}})
    res = []
    for sig in ["myCAF", "iCAF", "apCAF"]:
        a = df.loc[df["band"] == "adjacent", sig]
        d = df.loc[df["band"] == "distal", sig]
        u = mannwhitneyu(a, d, alternative="two-sided")
        rho, p = spearmanr(df.loc[df["band"] != "tumor", "dist_to_tumor"],
                           df.loc[df["band"] != "tumor", sig])
        med = df.groupby("band")[sig].median()
        res.append({"signature": sig,
                    "median_tumor": med.get("tumor", np.nan),
                    "median_adjacent": med.get("adjacent", np.nan),
                    "median_mid": med.get("mid", np.nan),
                    "median_distal": med.get("distal", np.nan),
                    "p_adjacent_vs_distal": u.pvalue,
                    "spearman_dist_rho": rho, "spearman_p": p})
        print(f"{sig}: medians tumor/adj/mid/dist = "
              f"{med.get('tumor',np.nan):.3f}/{med.get('adjacent',np.nan):.3f}/"
              f"{med.get('mid',np.nan):.3f}/{med.get('distal',np.nan):.3f}; "
              f"adj vs distal Wilcoxon p={u.pvalue:.2e}; "
              f"rho(score,dist)={rho:.3f} p={p:.2e}")
    pd.DataFrame(res).to_csv(os.path.join(OUT, "band_stats.csv"), index=False)
    df.to_csv(os.path.join(OUT, "spot_scores.csv"), index=False)
    print("wrote spot_scores.csv (%d rows)" % len(df))

    # ---- 6. figures -------------------------------------------------------
    xs, ys = pos["pxl_col_in_fullres"].to_numpy(), pos["pxl_row_in_fullres"].to_numpy()

    def scatter(vals, title, fname, tumor_outline=True, vmin=None, vmax=None,
                cmap="magma"):
        fig, ax = plt.subplots(figsize=(8, 7))
        if tumor_outline:
            t = is_tumor.to_numpy()
            ax.scatter(xs[~t], ys[~t], c="lightgray", s=6, lw=0)
        else:
            ax.scatter(xs, ys, c="lightgray", s=6, lw=0)
        sc = ax.scatter(xs, ys, c=vals, s=8, cmap=cmap, vmin=vmin, vmax=vmax,
                        lw=0)
        ax.set_aspect("equal"); ax.invert_yaxis()
        ax.set_title(title, fontsize=11)
        plt.colorbar(sc, ax=ax, shrink=0.8, label="rank score")
        ax.set_xticks([]); ax.set_yticks([])
        fig.tight_layout(); fig.savefig(os.path.join(OUT, fname), dpi=150)
        plt.close(fig)

    scatter(df["is_tumor"].astype(float),
            "Tumor spots (EPCAM/KRT7/KRT8/KRT18/MUC1/CEACAM5/CEACAM6/SOX9, top quartile)\n"
            "GSM8443449 PDAC-p1 — untreated human PDAC, GSE274103",
            "tumor_map.png", tumor_outline=False, vmin=0, vmax=1, cmap="Reds")
    scatter(df["myCAF"], "myCAF score (Elyada 2019 panel, rank/AUCell-like)\n"
            "GSM8443449 PDAC-p1 — untreated human PDAC", "mycaf_map.png")
    scatter(df["iCAF"], "iCAF score (Elyada 2019 panel, rank/AUCell-like)\n"
            "GSM8443449 PDAC-p1 — untreated human PDAC", "icaf_map.png")
    print("figures written to", OUT)


if __name__ == "__main__":
    main()
