#!/usr/bin/env python3
"""ANALYSIS A: ligand-receptor search for the indirect-reprogramming hypothesis.

Question: which ligands change in the NON-fibroblast compartments (tumor/
epithelial, immune) under RMC-6236 vs vehicle (KPC2 primary; MRTX1133 arm
secondary, KPC1) in a way that could drive the fibroblast iCAF shift seen
in the main GSE337208 CAF mining?

Method:
  1. Reuse GSE337208 matrices (data/) with the SAME conventions as
     gse337208_caf_mining/run.py: QC (n_genes>=200, totals>=500, pct_mt<=15),
     CP10k + log1p, deterministic k-means on the 16-gene fibroblast/
     epithelial/immune marker panel. Cluster summaries are cross-checked
     against the sibling output cluster_annotation.csv.
  2. Compartment labels from the annotated clusters:
       fibroblast: fibroblast_cluster==True
       epithelial : non-fibroblast, epi_mean >= 1.0
       immune     : non-fibroblast, epi_mean < 1.0, imm_mean >= 0.7
  3. Per-sample pseudobulk counts for the epithelial and immune
     compartments (KPC2 Veh GSM9850800 vs KPC2 RMC6236 GSM9850799;
     secondary: KPC1 Veh GSM9850794 vs KPC1 MRTX1133 GSM9850795).
  4. For candidate iCAF-inducing ligand genes, compute log2 fold change
     (treated vs vehicle) of CP10k-normalized log2(x+1) expression,
     per compartment. n=1 sample/arm: descriptive only, no p-values.
  5. Check cognate receptor expression in fibroblasts from the sibling
     fibroblast_pseudobulk_counts.tsv (KPC2 pair).

Ligand -> receptor map:
  Il1a/Il1b -> Il1r1 (+Il1r2 decoy), Il6 -> Il6ra + Il6st(gp130),
  Lif -> Lifr + Il6st, Osm -> Osmr + Il6st, Cxcl1 -> Cxcr2,
  Csf2 -> Csf2ra + Csf2rb, Csf3 -> Csf3r, Tnf -> Tnfrsf1a + Tnfrsf1b.

Re-runnable: python3 run.py   (writes output/ + figures)
"""
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread

HERE = os.path.dirname(os.path.abspath(__file__))
SIBLING = os.path.abspath(os.path.join(HERE, "..", "gse337208_caf_mining"))
DATA = os.path.join(SIBLING, "data")
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)

# gsm -> (file prefix, model, treatment)
SAMPLES = {
    "GSM9850794": ("343P", "KPC1", "Veh"),
    "GSM9850795": ("MRTX1133343P", "KPC1", "MRTX1133"),
    "GSM9850799": ("Hy-RMC6236", "KPC2", "RMC6236"),
    "GSM9850800": ("Hy-Veh", "KPC2", "Veh"),
}

FIBRO_MARKERS = ["Col1a1", "Col1a2", "Col3a1", "Dcn", "Pdpn", "Lum", "Sparc"]
EPI_MARKERS = ["Epcam", "Krt19", "Krt8", "Krt18"]
IMM_MARKERS = ["Ptprc", "Cd3e", "Cd68", "Cd79a", "Itgam"]
PANEL = {"fib": FIBRO_MARKERS, "epi": EPI_MARKERS, "imm": IMM_MARKERS}

LIGANDS = {
    "Il1a": ["Il1r1", "Il1r2"],
    "Il1b": ["Il1r1", "Il1r2"],
    "Il6": ["Il6ra", "Il6st"],
    "Lif": ["Lifr", "Il6st"],
    "Osm": ["Osmr", "Il6st"],
    "Cxcl1": ["Cxcr2"],
    "Csf2": ["Csf2ra", "Csf2rb"],
    "Csf3": ["Csf3r"],
    "Tnf": ["Tnfrsf1a", "Tnfrsf1b"],
}

CONTRASTS = [
    ("KPC2", "RMC6236", "GSM9850799", "GSM9850800"),   # primary
    ("KPC1", "MRTX1133", "GSM9850795", "GSM9850794"),   # secondary
]


def kmeans(X, k, seed=0, iters=100):
    """Deterministic k-means (numpy only), identical to sibling run.py."""
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
    assert mat.shape[0] == len(genes) and mat.shape[1] == len(barcodes)
    return genes, mat, model, trt


def compartment_masks(genes, mat):
    """Return dict compartment -> boolean mask over QC'd cells."""
    totals = np.asarray(mat.sum(axis=0)).ravel()
    n_genes = np.asarray((mat > 0).sum(axis=0)).ravel()
    is_mt = np.array([g.startswith("mt-") for g in genes])
    mt_counts = np.asarray(mat[is_mt, :].sum(axis=0)).ravel()
    pct_mt = np.where(totals > 0, 100 * mt_counts / totals, 0)
    keep = (n_genes >= 200) & (totals >= 500) & (pct_mt <= 15)
    mat = mat[:, keep]
    totals = totals[keep]
    scale = 1e4 / totals
    logn = (mat @ sparse.diags(scale)).tocsr()
    logn.data = np.log1p(logn.data)
    gidx = {g: i for i, g in enumerate(genes)}
    panel_genes = [g for grp in PANEL.values() for g in grp if g in gidx]
    Xp = np.asarray(logn[[gidx[g] for g in panel_genes], :].toarray()).T
    Xz = (Xp - Xp.mean(0)) / np.maximum(Xp.std(0), 1e-8)
    lab = kmeans(Xz, 6, seed=0)
    prog_of = {}
    for grp, gl in PANEL.items():
        cols = [panel_genes.index(g) for g in gl if g in panel_genes]
        prog_of[grp] = Xp[:, cols].mean(1) if cols else np.zeros(Xp.shape[0])
    is_fib = lab == (lambda: max(range(6), key=lambda i: (
        (prog_of["fib"][lab == i].mean() if (lab == i).any() else -1e9) -
        max(prog_of["epi"][lab == i].mean() if (lab == i).any() else -1e9,
            prog_of["imm"][lab == i].mean() if (lab == i).any() else -1e9))))()
    is_fib = is_fib & ~(prog_of["epi"] > prog_of["fib"])
    is_fib = is_fib & ~(prog_of["imm"] > prog_of["fib"])
    return {
        "fibroblast": is_fib,
        "epithelial": (~is_fib) & (prog_of["epi"] > prog_of["fib"]) &
                      (prog_of["epi"] > prog_of["imm"]),
        "immune": (~is_fib) & (prog_of["imm"] > prog_of["fib"]) &
                  (prog_of["imm"] > prog_of["epi"]),
        "mat_qc": mat, "genes": genes,
    }


def main():
    # 1. sanity: compartment labels vs sibling cluster_annotation.csv
    annot = pd.read_csv(os.path.join(SIBLING, "output",
                                     "cluster_annotation.csv"))
    pb = {}  # (gsm, compartment) -> pseudobulk count Series
    for gsm in SAMPLES:
        print(f"--- {gsm} ({SAMPLES[gsm][1]} {SAMPLES[gsm][2]})", flush=True)
        genes, mat, model, trt = load_sample(gsm)
        d = compartment_masks(genes, mat)
        qmat = d["mat_qc"]
        for comp in ("fibroblast", "epithelial", "immune"):
            m = d[comp]
            s = pd.Series(np.asarray(qmat[:, m].sum(axis=1)).ravel(),
                          index=genes).groupby(level=0).sum()
            pb[(gsm, comp)] = s
            print(f"  {comp}: {int(m.sum())} cells", flush=True)
        # cross-check fibroblast count vs sibling annotation
        sib_nfib = annot[(annot["sample"] == gsm) &
                         (annot["fibroblast_cluster"])]["n_cells"]
        sib_nfib = int(sib_nfib.sum()) if len(sib_nfib) else 0
        print(f"  sibling fibroblast-cluster total: {sib_nfib}", flush=True)

    genes_all = sorted({g for s in pb.values() for g in s.index})
    pbm = pd.concat(pb, axis=1).fillna(0).astype(int)
    pbm.columns = pd.MultiIndex.from_tuples(pbm.columns,
                                            names=["sample", "compartment"])

    def cp10k(col):
        return 1e4 * col / col.sum()

    ligand_rows = []
    for model, drug, g_t, g_v in CONTRASTS:
        for comp in ("epithelial", "immune"):
            t, v = cp10k(pbm[(g_t, comp)]), cp10k(pbm[(g_v, comp)])
            nt, nv = int(pbm[(g_t, comp)].sum()), int(pbm[(g_v, comp)].sum())
            for lig, recs in LIGANDS.items():
                xt = float(t.get(lig, 0.0))
                xv = float(v.get(lig, 0.0))
                l2fc = np.log2(xt + 1) - np.log2(xv + 1)
                ligand_rows.append({
                    "contrast": f"{model}_{drug}_vs_Veh",
                    "compartment": comp,
                    "ligand": lig, "receptors": ",".join(recs),
                    "cp10k_treated": round(xt, 2),
                    "cp10k_vehicle": round(xv, 2),
                    "log2FC_treated_vs_vehicle": round(float(l2fc), 3),
                    "direction": ("up" if l2fc > 0.5 else
                                  "down" if l2fc < -0.5 else "flat"),
                    "n_umi_treated": nt, "n_umi_vehicle": nv,
                })
    lig_df = pd.DataFrame(ligand_rows)
    lig_df.to_csv(os.path.join(OUT, "ligand_receptor_degs.csv"), index=False)
    print(lig_df.to_string(index=False))

    # 2. receptor expression in fibroblasts (KPC2 pair, sibling pseudobulk)
    fib_pb = pd.read_csv(
        os.path.join(SIBLING, "output", "fibroblast_pseudobulk_counts.tsv"),
        sep="\t").set_index("gene")
    rec_rows = []
    for lig, recs in LIGANDS.items():
        for rec in recs:
            row = {"ligand": lig, "receptor": rec}
            for gsm in ("GSM9850800", "GSM9850799"):
                if rec in fib_pb.index:
                    c = 1e4 * float(fib_pb.loc[rec, gsm]) / fib_pb[gsm].sum()
                    row[f"{gsm}_cp10k"] = round(c, 2)
                else:
                    row[f"{gsm}_cp10k"] = np.nan
            rec_rows.append(row)
    rec_df = pd.DataFrame(rec_rows)
    rec_df.to_csv(os.path.join(OUT, "receptor_expression_fibroblasts.csv"),
                  index=False)
    print(rec_df.to_string(index=False))

    # 3. figures
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pri = lig_df[lig_df["contrast"] == "KPC2_RMC6236_vs_Veh"]
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    for ax, comp in zip(axes, ("epithelial", "immune")):
        sub = pri[pri["compartment"] == comp]
        colors = ["#c0392b" if d == "up" else "#2980b9" if d == "down"
                  else "#7f8c8d" for d in sub["direction"]]
        ax.barh(sub["ligand"], sub["log2FC_treated_vs_vehicle"],
                color=colors, edgecolor="black")
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title(f"KPC2 RMC6236 vs Veh — {comp} ligands\n"
                     f"(log2FC of CP10k-normalized log counts, n=1/arm)")
        ax.set_xlabel("log2FC (RMC6236 vs vehicle)")
    axes[0].set_ylabel("ligand gene")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_ligand_log2fc.png"), dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    melt = rec_df.melt(id_vars=["ligand", "receptor"],
                       var_name="sample", value_name="cp10k")
    melt["arm"] = melt["sample"].map({"GSM9850800_cp10k": "Veh",
                                      "GSM9850799_cp10k": "RMC6236"})
    wide = melt.pivot_table(index=["ligand", "receptor"], columns="arm",
                            values="cp10k").fillna(0)
    im = ax.imshow(np.log10(wide.values + 1), cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(wide.columns)))
    ax.set_xticklabels(wide.columns)
    ax.set_yticks(range(len(wide)))
    ax.set_yticklabels([f"{l}:{r}" for l, r in wide.index], fontsize=8)
    ax.set_title("KPC2 fibroblast receptor expression\n"
                 "log10(CP10k + 1), fibroblast pseudobulk")
    fig.colorbar(im, ax=ax, label="log10(CP10k+1)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_receptor_expression.png"), dpi=150)
    plt.close(fig)

    print("DONE ->", OUT)


if __name__ == "__main__":
    main()
