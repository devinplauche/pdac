#!/usr/bin/env python3
"""
TCGA-PAAD CAF subtype signature scoring + clinical correlates.

Pipeline (re-runnable):
  1. Extract GDC STAR-Counts tarball from data/ and build a gene x sample
     counts matrix (gene_name, unstranded counts).
  2. Score the repo signature panel AS-IS (signatures_v1.yaml, human,
     zscore) via analyses/signature_panel/score_signatures.py.
  3. Naive immune/stromal scores (mean z of marker genes; crude by design).
  4. Merge with GDC clinical (OS + AJCC pathologic stage).
  5. Statistics: KM + manual log-rank (median split) for iCAF/myCAF;
     Kruskal-Wallis for stage; Spearman/Pearson vs immune score.
  6. CSVs + PNG figures in output/.

Writes large intermediates only to data/ (gitignored); output/ stays small.
"""
import json
import os
import subprocess
import sys
import tarfile

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "output")
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PANEL = os.path.join(REPO, "analyses", "signature_panel")
os.makedirs(OUT, exist_ok=True)

TAR = os.path.join(DATA, "star_counts.tar")
FILEMAP = os.path.join(DATA, "file_ids.json")          # gdc file_id list
CLIN_JSON = os.path.join(DATA, "gdc_clinical.json")     # gdc /cases dump
MATRIX = os.path.join(DATA, "counts_matrix.tsv.gz")    # built gene x sample
SCORES = os.path.join(OUT, "caf_scores.csv")
COVERAGE = os.path.join(OUT, "signature_coverage.csv")

IMMUNE_GENES = ["PTPRC", "CD3E", "CD8A", "CD4", "MS4A1", "CD68", "NKG7", "GZMB"]
STROMAL_GENES = ["COL1A1", "COL3A1", "DCN", "PDPN", "FAP", "VIM"]


# ---------------------------------------------------------------- matrix
def manifest_map():
    """requested file_id -> tar member uuid (from MANIFEST.txt, stream-safe)."""
    mp = {}
    with tarfile.open(TAR, "r|gz") as tf:
        for m in tf:
            if m.name == "MANIFEST.txt":
                for line in tf.extractfile(m).read().decode().splitlines()[1:]:
                    parts = line.split("\t")
                    if len(parts) >= 2:
                        tar_uuid = parts[1].split("/")[1].split(".")[0]
                        mp[parts[0]] = tar_uuid
                break
    return mp


def ensure_raw_files(mmap):
    """Extract recoverable TSVs; download any missing file_ids individually."""
    os.makedirs(os.path.join(DATA, "raw"), exist_ok=True)
    have = {f.split(".")[0] for f in os.listdir(os.path.join(DATA, "raw"))}
    # extract from tar (stream mode; tolerates truncation)
    try:
        with tarfile.open(TAR, "r|gz") as tf:
            for m in tf:
                if m.name.endswith(".tsv"):
                    uuid = os.path.basename(m.name).split(".")[0]
                    if uuid not in have:
                        fh = tf.extractfile(m)
                        with open(os.path.join(DATA, "raw", os.path.basename(m.name)), "wb") as o:
                            o.write(fh.read())
                        have.add(uuid)
    except Exception as e:
        print("tar stream ended early (expected):", e)
    # which requested ids are still missing?
    inv = {v: k for k, v in mmap.items()}
    missing = [inv[u] for u in set(mmap.values()) - have if u in inv]
    print(f"have {len(have)} files, downloading {len(missing)} missing individually")
    for fid in missing:
        out = os.path.join(DATA, "raw", fid + ".rna_seq.augmented_star_gene_counts.tsv")
        r = subprocess.run(["curl", "-s", "--max-time", "300", "--retry", "3",
                            f"https://api.gdc.cancer.gov/data/{fid}", "-o", out],
                           capture_output=True)
        if os.path.getsize(out) < 1000:
            print("FAILED download:", fid, open(out).read()[:200])
    return inv


def build_matrix():
    if os.path.exists(MATRIX):
        print("matrix exists, skipping rebuild")
        return pd.read_csv(MATRIX, sep="\t", index_col=0)
    mmap = manifest_map()
    print("manifest entries:", len(mmap))
    inv = ensure_raw_files(mmap)   # tar_uuid -> requested file_id
    meta = {h["file_id"]: (h["cases"][0]["submitter_id"],
                           h["cases"][0]["samples"][0]["sample_type"])
            for h in json.load(open(os.path.join(DATA, "gdc_star_files.json")))["data"]["hits"]}
    frames, skipped = [], []
    for f in sorted(os.listdir(os.path.join(DATA, "raw"))):
        uuid = f.split(".")[0]
        fid = inv.get(uuid, uuid)
        if fid not in meta:
            skipped.append(f); continue
        try:
            df = pd.read_csv(os.path.join(DATA, "raw", f), sep="\t", comment="#")
            s = df.set_index("gene_name")["unstranded"]
        except Exception as e:
            print("skipping unreadable", f, e); skipped.append(f); continue
        s.name = fid
        frames.append(s)
    print(f"loaded {len(frames)} files, skipped {len(skipped)}: {skipped[:5]}")
    mat = pd.concat(frames, axis=1)
    mat.index.name = "gene"
    mat = mat.groupby(level=0).mean()   # collapse duplicated symbols
    mat.to_csv(MATRIX, sep="\t", compression="gzip")
    print("matrix:", mat.shape)
    # persist the fid -> (case, sample_type) map for map_files_to_samples
    json.dump({fid: {"case_id": c, "sample_type": st} for fid, (c, st) in meta.items()
               if fid in {fr.name for fr in frames}},
              open(os.path.join(DATA, "fid_map.json"), "w"))
    return mat


def map_files_to_samples(mat):
    """file uuid -> (case submitter_id, sample_type) using persisted fid map."""
    mp = pd.DataFrame.from_dict(json.load(open(os.path.join(DATA, "fid_map.json"))),
                                orient="index")
    cols = [c for c in mat.columns if c in mp.index]
    print(f"mapped {len(cols)}/{mat.shape[1]} files to cases")
    return mp.loc[cols]


# ---------------------------------------------------------------- scoring
def run_panel(mat):
    mat.to_csv(os.path.join(DATA, "matrix_for_scoring.tsv"), sep="\t")
    cmd = [sys.executable, os.path.join(PANEL, "score_signatures.py"),
           "--input", os.path.join(DATA, "matrix_for_scoring.tsv"),
           "--signatures", os.path.join(PANEL, "signatures_v1.yaml"),
           "--organism", "human", "--method", "zscore",
           "--out", SCORES, "--coverage-out", COVERAGE]
    print(" ".join(cmd))
    subprocess.run(cmd, check=True)
    return pd.read_csv(SCORES, index_col=0)


def zscore_df(mat):
    logm = np.log2(mat + 1)
    mu = logm.mean(axis=1)
    sd = logm.std(axis=1, ddof=1).replace(0, np.nan)
    return logm.sub(mu, axis=0).div(sd, axis=0)


def naive_score(z, genes, label):
    present = [g for g in genes if g in z.index]
    missing = [g for g in genes if g not in z.index]
    print(f"{label}: {len(present)}/{len(genes)} genes present"
          + (f" (missing: {missing})" if missing else ""))
    return z.loc[present].mean(axis=0) if present else pd.Series(np.nan, index=z.columns)


# ---------------------------------------------------------------- clinical
def load_clinical():
    d = json.load(open(CLIN_JSON))["data"]["hits"]
    rows = []
    for h in d:
        dg = h["diagnoses"][0]
        demo = h["demographic"]
        dead = demo.get("vital_status") == "Dead"
        t = demo.get("days_to_death") if dead else dg.get("days_to_last_follow_up")
        rows.append({"case_id": h["submitter_id"],
                     "vital_status": demo.get("vital_status"),
                     "os_event": int(dead),
                     "os_days": t,
                     "stage": dg.get("ajcc_pathologic_stage"),
                     "histology": dg.get("primary_diagnosis"),
                     "age": dg.get("age_at_diagnosis")})
    return pd.DataFrame(rows).set_index("case_id")


# ---------------------------------------------------------------- stats
def km_curve(times, events):
    """Kaplan-Meier estimate; returns (t_grid, surv)."""
    order = np.argsort(times)
    t, e = times[order], events[order]
    surv, grid = [1.0], [0.0]
    s = 1.0
    for ti in np.unique(t[e == 1]):
        d = np.sum((t == ti) & (e == 1))
        n = np.sum(t >= ti)
        if n > 0:
            s *= (1 - d / n)
        surv.append(s)
        grid.append(ti)
    return np.array(grid), np.array(surv)


def logrank(t1, e1, t2, e2):
    """Two-group log-rank test (Mantel-Haenszel). Returns (stat, p)."""
    t1, e1, t2, e2 = map(np.asarray, (t1, e1, t2, e2))
    times = np.unique(np.concatenate([t1[e1 == 1], t2[e2 == 1]]))
    O, E, V = 0.0, 0.0, 0.0
    for ti in times:
        n1 = np.sum(t1 >= ti); n2 = np.sum(t2 >= ti)
        d1 = np.sum((t1 == ti) & (e1 == 1)); d2 = np.sum((t2 == ti) & (e2 == 1))
        n, d = n1 + n2, d1 + d2
        if n > 1 and d > 0:
            E += d * n1 / n
            O += d1
            V += d * (n1 / n) * (n2 / n) * (n - d) / (n - 1)
    stat = (O - E) ** 2 / V if V > 0 else 0.0
    return stat, float(stats.chi2.sf(stat, 1)), O, E, V


def main():
    mat = build_matrix()
    fmap = map_files_to_samples(mat)
    clin = load_clinical()

    # ---- panel scoring
    scores = run_panel(mat)          # rows=samples(file uuids)
    cov = pd.read_csv(COVERAGE, index_col=0)
    print(cov.to_string())

    z = zscore_df(mat)
    immune = naive_score(z, IMMUNE_GENES, "immune-naive")
    stromal = naive_score(z, STROMAL_GENES, "stromal-naive")

    df = pd.DataFrame({"case_id": fmap["case_id"], "sample_type": fmap["sample_type"]})
    for c in scores.columns:
        df[c] = scores[c].values
    df["immune_naive"] = immune.values
    df["stromal_naive"] = stromal.values
    df = df.join(clin, on="case_id", how="left")
    df["os_years"] = df["os_days"] / 365.25
    df.to_csv(os.path.join(OUT, "scores_clinical.csv"))

    tum = df[df["sample_type"] == "Primary Tumor"].copy()
    tum = tum.dropna(subset=["os_days"])
    print(f"primary tumors with OS: {len(tum)}")

    stats_rows = []

    # ---- (1) survival: median split, iCAF and myCAF (+ apCAF, SASP for completeness)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, sig in zip(axes, ["iCAF", "myCAF"]):
        vals = tum[sig].dropna()
        med = vals.median()
        hi = tum[tum[sig] > med]; lo = tum[tum[sig] <= med]
        st, p, O, E, V = logrank(hi["os_days"].values, hi["os_event"].values,
                                 lo["os_days"].values, lo["os_event"].values)
        g1, s1 = km_curve(hi["os_days"].values, hi["os_event"].values)
        g0, s0 = km_curve(lo["os_days"].values, lo["os_event"].values)
        ax.step(g1, s1, where="post", label=f"high (n={len(hi)})")
        ax.step(g0, s0, where="post", label=f"low (n={len(lo)})")
        ax.set_title(f"{sig} (median split)\nlog-rank p={p:.3g}")
        ax.set_xlabel("days"); ax.set_ylabel("survival")
        ax.legend(fontsize=8); ax.set_ylim(0, 1.02)
        # effect size: median survival per group (numeric summary appended below
        # for all four signatures to avoid duplicates)
    # all four signatures, numeric summary
    for sig in ["iCAF", "myCAF", "apCAF", "SASP"]:
        vals = tum[sig].dropna()
        med = vals.median()
        hi = tum[tum[sig] > med]; lo = tum[tum[sig] <= med]
        st, p, O, E, V = logrank(hi["os_days"].values, hi["os_event"].values,
                                 lo["os_days"].values, lo["os_event"].values)
        stats_rows.append({"test": f"logrank_OS_{sig}_median_split",
                           "n_high": len(hi), "n_low": len(lo),
                           "statistic": st, "p_value": p,
                           "median_OS_days_high": hi["os_days"][hi["os_event"] == 1].median(),
                           "median_OS_days_low": lo["os_days"][lo["os_event"] == 1].median()})
    fig.suptitle("TCGA-PAAD overall survival by CAF score (median split)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_survival_km.png"), dpi=150)

    # ---- (2) stage
    def stage_group(s):
        if pd.isna(s): return np.nan
        s = str(s)
        if "IV" in s: return "IV"
        if "III" in s: return "III"
        if "II" in s: return "II"
        if "I" in s or "0" in s: return "I"
        return np.nan
    tum["stage_grp"] = tum["stage"].map(stage_group)
    st2 = tum.dropna(subset=["stage_grp"])
    print("stage counts:\n", st2["stage_grp"].value_counts())
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=False)
    for ax, sig in zip(axes, ["iCAF", "myCAF"]):
        groups = [g[sig].dropna().values for _, g in st2.groupby("stage_grp")]
        labels = sorted(st2["stage_grp"].unique())
        kw = stats.kruskal(*groups)
        ax.boxplot(groups, labels=[f"{l}\n(n={len(g)})" for l, g in zip(labels, groups)])
        ax.set_title(f"{sig} by pathologic stage\nKruskal-Wallis p={kw.pvalue:.3g}")
        ax.set_ylabel("signature score (mean z)")
        stats_rows.append({"test": f"kruskal_stage_{sig}", "n": len(st2),
                           "statistic": kw.statistic, "p_value": float(kw.pvalue)})
    fig.suptitle("CAF scores vs pathologic stage (TCGA-PAAD, primary tumors)")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_stage_boxplot.png"), dpi=150)

    # ---- (3) immune / stromal correlation
    corrs = {}
    for sig in ["myCAF", "iCAF", "apCAF", "SASP"]:
        for ref, refname in [("immune_naive", "immune"), ("stromal_naive", "stromal")]:
            x = tum[[sig, ref]].dropna()
            r, p = stats.spearmanr(x[sig], x[ref])
            corrs[(sig, refname)] = (r, p, len(x))
            stats_rows.append({"test": f"spearman_{sig}_vs_{refname}",
                               "n": len(x), "statistic": r, "p_value": p})
            print(f"{sig} vs {refname}: spearman r={r:.3f} p={p:.3g} n={len(x)}")
    fig, ax = plt.subplots(figsize=(6, 5))
    sigs = ["myCAF", "iCAF", "apCAF", "SASP"]
    mat_c = np.array([[corrs[(s, "immune")][0], corrs[(s, "stromal")][0]] for s in sigs])
    im = ax.imshow(mat_c, cmap="coolwarm", vmin=-1, vmax=1)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["immune-naive", "stromal-naive"])
    ax.set_yticks(range(4)); ax.set_yticklabels(sigs)
    for i in range(4):
        for j in range(2):
            r, p, n = corrs[(sigs[i], ["immune", "stromal"][j])]
            ax.text(j, i, f"{r:.2f}\np={p:.1g}", ha="center", va="center", fontsize=9,
                    color="white" if abs(r) > 0.5 else "black")
    fig.colorbar(im, ax=ax, label="Spearman r")
    ax.set_title("CAF signatures vs naive immune/stromal scores")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig_immune_corr.png"), dpi=150)

    # ---- effect-size table
    pd.DataFrame(stats_rows).to_csv(os.path.join(OUT, "stats_summary.csv"), index=False)
    print("wrote output CSVs + PNGs")

    # ---- console summary for the report
    print("\n=== coverage ==="); print(cov.to_string())
    print("\n=== score correlations (tumors) ===")
    print(tum[["myCAF", "iCAF", "apCAF", "SASP"]].corr(method="spearman").round(2).to_string())


if __name__ == "__main__":
    main()
