#!/usr/bin/env python3
"""Cross-tissue validation of KRAS-inhibitor CAF reprogramming in GSE336609.

GSE336609 (Shen et al., J Hematol Oncol 2026): human appendiceal
adenocarcinoma PDX treated with KRAS inhibitors vs vehicle; bulk RNA-seq
with human (tumor) and mouse (stroma/host) reads extracted separately.
CRITICAL: raw/scRNA-seq data were NOT deposited (patient privacy) — GEO
holds only per-sample count matrices (treated vs vehicle contrast files).
So fibroblast-level per-cell scoring (as in the GSE337208 mining) is
impossible here; this analysis scores the bulk stromal (mouse-read) signal
with the repo's validated signature panel AS-IS (signatures_v1.yaml,
--organism mouse), and runs treated-vs-vehicle contrasts per model.

Design:
  AAP01 (KRAS G12D): MRTX1133 (G12D-selective) n=3 vs vehicle n=3
  AAP16 (KRAS G12V): RMC-6236 (pan-KRAS)   n=2 vs vehicle n=3
Mouse reads = host stroma (fibroblasts + endothelium + immune). Human reads
= tumor compartment (on-target PD sanity check: MAPK/RASON signatures).

Statistics (honest labels):
  - Per-sample rank-method scores (AUCell-like, repo scorer as-is).
  - Mann-Whitney U treated vs vehicle per model per signature (sample-level;
    underpowered by construction: min achievable two-sided p is 0.10 for
    3v3 and 0.20 for 2v3 — reported as exact p with this caveat).
  - Member-gene coherence: median log2FC per signature's genes + Wilcoxon
    signed-rank on member-gene log2FCs (genes as observations; labeled
    gene-level/descriptive, genes are not independent).
  - Compartment-composition check (fibroblast vs immune vs endothelial
    marker fractions) to ground whether the mouse signal is fibroblast-
    dominated or could reflect immune infiltration.

Re-runnable: python3 run.py  (downloads GEO files if missing; resolves the
~150 signature/marker symbols to matrix Ensembl ids via one mygene.info
batch query per organism and caches data/symbol_to_ensembl_{mouse,human}.tsv;
non-target rows keep Ensembl ids, which leaves the rank distribution used by
the repo scorer unchanged)
"""
import gzip
import json
import os
import subprocess
import sys
import time
import urllib.request

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, wilcoxon

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "output")
os.makedirs(DATA, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

SCORER = os.path.abspath(os.path.join(HERE, "..", "signature_panel",
                                       "score_signatures.py"))
SIGYAML = os.path.abspath(os.path.join(HERE, "..", "signature_panel",
                                        "signatures_v1.yaml"))

GEO_BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE336nnn/GSE336609/suppl"
GEO_FILES = [
    "GSE336609_AAP01_KRAS-G12D_AA_Human_Reads_MRTX1133_vs_Ctrl.csv.gz",
    "GSE336609_AAP01_KRAS-G12D_AA_Mouse_Reads_MRTX1133_vs_Ctrl.csv.gz",
    "GSE336609_AAP16_KRAS-G12V_AA_Human_Reads_RMC-6236_vs_Ctrl.csv.gz",
    "GSE336609_AAP16_KRAS-G12V_AA_Mouse_Reads_RMC-6236_vs_Ctrl.csv.gz",
]

# model -> (mouse file key, human file key, drug, trt prefix, veh prefix)
MODELS = {
    "AAP01": ("Mouse", "MRTX1133", ["AAP01_MRTX1_Mouse", "AAP01_MRTX2_Mouse",
                                     "AAP01_MRTX3_Mouse"],
              ["AAP01_VC1_Mouse", "AAP01_VC2_Mouse", "AAP01_VC3_Mouse"],
              "KRAS G12D", "G12D-selective"),
    "AAP16": ("Mouse", "RMC-6236", ["AAP16_RMC1_Mouse", "AAP16_RMC2_Mouse"],
              ["AAP16_VC1_Mouse", "AAP16_VC2_Mouse", "AAP16_VC3_Mouse"],
              "KRAS G12V", "pan-KRAS"),
}

COMPARTMENTS = {  # mouse symbols; host-stroma composition check
    "fibroblast": ["Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Pdpn",
                   "Sparc", "Fn1", "Vim", "Postn", "Tnc"],
    "myofibroblast_smc": ["Acta2", "Tagln", "Myh11"],
    "pericyte": ["Pdgfrb", "Rgs5", "Cspg4"],
    "immune": ["Ptprc", "Cd68", "Cd3e", "Cd4", "Cd8a", "Cd19",
               "Itgam", "Lyz2", "Cd14", "Nkg7", "Adgre1"],
    "endothelial": ["Pecam1", "Cdh5", "Vwf", "Kdr", "Tek"],
    "epithelial_host": ["Epcam", "Krt19", "Krt8", "Krt18", "Cdh1"],
}


def download_geo():
    for f in GEO_FILES:
        p = os.path.join(DATA, f)
        if os.path.exists(p) and os.path.getsize(p) > 0:
            print(f"  have {f}", flush=True)
            continue
        print(f"  downloading {f} ...", flush=True)
        urllib.request.urlretrieve(f"{GEO_BASE}/{f}", p)


def strip_version(eid):
    return eid.split(".")[0]


def hgnc_human_map(symbols):
    """Map human symbols -> Ensembl gene ids via the HGNC REST API.

    One GET per symbol (fast, reliable); mygene.info batch queries were
    intermittently refused for the human batches. Returns
    {symbol: [ensembl ids]}.
    """
    import urllib.request as urlreq
    out = {}
    symbols = sorted(set(symbols))
    for i, sym in enumerate(symbols):
        url = f"https://rest.genenames.org/fetch/symbol/{sym}"
        req = urlreq.Request(url, headers={"Accept": "application/json"})
        for attempt in range(6):
            try:
                with urlreq.urlopen(req, timeout=30) as fh:
                    d = json.loads(fh.read())
                break
            except Exception as e:
                if attempt == 5:
                    print(f"  HGNC FAILED for {sym}: {e}", flush=True)
                    d = None
                else:
                    time.sleep(2 + attempt * 2)
        if not d:
            continue
        docs = d.get("response", {}).get("docs", [])
        if docs:
            eg = docs[0].get("ensembl_gene_id")
            if eg:
                out[sym] = [strip_version(eg)]
        if (i + 1) % 20 == 0:
            print(f"  HGNC human: {i + 1}/{len(symbols)}", flush=True)
    return out


def mygene_symbol_to_ensembl(symbols, species, batch=12):
    """Map gene symbols -> Ensembl gene ids via mygene.info batch queries.

    POST /query treats each whitespace-separated token of q as an independent
    query against scopes=symbol and returns a per-token result list. Small
    batches: the API intermittently resets larger responses from this host.
    Returns {symbol: [ensembl ids]}.
    """
    import urllib.parse
    import urllib.request as urlreq
    symbols = sorted(set(symbols))
    out = {}

    def fetch(chunk):
        data = urllib.parse.urlencode({
            "q": " ".join(chunk),
            "scopes": "symbol",
            "species": species,
            "fields": "symbol,ensembl.gene",
            "size": len(chunk),
        }).encode()
        req = urlreq.Request("https://mygene.info/v3/query", data=data)
        res = None
        for attempt in range(8):
            try:
                with urlreq.urlopen(req, timeout=120) as fh:
                    res = json.loads(fh.read())
                break
            except Exception as e:
                print(f"  mygene retry {attempt + 1} ({species}): {e}",
                      flush=True)
                time.sleep(3 + attempt * 2)
        if res is None:
            raise RuntimeError(
                f"mygene.info batch query failed ({species}, "
                f"chunk {chunk[:3]}...)")
        return res

    for i in range(0, len(symbols), batch):
        chunk = symbols[i:i + batch]
        for r in fetch(chunk):
            if r.get("notfound"):
                continue
            sym = r.get("symbol")
            e = r.get("ensembl")
            ids = []
            if isinstance(e, dict):
                e = [e]
            for x in (e or []):
                g = x.get("gene")
                if g:
                    ids.append(strip_version(g))
            if sym and ids:
                out.setdefault(sym, [])
                for x in ids:
                    if x not in out[sym]:
                        out[sym].append(x)
        if (i // batch) % 5 == 0:
            print(f"  mygene {species}: {i + len(chunk)}/{len(symbols)}",
                  flush=True)
    return out


def target_symbols(comp):
    """All symbols we need resolved: signature genes + compartment markers."""
    syms_mouse, syms_human = set(), set()
    for s in comp.values():
        syms_mouse.update(s["mouse"])
        syms_human.update(s["human"])
    for genes in COMPARTMENTS.values():
        syms_mouse.update(genes)
    return syms_mouse, syms_human


def resolve_maps_if_missing(comp):
    """Resolve target symbols -> matrix Ensembl ids (one mygene batch each).

    We deliberately do NOT map all ~53k matrix rows: the repo scorer's rank
    method computes percentile ranks over the full gene distribution, so
    non-signature rows keep their (unique) Ensembl ids and only target genes
    are relabeled to symbols. Saved to data/symbol_to_ensembl_{org}.tsv.
    """
    syms_mouse, syms_human = target_symbols(comp)
    got = {}
    for org, species, syms in (("mouse", "mouse", syms_mouse),
                               ("human", "human", syms_human)):
        mp = os.path.join(DATA, f"symbol_to_ensembl_{org}.tsv")
        if os.path.exists(mp):
            d = pd.read_csv(mp, sep="\t")
            m = {}
            for _, r in d.iterrows():
                m.setdefault(r["symbol"], []).append(r["ensembl"])
            got[org] = m
            print(f"  have {mp} ({len(m)} symbols)", flush=True)
            continue
        print(f"  resolving {len(syms)} {org} symbols ...", flush=True)
        if org == "human":
            m = hgnc_human_map(syms)
        else:
            m = mygene_symbol_to_ensembl(syms, species)
        pd.DataFrame(
            [(s, e) for s, el in sorted(m.items()) for e in el],
            columns=["symbol", "ensembl"]).to_csv(mp, sep="\t", index=False)
        got[org] = m
        print(f"  {org}: {len(m)}/{len(syms)} symbols resolved", flush=True)
    return got


def load_matrix(compartment, sym2ens):
    """Load count matrix (genes x samples).

    Rows: target genes relabeled to official symbols via sym2ens
    ({symbol: [ensembl ids]}); all other rows keep version-stripped Ensembl
    ids (unique, so the rank distribution is unchanged). Duplicate symbols
    (multiple Ensembl ids -> one symbol) are summed.
    """
    org = "mouse" if compartment == "Mouse" else "human"
    ens2sym = {}
    for sym, el in sym2ens[org].items():
        for e in el:
            ens2sym[e] = sym
    frames = []
    for f in GEO_FILES:
        if f"_{compartment}_" not in f:
            continue
        model = "AAP01" if "AAP01" in f else "AAP16"
        with gzip.open(os.path.join(DATA, f), "rt") as fh:
            df = pd.read_csv(fh)
        ids = df.iloc[:, 0].map(strip_version)
        df = df.set_index(ids)
        df.index.name = "gene"
        df = df[[c for c in df.columns if c != df.columns[0]]]
        df = df.rename(index=ens2sym)  # only target genes relabeled
        before = len(df)
        df = df.groupby(level=0).sum().astype(int)  # collapse dups
        n_sym = sum(1 for i in df.index
                    if not str(i).startswith("ENS"))
        print(f"  {f}: {before} rows -> {len(df)} unique "
              f"({n_sym} symbol-labeled)", flush=True)
        df.columns = [f"{model}|{c}" for c in df.columns]
        frames.append(df)
    mat = pd.concat(frames, axis=1).fillna(0).astype(int)
    mat.index.name = "gene"
    return mat


def compartment_fractions(mat):
    rows = []
    lib = mat.sum(axis=0)
    for comp, genes in COMPARTMENTS.items():
        present = [g for g in genes if g in mat.index]
        frac = (mat.loc[present].sum(axis=0) / lib * 100) if present \
            else pd.Series(0.0, index=mat.columns)
        for s in mat.columns:
            rows.append({"sample": s, "compartment": comp,
                         "genes_found": f"{len(present)}/{len(genes)}",
                         "pct_of_mouse_reads": round(float(frac[s]), 2)})
    return pd.DataFrame(rows)


def score_with_repo(mat, org, tag):
    inp = os.path.join(DATA, f"{tag}_symbol_counts.tsv")
    mat.reset_index().rename(columns={"index": "gene"}).to_csv(
        inp, sep="\t", index=False)
    out = {}
    for method in ("rank", "zscore"):
        sp = os.path.join(OUT, f"{tag}_scores_{method}.csv")
        cp = os.path.join(OUT, f"{tag}_gene_coverage.csv")
        r = subprocess.run(
            [sys.executable, SCORER, "--input", inp, "--signatures", SIGYAML,
             "--organism", org, "--method", method, "--out", sp,
             "--coverage-out", cp] + (["--no-log"] if False else []),
            capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr, file=sys.stderr)
            sys.exit(1)
        out[method] = pd.read_csv(sp, index_col="sample")
    return out["rank"], out["zscore"]


def log2cpm(mat):
    lib = mat.sum(axis=0).replace(0, np.nan)
    return np.log2(mat.div(lib, axis=1) * 1e6 + 1)


def member_log2fc(lc, sig_genes, trt_cols, veh_cols):
    rows = []
    for g in sig_genes:
        if g not in lc.index:
            rows.append({"gene": g, "log2FC_trt_vs_veh": np.nan,
                         "detected": False})
            continue
        fc = float(lc.loc[g, trt_cols].median() -
                   lc.loc[g, veh_cols].median())
        rows.append({"gene": g, "log2FC_trt_vs_veh": round(fc, 3),
                     "detected": True})
    return pd.DataFrame(rows)


def main():
    import yaml
    print("== download ==", flush=True)
    download_geo()
    print("== gene map ==", flush=True)
    with open(SIGYAML) as f:
        comp = yaml.safe_load(f)["signatures"]
    sym2ens = resolve_maps_if_missing(comp)
    sig_genes_mouse = {s: comp[s]["mouse"] for s in comp}
    sig_genes_human = {s: comp[s]["human"] for s in comp}

    # ---------------- mouse (stroma) ----------------
    print("== mouse matrix ==", flush=True)
    mmat = load_matrix("Mouse", sym2ens)
    comp_frac = compartment_fractions(mmat)
    comp_frac.to_csv(os.path.join(OUT, "mouse_compartment_fractions.csv"),
                     index=False)
    print(comp_frac.pivot_table(index="sample", columns="compartment",
                                values="pct_of_mouse_reads").to_string())
    mlc = log2cpm(mmat)
    print("== scoring mouse (rank + zscore) ==", flush=True)
    mrank, mz = score_with_repo(mmat, "mouse", "mouse")

    # ---------------- human (tumor, PD check) ----------------
    print("== human matrix ==", flush=True)
    hmat = load_matrix("Human", sym2ens)
    hlc = log2cpm(hmat)
    print("== scoring human (rank + zscore) ==", flush=True)
    hrank, hz = score_with_repo(hmat, "human", "human")

    # ---------------- contrasts ----------------
    stat_rows, gene_fc_rows = [], []
    for model, (comp_key, drug, trt, veh, kras, sel) in MODELS.items():
        tcols = [f"{model}|{c}" for c in trt]
        vcols = [f"{model}|{c}" for c in veh]
        for sname in comp:
            a = mrank.loc[vcols, sname].dropna()
            b = mrank.loc[tcols, sname].dropna()
            u, p = mannwhitneyu(a, b, alternative="two-sided")
            gfc = member_log2fc(mlc, sig_genes_mouse[sname], tcols, vcols)
            det = gfc.dropna(subset=["log2FC_trt_vs_veh"])
            if len(det) >= 3:
                _, pw = wilcoxon(det["log2FC_trt_vs_veh"])
            else:
                pw = np.nan
            gfc["model"], gfc["signature"], gfc["compartment"] = \
                model, sname, "stroma(mouse)"
            gene_fc_rows.append(gfc)
            stat_rows.append({
                "model": model, "kras": kras, "drug": drug,
                "signature": sname,
                "compartment": "stroma(mouse)",
                "n_veh": len(a), "n_trt": len(b),
                "median_score_veh": round(float(a.median()), 4),
                "median_score_trt": round(float(b.median()), 4),
                "delta_trt_minus_veh": round(float(b.median() - a.median()),
                                             4),
                "MWU_p_sample_level": float(p),
                "member_genes_found":
                    f"{len(det)}/{len(sig_genes_mouse[sname])}",
                "median_member_log2FC":
                    round(float(det["log2FC_trt_vs_veh"].median()), 3)
                    if len(det) else np.nan,
                "signedrank_p_gene_level": float(pw)
                    if not np.isnan(pw) else np.nan,
            })
        # human PD signatures
        htcols = [c.replace("_Mouse", "_Human") for c in tcols]
        hvcols = [c.replace("_Mouse", "_Human") for c in vcols]
        for sname in ("MAPK_targets", "RASON_inhibition_PD"):
            a = hrank.loc[hvcols, sname].dropna()
            b = hrank.loc[htcols, sname].dropna()
            u, p = mannwhitneyu(a, b, alternative="two-sided")
            gfc = member_log2fc(hlc, sig_genes_human[sname], htcols, hvcols)
            det = gfc.dropna(subset=["log2FC_trt_vs_veh"])
            _, pw = wilcoxon(det["log2FC_trt_vs_veh"]) if len(det) >= 3 \
                else (np.nan, np.nan)
            gfc["model"], gfc["signature"], gfc["compartment"] = \
                model, sname, "tumor(human)"
            gene_fc_rows.append(gfc)
            stat_rows.append({
                "model": model, "kras": kras, "drug": drug,
                "signature": sname,
                "compartment": "tumor(human)",
                "n_veh": len(a), "n_trt": len(b),
                "median_score_veh": round(float(a.median()), 4),
                "median_score_trt": round(float(b.median()), 4),
                "delta_trt_minus_veh": round(float(b.median() - a.median()),
                                             4),
                "MWU_p_sample_level": float(p),
                "member_genes_found":
                    f"{len(det)}/{len(sig_genes_human[sname])}",
                "median_member_log2FC":
                    round(float(det["log2FC_trt_vs_veh"].median()), 3)
                    if len(det) else np.nan,
                "signedrank_p_gene_level": float(pw)
                    if not (isinstance(pw, float) and np.isnan(pw))
                    else np.nan,
            })
    stats = pd.DataFrame(stat_rows)
    stats.to_csv(os.path.join(OUT, "contrast_stats.csv"), index=False)
    pd.concat(gene_fc_rows, ignore_index=True).to_csv(
        os.path.join(OUT, "member_gene_log2fc.csv"), index=False)

    # attach meta to per-sample scores
    def with_meta(df):
        d = df.copy()
        d["model"] = d.index.map(lambda s: s.split("|")[0])
        d["sample"] = d.index.map(lambda s: s.split("|")[1])
        d["treatment"] = d["sample"].map(
            lambda s: "treated" if ("MRTX" in s or "_RMC" in s) else "vehicle")
        return d
    with_meta(mrank).to_csv(os.path.join(OUT, "mouse_scores_per_sample.csv"))
    with_meta(hrank).to_csv(os.path.join(OUT, "human_scores_per_sample.csv"))
    with_meta(mz).to_csv(os.path.join(OUT,
                                      "mouse_scores_per_sample_zscore.csv"))

    qc_rows = []
    for s in mmat.columns:
        lib = int(mmat[s].sum())
        ng = int((mmat[s] > 0).sum())
        qc_rows.append({"sample": s, "compartment": "stroma(mouse)",
                        "library_size": lib, "genes_detected": ng})
    for s in hmat.columns:
        lib = int(hmat[s].sum())
        ng = int((hmat[s] > 0).sum())
        qc_rows.append({"sample": s, "compartment": "tumor(human)",
                        "library_size": lib, "genes_detected": ng})
    pd.DataFrame(qc_rows).to_csv(os.path.join(OUT, "qc_table.csv"),
                                 index=False)

    print()
    print(stats[["model", "drug", "signature", "compartment",
                 "median_score_veh", "median_score_trt",
                 "delta_trt_minus_veh", "MWU_p_sample_level",
                 "median_member_log2FC",
                 "signedrank_p_gene_level"]].to_string(index=False))
    print()
    print("KEY: min achievable two-sided MWU p = 0.10 for 3v3 (AAP01), "
          "0.20 for 2v3 (AAP16).")

    # ---------------- figures ----------------
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    caf_sigs = ["myCAF", "iCAF", "apCAF", "SASP", "matrisome_core",
                "MAPK_targets", "RASON_inhibition_PD"]
    labels = {"myCAF": "myCAF", "iCAF": "iCAF", "apCAF": "apCAF",
              "SASP": "SASP", "matrisome_core": "matrisome",
              "MAPK_targets": "MAPK out", "RASON_inhibition_PD": "RAS(ON) PD"}

    # Fig 1: per-sample rank scores, stroma, both models
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.2), sharex=True)
    for ax, model in zip(axes, ["AAP01", "AAP16"]):
        drug = MODELS[model][1]
        tcols = [f"{model}|{c}" for c in MODELS[model][2]]
        vcols = [f"{model}|{c}" for c in MODELS[model][3]]
        y = np.arange(len(caf_sigs))
        for cols, color, lab, ms in ((vcols, "#7f7f7f", "vehicle", "o"),
                                     (tcols, "#d62728", drug, "D")):
            xs = mrank.loc[cols, caf_sigs].values
            for j in range(len(caf_sigs)):
                ax.scatter(xs[:, j], np.full(len(cols), y[j]) +
                           np.linspace(-0.12, 0.12, len(cols)),
                           color=color, s=64, marker=ms, zorder=3,
                           label=lab if j == 0 else None)
            ax.scatter(np.median(xs, axis=0), y, color="black", s=200,
                       marker="|", linewidths=3, zorder=4)
        ax.set_yticks(y)
        ax.set_yticklabels([labels[s] for s in caf_sigs])
        ax.set_title(f"{model} stroma (mouse reads)\n{drug} vs vehicle",
                     fontsize=11)
        ax.set_xlabel("per-sample rank score (AUCell-like)")
        ax.grid(axis="x", alpha=0.3)
    axes[0].legend(frameon=False, loc="lower right")
    fig.suptitle("KRAS inhibition vs vehicle — stromal CAF signatures "
                 "(GSE336609, per-sample rank scores; tick = median)",
                 fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig1_stromal_caf_scores.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)

    # Fig 2: member-gene log2FCs, stroma, both models
    gene_fc = pd.read_csv(os.path.join(OUT, "member_gene_log2fc.csv"))
    gf = gene_fc[gene_fc.compartment == "stroma(mouse)"]
    focus = ["myCAF", "iCAF", "apCAF"]
    fig, axes = plt.subplots(1, 3, figsize=(13, 5.5), sharex=True)
    for ax, sname in zip(axes, focus):
        sub = gf[gf.signature == sname].copy()
        sub = sub.dropna(subset=["log2FC_trt_vs_veh"])
        genes = list(dict.fromkeys(sub["gene"]))
        y = np.arange(len(genes))
        for model, color, mk in (("AAP01", "#1f77b4", "o"),
                                 ("AAP16", "#9467bd", "s")):
            m = sub[sub.model == model].set_index("gene")
            xs = [m.loc[g, "log2FC_trt_vs_veh"] if g in m.index else np.nan
                  for g in genes]
            ax.scatter(xs, y + (0.15 if model == "AAP01" else -0.15),
                       color=color, s=70, marker=mk, label=model)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yticks(y)
        ax.set_yticklabels(genes, fontsize=9)
        ax.set_title(sname, fontsize=11)
        ax.set_xlabel("log2FC treated vs vehicle\n(median, per model)")
        ax.grid(axis="x", alpha=0.3)
    axes[0].legend(frameon=False, title="model")
    fig.suptitle("Stromal CAF-signature member genes — treated vs vehicle "
                 "log2FC (GSE336609)", fontsize=12, y=1.02)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "fig2_member_gene_log2fc.png"), dpi=150,
                bbox_inches="tight")
    plt.close(fig)
    print("figures written")
    print("ALL DONE")


if __name__ == "__main__":
    main()
