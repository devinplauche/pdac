#!/usr/bin/env python3
"""RTK expression in PDAC CAFs vs malignant cells, from TISCH2.

Question: do cancer-associated fibroblasts express the RTKs (EGFR, MET,
ERBB2/HER2, ERBB3) that daraxonrasib induces on tumor cells? If yes, chronic
RAS(ON) inhibition may rewire CAFs the same way it rewires tumor cells --
a testable hypothesis for the resistant-CAF project.

Data: TISCH2 (http://tisch.comp-genomics.org), PAAD datasets, precomputed
per-cell-type mean expression queried via the site's own search-gene API.
No S3 access needed (the CELLxGENE Census S3 endpoint is unreachable from
this sandbox).

PAAD datasets (from TISCH2's embedded dataset map):
  PAAD_CRA001160, PAAD_GSE111672, PAAD_GSE141017, PAAD_GSE148673,
  PAAD_GSE154763, PAAD_GSE154778, PAAD_GSE158356, PAAD_GSE162708,
  PAAD_GSE165399
"""
import json
import os
import re
import time

import requests
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = "https://tisch.compbio.cn"
GENES = ["EGFR", "MET", "ERBB2", "ERBB3"]
DATASETS = [
    "PAAD_CRA001160", "PAAD_GSE111672", "PAAD_GSE141017", "PAAD_GSE148673",
    "PAAD_GSE154763", "PAAD_GSE154778", "PAAD_GSE158356", "PAAD_GSE162708",
    "PAAD_GSE165399",
]
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUT, exist_ok=True)


def session_with_csrf():
    s = requests.Session()
    # TISCH2's Django view gates on request.is_ajax(); jQuery sets this header.
    s.headers.update({
        "User-Agent": "Mozilla/5.0",
        "X-Requested-With": "XMLHttpRequest",
        "Referer": f"{BASE}/search-gene/",
    })
    r = s.get(f"{BASE}/search-gene/", timeout=60)
    r.raise_for_status()
    m = re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', r.text)
    return s, m.group(1)


def query_gene(s, token, gene):
    """POST gene search, return parsed heatmap JSON (dataset -> celltype -> value)."""
    payload = [
        ("genesearch", gene),
        ("annotation", "Celltype_curated"),   # major-lineage: Fibroblasts, Malignant, ...
        ("cancer", "PAAD"),
        ("Cell_lineage", "All lineage"),
        ("plottype", "heatmap"),
        ("csrfmiddlewaretoken", token),
    ]
    # jQuery serializes a JS array as dataset[]=...; Django reads getlist('dataset[]')
    for d in DATASETS:
        payload.append(("dataset[]", d))
    resp = s.post(f"{BASE}/search-gene/", data=payload, timeout=300)
    resp.raise_for_status()
    payload = resp.json()
    heatmap_url = payload.get("heatmap_file")
    if not heatmap_url:
        raise RuntimeError(f"no heatmap_file in response: {str(payload)[:300]}")
    if heatmap_url.startswith("/"):
        heatmap_url = BASE + heatmap_url
    hm = s.get(heatmap_url, timeout=120)
    hm.raise_for_status()
    return hm.json(), payload.get("selected_gene", gene)


def main():
    s, token = session_with_csrf()
    records = []
    meta = {"source": "TISCH2 search-gene API", "base": BASE,
            "annotation": "Celltype_curated (major-lineage)",
            "datasets_queried": DATASETS, "genes": GENES, "per_gene": {}}
    for gene in GENES:
        print(f"querying {gene} ...", flush=True)
        for attempt in range(3):
            try:
                data, selected = query_gene(s, token, gene)
                break
            except Exception as e:
                print(f"  attempt {attempt+1} failed: {type(e).__name__}: {str(e)[:150]}")
                time.sleep(10)
                s, token = session_with_csrf()
        else:
            raise RuntimeError(f"failed to query {gene}")
        meta["per_gene"][gene] = {"selected_gene": selected,
                                  "heatmap_keys": list(data.keys())[:10]}
        # Inspect structure on first gene
        if gene == GENES[0]:
            with open(os.path.join(OUT, "_heatmap_raw_EGFR.json"), "w") as f:
                json.dump(data, f)
            print("  heatmap JSON top-level keys:", list(data.keys()))
        # Parse: expected keys like 'x' (celltypes), 'y' (datasets), 'z' (values)
        # fall back to generic handling below
        parsed = parse_heatmap(data, gene)
        records.extend(parsed)
        print(f"  parsed {len(parsed)} dataset x celltype values")
        time.sleep(2)

    df = pd.DataFrame(records)
    df.to_csv(os.path.join(OUT, "rtk_by_celltype_all.csv"), index=False)

    # ---- summary: fibroblasts vs malignant, per gene ----
    focus = df[df["celltype"].isin(["Fibroblasts", "Myofibroblasts", "Malignant"])]
    piv = (focus.groupby(["gene", "celltype"])["value"].mean()
           .unstack("celltype"))
    piv.to_csv(os.path.join(OUT, "rtk_fibroblast_vs_malignant.csv"))

    # ---- plot ----
    genes = GENES
    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for ax, gene in zip(axes.flat, genes):
        sub = piv.loc[gene]
        cts = [c for c in ["Fibroblasts", "Malignant"] if c in sub.index]
        vals = [sub[c] for c in cts]
        bars = ax.bar(cts, vals, color=["#d62728", "#1f77b4"], edgecolor="black")
        ax.set_title(gene, fontweight="bold")
        ax.set_ylabel("mean expression (TISCH2, log-normalized)")
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() * 1.03,
                    f"{v:.2f}", ha="center", fontsize=10)
    fig.suptitle("RTK expression: fibroblasts vs malignant cells in human PDAC\n"
                 "TISCH2, 9 PAAD scRNA-seq datasets, major-lineage annotation "
                 "(mean across datasets)",
                 fontsize=12, fontweight="bold")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "rtk_fibroblast_vs_malignant.png"), dpi=150)

    with open(os.path.join(OUT, "meta.json"), "w") as f:
        json.dump(meta, f, indent=2)
    print(piv.to_string())
    print("\nDone. Outputs in", OUT)


def parse_heatmap(data, gene):
    """Extract (dataset, celltype, value) triples from TISCH2 heatmap JSON.

    Structure: {'celltype': [...], 'dataset': [...],
                'expression': [[celltype_idx, dataset_idx, value], ...]}
    """
    recs = []
    cts, dss = data.get("celltype", []), data.get("dataset", [])
    for row in data.get("expression", []):
        try:
            ci, di, val = int(row[0]), int(row[1]), float(row[2])
        except (TypeError, ValueError, IndexError):
            continue
        if 0 <= ci < len(cts) and 0 <= di < len(dss):
            recs.append({"gene": gene, "dataset": dss[di],
                         "celltype": cts[ci], "value": val})
    if not recs:
        raise RuntimeError(f"could not parse heatmap JSON for {gene}; "
                           f"keys={list(data.keys()) if isinstance(data, dict) else type(data)}")
    return recs


if __name__ == "__main__":
    main()
