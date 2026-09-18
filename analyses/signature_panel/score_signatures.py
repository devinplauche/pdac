#!/usr/bin/env python3
"""Reusable signature scorer.

Scores gene signatures (from a YAML compendium like signatures_v1.yaml)
on an expression matrix.

Input:  CSV/TSV with one row per gene, first column = gene symbol,
        remaining columns = samples (TPM, counts, or log-normalized values).
Output: CSV with one row per sample and one column per signature.

Methods:
  zscore (default): log2(x+1) transform, per-gene z-score across samples,
                    signature score = mean z of member genes present.
  rank:             per-sample fractional rank of each gene (0..1),
                    signature score = mean rank fraction of member genes.
                    (AUCell-like; robust to normalization, good for scRNA-seq.)

Gene coverage: the scorer reports, per signature, how many member genes were
found in the matrix. Scores from <3 genes are still computed but flagged.

Usage:
    python score_signatures.py --input expr.csv --signatures signatures_v1.yaml \\
        --organism mouse --out scores.csv
"""
import argparse
import sys

import numpy as np
import pandas as pd
import yaml


def load_signatures(path, organism):
    with open(path) as f:
        comp = yaml.safe_load(f)
    sigs = {}
    for name, s in comp["signatures"].items():
        genes = s.get(organism)
        if genes is None:
            # fallback: title-case heuristic, e.g. DUSP6 -> Dusp6
            genes = [g.capitalize() for g in s["human"]]
        sigs[name] = {"label": s.get("label", name), "genes": genes,
                      "direction": s.get("direction", "increases with program")}
    return sigs, comp.get("version")


def score_zscore(expr, genes):
    """expr: DataFrame genes x samples (log2-transformed already)."""
    present = [g for g in genes if g in expr.index]
    if not present:
        return np.nan, 0
    sub = expr.loc[present]
    z = (sub.sub(sub.mean(axis=1), axis=0)
             .div(sub.std(axis=1, ddof=1).replace(0, np.nan), axis=0))
    return z.mean(axis=0), len(present)


def score_rank(expr, genes):
    """Rank-based (AUCell-like): mean fractional rank of member genes per sample."""
    present = [g for g in genes if g in expr.index]
    if not present:
        return np.nan, 0
    ranks = expr.rank(axis=0, pct=True)
    return ranks.loc[present].mean(axis=0), len(present)


def main():
    ap = argparse.ArgumentParser(description="Score gene signatures on an expression matrix.")
    ap.add_argument("--input", required=True, help="CSV/TSV: rows=genes, cols=samples")
    ap.add_argument("--signatures", required=True, help="YAML signature compendium")
    ap.add_argument("--organism", default="mouse", choices=["mouse", "human"])
    ap.add_argument("--gene-col", default=None, help="gene symbol column name (default: first column)")
    ap.add_argument("--sep", default=None, help="delimiter (default: auto-detect)")
    ap.add_argument("--method", default="zscore", choices=["zscore", "rank"])
    ap.add_argument("--no-log", action="store_true",
                    help="skip log2(x+1) transform (input already logged)")
    ap.add_argument("--out", required=True, help="output scores CSV path")
    ap.add_argument("--coverage-out", default=None, help="optional per-signature coverage CSV path")
    args = ap.parse_args()

    sep = args.sep or ("\t" if args.input.endswith((".tsv", ".txt")) else ",")
    df = pd.read_csv(args.input, sep=sep)
    gene_col = args.gene_col or df.columns[0]
    expr = df.set_index(gene_col)
    expr = expr.apply(pd.to_numeric, errors="coerce")
    expr = expr.dropna(axis=1, how="all").fillna(0)

    if args.method == "zscore" and not args.no_log:
        expr = np.log2(expr + 1)

    sigs, version = load_signatures(args.signatures, args.organism)
    scores, coverage = {}, {}
    for name, s in sigs.items():
        fn = score_zscore if args.method == "zscore" else score_rank
        sc, n = fn(expr, s["genes"])
        scores[name] = sc if not np.isscalar(sc) else pd.Series(
            [sc] * expr.shape[1], index=expr.columns)
        coverage[name] = {"label": s["label"], "n_genes": len(s["genes"]),
                          "n_found": n, "direction": s["direction"]}

    out = pd.DataFrame(scores)
    out.index.name = "sample"
    out.to_csv(args.out)
    print(f"compendium v{version}; method={args.method}; "
          f"{out.shape[0]} samples x {out.shape[1]} signatures -> {args.out}")
    cov = pd.DataFrame(coverage).T
    print(cov[["label", "n_genes", "n_found"]].to_string())
    if args.coverage_out:
        cov.to_csv(args.coverage_out)
    missing = cov[cov["n_found"] < 3]
    if not missing.empty:
        print(f"WARNING: <3 genes found for: {list(missing.index)}", file=sys.stderr)


if __name__ == "__main__":
    main()
