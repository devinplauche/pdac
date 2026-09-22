#!/usr/bin/env python3
"""
GSE335219 MAF comparison: HPAC parental P18 vs daraxonrasib-resistant P28.

Aronchik et al., Nature Medicine 2026 — deposited data are MAFs only (2 samples).
No copy-number segments deposited, so no direct CNV verification is possible.
This script does the mutation-level comparison the deposited data support:
  1. Shared vs private high-confidence somatic-ish variants (PASS filter,
     protein-coding consequence, gnomAD AF < 1% or blank).
  2. KRAS: list every KRAS variant in both samples; check for secondary
     KRAS mutations private to the resistant sample (paper says none).
  3. MAPK / PI3K / RTK / MYC pathway mutations, split shared vs private.

Outputs -> ../output/
  variant_comparison_summary.csv
  private_mutations.csv
  kras_mutations.csv
  pathway_private_mutations.csv
  vaf_shared_vs_private.png

Re-runnable: python3 run.py  (expects ../data/*.maf.gz)
"""
import csv, gzip, os, sys
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)

PARENTAL = os.path.join(DATA, "GSM9808344_HPAC-Parental-P18-1mil-SB-10-14-24_markdup_tnscope.vcf.maf.gz")
RESIST = os.path.join(DATA, "GSM9808345_HPAC-0-03-R9945-0-1uM-9945-P28-SB-1mil-10-16-24_markdup_tnscope.vcf.maf.gz")

CODING = {
    "Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins",
    "In_Frame_Del", "In_Frame_Ins", "Splice_Site", "Translation_Start_Site",
    "Nonstop_Mutation", "De_novo_Start_InFrame", "De_novo_Start_OutOfFrame",
}

MAPK = {"KRAS", "NRAS", "HRAS", "BRAF", "RAF1", "ARAF", "MAP2K1", "MAP2K2",
        "MAPK1", "MAPK3", "NF1", "SOS1", "SHOC2", "PTPN11", "RIT1", "MRAS"}
PI3K = {"PIK3CA", "PIK3CB", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2",
        "PTEN", "AKT1", "AKT2", "AKT3", "MTOR", "TSC1", "TSC2", "RHEB"}
RTK = {"EGFR", "ERBB2", "ERBB3", "ERBB4", "MET", "FGFR1", "FGFR2", "FGFR3",
       "FGFR4", "PDGFRA", "PDGFRB", "KIT", "RET", "ALK", "ROS1", "NTRK1",
       "NTRK2", "NTRK3", "IGF1R", "KDR", "FLT1", "FLT4", "TEK", "AXL", "DDR2"}
MYC = {"MYC", "MYCN", "MYCL"}
PATHWAY = {g: "MAPK" for g in MAPK} | {g: "PI3K" for g in PI3K} | \
          {g: "RTK" for g in RTK} | {g: "MYC" for g in MYC}


def af_ok(af):
    if af in ("", ".", None):
        return True  # absent from gnomAD
    try:
        return float(af) < 0.01
    except ValueError:
        return True


def vaf(t_ref, t_alt, depth):
    try:
        r, a, d = int(t_ref), int(t_alt), int(depth)
    except (ValueError, TypeError):
        return None
    tot = r + a
    return (a / tot) if tot > 0 else None


def load(path, label):
    """Stream MAF, keep high-confidence somatic-ish coding variants."""
    variants = {}   # key -> dict
    n_total = n_pass = n_coding = 0
    with gzip.open(path, "rt") as f:
        header = None
        for line in f:
            if line.startswith("#"):
                continue
            row = line.rstrip("\n").split("\t")
            if header is None:
                header = row
                idx = {c: i for i, c in enumerate(header)}
                def g(r, c):
                    i = idx.get(c)
                    return r[i] if i is not None and i < len(r) else ""
                continue
            n_total += 1
            if g(row, "FILTER") != "PASS":
                continue
            n_pass += 1
            if g(row, "Variant_Classification") not in CODING:
                continue
            if not af_ok(g(row, "gnomAD_AF")):
                continue
            n_coding += 1
            ref = g(row, "Reference_Allele")
            alt = g(row, "Tumor_Seq_Allele2")
            if alt == ref:
                alt = g(row, "Tumor_Seq_Allele1")
            key = (g(row, "Chromosome"), g(row, "Start_Position"),
                   g(row, "End_Position"), ref, alt)
            variants[key] = {
                "label": label,
                "gene": g(row, "Hugo_Symbol"),
                "class": g(row, "Variant_Classification"),
                "ref": ref, "alt": alt,
                "hgvsp": g(row, "HGVSp_Short"),
                "hgvsc": g(row, "HGVSc"),
                "vaf": vaf(g(row, "t_ref_count"), g(row, "t_alt_count"),
                           g(row, "t_depth")),
                "depth": g(row, "t_depth"),
                "pathway": PATHWAY.get(g(row, "Hugo_Symbol"), ""),
            }
    print(f"{label}: total={n_total} pass={n_pass} coding_somatic={len(variants)}")
    return variants


def main():
    for p in (PARENTAL, RESIST):
        assert os.path.exists(p), f"missing {p} — extract GSE335219_RAW.tar into data/ first"

    par = load(PARENTAL, "parental_P18")
    res = load(RESIST, "resistant_P28")

    keys_par, keys_res = set(par), set(res)
    shared = keys_par & keys_res
    par_only = keys_par - keys_res
    res_only = keys_res - keys_par
    print(f"shared={len(shared)} parental_private={len(par_only)} resistant_private={len(res_only)}")

    # ---- summary CSV ----
    with open(os.path.join(OUT, "variant_comparison_summary.csv"), "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "value"])
        w.writerows([
            ["parental_P18_coding_somatic_variants", len(par)],
            ["resistant_P28_coding_somatic_variants", len(res)],
            ["shared_variants", len(shared)],
            ["parental_private_variants", len(par_only)],
            ["resistant_private_variants", len(res_only)],
        ])

    # ---- private mutations CSV ----
    def priv_rows(keys, src, status):
        rows = []
        for k in sorted(keys):
            v = src[k]
            rows.append({
                "sample_status": status,
                "gene": v["gene"],
                "chromosome": k[0], "start": k[1], "end": k[2],
                "ref": v["ref"], "alt": v["alt"],
                "variant_classification": v["class"],
                "hgvsp": v["hgvsp"], "hgvsc": v["hgvsc"],
                "vaf": round(v["vaf"], 4) if v["vaf"] is not None else "",
                "depth": v["depth"],
                "pathway": v["pathway"],
            })
        return rows

    priv = priv_rows(par_only, par, "parental_private") + \
           priv_rows(res_only, res, "resistant_private")
    with open(os.path.join(OUT, "private_mutations.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(priv[0].keys()))
        w.writeheader(); w.writerows(priv)

    # ---- KRAS table ----
    def kras_rows(src, status):
        rows = []
        for k, v in src.items():
            if v["gene"] == "KRAS":
                rows.append({"sample_status": status, "chromosome": k[0],
                             "start": k[1], "end": k[2], "ref": v["ref"],
                             "alt": v["alt"],
                             "variant_classification": v["class"],
                             "hgvsp": v["hgvsp"], "hgvsc": v["hgvsc"],
                             "vaf": round(v["vaf"], 4) if v["vaf"] is not None else "",
                             "depth": v["depth"]})
        return sorted(rows, key=lambda r: r["start"])

    kr = kras_rows(par, "parental_P18") + kras_rows(res, "resistant_P28")
    kras_keys_par = {k for k, v in par.items() if v["gene"] == "KRAS"}
    kras_keys_res = {k for k, v in res.items() if v["gene"] == "KRAS"}
    secondary_kras = kras_keys_res - kras_keys_par
    with open(os.path.join(OUT, "kras_mutations.csv"), "w", newline="") as f:
        cols = ["sample_status", "chromosome", "start", "end", "ref", "alt",
                "variant_classification", "hgvsp", "hgvsc", "vaf", "depth"]
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader(); w.writerows(kr)
    print(f"KRAS variants: parental={len(kras_keys_par)} resistant={len(kras_keys_res)} "
          f"resistant_private={len(secondary_kras)}")

    # ---- pathway private mutations ----
    pw = [r for r in priv if r["pathway"]]
    with open(os.path.join(OUT, "pathway_private_mutations.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(pw[0].keys()) if pw else ["note"])
        w.writeheader()
        if pw:
            w.writerows(pw)
        else:
            w.writerow({"note": "no private protein-coding PASS variants in MAPK/PI3K/RTK/MYC gene panels"})
    print("pathway private rows:", len(pw))

    # ---- figure: VAF of shared vs private ----
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        sv = [par[k]["vaf"] for k in shared if par[k]["vaf"] is not None]
        rv = [res[k]["vaf"] for k in res_only if res[k]["vaf"] is not None]
        pv = [par[k]["vaf"] for k in par_only if par[k]["vaf"] is not None]
        fig, ax = plt.subplots(figsize=(7, 4.5))
        ax.hist(sv, bins=40, alpha=0.6, label=f"shared (n={len(sv)})")
        ax.hist(pv, bins=40, alpha=0.6, label=f"parental-private (n={len(pv)})")
        ax.hist(rv, bins=40, alpha=0.6, label=f"resistant-private (n={len(rv)})")
        ax.set_xlabel("VAF"); ax.set_ylabel("variants")
        ax.set_title("GSE335219 HPAC P18 vs P28: VAF distributions (PASS coding somatic)")
        ax.legend(); fig.tight_layout()
        fig.savefig(os.path.join(OUT, "vaf_shared_vs_private.png"), dpi=150)
        print("wrote vaf_shared_vs_private.png")
    except ImportError:
        print("matplotlib unavailable, skipping figure")


if __name__ == "__main__":
    main()
