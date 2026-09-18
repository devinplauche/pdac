#!/usr/bin/env python3
"""Figures for the GSE337208 CAF mining analysis."""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "output")
SIGS = ["myCAF", "iCAF", "apCAF", "SASP", "matrisome_core",
        "MAPK_targets", "RASON_inhibition_PD"]
LABELS = {"myCAF": "myCAF", "iCAF": "iCAF", "apCAF": "apCAF", "SASP": "SASP",
          "matrisome_core": "matrisome", "MAPK_targets": "MAPK output",
          "RASON_inhibition_PD": "RAS(ON) PD"}

# --- Fig 1: cell-level score distributions, KPC2 Veh vs RMC6236 ---
cells = pd.read_csv(os.path.join(OUT, "kpc2_cell_signature_scores.csv"))
stats = pd.read_csv(os.path.join(OUT, "kpc2_contrast_stats.csv")
                    ).set_index("signature")

fig, axes = plt.subplots(1, 7, figsize=(14, 4.2), sharey=True)
for ax, s in zip(axes, SIGS):
    a = cells.loc[cells.treatment == "Veh", s].dropna()
    b = cells.loc[cells.treatment == "RMC6236", s].dropna()
    ax.boxplot([a, b], labels=["Veh\n(n=42)", "RMC6236\n(n=53)"],
               patch_artist=True,
               boxprops=dict(facecolor="#dbe9f6"),
               medianprops=dict(color="black"))
    p = stats.loc[s, "MWU_p_cell_level"]
    ax.set_title(f"{LABELS[s]}\nMWU p={p:.4f}", fontsize=10)
    ax.tick_params(labelsize=9)
axes[0].set_ylabel("signature score\n(mean fractional rank / cell)")
fig.suptitle("GSE337208 KPC2: fibroblast signature scores, RMC6236 vs vehicle\n"
             "cell-level (pseudoreplicated: 1 sample/arm)", fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig1_kpc2_cell_scores.png"), dpi=150)
plt.close(fig)

# --- Fig 2: pseudobulk z-scores across samples ---
pb = pd.read_csv(os.path.join(OUT, "pseudobulk_scores_with_meta.csv"),
                 index_col="sample")
order = ["GSM9850794", "GSM9850795", "GSM9850796", "GSM9850797",
         "GSM9850798", "GSM9850799", "GSM9850800"]
pb = pb.loc[order]
labels = [f"{r['model']} {r['treatment']}\n(n={int(r['n_fibroblasts'])} fib)"
          for _, r in pb.iterrows()]

fig, ax = plt.subplots(figsize=(9, 5.5))
v = pb[SIGS].values
im = ax.imshow(v, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5)
ax.set_xticks(range(len(SIGS)))
ax.set_xticklabels([LABELS[s] for s in SIGS], rotation=30, ha="right")
ax.set_yticks(range(len(pb)))
ax.set_yticklabels(labels, fontsize=9)
for i in range(v.shape[0]):
    for j in range(v.shape[1]):
        ax.text(j, i, f"{v[i, j]:.2f}", ha="center", va="center",
                fontsize=8, color="white" if abs(v[i, j]) > 0.9 else "black")
# model separators
for i in [1.5, 4.5]:
    ax.axhline(i, color="black", linewidth=1.5)
fig.colorbar(im, ax=ax, label="z-score across samples")
ax.set_title("Fibroblast pseudobulk signature z-scores (7 samples)\n"
             "z-scored across samples: compare within KPC model only",
             fontsize=11)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "fig2_pseudobulk_heatmap.png"), dpi=150)
plt.close(fig)
print("figures written")
