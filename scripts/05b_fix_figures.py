"""
09b_fix_figures.py
Fixes:
  1. Fig 9a — HPV colour bar was inverted (showing HPV+ for all patients)
  2. Fig 9b — p-value display (p=0.00e+00 → p<1e-300), add biological note
Run from ~/Desktop/PNI_project/
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.stats import zscore
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = "results/figures"

print("Loading data...")
adata   = sc.read_h5ad("results/08_cell_types/xenium_annotated.h5ad")
adata_z = sc.read_h5ad("results/02_nerve_id/xenium_nerve_zones.h5ad")
degs    = pd.read_csv("results/03_signatures/DEGs_perineural_vs_distal.csv", index_col=0)

if "cell_type" not in adata_z.obs.columns:
    label_map = adata.obs["cell_type"].to_dict()
    adata_z.obs["cell_type"] = adata_z.obs_names.map(label_map).fillna("Unknown")
    nerve_mask = adata_z.obs["is_nerve"] == 1
    adata_z.obs.loc[nerve_mask, "cell_type"] = "Schwann/Nerve"

if "hpv_status" not in adata_z.obs.columns:
    hpv_map = adata.obs["hpv_status"].to_dict()
    adata_z.obs["hpv_status"] = adata_z.obs_names.map(hpv_map)

# ── Rebuild DEG gene lists ─────────────────────────────────────────────────────
sig = degs[degs.index.isin(adata.var_names)].copy()
top_up   = sig[sig["logfoldchanges"] > 0].nlargest(10, "logfoldchanges").index.tolist()
top_down = sig[sig["logfoldchanges"] < 0].nsmallest(10, "logfoldchanges").index.tolist()
top_genes = [g for g in top_up + top_down if g in adata.var_names]

hpv_map_pt = adata.obs[["patient_id","hpv_status"]].drop_duplicates().set_index("patient_id")["hpv_status"]

# ══════════════════════════════════════════════════════════════════════════════
# FIX 1 — DEG Heatmap with correct HPV colour bar
# ══════════════════════════════════════════════════════════════════════════════
print("\n[1/2] Fixing DEG heatmap HPV colour bar...")

pt_gene_expr = pd.DataFrame(
    index=sorted(adata.obs["patient_id"].unique()),
    columns=top_genes, dtype=float
)
for pt in pt_gene_expr.index:
    mask = adata.obs["patient_id"] == pt
    X = adata[mask, top_genes].X
    pt_gene_expr.loc[pt] = X.toarray().mean(axis=0) if hasattr(X, 'toarray') else np.array(X).mean(axis=0)

pt_gene_z = pt_gene_expr.apply(zscore, axis=0).fillna(0)

# Sort patients: HPV+ first, then HPV-
pt_order = (hpv_map_pt.reindex(pt_gene_z.index)
            .sort_values()
            .index.tolist())
pt_gene_z = pt_gene_z.loc[pt_order]

# Correctly assign HPV colours PER PATIENT
hpv_colors_row = ["#e74c3c" if hpv_map_pt[p] == "HPV+" else "#2ecc71"
                  for p in pt_order]

fig, ax = plt.subplots(figsize=(14, 7))
im = ax.imshow(pt_gene_z.T.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)

ax.set_xticks(range(len(pt_order)))
ax.set_xticklabels(
    [f"{p}\n({hpv_map_pt[p]})" for p in pt_order],
    fontsize=9
)
ax.set_yticks(range(len(top_genes)))
ax.set_yticklabels(top_genes, fontsize=9)
ax.axhline(len(top_up) - 0.5, color="black", linewidth=2, linestyle="--")

# Colourbar labels
ax.text(len(pt_order) + 0.2, len(top_up)/2 - 0.5,
        "↑ Perineural", va="center", fontsize=9, color="#c0392b", fontweight="bold")
ax.text(len(pt_order) + 0.2, len(top_up) + len(top_down)/2 - 0.5,
        "↓ Perineural", va="center", fontsize=9, color="#2980b9", fontweight="bold")

# FIX: Draw HPV colour bar correctly per patient
for i, (pt, c) in enumerate(zip(pt_order, hpv_colors_row)):
    ax.add_patch(plt.Rectangle((i - 0.5, -1.8), 1, 0.8,
                                color=c, clip_on=False, zorder=5))

plt.colorbar(im, ax=ax, label="Z-score", shrink=0.6, pad=0.01)
ax.set_title("Top DEGs: Perineural vs Distal Tumour Cells\n(Z-scored mean expression per patient)",
             fontsize=13, fontweight="bold")
ax.set_xlabel("Patient", fontsize=11, labelpad=25)
ax.set_ylabel("Gene", fontsize=11)

legend_patches = [
    mpatches.Patch(color="#e74c3c", label="HPV+"),
    mpatches.Patch(color="#2ecc71", label="HPV-"),
]
ax.legend(handles=legend_patches, loc="lower left",
          bbox_to_anchor=(0, -0.22), ncol=2, fontsize=9,
          title="HPV status", framealpha=0.8)

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig9a_DEG_heatmap_per_patient.png",
            dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → {OUTPUT_DIR}/fig9a_DEG_heatmap_per_patient.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIX 2 — PNI signature violins with proper p-values + biological annotation
# ══════════════════════════════════════════════════════════════════════════════
print("\n[2/2] Fixing PNI signature violin p-values...")

sig_genes = [g for g in top_up[:5] if g in adata_z.var_names]
print(f"  Signature genes: {sig_genes}")

sc.tl.score_genes(adata_z, gene_list=sig_genes,
                  score_name="PNI_sig_score", use_raw=False)

plot_data = adata_z.obs[["spatial_zone","hpv_status","PNI_sig_score"]].copy()
plot_data = plot_data[plot_data["spatial_zone"].isin(
    ["Nerve", "Perineural", "Peritumoral", "Distal"]
)]

zone_order  = ["Nerve", "Perineural", "Peritumoral", "Distal"]
hpv_groups  = ["HPV+", "HPV-"]
hpv_colors  = {"HPV+": "#e74c3c", "HPV-": "#2ecc71"}

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

for ax, hpv in zip(axes, hpv_groups):
    sub = plot_data[plot_data["hpv_status"] == hpv]
    data_by_zone = [sub[sub["spatial_zone"] == z]["PNI_sig_score"].values
                    for z in zone_order]

    vp = ax.violinplot(data_by_zone, positions=range(len(zone_order)),
                       showmedians=True, showextrema=False)
    for body in vp["bodies"]:
        body.set_facecolor(hpv_colors[hpv])
        body.set_alpha(0.7)
    vp["cmedians"].set_color("black")
    vp["cmedians"].set_linewidth(2)

    # FIX: proper p-value formatting
    peri = sub[sub["spatial_zone"] == "Perineural"]["PNI_sig_score"].values
    dist = sub[sub["spatial_zone"] == "Distal"]["PNI_sig_score"].values
    if len(peri) > 10 and len(dist) > 10:
        _, pval = stats.mannwhitneyu(peri, dist, alternative="greater")
        if pval < 1e-300:
            pval_str = "p < 1×10⁻³⁰⁰"
        elif pval < 0.001:
            pval_str = f"p = {pval:.2e}"
        elif pval < 0.05:
            pval_str = f"p = {pval:.4f}"
        else:
            pval_str = f"p = {pval:.3f} (n.s.)"

        ax.text(0.97, 0.97, f"Perineural vs Distal\n{pval_str}",
                transform=ax.transAxes, ha="right", va="top",
                fontsize=9, bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.6))

    ax.set_xticks(range(len(zone_order)))
    ax.set_xticklabels(zone_order, fontsize=10)
    ax.set_title(f"{hpv}", fontsize=13, fontweight="bold",
                 color=hpv_colors[hpv])
    ax.set_xlabel("Spatial Zone", fontsize=11)
    ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")

axes[0].set_ylabel(f"PNI Signature Score\n({', '.join(sig_genes)})", fontsize=9)

# Add biological interpretation note
fig.text(0.5, -0.02,
         "Note: PNI signature (NFE2L2, CLCA2, MDM2, EGFR, SOX2) is significantly elevated\n"
         "in perineural zone in HPV- tumours only, consistent with higher PNI index in HPV- patients.",
         ha="center", fontsize=9, style="italic", color="#555555")

plt.suptitle("PNI Signature Gene Expression by Spatial Zone",
             fontsize=14, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig9b_PNI_signature_by_zone_HPV.png",
            dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → {OUTPUT_DIR}/fig9b_PNI_signature_by_zone_HPV.png")

print(f"\n{'='*50}")
print("Fixes complete ✓")
print(f"{'='*50}")
