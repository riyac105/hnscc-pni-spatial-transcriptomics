"""
09_additional_figures.py
Generates 4 additional figures:
  1. DEG heatmap (top 20) across all patients
  2. PNI signature score per zone per patient (violin) split by HPV
  3. Nerve density vs PNI index scatter
  4. Cell type zone enrichment with odds ratios + chi-square stats
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
from scipy.stats import chi2_contingency
import warnings
warnings.filterwarnings('ignore')

OUTPUT_DIR = "results/figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading data...")
adata      = sc.read_h5ad("results/08_cell_types/xenium_annotated.h5ad")
adata_z    = sc.read_h5ad("results/02_nerve_id/xenium_nerve_zones.h5ad")
degs       = pd.read_csv("results/03_signatures/DEGs_perineural_vs_distal.csv", index_col=0)

print(f"  Main object:  {adata.shape}")
print(f"  Zones object: {adata_z.shape}")
print(f"  DEGs: {len(degs)}")

# Transfer cell_type and hpv to zones object if not present
if "cell_type" not in adata_z.obs.columns:
    label_map = adata.obs["cell_type"].to_dict()
    adata_z.obs["cell_type"] = adata_z.obs_names.map(label_map).fillna("Unknown")
    nerve_mask = adata_z.obs["is_nerve"] == 1
    adata_z.obs.loc[nerve_mask, "cell_type"] = "Schwann/Nerve"

if "hpv_status" not in adata_z.obs.columns:
    hpv_map = adata.obs["hpv_status"].to_dict()
    adata_z.obs["hpv_status"] = adata_z.obs_names.map(hpv_map)

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — DEG Heatmap: top 20 genes × all patients
# ══════════════════════════════════════════════════════════════════════════════
print("\n[1/4] DEG heatmap across patients...")

# Pick top 10 up + top 10 down by LFC (confirmed in panel)
if "logfoldchanges" in degs.columns:
    lfc_col = "logfoldchanges"
elif "log2FoldChange" in degs.columns:
    lfc_col = "log2FoldChange"
else:
    lfc_col = degs.columns[0]
    print(f"  Using column: {lfc_col}")

print(f"  DEG columns: {degs.columns.tolist()}")

sig = degs[degs.index.isin(adata.var_names)].copy()
top_up   = sig[sig[lfc_col] > 0].nlargest(10, lfc_col).index.tolist()
top_down = sig[sig[lfc_col] < 0].nsmallest(10, lfc_col).index.tolist()
top_genes = top_up + top_down
top_genes = [g for g in top_genes if g in adata.var_names]
print(f"  Top genes selected: {top_genes}")

# Mean expression per patient per gene (log-normalised)
sc.pp.normalize_total(adata, target_sum=1e4, inplace=True) if "log1p" not in adata.uns else None
pt_gene_expr = pd.DataFrame(
    index=sorted(adata.obs["patient_id"].unique()),
    columns=top_genes,
    dtype=float
)
for pt in pt_gene_expr.index:
    mask = adata.obs["patient_id"] == pt
    pt_gene_expr.loc[pt] = adata[mask, top_genes].X.toarray().mean(axis=0) \
        if hasattr(adata[mask, top_genes].X, 'toarray') \
        else adata[mask, top_genes].X.mean(axis=0)

# Z-score across patients per gene
from scipy.stats import zscore
pt_gene_z = pt_gene_expr.apply(zscore, axis=0).fillna(0)

# Sort patients by HPV status
hpv_map_pt = adata.obs[["patient_id","hpv_status"]].drop_duplicates().set_index("patient_id")["hpv_status"]
pt_order = (hpv_map_pt.reindex(pt_gene_z.index)
            .sort_values()
            .index.tolist())
pt_gene_z = pt_gene_z.loc[pt_order]

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
ax.text(len(pt_order) + 0.1, len(top_up)/2 - 0.5, "↑ Perineural",
        va="center", fontsize=9, color="#c0392b", fontweight="bold")
ax.text(len(pt_order) + 0.1, len(top_up) + len(top_down)/2 - 0.5, "↓ Perineural",
        va="center", fontsize=9, color="#2980b9", fontweight="bold")

# HPV status colour bar on top
hpv_colors_row = ["#e74c3c" if hpv_map_pt[p] == "HPV+" else "#2ecc71" for p in pt_order]
for i, c in enumerate(hpv_colors_row):
    ax.add_patch(plt.Rectangle((i-0.5, -1.5), 1, 0.8, color=c, clip_on=False))

plt.colorbar(im, ax=ax, label="Z-score", shrink=0.6, pad=0.01)
ax.set_title("Top DEGs: Perineural vs Distal Tumour Cells\n(Z-scored mean expression per patient)",
             fontsize=13, fontweight="bold")
ax.set_xlabel("Patient", fontsize=11)
ax.set_ylabel("Gene", fontsize=11)

legend_patches = [
    mpatches.Patch(color="#e74c3c", label="HPV+"),
    mpatches.Patch(color="#2ecc71", label="HPV-"),
]
ax.legend(handles=legend_patches, loc="upper left",
          bbox_to_anchor=(0, -0.15), ncol=2, fontsize=9, title="HPV status")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig9a_DEG_heatmap_per_patient.png",
            dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → {OUTPUT_DIR}/fig9a_DEG_heatmap_per_patient.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — PNI signature gene scores by zone × HPV status
# ══════════════════════════════════════════════════════════════════════════════
print("\n[2/4] PNI signature scores by zone and HPV...")

# Use top 5 upregulated genes as PNI signature
sig_genes = [g for g in top_up[:5] if g in adata_z.var_names]
print(f"  Signature genes: {sig_genes}")

if sig_genes:
    sc.tl.score_genes(adata_z, gene_list=sig_genes,
                      score_name="PNI_sig_score", use_raw=False)

    plot_data = adata_z.obs[["spatial_zone", "hpv_status", "PNI_sig_score"]].copy()
    plot_data = plot_data[plot_data["spatial_zone"].isin(
        ["Nerve", "Perineural", "Peritumoral", "Distal"]
    )]

    zone_order = ["Nerve", "Perineural", "Peritumoral", "Distal"]
    hpv_groups = ["HPV+", "HPV-"]
    hpv_colors = {"HPV+": "#e74c3c", "HPV-": "#2ecc71"}

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

        # Stats: Perineural vs Distal
        peri = sub[sub["spatial_zone"] == "Perineural"]["PNI_sig_score"].values
        dist = sub[sub["spatial_zone"] == "Distal"]["PNI_sig_score"].values
        if len(peri) > 10 and len(dist) > 10:
            _, pval = stats.mannwhitneyu(peri, dist, alternative="greater")
            ax.text(0.5, 0.95, f"Perineural vs Distal\np={pval:.2e}",
                    transform=ax.transAxes, ha="center", va="top",
                    fontsize=9, bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

        ax.set_xticks(range(len(zone_order)))
        ax.set_xticklabels(zone_order, fontsize=10)
        ax.set_title(f"{hpv}", fontsize=13, fontweight="bold",
                     color=hpv_colors[hpv])
        ax.set_xlabel("Spatial Zone", fontsize=11)
        ax.axhline(0, color="gray", linewidth=0.8, linestyle="--")

    axes[0].set_ylabel(f"PNI Signature Score\n({', '.join(sig_genes)})", fontsize=10)
    plt.suptitle("PNI Signature Gene Expression by Spatial Zone",
                 fontsize=14, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/fig9b_PNI_signature_by_zone_HPV.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {OUTPUT_DIR}/fig9b_PNI_signature_by_zone_HPV.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Nerve density vs PNI index scatter
# ══════════════════════════════════════════════════════════════════════════════
print("\n[3/4] Nerve density vs PNI index scatter...")

by_pt = []
for pt in adata_z.obs["patient_id"].unique():
    mask = adata_z.obs["patient_id"] == pt
    sub  = adata_z.obs[mask]
    n_nerve   = (sub["is_nerve"] == 1).sum()
    n_pni_pos = (sub["PNI_positive"] == 1).sum()
    hpv = sub["hpv_status"].iloc[0]
    by_pt.append({
        "patient"    : pt,
        "hpv"        : hpv,
        "pct_nerve"  : 100 * n_nerve / len(sub),
        "pni_index"  : n_pni_pos / max(n_nerve, 1),
    })
df = pd.DataFrame(by_pt)

hpv_colors = {"HPV+": "#e74c3c", "HPV-": "#2ecc71"}
fig, ax = plt.subplots(figsize=(7, 6))

for _, row in df.iterrows():
    ax.scatter(row["pct_nerve"], row["pni_index"],
               color=hpv_colors[row["hpv"]], s=120,
               edgecolors="black", linewidth=0.8, zorder=5)
    ax.annotate(row["patient"],
                (row["pct_nerve"], row["pni_index"]),
                textcoords="offset points", xytext=(6, 4), fontsize=9)

# Pearson correlation
r, pval = stats.pearsonr(df["pct_nerve"], df["pni_index"])
m, b = np.polyfit(df["pct_nerve"], df["pni_index"], 1)
x_line = np.linspace(df["pct_nerve"].min(), df["pct_nerve"].max(), 100)
ax.plot(x_line, m*x_line + b, "k--", linewidth=1.5, alpha=0.6)
ax.text(0.05, 0.92,
        f"Pearson r = {r:.2f}\np = {pval:.3f}",
        transform=ax.transAxes, fontsize=10,
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5))

legend_patches = [mpatches.Patch(color=c, label=h) for h, c in hpv_colors.items()]
ax.legend(handles=legend_patches, fontsize=10, framealpha=0.8)
ax.set_xlabel("Nerve Cell Density (% of total cells)", fontsize=12)
ax.set_ylabel("PNI Index\n(perineural tumour cells per nerve cell)", fontsize=11)
ax.set_title("Nerve Density vs Perineural Invasion Index\nper Patient",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig9c_nerve_density_vs_PNI_index.png",
            dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → {OUTPUT_DIR}/fig9c_nerve_density_vs_PNI_index.png")


# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 4 — Cell type zone enrichment: odds ratios + chi-square
# ══════════════════════════════════════════════════════════════════════════════
print("\n[4/4] Cell type zone enrichment statistics...")

CT_COLORS = {
    'Tumour'        : '#e74c3c',
    'Tumour_Prolif' : '#ff5722',
    'Schwann/Nerve' : '#9b59b6',
    'T_cell'        : '#3498db',
    'B_cell'        : '#1abc9c',
    'DC_LAMP3'      : '#F39C12',
    'Macrophage'    : '#e67e22',
    'Fibroblast'    : '#795548',
    'Endothelial'   : '#e91e63',
    'Mast'          : '#607d8b',
}

focus_zone  = "Perineural"
ref_zone    = "Distal"
cell_types  = [ct for ct in CT_COLORS if ct != "Schwann/Nerve"]  # exclude nerve (trivially enriched in Nerve zone)

zone_ct = adata_z.obs.groupby(["spatial_zone","cell_type"]).size().unstack(fill_value=0)

results = []
for ct in cell_types:
    if ct not in zone_ct.columns:
        continue
    # 2×2 contingency: ct vs not-ct, perineural vs distal
    if focus_zone not in zone_ct.index or ref_zone not in zone_ct.index:
        continue
    a = zone_ct.loc[focus_zone, ct]                                    # ct in perineural
    b = zone_ct.loc[focus_zone].sum() - a                              # not-ct in perineural
    c = zone_ct.loc[ref_zone, ct]                                      # ct in distal
    d = zone_ct.loc[ref_zone].sum() - c                                # not-ct in distal

    if c == 0 or b == 0 or d == 0:
        continue

    # Odds ratio
    OR = (a / b) / (c / d)
    # 95% CI using log method
    log_or = np.log(OR)
    se = np.sqrt(1/a + 1/b + 1/c + 1/d) if a > 0 else np.nan
    ci_low  = np.exp(log_or - 1.96*se) if a > 0 else np.nan
    ci_high = np.exp(log_or + 1.96*se) if a > 0 else np.nan

    # Chi-square
    chi2, pval, _, _ = chi2_contingency([[a, b], [c, d]])

    # Bonferroni correction
    pval_adj = min(pval * len(cell_types), 1.0)

    results.append({
        "cell_type" : ct,
        "OR"        : OR,
        "CI_low"    : ci_low,
        "CI_high"   : ci_high,
        "chi2"      : chi2,
        "pval"      : pval,
        "pval_adj"  : pval_adj,
        "n_perineural": a,
        "n_distal"  : c,
        "enriched"  : OR > 1,
    })

res_df = pd.DataFrame(results).sort_values("OR", ascending=True)
print("\nCell type enrichment in Perineural vs Distal zone:")
print(res_df[["cell_type","OR","CI_low","CI_high","pval","pval_adj","n_perineural","n_distal"]].round(3).to_string())

# Save stats table
res_df.to_csv(f"{OUTPUT_DIR}/cell_type_enrichment_stats.csv", index=False)
print(f"\n  Saved table → {OUTPUT_DIR}/cell_type_enrichment_stats.csv")

# Forest plot
fig, ax = plt.subplots(figsize=(9, 7))
y_pos = range(len(res_df))

for i, (_, row) in enumerate(res_df.iterrows()):
    color = CT_COLORS.get(row["cell_type"], "#bdc3c7")
    sig   = "**" if row["pval_adj"] < 0.01 else ("*" if row["pval_adj"] < 0.05 else "")
    alpha = 1.0 if row["pval_adj"] < 0.05 else 0.45

    # CI line
    ax.plot([row["CI_low"], row["CI_high"]], [i, i],
            color=color, linewidth=2, alpha=alpha)
    # Point estimate
    ax.scatter(row["OR"], i, color=color, s=120,
               zorder=5, alpha=alpha, edgecolors="black", linewidth=0.5)
    # Label
    label = f"{row['cell_type']}  OR={row['OR']:.2f} {sig}"
    ax.text(row["CI_high"] + 0.05, i, label,
            va="center", fontsize=9, alpha=alpha)

ax.axvline(1.0, color="black", linewidth=1.5, linestyle="--", alpha=0.7)
ax.set_yticks([])
ax.set_xlabel("Odds Ratio (Perineural vs Distal)", fontsize=12)
ax.set_title(f"Cell Type Enrichment in Perineural Zone\n"
             f"(vs Distal, * p_adj<0.05, ** p_adj<0.01, Bonferroni corrected)",
             fontsize=12, fontweight="bold")
ax.set_xlim(left=0)

# Shade enriched vs depleted
ax.axvspan(0, 1, alpha=0.04, color="blue", label="Depleted in Perineural")
ax.axvspan(1, ax.get_xlim()[1], alpha=0.04, color="red", label="Enriched in Perineural")
ax.legend(fontsize=9, loc="lower right")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig9d_cell_type_enrichment_forest.png",
            dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → {OUTPUT_DIR}/fig9d_cell_type_enrichment_forest.png")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("Module 09 complete ✓")
print(f"{'='*60}")
for f in sorted(os.listdir(OUTPUT_DIR)):
    if "fig9" in f or "enrichment" in f:
        size = os.path.getsize(os.path.join(OUTPUT_DIR, f)) // 1024
        print(f"  {size:>5} KB — {f}")
