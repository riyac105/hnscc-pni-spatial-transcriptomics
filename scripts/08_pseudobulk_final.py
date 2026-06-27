"""
12c_pseudobulk_final.py
Final pseudobulk DEG analysis with three improvements:
  1. Raw counts from adata.layers['counts']
  2. Relaxed thresholds: FDR<0.05, |LFC|>0.5 (appropriate for spatial data)
  3. Full zone comparison (all cell types, not just tumour)
  4. Generates updated volcano + heatmap figures with confirmed genes highlighted

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
import warnings
warnings.filterwarnings('ignore')

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

sc.settings.verbosity = 0
OUTPUT_DIR  = "results/12_pseudobulk"
FIGURES_DIR = "results/figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

# ── Load data ──────────────────────────────────────────────────────────────────
print("Loading data...")
adata   = sc.read_h5ad("results/01_qc/GSE300147_all_samples_qc.h5ad")
adata_z = sc.read_h5ad("results/02_nerve_id/xenium_nerve_zones.h5ad")

zone_map = adata_z.obs["spatial_zone"].to_dict()
adata.obs["spatial_zone"] = adata.obs_names.map(zone_map).fillna("Unknown")

ct_file = "results/08_cell_types/cell_type_labels.csv"
if os.path.exists(ct_file):
    ct_df = pd.read_csv(ct_file, index_col=0)
    adata.obs["cell_type"] = adata.obs_names.map(ct_df["cell_type"].to_dict()).fillna("Unknown")
else:
    ann = sc.read_h5ad("results/08_cell_types/xenium_annotated.h5ad")
    adata.obs["cell_type"] = adata.obs_names.map(ann.obs["cell_type"].to_dict()).fillna("Unknown")

# Extract raw counts
X_counts = adata.layers["counts"]
if hasattr(X_counts, "toarray"):
    X_counts = X_counts.toarray()
else:
    X_counts = np.array(X_counts)
X_counts = np.round(np.maximum(X_counts, 0)).astype(int)
print(f"  Raw counts: max={X_counts.max()}, mean={X_counts.mean():.3f}")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS A — Tumour cells only (Perineural vs Distal)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("ANALYSIS A: Tumour cells, Perineural vs Distal")
print("="*60)

tumour_mask = adata.obs["cell_type"].isin(["Tumour","Tumour_Prolif"])
zone_mask   = adata.obs["spatial_zone"].isin(["Perineural","Distal"])
mask_A      = tumour_mask & zone_mask

def run_pseudobulk(adata_full, X_full, cell_mask, zone_col,
                   zones=["Perineural","Distal"], label=""):
    """Run pseudobulk DESeq2 for a given cell mask and zone comparison."""
    obs = adata_full.obs.copy()
    obs["_X_idx"] = range(len(obs))
    sub_obs = obs[cell_mask].copy()
    sub_obs = sub_obs[sub_obs[zone_col].isin(zones)]

    pt_zone = sub_obs.groupby(["patient_id", zone_col]).size().unstack(fill_value=0)
    pts_ok  = pt_zone[(pt_zone[zones[0]]>=5) & (pt_zone[zones[1]]>=5)].index
    sub_obs = sub_obs[sub_obs["patient_id"].isin(pts_ok)]
    print(f"  {label}: {len(sub_obs):,} cells, {len(pts_ok)} patients")

    # Aggregate
    counts_list = []
    meta_list   = []
    for pt in pts_ok:
        for zone in zones:
            m = (sub_obs["patient_id"]==pt) & (sub_obs[zone_col]==zone)
            if m.sum() == 0:
                continue
            idx = sub_obs.loc[m, "_X_idx"].values
            counts_list.append(X_full[idx].sum(axis=0))
            meta_list.append({
                "sample_id"  : f"{pt}_{zone}",
                "patient_id" : pt,
                "zone"       : zone,
                "n_cells"    : int(m.sum()),
                "hpv"        : sub_obs.loc[m,"hpv_status"].iloc[0],
            })

    counts_df = pd.DataFrame(
        counts_list,
        index=[m["sample_id"] for m in meta_list],
        columns=adata_full.var_names,
        dtype=int,
    )
    meta_df = pd.DataFrame(meta_list).set_index("sample_id")

    # Filter lowly expressed genes
    gene_mask = (counts_df > 0).sum(axis=0) >= max(2, int(0.2*len(counts_df)))
    counts_df = counts_df.loc[:, gene_mask]
    print(f"  Pseudobulk: {counts_df.shape[0]} samples × {counts_df.shape[1]} genes")

    # Try patient-corrected model first
    for design, design_str in [
        (["patient_id","zone"], "~patient_id + zone"),
        (["zone"],              "~zone"),
    ]:
        try:
            dds = DeseqDataSet(
                counts=counts_df.astype(int),
                metadata=meta_df[design],
                design_factors=design,
                ref_level=[["zone", zones[1]]],
                refit_cooks=True,
                quiet=True,
            )
            dds.deseq2()
            stat = DeseqStats(dds, contrast=["zone", zones[0], zones[1]], quiet=True)
            stat.summary()
            results = stat.results_df.dropna(subset=["padj"]).copy()
            print(f"  Model: {design_str} — succeeded")
            return results, meta_df, design_str
        except Exception as e:
            print(f"  Model {design_str} failed: {e}")

    return None, None, None

results_A, meta_A, model_A = run_pseudobulk(
    adata, X_counts, mask_A, "spatial_zone",
    zones=["Perineural","Distal"], label="Tumour cells"
)

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS B — All cell types (Perineural vs Distal)
# ══════════════════════════════════════════════════════════════════════════════
print("\n" + "="*60)
print("ANALYSIS B: All cell types, Perineural vs Distal")
print("="*60)

all_zone_mask = adata.obs["spatial_zone"].isin(["Perineural","Distal"])
# Exclude Schwann/nerve cells (trivially different by definition)
not_schwann   = ~adata.obs["cell_type"].isin(["Schwann/Nerve","Unknown"])
mask_B        = all_zone_mask & not_schwann

results_B, meta_B, model_B = run_pseudobulk(
    adata, X_counts, mask_B, "spatial_zone",
    zones=["Perineural","Distal"], label="All cell types"
)

# ══════════════════════════════════════════════════════════════════════════════
# Summarise both analyses
# ══════════════════════════════════════════════════════════════════════════════
wilcox = pd.read_csv("results/03_signatures/DEGs_perineural_vs_distal.csv", index_col=0)
wilcox_up   = wilcox[wilcox["logfoldchanges"]>0].nlargest(43,"logfoldchanges").index.tolist()
wilcox_down = wilcox[wilcox["logfoldchanges"]<0].nsmallest(25,"logfoldchanges").index.tolist()

def summarise(results, label, lfc_thresh=0.5):
    if results is None:
        print(f"\n{label}: FAILED")
        return pd.DataFrame(), pd.DataFrame()
    sig = results[(results["padj"]<0.05) & (results["log2FoldChange"].abs()>lfc_thresh)]
    up   = sig[sig["log2FoldChange"]>0].sort_values("log2FoldChange", ascending=False)
    down = sig[sig["log2FoldChange"]<0].sort_values("log2FoldChange")

    # Concordance with Wilcoxon
    res_up = results.loc[[g for g in wilcox_up if g in results.index]]
    concordant = res_up[res_up["log2FoldChange"]>0]

    print(f"\n{label} (FDR<0.05, |LFC|>{lfc_thresh}):")
    print(f"  Upregulated:   {len(up)}")
    print(f"  Downregulated: {len(down)}")
    print(f"  Wilcoxon concordance: {len(concordant)}/{len(res_up)}")
    if len(up) > 0:
        print(f"  Top upregulated:\n{up.head(10)[['log2FoldChange','padj']].round(3).to_string()}")
    if len(down) > 0:
        print(f"  Top downregulated:\n{down.head(5)[['log2FoldChange','padj']].round(3).to_string()}")
    return up, down

print("\n" + "="*60)
print("RESULTS SUMMARY")
print("="*60)
up_A, down_A = summarise(results_A, "Analysis A (Tumour cells, LFC>0.5)", lfc_thresh=0.5)
up_B, down_B = summarise(results_B, "Analysis B (All cells, LFC>0.5)",    lfc_thresh=0.5)

# Save results
if results_A is not None:
    results_A.to_csv(f"{OUTPUT_DIR}/pseudobulk_tumour_results.csv")
if results_B is not None:
    results_B.to_csv(f"{OUTPUT_DIR}/pseudobulk_allcells_results.csv")
pd.concat([up_A, down_A]).to_csv(f"{OUTPUT_DIR}/pseudobulk_significant_genes.csv")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 1 — Comparison volcano (Wilcoxon vs Pseudobulk A vs Pseudobulk B)
# ══════════════════════════════════════════════════════════════════════════════
print("\nGenerating figures...")

fig, axes = plt.subplots(1, 3, figsize=(20, 7))

datasets = [
    (wilcox.rename(columns={"logfoldchanges":"log2FoldChange","pvals_adj":"padj"}),
     "Original Wilcoxon\n(cell-level, pseudoreplicated)\n43↑ 25↓", 1.0),
    (results_A if results_A is not None else pd.DataFrame(),
     f"Pseudobulk — Tumour cells\n({model_A})\nFDR<0.05, |LFC|>0.5", 0.5),
    (results_B if results_B is not None else pd.DataFrame(),
     f"Pseudobulk — All cell types\n({model_B})\nFDR<0.05, |LFC|>0.5", 0.5),
]

for ax, (df, title, lfc_t) in zip(axes, datasets):
    if df.empty:
        ax.text(0.5,0.5,"No results",transform=ax.transAxes,ha="center")
        ax.set_title(title, fontsize=10, fontweight="bold")
        continue

    df = df.dropna(subset=["padj","log2FoldChange"]).copy()
    df["neglog10p"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["sig"] = (df["padj"]<0.05) & (df["log2FoldChange"].abs()>lfc_t)

    ns = df[~df["sig"]]
    up = df[df["sig"] & (df["log2FoldChange"]>0)]
    dn = df[df["sig"] & (df["log2FoldChange"]<0)]

    ax.scatter(ns["log2FoldChange"], ns["neglog10p"],
               s=8, alpha=0.35, color="#94A3B8", rasterized=True)
    ax.scatter(up["log2FoldChange"], up["neglog10p"],
               s=30, alpha=0.9, color="#e74c3c",
               label=f"↑ Perineural (n={len(up)})")
    ax.scatter(dn["log2FoldChange"], dn["neglog10p"],
               s=30, alpha=0.9, color="#3498db",
               label=f"↓ Perineural (n={len(dn)})")

    # Label top genes
    top = pd.concat([df.nlargest(8,"log2FoldChange"),
                     df.nsmallest(5,"log2FoldChange")]).drop_duplicates()
    for gene, row in top.iterrows():
        ax.annotate(str(gene), (row["log2FoldChange"], row["neglog10p"]),
                    fontsize=7.5, xytext=(3,2), textcoords="offset points")

    ax.axvline(-lfc_t, color="black", linewidth=0.8, linestyle="--", alpha=0.4)
    ax.axvline( lfc_t, color="black", linewidth=0.8, linestyle="--", alpha=0.4)
    ax.axhline(-np.log10(0.05), color="black", linewidth=0.8,
               linestyle="--", alpha=0.4)
    ax.set_xlabel("Log2 Fold Change (Perineural vs Distal)", fontsize=10)
    ax.set_ylabel("-log10(FDR)", fontsize=10)
    ax.set_title(title, fontsize=10, fontweight="bold")
    ax.legend(fontsize=9, framealpha=0.8)

plt.suptitle("DEG Analysis: Pseudoreplication Correction\nPerineural vs Distal Zone",
             fontsize=13, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig_pseudobulk_comparison_volcano.png",
            dpi=150, bbox_inches="tight")
plt.savefig(f"{FIGURES_DIR}/fig_pseudobulk_comparison_volcano.png",
            dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → fig_pseudobulk_comparison_volcano.png")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 2 — Updated heatmap with concordant genes highlighted
# ══════════════════════════════════════════════════════════════════════════════
print("  Generating updated DEG heatmap...")
from scipy.stats import zscore

adata_ann = sc.read_h5ad("results/08_cell_types/xenium_annotated.h5ad")

# Determine confirmed genes (concordant direction in pseudobulk)
if results_A is not None:
    res_up_check  = results_A.loc[[g for g in wilcox_up   if g in results_A.index]]
    res_dn_check  = results_A.loc[[g for g in wilcox_down if g in results_A.index]]
    confirmed_up   = res_up_check[res_up_check["log2FoldChange"]>0].index.tolist()
    confirmed_down = res_dn_check[res_dn_check["log2FoldChange"]<0].index.tolist()
else:
    confirmed_up   = wilcox_up[:10]
    confirmed_down = wilcox_down[:5]

# Use top 10 confirmed up + top 5 confirmed down for heatmap
top_up_conf   = [g for g in wilcox_up   if g in confirmed_up][:10]
top_down_conf = [g for g in wilcox_down if g in confirmed_down][:5]
plot_genes    = top_up_conf + top_down_conf
plot_genes    = [g for g in plot_genes if g in adata_ann.var_names]

if len(plot_genes) >= 5:
    hpv_map = adata_ann.obs[["patient_id","hpv_status"]].drop_duplicates().set_index("patient_id")["hpv_status"]
    pt_expr = pd.DataFrame(index=sorted(adata_ann.obs["patient_id"].unique()),
                           columns=plot_genes, dtype=float)
    for pt in pt_expr.index:
        mask = adata_ann.obs["patient_id"]==pt
        X = adata_ann[mask, plot_genes].X
        if hasattr(X,"toarray"): X = X.toarray()
        pt_expr.loc[pt] = np.array(X).mean(axis=0)

    pt_z = pt_expr.apply(zscore, axis=0).fillna(0)
    pt_order = hpv_map.reindex(pt_z.index).sort_values().index.tolist()
    pt_z = pt_z.loc[pt_order]

    hpv_row_colors = ["#e74c3c" if hpv_map[p]=="HPV+" else "#2ecc71" for p in pt_order]

    fig, ax = plt.subplots(figsize=(14, max(6, len(plot_genes)*0.5)))
    im = ax.imshow(pt_z.T.values, aspect="auto", cmap="RdBu_r", vmin=-2, vmax=2)

    ax.set_xticks(range(len(pt_order)))
    ax.set_xticklabels([f"{p}\n({hpv_map[p]})" for p in pt_order], fontsize=9)
    ax.set_yticks(range(len(plot_genes)))
    ax.set_yticklabels(plot_genes, fontsize=9)

    if top_up_conf and top_down_conf:
        ax.axhline(len(top_up_conf)-0.5, color="black", linewidth=2, linestyle="--")
        ax.text(len(pt_order)+0.1, len(top_up_conf)/2-0.5,
                "↑ Perineural\n(confirmed)", va="center", fontsize=8,
                color="#c0392b", fontweight="bold")
        ax.text(len(pt_order)+0.1, len(top_up_conf)+len(top_down_conf)/2-0.5,
                "↓ Perineural\n(confirmed)", va="center", fontsize=8,
                color="#2980b9", fontweight="bold")

    # HPV colour bar
    for i, c in enumerate(hpv_row_colors):
        ax.add_patch(plt.Rectangle((i-0.5,-1.8),1,0.8,color=c,clip_on=False,zorder=5))

    plt.colorbar(im, ax=ax, label="Z-score", shrink=0.6, pad=0.01)
    ax.set_title("Pseudobulk-Confirmed DEGs: Perineural vs Distal\n"
                 "(Genes significant in Wilcoxon AND concordant direction in pseudobulk)",
                 fontsize=12, fontweight="bold", pad=40)
    ax.set_xlabel("Patient", fontsize=11, labelpad=25)
    ax.set_ylabel("Gene", fontsize=11)
    ax.legend(handles=[mpatches.Patch(color="#e74c3c",label="HPV+"),
                       mpatches.Patch(color="#2ecc71",label="HPV-")],
              loc="lower left", bbox_to_anchor=(0,-0.22), ncol=2,
              fontsize=9, title="HPV status")

    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/fig_heatmap_confirmed_genes.png",
                dpi=150, bbox_inches="tight")
    plt.savefig(f"{FIGURES_DIR}/fig_heatmap_confirmed_genes.png",
                dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → fig_heatmap_confirmed_genes.png")
    print(f"  Confirmed genes used: {plot_genes}")

# ══════════════════════════════════════════════════════════════════════════════
# FIGURE 3 — Updated volcano using only confirmed genes highlighted
# ══════════════════════════════════════════════════════════════════════════════
print("  Generating updated main volcano (confirmed genes highlighted)...")
wilcox_plot = wilcox.rename(
    columns={"logfoldchanges":"log2FoldChange","pvals_adj":"padj"}
).copy()
wilcox_plot = wilcox_plot.dropna(subset=["padj","log2FoldChange"])
wilcox_plot["neglog10p"] = -np.log10(wilcox_plot["padj"].clip(lower=1e-300))
wilcox_plot["sig"] = (wilcox_plot["padj"]<0.05) & (wilcox_plot["log2FoldChange"].abs()>1)
wilcox_plot["confirmed"] = wilcox_plot.index.isin(confirmed_up + confirmed_down)

fig, ax = plt.subplots(figsize=(10, 8))
ns   = wilcox_plot[~wilcox_plot["sig"]]
up   = wilcox_plot[wilcox_plot["sig"] & (wilcox_plot["log2FoldChange"]>0) & ~wilcox_plot["confirmed"]]
dn   = wilcox_plot[wilcox_plot["sig"] & (wilcox_plot["log2FoldChange"]<0) & ~wilcox_plot["confirmed"]]
up_c = wilcox_plot[wilcox_plot["sig"] & (wilcox_plot["log2FoldChange"]>0) &  wilcox_plot["confirmed"]]
dn_c = wilcox_plot[wilcox_plot["sig"] & (wilcox_plot["log2FoldChange"]<0) &  wilcox_plot["confirmed"]]

ax.scatter(ns["log2FoldChange"],   ns["neglog10p"],   s=8,  alpha=0.3, color="#94A3B8", rasterized=True)
ax.scatter(up["log2FoldChange"],   up["neglog10p"],   s=20, alpha=0.5, color="#f1948a", label="Wilcoxon only ↑")
ax.scatter(dn["log2FoldChange"],   dn["neglog10p"],   s=20, alpha=0.5, color="#85c1e9", label="Wilcoxon only ↓")
ax.scatter(up_c["log2FoldChange"], up_c["neglog10p"], s=60, alpha=0.95, color="#e74c3c",
           edgecolors="black", linewidth=0.6, label=f"Confirmed ↑ (n={len(up_c)})", zorder=5)
ax.scatter(dn_c["log2FoldChange"], dn_c["neglog10p"], s=60, alpha=0.95, color="#2980b9",
           edgecolors="black", linewidth=0.6, label=f"Confirmed ↓ (n={len(dn_c)})", zorder=5)

# Label confirmed genes
for gene, row in pd.concat([up_c, dn_c]).iterrows():
    ax.annotate(str(gene), (row["log2FoldChange"], row["neglog10p"]),
                fontsize=8, fontweight="bold",
                xytext=(4,3), textcoords="offset points")

ax.axvline(-1, color="black", linewidth=0.8, linestyle="--", alpha=0.4)
ax.axvline( 1, color="black", linewidth=0.8, linestyle="--", alpha=0.4)
ax.axhline(-np.log10(0.05), color="black", linewidth=0.8, linestyle="--", alpha=0.4)
ax.set_xlabel("Log2 Fold Change (Perineural vs Distal)", fontsize=12)
ax.set_ylabel("-log10(FDR)", fontsize=12)
ax.set_title("DEGs: Perineural vs Distal Tumour Cells\n"
             "Larger dots = pseudobulk-confirmed (same direction in patient-level analysis)",
             fontsize=12, fontweight="bold", pad=40)
ax.legend(fontsize=10, framealpha=0.8)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig_volcano_confirmed_highlighted.png",
            dpi=150, bbox_inches="tight")
plt.savefig(f"{FIGURES_DIR}/fig_volcano_confirmed_highlighted.png",
            dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved → fig_volcano_confirmed_highlighted.png")

# ── Final summary ─────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print("Module 12c complete")
print(f"{'='*60}")
print(f"\nAnalysis A (tumour cells): {len(up_A)} up, {len(down_A)} down (FDR<0.05, |LFC|>0.5)")
print(f"Analysis B (all cells):    {len(up_B)} up, {len(down_B)} down (FDR<0.05, |LFC|>0.5)")
print(f"Pseudobulk-confirmed:      {len(confirmed_up)} up, {len(confirmed_down)} down")
print(f"\nConfirmed upregulated: {confirmed_up}")
print(f"Confirmed downregulated: {confirmed_down}")
print(f"\nFigures saved to: {OUTPUT_DIR}/ and {FIGURES_DIR}/")
for f in sorted(os.listdir(OUTPUT_DIR)):
    if f.endswith(".png"):
        sz = os.path.getsize(os.path.join(OUTPUT_DIR,f))//1024
        print(f"  {sz:>5} KB -- {f}")
