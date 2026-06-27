"""
12d_pseudobulk_HPVneg.py
Pseudobulk DEG analysis restricted to HPV- patients only.

Rationale: The PNI transcriptional signature is HPV--specific (violin plot).
Running DEGs on all 10 patients dilutes the HPV- signal with flat HPV+ results.
Testing HPV- only directly addresses the biological hypothesis.

Three analyses:
  A. HPV- tumour cells only (n=5 patients)
  B. HPV- all cell types (n=5 patients, more power)
  C. Bonus: HPV+ tumour cells only (expect few/no DEGs - negative control)

Run from ~/Desktop/PNI_project/
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

sc.settings.verbosity = 0
OUTPUT_DIR = "results/12_pseudobulk"
os.makedirs(OUTPUT_DIR, exist_ok=True)

HPV_NEG = ['P12','P15','P17','P20','P28']
HPV_POS = ['P1','P9','P13','P19','P23']

# ── Load ──────────────────────────────────────────────────────────────────────
print("Loading data...")
adata   = sc.read_h5ad("results/01_qc/GSE300147_all_samples_qc.h5ad")
adata_z = sc.read_h5ad("results/02_nerve_id/xenium_nerve_zones.h5ad")

zone_map = adata_z.obs["spatial_zone"].to_dict()
adata.obs["spatial_zone"] = adata.obs_names.map(zone_map).fillna("Unknown")

ct_file = "results/08_cell_types/cell_type_labels.csv"
ct_df   = pd.read_csv(ct_file, index_col=0)
adata.obs["cell_type"] = adata.obs_names.map(ct_df["cell_type"].to_dict()).fillna("Unknown")

# Raw counts
X_counts = adata.layers["counts"]
if hasattr(X_counts,"toarray"): X_counts = X_counts.toarray()
X_counts = np.round(np.maximum(np.array(X_counts), 0)).astype(int)
print(f"  {adata.shape}, raw counts max={X_counts.max()}")

def run_pseudobulk_subset(patient_list, cell_types, label, lfc_thresh=0.5):
    """Run pseudobulk DESeq2 for a subset of patients and cell types."""
    print(f"\n{'='*55}")
    print(f"{label}")
    print(f"  Patients: {patient_list}")
    print(f"  Cell types: {cell_types}")
    print(f"{'='*55}")

    mask = (adata.obs["patient_id"].isin(patient_list) &
            adata.obs["cell_type"].isin(cell_types) &
            adata.obs["spatial_zone"].isin(["Perineural","Distal"]))

    sub_obs = adata.obs[mask].copy()
    sub_obs["_idx"] = np.where(mask)[0]

    pt_zone = sub_obs.groupby(["patient_id","spatial_zone"]).size().unstack(fill_value=0)
    pts_ok  = pt_zone[(pt_zone.get("Perineural",0)>=5) &
                      (pt_zone.get("Distal",0)>=5)].index.tolist()
    sub_obs = sub_obs[sub_obs["patient_id"].isin(pts_ok)]
    print(f"  {len(sub_obs):,} cells, {len(pts_ok)} patients with data in both zones")

    if len(pts_ok) < 3:
        print("  Too few patients — skipping")
        return None

    # Aggregate
    counts_list, meta_list = [], []
    for pt in pts_ok:
        for zone in ["Perineural","Distal"]:
            m = (sub_obs["patient_id"]==pt) & (sub_obs["spatial_zone"]==zone)
            if m.sum() == 0: continue
            idx = sub_obs.loc[m,"_idx"].values
            counts_list.append(X_counts[idx].sum(axis=0))
            meta_list.append({"sample_id":f"{pt}_{zone}",
                               "patient_id":pt, "zone":zone,
                               "n_cells":int(m.sum())})

    counts_df = pd.DataFrame(counts_list,
                              index=[m["sample_id"] for m in meta_list],
                              columns=adata.var_names, dtype=int)
    meta_df   = pd.DataFrame(meta_list).set_index("sample_id")

    # Filter low-expressed genes
    gene_mask = (counts_df>0).sum(axis=0) >= max(2, int(0.2*len(counts_df)))
    counts_df = counts_df.loc[:, gene_mask]
    print(f"  Pseudobulk: {counts_df.shape[0]} samples × {counts_df.shape[1]} genes")
    print(f"  Count range: {counts_df.values.min()} – {counts_df.values.max()}")

    # Model selection
    design = ["patient_id","zone"] if len(pts_ok) >= 4 else ["zone"]
    design_str = "~patient_id + zone" if len(pts_ok)>=4 else "~zone (no patient correction, n<4)"
    print(f"  Model: {design_str}")

    try:
        dds = DeseqDataSet(
            counts=counts_df.astype(int),
            metadata=meta_df[design],
            design_factors=design,
            ref_level=[["zone","Distal"]],
            refit_cooks=True, quiet=True,
        )
        dds.deseq2()
        stat = DeseqStats(dds, contrast=["zone","Perineural","Distal"], quiet=True)
        stat.summary()
        results = stat.results_df.dropna(subset=["padj"]).copy()
    except Exception as e:
        print(f"  Model failed: {e}")
        return None

    sig = results[(results["padj"]<0.05) & (results["log2FoldChange"].abs()>lfc_thresh)]
    up   = sig[sig["log2FoldChange"]>0].sort_values("log2FoldChange", ascending=False)
    down = sig[sig["log2FoldChange"]<0].sort_values("log2FoldChange")

    print(f"\n  RESULTS (FDR<0.05, |LFC|>{lfc_thresh}):")
    print(f"    Upregulated:   {len(up)}")
    print(f"    Downregulated: {len(down)}")

    if len(up) > 0:
        print(f"\n  Top upregulated:")
        print(up.head(15)[["log2FoldChange","pvalue","padj"]].round(4).to_string())
    if len(down) > 0:
        print(f"\n  Top downregulated:")
        print(down.head(10)[["log2FoldChange","pvalue","padj"]].round(4).to_string())

    # Concordance with Wilcoxon
    wilcox = pd.read_csv("results/03_signatures/DEGs_perineural_vs_distal.csv", index_col=0)
    wilcox_up = wilcox[wilcox["logfoldchanges"]>0].nlargest(43,"logfoldchanges").index.tolist()
    res_up = results.loc[[g for g in wilcox_up if g in results.index]]
    concordant = res_up[res_up["log2FoldChange"]>0]
    print(f"\n  Wilcoxon concordance: {len(concordant)}/{len(res_up)} genes")

    results.to_csv(f"{OUTPUT_DIR}/pseudobulk_{label.replace(' ','_')}_results.csv")
    return results, up, down

# ── Run analyses ──────────────────────────────────────────────────────────────
res_A = run_pseudobulk_subset(
    HPV_NEG, ["Tumour","Tumour_Prolif"],
    "HPV-neg tumour cells", lfc_thresh=0.5
)
res_B = run_pseudobulk_subset(
    HPV_NEG, ["Tumour","Tumour_Prolif","Macrophage","T_cell","NK",
               "B_cell","Fibroblast","Endothelial","Mast"],
    "HPV-neg all cell types", lfc_thresh=0.5
)
res_C = run_pseudobulk_subset(
    HPV_POS, ["Tumour","Tumour_Prolif"],
    "HPV-pos tumour cells (negative control)", lfc_thresh=0.5
)

# ── Summary figure ────────────────────────────────────────────────────────────
print("\nGenerating summary comparison figure...")
fig, axes = plt.subplots(1, 3, figsize=(18, 7))
titles = [
    "HPV\u2212 Tumour Cells\n(n=5 patients)",
    "HPV\u2212 All Cell Types\n(n=5 patients)",
    "HPV+ Tumour Cells\n(negative control, n=5)",
]
colors = ["#2ecc71","#27ae60","#e74c3c"]

for ax, res_tuple, title, col in zip(axes,
    [res_A, res_B, res_C], titles, colors):
    if res_tuple is None:
        ax.text(0.5,0.5,"Insufficient data",transform=ax.transAxes,ha="center")
        ax.set_title(title, fontsize=11, fontweight="bold")
        continue

    results, up, down = res_tuple
    df = results.copy()
    df["neglog10p"] = -np.log10(df["padj"].clip(lower=1e-300))
    df["sig"] = (df["padj"]<0.05) & (df["log2FoldChange"].abs()>0.5)

    ns = df[~df["sig"]]
    u  = df[df["sig"] & (df["log2FoldChange"]>0)]
    d  = df[df["sig"] & (df["log2FoldChange"]<0)]

    ax.scatter(ns["log2FoldChange"], ns["neglog10p"],
               s=8, alpha=0.35, color="#94A3B8", rasterized=True)
    ax.scatter(u["log2FoldChange"],  u["neglog10p"],
               s=40, alpha=0.9, color=col, edgecolors="black",
               linewidth=0.5, label=f"Up (n={len(u)})")
    ax.scatter(d["log2FoldChange"],  d["neglog10p"],
               s=40, alpha=0.9, color="#2980b9", edgecolors="black",
               linewidth=0.5, label=f"Down (n={len(d)})")

    # Label top genes with adjustText to prevent overlap
    from adjustText import adjust_text
    top = pd.concat([df.nlargest(6,"log2FoldChange"),
                     df.nsmallest(4,"log2FoldChange")]).drop_duplicates()
    texts = []
    for gene, row in top.iterrows():
        t = ax.text(row["log2FoldChange"], row["neglog10p"],
                    str(gene), fontsize=8, fontweight="bold",
                    color="#1a1a1a")
        texts.append(t)
    adjust_text(texts, ax=ax,
                arrowprops=dict(arrowstyle="-", color="#888888",
                                lw=0.8, shrinkA=2, shrinkB=2),
                expand=(1.4, 1.6), force_text=(0.4, 0.6))

    ax.axvline(-0.5,color="black",linewidth=0.8,linestyle="--",alpha=0.4)
    ax.axvline( 0.5,color="black",linewidth=0.8,linestyle="--",alpha=0.4)
    ax.axhline(-np.log10(0.05),color="black",linewidth=0.8,linestyle="--",alpha=0.4)
    ax.set_xlabel("Log2 Fold Change (Perineural vs Distal)",fontsize=10)
    ax.set_ylabel("-log10(FDR)",fontsize=10)
    ax.set_title(title,fontsize=11,fontweight="bold")
    ax.legend(fontsize=9,framealpha=0.8)

plt.suptitle("Pseudobulk DEGs by HPV Status\nPerineural vs Distal Zone",
             fontsize=13,fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/fig_pseudobulk_HPVstratified.png",
            dpi=150,bbox_inches="tight")
plt.savefig("results/figures/fig_pseudobulk_HPVstratified.png",
            dpi=150,bbox_inches="tight")
plt.close()
print(f"  Saved → {OUTPUT_DIR}/fig_pseudobulk_HPVstratified.png")

print(f"\n{'='*55}")
print("Module 12d complete")
print(f"{'='*55}")
