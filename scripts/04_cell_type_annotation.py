"""
08_cell_type_annotation.py
Annotates Leiden clusters with cell type labels using known marker genes.
Transfers labels to the zones object and generates publication figures.
Run from ~/Desktop/PNI_project/

Fixes applied:
  1. Proliferating cluster renamed to Tumour_Prolif (proliferating tumour cells)
  2. Schwann/Nerve included in spatial map legend and plots
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

sc.settings.verbosity = 1

OUTPUT_DIR = "results/08_cell_types"
os.makedirs(OUTPUT_DIR, exist_ok=True)

MARKERS = {
    'Schwann/Nerve' : ['PMP22', 'EDNRB', 'PTN', 'LGI4'],
    'Tumour'        : ['EGFR', 'EPCAM', 'SOX2', 'KRT7', 'CLCA2', 'SERPINB3'],
    'T_cell'        : ['CD3D', 'CD3E', 'CD8A', 'CD4', 'FOXP3', 'GZMB'],
    'B_cell'        : ['CD19', 'MS4A1', 'CD79A'],
    'NK'            : ['NKG7', 'GNLY', 'KLRD1'],
    'Macrophage'    : ['CD68', 'CD163', 'MRC1', 'CD14', 'TREM2'],
    'DC'            : ['CD1C', 'LAMP3', 'CLEC10A'],
    'Fibroblast'    : ['ACTA2', 'PDGFRA', 'PDGFRB', 'THY1'],
    'Endothelial'   : ['PECAM1', 'CD34', 'VWF', 'EGFL7'],
    'Tumour_Prolif' : ['MKI67', 'TOP2A', 'CDK1', 'UBE2C'],
    'Mast'          : ['KIT', 'CPA3', 'MS4A2'],
}

CT_COLORS = {
    'Tumour'        : '#e74c3c',
    'Tumour_Prolif' : '#ff5722',
    'Schwann/Nerve' : '#9b59b6',
    'T_cell'        : '#3498db',
    'B_cell'        : '#1abc9c',
    'NK'            : '#2ecc71',
    'Macrophage'    : '#e67e22',
    'DC'            : '#f39c12',
    'Fibroblast'    : '#795548',
    'Endothelial'   : '#e91e63',
    'Mast'          : '#607d8b',
    'Unknown'       : '#bdc3c7',
}

# ==============================================================================
print("Loading QC object...")
adata = sc.read_h5ad("results/01_qc/GSE300147_all_samples_qc.h5ad")
print(f"  {adata.shape[0]:,} cells x {adata.shape[1]:,} genes")

print("\nScoring clusters against marker genes...")
for ct, genes in MARKERS.items():
    present = [g for g in genes if g in adata.var_names]
    if present:
        sc.tl.score_genes(adata, gene_list=present, score_name=f"score_{ct}", use_raw=False)
        print(f"  {ct:15s}: {len(present)}/{len(genes)} genes")

score_cols = [f"score_{ct}" for ct in MARKERS if f"score_{ct}" in adata.obs.columns]
cluster_scores = (
    adata.obs.groupby("leiden_r05")[score_cols]
    .mean()
    .rename(columns={f"score_{ct}": ct for ct in MARKERS})
)
print("\nCluster x cell type score matrix:")
print(cluster_scores.round(3).to_string())

print("\nAuto-assigning cell types to clusters...")
cluster_to_celltype = {}
for cluster in cluster_scores.index:
    scores = cluster_scores.loc[cluster]
    best_ct = scores.idxmax()
    best_score = scores.max()
    cluster_to_celltype[cluster] = best_ct if best_score > 0.05 else "Unknown"
    print(f"  Cluster {cluster:>2s} -> {cluster_to_celltype[cluster]:15s} (score={best_score:.3f})")

adata.obs["cell_type"] = adata.obs["leiden_r05"].map(cluster_to_celltype).fillna("Unknown")
print("\nCell type counts (QC object):")
print(adata.obs["cell_type"].value_counts().to_string())

# Dotplot
print("\nGenerating marker gene dotplot...")
fig, ax = plt.subplots(figsize=(20, 8))
sc.pl.dotplot(
    adata,
    var_names={ct: [g for g in genes if g in adata.var_names] for ct, genes in MARKERS.items()},
    groupby="leiden_r05", ax=ax, show=False,
    title="Marker gene expression per Leiden cluster", standard_scale="var",
)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/dotplot_clusters_vs_markers.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved -> {OUTPUT_DIR}/dotplot_clusters_vs_markers.png")

# UMAP
print("\nGenerating cell type UMAP...")
umap = adata.obsm["X_umap"]
cell_types = adata.obs["cell_type"].astype(str)
unique_cts = [ct for ct in CT_COLORS if ct in cell_types.unique()]

fig, axes = plt.subplots(1, 2, figsize=(18, 7))
for ct in unique_cts:
    mask = cell_types == ct
    axes[0].scatter(umap[mask,0], umap[mask,1], s=0.4, alpha=0.5,
                    color=CT_COLORS[ct], rasterized=True)
axes[0].set_title("Cell Type Annotation", fontsize=14, fontweight="bold")
axes[0].set_xlabel("UMAP1"); axes[0].set_ylabel("UMAP2")
axes[0].set_xticks([]); axes[0].set_yticks([])
axes[0].legend(handles=[mpatches.Patch(facecolor=CT_COLORS[ct], label=ct) for ct in unique_cts],
               fontsize=9, bbox_to_anchor=(1.01,1), loc="upper left", framealpha=0.8)

hpv = adata.obs["hpv_status"].astype(str)
for h, c in {"HPV+":"#e74c3c","HPV-":"#2ecc71"}.items():
    mask = hpv == h
    axes[1].scatter(umap[mask,0], umap[mask,1], s=0.4, alpha=0.4, color=c, label=h, rasterized=True)
axes[1].set_title("HPV Status", fontsize=14, fontweight="bold")
axes[1].set_xlabel("UMAP1"); axes[1].set_yticks([]); axes[1].set_xticks([])
axes[1].legend(markerscale=6, fontsize=11, framealpha=0.8)
plt.suptitle("Cell Type Annotation - 560,492 cells, 10 HNSCC patients", fontsize=15, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/umap_cell_types.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved -> {OUTPUT_DIR}/umap_cell_types.png")

# Composition per patient
print("\nGenerating cell type composition chart...")
comp = adata.obs.groupby(["patient_id","cell_type"]).size().unstack(fill_value=0)
comp_pct = comp.div(comp.sum(axis=1), axis=0) * 100
comp_pct = comp_pct[comp_pct.mean().sort_values(ascending=False).index]
meta = adata.obs[["patient_id","hpv_status"]].drop_duplicates().set_index("patient_id")
comp_pct = comp_pct.join(meta).sort_values("hpv_status")
hpv_labels = comp_pct.pop("hpv_status")
colors = [CT_COLORS.get(ct,"#bdc3c7") for ct in comp_pct.columns]
fig, ax = plt.subplots(figsize=(14, 6))
comp_pct.plot(kind="bar", stacked=True, color=colors, ax=ax, edgecolor="white", linewidth=0.3)
ax.set_xlabel("Patient", fontsize=12); ax.set_ylabel("% of cells", fontsize=12)
ax.set_title("Cell Type Composition per Patient", fontsize=14, fontweight="bold")
ax.set_xticklabels([f"{p}\n({hpv_labels[p]})" for p in comp_pct.index], rotation=45, ha="right", fontsize=10)
ax.legend(title="Cell type", bbox_to_anchor=(1.01,1), loc="upper left", fontsize=9, framealpha=0.8)
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/cell_type_composition_per_patient.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved -> {OUTPUT_DIR}/cell_type_composition_per_patient.png")

# Transfer to zones + nerve override
print("\nTransferring cell type labels to zones object...")
adata_zones = sc.read_h5ad("results/02_nerve_id/xenium_nerve_zones.h5ad")
label_map = adata.obs["cell_type"].to_dict()
adata_zones.obs["cell_type"] = adata_zones.obs_names.map(label_map).fillna("Unknown")
if "is_nerve" in adata_zones.obs.columns:
    nerve_mask = adata_zones.obs["is_nerve"] == 1
    adata_zones.obs.loc[nerve_mask, "cell_type"] = "Schwann/Nerve"
    print(f"  Overrode {nerve_mask.sum():,} cells as Schwann/Nerve from spatial scoring")
print("\nCell type counts (zones object):")
print(adata_zones.obs["cell_type"].value_counts().to_string())

# unique_cts for zones — includes Schwann/Nerve
unique_cts_zones = [ct for ct in CT_COLORS if ct in adata_zones.obs["cell_type"].unique()]

# Spatial map with Schwann/Nerve in legend
print("\nGenerating spatial cell type map...")
coords   = adata_zones.obsm["spatial"]
ct_zones = adata_zones.obs["cell_type"].astype(str)
pts_z    = adata_zones.obs["patient_id"].astype(str)
ncols = 4
nrows = int(np.ceil(pts_z.nunique() / ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(ncols*4, nrows*4))
axes = axes.flatten()
for i, pt in enumerate(sorted(pts_z.unique())):
    ax = axes[i]
    mask_pt = pts_z == pt
    hpv = adata_zones.obs.loc[mask_pt, "hpv_status"].iloc[0]
    # Draw Schwann/Nerve last so it's on top
    ct_order = [ct for ct in unique_cts_zones if ct != "Schwann/Nerve"] + ["Schwann/Nerve"]
    for ct in ct_order:
        mask = mask_pt & (ct_zones == ct)
        if mask.sum() > 0:
            ax.scatter(coords[mask,0], coords[mask,1],
                       s=1.0 if ct == "Schwann/Nerve" else 0.5,
                       alpha=0.8 if ct == "Schwann/Nerve" else 0.5,
                       color=CT_COLORS[ct], rasterized=True)
    ax.set_title(f"{pt} ({hpv})", fontsize=11, fontweight="bold")
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    ax.invert_yaxis()
for j in range(i+1, len(axes)):
    axes[j].set_visible(False)
legend_elements = [mpatches.Patch(facecolor=CT_COLORS[ct], label=ct) for ct in unique_cts_zones]
fig.legend(handles=legend_elements, fontsize=9, loc="lower right",
           framealpha=0.9, ncol=2, bbox_to_anchor=(0.99, 0.01))
plt.suptitle("Spatial Cell Type Distribution - All Patients", fontsize=15, fontweight="bold")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/spatial_cell_types_all_patients.png", dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved -> {OUTPUT_DIR}/spatial_cell_types_all_patients.png")

# Zone composition
print("\nCell type composition within PNI zones...")
if "spatial_zone" in adata_zones.obs.columns:
    zone_ct = adata_zones.obs.groupby(["spatial_zone","cell_type"]).size().unstack(fill_value=0)
    zone_ct_pct = zone_ct.div(zone_ct.sum(axis=1), axis=0) * 100
    zone_ct_pct = zone_ct_pct.reindex(["Nerve","Perineural","Peritumoral","Distal"]).dropna()
    col_order = ["Schwann/Nerve"] + [c for c in zone_ct_pct.mean().sort_values(ascending=False).index
                                      if c != "Schwann/Nerve" and c in zone_ct_pct.columns]
    zone_ct_pct = zone_ct_pct[[c for c in col_order if c in zone_ct_pct.columns]]
    colors = [CT_COLORS.get(ct,"#bdc3c7") for ct in zone_ct_pct.columns]
    fig, ax = plt.subplots(figsize=(10, 6))
    zone_ct_pct.plot(kind="bar", stacked=True, color=colors, ax=ax, edgecolor="white", linewidth=0.3)
    ax.set_xlabel("Spatial Zone", fontsize=12); ax.set_ylabel("% of cells", fontsize=12)
    ax.set_title("Cell Type Composition per PNI Zone", fontsize=14, fontweight="bold")
    ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=11)
    ax.legend(title="Cell type", bbox_to_anchor=(1.01,1), loc="upper left", fontsize=9, framealpha=0.8)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/cell_type_composition_by_zone.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved -> {OUTPUT_DIR}/cell_type_composition_by_zone.png")

# Save
print("\nSaving updated objects...")
adata.write_h5ad(f"{OUTPUT_DIR}/xenium_annotated.h5ad")
adata_zones.obs[["cell_type"]].to_csv(f"{OUTPUT_DIR}/cell_type_labels.csv")

print(f"\n{'='*60}")
print("Module 8 (Cell Type Annotation) complete")
print(f"{'='*60}")
for f in sorted(os.listdir(OUTPUT_DIR)):
    if f.endswith('.png'):
        size = os.path.getsize(os.path.join(OUTPUT_DIR, f)) // 1024
        print(f"  {size:>5} KB -- {f}")
