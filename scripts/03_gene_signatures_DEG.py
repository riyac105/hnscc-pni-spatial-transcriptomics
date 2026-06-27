"""
Module 3: PNI Gene Signatures & Differential Expression
========================================================
Project: Perineural Invasion (PNI) in HNSCC — Spatial Transcriptomics
Platform: 10x Xenium
Description:
    - DEG analysis: Perineural tumor cells vs. Distal tumor cells
    - Pathway enrichment (GSEA / ORA)
    - PNI gene signature scoring
    - Spatially variable gene analysis using Squidpy
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import os

# Optional — install if available: pip install gseapy decoupler
try:
    import gseapy as gp
    GSEAPY = True
except ImportError:
    GSEAPY = False
    print("gseapy not installed — pathway enrichment will be skipped. "
          "Install with: pip install gseapy")

try:
    import decoupler as dc
    DECOUPLER = True
except ImportError:
    DECOUPLER = False

# ── Paths ──────────────────────────────────────────────────────────────────────
INPUT_H5AD = "results/02_nerve_id/xenium_nerve_zones.h5ad"
OUTPUT_DIR = "results/03_signatures/"
os.makedirs(OUTPUT_DIR, exist_ok=True)
sc.settings.figdir = OUTPUT_DIR

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
adata = sc.read_h5ad(INPUT_H5AD)
print(f"Loaded: {adata.shape[0]:,} cells")

# Restore raw counts layer for DEG testing
adata.X = adata.layers["counts"].copy()

# ==============================================================================
# 2. SUBSET TO TUMOR CELLS
# ==============================================================================
tumor_adata = adata[adata.obs["is_tumor"] == 1].copy()
print(f"Tumor cells: {tumor_adata.n_obs:,}")

# Keep only Perineural and Distal zones for the main comparison
tumor_adata = tumor_adata[
    tumor_adata.obs["spatial_zone"].isin(["Perineural", "Distal"])
].copy()

n_peri = (tumor_adata.obs["spatial_zone"] == "Perineural").sum()
n_dist = (tumor_adata.obs["spatial_zone"] == "Distal").sum()
print(f"  Perineural tumor cells : {n_peri:,}")
print(f"  Distal tumor cells     : {n_dist:,}")

if n_peri < 10 or n_dist < 10:
    raise ValueError("Too few cells in one zone — check thresholds in Module 2.")

# ==============================================================================
# 3. DIFFERENTIAL GENE EXPRESSION (Wilcoxon rank-sum)
# ==============================================================================
print("\nRunning DEG analysis (Perineural vs Distal tumor cells)...")

sc.tl.rank_genes_groups(
    tumor_adata,
    groupby="spatial_zone",
    groups=["Perineural"],
    reference="Distal",
    method="wilcoxon",
    key_added="deg_perineural_vs_distal",
    use_raw=False
)

# Extract results
deg_df = sc.get.rank_genes_groups_df(
    tumor_adata,
    group="Perineural",
    key="deg_perineural_vs_distal",
    pval_cutoff=0.05,
)
deg_df = deg_df.sort_values("scores", ascending=False)

print(f"  Significant DEGs (FDR<0.05): {len(deg_df):,}")
print(f"  Top 10 upregulated in Perineural zone:")
print(deg_df.head(10)[["names", "scores", "logfoldchanges", "pvals_adj"]].to_string())

deg_df.to_csv(f"{OUTPUT_DIR}/DEGs_perineural_vs_distal.csv", index=False)

# ── Volcano plot ───────────────────────────────────────────────────────────────
all_deg = sc.get.rank_genes_groups_df(
    tumor_adata,
    group="Perineural",
    key="deg_perineural_vs_distal"
)
all_deg["-log10_padj"] = -np.log10(all_deg["pvals_adj"].clip(1e-300))

fig, ax = plt.subplots(figsize=(8, 6))
colors = np.where(
    (all_deg["pvals_adj"] < 0.05) & (all_deg["logfoldchanges"] > 0.5),
    "firebrick",
    np.where(
        (all_deg["pvals_adj"] < 0.05) & (all_deg["logfoldchanges"] < -0.5),
        "steelblue", "lightgray"
    )
)
ax.scatter(all_deg["logfoldchanges"], all_deg["-log10_padj"],
           c=colors, s=8, alpha=0.7, linewidths=0)
ax.axhline(-np.log10(0.05), color="black", ls="--", lw=0.8, alpha=0.5)
ax.axvline(0.5,  color="black", ls="--", lw=0.8, alpha=0.5)
ax.axvline(-0.5, color="black", ls="--", lw=0.8, alpha=0.5)

# Label top 15 genes
top_genes = deg_df.head(15)
for _, row in top_genes.iterrows():
    ax.annotate(row["names"],
                xy=(row["logfoldchanges"],
                    -np.log10(max(row["pvals_adj"], 1e-300))),
                fontsize=6, alpha=0.8)

ax.set_xlabel("Log2 fold change (Perineural / Distal)")
ax.set_ylabel("-log10 adjusted p-value")
ax.set_title("Perineural vs. Distal tumor cells — Volcano plot")
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/volcano_perineural_vs_distal.pdf", dpi=150)
plt.close()

# ── Dotplot top 30 DEGs ────────────────────────────────────────────────────────
top30 = deg_df["names"].head(30).tolist()
sc.pl.dotplot(
    tumor_adata,
    var_names=top30,
    groupby="spatial_zone",
    save="_top30_DEGs.pdf",
    show=False
)

# ==============================================================================
# 4. KNOWN PNI GENE SIGNATURE SCORING
# ==============================================================================
# Curated PNI signatures from literature:
#   - Amit et al. Nature Communications 2020 (HNSCC PNI)
#   - Liebig et al. Cancer 2009 (PNI pan-cancer)
#   - Boilly et al. Cancer Cell 2017 (neurotrophic signalling)

PNI_SIGNATURE = {
    "PNI_neurotrophin_signalling": [
        "NGF", "BDNF", "NTF3", "NTF4",        # Neurotrophins (ligands)
        "NTRK1", "NTRK2", "NTRK3", "NGFR",    # Receptors
        "ARTN", "GFRA1", "RET",                 # GFL / GDNF family
    ],
    "PNI_axon_guidance": [
        "SEMA3A", "SEMA3F", "SEMA4D",
        "NRP1", "NRP2",
        "PLXNA1", "PLXNB1",
        "ROBO1", "SLIT2",
        "EFNB1", "EPHB2",
    ],
    "PNI_invasion_EMT": [
        "MMP2", "MMP9", "MMP14",
        "TWIST1", "SNAI1", "SNAI2", "ZEB1",
        "VIM", "FN1", "CDH2",
    ],
    "PNI_immune_evasion": [
        "CD274",   # PD-L1
        "PDCD1LG2","CTLA4",
        "TGFB1", "TGFB2", "IL10",
        "IDO1",
    ],
}

panel_genes = set(adata.var_names)
for sig_name, gene_list in PNI_SIGNATURE.items():
    genes_in_panel = [g for g in gene_list if g in panel_genes]
    if genes_in_panel:
        sc.tl.score_genes(adata, gene_list=genes_in_panel,
                          score_name=sig_name)
        print(f"  Scored {sig_name}: "
              f"{len(genes_in_panel)}/{len(gene_list)} genes in panel")
    else:
        print(f"  SKIP {sig_name}: no genes in panel")

# Plot signature scores across spatial zones (all cells)
sig_scores = [s for s in PNI_SIGNATURE if s in adata.obs.columns]
if sig_scores:
    fig, axes = plt.subplots(1, len(sig_scores),
                              figsize=(5*len(sig_scores), 5))
    if len(sig_scores) == 1:
        axes = [axes]
    zone_order = ["Nerve", "Perineural", "Peritumoral", "Distal"]
    zone_colors = ["#e74c3c", "#f39c12", "#3498db", "#bdc3c7"]

    for ax, sig in zip(axes, sig_scores):
        data_by_zone = [
            adata.obs.loc[adata.obs["spatial_zone"] == z, sig].dropna().values
            for z in zone_order
        ]
        ax.boxplot(data_by_zone, labels=zone_order,
                   patch_artist=True,
                   boxprops=dict(facecolor="white"),
                   medianprops=dict(color="black", lw=2))
        for patch, color in zip(ax.patches, zone_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.7)
        ax.set_title(sig.replace("PNI_", "").replace("_", " "), fontsize=9)
        ax.set_ylabel("Module score")
        ax.tick_params(axis='x', rotation=30)

    plt.suptitle("PNI signature scores across spatial zones", fontsize=12)
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/PNI_signature_scores_by_zone.pdf", dpi=150)
    plt.close()

# Spatial maps of PNI scores
for sig in sig_scores:
    sq.pl.spatial_scatter(
        adata, color=sig, shape=None, size=0.5,
        cmap="RdYlBu_r",
        title=sig.replace("_", " "),
        save=f"spatial_{sig}.pdf"
    )

# ==============================================================================
# 5. SPATIALLY VARIABLE GENES (Squidpy Moran's I)
# ==============================================================================
print("\nComputing spatially variable genes (Moran's I)...")

# Build spatial graph if not already present
if "spatial_neighbors" not in adata.uns:
    sq.gr.spatial_neighbors(adata, coord_type="generic",
                             n_neighs=10, key_added="spatial")

# n_perms=0 avoids multiprocessing which breaks on Mac with nohup
# Moran I skipped — multiprocessing incompatibility on Mac
print("Skipping Moran I spatial autocorr")

# ==============================================================================
# 6. PATHWAY ENRICHMENT (ORA with gseapy)
# ==============================================================================
if GSEAPY:
    print("\nRunning pathway enrichment (ORA)...")

    up_genes   = deg_df[deg_df["logfoldchanges"] > 0.5]["names"].tolist()
    down_genes = deg_df[deg_df["logfoldchanges"] < -0.5]["names"].tolist()

    for direction, gene_list in [("upregulated", up_genes),
                                   ("downregulated", down_genes)]:
        if not gene_list:
            continue
        try:
            enr = gp.enrichr(
                gene_list=gene_list,
                gene_sets=["MSigDB_Hallmark_2020",
                           "KEGG_2021_Human",
                           "Reactome_2022"],
                organism="human",
                outdir=f"{OUTPUT_DIR}/enrichr_{direction}",
                cutoff=0.05
            )
            # Plot top results
            gp.barplot(enr.results,
                       title=f"ORA — {direction} in Perineural",
                       ofname=f"{OUTPUT_DIR}/ORA_{direction}_barplot.pdf")
            print(f"  ORA ({direction}): "
                  f"{len(enr.results):,} enriched terms")
        except Exception as e:
            print(f"  ORA failed ({direction}): {e}")
else:
    print("\ngseapy not available — skipping ORA.")
    print("  Re-run after: pip install gseapy")

# ==============================================================================
# 7. SAVE
# ==============================================================================
out_path = "results/03_signatures/xenium_signatures.h5ad"
adata.write_h5ad(out_path)
print(f"\nSaved → {out_path}")
print("Module 3 complete ✓")
