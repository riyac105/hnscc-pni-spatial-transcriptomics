"""
Module 2: Computational Nerve Identification & PNI Zone Definition
==================================================================
Project: Perineural Invasion (PNI) in HNSCC — Spatial Transcriptomics
Platform: 10x Xenium
Description:
    - Score each cell for Schwann/nerve identity using canonical marker genes
    - Classify cells into nerve, tumor, stromal, immune compartments
    - Define spatial PNI zones by proximity analysis (nerve-adjacent vs. distal)
    - Visualise nerve maps and invasion zones
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.spatial import cKDTree
import os

# ── Paths ──────────────────────────────────────────────────────────────────────
INPUT_H5AD  = "results/01_qc/GSE300147_all_samples_qc.h5ad"
OUTPUT_DIR  = "results/02_nerve_id/"
os.makedirs(OUTPUT_DIR, exist_ok=True)
sc.settings.figdir = OUTPUT_DIR

# ==============================================================================
# 1. LOAD DATA
# ==============================================================================
adata = sc.read_h5ad(INPUT_H5AD)
print(f"Loaded: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")

# ==============================================================================
# 2. CELL-TYPE GENE SIGNATURES
# ==============================================================================
# Marker gene sets — intersect with detected panel genes
# References:
#   Nerve/Schwann : Jessen & Mirsky (2005); PNI lit (Amit et al. 2020)
#   Tumor (HNSCC) : Puram et al. Cell 2017
#   Fibroblast    : Elyada et al. 2019; Cancer-associated fibroblast markers
#   Immune        : Pan-cancer immune markers


# ── UPDATED for GSE300147 panel (477 genes) ────────────────────────────────────
# Nerve markers confirmed present in panel audit:
#   MBP, MPZ, PRX, PLP1         — peripheral myelin proteins
#   S100B, SOX10, NGFR, GFAP    — Schwann/glial markers
#   CDH19                        — Schwann cell adhesion molecule
# NOTE: TUBB3, NEFM, NTRK1, NGF, RET etc. are ABSENT from this panel
# NOTE: Tumour markers sparse (KRT5/KRT14 present, EGFR/EPCAM absent)
#       Use KRT5+KRT14+TP63 as tumour signature

# Marker genes confirmed present in actual GSE300147 panel (399 genes)
MARKER_GENES = {
    # Nerve / Schwann — confirmed in panel
    "Schwann_nerve"          : ["PMP22", "EDNRB", "PTN", "LGI4"],

    # Tumour HNSCC — confirmed in panel
    "Tumor_HNSCC"            : ["EGFR", "EPCAM", "SOX2", "SERPINB3",
                                  "SERPINB2", "CLCA2", "GPRC5A", "KRT7",
                                  "KRT20"],

    # Stroma / CAF — confirmed in panel
    "Fibroblast_CAF"         : ["ACTA2", "PDGFRA", "PDGFRB", "PDPN",
                                  "THY1", "FBN1", "SFRP2", "SFRP4"],

    # Endothelial — confirmed in panel
    "Endothelial"            : ["PECAM1", "CD34", "VWF", "EGFL7",
                                  "CLEC14A", "RAMP2"],

    # Immune — well covered
    "Macrophage"             : ["CD68", "CD163", "MRC1", "CD14",
                                  "ADGRE1", "VSIG4", "TREM2"],
    "T_cell"                 : ["CD3D", "CD3E", "CD8A", "CD4",
                                  "FOXP3", "GZMB", "GZMA", "GZMK"],
    "B_cell"                 : ["CD19", "MS4A1", "CD79A", "BANK1"],
    "NK_cell"                : ["NKG7", "GNLY", "KLRD1", "KLRB1"],
    "DC"                     : ["CD1C", "CLEC10A", "LAMP3", "CD83",
                                  "LILRA4"],
    "Proliferating"          : ["MKI67", "TOP2A", "CDK1", "PCNA",
                                  "UBE2C", "CENPF"],
}

# Filter to genes present in the panel
panel_genes = set(adata.var_names)
filtered_markers = {}
for ct, genes in MARKER_GENES.items():
    detected = [g for g in genes if g in panel_genes]
    if detected:
        filtered_markers[ct] = detected
        print(f"  {ct}: {len(detected)}/{len(genes)} markers detected → {detected}")
    else:
        print(f"  WARNING — {ct}: NO markers detected in panel")

# ==============================================================================
# 3. SCORE EACH CELL TYPE WITH sc.tl.score_genes
# ==============================================================================
print("\nScoring cell types...")
for ct, genes in filtered_markers.items():
    sc.tl.score_genes(adata, gene_list=genes,
                      score_name=f"score_{ct}", use_raw=False)

# ── Aggregate nerve score ──────────────────────────────────────────────────────
# Use both Schwann categories confirmed in panel audit
nerve_keys = [k for k in filtered_markers if "Schwann" in k or "nerve" in k.lower()]
if nerve_keys:
    adata.obs["score_Nerve_composite"] = adata.obs[
        [f"score_{k}" for k in nerve_keys]
    ].mean(axis=1)

# ==============================================================================
# 4. CELL-TYPE ASSIGNMENT (SCORE-BASED)
# ==============================================================================
score_cols = [c for c in adata.obs.columns if c.startswith("score_")]
score_df   = adata.obs[score_cols].copy()
score_df.columns = [c.replace("score_", "") for c in score_df.columns]

# Assign cell type as the compartment with the highest score
adata.obs["cell_type_score"] = score_df.idxmax(axis=1)

# ── Hard threshold: flag high-confidence nerve cells ─────────────────────────
NERVE_THRESHOLD = 0.5   # targeting ~5-15% nerve cells per patient   # lowered for GSE300147 — fewer nerve markers in panel
                        # if too few nerve cells identified, try 0.15
                        # if too many false positives, try 0.25
if "score_Nerve_composite" in adata.obs.columns:
    adata.obs["is_nerve"] = (
        adata.obs["score_Nerve_composite"] > NERVE_THRESHOLD
    ).astype(int)
else:
    # Fallback: use best individual nerve score
    best_nerve = score_df[
        [c for c in score_df.columns if "Schwann" in c or "Nerve" in c]
    ].max(axis=1)
    adata.obs["is_nerve"] = (best_nerve > NERVE_THRESHOLD).astype(int)

n_nerve = adata.obs["is_nerve"].sum()
print(f"\n  Nerve cells identified: {n_nerve:,} "
      f"({100*n_nerve/adata.n_obs:.1f}% of all cells)")

# ── Flag tumor cells ──────────────────────────────────────────────────────────
TUMOR_THRESHOLD = 0.2
adata.obs["is_tumor"] = (
    adata.obs["score_Tumor_HNSCC"] > TUMOR_THRESHOLD
).astype(int)

# ==============================================================================
# 5. SPATIAL PNI ZONE DEFINITION VIA KD-TREE PROXIMITY
# ==============================================================================
print("\nDefining PNI spatial zones...")

coords      = adata.obsm["spatial"]                  # (N, 2) array
nerve_mask  = adata.obs["is_nerve"].values == 1
nerve_coords = coords[nerve_mask]

if nerve_coords.shape[0] == 0:
    raise ValueError("No nerve cells identified — lower NERVE_THRESHOLD or "
                     "verify marker genes are in your panel.")

# Build KD-tree on nerve cell centroids
tree = cKDTree(nerve_coords)

# Query distance from every cell to its nearest nerve cell
dist_to_nerve, _ = tree.query(coords, k=1)
adata.obs["dist_to_nerve_um"] = dist_to_nerve   # units = microns (Xenium default)

# ── Define zones (adjust radii to match your tissue biology) ──────────────────
# Xenium coordinates are in microns
PERINEURAL_RADIUS_UM  = 7     # ~25th percentile of cell-nerve distances    # ≤50 µm  → perineural zone
PERITUMORAL_RADIUS_UM = 20    # ~75th percentile of cell-nerve distances   # 50–200 µm → peritumoral zone

def assign_zone(row):
    if row["is_nerve"] == 1:
        return "Nerve"
    d = row["dist_to_nerve_um"]
    if d <= PERINEURAL_RADIUS_UM:
        return "Perineural"
    elif d <= PERITUMORAL_RADIUS_UM:
        return "Peritumoral"
    else:
        return "Distal"

adata.obs["spatial_zone"] = adata.obs.apply(assign_zone, axis=1)

zone_counts = adata.obs["spatial_zone"].value_counts()
print("\n  Zone breakdown:")
print(zone_counts.to_string())

# ── PNI label: tumor cell in perineural zone ──────────────────────────────────
adata.obs["PNI_positive"] = (
    (adata.obs["is_tumor"] == 1) &
    (adata.obs["spatial_zone"] == "Perineural")
).astype(int)

n_pni = adata.obs["PNI_positive"].sum()
print(f"\n  PNI-positive tumor cells: {n_pni:,}")

# ==============================================================================
# 6. VISUALISATION
# ==============================================================================
print("\nGenerating spatial maps...")

# ── Nerve score map ────────────────────────────────────────────────────────────
sq.pl.spatial_scatter(
    adata,
    color="score_Nerve_composite" if "score_Nerve_composite" in adata.obs
          else f"score_{nerve_keys[0]}",
    shape=None, size=0.5,
    cmap="magma",
    title="Nerve composite score",
    save="nerve_score_spatial.pdf"
)

# ── Spatial zone map ──────────────────────────────────────────────────────────
zone_palette = {
    "Nerve"       : "#e74c3c",
    "Perineural"  : "#f39c12",
    "Peritumoral" : "#3498db",
    "Distal"      : "#bdc3c7",
}
# Assign palette as categorical colours directly on adata
import pandas as pd
adata.obs["spatial_zone"] = pd.Categorical(
    adata.obs["spatial_zone"],
    categories=["Nerve", "Perineural", "Peritumoral", "Distal"]
)
adata.uns["spatial_zone_colors"] = ["#e74c3c", "#f39c12", "#3498db", "#bdc3c7"]
sq.pl.spatial_scatter(
    adata, color="spatial_zone",
    shape=None, size=0.5,
    title="Spatial PNI zones",
    save="spatial_zones.pdf"
)

# ── Distance to nerve histogram ────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(adata.obs.loc[adata.obs["is_tumor"]==1, "dist_to_nerve_um"],
        bins=100, color="steelblue", edgecolor="none", alpha=0.8,
        label="Tumor cells")
ax.axvline(PERINEURAL_RADIUS_UM,  color="red",    ls="--", label=f"Perineural ({PERINEURAL_RADIUS_UM} µm)")
ax.axvline(PERITUMORAL_RADIUS_UM, color="orange", ls="--", label=f"Peritumoral ({PERITUMORAL_RADIUS_UM} µm)")
ax.set_xlabel("Distance to nearest nerve cell (µm)")
ax.set_ylabel("# Tumor cells")
ax.set_title("Tumor cell distance to nerve")
ax.legend()
plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/dist_to_nerve_histogram.pdf", dpi=150)
plt.close()

# ── Cell type map ──────────────────────────────────────────────────────────────
sq.pl.spatial_scatter(
    adata, color="cell_type_score",
    shape=None, size=0.5,
    title="Cell type (score-based)",
    save="cell_type_spatial.pdf"
)

# ── UMAP coloured by zone ──────────────────────────────────────────────────────
sc.pl.umap(adata, color=["spatial_zone", "is_nerve", "PNI_positive"],
           ncols=3, save="_zones.pdf", show=False)

# ==============================================================================
# 7. SAVE
# ==============================================================================
out_path = "results/02_nerve_id/xenium_nerve_zones.h5ad"
adata.write_h5ad(out_path)
print(f"\nSaved → {out_path}")
print("Module 2 complete ✓")
