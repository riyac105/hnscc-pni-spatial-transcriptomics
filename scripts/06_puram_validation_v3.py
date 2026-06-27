"""
10_puram_validation_v3.py
Full Puram et al. 2017 Cell correlation validation
380 shared genes available — sufficient for robust correlation

Run from ~/Desktop/PNI_project/
"""
import os, gzip
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import pearsonr, spearmanr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 0
OUTPUT_DIR = "results/10_puram_validation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ── Load Xenium data ───────────────────────────────────────────────────────────
print("Loading Xenium data...")
adata = sc.read_h5ad("results/08_cell_types/xenium_annotated.h5ad")
print(f"  {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
print(f"  Cell types: {list(adata.obs['cell_type'].unique())}")

# ── Load Puram reference ───────────────────────────────────────────────────────
print("\nLoading Puram et al. 2017 reference...")
raw = pd.read_csv("data/GSE103322_HNSCC_data.txt.gz",
                  sep="\t", index_col=0, compression="gzip",
                  header=None, low_memory=False)
print(f"  Raw shape: {raw.shape}")

# Find where numeric data starts
data_start = 0
for i in range(len(raw.index)):
    try:
        float(raw.iloc[i, 0])
        data_start = i
        break
    except (ValueError, TypeError):
        pass

meta = raw.iloc[:data_start].copy()
expr = raw.iloc[data_start:].copy()

# Strip quotes from gene names
expr.index = expr.index.str.strip("'\"")
expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)
print(f"  Expression: {expr.shape[0]} genes × {expr.shape[1]} cells")

# ── Extract Puram cell type labels ────────────────────────────────────────────
print("\nExtracting Puram cell type labels...")

# data_start=1 means metadata rows got included in expr — re-extract from raw
raw_clean = raw.copy()
raw_clean.index = raw_clean.index.astype(str).str.strip("'\"").str.strip()

META_ROWS = ["processed by Maxima enzyme", "Lymph node",
             "classified  as cancer cell", "classified as non-cancer cells",
             "non-cancer cell type"]

ct_labels = None
cancer_flag = None
for row_name in META_ROWS:
    if row_name in raw_clean.index:
        row_data = raw_clean.loc[row_name]
        if isinstance(row_data, pd.DataFrame):
            row_data = row_data.iloc[0]  # take first if duplicate
        if "non-cancer cell type" in row_name:
            ct_labels = row_data.astype(str)
            print(f"  Found: '{row_name}' — {ct_labels.value_counts().head(6).to_dict()}")
        if "classified" in row_name and "non" not in row_name:
            cancer_flag = row_data
            print(f"  Found: '{row_name}'")

# Also re-extract expr without metadata rows
expr = raw_clean.loc[[i for i in raw_clean.index
                      if i not in META_ROWS and i not in ["nan","NaN",""]]]
expr = expr.apply(pd.to_numeric, errors="coerce").fillna(0)
print(f"  Clean expression: {expr.shape[0]} genes × {expr.shape[1]} cells")

final_ct = ct_labels.copy()
if cancer_flag is not None:
    cancer_mask = cancer_flag.astype(str).str.strip() == "1"
    final_ct[cancer_mask] = "Malignant"

# Map to your cell type names
puram_map = {
    "Malignant"       : "Tumour",
    "T cell"          : "T_cell",
    "T cells"         : "T_cell",
    "B cell"          : "B_cell",
    "B cells"         : "B_cell",
    "Macrophage"      : "Macrophage",
    "Macrophages"     : "Macrophage",
    "Fibroblast"      : "Fibroblast",
    "Fibroblasts"     : "Fibroblast",
    "Endothelial"     : "Endothelial",
    "Endothelial cell": "Endothelial",
    "Mast cell"       : "Mast",
    "Mast"            : "Mast",
    "Dendritic"       : "DC_LAMP3",
    "DC"              : "DC_LAMP3",
    "NK cell"         : "DC_LAMP3",  # Puram NK may overlap with DC_LAMP3
    "NK"              : "DC_LAMP3",
    "Myofibroblast"   : "Fibroblast",
}
final_ct_mapped = final_ct.map(
    lambda x: puram_map.get(x.strip(), None)
)

print(f"  Puram cell type distribution:")
print(final_ct_mapped.value_counts(dropna=False).to_string())

# ── Find shared genes ──────────────────────────────────────────────────────────
shared = [g for g in adata.var_names if g in expr.index]
print(f"\n  Shared genes: {len(shared)}")

# ── Compute mean profiles ──────────────────────────────────────────────────────
print("\nComputing mean expression profiles...")

# Log-normalise Puram if needed
if expr.values.max() > 100:
    expr_log = np.log2(expr + 1)
else:
    expr_log = expr.copy()

# Puram profiles per cell type (shared genes only)
puram_profiles = {}
for ct in final_ct_mapped.dropna().unique():
    mask = final_ct_mapped == ct
    if mask.sum() >= 5:
        puram_profiles[ct] = expr_log.loc[shared, mask].mean(axis=1)
        print(f"  Puram {ct}: {mask.sum()} cells")

puram_df = pd.DataFrame(puram_profiles)  # genes × cell_types

# Your profiles per cell type (shared genes only)
your_profiles = {}
for ct in adata.obs["cell_type"].unique():
    mask = adata.obs["cell_type"] == ct
    X = adata[mask, shared].X
    if hasattr(X, "toarray"): X = X.toarray()
    your_profiles[ct] = np.array(X).mean(axis=0)
    print(f"  Your {ct}: {mask.sum()} cells")

your_df = pd.DataFrame(your_profiles, index=shared)  # genes × cell_types

# ── Compute Pearson correlation matrix ────────────────────────────────────────
print("\nComputing Pearson correlations...")
your_cts  = list(your_df.columns)
puram_cts = list(puram_df.columns)

corr_matrix = pd.DataFrame(index=your_cts, columns=puram_cts, dtype=float)
for yct in your_cts:
    for pct in puram_cts:
        r, p = pearsonr(your_df[yct].values,
                        puram_df[pct].values)
        corr_matrix.loc[yct, pct] = round(r, 3)

print("\nCorrelation matrix (your annotation × Puram reference):")
print(corr_matrix.to_string())
corr_matrix.to_csv(f"{OUTPUT_DIR}/puram_correlation_matrix_v3.csv")

# ── Check diagonal dominance ──────────────────────────────────────────────────
print("\nDiagonal dominance check:")
expected_pairs = {
    "Tumour"     : "Tumour",
    "Tumour_Prolif": "Tumour",
    "T_cell"     : "T_cell",
    "B_cell"     : "B_cell",
    "Macrophage" : "Macrophage",
    "DC_LAMP3"   : "DC_LAMP3",
    "Endothelial": "Endothelial",
    "Fibroblast" : "Fibroblast",
    "Mast"       : "Mast",
}
concordant = 0
for your_ct, expected_puram in expected_pairs.items():
    if your_ct not in corr_matrix.index:
        continue
    if expected_puram not in corr_matrix.columns:
        # Find best match
        row = corr_matrix.loc[your_ct].dropna()
        best_match = row.idxmax()
        best_r = row.max()
        print(f"  {your_ct:20s} expected={expected_puram:15s} "
              f"best_available={best_match} (r={best_r:.3f})")
        continue
    row = corr_matrix.loc[your_ct].dropna()
    best_match = row.idxmax()
    best_r = row.max()
    expected_r = corr_matrix.loc[your_ct, expected_puram]
    is_correct = best_match == expected_puram
    status = "✓" if is_correct else "✗"
    if is_correct: concordant += 1
    print(f"  {status} {your_ct:20s} → best={best_match:15s} "
          f"(r={best_r:.3f})  expected_r={expected_r:.3f}")

# ── Figure: Correlation heatmap ───────────────────────────────────────────────
print("\nGenerating correlation heatmap...")

# Sort for visual clarity
row_order = [ct for ct in ["Tumour","Tumour_Prolif","T_cell","B_cell",
             "Macrophage","DC_LAMP3","Endothelial","Fibroblast","Mast"]
             if ct in corr_matrix.index]
col_order = [ct for ct in ["Tumour","T_cell","B_cell","Macrophage",
             "DC_LAMP3","Endothelial","Fibroblast","Mast"]
             if ct in corr_matrix.columns]

corr_plot = corr_matrix.loc[row_order, col_order].astype(float)

fig, ax = plt.subplots(figsize=(max(10, len(col_order)*1.4),
                                max(7,  len(row_order)*0.9)))
fig.patch.set_facecolor("white")

sns.heatmap(corr_plot, annot=True, fmt=".3f", cmap="RdBu_r",
            center=0, vmin=-0.5, vmax=1.0,
            linewidths=0.8, linecolor="#EEEEEE",
            ax=ax, square=True,
            cbar_kws={"label": "Pearson r\n(380 shared genes)",
                      "shrink": 0.8})

# Highlight diagonal / expected matches
for i, your_ct in enumerate(row_order):
    for j, puram_ct in enumerate(col_order):
        if expected_pairs.get(your_ct) == puram_ct:
            ax.add_patch(plt.Rectangle((j, i), 1, 1,
                fill=False, edgecolor="lime",
                linewidth=3, zorder=5))

ax.set_title(
    "Cell Type Validation — Puram et al. 2017 (Cell) HNSCC Reference\n"
    f"Pearson correlation across {len(shared)} shared genes  |  "
    "Green boxes = expected concordant pairs\n"
    "Rows = your Xenium annotations  |  Columns = Puram scRNA-seq reference",
    fontsize=11, fontweight="bold", pad=12
)
ax.set_xlabel("Puram et al. 2017 cell type", fontsize=11)
ax.set_ylabel("Your annotation", fontsize=11)
ax.tick_params(axis='x', rotation=45, labelsize=10)
ax.tick_params(axis='y', rotation=0, labelsize=10)

plt.tight_layout()
out = f"{OUTPUT_DIR}/fig_puram_correlation_v3.png"
plt.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
plt.close()
print(f"  Saved → {out}")

print(f"\n{'='*60}")
print("PURAM VALIDATION SUMMARY")
print(f"{'='*60}")
print(f"Shared genes: {len(shared)}/399 ({len(shared)/399*100:.0f}% of Xenium panel)")
print(f"Puram cell types available: {puram_cts}")
print(f"Concordant annotations: {concordant}/{len([p for p in expected_pairs.values() if p in corr_matrix.columns])}")
