"""
13_neighbourhood_analysis.py  v2 — Mac-compatible
Custom neighbourhood enrichment without squidpy multiprocessing.

Analyses:
  1. Cell type neighbourhood enrichment around nerves (manual permutation)
  2. NK cell distance to nearest nerve — HPV+ vs HPV-
  3. Tumour cell proximity to nerves — HPV+ vs HPV-
  4. Neighbourhood composition heatmap

Run from ~/Desktop/PNI_project/
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from scipy import stats
from scipy.spatial import cKDTree
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 0
OUTPUT_DIR  = "results/13_neighbourhood"
FIGURES_DIR = "results/figures"
os.makedirs(OUTPUT_DIR, exist_ok=True)
os.makedirs(FIGURES_DIR, exist_ok=True)

RADIUS   = 50    # µm neighbourhood radius
N_PERMS  = 200   # permutations for enrichment score

# ── Load ──────────────────────────────────────────────────────────────────────
print("Loading data...")
adata = sc.read_h5ad("results/08_cell_types/xenium_annotated.h5ad")
adata_z = sc.read_h5ad("results/02_nerve_id/xenium_nerve_zones.h5ad")
zone_map = adata_z.obs["spatial_zone"].to_dict()
adata.obs["spatial_zone"] = adata.obs_names.map(zone_map).fillna("Unknown")

print(f"  {adata.shape[0]:,} cells")
print(f"  Cell types: {adata.obs['cell_type'].value_counts().to_dict()}")

coords   = adata.obsm["spatial"]
ct_array = adata.obs["cell_type"].values
hpv_arr  = adata.obs["hpv_status"].values
cell_types = [ct for ct in adata.obs["cell_type"].unique()
              if ct not in ["Unknown"]]

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 1 — Manual neighbourhood enrichment (nerve-centric)
# For each nerve cell, find all cells within RADIUS µm
# Count cell types observed vs expected by random
# ══════════════════════════════════════════════════════════════════════════════
print(f"\nAnalysis 1: Neighbourhood enrichment around nerve cells (r={RADIUS}µm)...")

nerve_mask = adata.obs["spatial_zone"].values == "Nerve"
nerve_coords  = coords[nerve_mask]
other_coords  = coords[~nerve_mask]
other_ct      = ct_array[~nerve_mask]
other_hpv     = hpv_arr[~nerve_mask]

print(f"  Nerve cells: {nerve_mask.sum():,}")
print(f"  Other cells: {(~nerve_mask).sum():,}")

# Build KD-tree for non-nerve cells
tree = cKDTree(other_coords)

# Find neighbours of each nerve cell
print(f"  Querying neighbours within {RADIUS}µm...")
neighbour_cts = []
for nc in nerve_coords:
    idxs = tree.query_ball_point(nc, r=RADIUS)
    if idxs:
        neighbour_cts.extend(other_ct[idxs])

from collections import Counter
observed_counts = Counter(neighbour_cts)
total_neighbours = sum(observed_counts.values())
print(f"  Total neighbour observations: {total_neighbours:,}")

# Expected: global cell type proportions × total neighbours
global_props = adata.obs["cell_type"].value_counts(normalize=True)
expected_counts = {ct: global_props.get(ct, 0) * total_neighbours
                   for ct in cell_types}

# Enrichment score = log2(observed/expected)
enrichment = {}
for ct in cell_types:
    obs = observed_counts.get(ct, 0)
    exp = expected_counts.get(ct, 1)
    enrichment[ct] = np.log2((obs + 1) / (exp + 1))

# Permutation test — shuffle cell type labels
print(f"  Running {N_PERMS} permutations...")
perm_enrichments = {ct: [] for ct in cell_types}
all_ct = other_ct.copy()

for perm_i in range(N_PERMS):
    if perm_i % 50 == 0:
        print(f"    Permutation {perm_i}/{N_PERMS}")
    shuffled = np.random.permutation(all_ct)
    perm_obs = Counter()
    for nc in nerve_coords:
        idxs = tree.query_ball_point(nc, r=RADIUS)
        if idxs:
            perm_obs.update(shuffled[idxs])
    perm_total = sum(perm_obs.values())
    for ct in cell_types:
        obs_p = perm_obs.get(ct, 0)
        exp_p = global_props.get(ct, 0) * perm_total
        perm_enrichments[ct].append(np.log2((obs_p + 1) / (exp_p + 1)))

# Z-scores and p-values from permutations
results = []
for ct in cell_types:
    perms = np.array(perm_enrichments[ct])
    obs   = enrichment[ct]
    z     = (obs - perms.mean()) / (perms.std() + 1e-10)
    p     = (np.abs(perms) >= np.abs(obs)).mean()
    results.append({
        "cell_type" : ct,
        "enrichment": round(obs, 3),
        "z_score"   : round(z, 3),
        "p_value"   : round(p, 4),
        "n_observed": observed_counts.get(ct, 0),
        "n_expected": round(expected_counts.get(ct, 0), 1),
        "direction" : "Enriched" if obs > 0 else "Depleted",
    })

res_df = pd.DataFrame(results).sort_values("z_score", ascending=False)
res_df.to_csv(f"{OUTPUT_DIR}/neighbourhood_enrichment.csv", index=False)
print("\nNeighbourhood enrichment results:")
print(res_df[["cell_type","enrichment","z_score","p_value","direction"]].to_string(index=False))

# ── Figure 1: Neighbourhood enrichment bar plot ────────────────────────────────
fig, ax = plt.subplots(figsize=(12, 7))
colors = ["#e74c3c" if v > 0 else "#3498db" for v in res_df["enrichment"]]
bars = ax.barh(range(len(res_df)), res_df["enrichment"],
               color=colors, alpha=0.85, edgecolor="black", linewidth=0.6)
ax.set_yticks(range(len(res_df)))
ax.set_yticklabels(res_df["cell_type"].tolist(), fontsize=11)

# Labels outside bars with enough offset to avoid overlap
xmax = res_df["enrichment"].abs().max()
for i, (_, row) in enumerate(res_df.iterrows()):
    p    = row["p_value"]
    sig  = "**" if p < 0.01 else ("*" if p < 0.05 else "")
    x    = row["enrichment"]
    # Fix sign: use actual enrichment direction for label
    direction_label = "Enriched" if x > 0 else "Depleted"
    label = f"z={row['z_score']:.1f}  {sig}"
    offset = xmax * 0.05
    ha = "left" if x >= 0 else "right"
    xpos = x + offset if x >= 0 else x - offset
    ax.text(xpos, i, label, va="center", ha=ha, fontsize=9,
            fontweight="bold" if sig else "normal")

ax.axvline(0, color="black", linewidth=1.2, linestyle="--", alpha=0.5)
ax.set_xlim(-xmax*1.6, xmax*1.6)
ax.set_xlabel("Log2 Enrichment (observed/expected)\nnear Schwann/Nerve cells  (r=50\u00b5m)", fontsize=11)
ax.set_title("Cell Type Neighbourhood Enrichment Around Nerve Cells\n"
             f"(radius=50\u00b5m, n=200 permutations)  **p<0.01, *p<0.05",
             fontsize=12, fontweight="bold")
ax.spines[["top","right"]].set_visible(False)
plt.tight_layout()
for path in [f"{OUTPUT_DIR}/fig_nhood_enrichment.png",
             f"{FIGURES_DIR}/fig_nhood_enrichment.png"]:
    plt.savefig(path, dpi=150, bbox_inches="tight")
plt.close()
print(f"\n  Saved fig_nhood_enrichment.png")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 2 — NK cell distance to nearest nerve by HPV status
# ══════════════════════════════════════════════════════════════════════════════
print("\nAnalysis 2: NK cell distance to nearest nerve by HPV status...")

nerve_tree = cKDTree(nerve_coords)

for ct_name, ct_label in [("DC_LAMP3","DC_LAMP3 cells"), ("T_cell","T cells"),
                            ("Tumour","Tumour cells")]:
    ct_mask   = ct_array == ct_name
    ct_coords = coords[ct_mask]
    ct_hpv    = hpv_arr[ct_mask]

    if ct_mask.sum() == 0:
        continue

    dists, _ = nerve_tree.query(ct_coords, k=1)

    hpv_pos = dists[ct_hpv == "HPV+"]
    hpv_neg = dists[ct_hpv == "HPV-"]

    stat, p = stats.mannwhitneyu(hpv_pos, hpv_neg, alternative="two-sided")
    pstr = f"p={p:.4f}" if p >= 0.0001 else "p<0.0001"
    print(f"  {ct_name}: HPV+ median={np.median(hpv_pos):.1f}µm "
          f"vs HPV- median={np.median(hpv_neg):.1f}µm  {pstr}")

# ── Figure 2: Distance distributions ──────────────────────────────────────────
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
for ax, (ct_name, ct_label) in zip(axes, [
    ("DC_LAMP3","DC_LAMP3 cells"), ("T_cell","T cells"), ("Tumour","Tumour cells")
]):
    ct_mask   = ct_array == ct_name
    ct_coords = coords[ct_mask]
    ct_hpv    = hpv_arr[ct_mask]
    if ct_mask.sum() == 0:
        ax.text(0.5,0.5,"No data",transform=ax.transAxes,ha="center")
        continue

    dists, _ = nerve_tree.query(ct_coords, k=1)
    hpv_pos  = dists[ct_hpv == "HPV+"]
    hpv_neg  = dists[ct_hpv == "HPV-"]

    bins = np.linspace(0, min(500, np.percentile(dists, 95)), 40)
    ax.hist(hpv_pos, bins=bins, alpha=0.6, color="#e74c3c",
            density=True, label=f"HPV+ (n={len(hpv_pos):,})")
    ax.hist(hpv_neg, bins=bins, alpha=0.6, color="#2ecc71",
            density=True, label=f"HPV\u2212 (n={len(hpv_neg):,})")

    stat, p = stats.mannwhitneyu(hpv_pos, hpv_neg, alternative="two-sided")
    pstr = f"p={p:.4f}" if p >= 0.0001 else "p<0.0001"
    ax.text(0.97, 0.97, f"MWU {pstr}", transform=ax.transAxes,
            ha="right", va="top", fontsize=9,
            bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.7))
    ax.set_xlabel("Distance to nearest nerve (\u00b5m)", fontsize=10)
    ax.set_ylabel("Density", fontsize=10)
    ax.set_title(ct_label, fontsize=11, fontweight="bold")
    ax.legend(fontsize=9)
    ax.spines[["top","right"]].set_visible(False)

plt.suptitle("Cell-to-Nerve Distance Distributions by HPV Status",
             fontsize=12, fontweight="bold")
plt.tight_layout()
for path in [f"{OUTPUT_DIR}/fig_nerve_distance_by_HPV.png",
             f"{FIGURES_DIR}/fig_nerve_distance_by_HPV.png"]:
    plt.savefig(path, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved fig_nerve_distance_by_HPV.png")

# ══════════════════════════════════════════════════════════════════════════════
# ANALYSIS 3 — Neighbourhood composition heatmap
# What cell types co-occur within RADIUS µm of each cell type?
# ══════════════════════════════════════════════════════════════════════════════
print("\nAnalysis 3: Neighbourhood composition heatmap...")

# For each cell type, what is the mean composition of its neighbourhood?
all_tree = cKDTree(coords)
n_ct     = len(cell_types)
ct_to_idx = {ct: i for i, ct in enumerate(cell_types)}
comp_matrix = np.zeros((n_ct, n_ct))

# Sample max 2000 cells per cell type for speed
np.random.seed(42)
for i, source_ct in enumerate(cell_types):
    src_mask  = ct_array == source_ct
    src_idx   = np.where(src_mask)[0]
    sample    = src_idx if len(src_idx) <= 2000 else np.random.choice(src_idx, 2000, replace=False)
    neighbour_counts = np.zeros(n_ct)
    for idx in sample:
        neighbours = all_tree.query_ball_point(coords[idx], r=RADIUS)
        neighbours = [n for n in neighbours if n != idx]
        for n_idx in neighbours:
            n_ct_name = ct_array[n_idx]
            if n_ct_name in ct_to_idx:
                neighbour_counts[ct_to_idx[n_ct_name]] += 1
    total = neighbour_counts.sum()
    if total > 0:
        comp_matrix[i] = neighbour_counts / total

# Subtract global expected (global props)
for j, ct in enumerate(cell_types):
    comp_matrix[:, j] -= global_props.get(ct, 0)

comp_df = pd.DataFrame(comp_matrix, index=cell_types, columns=cell_types)
comp_df.to_csv(f"{OUTPUT_DIR}/neighbourhood_composition.csv")

# Plot
fig, ax = plt.subplots(figsize=(11, 9))
vmax = np.abs(comp_matrix).max()
im = ax.imshow(comp_matrix, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")

ax.set_xticks(range(n_ct)); ax.set_xticklabels(cell_types, rotation=45, ha="right", fontsize=9)
ax.set_yticks(range(n_ct)); ax.set_yticklabels(cell_types, fontsize=9)

for i in range(n_ct):
    for j in range(n_ct):
        val = comp_matrix[i, j]
        ax.text(j, i, f"{val:.3f}", ha="center", va="center",
                fontsize=7, color="white" if abs(val) > vmax*0.6 else "black")

plt.colorbar(im, ax=ax, label="Co-occurrence enrichment\n(observed − expected proportion)",
             shrink=0.8)
ax.set_title(f"Spatial Neighbourhood Composition\n"
             f"(enrichment of column cell type in neighbourhood of row cell type, r={RADIUS}\u00b5m)",
             fontsize=11, fontweight="bold")
ax.set_xlabel("Neighbour cell type", fontsize=10)
ax.set_ylabel("Source cell type", fontsize=10)
plt.tight_layout()
for path in [f"{OUTPUT_DIR}/fig_neighbourhood_composition.png",
             f"{FIGURES_DIR}/fig_neighbourhood_composition.png"]:
    plt.savefig(path, dpi=150, bbox_inches="tight")
plt.close()
print(f"  Saved fig_neighbourhood_composition.png")

# ── Summary ───────────────────────────────────────────────────────────────────
print(f"\n{'='*55}")
print("Module 13 complete")
print(f"{'='*55}")
print("\nKey neighbourhood enrichment findings:")
for _, row in res_df.iterrows():
    sig = "**" if row['p_value'] < 0.01 else ("*" if row['p_value'] < 0.05 else "n.s.")
    print(f"  {row['cell_type']:20s} z={row['z_score']:+.2f}  "
          f"({row['direction']})  {sig}")
print(f"\nFigures saved to: {OUTPUT_DIR}/")
