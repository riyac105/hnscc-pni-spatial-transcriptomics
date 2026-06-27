"""
11_replicate_concordance.py (v3 — fixed folder paths)
Validates PNI metric reproducibility across technical replicates.

Folder structure is: data/GSE300147/GSM9054472/Supp_GSM9054472_Patient_9__run_1/
Script now searches for Supp_* subdirectory automatically.

Run from ~/Desktop/PNI_project/
"""

import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
warnings.filterwarnings('ignore')

sc.settings.verbosity = 0
OUTPUT_DIR = "results/11_replicate_concordance"
os.makedirs(OUTPUT_DIR, exist_ok=True)

NERVE_THRESHOLD   = 0.8
PERINEURAL_RADIUS = 7.0
NERVE_MARKERS     = ["PMP22", "EDNRB", "PTN", "LGI4"]
TUMOUR_MARKERS    = ["EGFR", "EPCAM", "SOX2", "KRT7", "CLCA2", "SERPINB3"]

REPLICATES = {
    "P9"  : {"run1": "GSM9054472", "run2": "GSM9054473", "hpv": "HPV+"},
    "P12" : {"run1": "GSM9054474", "run2": "GSM9054475", "hpv": "HPV-"},
    "P13" : {"run1": "GSM9054476", "run2": "GSM9054477", "hpv": "HPV+"},
    "P17" : {"run1": "GSM9054479", "run2": "GSM9054480", "hpv": "HPV-"},
    "P19" : {"run1": "GSM9054481", "run2": "GSM9054482", "hpv": "HPV+"},
    "P20" : {"run1": "GSM9054483", "run2": "GSM9054484", "hpv": "HPV-"},
    "P28" : {"run1": "GSM9054486", "run2": "GSM9054487", "hpv": "HPV-"},
}

def find_xenium_dirs(gsm_id):
    """Find the Xenium data directories for a GSM ID.
    Handles folder structure: data/GSE300147/GSM_ID/Supp_GSM_ID_*/
    Returns (matrix_path, coords_path) or (None, None)
    """
    base = f"data/GSE300147/{gsm_id}"
    if not os.path.exists(base):
        return None, None

    # Look for Supp_* subdirectory
    supp_dir = None
    for entry in os.listdir(base):
        if entry.startswith("Supp_") and os.path.isdir(os.path.join(base, entry)):
            supp_dir = os.path.join(base, entry)
            break

    if supp_dir is None:
        print(f"    No Supp_* directory found in {base}")
        print(f"    Contents: {os.listdir(base)}")
        return None, None

    print(f"    Found: {supp_dir}")

    # Find matrix.mtx.gz — search recursively
    matrix_path = None
    coords_path = None
    for root, dirs, files in os.walk(supp_dir):
        if "matrix.mtx.gz" in files or "matrix.mtx" in files:
            matrix_path = root
        for fname in ["cells.csv.gz", "cells.csv", "cell_metadata.csv"]:
            if fname in files:
                coords_path = os.path.join(root, fname)

    if matrix_path is None:
        # Sometimes the matrix is directly in the supp dir
        # Check for barcodes.tsv.gz which indicates a 10x matrix
        for root, dirs, files in os.walk(supp_dir):
            if "barcodes.tsv.gz" in files or "barcodes.tsv" in files:
                matrix_path = root
                break

    return matrix_path, coords_path

def compute_pni_metrics(gsm_id, patient, run_label):
    print(f"  Processing {gsm_id} ({patient} {run_label})...")
    matrix_path, coords_path = find_xenium_dirs(gsm_id)

    if matrix_path is None:
        print(f"    No matrix found — skipping")
        return None

    try:
        adata = sc.read_10x_mtx(matrix_path, var_names="gene_symbols",
                                  cache=False, gex_only=True)
        print(f"    Loaded: {adata.shape}")
    except Exception as e:
        print(f"    Load error: {e}")
        return None

    # QC
    sc.pp.filter_cells(adata, min_genes=3)
    sc.pp.filter_cells(adata, min_counts=5)

    # Coordinates
    if coords_path is None:
        # Search in base directory
        base = f"data/GSE300147/{gsm_id}"
        for root, dirs, files in os.walk(base):
            for fname in ["cells.csv.gz","cells.csv"]:
                if fname in files:
                    coords_path = os.path.join(root, fname)
                    break

    if coords_path is None:
        print(f"    No coordinates file found")
        return None

    coords_df = pd.read_csv(coords_path)
    x_col = next((c for c in ["x_centroid","cell_x","x"] if c in coords_df.columns), None)
    y_col = next((c for c in ["y_centroid","cell_y","y"] if c in coords_df.columns), None)
    id_col = coords_df.columns[0]

    if x_col is None or y_col is None:
        print(f"    Coordinate columns not found. Available: {coords_df.columns.tolist()}")
        return None

    coords_df = coords_df.set_index(id_col)
    common = adata.obs_names.intersection(coords_df.index.astype(str))
    if len(common) == 0:
        coords_df.index = coords_df.index.astype(str)
        common = adata.obs_names.intersection(coords_df.index)
    if len(common) == 0:
        print(f"    Cannot match cell IDs between matrix and coords")
        return None

    adata = adata[common]
    adata.obsm["spatial"] = coords_df.loc[common, [x_col, y_col]].values

    # Score
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    nerve_genes  = [g for g in NERVE_MARKERS  if g in adata.var_names]
    tumour_genes = [g for g in TUMOUR_MARKERS if g in adata.var_names]

    if len(nerve_genes) < 2:
        print(f"    Too few nerve genes: {nerve_genes}")
        return None

    sc.tl.score_genes(adata, nerve_genes,  score_name="nerve_score",  use_raw=False)
    sc.tl.score_genes(adata, tumour_genes, score_name="tumour_score", use_raw=False)

    adata.obs["is_nerve"]  = (adata.obs["nerve_score"]  > NERVE_THRESHOLD).astype(int)
    adata.obs["is_tumour"] = (adata.obs["tumour_score"] > 0.1).astype(int)

    coords   = adata.obsm["spatial"]
    nerve_c  = coords[adata.obs["is_nerve"]==1]
    if len(nerve_c) < 5:
        print(f"    Too few nerve cells ({len(nerve_c)})")
        return None

    tree = cKDTree(nerve_c)
    dist, _ = tree.query(coords, k=1)
    adata.obs["dist_to_nerve"] = dist
    adata.obs["zone"] = np.where(
        adata.obs["is_nerve"]==1, "Nerve",
        np.where(adata.obs["dist_to_nerve"]<=PERINEURAL_RADIUS, "Perineural", "Other")
    )

    n_nerve   = int((adata.obs["is_nerve"]==1).sum())
    n_pni     = int(((adata.obs["zone"]=="Perineural") &
                     (adata.obs["is_tumour"]==1)).sum())
    n_cells   = len(adata)
    pni_index = n_pni / max(n_nerve, 1)

    print(f"    {n_cells:,} cells | nerve={100*n_nerve/n_cells:.1f}% | PNI_index={pni_index:.3f}")
    return {
        "patient": patient, "run": run_label, "gsm": gsm_id,
        "hpv": REPLICATES[patient]["hpv"],
        "n_cells": n_cells,
        "pct_nerve": 100*n_nerve/n_cells,
        "pct_perineural": 100*(adata.obs["zone"]=="Perineural").sum()/n_cells,
        "pni_index": pni_index,
    }

# ── Run analysis ──────────────────────────────────────────────────────────────
print("Computing PNI metrics for all replicates...")
records = []
for patient, info in REPLICATES.items():
    print(f"\n{patient} ({info['hpv']}):")
    for run_label, gsm_id in [("Run1", info["run1"]), ("Run2", info["run2"])]:
        result = compute_pni_metrics(gsm_id, patient, run_label)
        if result:
            records.append(result)

if len(records) < 4:
    print(f"\nOnly {len(records)} samples loaded — falling back to h5ad consistency analysis")
    adata_z = sc.read_h5ad("results/02_nerve_id/xenium_nerve_zones.h5ad")
    fb = []
    for pt in adata_z.obs["patient_id"].unique():
        mask = adata_z.obs["patient_id"]==pt
        sub  = adata_z.obs[mask]
        n_n  = (sub["is_nerve"]==1).sum()
        n_p  = (sub["PNI_positive"]==1).sum() if "PNI_positive" in sub.columns else 0
        fb.append({
            "patient": pt, "hpv": sub["hpv_status"].iloc[0],
            "n_cells": len(sub),
            "pct_nerve": 100*n_n/len(sub),
            "pct_perineural": 100*(sub["spatial_zone"]=="Perineural").sum()/len(sub),
            "pni_index": n_p/max(n_n,1),
        })
    df_fb = pd.DataFrame(fb)
    print(df_fb.to_string(index=False))
    df_fb.to_csv(f"{OUTPUT_DIR}/patient_metrics.csv", index=False)

    hpv_colors = {"HPV+":"#e74c3c","HPV-":"#2ecc71"}
    fig, axes = plt.subplots(1,3,figsize=(15,5))
    for ax,(col,label) in zip(axes,[
        ("pni_index","PNI Index"),
        ("pct_nerve","Nerve Density (%)"),
        ("pct_perineural","Perineural Zone (%)")
    ]):
        for hpv, grp in df_fb.groupby("hpv"):
            xs = [list(hpv_colors.keys()).index(hpv)] * len(grp)
            ax.scatter(xs, grp[col], color=hpv_colors[hpv],
                       s=120, edgecolors="black", linewidth=0.8, zorder=5)
            for _, row in grp.iterrows():
                ax.annotate(row["patient"], (xs[0], row[col]),
                            textcoords="offset points", xytext=(8,0), fontsize=8)
            ax.hlines(grp[col].mean(),
                      list(hpv_colors.keys()).index(hpv)-0.25,
                      list(hpv_colors.keys()).index(hpv)+0.25,
                      colors=hpv_colors[hpv], linewidth=2.5, linestyle="--")
        ax.set_xticks([0,1]); ax.set_xticklabels(["HPV+","HPV-"],fontsize=11)
        ax.set_ylabel(label,fontsize=11); ax.set_title(label,fontsize=11,fontweight="bold")

    plt.suptitle("Per-Patient PNI Metrics by HPV Status",fontsize=13,fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/fig_patient_consistency.png",dpi=150,bbox_inches="tight")
    plt.close()
    print(f"  Saved → {OUTPUT_DIR}/fig_patient_consistency.png")

else:
    df = pd.DataFrame(records)
    df.to_csv(f"{OUTPUT_DIR}/replicate_metrics.csv", index=False)
    print(f"\nLoaded {len(records)} samples:")
    print(df.to_string(index=False))

    run1 = df[df["run"]=="Run1"].set_index("patient")
    run2 = df[df["run"]=="Run2"].set_index("patient")
    common = run1.index.intersection(run2.index)
    print(f"\nMatched pairs: {list(common)}")

    metrics = [("pni_index","PNI Index"),
               ("pct_nerve","Nerve Density (%)"),
               ("pct_perineural","Perineural Zone (%)")]
    hpv_colors = {"HPV+":"#e74c3c","HPV-":"#2ecc71"}
    fig, axes = plt.subplots(1,3,figsize=(16,5))
    for ax,(col,label) in zip(axes,metrics):
        x = run1.loc[common,col].values.astype(float)
        y = run2.loc[common,col].values.astype(float)
        for i,pt in enumerate(common):
            hpv = run1.loc[pt,"hpv"]
            ax.scatter(x[i],y[i],color=hpv_colors[hpv],s=120,
                       edgecolors="black",linewidth=0.8,zorder=5)
            ax.annotate(pt,(x[i],y[i]),textcoords="offset points",
                        xytext=(6,4),fontsize=9)
        if len(x)>=3:
            r,p = stats.pearsonr(x,y)
            mn = min(x.min(),y.min())*0.9; mx = max(x.max(),y.max())*1.1
            ax.plot([mn,mx],[mn,mx],"k--",linewidth=1.5,alpha=0.5)
            ax.set_xlim(mn,mx); ax.set_ylim(mn,mx)
            ax.text(0.05,0.92,f"r={r:.2f}\np={p:.3f}",
                    transform=ax.transAxes,fontsize=10,
                    bbox=dict(boxstyle="round",facecolor="wheat",alpha=0.6))
            print(f"  {label}: r={r:.2f}, p={p:.3f}")
        ax.set_xlabel("Run 1",fontsize=11); ax.set_ylabel("Run 2",fontsize=11)
        ax.set_title(label,fontsize=12,fontweight="bold"); ax.set_aspect("equal")

    legend_patches=[mpatches.Patch(color=c,label=h) for h,c in hpv_colors.items()]
    fig.legend(handles=legend_patches,fontsize=10,loc="lower center",
               ncol=2,bbox_to_anchor=(0.5,-0.05))
    plt.suptitle("Technical Replicate Concordance (Run 1 vs Run 2)",
                 fontsize=14,fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/fig_replicate_concordance.png",dpi=150,bbox_inches="tight")
    plt.close()
    print(f"  Saved → {OUTPUT_DIR}/fig_replicate_concordance.png")

print(f"\n{'='*60}")
print("Module 11 complete")
print(f"{'='*60}")
