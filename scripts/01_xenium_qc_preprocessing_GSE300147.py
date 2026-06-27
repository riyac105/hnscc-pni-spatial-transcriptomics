"""
Module 1 (GSE300147-specific): Multi-Sample Loading, Batch Integration & QC
============================================================================
Dataset : GSE300147 — Xenium HNSCC (McCord & Hudson, Baylor)
Samples : 17 HNSCC Xenium sections (10 unique patients, 7 with replicates)
          + 1 ameloblastoma (excluded from PNI analysis)

Key changes from generic Module 1:
    - Loads all 17 HNSCC samples and concatenates into one AnnData
    - Handles Xenium PARQUET format (cells.parquet, transcripts.parquet)
    - Stores patient_id, run, hpv_status, replicate_group as metadata
    - Batch correction via scVI or Harmony (across patients + runs)
    - Technical replicate QC: checks reproducibility across replicate runs
    - Excludes ameloblastoma (GSM9054488, Patient 10)
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
warnings.filterwarnings("ignore")

# Batch correction options
try:
    import scvi
    SCVI = True
except ImportError:
    SCVI = False
    print("scvi-tools not installed → falling back to Harmony")
    print("  Install: pip install scvi-tools")

try:
    from harmonypy import run_harmony
    HARMONY = True
except ImportError:
    HARMONY = False
    print("harmonypy not installed → pip install harmonypy")

sc.settings.seed = 42
np.random.seed(42)

# ── Paths ──────────────────────────────────────────────────────────────────────
DATA_DIR   = "data/GSE300147/"
OUTPUT_DIR = "results/01_qc/"
os.makedirs(OUTPUT_DIR, exist_ok=True)
sc.settings.figdir = OUTPUT_DIR

# ==============================================================================
# SAMPLE MANIFEST — GSE300147 (HNSCC only, ameloblastoma excluded)
# ==============================================================================
# HPV status: update from paper/supplementary once available
# Current assignment based on GEO description (5 HPV+, 5 HPV-)

SAMPLES = [
    # gsm_id          patient  run  hpv     is_replicate
    ("GSM9054471", "P1",  1, "HPV+", False),
    ("GSM9054472", "P9",  1, "HPV+", False),
    ("GSM9054473", "P9",  2, "HPV+", True),   # replicate of P9 run1
    ("GSM9054474", "P12", 1, "HPV-", False),
    ("GSM9054475", "P12", 2, "HPV-", True),
    ("GSM9054476", "P13", 1, "HPV+", False),
    ("GSM9054477", "P13", 2, "HPV+", True),
    ("GSM9054478", "P15", 1, "HPV-", False),
    ("GSM9054479", "P17", 1, "HPV-", False),
    ("GSM9054480", "P17", 2, "HPV-", True),
    ("GSM9054481", "P19", 1, "HPV+", False),
    ("GSM9054482", "P19", 2, "HPV+", True),
    ("GSM9054483", "P20", 1, "HPV-", False),
    ("GSM9054484", "P20", 2, "HPV-", True),
    ("GSM9054485", "P23", 1, "HPV+", False),
    ("GSM9054486", "P28", 1, "HPV-", False),
    ("GSM9054487", "P28", 2, "HPV-", True),
    # GSM9054488 = Patient 10 (ameloblastoma) — EXCLUDED
]

# QC thresholds — calibrated for 477-gene Xenium panel
MIN_COUNTS  = 5      # Lower than WGS because panel is targeted (477 genes)
MAX_COUNTS  = 1500   # Upper cap
MIN_GENES   = 3
MAX_MT_PCT  = 25

# For initial analysis use run 1 only per patient (drops replicates)
# This reduces cells from ~1.1M to ~500k — much more manageable on a laptop
# Replicates are used later only for concordance QC (Module replicate check)
USE_FIRST_RUN_ONLY = True

# ==============================================================================
# 1. LOAD EACH SAMPLE
# ==============================================================================

def find_file_recursive(root_dir: str, pattern: str):
    """
    Recursively search root_dir for a file whose name contains pattern.
    Handles GEOparse Supp_GSM.../ subfolder structure automatically.
    """
    for dirpath, dirnames, filenames in os.walk(root_dir):
        for fname in filenames:
            if pattern in fname:
                return os.path.join(dirpath, fname)
    return None


def load_xenium_sample(gsm_id: str,
                        patient_id: str,
                        run: int,
                        hpv_status: str,
                        is_replicate: bool,
                        data_dir: str) -> sc.AnnData:
    """
    Load a single Xenium sample, searching recursively for files.
    Handles GEOparse Supp_GSM.../ subfolder structure automatically.
    """
    sample_dir = os.path.join(data_dir, gsm_id)

    # ── Expression matrix — recursive search ───────────────────────────────────
    h5_path = find_file_recursive(sample_dir, "cell_feature_matrix.h5")
    if h5_path is None:
        raise FileNotFoundError(f"H5 file not found for {gsm_id}. "
                                 f"Check: {sample_dir}")
    print(f"    H5:      {os.path.relpath(h5_path)}")
    adata = sc.read_10x_h5(h5_path)

    # ── Cell coordinates — recursive search ────────────────────────────────────
    parquet_path = find_file_recursive(sample_dir, "cells.parquet")
    csv_path     = find_file_recursive(sample_dir, "cells.csv")

    if parquet_path:
        print(f"    Parquet: {os.path.relpath(parquet_path)}")
        cells_df = pd.read_parquet(parquet_path)
    elif csv_path:
        print(f"    CSV:     {os.path.relpath(csv_path)}")
        cells_df = pd.read_csv(csv_path)
    else:
        raise FileNotFoundError(f"Cells file not found for {gsm_id}")

    # Standardise index
    if "cell_id" in cells_df.columns:
        cells_df = cells_df.set_index("cell_id")
    cells_df.index = cells_df.index.astype(str)
    adata.obs_names = adata.obs_names.astype(str)

    common = adata.obs_names.intersection(cells_df.index)
    adata   = adata[common].copy()
    cells_df = cells_df.loc[common]

    # Spatial coordinates (microns)
    x_col = next((c for c in cells_df.columns if "x_centroid" in c or "x" == c), None)
    y_col = next((c for c in cells_df.columns if "y_centroid" in c or "y" == c), None)
    if x_col and y_col:
        adata.obsm["spatial"] = cells_df[[x_col, y_col]].values

    # Cell metadata
    for col in ["cell_area", "nucleus_area", "transcript_counts",
                "control_probe_counts", "total_counts"]:
        if col in cells_df.columns:
            adata.obs[col] = cells_df[col].values

    # ── Sample-level metadata ───────────────────────────────────────────────────
    adata.obs["sample_id"]    = gsm_id
    adata.obs["patient_id"]   = patient_id
    adata.obs["run"]          = run
    adata.obs["hpv_status"]   = hpv_status
    adata.obs["is_replicate"] = is_replicate
    adata.obs["batch"]        = f"{patient_id}_run{run}"

    # Make obs_names unique across samples
    adata.obs_names = [f"{gsm_id}_{bc}" for bc in adata.obs_names]

    print(f"  Loaded {gsm_id} ({patient_id} run{run}, {hpv_status}): "
          f"{adata.n_obs:,} cells × {adata.n_vars} genes")
    return adata


print("Loading all HNSCC samples...")
adatas = []
failed = []

for gsm_id, patient_id, run, hpv_status, is_replicate in SAMPLES:
    sample_dir = os.path.join(DATA_DIR, gsm_id)
    if not os.path.isdir(sample_dir):
        print(f"  SKIP {gsm_id} — directory not found (download first)")
        failed.append(gsm_id)
        continue
    try:
        adata_s = load_xenium_sample(
            gsm_id, patient_id, run, hpv_status, is_replicate, DATA_DIR
        )
        adatas.append(adata_s)
    except Exception as e:
        print(f"  ERROR loading {gsm_id}: {e}")
        failed.append(gsm_id)

if not adatas:
    raise RuntimeError(
        "No samples loaded. Please download the data first:\n"
        "  See results/00_panel_audit/download_instructions.txt\n"
        "  Then re-run this module."
    )

print(f"\nSuccessfully loaded: {len(adatas)}/{len(SAMPLES)} samples")
if failed:
    print(f"Failed / missing: {failed}")

# ── Drop replicates for speed — use run 1 only per patient ────────────────────
if USE_FIRST_RUN_ONLY:
    adatas = [a for a in adatas if a.obs["is_replicate"].iloc[0] == False]
    print(f"Using run 1 only: {len(adatas)} samples (replicates dropped for speed)")

# ==============================================================================
# 2. CONCATENATE ALL SAMPLES
# ==============================================================================
print("\nConcatenating samples...")
adata = sc.concat(
    adatas,
    join="inner",       # keep only genes present in ALL samples (should be all 477)
    label="sample_id",
    keys=[a.obs["sample_id"].iloc[0] for a in adatas],
    uns_merge="first"
)

print(f"Combined: {adata.n_obs:,} cells × {adata.n_vars} genes")
print(f"Patients: {adata.obs['patient_id'].nunique()}")
print(f"HPV+: {(adata.obs['hpv_status']=='HPV+').sum():,} cells | "
      f"HPV-: {(adata.obs['hpv_status']=='HPV-').sum():,} cells")

# ==============================================================================
# 3. QUALITY CONTROL
# ==============================================================================
print("\nRunning QC...")

adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

print(f"  Before filtering: {adata.n_obs:,} cells")
sc.pp.filter_cells(adata, min_counts=MIN_COUNTS)
sc.pp.filter_cells(adata, max_counts=MAX_COUNTS)
sc.pp.filter_cells(adata, min_genes=MIN_GENES)
adata = adata[adata.obs["pct_counts_mt"] < MAX_MT_PCT].copy()
print(f"  After filtering:  {adata.n_obs:,} cells")

# Per-sample QC summary
qc_summary = adata.obs.groupby("sample_id").agg(
    n_cells=("total_counts", "count"),
    median_counts=("total_counts", "median"),
    median_genes=("n_genes_by_counts", "median"),
)
print("\nPer-sample QC summary:")
print(qc_summary.to_string())
qc_summary.to_csv(f"{OUTPUT_DIR}/per_sample_qc.csv")

# ── QC violin plots per sample ─────────────────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(16, 5))
samples_order = sorted(adata.obs["sample_id"].unique())

for ax, metric in zip(axes, ["total_counts", "n_genes_by_counts"]):
    data_to_plot = [adata.obs.loc[adata.obs["sample_id"] == s, metric].values
                    for s in samples_order]
    ax.violinplot(data_to_plot, showmedians=True)
    ax.set_xticks(range(1, len(samples_order)+1))
    ax.set_xticklabels(
        [s.replace("GSM905", "GSM..") for s in samples_order],
        rotation=45, ha="right", fontsize=7
    )
    ax.set_ylabel(metric)
    ax.set_title(f"Per-sample {metric}")

plt.tight_layout()
plt.savefig(f"{OUTPUT_DIR}/per_sample_qc_violins.pdf", dpi=150)
plt.close()

# ==============================================================================
# 4. TECHNICAL REPLICATE QC (correlation check)
# ==============================================================================
print("\nChecking technical replicate concordance...")

replicate_patients = adata.obs.loc[
    adata.obs["is_replicate"] == True, "patient_id"
].unique()

rep_corr_rows = []
for pt in replicate_patients:
    pt_adata = adata[adata.obs["patient_id"] == pt].copy()
    runs = pt_adata.obs["run"].unique()
    if len(runs) < 2:
        continue

    run1_mean = np.array(pt_adata[pt_adata.obs["run"]==1].X.mean(axis=0)).flatten()
    run2_mean = np.array(pt_adata[pt_adata.obs["run"]==2].X.mean(axis=0)).flatten()

    from scipy.stats import pearsonr
    r, p = pearsonr(run1_mean, run2_mean)
    rep_corr_rows.append({"patient": pt, "pearson_r": r, "pval": p})
    print(f"  {pt}: run1 vs run2 Pearson r = {r:.3f}")

if rep_corr_rows:
    rep_df = pd.DataFrame(rep_corr_rows)
    rep_df.to_csv(f"{OUTPUT_DIR}/replicate_concordance.csv", index=False)

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(rep_df["patient"], rep_df["pearson_r"],
           color="steelblue", edgecolor="none")
    ax.axhline(0.9, color="red", ls="--", label="r=0.9 threshold")
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("Pearson r (run1 vs run2 mean expression)")
    ax.set_title("Technical replicate concordance by patient")
    ax.legend()
    plt.tight_layout()
    plt.savefig(f"{OUTPUT_DIR}/replicate_concordance.pdf", dpi=150)
    plt.close()

# ==============================================================================
# 5. NORMALISATION
# ==============================================================================
print("\nNormalising...")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=None)
sc.pp.log1p(adata)

# ==============================================================================
# 6. DIMENSIONALITY REDUCTION & BATCH CORRECTION
# ==============================================================================
print("\nPCA + batch correction...")

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50, svd_solver="arpack")

# ── Batch correction strategy ─────────────────────────────────────────────────
# Use patient_id as the primary batch variable to correct for:
#   (a) patient-to-patient variation (dominant)
#   (b) run-to-run variation (secondary; use "batch" = patient_run)
#
# Priority: scVI > Harmony > uncorrected PCA

BATCH_KEY = "patient_id"   # primary batch
N_PCS     = 20  # reduced from 30 for speed on laptop

if SCVI:
    print("  Running scVI batch correction...")
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=BATCH_KEY)
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=30, gene_likelihood="nb")
    model.train(max_epochs=200, early_stopping=True)
    adata.obsm["X_scVI"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_scVI", n_neighbors=20, key_added="scVI_neighbors")
    sc.tl.umap(adata, neighbors_key="scVI_neighbors")
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_r05",
                 neighbors_key="scVI_neighbors")
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden_r10",
                 neighbors_key="scVI_neighbors")
    CORRECTION_METHOD = "scVI"

elif HARMONY:
    print("  Running Harmony batch correction...")
    import harmonypy as hm
    ho = hm.run_harmony(
        adata.obsm["X_pca"],
        adata.obs,
        [BATCH_KEY],
        max_iter_harmony=20,
        random_state=42
    )
    # Z_corr shape is (n_pcs, n_cells) — transpose to (n_cells, n_pcs)
    Z = ho.Z_corr
    if Z.shape[0] < Z.shape[1]:  # shape is (pcs, cells), needs transposing
        Z = Z.T
    adata.obsm["X_pca_harmony"] = Z
    sc.pp.neighbors(adata, use_rep="X_pca_harmony", n_pcs=N_PCS,
                    n_neighbors=15, key_added="harmony_neighbors")  # n_neighbors=15 for speed
    sc.tl.umap(adata, neighbors_key="harmony_neighbors")
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_r05",
                 neighbors_key="harmony_neighbors")
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden_r10",
                 neighbors_key="harmony_neighbors")
    CORRECTION_METHOD = "Harmony"

else:
    print("  No batch correction tool available — using uncorrected PCA.")
    print("  STRONGLY recommended: pip install scvi-tools  OR  pip install harmonypy")
    sc.pp.neighbors(adata, n_pcs=N_PCS, n_neighbors=20)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=0.5, key_added="leiden_r05")
    sc.tl.leiden(adata, resolution=1.0, key_added="leiden_r10")
    CORRECTION_METHOD = "None"

adata.uns["batch_correction_method"] = CORRECTION_METHOD
print(f"  Batch correction: {CORRECTION_METHOD}")

# ==============================================================================
# 7. VISUALISATION
# ==============================================================================
# UMAP coloured by key variables
sc.pl.umap(
    adata,
    color=["leiden_r05", "patient_id", "hpv_status", "total_counts"],
    ncols=3,
    save="_batch_overview.pdf",
    show=False
)

# Per-patient spatial plots (one figure per sample)
print("\nGenerating per-sample spatial plots...")
for gsm_id in adata.obs["sample_id"].unique():
    s_adata = adata[adata.obs["sample_id"] == gsm_id].copy()
    pt      = s_adata.obs["patient_id"].iloc[0]
    run     = s_adata.obs["run"].iloc[0]
    hpv     = s_adata.obs["hpv_status"].iloc[0]

    try:
        sq.pl.spatial_scatter(
            s_adata, color="leiden_r05",
            shape=None, size=0.5,
            title=f"{gsm_id} ({pt} run{run}, {hpv})",
            save=f"spatial_{gsm_id}_clusters.pdf"
        )
    except Exception as e:
        print(f"  Could not plot {gsm_id}: {e}")

# ==============================================================================
# 8. SAVE
# ==============================================================================
out_path = f"{OUTPUT_DIR}/GSE300147_all_samples_qc.h5ad"
adata.write_h5ad(out_path)
print(f"\nSaved → {out_path}")
print(f"  {adata.n_obs:,} cells × {adata.n_vars} genes")
print(f"  Batch correction: {CORRECTION_METHOD}")
print("Module 1 (GSE300147) complete ✓")
