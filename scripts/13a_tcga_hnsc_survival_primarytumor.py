"""
TCGA primary-tumor correction + REMARK reporting
-------------------------------------------------
Purpose:
1) Restrict TCGA RNA-seq to Primary Solid Tumor samples (sample type 01)
   BEFORE truncating barcodes to 12-character patient IDs.
2) Recompute the MDM2/PPARG/NFE2L2 score exactly as in the final Figure 6 workflow:
   gene-wise z-normalization across the TCGA primary-tumor cohort, then mean score.
3) Apply ONE median cutoff across the full scored TCGA cohort, retained after HPV stratification.
4) Re-run the 7-year HPV-stratified Kaplan-Meier analysis.
5) Print REMARK-ready analyzable N / event counts and univariate Cox effect estimates.
6) Save a corrected scored patient-level CSV and corrected KM figure without overwriting old files.

Run from:
    ~/Desktop/PNI_project/
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings("ignore")

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import mygene

OUTPUT_DIR = "results/06_TCGA_validation"
MAX_DAYS = 2555
TIMEPOINTS = [0, 500, 1000, 1500, 2000, 2555]
os.makedirs(OUTPUT_DIR, exist_ok=True)

def censor(df, tc, ec, mx):
    df = df.copy()
    beyond = df[tc] > mx
    df.loc[beyond, ec] = 0
    df.loc[beyond, tc] = mx
    return df

# -------------------------------------------------------------------------
# 1. Load patient-level clinical cohort and raw RNA-seq matrix
# -------------------------------------------------------------------------
print("Loading patient-level clinical data...")
clin = pd.read_csv(f"{OUTPUT_DIR}/TCGA_HPV_PNI_merged_clean.csv")

print("Loading TCGA RNA-seq counts...")
expr = pd.read_csv(
    "data/TCGA/TCGA-HNSC.star_counts.tsv.gz",
    sep="\t",
    index_col=0,
    compression="gzip"
)
expr.index = expr.index.astype(str).str.split(".").str[0]

# -------------------------------------------------------------------------
# 2. Restrict to PRIMARY SOLID TUMOR (01) BEFORE collapsing to patient ID
# -------------------------------------------------------------------------
original_cols = pd.Index(expr.columns.astype(str))
sample_types = original_cols.to_series(index=original_cols).str.split("-").str[3].str[:2]

primary_cols = sample_types[sample_types == "01"].index.tolist()
expr = expr.loc[:, primary_cols].copy()

# Now collapse barcode to 12-character patient ID.
expr.columns = [str(c)[:12] for c in expr.columns]

if len(expr.columns) != len(set(expr.columns)):
    dupes = pd.Series(expr.columns).value_counts()
    dupes = dupes[dupes > 1]
    raise RuntimeError(
        "Unexpected duplicate patient IDs remain among primary-tumor samples:\n"
        + dupes.to_string()
    )

print(f"Primary-tumor RNA-seq columns: {expr.shape[1]}")
print(f"Unique primary-tumor patients: {len(set(expr.columns))}")

# -------------------------------------------------------------------------
# 3. Recompute the 3-gene score
# -------------------------------------------------------------------------
genes = ["MDM2", "PPARG", "NFE2L2"]
mg = mygene.MyGeneInfo()
res = mg.querymany(
    genes,
    scopes="symbol",
    fields="ensembl.gene",
    species="human",
    returnall=False,
    verbose=False
)

sym_to_ens = {}
for r in res:
    if "ensembl" in r:
        ens = r["ensembl"]
        if isinstance(ens, list):
            ens = ens[0]
        if isinstance(ens, dict):
            sym_to_ens[r["query"]] = ens.get("gene", "")

found = {symbol: ens for symbol, ens in sym_to_ens.items() if ens in expr.index}
print("Mapped genes:", found)

if len(found) != 3:
    raise RuntimeError(f"Expected all 3 genes to map; found {len(found)}: {found}")

sub = expr.loc[list(found.values())].copy()
reverse_map = {v: k for k, v in found.items()}
sub.index = [reverse_map[i] for i in sub.index]

# Gene-wise z-normalization across the primary-tumor TCGA cohort
ez = sub.T
ez = (ez - ez.mean()) / ez.std()

score = ez.mean(axis=1)
score.name = "score"
score.index.name = "id12"

# -------------------------------------------------------------------------
# 4. Merge ONE expression score per patient into clinical data
# -------------------------------------------------------------------------
scored = clin.merge(score.reset_index(), on="id12", how="inner")

# One median across the full scored cohort, retained after HPV stratification
median_score = scored["score"].median()
scored["grp"] = (scored["score"] > median_score).map({True: "High", False: "Low"})

out_csv = f"{OUTPUT_DIR}/TCGA_threegene_primarytumor_scored.csv"
scored.to_csv(out_csv, index=False)

print(f"\nScored patient-level rows: {len(scored)}")
print("HPV counts after primary-tumor merge:")
print(scored["HPV"].value_counts(dropna=False).to_string())

# -------------------------------------------------------------------------
# 5. REMARK reporting + univariate effect estimates
# -------------------------------------------------------------------------
print("\n" + "=" * 78)
print("TCGA PRIMARY-TUMOR CORRECTED — REMARK REPORTING")
print("=" * 78)

for hpv in ["HPV-", "HPV+"]:
    dat = scored[scored["HPV"] == hpv].dropna(
        subset=["OS.time", "OS.event", "score", "grp"]
    ).copy()
    dat = censor(dat, "OS.time", "OS.event", MAX_DAYS)

    hi = dat[dat["grp"] == "High"].copy()
    lo = dat[dat["grp"] == "Low"].copy()

    lr = logrank_test(
        hi["OS.time"], lo["OS.time"],
        event_observed_A=hi["OS.event"],
        event_observed_B=lo["OS.event"]
    )

    print(f"\n{hpv}")
    print(f"  Analyzable N={len(dat)}")
    print(f"  7-year events={int(dat['OS.event'].sum())}")
    print(f"  Signature-high: N={len(hi)}, events={int(hi['OS.event'].sum())}")
    print(f"  Signature-low:  N={len(lo)}, events={int(lo['OS.event'].sum())}")
    print(f"  Log-rank p={lr.p_value:.8g}")

    # High vs low Cox model
    cox_group = dat[["OS.time", "OS.event", "grp"]].copy()
    cox_group["signature_high"] = (cox_group["grp"] == "High").astype(int)
    cox_group = cox_group[["OS.time", "OS.event", "signature_high"]].dropna()

    cph_group = CoxPHFitter()
    cph_group.fit(cox_group, duration_col="OS.time", event_col="OS.event")
    r1 = cph_group.summary.loc["signature_high"]

    print("  Cox high vs low:")
    print(f"    HR={r1['exp(coef)']:.3f}")
    print(
        f"    95% CI={r1['exp(coef) lower 95%']:.3f}–"
        f"{r1['exp(coef) upper 95%']:.3f}"
    )
    print(f"    p={r1['p']:.8g}")

    # Continuous score Cox model
    cox_cont = dat[["OS.time", "OS.event", "score"]].dropna().copy()
    cph_cont = CoxPHFitter()
    cph_cont.fit(cox_cont, duration_col="OS.time", event_col="OS.event")
    r2 = cph_cont.summary.loc["score"]

    print("  Cox continuous score:")
    print(f"    HR={r2['exp(coef)']:.3f}")
    print(
        f"    95% CI={r2['exp(coef) lower 95%']:.3f}–"
        f"{r2['exp(coef) upper 95%']:.3f}"
    )
    print(f"    p={r2['p']:.8g}")

# -------------------------------------------------------------------------
# 6. Save corrected TCGA KM figure
# -------------------------------------------------------------------------
fig = plt.figure(figsize=(14, 8))
fig.patch.set_facecolor("white")
gs = gridspec.GridSpec(
    1, 2, figure=fig,
    left=0.08, right=0.97, top=0.84, bottom=0.06, wspace=0.25
)

for ci, (hpv, subtitle) in enumerate([
    ("HPV-", "HPV− Patients"),
    ("HPV+", "HPV+ Patients — Negative Control")
]):
    gs_in = gridspec.GridSpecFromSubplotSpec(
        2, 1, subplot_spec=gs[ci],
        height_ratios=[4.5, 1], hspace=0.18
    )
    ax_km = fig.add_subplot(gs_in[0])
    ax_risk = fig.add_subplot(gs_in[1])

    dat = scored[scored["HPV"] == hpv].dropna(
        subset=["OS.time", "OS.event", "score", "grp"]
    ).copy()
    dat = censor(dat, "OS.time", "OS.event", MAX_DAYS)
    hi = dat[dat["grp"] == "High"]
    lo = dat[dat["grp"] == "Low"]

    KaplanMeierFitter().fit(
        hi["OS.time"], hi["OS.event"],
        label=f"Signature-high (n={len(hi)})"
    ).plot_survival_function(ax=ax_km, ci_show=True, linewidth=2)

    KaplanMeierFitter().fit(
        lo["OS.time"], lo["OS.event"],
        label=f"Signature-low (n={len(lo)})"
    ).plot_survival_function(ax=ax_km, ci_show=True, linewidth=2)

    lr = logrank_test(
        hi["OS.time"], lo["OS.time"],
        event_observed_A=hi["OS.event"],
        event_observed_B=lo["OS.event"]
    )

    ax_km.text(
        0.05, 0.08,
        f"Log-rank p = {lr.p_value:.4g}",
        transform=ax_km.transAxes,
        fontsize=10
    )
    ax_km.set_title(f"{subtitle} (n={len(dat)})", fontweight="bold")
    ax_km.set_ylabel("Overall Survival Probability")
    ax_km.set_xlabel("")
    ax_km.set_xlim(-50, MAX_DAYS + 80)
    ax_km.set_ylim(0, 1.05)
    ax_km.set_xticks(TIMEPOINTS)
    ax_km.tick_params(labelbottom=False)
    ax_km.spines[["top", "right"]].set_visible(False)

    hi_risk = [int((hi["OS.time"] >= t).sum()) for t in TIMEPOINTS]
    lo_risk = [int((lo["OS.time"] >= t).sum()) for t in TIMEPOINTS]

    ax_risk.set_xlim(ax_km.get_xlim())
    ax_risk.set_ylim(0, 1)
    ax_risk.axis("off")

    xlim = ax_km.get_xlim()
    xrange = xlim[1] - xlim[0]
    for t, nh, nl in zip(TIMEPOINTS, hi_risk, lo_risk):
        xpos = (t - xlim[0]) / xrange
        ax_risk.text(xpos, 0.72, str(nh), ha="center", va="center",
                     transform=ax_risk.transAxes, fontsize=9)
        ax_risk.text(xpos, 0.25, str(nl), ha="center", va="center",
                     transform=ax_risk.transAxes, fontsize=9)

    ax_risk.text(0.01, 0.72, "High", ha="left", va="center",
                 transform=ax_risk.transAxes, fontweight="bold", fontsize=9)
    ax_risk.text(0.01, 0.25, "Low", ha="left", va="center",
                 transform=ax_risk.transAxes, fontweight="bold", fontsize=9)

fig.suptitle(
    "TCGA-HNSC: NFE2L2/MDM2/PPARG Signature — Primary-Tumor RNA-seq, 7-Year OS",
    fontsize=13, fontweight="bold", y=0.97
)

out_fig = f"{OUTPUT_DIR}/fig_KM_threegene_primarytumor_corrected.png"
plt.savefig(out_fig, dpi=150, bbox_inches="tight", facecolor="white")
plt.close()

print("\nSaved corrected patient-level score file:")
print(" ", out_csv)
print("Saved corrected KM figure:")
print(" ", out_fig)
