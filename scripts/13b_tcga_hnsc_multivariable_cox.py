"""
Corrected TCGA multivariable Cox models for the primary-tumor analysis
---------------------------------------------------------------------
Reads:
  results/06_TCGA_validation/TCGA_threegene_primarytumor_scored.csv

Primary model:
  HPV-negative patients only
  7-year overall survival endpoint (2555 days)
  continuous NFE2L2/MDM2/PPARG composite score
  adjusted for age + AJCC overall stage

Stage handling:
  STAGE IVA / IVB / IVC are collapsed to STAGE IV.
  Stage is modeled categorically with STAGE I as the reference group.

Missing covariates:
  complete-case analysis.

Sensitivity model:
  primary model + radiation therapy (Yes/No).

Run from:
  ~/Desktop/PNI_project/
"""

import pandas as pd
import numpy as np
from lifelines import CoxPHFitter
from lifelines.statistics import proportional_hazard_test

INPUT = "results/06_TCGA_validation/TCGA_threegene_primarytumor_scored.csv"
MAX_DAYS = 2555

AGE = "Diagnosis Age"
STAGE = "Neoplasm Disease Stage American Joint Committee on Cancer Code"
RAD = "Radiation Therapy"

def censor_7yr(df):
    df = df.copy()
    beyond = df["OS.time"] > MAX_DAYS
    df.loc[beyond, "OS.event"] = 0
    df.loc[beyond, "OS.time"] = MAX_DAYS
    return df

def broad_stage(x):
    if pd.isna(x):
        return np.nan
    x = str(x).strip().upper()
    if x == "STAGE I":
        return "I"
    if x == "STAGE II":
        return "II"
    if x == "STAGE III":
        return "III"
    if x.startswith("STAGE IV"):
        return "IV"
    return np.nan

def fit_model(df, include_radiation=False):
    work = df.copy()
    work["stage4"] = work[STAGE].apply(broad_stage)

    needed = ["OS.time", "OS.event", "score", AGE, "stage4"]
    if include_radiation:
        needed.append(RAD)

    work = work[needed].dropna().copy()

    # Apply 7-year administrative censoring AFTER choosing complete cases.
    work = censor_7yr(work)

    # Encode stage categorically, reference = Stage I.
    stage_cat = pd.Categorical(
        work["stage4"],
        categories=["I", "II", "III", "IV"],
        ordered=False
    )
    dummies = pd.get_dummies(stage_cat, prefix="stage", drop_first=True, dtype=int)
    dummies.index = work.index

    model = pd.concat(
        [
            work[["OS.time", "OS.event", "score", AGE]].rename(
                columns={AGE: "age"}
            ),
            dummies
        ],
        axis=1
    )

    if include_radiation:
        rad_map = {"Yes": 1, "No": 0}
        model["radiation_yes"] = work[RAD].map(rad_map)
        model = model.dropna()

    cph = CoxPHFitter()
    cph.fit(model, duration_col="OS.time", event_col="OS.event")

    ph = proportional_hazard_test(cph, model, time_transform="rank")

    return work, model, cph, ph

# ---------------------------------------------------------------------
# Load corrected patient-level primary-tumor cohort
# ---------------------------------------------------------------------
d = pd.read_csv(INPUT)
d = d[d["HPV"] == "HPV-"].copy()

print("=" * 80)
print("TCGA CORRECTED MULTIVARIABLE COX — PRIMARY TUMOR, HPV-")
print("=" * 80)
print(f"Starting HPV- cohort: N={len(d)}")

print("\nStage recoding:")
stage_check = pd.DataFrame({
    "original_stage": d[STAGE],
    "stage4": d[STAGE].apply(broad_stage)
})
print(stage_check["stage4"].value_counts(dropna=False).to_string())

# ---------------------------------------------------------------------
# Primary model: score + age + stage
# ---------------------------------------------------------------------
work1, model1, cph1, ph1 = fit_model(d, include_radiation=False)

print("\n" + "-" * 80)
print("PRIMARY MODEL: score + age + categorical stage")
print("-" * 80)
print(f"Complete-case N={len(model1)}")
print(f"7-year events={int(model1['OS.event'].sum())}")

r = cph1.summary.loc["score"]
print("\nThree-gene score:")
print(f"  adjusted HR={r['exp(coef)']:.3f}")
print(
    f"  95% CI={r['exp(coef) lower 95%']:.3f}–"
    f"{r['exp(coef) upper 95%']:.3f}"
)
print(f"  p={r['p']:.8g}")

print("\nFull primary model:")
show_cols = [
    "coef", "exp(coef)", "exp(coef) lower 95%",
    "exp(coef) upper 95%", "p"
]
print(cph1.summary[show_cols].to_string())

print("\nProportional-hazards test p-values:")
print(ph1.summary["p"].to_string())

# ---------------------------------------------------------------------
# Sensitivity model: + radiation
# ---------------------------------------------------------------------
work2, model2, cph2, ph2 = fit_model(d, include_radiation=True)

print("\n" + "-" * 80)
print("SENSITIVITY MODEL: score + age + categorical stage + radiation")
print("-" * 80)
print(f"Complete-case N={len(model2)}")
print(f"7-year events={int(model2['OS.event'].sum())}")

r2 = cph2.summary.loc["score"]
print("\nThree-gene score:")
print(f"  adjusted HR={r2['exp(coef)']:.3f}")
print(
    f"  95% CI={r2['exp(coef) lower 95%']:.3f}–"
    f"{r2['exp(coef) upper 95%']:.3f}"
)
print(f"  p={r2['p']:.8g}")

print("\nFull sensitivity model:")
print(cph2.summary[show_cols].to_string())

print("\nProportional-hazards test p-values:")
print(ph2.summary["p"].to_string())

print("\n" + "=" * 80)
print("END")
print("=" * 80)
