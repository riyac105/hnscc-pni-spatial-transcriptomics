"""
gse65858_final.py — Clean final version
Independent validation of MDM2/PPARG/NFE2L2 in GSE65858
Run from ~/Desktop/PNI_project/
"""
import os, gzip, numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import warnings; warnings.filterwarnings('ignore')

OUTPUT_DIR = "results/06_TCGA_validation"
SOFT_FILE  = "data/GSE65858/GSE65858_family.soft.gz"
MAX_DAYS   = 1825
TIMEPOINTS = [0, 365, 730, 1095, 1460, 1825]
os.makedirs(OUTPUT_DIR, exist_ok=True)

def censor(df, tc, ec, mx):
    df = df.copy(); b = df[tc] > mx
    df.loc[b, ec] = 0; df.loc[b, tc] = mx
    return df

# ── Parse expression matrix ───────────────────────────────────────────────────
print("Parsing SOFT file...")
samples = {}; current = None; in_table = False; table_rows = []; col_headers = []
plat_rows = []; plat_headers = []; in_plat = False

with gzip.open(SOFT_FILE, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.rstrip("\n")
        if line.startswith("^SAMPLE"):
            if current and table_rows:
                samples[current]["_table"] = pd.DataFrame(table_rows, columns=col_headers)
            current = line.split("=")[1].strip()
            samples[current] = {}
            in_table = False; table_rows = []; col_headers = []
        elif current:
            if line.startswith("!sample_table_begin"):
                in_table = True
            elif line.startswith("!sample_table_end"):
                in_table = False
                if table_rows:
                    samples[current]["_table"] = pd.DataFrame(table_rows, columns=col_headers)
                table_rows = []; col_headers = []
            elif in_table:
                parts = line.split("\t")
                if not col_headers: col_headers = parts
                else: table_rows.append(parts)
        if line.startswith("!platform_table_begin"): in_plat = True
        elif line.startswith("!platform_table_end"): in_plat = False
        elif in_plat:
            parts = line.split("\t")
            if not plat_headers: plat_headers = parts
            else: plat_rows.append(parts)

if current and table_rows:
    samples[current]["_table"] = pd.DataFrame(table_rows, columns=col_headers)
print(f"  Parsed {len(samples)} samples")

# Probe → gene map
probe_map = {}
if plat_rows and plat_headers:
    plat_df = pd.DataFrame(plat_rows, columns=plat_headers)
    sym_col = next((c for c in plat_df.columns if "symbol" in c.lower()), None)
    if sym_col:
        probe_map = plat_df.set_index(plat_df.columns[0])[sym_col].to_dict()
        print(f"  Probe map: {len(probe_map)} probes")

# Expression matrix
expr_dict = {}
for gsm, meta in samples.items():
    if "_table" in meta:
        tbl = meta["_table"]
        if "ID_REF" in tbl.columns and "VALUE" in tbl.columns:
            expr_dict[gsm] = pd.to_numeric(tbl.set_index("ID_REF")["VALUE"], errors="coerce")

expr_df = pd.DataFrame(expr_dict).T
if probe_map:
    expr_df.columns = [probe_map.get(str(c), str(c)) for c in expr_df.columns]
    expr_df = expr_df.T.groupby(level=0).mean().T
print(f"  Expression: {expr_df.shape}")

# ── Clinical data — direct re-parse ──────────────────────────────────────────
print("Parsing clinical data...")
clin_direct = {}
current = None
with gzip.open(SOFT_FILE, "rt", encoding="utf-8", errors="replace") as f:
    for line in f:
        line = line.strip()
        if line.startswith("^SAMPLE"):
            current = line.split("=")[1].strip()
            clin_direct[current] = {}
        elif current and "!Sample_characteristics_ch1" in line:
            val = line.split("= ", 1)[-1].strip()
            if ": " in val:
                k, v = val.split(": ", 1)
                clin_direct[current][k.strip().lower()] = v.strip()

records = []
for gsm, cd in clin_direct.items():
    os_raw  = cd.get("os", None)
    os_days = float(os_raw) if os_raw and os_raw != "NA" else None
    records.append({
        "gsm"     : gsm,
        "OS.days" : os_days,
        "OS.event": cd.get("os_event", None),
        "HPV_raw" : cd.get("hpv16_dna_rna", cd.get("hpv_dna", None)),
    })

clin = pd.DataFrame(records)

def parse_event(x):
    x = str(x).upper().strip()
    if x == "TRUE": return 1
    if x == "FALSE": return 0
    return np.nan

def parse_hpv(x):
    x = str(x).upper().strip()
    if "RNA+" in x: return "HPV+"
    if x == "DNA-" or "RNA-" in x: return "HPV-"
    return np.nan

clin["OS.event"] = clin["OS.event"].apply(parse_event)
clin["HPV"]      = clin["HPV_raw"].apply(parse_hpv)

print(f"  OS valid: {clin['OS.days'].notna().sum()}")
print(f"  HPV: {clin['HPV'].value_counts().to_dict()}")

# ── Score MDM2/PPARG/NFE2L2 ───────────────────────────────────────────────────
print("Scoring signature...")
sig_genes = ["MDM2","PPARG","NFE2L2"]
avail     = [g for g in sig_genes if g in expr_df.columns]
print(f"  Available: {avail}")

clin_idx = clin.set_index("gsm")
common   = expr_df.index.intersection(clin_idx.index)
expr_sub = expr_df.loc[common, avail]
clin_sub = clin_idx.loc[common].copy()

ez = (expr_sub - expr_sub.mean()) / expr_sub.std()
clin_sub["PNI_score"] = ez.mean(axis=1)

data = clin_sub[["OS.days","OS.event","HPV","PNI_score"]].dropna()
med  = data["PNI_score"].median()
data["grp"] = (data["PNI_score"] > med).map({True:"High", False:"Low"})
print(f"  Complete: {len(data)} | HPV: {data['HPV'].value_counts().to_dict()}")

# ── REMARK reporting only (does not alter analysis) ────────────────────────────
print("\n" + "=" * 72)
print("GSE65858 REMARK REPORTING")
print("=" * 72)

print(f"Source cohort: N={len(clin)}")
print(f"Valid OS time: {clin['OS.days'].notna().sum()}")
print(f"Valid OS event: {clin['OS.event'].notna().sum()}")
print(f"Valid HPV status: {clin['HPV'].notna().sum()}")
print(f"Valid three-gene score: {clin_sub['PNI_score'].notna().sum()}")
print(f"Complete analyzable cohort: N={len(data)}")
print(f"Excluded from survival analysis: N={len(clin) - len(data)}")

print("\nMissingness in source cohort:")
print(f"  Missing OS time: {clin['OS.days'].isna().sum()}")
print(f"  Missing OS event: {clin['OS.event'].isna().sum()}")
print(f"  Missing/unclassifiable HPV status: {clin['HPV'].isna().sum()}")
print(f"  Missing three-gene score among matched samples: {clin_sub['PNI_score'].isna().sum()}")

all5 = censor(data, "OS.days", "OS.event", MAX_DAYS)
print(f"\nAll analyzable patients: N={len(all5)}, 5-year events={int(all5['OS.event'].sum())}")

for hpv in ["HPV-", "HPV+"]:
    sub_r = data[data["HPV"] == hpv].copy()
    sub_r = censor(sub_r, "OS.days", "OS.event", MAX_DAYS)
    hi_r = sub_r[sub_r["grp"] == "High"]
    lo_r = sub_r[sub_r["grp"] == "Low"]

    print(f"\n{hpv}:")
    print(f"  N={len(sub_r)}")
    print(f"  5-year events={int(sub_r['OS.event'].sum())}")
    print(f"  Signature-high: N={len(hi_r)}, events={int(hi_r['OS.event'].sum())}")
    print(f"  Signature-low:  N={len(lo_r)}, events={int(lo_r['OS.event'].sum())}")

print("=" * 72)


# ── KM figure ─────────────────────────────────────────────────────────────────
print("Generating figure...")
fig = plt.figure(figsize=(14, 8)); fig.patch.set_facecolor("white")
gs  = gridspec.GridSpec(1, 2, figure=fig,
      left=0.08, right=0.97, top=0.82, bottom=0.04, wspace=0.30)

for ci, (hpv, col_hi, col_lo, subtitle) in enumerate([
    ("HPV-","#c0392b","#27ae60","HPV\u2212 Patients Only"),
    ("HPV+","#e74c3c","#2ecc71","HPV+ Patients \u2014 Negative Control"),
]):
    gs_in   = gridspec.GridSpecFromSubplotSpec(2, 1,
        subplot_spec=gs[ci], height_ratios=[4.5, 1], hspace=0.25)
    ax_km   = fig.add_subplot(gs_in[0])
    ax_risk = fig.add_subplot(gs_in[1])

    sub = data[data["HPV"]==hpv].copy()
    sub = censor(sub, "OS.days", "OS.event", MAX_DAYS)
    hi  = sub[sub["grp"]=="High"]; lo = sub[sub["grp"]=="Low"]

    s = sub.rename(columns={"OS.days":"Days","OS.event":"Event"})
    sh = s[s["grp"]=="High"]; sl = s[s["grp"]=="Low"]

    KaplanMeierFitter().fit(sh["Days"], sh["Event"],
        label=f"High PNI (n={len(sh)})").plot_survival_function(
        ax=ax_km, color=col_hi, ci_show=True, ci_alpha=0.12, linewidth=2.2)
    KaplanMeierFitter().fit(sl["Days"], sl["Event"],
        label=f"Low PNI  (n={len(sl)})").plot_survival_function(
        ax=ax_km, color=col_lo, ci_show=True, ci_alpha=0.12, linewidth=2.2)
    ax_km.set_xlabel("Time (days)", fontsize=10)

    r    = logrank_test(hi["OS.days"], lo["OS.days"],
                        event_observed_A=hi["OS.event"],
                        event_observed_B=lo["OS.event"])
    p    = r.p_value
    pstr = f"p = {p:.4f}" if p >= 0.0001 else "p < 0.0001"
    sig  = "\u2731" if p < 0.05 else "n.s."
    ax_km.text(0.05, 0.08, f"Log-rank {pstr}   {sig}",
        transform=ax_km.transAxes, fontsize=11, fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#FFFBE6",
                  edgecolor="#C8A800", linewidth=1.3, alpha=0.97))
    ax_km.set_title(f"{subtitle}\n(n={len(sub)})", fontsize=12, fontweight="bold", pad=8)
    ax_km.set_ylabel("5-Year Survival Probability", fontsize=10)
    ax_km.set_ylim(0, 1.05); ax_km.set_xlim(-30, MAX_DAYS+50)
    ax_km.set_xticks(TIMEPOINTS)
    ax_km.spines[["top","right"]].set_visible(False)
    ax_km.grid(axis="y", alpha=0.2, linestyle="--")
    ax_km.legend(fontsize=9.5, loc="upper right")

    # Risk table
    hi_n = [int((hi["OS.days"]>=t).sum()) for t in TIMEPOINTS]
    lo_n = [int((lo["OS.days"]>=t).sum()) for t in TIMEPOINTS]
    xlim = ax_km.get_xlim(); xr = xlim[1] - xlim[0]
    ax_risk.set_xlim(*xlim); ax_risk.set_ylim(0, 1); ax_risk.axis("off")
    for lbl, y, col, vals in [("High",0.75,col_hi,hi_n),("Low",0.25,col_lo,lo_n)]:
        ax_risk.text(-0.02, y, lbl, color=col, fontsize=9, fontweight="bold",
            ha="right", va="center", transform=ax_risk.transAxes)
        for t, v in zip(TIMEPOINTS, vals):
            xpos = (t - xlim[0]) / xr
            ax_risk.text(xpos, y, str(v), color=col, fontsize=9,
                fontweight="bold", ha="center", va="center",
                transform=ax_risk.transAxes)
    ax_risk.axhline(0.5, color="#DDDDDD", linewidth=0.8)
    print(f"  {hpv}: {pstr} {sig} | hi={hi_n} lo={lo_n}")

fig.text(0.5, 0.97,
    "GSE65858 Independent Validation \u2014 MDM2/PPARG/NFE2L2 PNI Signature",
    ha="center", va="top", fontsize=13, fontweight="bold")
fig.text(0.5, 0.91,
    "University Hospital Leipzig (n=270 HNSCC)  |  Illumina microarray  |  "
    "5-year OS  |  Independent of TCGA gene selection",
    ha="center", va="top", fontsize=9.5, color="#444444")

out = f"{OUTPUT_DIR}/fig_KM_GSE65858_validation.png"
plt.savefig(out, dpi=150, bbox_inches="tight", facecolor="white")
plt.close()
print(f"\nSaved \u2192 {out}")
