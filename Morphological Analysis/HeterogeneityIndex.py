# ==========================================================
# CLONE HETEROGENEITY INDEX — CLEAN, TALL FIGURE (+ CSV OUTPUT)
# ==========================================================
# - Reads: Morphotype_Integrated_Full_FINAL_with_sholl.xlsx
# - Uses only rows with has_swc == TRUE (reconstructed neurons)
# - Morphotypes = unique values in "CellType_With_Complexity" per CloneKey
# - Heterogeneity index per clone:
#       M = # of distinct morphotypes in the clone
#       N = # of reconstructed neurons in the clone
#       if N <= 1 or M <= 1: HI = 0
#       else: HI = (M - 1) / (N - 1)
#
# - Outputs:
#   1) Console table with [CloneKey, N, M, HI]
#   2) One figure with 3 panels:
#        (A) Bar: clone vs heterogeneity index (with n on top)
#        (B) Histogram of indices (mean & median)
#        (C) Scatter: clone size vs heterogeneity index
#   3) Saved as SVG and JPEG in the same folder as the Excel.
#   4) NEW: CSV with ALL reconstructed rows + per-clone columns:
#        - N_morphotypes_in_clone
#        - HeterogeneityIndex_in_clone
#        - N_reconstructed_in_clone
# ==========================================================

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

########################################################################################################################################################################################################
# User input
EXCEL_PATH = r"C:\Users\fdasilve\Morphological Analysis\Morphotype_Integrated_Full_FINAL_with_sholl.xlsx" # Path of Cell_morphotypes xlsx file

CLONE_COL = "CloneKey"
MORPHO_COL = "CellType_With_Complexity"   # or "CellType"
HAS_SWC_COL = "has_swc"

OUT_DIR = "C:/Users/fdasilve\Morphological Analysis\Clone_Richness_Permutation"
OUT_BASENAME = "clone_heterogeneity_index_"

DPI = 300  # high resolution

# Load Data
print("[INFO] Loading Excel:", EXCEL_PATH)
df = pd.read_excel(EXCEL_PATH)

expected_cols = {CLONE_COL, MORPHO_COL, HAS_SWC_COL}
missing = expected_cols - set(df.columns)
if missing:
    raise ValueError(f"Missing expected columns in Excel: {missing}")

# ---------- FILTER TO RECONSTRUCTED NEURONS ----------
has_swc_mask = df[HAS_SWC_COL].astype(str).str.upper().isin(["TRUE", "1", "YES"])
df_recon = df[has_swc_mask].copy()

print(f"[INFO] Total rows in table       : {len(df)}")
print(f"[INFO] Reconstructed (has_swc)   : {len(df_recon)}")

# keep only rows that can be assigned to a clone + morphotype
df_recon = df_recon.dropna(subset=[CLONE_COL, MORPHO_COL]).copy()

# ---------- COMPUTE HETEROGENEITY PER CLONE ----------
records = []
for clone_key, sub in df_recon.groupby(CLONE_COL):
    N = len(sub)                              # reconstructed neurons
    M = sub[MORPHO_COL].nunique(dropna=True)  # distinct morphotypes

    if N <= 1 or M <= 1:
        hi = 0.0
    else:
        hi = (M - 1.0) / (N - 1.0)

    records.append({
        CLONE_COL: clone_key,
        "N_reconstructed": N,
        "N_morphotypes": M,
        "HeterogeneityIndex": hi
    })

df_hi = pd.DataFrame(records)
if df_hi.empty:
    raise RuntimeError("No reconstructed cells found after filtering; cannot plot.")

# Sort by heterogeneity index (descending)
df_hi = df_hi.sort_values(by="HeterogeneityIndex", ascending=False).reset_index(drop=True)

print("\n[SUMMARY] Clone heterogeneity:")
print(df_hi[[CLONE_COL, "N_reconstructed", "N_morphotypes", "HeterogeneityIndex"]].to_string(index=False))

df_hi.to_csv(
    os.path.join(OUT_DIR, OUT_BASENAME + "per_clone_summary.csv"),
    index=False
)


# ---------- NEW: MERGE BACK TO ROW-LEVEL TABLE & SAVE CSV ----------
# This produces a CSV with ALL ORIGINAL COLUMNS (for reconstructed rows)
# plus the per-clone summary columns repeated on each row of that clone.
df_recon_with_hi = df_recon.merge(
    df_hi[[CLONE_COL, "N_reconstructed", "N_morphotypes", "HeterogeneityIndex"]],
    on=CLONE_COL,
    how="left"
)

# Optional: clearer names for the added columns (keep or remove as you prefer)
df_recon_with_hi = df_recon_with_hi.rename(columns={
    "N_reconstructed": "N_reconstructed_in_clone",
    "N_morphotypes": "N_morphotypes_in_clone",
    "HeterogeneityIndex": "HeterogeneityIndex_in_clone"
})

csv_path = os.path.join(OUT_DIR, OUT_BASENAME + "reconstructed_rows_with_HI.csv")
df_recon_with_hi.to_csv(csv_path, index=False)
print("\n[INFO] Saved reconstructed rows + per-clone HI CSV:")
print("  CSV :", csv_path)

# ---------- PLOTTING STYLE HELPERS ----------
plt.style.use("default")
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial"],
    "axes.linewidth": 1.0,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "xtick.labelsize": 8,
    "ytick.labelsize": 9,
})

BAR_COLOR = "#4C72B0"
HIST_COLOR = "#55A868"
SCATTER_COLOR = "#C44E52"

def clean_axis(ax):
    """Remove top/right spines and add light horizontal grid."""
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="y", linestyle="--", linewidth=0.5, alpha=0.4)

# ---------- DATA ARRAYS ----------
clone_labels = df_hi[CLONE_COL].astype(str).values
hi_values = df_hi["HeterogeneityIndex"].values
clone_sizes = df_hi["N_reconstructed"].values

# ---------- FIGURE LAYOUT (taller so nothing is cut) ----------
fig = plt.figure(figsize=(18, 8), dpi=DPI)
gs = fig.add_gridspec(1, 3, width_ratios=[2.4, 1.1, 1.3])
ax_bar = fig.add_subplot(gs[0, 0])
ax_hist = fig.add_subplot(gs[0, 1])
ax_scatter = fig.add_subplot(gs[0, 2])

# --- (A) BAR PLOT: CLONE vs HETEROGENEITY ---
x_pos = np.arange(len(df_hi))
ax_bar.bar(
    x_pos, hi_values,
    color=BAR_COLOR, edgecolor="black", linewidth=0.6, alpha=0.9
)

# annotate N on top of each bar
for x, h, n_cells in zip(x_pos, hi_values, clone_sizes):
    ax_bar.text(x, h + 0.03, f"n={n_cells}", ha="center", va="bottom", fontsize=7)

ax_bar.set_xticks(x_pos)
ax_bar.set_xticklabels(clone_labels, rotation=90, ha="center")
ax_bar.set_ylabel("Heterogeneity index")
ax_bar.set_xlabel("Clone")
ax_bar.set_title("Heterogeneity per clone", pad=12)
ax_bar.set_ylim(0, 1.1)
clean_axis(ax_bar)

# --- (B) HISTOGRAM: DISTRIBUTION OF INDICES ---
bins = np.linspace(0, 1, 11)
ax_hist.hist(hi_values, bins=bins, color=HIST_COLOR, edgecolor="black", alpha=0.85)

mean_hi = float(np.mean(hi_values))
median_hi = float(np.median(hi_values))
ax_hist.axvline(mean_hi, linestyle="--", linewidth=1.2, label=f"Mean = {mean_hi:.2f}")
ax_hist.axvline(median_hi, linestyle=":", linewidth=1.2, label=f"Median = {median_hi:.2f}")

ax_hist.set_xlabel("Heterogeneity index")
ax_hist.set_ylabel("Number of clones")
ax_hist.set_title("Distribution of heterogeneity indices", pad=12)
ax_hist.legend(frameon=False, fontsize=8)
clean_axis(ax_hist)

# --- (C) SCATTER: CLONE SIZE vs HETEROGENEITY ---
ax_scatter.scatter(
    clone_sizes, hi_values,
    s=45, color=SCATTER_COLOR, edgecolor="black", alpha=0.85
)

for _, row in df_hi.iterrows():
    ax_scatter.text(
        row["N_reconstructed"] + 0.05,
        row["HeterogeneityIndex"] + 0.01,
        str(row[CLONE_COL]),
        fontsize=7, ha="left", va="bottom"
    )

ax_scatter.set_xlabel("Number of reconstructed neurons per clone")
ax_scatter.set_ylabel("Heterogeneity index")
ax_scatter.set_title("Clone size vs heterogeneity", pad=12)
clean_axis(ax_scatter)

# ---- GLOBAL TITLE (with extra top space so it never gets cut) ----
fig.suptitle("Clone heterogeneity based on morphotypes", fontsize=16, y=1.05)

# Tight layout with reserved space for the suptitle
fig.tight_layout(rect=[0, 0, 1, 0.97])

# ---------- SAVE FIGURE ----------
svg_path = os.path.join(OUT_DIR, OUT_BASENAME + ".svg")
jpg_path = os.path.join(OUT_DIR, OUT_BASENAME + ".jpg")

fig.savefig(svg_path, format="svg")
fig.savefig(jpg_path, format="jpg", dpi=DPI)

print("\n[INFO] Saved figure as:")
print("  SVG :", svg_path)
print("  JPEG:", jpg_path)

plt.close(fig)
