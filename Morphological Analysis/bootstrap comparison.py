"""
================================================================================
CLONE MORPHOTYPE RICHNESS — PERMUTATION / BOOTSTRAP ANALYSIS
RICHNESS (ΔR) + HETEROGENEITY INDEX (H, ΔH) • REAL + SIMULATED DATASETS
================================================================================
Purpose:
Quantify whether each clone’s morphotype diversity is higher or lower than expected
by chance, using a permutation-based null model.

What this pipeline does:
- Per-clone permutation test of morphotype richness (R = # unique morphotypes)
- Computes ΔR = observed R − mean null R (plus z-score, two-tailed p-values, CI)
- Computes heterogeneity index:
      H = (R − 1) / (N − 1)   in [0, 1]
  and ΔH = observed H − mean null H
- Produces publication-ready summary figures for:
  • Real dataset alone + comparison across all datasets
  • Homogeneous vs heterogeneous simulated datasets
  • Random vs clonal paired comparisons (richness and H)

CATEGORY POLICY (IMPORTANT):
- Uses "Morphotype" as the preferred category column
- Falls back ONLY to "CellType_With_Complexity" (for simulated/legacy tables)
- "CellType" is NOT accepted (geometry-only; not a morphotype category here)

Input requirements:
- Clone identity column: "CloneKey"
- Category column: "Morphotype" OR "CellType_With_Complexity"
- Optional: "has_swc" boolean to restrict to reconstructed cells

Outputs (saved to OUTPUT_DIR):
- Per-clone results CSV + summary CSV + text report (per dataset)
- Figure 1: Richness (Real + All)
- Figure 2: Richness (Simulated: Homo vs Het)
- Figure 3: Random vs Clonal richness (paired)
- Figure 4: Random vs Clonal richness (paired, colored by clone size)
- Figure 5: Random vs Clonal heterogeneity H (paired)
- Figure 6: Heterogeneity (Real + All)  [NEW]

================================================================================
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
from scipy import stats  # for paired t-tests
from matplotlib.cm import get_cmap, ScalarMappable
from matplotlib.colors import Normalize

warnings.filterwarnings("ignore")

# =============================================================================
# CONFIGURATION
# =============================================================================

# Real dataset can be CSV (your MorphotypeClassifier output) OR Excel.
REAL_PATH = r"Z:\People\Francisco\Code_testing_folder\Output\cell_analysis_detailed_with_complexity.csv"

# Simulated datasets (Excel)
HOMO_PATH = r"Z:\People\Francisco\Code_testing_folder\Input\Homogeneous Simulated Dataset.xlsx"
HET_PATH  = r"Z:\People\Francisco\Code_testing_folder\Input\Heterogeneous Simulated Dataset.xlsx"

OUTPUT_DIR = r"Z:\People\Francisco\Code_testing_folder\Output\Clone_Richness_Permutation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.style.use("seaborn-v0_8-whitegrid")
sns.set_context("paper", font_scale=1.2)

DPI = 150

# Core column names
CLONE_COL = "CloneKey"

# UPDATED TERMINOLOGY:
# - preferred morphotype column name in your per-cell CSV
CATEGORY_COL_PREFERRED = "Morphotype"
# - fallback for older/simulated tables
CATEGORY_COL_FALLBACK  = "CellType_With_Complexity"

# NOTE: "CellType" is intentionally NOT allowed here (geometry-only).
DISALLOWED_GEOMETRY_COL = "CellType"

HAS_SWC_COL = "has_swc"

# Permutation settings
N_PERMUTATIONS = 10000
MIN_CLONE_SIZE = 2
SAMPLE_WITH_REPLACEMENT = True

# Colors for datasets
COLOR_REAL = "grey"
COLOR_HOMO = "red"
COLOR_HET  = "blue"

# Plot readability helpers (overlap reduction)
TICK_FONTSIZE = 7
LEGEND_FONTSIZE = 9


# =============================================================================
# IO HELPERS
# =============================================================================

def load_table(path):
    """Load CSV or Excel transparently."""
    if str(path).lower().endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_excel(path)


def pick_category_column(df):
    """
    Prefer Morphotype; fallback to CellType_With_Complexity.

    IMPORTANT:
    - We intentionally do NOT accept 'CellType' as a fallback, because it encodes
      geometry only (not morphotype categories for this analysis).
    """
    if CATEGORY_COL_PREFERRED in df.columns:
        return CATEGORY_COL_PREFERRED
    if CATEGORY_COL_FALLBACK in df.columns:
        return CATEGORY_COL_FALLBACK
    return None


# =============================================================================
# BASIC HELPERS
# =============================================================================

def summarize_dataset(df, category_col, dataset_name):
    print(f"[INFO] Dataset summary for richness analysis — {dataset_name}")
    print(f"  Total cells (after filtering): {len(df)}")
    print(f"  Total clones: {df[CLONE_COL].nunique()}")
    print(f"  Clones with ≥{MIN_CLONE_SIZE} cells: "
          f"{(df.groupby(CLONE_COL).size() >= MIN_CLONE_SIZE).sum()}")
    print(f"  Unique categories ({category_col}): {df[category_col].nunique()}")
    print()


def generate_null_richness_distribution(global_labels, clone_size, n_perm, sample_with_replacement):
    """
    Null distribution of richness by sampling clone_size labels from global pool.
    """
    N = len(global_labels)
    richness_values = np.empty(n_perm, dtype=float)

    for i in range(n_perm):
        if sample_with_replacement:
            idx = np.random.randint(0, N, size=clone_size)
        else:
            idx = np.random.choice(N, size=clone_size, replace=False)
        sample = global_labels[idx]
        richness_values[i] = len(np.unique(sample))

    return richness_values


def heterogeneity_index(R, N):
    """
    Normalized richness / heterogeneity index in [0, 1]:
        H = (R - 1) / (N - 1)
    Works for scalar or vector R.
    """
    R = np.asarray(R, dtype=float)
    N = np.asarray(N, dtype=float)
    H = np.full_like(R, np.nan, dtype=float)
    mask = (N > 1)
    H[mask] = (R[mask] - 1.0) / (N[mask] - 1.0)
    return H


def analyze_clone_richness(df, category_col, dataset_name):
    print("======================================================================")
    print(f"RUNNING CLONE RICHNESS PERMUTATION ANALYSIS — {dataset_name}")
    print("======================================================================")

    clone_sizes = df.groupby(CLONE_COL).size()
    valid_clones = clone_sizes[clone_sizes >= MIN_CLONE_SIZE].index
    df_valid = df[df[CLONE_COL].isin(valid_clones)].copy()

    global_labels = df_valid[category_col].values
    global_K = df_valid[category_col].nunique()

    print(f"[INFO] Using category column: {category_col}")
    print(f"[INFO] Clone ID column:       {CLONE_COL}")
    print(f"[INFO] Clones analyzed:       {len(valid_clones)}")
    print(f"[INFO] Total cells in analyzed clones: {len(df_valid)}")
    print(f"[INFO] Number of categories in this label: {global_K}")
    print(f"[INFO] Permutations per clone: {N_PERMUTATIONS}")
    print()

    results = []

    for clone in valid_clones:
        clone_df = df_valid[df_valid[CLONE_COL] == clone]
        clone_size = len(clone_df)

        observed_richness = clone_df[category_col].nunique()

        null_richness = generate_null_richness_distribution(
            global_labels,
            clone_size,
            N_PERMUTATIONS,
            SAMPLE_WITH_REPLACEMENT
        )

        mean_null = null_richness.mean()
        sd_null = null_richness.std(ddof=1)

        # heterogeneity index (H)
        observed_H = heterogeneity_index(observed_richness, clone_size)
        null_H = heterogeneity_index(null_richness, clone_size)

        mean_null_H = np.nanmean(null_H)
        sd_null_H = np.nanstd(null_H, ddof=1)

        # p-values for richness
        p_lower = (np.sum(null_richness <= observed_richness) + 1) / (N_PERMUTATIONS + 1)
        p_upper = (np.sum(null_richness >= observed_richness) + 1) / (N_PERMUTATIONS + 1)
        p_two = 2 * min(p_lower, p_upper)
        p_two = min(p_two, 1.0)

        # p-values for H
        pH_lower = (np.sum(null_H <= observed_H) + 1) / (N_PERMUTATIONS + 1)
        pH_upper = (np.sum(null_H >= observed_H) + 1) / (N_PERMUTATIONS + 1)
        pH_two = 2 * min(pH_lower, pH_upper)
        pH_two = min(pH_two, 1.0)

        # z-score for richness
        z_score = (observed_richness - mean_null) / sd_null if sd_null > 0 else np.nan

        # scaled index for richness
        r_max = min(clone_size, global_K)
        scaled_index = ((observed_richness - mean_null) / (r_max - mean_null)) if (r_max > mean_null) else np.nan

        # quantiles
        q_lower, q_median, q_upper = np.percentile(null_richness, [2.5, 50, 97.5])
        qH_lower, qH_median, qH_upper = np.nanpercentile(null_H, [2.5, 50, 97.5])

        results.append({
            "CloneKey": clone,
            "n_cells": clone_size,
            "observed_richness": observed_richness,
            "mean_null_richness": mean_null,
            "sd_null_richness": sd_null,
            "q2.5_null": q_lower,
            "q50_null": q_median,
            "q97.5_null": q_upper,
            "p_lower": p_lower,
            "p_upper": p_upper,
            "p_two_tailed": p_two,
            "z_score": z_score,
            "scaled_index": scaled_index,
            "r_max": r_max,

            "observed_H": float(observed_H),
            "mean_null_H": float(mean_null_H),
            "sd_null_H": float(sd_null_H),
            "q2.5_null_H": float(qH_lower),
            "q50_null_H": float(qH_median),
            "q97.5_null_H": float(qH_upper),
            "pH_lower": float(pH_lower),
            "pH_upper": float(pH_upper),
            "pH_two_tailed": float(pH_two),
        })

    per_clone_df = pd.DataFrame(results)
    per_clone_df["delta_richness"] = per_clone_df["observed_richness"] - per_clone_df["mean_null_richness"]
    per_clone_df["delta_H"] = per_clone_df["observed_H"] - per_clone_df["mean_null_H"]

    print("------------------------------------------------------------------")
    print(f"[DIAGNOSTIC] Number of clones in per_clone_df ({dataset_name}):", len(per_clone_df))
    print("[DIAGNOSTIC] Number of unique (expected, observed) pairs:")
    print(per_clone_df.groupby(["mean_null_richness", "observed_richness"]).size())
    print("------------------------------------------------------------------\n")

    n_clones = len(per_clone_df)
    n_low = (per_clone_df["p_lower"] < 0.05).sum()
    n_high = (per_clone_df["p_upper"] < 0.05).sum()

    print("======================================================================")
    print(f"RESULTS SUMMARY (Richness vs Chance) — {dataset_name}")
    print("======================================================================")
    print(f"Total clones analyzed: {n_clones}")
    print(f"Clones significantly MORE homogeneous (p_lower < 0.05): {n_low}")
    print(f"Clones significantly MORE heterogeneous (p_upper < 0.05): {n_high}")
    print()
    print("Mean observed richness across clones:  {:.2f}".format(per_clone_df["observed_richness"].mean()))
    print("Mean expected richness across clones:  {:.2f}".format(per_clone_df["mean_null_richness"].mean()))
    print("Mean observed - expected difference:   {:.2f}".format(per_clone_df["delta_richness"].mean()))
    print("======================================================================\n")

    global_info = {
        "n_clones": n_clones,
        "n_low": n_low,
        "n_high": n_high,
        "category_col": category_col,
        "dataset_name": dataset_name
    }

    return per_clone_df, global_info


def compute_regression(x, y):
    x = np.asarray(x)
    y = np.asarray(y)
    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]
    if len(x) < 2:
        return np.nan, np.nan, np.nan
    slope, intercept = np.polyfit(x, y, 1)
    r = np.corrcoef(x, y)[0, 1]
    return slope, intercept, r


def add_regression_line(ax, x, y, color, linestyle, label_prefix):
    slope, intercept, r = compute_regression(x, y)
    if not np.isfinite(slope):
        return np.nan, None
    xx = np.linspace(np.nanmin(x), np.nanmax(x), 100)
    yy = slope * xx + intercept
    label = f"{label_prefix} (r={r:.2f})"
    ax.plot(xx, yy, color=color, linestyle=linestyle, linewidth=2.0, label=label)
    return r, label


def save_results_and_report(per_clone_df, global_info, category_col, dataset_tag, dataset_name):
    csv_path = os.path.join(
        OUTPUT_DIR,
        f"clone_richness_permutation_results_{dataset_tag}_{category_col}.csv"
    )
    per_clone_df.to_csv(csv_path, index=False)
    print(f"[INFO] Saved per-clone results to: {csv_path}")

    summary_cols = [
        "CloneKey", "n_cells",
        "observed_richness", "mean_null_richness", "q2.5_null", "q97.5_null", "p_two_tailed",
        "observed_H", "mean_null_H", "q2.5_null_H", "q97.5_null_H", "pH_two_tailed"
    ]
    summary_df = per_clone_df.loc[:, [c for c in summary_cols if c in per_clone_df.columns]].copy()
    summary_path = os.path.join(
        OUTPUT_DIR,
        f"clone_richness_H_summary_{dataset_tag}_{category_col}.csv"
    )
    summary_df.to_csv(summary_path, index=False)
    print(f"[INFO] Saved summary table to: {summary_path}")

    report_path = os.path.join(
        OUTPUT_DIR,
        f"clone_richness_permutation_report_{dataset_tag}_{category_col}.txt"
    )

    with open(report_path, "w", encoding="utf-8") as f:
        f.write("=" * 80 + "\n")
        f.write(f"CLONE MORPHOTYPE RICHNESS — PERMUTATION TEST REPORT ({dataset_name})\n")
        f.write("=" * 80 + "\n\n")
        f.write(f"Category analyzed: {category_col}\n")
        f.write(f"Total clones analyzed: {global_info['n_clones']}\n")
        f.write(f"Clones more homogeneous than expected (p_lower < 0.05): {global_info['n_low']}\n")
        f.write(f"Clones more heterogeneous than expected (p_upper < 0.05): {global_info['n_high']}\n\n")

        mean_obs = per_clone_df["observed_richness"].mean()
        mean_exp = per_clone_df["mean_null_richness"].mean()
        mean_delta = per_clone_df["delta_richness"].mean()

        f.write("SUMMARY ACROSS CLONES:\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Mean observed richness:         {mean_obs:.3f}\n")
        f.write(f"  Mean expected richness:         {mean_exp:.3f}\n")
        f.write(f"  Mean (observed - expected):     {mean_delta:.3f}\n\n")

        f.write("=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

    print(f"[INFO] Saved text report to: {report_path}")


# =============================================================================
# FIGURES
# =============================================================================

def make_figure1(real_entry, all_entries, category_col):
    """
    Figure 1 (Richness):
      A) Real: Observed vs Expected richness + dotted trend lines for all datasets
      B) Real: ΔR histogram
      C) Real: ΔR bar plot
      D) All:  Obs vs Exp (richness) + per-dataset regression
    """
    real_df = real_entry["per_clone_df"]

    fig, axes = plt.subplots(1, 4, figsize=(32, 8), dpi=DPI)
    ax1, ax2, ax3, ax4 = axes

    # Panel A
    x_real = real_df["mean_null_richness"].values
    y_real = real_df["observed_richness"].values
    n_cells_real = real_df["n_cells"].values
    sizes_real = 40 + 10 * np.log2(n_cells_real)

    rng = np.random.default_rng(0)
    x_real_jit = x_real + rng.normal(scale=0.03, size=len(x_real))

    scatter1 = ax1.scatter(
        x_real_jit, y_real,
        s=sizes_real,
        c=n_cells_real,
        cmap="viridis",
        edgecolor="black",
        alpha=0.85
    )
    max_val_real = max(x_real.max(), y_real.max()) + 0.5
    ax1.plot([0, max_val_real], [0, max_val_real], "k--", alpha=0.7, label="y = x")

    for entry in all_entries:
        df = entry["per_clone_df"]
        add_regression_line(
            ax1,
            df["mean_null_richness"].values,
            df["observed_richness"].values,
            color=entry["color"],
            linestyle=":",
            label_prefix=entry["dataset_name"].replace(" Dataset", "")
        )

    ax1.set_xlabel("Expected richness (null model)")
    ax1.set_ylabel("Observed richness")
    ax1.set_title(
        "A) Real Data — Observed vs Expected Richness\n"
        "Points: real data (colored by clone size)\n"
        "Dotted lines: tendency per dataset"
    )
    cbar = plt.colorbar(scatter1, ax=ax1, fraction=0.05, pad=0.04)
    cbar.set_label("Clone size (n cells)")
    ax1.legend(loc="upper left", fontsize=LEGEND_FONTSIZE)

    # Panel B
    delta_real = real_df["delta_richness"]
    sns.histplot(delta_real, bins=10, kde=True, ax=ax2, color="skyblue", edgecolor="black")
    ax2.axvline(0, color="gray", linestyle="--")
    ax2.set_xlabel("ΔR = Observed - Expected richness")
    ax2.set_ylabel("Number of clones")
    ax2.set_title("B) Real Data — Distribution of ΔR")

    # Panel C
    df_bar = real_df.sort_values("delta_richness")
    x_idx = np.arange(len(df_bar))
    ax3.bar(x_idx, df_bar["delta_richness"], color="lightsteelblue", edgecolor="black")
    ax3.axhline(0, color="gray", linestyle="--", linewidth=1)

    xlabels = [f"{ck}\n(n={int(n)})" for ck, n in zip(df_bar["CloneKey"], df_bar["n_cells"])]
    ax3.set_xticks(x_idx)
    ax3.set_xticklabels(xlabels, rotation=90, ha="right", fontsize=TICK_FONTSIZE)
    ax3.set_ylabel("ΔR")
    ax3.set_title("C) Real Data — Clone-specific ΔR\n(sorted by ΔR)")

    # Panel D
    rng = np.random.default_rng(1)
    for entry in all_entries:
        df = entry["per_clone_df"]
        x = df["mean_null_richness"].values
        y = df["observed_richness"].values
        sizes = 40 + 10 * np.log2(df["n_cells"].values)
        x_jit = x + rng.normal(scale=0.03, size=len(x))

        ax4.scatter(x_jit, y, s=sizes, color=entry["color"], edgecolor="black", alpha=0.85)
        add_regression_line(ax4, x, y, color=entry["color"], linestyle="-", label_prefix=entry["dataset_name"])

    x_all = np.concatenate([e["per_clone_df"]["mean_null_richness"].values for e in all_entries])
    y_all = np.concatenate([e["per_clone_df"]["observed_richness"].values for e in all_entries])
    max_val_all = max(x_all.max(), y_all.max()) + 0.5
    ax4.plot([0, max_val_all], [0, max_val_all], "k--", alpha=0.7, label="y = x")

    ax4.set_xlabel("Expected richness (null model)")
    ax4.set_ylabel("Observed richness")
    ax4.set_title("D) All Datasets — Observed vs Expected\nSolid lines: regression per dataset (with r)")
    ax4.legend(loc="upper left", fontsize=LEGEND_FONTSIZE)

    fig.suptitle(
        f"Figure 1 — Clone Morphotype Richness vs Null Model\nCategory: {category_col}",
        fontsize=16, y=0.98
    )

    # overlap protection (Panel C labels)
    fig.subplots_adjust(left=0.05, right=0.98, top=0.85, bottom=0.32, wspace=0.35)

    out_path = os.path.join(OUTPUT_DIR, f"Figure1_clone_richness_real_and_all_{category_col}.svg")
    fig.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 1 to: {out_path}")


def make_figure1_heterogeneity(real_entry, all_entries, category_col):
    """
    Figure 6 (Heterogeneity H): same layout as Figure 1 but with H and ΔH.
    """
    real_df = real_entry["per_clone_df"].copy()

    fig, axes = plt.subplots(1, 4, figsize=(32, 8), dpi=DPI)
    ax1, ax2, ax3, ax4 = axes

    # Panel A
    x_real = real_df["mean_null_H"].values
    y_real = real_df["observed_H"].values
    n_cells_real = real_df["n_cells"].values
    sizes_real = 40 + 10 * np.log2(n_cells_real)

    rng = np.random.default_rng(10)
    x_real_jit = x_real + rng.normal(scale=0.01, size=len(x_real))

    scatter1 = ax1.scatter(
        x_real_jit, y_real,
        s=sizes_real,
        c=n_cells_real,
        cmap="viridis",
        edgecolor="black",
        alpha=0.85
    )
    max_val_real = max(np.nanmax(x_real), np.nanmax(y_real)) + 0.05
    ax1.plot([0, max_val_real], [0, max_val_real], "k--", alpha=0.7, label="y = x")

    for entry in all_entries:
        df = entry["per_clone_df"]
        add_regression_line(
            ax1,
            df["mean_null_H"].values,
            df["observed_H"].values,
            color=entry["color"],
            linestyle=":",
            label_prefix=entry["dataset_name"].replace(" Dataset", "")
        )

    ax1.set_xlabel("Expected heterogeneity H (null model)")
    ax1.set_ylabel("Observed heterogeneity H")
    ax1.set_title(
        "A) Real Data — Observed vs Expected H\n"
        "Points: real data (colored by clone size)\n"
        "Dotted lines: tendency per dataset"
    )
    cbar = plt.colorbar(scatter1, ax=ax1, fraction=0.05, pad=0.04)
    cbar.set_label("Clone size (n cells)")
    ax1.legend(loc="upper left", fontsize=LEGEND_FONTSIZE)

    # Panel B
    delta_H = real_df["delta_H"]
    sns.histplot(delta_H, bins=10, kde=True, ax=ax2, color="skyblue", edgecolor="black")
    ax2.axvline(0, color="gray", linestyle="--")
    ax2.set_xlabel("ΔH = Observed - Expected")
    ax2.set_ylabel("Number of clones")
    ax2.set_title("B) Real Data — Distribution of ΔH")

    # Panel C
    df_bar = real_df.sort_values("delta_H")
    x_idx = np.arange(len(df_bar))
    ax3.bar(x_idx, df_bar["delta_H"], color="lightsteelblue", edgecolor="black")
    ax3.axhline(0, color="gray", linestyle="--", linewidth=1)

    xlabels = [f"{ck}\n(n={int(n)})" for ck, n in zip(df_bar["CloneKey"], df_bar["n_cells"])]
    ax3.set_xticks(x_idx)
    ax3.set_xticklabels(xlabels, rotation=90, ha="right", fontsize=TICK_FONTSIZE)
    ax3.set_ylabel("ΔH")
    ax3.set_title("C) Real Data — Clone-specific ΔH\n(sorted by ΔH)")

    # Panel D
    rng = np.random.default_rng(11)
    for entry in all_entries:
        df = entry["per_clone_df"]
        x = df["mean_null_H"].values
        y = df["observed_H"].values
        sizes = 40 + 10 * np.log2(df["n_cells"].values)
        x_jit = x + rng.normal(scale=0.01, size=len(x))

        ax4.scatter(x_jit, y, s=sizes, color=entry["color"], edgecolor="black", alpha=0.85)
        add_regression_line(ax4, x, y, color=entry["color"], linestyle="-", label_prefix=entry["dataset_name"])

    x_all = np.concatenate([e["per_clone_df"]["mean_null_H"].values for e in all_entries])
    y_all = np.concatenate([e["per_clone_df"]["observed_H"].values for e in all_entries])
    max_val_all = max(np.nanmax(x_all), np.nanmax(y_all)) + 0.05
    ax4.plot([0, max_val_all], [0, max_val_all], "k--", alpha=0.7, label="y = x")

    ax4.set_xlabel("Expected heterogeneity H (null model)")
    ax4.set_ylabel("Observed heterogeneity H")
    ax4.set_title("D) All Datasets — Observed vs Expected H\nSolid lines: regression per dataset (with r)")
    ax4.legend(loc="upper left", fontsize=LEGEND_FONTSIZE)

    fig.suptitle(
        f"Figure 6 — Clone Heterogeneity Index vs Null Model\nCategory: {category_col}",
        fontsize=16, y=0.98
    )

    fig.subplots_adjust(left=0.05, right=0.98, top=0.85, bottom=0.32, wspace=0.35)

    out_path = os.path.join(OUTPUT_DIR, f"Figure6_clone_heterogeneity_real_and_all_{category_col}.svg")
    fig.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 6 to: {out_path}")


def _plot_three_panels_for_dataset(axA, axB, axC, df, dataset_label_prefix, dataset_name):
    """
    Helper: Obs vs Exp richness, ΔR histogram, ΔR bar plot for one dataset.
    """
    x = df["mean_null_richness"].values
    y = df["observed_richness"].values
    n_cells = df["n_cells"].values
    sizes = 40 + 10 * np.log2(n_cells)

    r_val = np.corrcoef(x, y)[0, 1] if len(x) > 1 else np.nan
    rng = np.random.default_rng(0)
    x_jit = x + rng.normal(scale=0.03, size=len(x))

    sc = axA.scatter(x_jit, y, s=sizes, c=n_cells, cmap="viridis",
                     edgecolor="black", alpha=0.85)
    max_val = max(x.max(), y.max()) + 0.5
    axA.plot([0, max_val], [0, max_val], "k--", alpha=0.7, label="y = x")
    axA.set_xlabel("Expected richness (null model)")
    axA.set_ylabel("Observed richness")
    axA.set_title(f"{dataset_label_prefix}) {dataset_name} — Obs vs Exp\n(r = {r_val:.2f})")
    cbar = plt.colorbar(sc, ax=axA, fraction=0.05, pad=0.04)
    cbar.set_label("Clone size (n cells)")
    axA.legend(loc="upper left", fontsize=LEGEND_FONTSIZE)

    delta = df["delta_richness"]
    sns.histplot(delta, bins=10, kde=True, ax=axB, color="skyblue", edgecolor="black")
    axB.axvline(0, color="gray", linestyle="--")
    axB.set_xlabel("ΔR = Observed - Expected richness")
    axB.set_ylabel("Number of clones")
    axB.set_title(f"{chr(ord(dataset_label_prefix)+1)}) {dataset_name} — ΔR Distribution")

    df_bar = df.sort_values("delta_richness")
    x_idx = np.arange(len(df_bar))
    axC.bar(x_idx, df_bar["delta_richness"], color="lightsteelblue", edgecolor="black")
    axC.axhline(0, color="gray", linestyle="--", linewidth=1)

    labels = [f"{ck}\n(n={int(n)})" for ck, n in zip(df_bar["CloneKey"], df_bar["n_cells"])]
    axC.set_xticks(x_idx)
    axC.set_xticklabels(labels, rotation=90, ha="right", fontsize=TICK_FONTSIZE)
    axC.set_ylabel("ΔR")
    axC.set_title(f"{chr(ord(dataset_label_prefix)+2)}) {dataset_name} — Clone-specific ΔR")


def make_figure2(homo_entry, het_entry, category_col):
    if homo_entry is None or het_entry is None:
        print("[WARN] Cannot make Figure 2 — missing homo or het dataset.")
        return

    homo_df = homo_entry["per_clone_df"]
    het_df = het_entry["per_clone_df"]

    fig, axes = plt.subplots(2, 3, figsize=(28, 12), dpi=DPI)
    axes = axes.flatten()

    _plot_three_panels_for_dataset(axes[0], axes[1], axes[2], homo_df, "A", "Homogeneous Simulated")
    _plot_three_panels_for_dataset(axes[3], axes[4], axes[5], het_df,  "D", "Heterogeneous Simulated")

    fig.suptitle(
        f"Figure 2 — Simulated Datasets: Homogeneous vs Heterogeneous\nCategory: {category_col}",
        fontsize=16, y=0.98
    )

    fig.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.28, hspace=0.35, wspace=0.25)

    out_path = os.path.join(OUTPUT_DIR, f"Figure2_clone_richness_simulated_{category_col}.svg")
    fig.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 2 to: {out_path}")


def _paired_panel_richness(ax, df, nice_name, panel_label):
    expected = df["mean_null_richness"].values
    observed = df["observed_richness"].values

    t_stat, p_val = stats.ttest_rel(observed, expected, nan_policy="omit")
    mean_delta = np.nanmean(observed - expected)

    tmp = df.copy()
    tmp["expected"] = expected
    tmp["observed"] = observed
    tmp = tmp.sort_values("expected").reset_index(drop=True)

    x_clonal = 0
    x_random = 1

    for _, row in tmp.iterrows():
        ax.plot([x_clonal, x_random], [row["observed"], row["expected"]],
                color="0.7", alpha=0.7, linewidth=1.2)
        ax.scatter([x_clonal], [row["observed"]], color="tab:red", s=25, zorder=3)
        ax.scatter([x_random], [row["expected"]], color="tab:blue", s=25, zorder=3)

    ax.set_xlim(-0.2, 1.2)
    ax.set_xticks([x_clonal, x_random])
    ax.set_xticklabels(["Clonal", "Random"])
    ax.set_ylabel("Richness (# morphotypes)")
    ax.set_title(f"{panel_label} {nice_name} — Random vs Clonal richness\n"
                 f"t = {t_stat:.2f}, p = {p_val:.3g}, mean Δ = {mean_delta:.2f}")
    ax.grid(axis="y", alpha=0.3)


def make_figure3_intra_vs_inter(all_entries, category_col):
    mapping = {}
    for entry in all_entries:
        if entry["dataset_tag"] == "real":
            mapping["A"] = ("Real Data", entry)
        elif entry["dataset_tag"] == "homo":
            mapping["B"] = ("Homogeneous Simulated", entry)
        elif entry["dataset_tag"] == "het":
            mapping["C"] = ("Heterogeneous Simulated", entry)

    ordered_labels = [k for k in ["A", "B", "C"] if k in mapping]
    if not ordered_labels:
        print("[WARN] No datasets available for Figure 3.")
        return

    fig, axes = plt.subplots(1, len(ordered_labels), figsize=(6 * len(ordered_labels), 7), dpi=DPI)
    if len(ordered_labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, ordered_labels):
        nice_name, entry = mapping[label]
        _paired_panel_richness(ax, entry["per_clone_df"], nice_name, f"{label})")

    fig.suptitle(
        f"Figure 3 — Random vs Clonal Richness (All Datasets)\nCategory: {category_col}",
        fontsize=16, y=0.98
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.12, wspace=0.30)

    out_path = os.path.join(OUTPUT_DIR, f"Figure3_random_vs_clonal_richness_{category_col}.svg")
    fig.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 3 to: {out_path}")


def _paired_panel_richness_colored(ax, df, nice_name, panel_label):
    expected = df["mean_null_richness"].values
    observed = df["observed_richness"].values
    n_cells = df["n_cells"].values

    t_stat, p_val = stats.ttest_rel(observed, expected, nan_policy="omit")
    mean_delta = np.nanmean(observed - expected)

    tmp = df.copy()
    tmp["expected"] = expected
    tmp["observed"] = observed
    tmp = tmp.sort_values("expected").reset_index(drop=True)

    x_random = 0
    x_clonal = 1

    cmap = get_cmap("viridis")
    norm = Normalize(vmin=n_cells.min(), vmax=n_cells.max())

    for _, row in tmp.iterrows():
        color_point = cmap(norm(row["n_cells"]))
        ax.plot([x_random, x_clonal], [row["expected"], row["observed"]],
                color="black", alpha=0.45, linewidth=1.2)
        ax.scatter([x_random], [row["expected"]], color=color_point, s=30, zorder=3)
        ax.scatter([x_clonal], [row["observed"]], color=color_point, s=30, zorder=3)

    ax.set_xlim(-0.2, 1.2)
    ax.set_xticks([x_random, x_clonal])
    ax.set_xticklabels(["Random", "Clonal"])
    ax.set_ylabel("Richness (# morphotypes)")
    ax.set_title(f"{panel_label} {nice_name} — Random vs Clonal\n"
                 f"t = {t_stat:.2f}, p = {p_val:.3g}, mean Δ = {mean_delta:.2f}")
    ax.grid(axis="y", alpha=0.3)

    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.05, pad=0.04)
    cbar.set_label("Clone size (n cells)")


def make_figure4_intra_vs_inter_colored(all_entries, category_col):
    mapping = {}
    for entry in all_entries:
        if entry["dataset_tag"] == "real":
            mapping["A"] = ("Real Data", entry)
        elif entry["dataset_tag"] == "homo":
            mapping["B"] = ("Homogeneous Simulated", entry)
        elif entry["dataset_tag"] == "het":
            mapping["C"] = ("Heterogeneous Simulated", entry)

    ordered_labels = [k for k in ["A", "B", "C"] if k in mapping]
    if not ordered_labels:
        print("[WARN] No datasets available for Figure 4.")
        return

    fig, axes = plt.subplots(1, len(ordered_labels), figsize=(6 * len(ordered_labels), 7), dpi=DPI)
    if len(ordered_labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, ordered_labels):
        nice_name, entry = mapping[label]
        _paired_panel_richness_colored(ax, entry["per_clone_df"], nice_name, f"{label})")

    fig.suptitle(
        f"Figure 4 — Random vs Clonal Richness (Colour-coded by clone size)\nCategory: {category_col}",
        fontsize=16, y=0.98
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.12, wspace=0.30)

    out_path = os.path.join(OUTPUT_DIR, f"Figure4_random_vs_clonal_colored_{category_col}.svg")
    fig.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 4 to: {out_path}")


def _paired_panel_heterogeneity(ax, df, nice_name, panel_label):
    N = df["n_cells"].values
    H_clonal = df["observed_H"].values
    H_random = df["mean_null_H"].values

    t_stat, p_val = stats.ttest_rel(H_clonal, H_random, nan_policy="omit")
    mean_delta = np.nanmean(H_clonal - H_random)

    tmp = df.copy()
    tmp["H_random"] = H_random
    tmp["H_clonal"] = H_clonal
    tmp = tmp.sort_values("H_random").reset_index(drop=True)

    x_random = 0
    x_clonal = 1

    for _, row in tmp.iterrows():
        ax.plot([x_random, x_clonal], [row["H_random"], row["H_clonal"]],
                color="0.7", alpha=0.7, linewidth=1.2)
        ax.scatter([x_random], [row["H_random"]], color="tab:blue", s=25, zorder=3)
        ax.scatter([x_clonal], [row["H_clonal"]], color="tab:red", s=25, zorder=3)

    ax.set_xlim(-0.2, 1.2)
    ax.set_xticks([x_random, x_clonal])
    ax.set_xticklabels(["Random", "Clonal"])
    ax.set_ylabel("Heterogeneity index H (0–1)")
    ax.set_title(f"{panel_label} {nice_name} — Random vs Clonal H\n"
                 f"t = {t_stat:.2f}, p = {p_val:.3g}, mean ΔH = {mean_delta:.2f}")
    ax.grid(axis="y", alpha=0.3)


def make_figure5_heterogeneity(all_entries, category_col):
    mapping = {}
    for entry in all_entries:
        if entry["dataset_tag"] == "real":
            mapping["A"] = ("Real Data", entry)
        elif entry["dataset_tag"] == "homo":
            mapping["B"] = ("Homogeneous Simulated", entry)
        elif entry["dataset_tag"] == "het":
            mapping["C"] = ("Heterogeneous Simulated", entry)

    ordered_labels = [k for k in ["A", "B", "C"] if k in mapping]
    if not ordered_labels:
        print("[WARN] No datasets available for Figure 5.")
        return

    fig, axes = plt.subplots(1, len(ordered_labels), figsize=(6 * len(ordered_labels), 7), dpi=DPI)
    if len(ordered_labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, ordered_labels):
        nice_name, entry = mapping[label]
        _paired_panel_heterogeneity(ax, entry["per_clone_df"], nice_name, f"{label})")

    fig.suptitle(
        f"Figure 5 — Random vs Clonal Heterogeneity Index (All Datasets)\nCategory: {category_col}",
        fontsize=16, y=0.98
    )
    fig.subplots_adjust(left=0.08, right=0.98, top=0.88, bottom=0.12, wspace=0.30)

    out_path = os.path.join(OUTPUT_DIR, f"Figure5_random_vs_clonal_heterogeneity_{category_col}.svg")
    fig.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 5 to: {out_path}")


# =============================================================================
# RUN A SINGLE DATASET (CSV or Excel)
# =============================================================================

def run_single_dataset(path, dataset_name, dataset_tag):
    print("======================================================================")
    print(f"PROCESSING DATASET: {dataset_name}")
    print("======================================================================")

    try:
        print("[INFO] Loading:", path)
        df = load_table(path)
        print("[INFO] Data loaded. Shape:", df.shape)
    except FileNotFoundError:
        print("[ERROR] Could not find file at:", path)
        return None
    except Exception as e:
        print("[ERROR] Problem loading file:", e)
        return None

    # has_swc handling
    if HAS_SWC_COL in df.columns:
        initial_n = len(df)
        df = df[df[HAS_SWC_COL] == True].copy()
        print(f"[INFO] Filtered to reconstructed cells using '{HAS_SWC_COL}': {len(df)} of {initial_n}")
    else:
        print("[INFO] 'has_swc' not found -> assuming all rows are reconstructed (has_swc=True).")
        df[HAS_SWC_COL] = True

    # choose category column (Morphotype preferred; fallback only to CellType_With_Complexity)
    category_col = pick_category_column(df)
    if category_col is None:
        print("[ERROR] Missing morphotype column. Need 'Morphotype' or 'CellType_With_Complexity'.")
        if DISALLOWED_GEOMETRY_COL in df.columns:
            print(f"       NOTE: '{DISALLOWED_GEOMETRY_COL}' is present, but it encodes geometry only and is NOT used here.")
        print("       Available columns:", df.columns.tolist())
        return None

    # required cols
    required_cols = [CLONE_COL, category_col]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print("[ERROR] Missing required columns:", missing)
        print("       Available columns:", df.columns.tolist())
        return None

    df_clean = df.dropna(subset=[CLONE_COL, category_col]).copy()

    summarize_dataset(df_clean, category_col, dataset_name)

    per_clone_df, global_info = analyze_clone_richness(df_clean, category_col, dataset_name)

    # Save per-dataset outputs
    save_results_and_report(per_clone_df, global_info, category_col, dataset_tag, dataset_name)

    return {
        "per_clone_df": per_clone_df,
        "global_info": global_info,
        "dataset_name": dataset_name,
        "dataset_tag": dataset_tag
    }


# =============================================================================
# MAIN
# =============================================================================

def main():
    dataset_configs = [
        {"path": REAL_PATH, "dataset_name": "Real Data", "dataset_tag": "real", "color": COLOR_REAL},
        {"path": HOMO_PATH, "dataset_name": "Homogeneous Simulated Dataset", "dataset_tag": "homo", "color": COLOR_HOMO},
        {"path": HET_PATH,  "dataset_name": "Heterogeneous Simulated Dataset", "dataset_tag": "het",  "color": COLOR_HET},
    ]

    all_entries = []
    for cfg in dataset_configs:
        entry = run_single_dataset(cfg["path"], cfg["dataset_name"], cfg["dataset_tag"])
        if entry is not None:
            entry["color"] = cfg["color"]
            all_entries.append(entry)

    if not all_entries:
        print("[ERROR] No datasets were successfully processed.")
        return

    real_entry = next((e for e in all_entries if e["dataset_tag"] == "real"), None)
    homo_entry = next((e for e in all_entries if e["dataset_tag"] == "homo"), None)
    het_entry  = next((e for e in all_entries if e["dataset_tag"] == "het"), None)

    # Determine label used for outputs (use preferred if real has it; otherwise fallback)
    if real_entry is not None:
        category_used = real_entry["global_info"]["category_col"]
    else:
        category_used = CATEGORY_COL_PREFERRED

    # Figure 1 (richness) + Figure 6 (heterogeneity)
    if real_entry is not None:
        make_figure1(real_entry, all_entries, category_used)
        make_figure1_heterogeneity(real_entry, all_entries, category_used)
    else:
        print("[WARN] Real dataset missing — cannot make Figure 1 / Figure 6.")

    # Figure 2
    make_figure2(homo_entry, het_entry, category_used)

    # Figure 3
    make_figure3_intra_vs_inter(all_entries, category_used)

    # Figure 4
    make_figure4_intra_vs_inter_colored(all_entries, category_used)

    # Figure 5
    make_figure5_heterogeneity(all_entries, category_used)

    print("\n======================================================================")
    print("ANALYSIS COMPLETE FOR ALL DATASETS!")
    print("======================================================================")
    print(f"All results saved to: {OUTPUT_DIR}")
    print("  - clone_richness_permutation_results_*_*.csv")
    print("  - clone_richness_H_summary_*_*.csv")
    print("  - clone_richness_permutation_report_*_*.txt")
    print("  - Figure1_clone_richness_real_and_all_*.svg")
    print("  - Figure2_clone_richness_simulated_*.svg")
    print("  - Figure3_random_vs_clonal_richness_*.svg")
    print("  - Figure4_random_vs_clonal_colored_*.svg")
    print("  - Figure5_random_vs_clonal_heterogeneity_*.svg")
    print("  - Figure6_clone_heterogeneity_real_and_all_*.svg")
    print("======================================================================")


if __name__ == "__main__":
    main()
