"""
================================================================================
CLONE MORPHOTYPE RICHNESS — PERMUTATION / BOOTSTRAP ANALYSIS
Multi-dataset version with 3 main figures + additional variants

  Figure 1 (Real + All Datasets)
    A) Real data: Observed vs Expected
       - Points: real dataset only, colored by clone size
       - Overlaid: regression / tendency lines for all 3 datasets
         (grey = real, red = homo, blue = hetero; dotted)
    B) Real data: Histogram of ΔR = observed - expected
    C) Real data: Clone-specific ΔR bar plot
    D) All datasets: Scatter of Obs vs Exp for all 3 datasets
       + separate regression line and r for EACH dataset (in legend)

  Figure 2 (Simulated datasets)
    A–C) Homogeneous simulated dataset: same 3 panels as real data (but for homo)
    D–F) Heterogeneous simulated dataset: same 3 panels (but for hetero)

  Figure 3 (Intra vs inter comparison for ALL three datasets, richness)
    A) Paired Random vs Clonal richness — Real dataset
    B) Paired Random vs Clonal richness — Homogeneous simulated
    C) Paired Random vs Clonal richness — Heterogeneous simulated

  Figure 4 (Same as Figure 3, but coloured by clone size)
    A–C) As above, but datapoints colour-coded by clone size (with colorbar)
         and LINES IN BLACK

  Figure 5 (Heterogeneity index instead of raw richness)
    A–C) Paired Random vs Clonal H (per-clone heterogeneity index), all datasets
         H = (R - 1) / (N - 1)  (normalized richness 0–1)

  Supplementary figure (for ΔR) — NOT CALLED FROM main()
    S1) Real data only:
        Left: ΔR vs clone size scatter
        Right: Histogram of ΔR

Datasets must have:
  - Clone identity column: "CloneKey"
  - Morphotype category column: "CellType_With_Complexity"
  - Optional: "has_swc" boolean, to restrict to reconstructed cells
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

########################################################################################################################################################################################################
# User input
REAL_EXCEL_PATH = r"Z:\People\Francisco\Morphological Analysis\Data Tables\Morphotype_Integrated_Full_FINAL_with_sholl.xlsx" # Path of Cell_morphotypes xlsx file
HOMO_EXCEL_PATH = r"Z:\People\Francisco\Morphological Analysis\Data Tables\Homogeneous Simulated Dataset.xlsx" # Path of Homogeneous Simulated Dataset file
HET_EXCEL_PATH  = r"Z:\People\Francisco\Morphological Analysis\Data Tables\Heterogeneous Simulated Dataset.xlsx" # Path of Heterogeneous Simulated Dataset file

OUTPUT_DIR = r"Z:\People\Francisco\Morphological Analysis\Clone_Richness_Permutation"
os.makedirs(OUTPUT_DIR, exist_ok=True)

plt.style.use("seaborn-v0_8-whitegrid")
sns.set_context("paper", font_scale=1.2)
DPI = 150

# Core column names
CATEGORY_COL = "CellType_With_Complexity"  # morphotype category
CLONE_COL = "CloneKey"                     # clone identity

# Permutation settings
N_PERMUTATIONS = 10000          # random pseudo-clones per real clone
MIN_CLONE_SIZE = 2              # only analyze clones with at least this many cells
SAMPLE_WITH_REPLACEMENT = True  # True = bootstrap-style sampling

# Colors for datasets
COLOR_REAL = "grey"
COLOR_HOMO = "red"
COLOR_HET  = "blue"

########################################################################################################################################################################################################
# Helpers

def summarize_dataset(df, category_col, dataset_name):
    """Print basic dataset summary."""
    print(f"[INFO] Dataset summary for richness analysis — {dataset_name}")
    print(f"  Total cells (after filtering): {len(df)}")
    print(f"  Total clones: {df[CLONE_COL].nunique()}")
    print(f"  Clones with ≥{MIN_CLONE_SIZE} cells: "
          f"{(df.groupby(CLONE_COL).size() >= MIN_CLONE_SIZE).sum()}")
    print(f"  Unique categories ({category_col}): {df[category_col].nunique()}")
    print()


def generate_null_richness_distribution(global_labels, clone_size,
                                        n_perm, sample_with_replacement):
    """
    Generate null distribution of richness (number of unique labels)
    by randomly sampling 'clone_size' labels from the global pool.
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


def analyze_clone_richness(df, category_col, dataset_name):
    """
    Core analysis:
      - For each clone with >= MIN_CLONE_SIZE:
          * compute observed richness
          * compute null richness distribution
          * derive p-values, z-score, etc.

    Returns:
      per_clone_df (DataFrame with all stats)
      global_info  (dict with summary info)
    """
    print("======================================================================")
    print(f"RUNNING CLONE RICHNESS PERMUTATION ANALYSIS — {dataset_name}")
    print("======================================================================")

    # Keep only clones with at least MIN_CLONE_SIZE cells
    clone_sizes = df.groupby(CLONE_COL).size()
    valid_clones = clone_sizes[clone_sizes >= MIN_CLONE_SIZE].index

    df_valid = df[df[CLONE_COL].isin(valid_clones)].copy()

    # Global pool of labels
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

    # Loop over clones
    for clone in valid_clones:
        clone_df = df_valid[df_valid[CLONE_COL] == clone]
        clone_size = len(clone_df)

        # Observed richness
        observed_richness = clone_df[category_col].nunique()

        # Null distribution
        null_richness = generate_null_richness_distribution(
            global_labels,
            clone_size,
            N_PERMUTATIONS,
            SAMPLE_WITH_REPLACEMENT
        )

        mean_null = null_richness.mean()
        sd_null = null_richness.std(ddof=1)
        
        # --- Heterogeneity index (normalized richness) ---
        observed_H = heterogeneity_index(observed_richness, clone_size)  # scalar
        null_H = heterogeneity_index(null_richness, clone_size)          # vector

        mean_null_H = np.nanmean(null_H)
        sd_null_H = np.nanstd(null_H, ddof=1)

        # Clone-wise Monte-Carlo p-values for H (with +1 correction)
        pH_lower = (np.sum(null_H <= observed_H) + 1) / (N_PERMUTATIONS + 1)
        pH_upper = (np.sum(null_H >= observed_H) + 1) / (N_PERMUTATIONS + 1)
        pH_two = 2 * min(pH_lower, pH_upper)
        pH_two = min(pH_two, 1.0)

        # 95% CI for H under the null
        qH_lower, qH_median, qH_upper = np.nanpercentile(null_H, [2.5, 50, 97.5])

        # p-values (lower tail: more homogeneous)
        p_lower = (np.sum(null_richness <= observed_richness) + 1) / (N_PERMUTATIONS + 1)
        p_upper = (np.sum(null_richness >= observed_richness) + 1) / (N_PERMUTATIONS + 1)
        p_two = 2 * min(p_lower, p_upper)
        p_two = min(p_two, 1.0)

        # z-score
        if sd_null > 0:
            z_score = (observed_richness - mean_null) / sd_null
        else:
            z_score = np.nan

        # Maximum possible richness for this clone
        r_max = min(clone_size, global_K)
        if r_max > mean_null:
            scaled_index = (observed_richness - mean_null) / (r_max - mean_null)
        else:
            scaled_index = np.nan

        # quantiles
        q_lower, q_median, q_upper = np.percentile(null_richness, [2.5, 50, 97.5])

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

    # Add richness deviation column ΔR = observed - expected
    per_clone_df["delta_richness"] = (
        per_clone_df["observed_richness"] - per_clone_df["mean_null_richness"]
    )

    # Diagnostics about dot overlap
    print("------------------------------------------------------------------")
    print(f"[DIAGNOSTIC] Number of clones in per_clone_df ({dataset_name}):", len(per_clone_df))
    print("[DIAGNOSTIC] Number of unique (expected, observed) pairs:")
    print(per_clone_df.groupby(["mean_null_richness", "observed_richness"]).size())
    print("------------------------------------------------------------------\n")

    # Global summary
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
    print("Mean observed richness across clones:  {:.2f}".format(
        per_clone_df["observed_richness"].mean()))
    print("Mean expected richness across clones:  {:.2f}".format(
        per_clone_df["mean_null_richness"].mean()))
    print("Mean observed - expected difference:   {:.2f}".format(
        per_clone_df["delta_richness"].mean()))
    print("======================================================================\n")

    global_info = {
        "n_clones": n_clones,
        "n_low": n_low,
        "n_high": n_high,
        "category_col": category_col,
        "dataset_name": dataset_name
    }

    return per_clone_df, global_info


def save_results_and_report(per_clone_df, global_info, category_col,
                            dataset_tag, dataset_name):
    """Save per-clone CSV and a short report for one dataset."""
    csv_path = os.path.join(
        OUTPUT_DIR,
        f"clone_richness_permutation_results_{dataset_tag}_{category_col}.csv"
    )
    per_clone_df.to_csv(csv_path, index=False)
    print(f"[INFO] Saved per-clone results to: {csv_path}")

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
        f.write(f"Clones more homogeneous than expected (p_lower < 0.05): "
                f"{global_info['n_low']}\n")
        f.write(f"Clones more heterogeneous than expected (p_upper < 0.05): "
                f"{global_info['n_high']}\n\n")

        mean_obs = per_clone_df["observed_richness"].mean()
        mean_exp = per_clone_df["mean_null_richness"].mean()
        mean_delta = per_clone_df["delta_richness"].mean()

        f.write("SUMMARY ACROSS CLONES:\n")
        f.write("-" * 40 + "\n")
        f.write(f"  Mean observed richness:         {mean_obs:.3f}\n")
        f.write(f"  Mean expected richness:         {mean_exp:.3f}\n")
        f.write(f"  Mean (observed - expected):     {mean_delta:.3f}\n\n")

        f.write("INTERPRETATION:\n")
        f.write("-" * 40 + "\n")
        f.write(
            "  If clones were strongly lineage-restricted (few morphotypes),\n"
            "  we would see many clones with significantly lower richness than\n"
            "  expected by chance (p_lower < 0.05). Conversely, strong multi-lineage\n"
            "  behavior would produce clones with consistently higher richness (p_upper < 0.05).\n\n"
            "  Here, the observed richness is compared against a null model in which\n"
            "  morphotypes are randomly assigned to cells according to the global\n"
            "  frequency distribution in this dataset.\n"
        )

        f.write("\nDETAILS PER CLONE (first few rows):\n")
        f.write("-" * 40 + "\n")
        f.write(per_clone_df.head().to_string(index=False))
        f.write("\n\n")

        f.write("=" * 80 + "\n")
        f.write("END OF REPORT\n")
        f.write("=" * 80 + "\n")

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
    print(f"[INFO] Saved text report to: {report_path}")


# Helpers for regression lines and het index
def compute_regression(x, y):
    """Return slope, intercept, Pearson r for given x,y arrays (no NaNs)."""
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
    """
    Compute and plot regression line on given axes.
    Returns the Pearson r and line label (with r embedded).
    """
    slope, intercept, r = compute_regression(x, y)
    if not np.isfinite(slope):
        return np.nan, None

    x_min = x.min()
    x_max = x.max()
    xx = np.linspace(x_min, x_max, 100)
    yy = slope * xx + intercept

    label = f"{label_prefix} (r={r:.2f})"
    ax.plot(xx, yy, color=color, linestyle=linestyle, linewidth=2.0, label=label)
    return r, label


def heterogeneity_index(R, N):
    """
    Normalized richness / heterogeneity index in [0, 1]:

        H = (R - 1) / (N - 1)

    where R = # morphotypes (richness), N = clone size (# reconstructed neurons).
    Works for float R as well (e.g. expected richness).
    """
    R = np.asarray(R, dtype=float)
    N = np.asarray(N, dtype=float)

    H = np.full_like(R, np.nan, dtype=float)
    # N >= 2 because we enforce MIN_CLONE_SIZE = 2, but keep mask just in case
    mask = (N > 1)
    H[mask] = (R[mask] - 1.0) / (N[mask] - 1.0)
    return H


# FIGURE 1 — REAL DATA + ALL DATASETS

def make_figure1(real_entry, all_entries, category_col):
    """
    Figure 1:
      A) Real data: Obs vs Exp (clone size color) + dotted regression lines
         for all datasets (real/homo/het), color-coded.
      B) Real data: ΔR histogram
      C) Real data: ΔR bar plot
      D) All datasets: Obs vs Exp with points + separate regression line and r
         for EACH dataset in the legend.
    """
    real_df = real_entry["per_clone_df"]

    fig, axes = plt.subplots(1, 4, figsize=(28, 6), dpi=DPI)
    ax1, ax2, ax3, ax4 = axes

    # ---------- Panel A: Real data points + dotted regression lines ----------
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
    ax1.plot([0, max_val_real], [0, max_val_real],
             "k--", alpha=0.7,
             label="y = x (observed = expected)")

    for entry in all_entries:
        df = entry["per_clone_df"]
        x = df["mean_null_richness"].values
        y = df["observed_richness"].values
        color = entry["color"]
        ds_name_short = entry["dataset_name"].replace(" Dataset", "")
        r_val, _ = add_regression_line(
            ax1, x, y,
            color=color,
            linestyle=":",
            label_prefix=ds_name_short
        )
        entry["regression_r"] = r_val

    ax1.set_xlabel("Expected richness (null model)")
    ax1.set_ylabel("Observed richness")
    ax1.set_title(
        "A) Real Data — Observed vs Expected Richness\n"
        "Points: real data (colored by clone size)\n"
        "Dotted lines: tendency per dataset"
    )
    cbar = plt.colorbar(scatter1, ax=ax1)
    cbar.set_label("Clone size (n cells)")
    ax1.legend(loc="upper left", fontsize=9)

    # ---------- Panel B: ΔR histogram (real) ----------
    delta_real = real_df["delta_richness"]
    sns.histplot(delta_real, bins=10, kde=True, ax=ax2,
                 color="skyblue", edgecolor="black")
    ax2.axvline(0, color="gray", linestyle="--")
    ax2.set_xlabel("Observed - Expected richness (ΔR)")
    ax2.set_ylabel("Number of clones")
    ax2.set_title("B) Real Data — Distribution of ΔR\n(ΔR = observed - expected)")

    # ---------- Panel C: Clone-specific ΔR bar plot (real) ----------
    df_bar = real_df.sort_values("delta_richness")
    x_idx = np.arange(len(df_bar))
    ax3.bar(x_idx, df_bar["delta_richness"],
            color="lightsteelblue", edgecolor="black")
    ax3.axhline(0, color="gray", linestyle="--", linewidth=1)

    xlabels = [
        f"{ck}\n(n={int(n)})"
        for ck, n in zip(df_bar["CloneKey"], df_bar["n_cells"])
    ]
    ax3.set_xticks(x_idx)
    ax3.set_xticklabels(xlabels, rotation=90, ha="right")
    ax3.set_ylabel("Observed - Expected richness (ΔR)")
    ax3.set_title("C) Real Data — Clone-specific ΔR\n(sorted by ΔR)")

    # ---------- Panel D: All datasets Obs vs Exp + per-dataset regression ----------
    rng = np.random.default_rng(1)
    for entry in all_entries:
        df = entry["per_clone_df"]
        color = entry["color"]
        name = entry["dataset_name"]
        x = df["mean_null_richness"].values
        y = df["observed_richness"].values
        sizes = 40 + 10 * np.log2(df["n_cells"].values)

        x_jit = x + rng.normal(scale=0.03, size=len(x))

        ax4.scatter(
            x_jit, y,
            s=sizes,
            color=color,
            edgecolor="black",
            alpha=0.85
        )

        r_val, _ = add_regression_line(
            ax4, x, y,
            color=color,
            linestyle="-",
            label_prefix=name
        )
        entry["regression_r"] = r_val

    x_all = np.concatenate([e["per_clone_df"]["mean_null_richness"].values for e in all_entries])
    y_all = np.concatenate([e["per_clone_df"]["observed_richness"].values for e in all_entries])
    max_val_all = max(x_all.max(), y_all.max()) + 0.5
    ax4.plot([0, max_val_all], [0, max_val_all],
             "k--", alpha=0.7,
             label="y = x")

    ax4.set_xlabel("Expected richness (null model)")
    ax4.set_ylabel("Observed richness")
    ax4.set_title(
        "D) All Datasets — Observed vs Expected\n"
        "Solid lines: regression per dataset (with r)"
    )
    ax4.legend(loc="upper left", fontsize=9)

    plt.suptitle(
        f"Figure 1 — Clone Morphotype Richness vs Null Model\nCategory: {category_col}",
        fontsize=16, y=1.04
    )
    plt.tight_layout()
    out_path = os.path.join(
        OUTPUT_DIR,
        f"Figure1_clone_richness_real_and_all_{category_col}.svg"
    )
    plt.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 1 to: {out_path}")


# FIGURE 2 — SIMULATED DATASETS ONLY

def _plot_three_panels_for_dataset(axA, axB, axC, df, dataset_label_prefix, dataset_name):
    """
    Helper: draw Obs vs Exp, ΔR histogram, and ΔR bar plot on given axes.
    Obs vs Exp panel uses clone-size color, *without* multi-dataset lines.
    """
    # Panel A: Obs vs Exp
    x = df["mean_null_richness"].values
    y = df["observed_richness"].values
    n_cells = df["n_cells"].values
    sizes = 40 + 10 * np.log2(n_cells)

    r_val = np.corrcoef(x, y)[0, 1] if len(x) > 1 else np.nan

    rng = np.random.default_rng(0)
    x_jit = x + rng.normal(scale=0.03, size=len(x))

    sc = axA.scatter(
        x_jit, y,
        s=sizes,
        c=n_cells,
        cmap="viridis",
        edgecolor="black",
        alpha=0.85
    )
    max_val = max(x.max(), y.max()) + 0.5
    axA.plot([0, max_val], [0, max_val], "k--", alpha=0.7,
             label="y = x (observed = expected)")
    axA.set_xlabel("Expected richness (null model)")
    axA.set_ylabel("Observed richness")
    axA.set_title(f"{dataset_label_prefix}) {dataset_name} — Obs vs Exp\n(r = {r_val:.2f})")
    cbar = plt.colorbar(sc, ax=axA)
    cbar.set_label("Clone size (n cells)")
    axA.legend(loc="upper left")

    # Panel B: ΔR histogram
    delta = df["delta_richness"]
    sns.histplot(delta, bins=10, kde=True, ax=axB,
                 color="skyblue", edgecolor="black")
    axB.axvline(0, color="gray", linestyle="--")
    axB.set_xlabel("Observed - Expected richness (ΔR)")
    axB.set_ylabel("Number of clones")
    axB.set_title(f"{chr(ord(dataset_label_prefix)+1)}) {dataset_name} — ΔR Distribution")

    # Panel C: ΔR bar plot
    df_bar = df.sort_values("delta_richness")
    x_idx = np.arange(len(df_bar))
    axC.bar(x_idx, df_bar["delta_richness"],
            color="lightsteelblue", edgecolor="black")
    axC.axhline(0, color="gray", linestyle="--", linewidth=1)

    labels = [
        f"{ck}\n(n={int(n)})" for ck, n in
        zip(df_bar["CloneKey"], df_bar["n_cells"])
    ]
    axC.set_xticks(x_idx)
    axC.set_xticklabels(labels, rotation=90, ha="right")
    axC.set_ylabel("Observed - Expected richness (ΔR)")
    axC.set_title(f"{chr(ord(dataset_label_prefix)+2)}) {dataset_name} — Clone-specific ΔR")


def make_figure2(homo_entry, het_entry, category_col):
    """
    Figure 2:
      A–C: Homo dataset (3 panels as in Figure 1 for real data)
      D–F: Hetero dataset (3 panels as in Figure 1 for real data)
    """
    if homo_entry is None or het_entry is None:
        print("[WARN] Cannot make Figure 2 — missing homo or het dataset.")
        return

    homo_df = homo_entry["per_clone_df"]
    het_df = het_entry["per_clone_df"]

    fig, axes = plt.subplots(2, 3, figsize=(24, 10), dpi=DPI)
    axes = axes.flatten()

    # Top row A–C: Homo
    _plot_three_panels_for_dataset(
        axes[0], axes[1], axes[2],
        homo_df,
        "A",
        "Homogeneous Simulated"
    )

    # Bottom row D–F: Hetero
    _plot_three_panels_for_dataset(
        axes[3], axes[4], axes[5],
        het_df,
        "D",
        "Heterogeneous Simulated"
    )

    plt.suptitle(
        f"Figure 2 — Simulated Datasets: Homogeneous vs Heterogeneous\nCategory: {category_col}",
        fontsize=16, y=1.03
    )
    plt.tight_layout()
    out_path = os.path.join(
        OUTPUT_DIR,
        f"Figure2_clone_richness_simulated_{category_col}.svg"
    )
    plt.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 2 to: {out_path}")


# FIGURE 3 — INTRA vs INTER (RANDOM vs CLONAL) — RICHNESS

def _paired_panel_richness(ax, df, nice_name, panel_label):
    """
    Draw a clean paired plot for richness: Random (expected) vs Clonal (observed).
    """
    expected = df["mean_null_richness"].values  # Random
    observed = df["observed_richness"].values   # Clonal

    # Paired t-test
    t_stat, p_val = stats.ttest_rel(observed, expected, nan_policy="omit")
    delta = observed - expected
    mean_delta = np.nanmean(delta)

    # Sort clones by expected for nicer layout
    tmp = df.copy()
    tmp["expected"] = expected
    tmp["observed"] = observed
    tmp = tmp.sort_values("expected").reset_index(drop=True)

    x_random = 1
    x_clonal = 0

    for _, row in tmp.iterrows():
        ax.plot([x_clonal, x_random], [row["observed"], row["expected"]],
        color="0.7", alpha=0.7, linewidth=1.5)

        ax.scatter([x_clonal], [row["observed"]], color="tab:red", s=25, zorder=3)
        ax.scatter([x_random], [row["expected"]], color="tab:blue", s=25, zorder=3)

        tmp = df.copy()
        tmp["expected"] = expected
        tmp["observed"] = observed
        tmp = tmp.sort_values("expected").reset_index(drop=True)
        
        p_txt = row.get("pH_two_tailed", np.nan)
        if np.isfinite(p_txt):
            ax.text(
                x_clonal + 0.05, row["observed"],
                f"p={p_txt:.3g}",
                va="center", ha="left",
                fontsize=8, color="black"
    )

    ax.set_xlim(-0.2, 1.2)
    ax.set_xticks([x_clonal, x_random])
    ax.set_xticklabels(["Clonal", "Random"])
    ax.set_ylabel("Richness (# morphotypes)")
    ax.set_title(
        f"{panel_label} {nice_name} — Random vs Clonal richness\n"
        f"t = {t_stat:.2f}, p = {p_val:.3f}, mean Δ = {mean_delta:.2f}"
    )
    ax.grid(axis="y", alpha=0.3)


def make_figure3_intra_vs_inter(all_entries, category_col):
    """
    Figure 3:
      A) Random vs Clonal richness — Real
      B) Random vs Clonal richness — Homogeneous simulated
      C) Random vs Clonal richness — Heterogeneous simulated
    """
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

    fig, axes = plt.subplots(1, len(ordered_labels),
                             figsize=(5 * len(ordered_labels), 6),
                             dpi=DPI)
    if len(ordered_labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, ordered_labels):
        nice_name, entry = mapping[label]
        df = entry["per_clone_df"]
        _paired_panel_richness(ax, df, nice_name, f"{label})")

    plt.suptitle(
        f"Figure 3 — Random vs Clonal Richness (All Datasets)\nCategory: {category_col}",
        fontsize=16,
        y=1.04
    )
    plt.tight_layout()
    out_path = os.path.join(
        OUTPUT_DIR,
        f"Figure3_random_vs_clonal_richness_{category_col}.svg"
    )
    plt.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 3 to: {out_path}")


# FIGURE 4 — SAME PLOTS BUT COLOURED BY CLONE SIZE (LINES BLACK)

def _paired_panel_richness_colored(ax, df, nice_name, panel_label):
    """
    Same as _paired_panel_richness but colour-code each clone's datapoints
    by clone size, with BLACK connecting lines (as requested).
    """
    expected = df["mean_null_richness"].values  # Random
    observed = df["observed_richness"].values   # Clonal
    n_cells = df["n_cells"].values

    t_stat, p_val = stats.ttest_rel(observed, expected, nan_policy="omit")
    delta = observed - expected
    mean_delta = np.nanmean(delta)

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
        # black line between Random and Clonal
        ax.plot(
            [x_random, x_clonal],
            [row["expected"], row["observed"]],
            color="black",
            alpha=0.5,
            linewidth=1.5
        )
        # coloured points at Random and Clonal positions
        ax.scatter(
            [x_random], [row["expected"]],
            color=color_point,
            s=30,
            zorder=3
        )
        ax.scatter(
            [x_clonal], [row["observed"]],
            color=color_point,
            s=30,
            zorder=3
        )

    ax.set_xlim(-0.2, 1.2)
    ax.set_xticks([x_random, x_clonal])
    ax.set_xticklabels(["Random", "Clonal"])
    ax.set_ylabel("Richness (# morphotypes)")
    ax.set_title(
        f"{panel_label} {nice_name} — Random vs Clonal\n"
        f"t = {t_stat:.2f}, p = {p_val:.3f}, mean Δ = {mean_delta:.2f}"
    )
    ax.grid(axis="y", alpha=0.3)

    # Colorbar for clone size
    sm = ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label("Clone size (n cells)")


def make_figure4_intra_vs_inter_colored(all_entries, category_col):
    """
    Figure 4:
      Same as Figure 3 but clones colour-coded by size.
      Datapoints are coloured by clone size, lines are black.
    """
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

    fig, axes = plt.subplots(1, len(ordered_labels),
                             figsize=(5 * len(ordered_labels), 6),
                             dpi=DPI)
    if len(ordered_labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, ordered_labels):
        nice_name, entry = mapping[label]
        df = entry["per_clone_df"]
        _paired_panel_richness_colored(ax, df, nice_name, f"{label})")

    plt.suptitle(
        f"Figure 4 — Random vs Clonal Richness (Colour-coded by clone size)\nCategory: {category_col}",
        fontsize=16,
        y=1.04
    )
    plt.tight_layout()
    out_path = os.path.join(
        OUTPUT_DIR,
        f"Figure4_random_vs_clonal_colored_{category_col}.svg"
    )
    plt.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 4 to: {out_path}")


# FIGURE 5 — RANDOM vs CLONAL HETEROGENEITY INDEX (NORMALIZED RICHNESS)

def _paired_panel_heterogeneity(ax, df, nice_name, panel_label):
    """
    Paired plot using heterogeneity index (normalized richness 0–1)

    Layout:
      LEFT  = Random
      RIGHT = Clonal

    Features:
      - Per-clone Monte Carlo p-values (pH_two_tailed) annotated to the right
        of the CLONAL datapoint
      - Thick mean line on the CLONAL column
    """
    N = df["n_cells"].values
    R_clonal = df["observed_richness"].values
    R_random = df["mean_null_richness"].values

    H_clonal = heterogeneity_index(R_clonal, N)
    H_random = heterogeneity_index(R_random, N)

    # Global paired t-test summary (kept for figure title)
    t_stat, p_val = stats.ttest_rel(H_clonal, H_random, nan_policy="omit")
    delta = H_clonal - H_random
    mean_delta = np.nanmean(delta)

    tmp = df.copy()
    tmp["H_random"] = H_random
    tmp["H_clonal"] = H_clonal

    # Clone-wise Monte Carlo p-values (if present)
    if "pH_two_tailed" not in tmp.columns:
        tmp["pH_two_tailed"] = np.nan

    # Sort by random H for cleaner stacking
    tmp = tmp.sort_values("H_random").reset_index(drop=True)

    # --- AXIS POSITIONS ---
    x_random = 0
    x_clonal = 1

    mean_H_random = np.nanmean(H_random)
    mean_H_clonal = np.nanmean(H_clonal)

    for _, row in tmp.iterrows():
        # connecting line
        ax.plot(
            [x_random, x_clonal],
            [row["H_random"], row["H_clonal"]],
            color="0.7",
            alpha=0.7,
            linewidth=1.5
        )

        # points
        ax.scatter([x_random], [row["H_random"]],
                   color="tab:blue", s=25, zorder=3)
        ax.scatter([x_clonal], [row["H_clonal"]],
                   color="tab:red", s=25, zorder=3)

        # p-value annotation (to the RIGHT of CLONAL point)
        p_txt = row["pH_two_tailed"]
        if np.isfinite(p_txt):
            ax.text(
                x_clonal + 0.06, row["H_clonal"],
                f"p={p_txt:.3g}",
                va="center", ha="left",
                fontsize=8, color="black"
            )

    # --- Thick mean paired line (Random → Clonal) ---
    ax.plot(
        [x_random, x_clonal],
        [mean_H_random, mean_H_clonal],
        color="black",
        linewidth=4.0,
        zorder=6
    )

    # Optional: emphasise mean points
    ax.scatter(
        [x_random, x_clonal],
        [mean_H_random, mean_H_clonal],
        color="black",
        s=60,
        zorder=7
    )

    ax.set_xlim(-0.2, 1.35)
    ax.set_xticks([x_random, x_clonal])
    ax.set_xticklabels(["Random", "Clonal"])
    ax.set_ylabel("Heterogeneity index H\n(normalized richness 0–1)")
    ax.set_title(
        f"{panel_label} {nice_name} — Random vs Clonal H\n"
        f"t = {t_stat:.2f}, p = {p_val:.3g}, mean ΔH = {mean_delta:.2f}"
    )
    ax.grid(axis="y", alpha=0.3)



def make_figure5_heterogeneity(all_entries, category_col):
    """
    Figure 5:
      A–C: Random vs Clonal heterogeneity index for all datasets.
    """
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

    fig, axes = plt.subplots(1, len(ordered_labels),
                             figsize=(5 * len(ordered_labels), 6),
                             dpi=DPI)
    if len(ordered_labels) == 1:
        axes = [axes]

    for ax, label in zip(axes, ordered_labels):
        nice_name, entry = mapping[label]
        df = entry["per_clone_df"]
        _paired_panel_heterogeneity(ax, df, nice_name, f"{label})")

    plt.suptitle(
        f"Figure 5 — Random vs Clonal Heterogeneity Index (All Datasets)\nCategory: {category_col}",
        fontsize=16,
        y=1.04
    )
    plt.tight_layout()
    out_path = os.path.join(
        OUTPUT_DIR,
        f"Figure5_random_vs_clonal_heterogeneity_{category_col}.svg"
    )
    plt.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Figure 5 to: {out_path}")


# SUPPLEMENTARY FIGURE — ΔR vs CLONE SIZE + HISTOGRAM (REAL)
# (NOT CALLED FROM main()

def make_supp_figure_delta_real(real_entry, category_col):
    """
    Supplementary: For real dataset, show
      Left: ΔR vs clone size
      Right: Histogram of ΔR
    """
    df = real_entry["per_clone_df"]
    delta = df["delta_richness"].values
    n_cells = df["n_cells"].values

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=DPI)

    # Left: scatter ΔR vs clone size
    sc = ax1.scatter(
        n_cells,
        delta,
        s=40 + 10 * np.log2(n_cells),
        c=delta,
        cmap="coolwarm",
        edgecolor="black",
        alpha=0.85
    )
    ax1.axhline(0, color="gray", linestyle="--", alpha=0.7)
    ax1.set_xlabel("Clone size (n cells)")
    ax1.set_ylabel("ΔR = Clonal - Random richness")
    ax1.set_title("ΔR vs clone size (Real Data)")
    cbar = plt.colorbar(sc, ax=ax1)
    cbar.set_label("ΔR")

    # Right: histogram of ΔR
    sns.histplot(delta, bins=10, kde=True, ax=ax2,
                 color="skyblue", edgecolor="black")
    ax2.axvline(0, color="gray", linestyle="--", alpha=0.7)
    ax2.set_xlabel("ΔR = Clonal - Random richness")
    ax2.set_ylabel("Number of clones")
    ax2.set_title("Distribution of ΔR (Real Data)")

    plt.suptitle(
        f"Supplementary — ΔR vs clone size and distribution (Real Dataset)\nCategory: {category_col}",
        fontsize=14,
        y=1.03
    )
    plt.tight_layout()
    out_path = os.path.join(
        OUTPUT_DIR,
        f"Supp_deltaR_real_{category_col}.svg"
    )
    plt.savefig(out_path, bbox_inches="tight")
    plt.show()
    print(f"[INFO] Saved Supplementary ΔR figure to: {out_path}")


# RUN A SINGLE DATASET
def run_single_dataset(excel_path, dataset_name, dataset_tag):
    """Load, clean, analyze, save, and return per-clone results for one dataset."""
    print("======================================================================")
    print(f"PROCESSING DATASET: {dataset_name}")
    print("======================================================================")

    try:
        print("[INFO] Loading Excel:", excel_path)
        df = pd.read_excel(excel_path)
        print("[INFO] Data loaded. Shape:", df.shape)
    except FileNotFoundError:
        print("[ERROR] Could not find file at:", excel_path)
        return None
    except Exception as e:
        print("[ERROR] Problem loading Excel:", e)
        return None

    # Filter to reconstructed cells if has_swc exists
    if "has_swc" in df.columns:
        initial_n = len(df)
        df = df[df["has_swc"] == True].copy()
        print(f"[INFO] Filtered to reconstructed cells using 'has_swc': "
              f"{len(df)} of {initial_n}")

    # Check required cols
    required_cols = [CLONE_COL, CATEGORY_COL]
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print("[ERROR] Missing required columns:", missing)
        print("       Available columns:", df.columns.tolist())
        return None

    # Drop NaNs
    df_clean = df.dropna(subset=[CLONE_COL, CATEGORY_COL]).copy()
    summarize_dataset(df_clean, CATEGORY_COL, dataset_name)

    # Analysis
    per_clone_df, global_info = analyze_clone_richness(df_clean, CATEGORY_COL, dataset_name)

    # Save per-dataset outputs
    save_results_and_report(per_clone_df, global_info, CATEGORY_COL, dataset_tag, dataset_name)

    return {
        "per_clone_df": per_clone_df,
        "global_info": global_info,
        "dataset_name": dataset_name,
        "dataset_tag": dataset_tag
    }

# MAIN

def main():
    # Dataset configuration (order: real, homo, het)
    dataset_configs = [
        {
            "excel_path": REAL_EXCEL_PATH,
            "dataset_name": "Real Data",
            "dataset_tag": "real",
            "color": COLOR_REAL
        },
        {
            "excel_path": HOMO_EXCEL_PATH,
            "dataset_name": "Homogeneous Simulated Dataset",
            "dataset_tag": "homo",
            "color": COLOR_HOMO
        },
        {
            "excel_path": HET_EXCEL_PATH,
            "dataset_name": "Heterogeneous Simulated Dataset",
            "dataset_tag": "het",
            "color": COLOR_HET
        }
    ]

    all_entries = []

    for cfg in dataset_configs:
        entry = run_single_dataset(
            cfg["excel_path"],
            cfg["dataset_name"],
            cfg["dataset_tag"]
        )
        if entry is not None:
            entry["color"] = cfg["color"]
            all_entries.append(entry)

    if not all_entries:
        print("[ERROR] No datasets were successfully processed.")
        return

    # Identify real/homo/het explicitly
    real_entry = next((e for e in all_entries if e["dataset_tag"] == "real"), None)
    homo_entry = next((e for e in all_entries if e["dataset_tag"] == "homo"), None)
    het_entry  = next((e for e in all_entries if e["dataset_tag"] == "het"), None)

    # Figure 1 — real + all datasets
    if real_entry is not None:
        make_figure1(real_entry, all_entries, CATEGORY_COL)
    else:
        print("[WARN] Real dataset missing — cannot make Figure 1.")

    # Figure 2 — simulated datasets
    make_figure2(homo_entry, het_entry, CATEGORY_COL)

    # Figure 3 — Random vs Clonal (richness)
    make_figure3_intra_vs_inter(all_entries, CATEGORY_COL)

    # Figure 4 — Random vs Clonal (richness, coloured by clone size, lines black)
    make_figure4_intra_vs_inter_colored(all_entries, CATEGORY_COL)

    # Figure 5 — Random vs Clonal (heterogeneity index, normalized richness)
    make_figure5_heterogeneity(all_entries, CATEGORY_COL)

    # Supplementary ΔR figure (real)

    print("\n======================================================================")
    print("ANALYSIS COMPLETE FOR ALL DATASETS!")
    print("======================================================================")
    print(f"All results saved to: {OUTPUT_DIR}")
    print("  - clone_richness_permutation_results_*_*.csv")
    print("  - clone_richness_permutation_report_*_*.txt")
    print("  - Figure1_clone_richness_real_and_all_*.svg")
    print("  - Figure2_clone_richness_simulated_*.svg")
    print("  - Figure3_random_vs_clonal_richness_*.svg")
    print("  - Figure4_random_vs_clonal_colored_*.svg")
    print("  - Figure5_random_vs_clonal_heterogeneity_*.svg")
    print("  # Supplementary figure function is defined but not called.")
    print("======================================================================")


if __name__ == "__main__":
    main()
