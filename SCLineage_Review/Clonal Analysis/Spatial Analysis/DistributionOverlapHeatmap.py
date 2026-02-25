#!/usr/bin/env python3
"""
Quantify overlap between AB distributions using KDE HDR regions and plot heatmaps.

For each reference file and each query file:
- Build a 2D KDE on a fixed grid.
- Find the density threshold corresponding to a target percentile (e.g. 70%).
- Define the HDR region: { (x,y) : KDE(x,y) >= threshold }.
- Compute overlap metrics between the two HDR regions.

Outputs
-------
1) Long-form CSV with columns:
   RefFile, QueryFile, Percentile,
   Area_ref, Area_query, Area_intersection, Area_union,
   Jaccard, Dice, Coverage_ref, Coverage_query

2) Wide CSV with Jaccard index matrix:
   rows = QueryFile, columns = RefFile

3) Heatmaps (PNG + SVG) for chosen metrics (e.g. Jaccard, Dice, Coverage_ref, Coverage_query).
"""

import os
import glob
import datetime
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.spatial.distance import pdist
from scipy.cluster.hierarchy import linkage, leaves_list

import matplotlib as mpl
import matplotlib.pyplot as plt

########################################################################################################################################################################################################
### User input

REF_DIR   = r"Z:\People\Sophie\10 AB Spatial Distribution Map\20251114 AB Spatial Map\20251114 Lumbar Matt and Giulia Counting\AB Files" # Path of the folder ("AB Files") containing the AB distributions (cell coordinates).
QUERY_DIR = r"Z:\People\Sophie\10 AB Spatial Distribution Map\20251114 AB Spatial Map\20251114 Lumbar Matt and Giulia Counting\AB Files" # Path of the folder containing the distributions you want to compare. If you want to look at the overlap between AB distributions, the path should be the same as the REF_DIR.

now = datetime.datetime.now()
curr_date = now.strftime("%Y%m%d_%H%M%S")
OUTPUT_DIR = rf"Z:\People\Sophie\10 AB Spatial Distribution Map\20260113 Clone Spatial Analysis\AB vs AB distribution comparison\Output_HDROverlap_{curr_date}"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Thoracic vs Lumbar scaling (if needed): set to 1/1.49 for T (as in MATLAB), else 1.0
X_SCALE = 1.0  # e.g., 1.0 for L, or (1/1.49) for T

# If the last N rows are landmarks to drop:
DROP_LAST_N = 0  # set 0 if your CSVs have no trailing landmarks

# Hemisphere handling (similar to your other script)
HEMISPHERES       = 'single'   # 'single' | 'LR' | 'random' | 'assigned'
SINGLE_HEMISPHERE = 'right'    # used only if HEMISPHERES='single'

# For HEMISPHERES='assigned' (optional); can be left empty if unused
LEFT_FILES  = []  # list of filename patterns for "forced left"
RIGHT_FILES = []  # list of filename patterns for "forced right"

AB_ORDER = [
    'dI5 (Lbx1+Tlx3+)', 'dI5 (Lmx1b+)', 'dI5 (Tlx3+)',
    'dI1 (Barhl1+)', 'dI2 (FoxD3+En1-)', 'dI3 (dorsal Isl+)',
    'dI3 (Isl+Tlx3+Lbx1-)',
    'dI6 (Dmrt3+)', 'V0 (Evx1+)', 'V1 (En1+)', 'V1 (En1+FoxD3+)',
    'V2 (Chx10+)', 'MN (Hb9+)', 'MN (ventral Isl+)', 'V3 (Nkx2.2+)'
]  # desired AB order, top→bottom on heatmap

APPEND_UNLISTED_ABS = True  # if True, ABs not in AB_ORDER are appended at the end

# KDE / grid
N_BINS    = 256
X_LIMS    = (0.0, 650.0)    # half-cord width used in your other script
Y_LIMS    = (-400.0, 500.0)

# Set X_GRID based on hemisphere mode
if HEMISPHERES == 'single':
    X_GRID = (0.0, X_LIMS[1]) if SINGLE_HEMISPHERE == 'right' else (-X_LIMS[1], 0.0)
else:
    X_GRID = (-X_LIMS[1], X_LIMS[1])

Y_GRID   = (Y_LIMS[0], Y_LIMS[1])
BW_METHOD = "scott"

# Target HDR percentile (e.g. 0.70 = 70% highest-density region)
TARGET_PERCENTILE = 0.70

MIN_KDE_POINTS = 5

# Heatmap appearance
FIG_SIZE_PX = (1500, 800)   # width, height in pixels
ANNOTATE_VALUES = False     # usually False for big matrices
CLUSTER_ROWS = True         # hierarchical cluster queries (rows)
CLUSTER_COLS = True         # hierarchical cluster refs (cols)

# Which metrics to plot as heatmaps
HEATMAP_METRICS = ["Area_intersection", "Jaccard", "Dice", "Coverage_ref", "Coverage_query"]

# --------------------------
# Matplotlib style (similar to your heatmap script)
# --------------------------
mpl.rcParams['svg.fonttype'] = 'none'  # keep text as text in SVG
mpl.rcParams['pdf.fonttype'] = 42      # TrueType
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']


# --------------------------
# Helpers
# --------------------------

def ensure_output_dir(path: str):
    os.makedirs(path, exist_ok=True)

def file_key_three_underscores(path: str) -> str:
    base = os.path.basename(path)
    name, _ = os.path.splitext(base)
    parts = name.split("_")
    return name if len(parts) < 3 else "_".join(parts[:3])

def file_matches_any(name: str, candidates: set[str]) -> bool:
    base = os.path.basename(name)
    stem, _ = os.path.splitext(base)
    pref3 = file_key_three_underscores(base)
    return any(k in (base, stem, pref3) for k in candidates)

def pick_xy_columns(df: pd.DataFrame):
    lower = {c.lower(): c for c in df.columns}
    for a, b in (("x", "y"), ("coord_x", "coord_y"), ("pos_x", "pos_y")):
        if a in lower and b in lower:
            return lower[a], lower[b]
    numeric = df.select_dtypes(include=["number"]).columns.tolist()
    if len(numeric) >= 2:
        return numeric[0], numeric[1]
    cols = df.columns.tolist()
    if len(cols) >= 2:
        return cols[0], cols[1]
    raise ValueError("Could not infer X/Y columns.")

def read_xy_points(csv_path: str, hemisphere_policy: str):
    df = pd.read_csv(csv_path)
    if DROP_LAST_N > 0 and len(df) > DROP_LAST_N:
        df = df.iloc[:-DROP_LAST_N, :].reset_index(drop=True)

    xcol, ycol = pick_xy_columns(df)
    x = df[xcol].astype(float).to_numpy()
    y = df[ycol].astype(float).to_numpy()

    # Thoracic/Lumbar scale (if used)
    x = x * X_SCALE

    if hemisphere_policy == 'single':
        if SINGLE_HEMISPHERE == 'right':
            x = np.abs(x)         # mirror everything to x >= 0
        else:
            x = -np.abs(x)        # mirror everything to x <= 0

    if hemisphere_policy == 'LR':
        hemi_col = next((c for c in df.columns if c.lower() == "hemisphere"), None)
        if hemi_col is None:
            raise ValueError(f"{os.path.basename(csv_path)} lacks 'Hemisphere' column for HEMISPHERES='LR'.")
        hemi = df[hemi_col].astype(str).str.strip().str.lower()
        sign = hemi.map(lambda h: -1.0 if h.startswith("l") else 1.0)
        x = x * sign.to_numpy(dtype=float)

    elif hemisphere_policy == 'random':
        x = x * (1.0 if np.random.rand() < 0.5 else -1.0)

    elif hemisphere_policy == 'assigned':
        force_left  = file_matches_any(csv_path, set(LEFT_FILES))
        force_right = file_matches_any(csv_path, set(RIGHT_FILES))
        if force_left and force_right:
            # If ambiguous, prefer right
            force_left = False
        if force_left:
            x = -x
        elif force_right:
            x = +x
        # else leave as-is

    return x, y

def kde2d_on_fixed_grid(x, y, x_grid, y_grid, n_bins=256, bw_method="scott"):
    gx = np.linspace(x_grid[0], x_grid[1], n_bins)
    gy = np.linspace(y_grid[0], y_grid[1], n_bins)
    GX, GY = np.meshgrid(gx, gy)
    kde = gaussian_kde(np.vstack([x, y]), bw_method=bw_method)
    Z = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)

    # Normalize to integrate to ~1 over the grid (uniform cell areas)
    dx = (gx[-1] - gx[0]) / (n_bins - 1)
    dy = (gy[-1] - gy[0]) / (n_bins - 1)
    Z = Z / (Z.sum() * dx * dy + 1e-16)
    return gx, gy, Z

def kde_percentile_levels(Z, percentiles):
    """
    Given a KDE grid Z and a list of percentiles (ascending, e.g. [0.7]),
    return density thresholds (same length as percentiles), where each threshold
    corresponds to the *highest-density region* containing that fraction of mass.
    """
    z = Z.ravel()
    z = z[np.isfinite(z)]
    z = z[z > 0]
    if z.size == 0:
        return None
    zs = np.sort(z)[::-1]            # descending density
    cdf = np.cumsum(zs) / zs.sum()   # cumulative mass from highest to lowest density

    levels = []
    for p in percentiles:
        idx = np.searchsorted(cdf, p)
        levels.append(float(zs[min(idx, zs.size - 1)]))
    return np.array(levels)

def order_labels_by_similarity(matrix: pd.DataFrame,
                               axis: int = 0,
                               metric: str = "euclidean",
                               method: str = "average") -> list[str]:
    """
    Order labels (rows or columns) so that similar profiles are adjacent,
    using hierarchical clustering.

    axis = 0 -> cluster rows, axis = 1 -> cluster columns.
    """
    if axis == 0:
        data = matrix.to_numpy(dtype=float)       # rows × cols
        labels = matrix.index
    else:
        data = matrix.to_numpy(dtype=float).T     # cols × rows
        labels = matrix.columns

    # Handle NaNs: replace with column means.
    col_means = np.nanmean(data, axis=0)
    col_means = np.where(np.isnan(col_means), 0.0, col_means)
    nan_rows, nan_cols = np.where(np.isnan(data))
    data[nan_rows, nan_cols] = col_means[nan_cols]

    if data.shape[0] <= 1:
        return list(labels)

    dist_vec = pdist(data, metric=metric)
    Z = linkage(dist_vec, method=method)
    leaf_order = leaves_list(Z)
    return labels[leaf_order].tolist()

def plot_overlap_heatmap(matrix: pd.DataFrame,
                         metric_name: str,
                         out_dir: str,
                         vmin: float = 0.0,
                         vmax: float = 1.0,
                         annotate: bool = False):
    """
    Plot a heatmap for a given overlap metric matrix (QueryFile × RefFile),
    in a style similar to your CloneDistributionHeatmap.
    """
    ensure_output_dir(out_dir)

    dpi = 100
    fig_w = max(FIG_SIZE_PX[0] / dpi, 7)
    fig_h = max(FIG_SIZE_PX[1] / dpi, 4)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    data = matrix.to_numpy(dtype=float)
    masked = np.ma.masked_invalid(data)
    cmap = plt.cm.viridis  # standard viridis, low=purple, high=yellow
    cmap = cmap.copy()
    cmap.set_bad(color="#d9d9d9")  # light gray for NaN

    im = ax.imshow(masked, aspect="auto", interpolation="nearest",
                   cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_xticklabels(matrix.columns, rotation=90, fontsize=9)
    ax.set_yticklabels(matrix.index, fontsize=9)

    ax.set_xlabel("Reference distributions (RefFile)", fontsize=11)
    ax.set_ylabel("Query distributions (QueryFile)", fontsize=11)
    ax.set_title(f"{metric_name} overlap of {int(TARGET_PERCENTILE*100)}% HDR regions")

    # Minor gridlines
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    if annotate and matrix.size <= 3600:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix.iat[i, j]
                if not np.isnan(val):
                    ax.text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label(metric_name)

    fig.tight_layout()
    base = os.path.join(out_dir, f"HDR{int(TARGET_PERCENTILE*100)}_{metric_name}_heatmap")
    fig.savefig(base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(base + ".svg", format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)


# --------------------------
# Main
# --------------------------

def main():
    ensure_output_dir(OUTPUT_DIR)

    # Precompute grid cell area
    gx = np.linspace(X_GRID[0], X_GRID[1], N_BINS)
    gy = np.linspace(Y_GRID[0], Y_GRID[1], N_BINS)
    dx = (gx[-1] - gx[0]) / (N_BINS - 1)
    dy = (gy[-1] - gy[0]) / (N_BINS - 1)
    cell_area = dx * dy

    percentile_list = [TARGET_PERCENTILE]

    # --- Collect reference files (build KDEs + HDR masks) ---
    ref_paths = sorted(glob.glob(os.path.join(REF_DIR, "*.csv")))
    if not ref_paths:
        raise FileNotFoundError(f"No reference CSVs found in {REF_DIR}")

    refs = []  # list of dicts: {name, Z, thr, mask, area}
    for rp in ref_paths:
        base = os.path.basename(rp)
        try:
            xr, yr = read_xy_points(rp, HEMISPHERES)
        except Exception as e:
            print(f"[REF SKIP READ] {base}: {e}")
            continue

        valid = np.isfinite(xr) & np.isfinite(yr)
        xr = xr[valid]
        yr = yr[valid]

        if xr.size < MIN_KDE_POINTS:
            print(f"[REF SKIP KDE] {base}: too few points ({xr.size}).")
            continue

        gx, gy, Z = kde2d_on_fixed_grid(xr, yr, X_GRID, Y_GRID,
                                        n_bins=N_BINS,
                                        bw_method=BW_METHOD)
        levels = kde_percentile_levels(Z, percentile_list)
        if levels is None:
            print(f"[REF SKIP KDE] {base}: KDE levels undefined.")
            continue
        thr = levels[0]

        mask = Z >= thr
        area = mask.sum() * cell_area

        refs.append(dict(
            name=os.path.splitext(base)[0],
            Z=Z,
            thr=thr,
            mask=mask,
            area=area,
        ))

    if not refs:
        raise RuntimeError("No valid reference KDEs could be built.")

    # --- Collect query files (build KDEs + HDR masks) ---
    q_paths = sorted(glob.glob(os.path.join(QUERY_DIR, "*.csv")))
    if not q_paths:
        raise FileNotFoundError(f"No query CSVs found in {QUERY_DIR}")

    queries = []  # list of dicts: {name, Z, thr, mask, area}
    for qp in q_paths:
        base = os.path.basename(qp)
        try:
            xq, yq = read_xy_points(qp, HEMISPHERES)
        except Exception as e:
            print(f"[QUERY SKIP READ] {base}: {e}")
            continue

        valid = np.isfinite(xq) & np.isfinite(yq)
        xq = xq[valid]
        yq = yq[valid]

        if xq.size < MIN_KDE_POINTS:
            print(f"[QUERY SKIP KDE] {base}: too few points ({xq.size}).")
            continue

        gx, gy, Z = kde2d_on_fixed_grid(xq, yq, X_GRID, Y_GRID,
                                        n_bins=N_BINS,
                                        bw_method=BW_METHOD)
        levels = kde_percentile_levels(Z, percentile_list)
        if levels is None:
            print(f"[QUERY SKIP KDE] {base}: KDE levels undefined.")
            continue
        thr = levels[0]

        mask = Z >= thr
        area = mask.sum() * cell_area

        queries.append(dict(
            name=os.path.splitext(base)[0],
            Z=Z,
            thr=thr,
            mask=mask,
            area=area,
        ))

    if not queries:
        raise RuntimeError("No valid query KDEs could be built.")

    # --- Compute overlap metrics for each ref/query pair ---
    rows = []
    for ref in refs:
        for q in queries:
            m_ref = ref["mask"]
            m_q   = q["mask"]

            inter = (m_ref & m_q).sum() * cell_area
            union = (m_ref | m_q).sum() * cell_area

            area_ref   = ref["area"]
            area_query = q["area"]

            if union > 0:
                jaccard = inter / union
            else:
                jaccard = np.nan

            denom_dice = area_ref + area_query
            if denom_dice > 0:
                dice = 2.0 * inter / denom_dice
            else:
                dice = np.nan

            cov_ref   = inter / area_ref   if area_ref   > 0 else np.nan
            cov_query = inter / area_query if area_query > 0 else np.nan

            rows.append(dict(
                RefFile            = ref["name"],
                QueryFile          = q["name"],
                Percentile         = TARGET_PERCENTILE,
                Area_ref           = area_ref,
                Area_query         = area_query,
                Area_intersection  = inter,
                Area_union         = union,
                Jaccard            = jaccard,
                Dice               = dice,
                Coverage_ref       = cov_ref,
                Coverage_query     = cov_query,
            ))

    df_long = pd.DataFrame(rows)
    long_csv = os.path.join(OUTPUT_DIR, f"overlap_HDR{int(TARGET_PERCENTILE*100)}_long_{curr_date}.csv")
    df_long.to_csv(long_csv, index=False)
    print(f"[SAVED] {long_csv}")

    # --- Build matrices for metrics of interest (queries × refs) ---
    metrics = ["Jaccard", "Dice", "Coverage_ref", "Coverage_query",
               "Area_ref", "Area_query", "Area_intersection", "Area_union"]

    matrices = {}
    for m in metrics:
        mat = df_long.pivot(index="QueryFile", columns="RefFile", values=m)
        matrices[m] = mat

    # Save Jaccard matrix as CSV (as before)
    jaccard_mat = matrices["Jaccard"]
    mat_csv = os.path.join(OUTPUT_DIR, f"overlap_HDR{int(TARGET_PERCENTILE*100)}_Jaccard_matrix_{curr_date}.csv")
    jaccard_mat.to_csv(mat_csv)
    print(f"[SAVED] {mat_csv}")

 # --- Ordering of rows/columns (optional) ---
    # We apply the same row/column ordering to all metric matrices.
    base_mat = matrices["Jaccard"]

    if AB_ORDER:
        # Use the user-specified AB_ORDER where possible
        # Rows (QueryFile)
        listed_rows = [a for a in AB_ORDER if a in base_mat.index]
        if APPEND_UNLISTED_ABS:
            extra_rows = [a for a in base_mat.index if a not in listed_rows]
            row_order = listed_rows + extra_rows
        else:
            row_order = listed_rows

        # Columns (RefFile)
        listed_cols = [a for a in AB_ORDER if a in base_mat.columns]
        if APPEND_UNLISTED_ABS:
            extra_cols = [a for a in base_mat.columns if a not in listed_cols]
            col_order = listed_cols + extra_cols
        else:
            col_order = listed_cols
    else:
        # Fall back to clustering or original ordering
        if CLUSTER_ROWS and base_mat.shape[0] > 1:
            row_order = order_labels_by_similarity(base_mat, axis=0)
        else:
            row_order = base_mat.index.tolist()

        if CLUSTER_COLS and base_mat.shape[1] > 1:
            col_order = order_labels_by_similarity(base_mat, axis=1)
        else:
            col_order = base_mat.columns.tolist()

    # Apply ordering
    for m in matrices:
        matrices[m] = matrices[m].reindex(index=row_order, columns=col_order)

    # --- Plot heatmaps for selected metrics ---
    for metric_name in HEATMAP_METRICS:
        if metric_name not in matrices:
            continue
        mat = matrices[metric_name]

        # For Jaccard/Dice/Coverage* we know the natural range is [0,1]
        if metric_name in ("Jaccard", "Dice", "Coverage_ref", "Coverage_query"):
            vmin, vmax = 0.0, 0.7
        else:
            # For raw areas, auto-scale
            vmin, vmax = np.nanmin(mat.to_numpy()), np.nanmax(mat.to_numpy())
            if not np.isfinite(vmin): vmin = 0.0
            if not np.isfinite(vmax): vmax = 1.0

        plot_overlap_heatmap(mat, metric_name, OUTPUT_DIR,
                             vmin=vmin, vmax=vmax,
                             annotate=ANNOTATE_VALUES)


if __name__ == "__main__":
    main()
