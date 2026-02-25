"""
Summarise per-clone AB overlap from per-cell KDE percentile calls and
generate heatmaps, per-cell membership tables, and mixing statistics.

Overview
--------
This script takes as input a *per-cell* AB membership table (the output of the
"query_cells_percentile_membership" script) and collapses it to a
*per-clone × AB* matrix of minimum KDE percentiles. It then:

1. Builds an AB (rows) × Clone (columns) matrix where each entry is the
   **smallest** percentile (best overlap) reached by any cell in that clone
   for that AB distribution.
2. Optionally reorders ABs and clones (hierarchical clustering or overlap
   score) to make structure in the matrix visually clearer.
3. Produces:
   - Heatmaps of the minimum percentiles per clone × AB,
   - A per-cell AB membership “table” (black/white) for all cells,
   - Individual per-clone cell×AB membership plots,
   - CSVs with the AB×clone matrix and per-clone mixing metrics,
   - Stacked bar plots and histograms describing cardinal class mixing,
   - A scatterplot of clone size vs. degree of cardinal class mixing.

The script is meant to be run after KDE-based scoring of individual cells
against reference AB distributions, and uses the same AB labels.

Expected input
--------------
INPUT_PATH
    Path to a CSV file produced by the KDE scoring script, typically named
    e.g. `query_cells_percentile_membership_YYYYMMDD_HHMMSS.csv`.

This input is expected to contain (at least) the following columns:

- Metadata / identifiers:
  - `QueryFile` (preferred) or `Clone`:
        clone identifier for each cell.
  - `PointIndex`:
        per-cell index within each query file (used to count cells/clone).
  - `TM` (optional):
        tamoxifen timepoint; used to colour the TM strip in the second heatmap.
  - `ID` (optional):
        per-cell ID; used in per-cell plots.

- Per-cell coordinates (X, Y, etc.):
  - Present but not used directly here, other than for diagnostics.

- AB membership columns:
  - All AB distribution columns are assumed to start at a fixed position:
        columns 9..end (0-based index), i.e. from the 10th column onward.
  - Each AB column contains strings such as:
        - `"Outside"`   → cell lies outside all query percentiles, or
        - `">=44%"`     → cell lies within the 44% HDR or better.
  - These are converted to integers representing percentiles:
        - `">=44%"` → 44
        - `"Outside"` → OUTSIDE_SENTINEL (default 100), or NaN if IGNORE_OUTSIDE=True.

Key user settings
-----------------
INPUT_PATH
    Per-cell AB membership CSV (from KDE script).

OUTPUT_DIR
    Directory where all outputs (CSVs, PNGs, SVGs) will be written.
    A timestamped subfolder is appended when the script runs.

OUTSIDE_SENTINEL
    Numeric value used when parsing `"Outside"` if IGNORE_OUTSIDE=False.
    Larger values correspond to “worse” overlap, so 100 means
    “as bad as it gets”.

IGNORE_OUTSIDE
    If True, `"Outside"` is treated as NaN and excluded from minima;
    if False, `"Outside"` is treated as OUTSIDE_SENTINEL.

EXCLUDE_SINGLE_CELL_CLONES
    If True, clones with only a single cell (based on `PointIndex` or row
    counts) are removed from all downstream analyses.

AB_ORDER
    Explicit top-to-bottom order for AB labels in the heatmaps. ABs not
    listed here can be appended alphabetically if APPEND_UNLISTED_ABS=True.

APPEND_UNLISTED_ABS
    If True, ABs not in AB_ORDER are appended (alphabetically) below the
    specified list. If False, only ABs in AB_ORDER are shown.

SORT_CLONES
    How columns (clones) are ordered in the AB×clone matrix:
      - "by_similarity": hierarchical clustering on clone profiles
                         (default here),
      - "by_overlap":    sort by an overlap score derived from the matrix
                         (see CLONE_OVERLAP_METHOD),
      - True:            alphabetical by clone name,
      - list:            custom canonical order of clones,
      - anything else:   leave input order unchanged.

CLONE_OVERLAP_METHOD
    When SORT_CLONES == "by_overlap", defines how overlap per clone is
    summarised from the AB×clone matrix:
      - "mean" (default), "median", "sum", or "max".
    Larger scores correspond to worse (higher) percentiles; smaller scores
    correspond to better (lower) percentiles.

MIX_THRESHOLD
    Percentile cutoff used to define whether a given AB (or cardinal class)
    is considered “present” in a clone. Typically 70:
      - values < MIX_THRESHOLD → AB present,
      - values ≥ MIX_THRESHOLD or NaN → considered absent for mixing metrics.

CELL_MEMBERSHIP_THRESHOLD
    Threshold for per-cell membership plots. For each cell×AB:
      - value ≤ CELL_MEMBERSHIP_THRESHOLD → black (present),
      - value > threshold or NaN → white/transparent.

LINEAGE_COMP, DV_GROUPS
    Lists of AB label groups defining lineage or dorsal/ventral groupings.
    Used to classify clones and build histograms of AB combinations.

Processing steps
----------------
1) Read and filter clones
   - The input CSV is read with automatic delimiter detection.
   - The script chooses the clone identifier column:
       - `clone_col = "QueryFile"` if available,
       - otherwise `clone_col = "Clone"`.
   - If EXCLUDE_SINGLE_CELL_CLONES is True, clones with only one unique
     `PointIndex` (or one row) are removed.

2) Identify AB columns and parse values
   - AB columns are taken as all columns from position 9 to the end
     (`df.columns[9:]`).
   - Each AB cell is parsed by `parse_ab_value()`:
       - `"Outside"` → NaN (if IGNORE_OUTSIDE) or OUTSIDE_SENTINEL,
       - strings like `">=44%"` → integer 44,
       - anything unparsable → NaN.
   - A numeric copy `df_num` is created with these parsed AB values.

3) Per-clone minimum percentile matrix (AB × Clone)
   - For each clone, the script computes the **minimum** percentile across
     all its cells for each AB column:
       grouped_min = df_num.groupby(clone_col)[ab_cols].min()
   - This is transposed to an AB (rows) × Clone (columns) matrix:
       matrix = grouped_min.T
   - Optionally, AB labels can be simplified or collapsed:
       - `simplify_ab_label()` can strip prefixes/suffixes (disabled here
         by KEEP_FULL_AB_LABELS=True).
       - If multiple AB columns collapse to the same label, they are merged
         by taking the minimum across duplicates.

4) Apply AB and clone ordering
   - AB rows:
       - If AB_ORDER is set:
           - Keep ABs in AB_ORDER that actually appear in the matrix,
           - Optionally append any unlisted ABs alphabetically at the end.
       - Otherwise: ABs are sorted alphabetically.
   - Clone columns:
       - "by_similarity": clones are ordered via hierarchical clustering
         on their AB profiles (NaNs are filled with per-AB means before
         distance computation).
       - "by_overlap": clones are ordered by an overlap score computed
         by `clone_scores_by_overlap()`.
       - True/list/other: handled as described under SORT_CLONES.

5) Per-cell AB membership table (all clones)
   - For each row in `df_num`, a cell label "Clone | ID" (or "Clone | idx=…")
     is created.
   - AB values are thresholded to a binary matrix:
       - value ≤ CELL_MEMBERSHIP_THRESHOLD → 1,
       - value > threshold → 0,
       - NaNs preserved.
   - This is arranged as an AB (rows) × Cell (columns) DataFrame and
     optionally given the same AB ordering as the main heatmap.
   - The table is:
       - Saved as CSV (`CELL_TABLE_CSV`),
       - Visualised as an image where:
           - black squares mark cell×AB pairs that satisfy the threshold,
           - NaNs are transparent,
           - alternating vertical bands mark different clones.

6) Per-clone mini cell×AB membership plots
   - For each clone separately, an AB×Cell binary matrix (as above) is built
     and plotted.
   - Each plot shows only the cells from that clone and uses the same AB
     ordering as the global heatmap.
   - Saved as PNG and SVG with filenames
     `cell_AB_membership_<clone>.png/svg`.

7) Cardinal class mixing metrics (per clone)
   - Using the AB×Clone minimum percentile matrix:
       - AB labels are mapped to **cardinal classes** via their first
         three characters (e.g. "dI1", "dI2", "V1", "MN", "V3", …).
       - For each clone:
           - ABs with value < MIX_THRESHOLD are considered “present”.
           - The **DegreeMix** is the number of *distinct* cardinal classes
             with at least one present AB.
           - IsMixed is True if DegreeMix > 1.
           - CloneSize is the number of cells in that clone (counting
             unique `PointIndex` where available).
   - These are collected into a per-clone table and saved as CSV
     (`MIXING_CSV`).

8) D/V mixing flag and clone-level heatmap CSV
   - A clone×AB version of the matrix is assembled (rows = clones, cols = ABs).
   - Using DV_GROUPS:
       - Group 0: “other cardinal classes”,
       - Group 1: “dI5” class.
   - For each clone, the script tests whether it has:
       - at least one dI5 AB present (< MIX_THRESHOLD), and
       - at least one non-dI5 AB present.
   - A "D/V mixing" flag ("yes"/"no") is added as the first column of the
     clone×AB table, then the table is saved as
     `query_cells_percentile_membership_clone_heatmap_DVmix.csv`.

9) Main heatmaps (AB × Clone minimum percentiles)
   - A masked array of the AB×Clone matrix (0–100, NaN allowed) is plotted
     with an inverted viridis colormap (yellow = low/minimal percentiles,
     purple = high/worse percentiles; NaNs in light gray).
   - Two versions are created:
       1) **Full clone labels**:
          - X-axis tick labels show "CloneName (n=###)".
       2) **Clone size only**:
          - X-axis tick labels show just the clone size (number of cells),
          - A coloured strip below the heatmap encodes TM per clone
            (E9.5/E10.5/E11.5 colours via `cap.ColPals`).
   - Both versions include:
       - AB labels on the y-axis,
       - minor grid lines to delineate cells,
       - a colourbar labelled "Minimum KDE percentile (0–100)".
   - Saved as PNG and SVG: `HEATMAP_PNG`, `HEATMAP_N_ONLY_PNG`.

10) Proportion of mixed cardinal classes
    - Using the per-clone IsMixed flag:
        - A stacked bar plot shows the proportion of clones that are
          "Mixed" (DegreeMix > 1) vs "Not mixed".
        - SVG is generated via `cap.snsStackedBarPlot`, and a PNG version
          is drawn with plain matplotlib.
    - Counts, proportions, and per-clone flags for dI5 + other mixing
      (below) are also exported as CSVs.

11) Proportion of clones with dI5 + any other cardinal class
    - Using the same MIX_THRESHOLD and cardinal class mapping:
        - For each clone, the script checks whether "dI5" is present and at
          least one other cardinal class is present at < MIX_THRESHOLD.
        - A stacked bar plot shows the proportion of clones with
          "dI5 + other" vs "No dI5 + other".
        - SVG via `cap.snsStackedBarPlot`, PNG via matplotlib.
    - CSVs are written for:
        - counts per category,
        - proportions,
        - a per-clone flag table.

12) Histograms of AB/lineage combinations (mixed clones only)
    - The helper `plot_lineage_histogram_from_groups()`:
        - Takes the AB×Clone matrix and a list of AB-label groups
          (LINEAGE_COMP or DV_GROUPS),
        - First thresholds presence at mix_threshold (< MIX_THRESHOLD),
        - Restricts to clones that:
            - have at least `min_abs_mixed` ABs present, and
            - have at least `min_classes_mixed` cardinal classes present
              (i.e. “mixed clones”).
        - Classifies each clone as belonging to a single AB-group subset
          or "Other".
        - Plots a histogram of counts per category and saves PNG + SVG.
    - Called twice:
        1) For LINEAGE_COMP with labels like
           "dI1/dI2/dI3", "dI5", "dI6/V0/V1", "V2/MN", "MN/V3".
        2) For DV_GROUPS with labels "dI1/dI2/dI3/dI6/V0/V1/V2/MN/V3" vs "dI5".

13) Scatter plot: clone size vs degree of mixing
    - A scatter plot is generated with:
        - x = CloneSize (cells),
        - y = DegreeMix (# cardinal classes with min percentile < MIX_THRESHOLD).
    - SVG created via `cap.snsScatterPlot`, PNG via seaborn/matplotlib.

Outputs
-------
All outputs are written to OUTPUT_DIR. Main files include:

- AB×Clone minimum percentile matrix:
    - `MATRIX_CSV`:
        AB (rows) × Clone (columns), values = min percentile (0–100, NaN).
- Per-clone mixing summary:
    - `MIXING_CSV`:
        per-clone table with DegreeMix, IsMixed, CloneSize.
- Clone-level heatmap table with D/V mixing:
    - `query_cells_percentile_membership_clone_heatmap_DVmix.csv`
- Global per-cell AB membership:
    - `CELL_TABLE_CSV`:
        AB (rows) × Cell (cols), entries 0/1/NaN after thresholding.
    - `CELL_TABLE_PNG`, `CELL_TABLE_SVG`:
        graphical black/white table.
- Per-clone cell×AB membership plots:
    - `cell_AB_membership_<Clone>.png/svg` for each clone.
- Heatmaps of min KDE percentiles:
    - `HEATMAP_PNG`, `HEATMAP_PNG.svg`
    - `HEATMAP_N_ONLY_PNG`, `HEATMAP_N_ONLY_PNG.svg`
- Proportion plots and data:
    - Mixed vs Not mixed cardinal classes:
        `proportion_mixed_cardinal_classes*.svg/png` + CSVs.
    - dI5 + any other cardinal class:
        `proportion_dI5_plus_other*.svg/png` + CSVs.
- Lineage and D/V histograms:
    - `lineage_subdivision_hist*.svg/png`
    - `dv_groups_hist*.svg/png`
- Clone size vs mixing scatter:
    - `clone_size_vs_degree_mixing*.svg/png`

Notes and assumptions
---------------------
- AB columns are assumed to start at a fixed position (column index 9).
  If the structure of the input CSV changes, this may need to be updated.
- "Outside" is treated as a worst-overlap state via OUTSIDE_SENTINEL, or aswith _tm_cel
  missing (NaN) if IGNORE_OUTSIDE=True.
- Lower percentiles correspond to **better** overlap (denser regions); all
  minima and thresholds are interpreted with that convention.
- Mixing metrics depend on both MIX_THRESHOLD and the mapping from AB labels
  to cardinal classes (first 3 characters); changing label conventions will
  change the grouping and results.
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import sys
sys.path.append(r"Z:\People\Sophie\9 Mouse SC lineage paper\20260126 Data Analysis and Plots\Code Versions\ClonalAnalysis") # File path of ClonalAnalysisPlots
import ClonalAnalysisPlots as cap
import datetime
import seaborn as sns

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import pdist


########################################################################################################################################################################################################
### User input
# Find current date and time, which will be used in file naming
now = datetime.datetime.now()
curr_date = now.strftime("%Y%m%d_%H%M%S")
    
INPUT_PATH  = r"Z:\People\Sophie\10 AB Spatial Distribution Map\20260113 Clone Spatial Analysis\AB vs Clone distribution comparison\Output_Neurons only 80% complete_20260115_155549\query_cells_percentile_membership_20260115_155928.csv" # Path of query_cells_percentile_membership csv file (output of CompareDistributions)
OUTPUT_DIR  = r"Z:\People\Sophie\10 AB Spatial Distribution Map\20260113 Clone Spatial Analysis\AB vs Clone distribution comparison\Output_Neurons only 80% complete_20260115_155549" # Desired output folder
# Output file naming
HEATMAP_PNG = "query_cells_percentile_membership_heatmap.png"
MATRIX_CSV  = "query_cells_percentile_membership_matrix.csv"
MIXING_CSV  = "query_cells_percentile_membership_clone_mixing.csv"
HEATMAP_N_ONLY_PNG = "query_cells_percentile_membership_heatmap_n_only.png"

# How to treat "Outside"
OUTSIDE_SENTINEL = 100      # value to use when not ignoring "Outside"
IGNORE_OUTSIDE   = False    # if True, treat "Outside" as NaN and ignore in mins

# Figure size in pixels (width, height)
FIG_SIZE_PX = (2000, 800)

# Whether to annotate values (may clutter for large matrices)
ANNOTATE_VALUES = False

# Clone sorting:
# - Set to True to sort alphabetically
# - Set to a list for a custom order
# - Set to "by_overlap" to sort by overall overlap (largest→smallest)
SORT_CLONES = "by_similarity"

#co.CUSTOM_CLONE_ORDER

# Overall-overlap method when SORT_CLONES == "by_overlap":
#   "mean" (default), "median", "sum", or "max"
CLONE_OVERLAP_METHOD = "mean"

# AB label simplification like "Lumbar_Barhl1_concatenated" -> "Barhl1"
AB_LABEL_SIMPLIFIER = re.compile(r'^(?:.*?_)?([^_]+?)(?:_concatenated)?$')
KEEP_FULL_AB_LABELS = True  # set True to keep full column names

# ABs
AB_ORDER =  ['dI5 (Lbx1+Tlx3+)', 'dI5 (Lmx1b+)', 'dI5 (Tlx3+)','dI1 (Barhl1+)', 'dI2 (FoxD3+En1-)', 'dI3 (dorsal Isl+)', 'dI3 (Isl+Tlx3+Lbx1-)', 
            'dI6 (Dmrt3+)', 'V0 (Evx1+)', 'V1 (En1+)', 'V1 (En1+FoxD3+)', 'V2 (Chx10+)', 'MN (Hb9+)', 'MN (ventral Isl+)', 'V3 (Nkx2.2+)'] # Top to bottom order of AN display on y axis of heatmap

# ['dI5 (Lbx1+Tlx3+)', 'dI5 (Lmx1b+)', 'dI5 (Tlx3+)','dI1 (Barhl1+)', 'dI2 (FoxD3+En1-)', 'dI3 (dorsal Isl+)', 'dI3 (Isl+Tlx3+Lbx1-)', 'dI6 (Dmrt3+)', 'V0 (Evx1+)', 'V1 (En1+)', 'V1 (En1+FoxD3+)', 'V2 (Chx10+)', 'V3 (Nkx2.2+)', 'MN (Hb9+)', 'MN (ventral Isl+)']

APPEND_UNLISTED_ABS = True   # if True, append ABs not in AB_ORDER (alphabetically). If False, show only AB_ORDER.
LINEAGE_COMP = [['dI1 (Barhl1+)', 'dI2 (FoxD3+En1-)', 'dI3 (dorsal Isl+)', 'dI3 (Isl+Tlx3+Lbx1-)'], 
                ['dI5 (Lbx1+Tlx3+)', 'dI5 (Lmx1b+)', 'dI5 (Tlx3+)'], 
                ['dI6 (Dmrt3+)', 'V0 (Evx1+)', 'V1 (En1+)', 'V1 (En1+FoxD3+)'], 
                ['V2 (Chx10+)','MN (Hb9+)', 'MN (ventral Isl+)'],
                ['V3 (Nkx2.2+)']]
DV_GROUPS = [['dI1 (Barhl1+)', 'dI2 (FoxD3+En1-)', 'dI3 (dorsal Isl+)', 'dI3 (Isl+Tlx3+Lbx1-)', 'dI6 (Dmrt3+)', 'V0 (Evx1+)', 'V1 (En1+)', 'V1 (En1+FoxD3+)',
              'V2 (Chx10+)','MN (Hb9+)', 'MN (ventral Isl+)', 'V3 (Nkx2.2+)'],
            ['dI5 (Lbx1+Tlx3+)', 'dI5 (Lmx1b+)', 'dI5 (Tlx3+)']]

# Do you want to exclude clones with only one cell?
EXCLUDE_SINGLE_CELL_CLONES = True 

PLOT_GLOBAL_CELL_TABLE = False 

# Threshold for an AB (or cardinal class) to be considered as “present” in a clone. Corresponds to a percentile
MIX_THRESHOLD = 70

# per-cell AB membership table output
CELL_TABLE_PNG       = "cell_AB_membership_table.png"
CELL_TABLE_SVG       = "cell_AB_membership_table.svg"
CELL_TABLE_CSV       = "cell_AB_membership_table_binary.csv"
CELL_MEMBERSHIP_THRESHOLD = 70   # black if percentile <= this value

# -------------------------
# Predicted identity settings
# -------------------------
PREDICTED_ID_COLNAME = "Predicted identity (single)"   # name of new column in df_num
PREDICT_TIE_SEED = 123                        # reproducible random tie-breaking

# -------------------------
# Per-cell *colour* heatmap settings
# -------------------------
PLOT_GLOBAL_CELL_TABLE_COLOURED = True   # new: coloured per-cell heatmap for all cells
PLOT_PER_CLONE_CELL_TABLE_COLOURED = True  # new: coloured per-cell heatmap per clone

CELL_COLOUR_VMIN = 0
CELL_COLOUR_VMAX = 100

# =========================
# ---- IMPLEMENTATION -----
# =========================

# Keep text as text (not paths) in SVG/PDF/PS
mpl.rcParams['svg.fonttype'] = 'none'  # <-- critical for Illustrator
mpl.rcParams['pdf.fonttype'] = 42      # TrueType; helps if you also export PDF
mpl.rcParams['ps.fonttype']  = 42

# Pick a common font that Illustrator will have
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']  # or 'Helvetica', etc.

def sanitize_filename(name: str) -> str:
    """Make a string safe to use as a filename."""
    return re.sub(r"[^\w\-\.]+", "_", str(name))

def canonical_clone_name(name: str) -> str:
    """
    Strip the trailing timestamp from clone names of the form
    'S017_a1_c3_20251119_091750' by removing the last 16 characters.

    If the name is shorter than 16 chars, return it unchanged.
    """
    if not isinstance(name, str):
        name = str(name)
    return name[:-16] if len(name) > 16 else name

def simplify_ab_label(col_name: str) -> str:
    if KEEP_FULL_AB_LABELS:
        return col_name
    #m = AB_LABEL_SIMPLIFIER.match(col_name)
    return m.group(1) if m else col_name

def ensure_output_dir(path: str):
    os.makedirs(path, exist_ok=True)

def read_table_auto(path: str) -> pd.DataFrame:
    # Auto-detect comma vs tab
    return pd.read_csv(path, sep=None, engine="python")

def parse_ab_value(x):
    """Parse an AB cell value: 'Outside' or '>=44%'. Return int or NaN."""
    if isinstance(x, str):
        xs = x.strip()
        if xs.lower() == "outside":
            return np.nan if IGNORE_OUTSIDE else OUTSIDE_SENTINEL
        if xs.startswith(">=") and xs.endswith("%"):
            inside = xs[2:-1].strip()
            try:
                return int(float(inside))
            except ValueError:
                return np.nan
    try:
        return int(float(x))
    except Exception:
        return np.nan
def choose_min_identity_random_tie(row_vals: pd.Series, rng: np.random.Generator) -> str | float:
    """
    Given a Series of AB percentiles for one cell (index = AB labels),
    return the AB label corresponding to the minimum percentile.
    If multiple ABs tie for the minimum, break ties randomly.
    Returns np.nan if all values are NaN.
    """
    # drop NaNs
    s = row_vals.dropna()
    if s.empty:
        return np.nan

    min_val = s.min()
    tied = s.index[s == min_val].tolist()
    if len(tied) == 1:
        return tied[0]
    return rng.choice(tied)


def plot_cell_heatmap_coloured(
    values: np.ndarray,
    row_labels: list[str],
    col_labels: list[str],
    out_png: str,
    out_svg: str,
    title: str,
    vmin: float = 0,
    vmax: float = 100,
    cmap_name: str = "viridis_r",
    nan_color: str = "#d9d9d9",
    figsize: tuple[float, float] | None = None,
    x_label: str = "Cells",
    y_label: str = "AB distribution",
    x_tick_fontsize: int = 6,
    y_tick_fontsize: int = 8,
):
    """
    Coloured heatmap for AB×Cell matrices where values are percentiles (0..100) and NaNs allowed.
    Matches the per-clone heatmap style: viridis_r, NaN in light grey.
    """
    masked = np.ma.masked_invalid(values)

    cmap = getattr(plt.cm, cmap_name).copy()
    cmap.set_bad(color=nan_color)

    if figsize is None:
        # heuristic sizing
        fig_w = max(6, 0.18 * len(col_labels))
        fig_h = max(4, 0.30 * len(row_labels))
        figsize = (fig_w, fig_h)

    fig, ax = plt.subplots(figsize=figsize, dpi=150)
    im = ax.imshow(masked, aspect="auto", interpolation="nearest", cmap=cmap, vmin=vmin, vmax=vmax)

    ax.set_yticks(np.arange(len(row_labels)))
    ax.set_yticklabels(row_labels, fontsize=y_tick_fontsize)

    ax.set_xticks(np.arange(len(col_labels)))
    ax.set_xticklabels(col_labels, rotation=90, fontsize=x_tick_fontsize)

    ax.set_xlabel(x_label, fontweight="bold")
    ax.set_ylabel(y_label, fontweight="bold")
    ax.set_title(title)

    # grid
    ax.set_xticks(np.arange(-0.5, len(col_labels), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(row_labels), 1), minor=True)
    ax.grid(which="minor", color="lightgrey", linestyle="-", linewidth=0.3)
    ax.tick_params(which="minor", bottom=False, left=False)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("KDE percentile (0–100)")

    fig.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)


def order_clones_by_similarity(matrix: pd.DataFrame,
                               metric: str = "euclidean",
                               method: str = "average") -> list[str]:
    """
    Return a list of clone (column) names ordered so that clones with similar
    AB profiles are adjacent, using hierarchical clustering.

    - matrix: AB (rows) × Clone (cols), numeric (0–100, NaNs allowed)
    """

    # Work on a copy of the underlying numeric data (Clone × AB)
    data = matrix.to_numpy(dtype=float).T  # shape: n_clones × n_ABs

    # Handle NaNs: replace with column means (per AB).
    # If an entire AB is NaN, np.nanmean -> NaN; replace those with 100.
    col_means = np.nanmean(data, axis=0)
    col_means = np.where(np.isnan(col_means), 100.0, col_means)
    nan_rows, nan_cols = np.where(np.isnan(data))
    data[nan_rows, nan_cols] = col_means[nan_cols]

    # Pairwise distances between clones
    dist_vec = pdist(data, metric=metric)

    # Hierarchical clustering
    Z = linkage(dist_vec, method=method)

    # Order of leaves in the dendrogram → column order
    leaf_order = leaves_list(Z)
    ordered_cols = matrix.columns[leaf_order].tolist()
    return ordered_cols

def clone_scores_by_overlap(matrix: pd.DataFrame, method: str = "mean") -> pd.Series:
    """
    Compute an overall-overlap score per clone (column) from the AB×Clone matrix.
    NaNs are ignored. Higher = more overlap.
    """
    method = method.lower()
    if method == "mean":
        scores = matrix.mean(axis=0, skipna=True)
    elif method == "median":
        scores = matrix.median(axis=0, skipna=True)
    elif method == "sum":
        scores = matrix.sum(axis=0, skipna=True)
    elif method == "max":
        scores = matrix.max(axis=0, skipna=True)
    else:
        raise ValueError(f"Unknown CLONE_OVERLAP_METHOD: {method}")
    # For clones that are all-NaN, push to the end
    scores = scores.fillna(-np.inf)
    return scores

def ab_to_cardinal_series(labels: list[str]) -> pd.Series:
    """Map each AB label to its cardinal class (first 3 characters)."""
    # assumes labels are strings with at least 3 chars
    return pd.Series({lab: str(lab)[:3] for lab in labels})

def min_by_clone_and_ab(df_num: pd.DataFrame, clone_col: str, ab_cols: list[str]) -> pd.DataFrame:
    """Return AB×Clone matrix of minimum percentile per clone"""
    grouped_min = df_num.groupby(clone_col, dropna=False)[ab_cols].min(numeric_only=True)
    return grouped_min.T  # AB rows × clone columns

def plot_lineage_histogram_from_groups(
    matrix: pd.DataFrame,
    mix_threshold: int,
    lineage_groups: list[list[str]],
    labels: list[str] | None,
    out_dir: str,
    base_name: str = "clone_lineage_category_hist",
    mixed_only: bool = True,
    min_abs_mixed: int = 2,            # ≥2 ABs present to consider “mixed” at AB level
    min_classes_mixed: int = 2         # ≥2 cardinal classes present to be “mixed” at class level
):
    """
    Classify each clone by whether its AB presence (AB min < mix_threshold) is a subset
    of exactly one of the provided lineage_groups (list of AB-label lists).
    If 0 or ≥2 groups match, classify as 'Other'.

    By default, only clones that are 'mixed' across cardinal classes (>= min_classes_mixed)
    are included in the histogram.

    Saves both PNG and SVG to out_dir/base_name.*
    Returns the counts DataFrame.
    """

    # 0/1 AB presence for each clone
    presence_ab = (matrix < mix_threshold).astype(bool).T  # Clone × AB

    # Keep only clones with ≥2 ABs present (still at AB level)
    presence_ab = presence_ab[presence_ab.sum(axis=1) >= min_abs_mixed]

    # Compute cardinal-class presence to decide “mixed clones”
    cardinal_of = ab_to_cardinal_series(matrix.index)            # AB -> first 3 chars
    class_presence = presence_ab.groupby(cardinal_of, axis=1).any()
    mixed_mask = class_presence.sum(axis=1) >= min_classes_mixed

    if mixed_only:
        presence_used = presence_ab.loc[mixed_mask]
    else:
        presence_used = presence_ab

    # Prepare lineage sets & labels
    lineage_sets = [set(group) for group in lineage_groups]
    if labels is None:
        labels = [f"Group {i+1}" for i in range(len(lineage_groups))]
    order = labels + ["Other"]

    def classify_clone_unique(ab_on: set[str]) -> str:
        """Return the single lineage label if ab_on ⊆ exactly one lineage set; else 'Other'."""
        if not ab_on:
            return "Other"
        matches = [lab for lab, s in zip(labels, lineage_sets) if ab_on.issubset(s)]
        return matches[0] if len(matches) == 1 else "Other"

    # Classify each clone
    rows = []
    for clone, row in presence_used.iterrows():
        present_abs = set(presence_used.columns[row.values.astype(bool)])
        rows.append((clone, classify_clone_unique(present_abs)))

    cat_df = pd.DataFrame(rows, columns=["Clone", "Category"])

    # Counts in fixed order
    counts = (
        cat_df["Category"]
             .value_counts()
             .reindex(order, fill_value=0)
             .rename_axis("Category")
             .reset_index(name="Frequency")
    )

    # Plot
    os.makedirs(out_dir, exist_ok=True)
    out_base = os.path.join(out_dir, base_name)
    fig, ax = plt.subplots(figsize=(max(3, 1.0*len(order)), 5), constrained_layout=True)
    sns.barplot(data=counts, x="Category", y="Frequency", order=order, edgecolor="black", ax=ax)
    ax.set_xlabel("Lineage category (subset match only)", fontweight="bold")
    ax.set_ylabel("Number of clones" + (" (mixed only)" if mixed_only else ""), fontweight="bold")
    ax.set_xticklabels(order, rotation=15, ha='right')
    fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_base + ".svg", bbox_inches="tight")
    plt.close(fig)

    return counts

def clone_to_animal(clone: str) -> str:
    """
    Extract animal ID from clone name.
    Assumes clone names like: S017_a1_c3_YYYYMMDD_HHMMSS
    → animal = S017
    """
    return str(clone).split("_")[0]


def main():
    ensure_output_dir(OUTPUT_DIR + '/Plotting_'+ str(curr_date))

    print(f"Reading: {INPUT_PATH}")
    df = read_table_auto(INPUT_PATH)

    # Identify clone column
    clone_col = "QueryFile" if "QueryFile" in df.columns else ("Clone" if "Clone" in df.columns else None)
    if clone_col is None:
        raise ValueError("Expected a 'QueryFile' (preferred) or 'Clone' column to identify clones.")

    if EXCLUDE_SINGLE_CELL_CLONES:
        # Count cells per clone. Prefer unique PointIndex if present.
        if "PointIndex" in df.columns:
            counts = df.groupby(clone_col)["PointIndex"].nunique(dropna=True)
        else:
            counts = df.groupby(clone_col).size()

        keep_clones = counts[counts > 1].index
        df = df[df[clone_col].isin(keep_clones)]
    
    # AB columns are strictly columns 9 to end by position
    if df.shape[1] < 6:
        raise ValueError("Not enough columns: need at least 5 metadata columns plus AB columns.")
    ab_cols = list(df.columns[9:])  # 0-based

    # Parse AB columns
    df_num = df.copy()
    for c in ab_cols:
        df_num[c] = df_num[c].apply(parse_ab_value)
    
    # ------------------------------------------------------------
    # Add a single predicted identity per cell (lowest percentile)
    # Random tie-breaking for equal minima
    # ------------------------------------------------------------
    rng_pred = np.random.default_rng(PREDICT_TIE_SEED)

    # Row-wise: choose AB label with minimum percentile (ties random)
    df_num[PREDICTED_ID_COLNAME] = df_num[ab_cols].apply(
        lambda r: choose_min_identity_random_tie(r, rng_pred),
        axis=1
    )

    # Save a per-cell CSV with predicted identity
    pred_csv = os.path.join(OUTPUT_DIR, "per_cell_predicted_identity.csv")
    df_num.to_csv(pred_csv, index=False)
    print(f"[SAVED] {pred_csv}")
    
    # Identify a per-cell ID column ONCE (used for global + per-clone plots)
    if "ID" in df_num.columns:
        id_col = "ID"
    elif "PointIndex" in df_num.columns:
        id_col = "PointIndex"
    else:
        id_col = None  # fall back to row index if necessary
    
    # Group by clone and take the minimum across cells per AB
    grouped_min = df_num.groupby(clone_col, dropna=False)[ab_cols].min(numeric_only=True)
    if PLOT_GLOBAL_CELL_TABLE:
        ### Per-cell table
        # Identify a per-cell ID column
        if "ID" in df_num.columns:
            id_col = "ID"
        elif "PointIndex" in df_num.columns:
            id_col = "PointIndex"
        else:
            id_col = None  # fall back to row index if necessary

        # Build human-readable cell labels: "Animal | ID"
        if id_col is not None:
            cell_labels = (
                df_num[clone_col].astype(str)
                + " | "
                + df_num[id_col].astype(str)
            )
        else:
            cell_labels = (
                df_num[clone_col].astype(str)
                + " | idx="
                + df_num.index.astype(str)
            )

        # Numeric matrix of AB values per cell
        cell_values = df_num[ab_cols].to_numpy(dtype=float)

        # Threshold into 0/1 (<= threshold == 1/black), keep NaNs as NaNs
        # Numeric matrix of AB values per cell
        cell_values = df_num[ab_cols].to_numpy(dtype=float)

        # --- NEW: memory-safe image for plotting (uint8 + mask) ---
        cell_present = (cell_values <= CELL_MEMBERSHIP_THRESHOLD) & ~np.isnan(cell_values)

        # uint8 image: 0/1 only (tiny memory)
        img = cell_present.astype(np.uint8)

        # mask NaNs so they become fully transparent
        mask = np.isnan(cell_values)

        # AB × Cell image for imshow
        masked_img = np.ma.array(img.T, mask=mask.T)

        # AB (rows) × Cell (cols)
        # Reconstruct a DataFrame version (0 / 1 / NaN) for CSVs and labels
        cell_matrix = pd.DataFrame(
            img.T.astype(float),   # 0/1
            index=ab_cols,
            columns=cell_labels
        )
        cell_matrix[mask.T] = np.nan

        # Optional: apply same AB ordering logic as main heatmap
        if AB_ORDER:
            listed_present = [a for a in AB_ORDER if a in cell_matrix.index]
            if APPEND_UNLISTED_ABS:
                unlisted = sorted([a for a in cell_matrix.index if a not in AB_ORDER])
                cell_row_order = listed_present + unlisted
            else:
                cell_row_order = listed_present
            cell_matrix = cell_matrix.loc[cell_row_order]
        else:
            cell_matrix = cell_matrix.sort_index(axis=0)

        
        # Clone label for each column (same order as original rows / columns)
        cell_clone_order = df_num[clone_col].astype(str).tolist()
        
        # Save the binary 0/1 table (NaN where no data / Outside if IGNORE_OUTSIDE=True)
        ensure_output_dir(OUTPUT_DIR)
        cell_table_csv_path = os.path.join(OUTPUT_DIR, CELL_TABLE_CSV)
        cell_matrix.to_csv(cell_table_csv_path)

        # Plot as black/white "graphical table"
        from matplotlib.colors import ListedColormap

        dpi = 100
        # scale figure size a bit with number of cells, but keep sane minimums
        fig_w = max(6.5, 0.20 * cell_matrix.shape[1])
        fig_h = max(4, 0.35 * cell_matrix.shape[0])

        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

        # alternating light background per clone block
        # two light colours that will alternate between contiguous clone segments
        clone_bg_colors = ["#e2e4e6", "#ffffff"]

        if len(cell_clone_order) == cell_matrix.shape[1]:
            current_clone = cell_clone_order[0]
            current_start = 0
            color_idx = 0

            # iterate plus a sentinel at the end to close the last span
            for col_idx, clone_name in enumerate(cell_clone_order + [None]):
                if col_idx == len(cell_clone_order) or clone_name != current_clone:
                    # shade columns [current_start .. col_idx-1]
                    ax.axvspan(
                        current_start - 0.5,
                        col_idx - 0.5,
                        facecolor=clone_bg_colors[color_idx % len(clone_bg_colors)],
                        alpha=1.0,
                        zorder=0,
                    )
                    current_start = col_idx
                    current_clone = clone_name
                    color_idx += 1
        # -----------------------------------------------------------

        masked_cells = _np.ma.masked_invalid(cell_matrix.to_numpy(dtype=float))

        # 0 -> semi-transparent white, 1 -> solid black
        bw_cmap = ListedColormap([
            (1.0, 1.0, 1.0, 0.25),  # value 0: white, 25% opaque
            (0.0, 0.0, 0.0, 1.0),   # value 1: black, 100% opaque
        ])
        # NaNs fully transparent (just show background band)
        bw_cmap.set_bad((0.0, 0.0, 0.0, 0.0))

        im_cells = ax.imshow(
            masked_img,
            aspect="auto",
            interpolation="nearest",
            resample=False,
            cmap=bw_cmap,
            vmin=0,
            vmax=1,
            zorder=1,
        )

        ax.set_yticks(_np.arange(cell_matrix.shape[0]))
        ax.set_yticklabels(cell_matrix.index, fontsize=8)

        ax.set_xticks(_np.arange(cell_matrix.shape[1]))
        ax.set_xticklabels(cell_matrix.columns, rotation=90, fontsize=6)

        ax.set_xlabel("Cell (QueryFile | ID)", fontweight="bold")
        ax.set_ylabel("AB distribution", fontweight="bold")
        ax.set_title(
            f"Cell × AB membership (black = percentile ≤ {CELL_MEMBERSHIP_THRESHOLD}%)",
            fontsize=10,
        )

        # Minor grid to emphasize boxes
        ax.set_xticks(_np.arange(-0.5, cell_matrix.shape[1], 1), minor=True)
        ax.set_yticks(_np.arange(-0.5, cell_matrix.shape[0], 1), minor=True)
        ax.grid(which="minor", color="lightgrey", linestyle="-", linewidth=0.3)
        ax.tick_params(which="minor", bottom=False, left=False)

        fig.tight_layout()
        cell_table_png_path = os.path.join(OUTPUT_DIR, CELL_TABLE_PNG)
        #fig.savefig(cell_table_png_path, dpi=300)
        #fig.savefig(os.path.join(OUTPUT_DIR, CELL_TABLE_SVG))
        #plt.show()
        plt.close(fig)
        # end of per-cell table
    
    ### Per-clone mini cell_AB_membership plots

    # Work clone-by-clone using df_num, so AB parsing is consistent
    all_clones = df_num[clone_col].astype(str).unique()

    for clone_name in all_clones:
        df_c = df_num[df_num[clone_col].astype(str) == clone_name].copy()
        if df_c.empty:
            continue

        # Build human-readable labels: "Clone | ID"
        if id_col is not None:
            col_labels_c = (
                df_c[clone_col].astype(str)
                + " | "
                + df_c[id_col].astype(str)
            )
        else:
            col_labels_c = (
                df_c[clone_col].astype(str)
                + " | idx="
                + df_c.index.astype(str)
            )

        # AB values for this clone only
        values_c = df_c[ab_cols].to_numpy(dtype=float)
        
        # Per-cell coloured heatmap
        if PLOT_PER_CLONE_CELL_TABLE_COLOURED:
            # Apply same AB ordering as global table
            row_order_c = ab_cols
            if AB_ORDER:
                listed_present = [a for a in AB_ORDER if a in ab_cols]
                if APPEND_UNLISTED_ABS:
                    unlisted = sorted([a for a in ab_cols if a not in AB_ORDER])
                    row_order_c = listed_present + unlisted
                else:
                    row_order_c = listed_present

            # Build AB×Cell percentiles matrix
            df_vals_c = pd.DataFrame(values_c, columns=ab_cols, index=col_labels_c)
            df_vals_c = df_vals_c[row_order_c]
            mat_c = df_vals_c.to_numpy(dtype=float).T  # AB×Cell

            clone_slug = sanitize_filename(clone_name)
            out_png = os.path.join(OUTPUT_DIR, f"cell_AB_membership_COLOURED_{clone_slug}.png")
            out_svg = os.path.join(OUTPUT_DIR, f"cell_AB_membership_COLOURED_{clone_slug}.svg")

            plot_cell_heatmap_coloured(
                values=mat_c,
                row_labels=row_order_c,
                col_labels=list(df_vals_c.index),
                out_png=out_png,
                out_svg=out_svg,
                title=f"{clone_name} – Cell × AB membership (COLOURED)",
                vmin=CELL_COLOUR_VMIN,
                vmax=CELL_COLOUR_VMAX,
                cmap_name="viridis_r",
                nan_color="#d9d9d9",
                x_label="Cell (QueryFile | ID)",
                y_label="AB distribution",
                x_tick_fontsize=6,
                y_tick_fontsize=8,
            )

        # Threshold into 0/1 (<= threshold == 1/black), keep NaNs as NaNs
        binary_c = (values_c <= CELL_MEMBERSHIP_THRESHOLD).astype(float)
        binary_c[np.isnan(values_c)] = np.nan

        # AB (rows) × Cell (cols) for this clone
        clone_cell_matrix = pd.DataFrame(
            binary_c.T,
            index=ab_cols,
            columns=col_labels_c
        )

        # Apply same AB ordering as global table
        if AB_ORDER:
            listed_present = [a for a in AB_ORDER if a in clone_cell_matrix.index]
            if APPEND_UNLISTED_ABS:
                unlisted = sorted([a for a in clone_cell_matrix.index if a not in AB_ORDER])
                row_order_c = listed_present + unlisted
            else:
                row_order_c = listed_present
            clone_cell_matrix = clone_cell_matrix.loc[row_order_c]
        else:
            clone_cell_matrix = clone_cell_matrix.sort_index(axis=0)

        # Plot for this clone
        dpi = 100
        # narrower figure: width depends on #cells in this clone
        fig_w_c = max(3, 0.25 * clone_cell_matrix.shape[1])
        fig_h_c = max(1.5, 0.35 * clone_cell_matrix.shape[0])

        fig_c, ax_c = plt.subplots(figsize=(fig_w_c, fig_h_c), dpi=dpi)

        # For a single clone we don't need alternating backgrounds, just a light base
        from matplotlib.colors import ListedColormap
        masked_c = np.ma.masked_invalid(clone_cell_matrix.to_numpy(dtype=float))

        bw_cmap_c = ListedColormap([
            (1.0, 1.0, 1.0, 0.25),  # 0 -> translucent white
            (0.0, 0.0, 0.0, 1.0),   # 1 -> solid black
        ])
        bw_cmap_c.set_bad((0.0, 0.0, 0.0, 0.0))  # NaNs transparent

        im_c = ax_c.imshow(
            masked_c,
            aspect="auto",
            interpolation="nearest",
            cmap=bw_cmap_c,
            vmin=0,
            vmax=1,
            zorder=1,
        )

        ax_c.set_yticks(np.arange(clone_cell_matrix.shape[0]))
        ax_c.set_yticklabels(clone_cell_matrix.index, fontsize=8)

        ax_c.set_xticks(np.arange(clone_cell_matrix.shape[1]))
        ax_c.set_xticklabels(clone_cell_matrix.columns, rotation=90, fontsize=6)

        ax_c.set_xlabel("Cell (QueryFile | ID)", fontweight="bold")
        ax_c.set_ylabel("AB distribution", fontweight="bold")
        ax_c.set_title(
            f"{clone_name} – cell × AB membership\n(black = percentile ≤ {CELL_MEMBERSHIP_THRESHOLD}%)",
            fontsize=10,
        )

        # Grid lines
        ax_c.set_xticks(np.arange(-0.5, clone_cell_matrix.shape[1], 1), minor=True)
        ax_c.set_yticks(np.arange(-0.5, clone_cell_matrix.shape[0], 1), minor=True)
        ax_c.grid(which="minor", color="lightgrey", linestyle="-", linewidth=0.3)
        ax_c.tick_params(which="minor", bottom=False, left=False)

        fig_c.tight_layout()

        clone_slug = sanitize_filename(clone_name)
        out_png = os.path.join(OUTPUT_DIR, f"cell_AB_membership_{clone_slug}.png")
        out_svg = os.path.join(OUTPUT_DIR, f"cell_AB_membership_{clone_slug}.svg")

        fig_c.savefig(out_png, bbox_inches="tight", dpi=300)
        fig_c.savefig(out_svg, bbox_inches="tight")
        plt.close(fig_c)
        # end per clone plots

    # Group by clone and take the minimum across cells per AB
    grouped_min = df_num.groupby(clone_col, dropna=False)[ab_cols].min(numeric_only=True)

    # Transpose to AB (rows) x Clone (cols)
    matrix = grouped_min.T

    # Apply simplified AB labels
    simplified_index = [simplify_ab_label(c) for c in matrix.index]
    matrix.index = simplified_index

    # Collapse duplicate AB labels by taking the minimum across duplicates
    if matrix.index.has_duplicates:
        matrix = matrix.groupby(matrix.index, sort=False).min()

    # === Apply custom AB order (rows) ===
    if AB_ORDER:
        listed_present = [a for a in AB_ORDER if a in matrix.index]
        if APPEND_UNLISTED_ABS:
            unlisted = sorted([a for a in matrix.index if a not in AB_ORDER])
            new_order_rows = listed_present + unlisted
        else:
            new_order_rows = listed_present
        matrix = matrix.loc[new_order_rows]
    else:
        matrix = matrix.sort_index(axis=0)

    # Sort clones in order of least to most overlapping ABs
    if SORT_CLONES == "by_overlap":
        scores = clone_scores_by_overlap(matrix, method=CLONE_OVERLAP_METHOD)
        # Stable sort (mergesort) so ties preserve original/alpha order
        ordered_cols = scores.sort_values(ascending=True, kind="mergesort").index.tolist()
        matrix = matrix[ordered_cols]
    elif SORT_CLONES == "by_similarity":
        # Use hierarchical clustering on clone profiles
        ordered_cols = order_clones_by_similarity(matrix,
                                                  metric="euclidean",
                                                  method="average")
        matrix = matrix[ordered_cols]
    elif SORT_CLONES is True:
        matrix = matrix.sort_index(axis=1)
    elif isinstance(SORT_CLONES, list):
        # Map from canonical (trimmed) clone name -> list of full column names
        canonical_to_full = {}
        for col in matrix.columns:
            can = canonical_clone_name(col)
            canonical_to_full.setdefault(can, []).append(col)

        # Build ordered list of full column names from CUSTOM_CLONE_ORDER
        ordered_full = []
        seen = set()
        for can in SORT_CLONES:
            if can in canonical_to_full:
                for full in canonical_to_full[can]:
                    if full not in seen:
                        ordered_full.append(full)
                        seen.add(full)

        # Any clones not mentioned in CUSTOM_CLONE_ORDER go at the end
        extras = [c for c in matrix.columns if c not in seen]
        matrix = matrix[ordered_full + extras]
    # else: leave as-is

        ###########################################################################
    # Predicted-identity presence heatmap (AB × Clone; black if any cell predicted)
    ###########################################################################
    # We use df_num (per-cell table), which already has PREDICTED_ID_COLNAME

    # Ensure AB row order matches your final AB order in `matrix`
    ab_order_final = list(matrix.index)

    # Ensure clone order matches `matrix.columns`
    clone_order_final = list(matrix.columns)

    # Build boolean presence: for each clone, which predicted identities exist?
    # Start with all False
    pred_presence = pd.DataFrame(False, index=ab_order_final, columns=clone_order_final)

    # For safety, filter to only clones we kept in analysis
    df_pred = df_num[df_num[clone_col].isin(clone_order_final)].copy()

    # Count presence (any cell in clone predicted as class)
    for clone, g in df_pred.groupby(clone_col, dropna=False):
        if clone not in pred_presence.columns:
            continue
        preds = g[PREDICTED_ID_COLNAME].dropna().astype(str).tolist()
        preds_set = set(preds)
        for ab in preds_set:
            if ab in pred_presence.index:
                pred_presence.loc[ab, clone] = True

    # Save as CSV (True/False)
    pred_presence_csv = os.path.join(OUTPUT_DIR, "predicted_identity_presence_ABxClone.csv")
    pred_presence.to_csv(pred_presence_csv)
    print(f"[SAVED] {pred_presence_csv}")

    # Plot: black boxes for True, transparent/white for False
    from matplotlib.colors import ListedColormap
    pred_img = pred_presence.to_numpy(dtype=int)  # 0/1

    bw_cmap = ListedColormap([
        (1, 1, 1, 0.0),  # 0 -> transparent
        (0, 0, 0, 1.0),  # 1 -> black
    ])

    fig_w = max(7, 0.22 * len(clone_order_final))
    fig_h = max(4, 0.35 * len(ab_order_final))

    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)
    ax.imshow(pred_img, aspect="auto", interpolation="nearest", cmap=bw_cmap, vmin=0, vmax=1)

    ax.set_xticks(np.arange(len(clone_order_final)))
    ax.set_yticks(np.arange(len(ab_order_final)))
    ax.set_xticklabels(clone_order_final, rotation=90, fontsize=7)
    ax.set_yticklabels(ab_order_final, fontsize=8)

    ax.set_title("Predicted identity presence per clone\n(black = ≥1 cell predicted as class)")
    ax.set_xlabel("Clone", fontweight="bold")
    ax.set_ylabel("Class / AB", fontweight="bold")

    ax.set_xticks(np.arange(-0.5, len(clone_order_final), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(ab_order_final), 1), minor=True)
    ax.grid(which="minor", color="lightgrey", linestyle="-", linewidth=0.3)
    ax.tick_params(which="minor", bottom=False, left=False)

    fig.tight_layout()
    out_base = os.path.join(OUTPUT_DIR, "predicted_identity_presence_heatmap")
    fig.savefig(out_base + ".png", dpi=300, bbox_inches="tight")
    fig.savefig(out_base + ".svg", bbox_inches="tight")
    plt.close(fig)

    print(f"[SAVED] {out_base}.png")
    print(f"[SAVED] {out_base}.svg")


    ### Calculate metrics of cardinal class mixing
    # matrix: AB(rows) × Clone(columns), values are min percentile (0..100), NaN allowed.
    AB_labels = matrix.index.tolist()
    cardinal_of = ab_to_cardinal_series(AB_labels)  # e.g., 'dI1', 'V1', ...

    # Degree of mixing per clone = number of 'distinct cardinal classes' that have value < threshold
    degree_mixing = {}
    present_mask = (matrix < MIX_THRESHOLD)  # AB×Clone boolean mask
    for clone in matrix.columns:
        # which ABs are present by threshold?
        present_ABs = [ab for ab, present in present_mask[clone].items() if bool(present)]
        # distinct cardinal classes among those ABs:
        degree_mixing[clone] = len(set(cardinal_of.loc[present_ABs].tolist()))

    degree_mixing = pd.Series(degree_mixing, name="DegreeMix")

    # "Mixed cardinal classes" flag: True if >1 cardinal class present
    is_mixed = degree_mixing > 1

    # Clone size = number of cells contributing to each clone (rows per clone in input)
    if "PointIndex" in df.columns:
        clone_size = df.groupby(clone_col)["PointIndex"].nunique(dropna=True)
    else:
        clone_size = df.groupby(clone_col).size()
    clone_size = clone_size.reindex(matrix.columns).fillna(0).astype(int)
    clone_size.name = "CloneSize"

    # Pack per-clone table
    per_clone = pd.DataFrame({
        "Clone": matrix.columns,
        "DegreeMix": degree_mixing.values,
        "IsMixed": is_mixed.values,
        "CloneSize": clone_size.values
    }).set_index("Clone")

    # Save files
    matrix_csv_path = os.path.join(OUTPUT_DIR, MATRIX_CSV)
    matrix.to_csv(matrix_csv_path)
    per_clone_csv_path = os.path.join(OUTPUT_DIR, MIXING_CSV)
    per_clone.to_csv(per_clone_csv_path)
    
    ###########################################################################
    # Per-clone heatmap data + D/V mixing flag
    ###########################################################################
    # Clone × AB matrix of minimum percentiles (same values as heatmap)
    clone_heatmap = matrix.T.copy()   # rows = clones, cols = ABs

    # Presence mask at MIX_THRESHOLD (AB × Clone)
    presence = (matrix < MIX_THRESHOLD)  # True = clone reaches that AB at < threshold

    # Define dorsal/ventral groups from DV_GROUPS:
    #   DV_GROUPS[0] = "other cardinal classes"
    #   DV_GROUPS[1] = "dI5" class
    other_ABs = set(DV_GROUPS[0])
    dI5_ABs   = set(DV_GROUPS[1])

    dv_mixing_flag = {}
    for clone in matrix.columns:
        # Does this clone have any dI5 AB "present"?
        has_dI5 = any(
            (ab in presence.index) and bool(presence.loc[ab, clone])
            for ab in dI5_ABs
        )
        # Does this clone have any non-dI5 (other) AB "present"?
        has_other = any(
            (ab in presence.index) and bool(presence.loc[ab, clone])
            for ab in other_ABs
        )
        dv_mixing_flag[clone] = "yes" if (has_dI5 and has_other) else "no"

    # Add D/V mixing column to the per-clone heatmap table
    clone_heatmap.insert(
        0,
        "D/V mixing",
        [dv_mixing_flag[c] for c in clone_heatmap.index]
    )

    # (Optional) also include clone size if you like:
    # clone_heatmap.insert(
    #     0,
    #     "CloneSize",
    #     clone_size.reindex(clone_heatmap.index).astype(int).tolist()
    # )

    # Save CSV: one row per clone, AB percentile columns, plus D/V mixing flag
    clone_heatmap_csv = os.path.join(
        OUTPUT_DIR,
        "query_cells_percentile_membership_clone_heatmap_DVmix.csv"
    )
    clone_heatmap.to_csv(clone_heatmap_csv)
    print(f"[SAVED] {clone_heatmap_csv}")
    
    ##############################################################################################################################################################################################################
    ### Plot heat map with inverted viridis (yellow = small values)
    ##############################################################################################################################################################################################################
    data = matrix.to_numpy(dtype=float)

    dpi = 100
    fig_w = max(FIG_SIZE_PX[0] / dpi, 7)
    fig_h = max(FIG_SIZE_PX[1] / dpi, 4)

    masked = np.ma.masked_invalid(data)
    cmap = plt.cm.viridis_r  # inverted: yellow ~ low, purple ~ high
    cmap = cmap.copy()
    cmap.set_bad(color="#d9d9d9")  # light gray for NaN

    # --- Clone sizes (used in both plots) ---
    if "PointIndex" in df.columns:
        clone_counts = df.groupby(clone_col)["PointIndex"].nunique(dropna=True)
    else:
        clone_counts = df.groupby(clone_col).size()
    clone_counts = clone_counts.reindex(matrix.columns).fillna(0).astype(int)

    # --- TM per clone (used for colouring tick labels in 2nd plot) ---
    if "TM" in df.columns:
        tm_per_clone = df.groupby(clone_col)["TM"].first().reindex(matrix.columns)
    else:
        tm_per_clone = pd.Series(index=matrix.columns, dtype=object)

    # Map TM -> colour (assumes ColPals is accessible via cap.ColPals)
    TM_COLOURS = {
        "E9.5":  cap.ColPals.light_orange,
        "E10.5": cap.ColPals.med_orange,
        "E11.5": cap.ColPals.dark_orange,
    }

    # ============================================================
    # 1) ORIGINAL HEATMAP (full clone label + n= on x axis)
    # ============================================================
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)
    im = ax.imshow(
        masked,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=100,
    )

    ax.set_xticks(np.arange(matrix.shape[1]))
    ax.set_yticks(np.arange(matrix.shape[0]))
    ax.set_yticklabels(matrix.index, fontsize=9)

    ax.set_xlabel("Clone (QueryFile)")
    ax.set_ylabel("AB distribution")
    ax.set_title("Minimum KDE Percentile per Clone × AB (yellow = smaller)")

    # Minor grid
    ax.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax.grid(which="minor", color="white", linestyle="-", linewidth=0.5)
    ax.tick_params(which="minor", bottom=False, left=False)

    # --- multi-line labels: "CloneName (n=###)" ---
    xticklabels = [f"{c} (n={clone_counts.get(c, 0)})" for c in matrix.columns]
    ax.set_xticklabels(xticklabels, rotation=90, fontsize=9)

    if ANNOTATE_VALUES and matrix.size <= 3600:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix.iat[i, j]
                if not np.isnan(val):
                    ax.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=7)

    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Minimum KDE percentile (0–100)")

    fig.tight_layout()
    heatmap_path = os.path.join(OUTPUT_DIR, HEATMAP_PNG)
    fig.savefig(heatmap_path, bbox_inches="tight")
    fig.savefig(str(heatmap_path) + ".svg", format="svg", dpi=300)
    plt.close(fig)

    # ============================================================
    # 2) SECOND HEATMAP – x axis shows ONLY clone size "n = x"
    #    with a TM-coloured strip directly under the heatmap
    # ============================================================
    import matplotlib.patches as mpatches
    from matplotlib.transforms import blended_transform_factory

    fig2, ax2 = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi, constrained_layout=True)

    im2 = ax2.imshow(
        masked,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=100,
    )

    ax2.set_xticks(np.arange(matrix.shape[1]))
    ax2.set_yticks(np.arange(matrix.shape[0]))
    ax2.set_yticklabels(matrix.index, fontsize=9)

    ax2.set_xlabel("Clone size (cells)")
    ax2.set_ylabel("AB distribution")
    ax2.set_title(
        "Minimum KDE Percentile per Clone × AB\n"
        "(x-axis = clone size, TM strip under heatmap)"
    )

    # Minor grid
    ax2.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax2.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax2.grid(which="minor", color="white", linestyle="-", linewidth=0.5)
    ax2.tick_params(which="minor", bottom=False, left=False)

    # --- labels: ONLY "n = x", in black ---
    n_only_labels = [f"{clone_counts.get(c, 0)}" for c in matrix.columns]
    ax2.set_xticklabels(n_only_labels, rotation=90, fontsize=9, color="black")

    if ANNOTATE_VALUES and matrix.size <= 3600:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix.iat[i, j]
                if not np.isnan(val):
                    ax2.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=7)

    # Make room at the bottom for strip + labels
    fig2.subplots_adjust(bottom=0.30)

    # ============================================================
    # 2) SECOND HEATMAP – TM strip + clone-size-bin strip
    # ============================================================
    import matplotlib.patches as mpatches
    from matplotlib.transforms import blended_transform_factory

    fig2, ax2 = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    im2 = ax2.imshow(
        masked,
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=100,
    )

    ax2.set_xticks(np.arange(matrix.shape[1]))
    ax2.set_yticks(np.arange(matrix.shape[0]))
    ax2.set_yticklabels(matrix.index, fontsize=9)

    ax2.set_xlabel("Clones")
    ax2.set_ylabel("AB distribution")
    ax2.set_title(
        "Minimum KDE Percentile per Clone × AB\n"
        "(TM strip + clone-size-bin strip under heatmap)"
    )

    # Minor grid
    ax2.set_xticks(np.arange(-0.5, matrix.shape[1], 1), minor=True)
    ax2.set_yticks(np.arange(-0.5, matrix.shape[0], 1), minor=True)
    ax2.grid(which="minor", color="white", linestyle="-", linewidth=0.5)
    ax2.tick_params(which="minor", bottom=False, left=False)

    # --- Option A: remove x tick labels entirely (recommended if you use strips) ---
    ax2.set_xticklabels([""] * matrix.shape[1])

    # (If you prefer keeping the clone sizes as labels, comment the line above and use this:)
    # n_only_labels = [f"{clone_counts.get(c, 0)}" for c in matrix.columns]
    # ax2.set_xticklabels(n_only_labels, rotation=90, fontsize=9, color="black")

    if ANNOTATE_VALUES and matrix.size <= 3600:
        for i in range(matrix.shape[0]):
            for j in range(matrix.shape[1]):
                val = matrix.iat[i, j]
                if not np.isnan(val):
                    ax2.text(j, i, f"{val:.0f}", ha="center", va="center", fontsize=7)

    # Make room at the bottom for BOTH strips
    fig2.subplots_adjust(bottom=0.32)

    # ----------------------------------------------------------------
    # Two colour strips: TM strip (upper) + size-bin strip (lower)
    # ----------------------------------------------------------------
    m = matrix.shape[0]
    cell_height_ax = 1.0 / max(1, m)

    # strip heights (> half a heatmap box)
    strip_height_ax = 0.5 * cell_height_ax

    # vertical positions in axes coords (negative = below heatmap)
    # TM strip sits closer to the heatmap, size strip below it
    tm_y0   = -strip_height_ax
    size_y0 = -2.05 * strip_height_ax  # gap-ish between them

    trans = blended_transform_factory(ax2.transData, ax2.transAxes)

    # --- TM colours (as before) ---
    TM_COLOURS = {
        "E9.5":  cap.ColPals.light_orange,
        "E10.5": cap.ColPals.med_orange,
        "E11.5": cap.ColPals.dark_orange,
    }

    # --- Clone size binning + greys ---
    # bins: 0-5, 6-10, 11-20, 21-40, 41-70
    def size_bin(n: int) -> int:
        if n <= 5:
            return 0
        elif n <= 10:
            return 1
        elif n <= 20:
            return 2
        elif n <= 40:
            return 3
        else:
            return 4  # 41-70 (and anything above)

    SIZE_BIN_LABELS = ["0–5", "6–10", "11–20", "21–40", "41–70"]
    SIZE_GREYS = ["#f0f0f0", "#d9d9d9", "#bdbdbd", "#969696", "#636363"]  # light -> dark

    # Draw strips column-by-column
    for j, clone in enumerate(matrix.columns):
        # --- TM strip ---
        tm_val = str(tm_per_clone.get(clone, ""))
        tm_colour = TM_COLOURS.get(tm_val, "lightgrey")

        rect_tm = mpatches.Rectangle(
            (j - 0.5, tm_y0),
            1.0,
            strip_height_ax,
            transform=trans,
            clip_on=False,
            facecolor=tm_colour,
            edgecolor="none",
        )
        ax2.add_patch(rect_tm)

        # --- size strip ---
        n = int(clone_counts.get(clone, 0))
        b = size_bin(n)
        size_colour = SIZE_GREYS[b]

        rect_sz = mpatches.Rectangle(
            (j - 0.5, size_y0),
            1.0,
            strip_height_ax,
            transform=trans,
            clip_on=False,
            facecolor=size_colour,
            edgecolor="none",
        )
        ax2.add_patch(rect_sz)

    # Optional: legends for the strips
    tm_handles = [
        mpatches.Patch(color=TM_COLOURS["E9.5"],  label="E9.5"),
        mpatches.Patch(color=TM_COLOURS["E10.5"], label="E10.5"),
        mpatches.Patch(color=TM_COLOURS["E11.5"], label="E11.5"),
    ]
    sz_handles = [
        mpatches.Patch(color=SIZE_GREYS[i], label=SIZE_BIN_LABELS[i])
        for i in range(len(SIZE_BIN_LABELS))
    ]

    ax2.legend(
        handles=tm_handles + sz_handles,
        frameon=False,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        title="TM / Clone size bins"
    )

    # Colourbar for the heatmap
    cbar2 = fig2.colorbar(im2, ax=ax2, fraction=0.046, pad=0.04)
    cbar2.set_label("Minimum KDE percentile (0–100)")

    heatmap_path_n = os.path.join(OUTPUT_DIR, HEATMAP_N_ONLY_PNG)
    fig2.savefig(heatmap_path_n)
    fig2.savefig(str(heatmap_path_n) + ".svg", format="svg", dpi=300)
    plt.close(fig2)
    
    ##############################################################################################################################################################################################################
    ### Proportion of clones with mixed cardinal classes
    # Stacked bar: proportion Mixed vs Not mixed
    stack_df = (
        per_clone["IsMixed"]
        .value_counts()
        .rename_axis("MixedFlag")
        .reset_index(name="Count")
    )

    # Convert to proportions and to a “long” frame compatible with snsStackedBarPlot
    total = stack_df["Count"].sum()
    plot_df = pd.DataFrame({
        "X": ["All clones", "All clones"],             # one category
        "Status": ["Not mixed", "Mixed"],
        "Weight": [
            float(stack_df.loc[stack_df["MixedFlag"] == False, "Count"].sum()) / total if total else 0.0,
            float(stack_df.loc[stack_df["MixedFlag"] == True, "Count"].sum()) / total if total else 0.0
        ]
    })

    stack_base = os.path.join(OUTPUT_DIR, "proportion_mixed_cardinal_classes")

    # Use your helper (saves .svg)
    cap.snsStackedBarPlot(
        data=plot_df,
        x="X",
        hue="Status",
        weights="Weight",
        x_axis_title="",
        y_axis_title="Proportion of clones",
        x_axis_ticks=["All clones"],
        output_path=stack_base,
        palette=("lightgrey","dimgray"),
        fig_size=(3.5,5),
        hue_order=["Not mixed", "Mixed"]
    )

    # Also export PNG (replot quickly using matplotlib so you get both formats)
    import matplotlib.pyplot as _plt
    _fig, _ax = _plt.subplots(figsize=(3.5,5))
    bottom = 0.0
    for label, val in zip(plot_df["Status"], plot_df["Weight"]):
        _ax.bar([0], [val], bottom=bottom, label=label)
        bottom += val
    _ax.set_xticks([0]); _ax.set_xticklabels(["All clones"], fontweight="bold")
    _ax.set_ylabel("Proportion of clones", fontweight="bold")
    _ax.legend(frameon=False)
    _fig.tight_layout()
    _fig.savefig(stack_base + ".png", dpi=300, bbox_inches="tight")
    _plt.close(_fig)

    ################################################################################################
    ### dI5 vs other clone-category plots (ALL-TM one bar + per-TM bars) + saved tables
    ################################################################################################

    tm_order  = ["E9.5", "E10.5", "E11.5"]
    cat_order = ["dI5 + other", "dI5-only", "other-only"]

    # ---------------------------
    # 1) Clone -> TM mapping (from THIS script's input table)
    # ---------------------------
    if "TM" not in df_num.columns:
        raise ValueError("TM column not found in input table; cannot make per-TM plots.")

    clone_tm = (
        df_num[[clone_col, "TM"]]
        .copy()
    )
    clone_tm[clone_col] = clone_tm[clone_col].astype(str).str.strip()
    clone_tm["TM"]      = clone_tm["TM"].astype(str).str.strip()

    # One TM per clone (take first)
    clone_tm = clone_tm.dropna(subset=[clone_col, "TM"]).drop_duplicates(subset=[clone_col])
    clone_tm = clone_tm[clone_tm["TM"].isin(tm_order)].copy()
    clone_tm["TM"] = pd.Categorical(clone_tm["TM"], categories=tm_order, ordered=True)

    # ---------------------------
    # 2) Determine clone categories using *CELL-LEVEL EXCLUSIVE* membership
    # ---------------------------

    # Identify di5 vs non-di5 AB columns from the per-cell table
    dI5_cols = [c for c in ab_cols if str(c).startswith("dI5")]
    other_cols = [c for c in ab_cols if c not in dI5_cols]

    if len(dI5_cols) == 0:
        raise ValueError("No dI5 AB columns found (expected columns starting with 'dI5').")

    df_cat = df_num[[clone_col, "TM"] + (["PointIndex"] if "PointIndex" in df_num.columns else []) + ab_cols].copy()
    df_cat[clone_col] = df_cat[clone_col].astype(str).str.strip()
    df_cat["TM"] = df_cat["TM"].astype(str).str.strip()

    # Cell-level membership flags (≤ threshold anywhere within group)
    df_cat["Has_dI5"] = (df_cat[dI5_cols].min(axis=1, skipna=True) <= MIX_THRESHOLD)
    df_cat["Has_other"] = (df_cat[other_cols].min(axis=1, skipna=True) <= MIX_THRESHOLD)

    # Define a stable per-cell key (prefer PointIndex)
    if "PointIndex" in df_cat.columns:
        cell_key = "PointIndex"
    else:
        df_cat = df_cat.reset_index().rename(columns={"index": "RowIndex"})
        cell_key = "RowIndex"

    # Collapse to ONE row per cell within clone (in case of duplicates)
    cell_level = (
        df_cat
        .groupby([clone_col, cell_key], dropna=False)
        .agg(
            TM=("TM", "first"),
            Has_dI5=("Has_dI5", "max"),
            Has_other=("Has_other", "max"),
        )
        .reset_index()
    )

    # Exclusive cell types
    cell_level["dI5_only_cell"] = cell_level["Has_dI5"] & (~cell_level["Has_other"])
    cell_level["other_only_cell"] = cell_level["Has_other"] & (~cell_level["Has_dI5"])
    cell_level["mixed_cell"] = cell_level["Has_dI5"] & (cell_level["Has_other"])

    # Per-clone counts of each cell type
    clone_summary = (
        cell_level
        .groupby(clone_col, dropna=False)
        .agg(
            TM=("TM", "first"),
            n_dI5_only=("dI5_only_cell", "sum"),
            n_other_only=("other_only_cell", "sum"),
            n_mixed=("mixed_cell", "sum"),
            n_any=("Has_dI5", "sum"),  # not used for logic; just a sanity count
        )
        .reset_index()
        .rename(columns={clone_col: "Clone"})
    )

    # Clone category logic (STRICT "only" bins)
    category = []
    for _, r in clone_summary.iterrows():
        n_di5 = int(r["n_dI5_only"])
        n_oth = int(r["n_other_only"])
        n_mix = int(r["n_mixed"])

        if (n_di5 > 0) and (n_oth == 0) and (n_mix == 0):
            cat = "dI5-only"
        elif (n_oth > 0) and (n_di5 == 0) and (n_mix == 0):
            cat = "other-only"
        elif (n_di5 > 0) and (n_oth > 0):
            # has both exclusive cell types
            cat = "dI5 + other"
        else:
            # e.g. only mixed cells, or dI5_only + mixed, or other_only + mixed, or nothing
            cat = "neither"
        category.append(cat)

    clone_summary["Category"] = category

    # Keep only valid TMs and set order (same as your plotting code)
    tm_order  = ["E9.5", "E10.5", "E11.5"]
    clone_summary = clone_summary[clone_summary["TM"].isin(tm_order)].copy()
    clone_summary["TM"] = pd.Categorical(clone_summary["TM"], categories=tm_order, ordered=True)

    # This is what the stacked-bar plots downstream should use:
    per_clone_cat_df = clone_summary.copy()
    per_clone_cat_df["Clone"] = per_clone_cat_df["Clone"].astype(str).str.strip()

    # Save a diagnostic table so you can inspect “neither” / mixed cases
    dI5_base = os.path.join(OUTPUT_DIR, "proportion_dI5_categories")
    per_clone_cat_df.to_csv(dI5_base + "_per_clone_category_with_TM_diagnostic.csv", index=False)

    # merge TM (clone_col may be "QueryFile" or "Clone")
    clone_tm_renamed = clone_tm.rename(columns={clone_col: "Clone"})
    per_clone_cat_df = per_clone_cat_df.merge(clone_tm_renamed, on="Clone", how="left")


    dI5_base = os.path.join(OUTPUT_DIR, "proportion_dI5_categories")

    # Save per-clone category table (with TM)
    per_clone_cat_df.to_csv(dI5_base + "_per_clone_category_with_TM.csv", index=False)

    # ---------------------------
    # Helper to make counts/proportions tables
    # ---------------------------
    def make_stacked_tables(df_in, x_col, drop_neither=True):
        d = df_in.copy()
        if drop_neither:
            d = d[d["Category"] != "neither"].copy()

        # x order
        if hasattr(d[x_col], "cat"):
            x_levels = list(d[x_col].cat.categories)
        else:
            x_levels = sorted(d[x_col].dropna().unique())

        counts = (
            d.groupby([x_col, "Category"])
            .size()
            .reindex(
                pd.MultiIndex.from_product([x_levels, cat_order], names=[x_col, "Category"]),
                fill_value=0
            )
            .reset_index(name="Count")
        )

        totals = counts.groupby(x_col)["Count"].transform("sum")
        counts["Total"] = totals
        counts["Weight"] = np.where(totals.to_numpy() > 0, counts["Count"] / totals, 0.0)

        # --- NEW: binomial SE of the proportion (always >= 0) ---
        counts["SE"] = np.sqrt(np.where(
            totals.to_numpy() > 0,
            counts["Weight"] * (1.0 - counts["Weight"]) / totals,
            0.0
        )).astype(float)

        # long df for plotting
        plot_df = (
            counts.rename(columns={x_col: "X", "Category": "Status"})
                [["X", "Status", "Weight", "SE"]]
        )
        return plot_df, counts

    # ---------------------------
    # 3) ONE-BAR plot (all TMs pooled)
    # ---------------------------
    per_clone_all = per_clone_cat_df.copy()
    per_clone_all["All"] = "All clones"

    plot_df_all, counts_all = make_stacked_tables(per_clone_all, x_col="All", drop_neither=True)
    dI5_all_base = dI5_base + "_ALLTM"

    # ---- Add clone + animal counts (ALL-TM) ----
    n_clones_all = (
        per_clone_cat_df[per_clone_cat_df["Category"] != "neither"]["Clone"]
        .nunique()
    )

    n_animals_all = (
        per_clone_cat_df[per_clone_cat_df["Category"] != "neither"]["Clone"]
        .apply(clone_to_animal)
        .nunique()
    )

    counts_all["NumClones"]  = n_clones_all
    counts_all["NumAnimals"] = n_animals_all
    
    cap.snsStackedBarPlot(
        data=plot_df_all,
        x="X",
        hue="Status",
        weights="Weight",
        x_axis_title="",
        y_axis_title=f"Proportion of clones (≤ {MIX_THRESHOLD}%)",
        x_axis_ticks=["All clones"],
        output_path=dI5_all_base,
        palette=("lightgrey", "darkgrey", "dimgray"),
        fig_size=(0.6, 5),
        hue_order=cat_order,
    )


    counts_all.to_csv(dI5_all_base + "_counts.csv", index=False)
    plot_df_all.to_csv(dI5_all_base + "_proportions.csv", index=False)

    # ---------------------------
    # 4) PER-TM plot (three bars)
    # ---------------------------
    plot_df_tm, counts_tm = make_stacked_tables(per_clone_cat_df, x_col="TM_x", drop_neither=True)
    dI5_tm_base = dI5_base + "_by_TM"
    
    # ---- Add clone + animal counts PER TM ----
    clone_counts_tm = (
        per_clone_cat_df[per_clone_cat_df["Category"] != "neither"]
        .groupby("TM_x")["Clone"]
        .nunique()
    )

    animal_counts_tm = (
        per_clone_cat_df[per_clone_cat_df["Category"] != "neither"]
        .assign(Animal=lambda d: d["Clone"].apply(clone_to_animal))
        .groupby("TM_x")["Animal"]
        .nunique()
    )

    counts_tm["NumClones"] = counts_tm["TM_x"].map(clone_counts_tm)
    counts_tm["NumAnimals"] = counts_tm["TM_x"].map(animal_counts_tm)

    cap.snsStackedBarPlot(
        data=plot_df_tm,
        x="X",
        hue="Status",
        weights="Weight",
        x_axis_title="TM",
        y_axis_title=f"Proportion of clones (≤ {MIX_THRESHOLD}%)",
        x_axis_ticks=tm_order,
        output_path=dI5_tm_base,
        palette=("lightgrey", "darkgrey", "dimgray"),
        fig_size=(4.5, 5),
        hue_order=cat_order,
    )

    counts_tm.to_csv(dI5_tm_base + "_counts.csv", index=False)
    plot_df_tm.to_csv(dI5_tm_base + "_proportions.csv", index=False)

    ################################################################################################
    ### Purity plot for mixed clones (dI5 + other): per-clone fraction of cells that are dI5 vs other
    ################################################################################################

    # --- Identify mixed clones from your existing per-clone category table ---
    mixed_clones = (
        per_clone_cat_df.loc[per_clone_cat_df["Category"] == "dI5 + other", "Clone"]
        .astype(str)
        .tolist()
    )

    purity_out_base = os.path.join(OUTPUT_DIR, f"purity_mixed_dI5_other_per_clone_thr{MIX_THRESHOLD}")

    if len(mixed_clones) == 0:
        print("[purity] No mixed (dI5 + other) clones found; skipping purity plot.")
    else:
        # --- Build a per-cell component label based on predicted identity ---
        df_purity = df_num[df_num[clone_col].astype(str).isin(mixed_clones)].copy()

        if PREDICTED_ID_COLNAME not in df_purity.columns:
            raise ValueError(
                f"[purity] Expected '{PREDICTED_ID_COLNAME}' in df_num. "
                "Make sure predicted identity is computed before this block."
            )

        # Predicted identity is an AB label; first 3 chars = cardinal class
        df_purity["PredCardinal"] = df_purity[PREDICTED_ID_COLNAME].astype(str).str[:3]
        df_purity["Component"] = np.where(df_purity["PredCardinal"] == "dI5", "dI5", "other")

        # Count CELLS per clone × component (use unique PointIndex if available)
        if "PointIndex" in df_purity.columns:
            counts = (
                df_purity.dropna(subset=["PointIndex"])
                         .groupby([clone_col, "Component"])["PointIndex"]
                         .nunique()
                         .rename("CellCount")
                         .reset_index()
                         .rename(columns={clone_col: "Clone"})
            )
        else:
            counts = (
                df_purity.groupby([clone_col, "Component"])
                         .size()
                         .rename("CellCount")
                         .reset_index()
                         .rename(columns={clone_col: "Clone"})
            )

        # Wide per-clone table (ensure both columns exist)
        wide = (
            counts.pivot(index="Clone", columns="Component", values="CellCount")
                  .fillna(0)
                  .reset_index()
        )
        if "dI5" not in wide.columns:
            wide["dI5"] = 0
        if "other" not in wide.columns:
            wide["other"] = 0

        wide["TotalCells"] = wide["dI5"] + wide["other"]

        # Drop any weird clones with 0 counted cells
        wide = wide[wide["TotalCells"] > 0].copy()

        # Fractions (each clone sums to 1.0)
        wide["Frac_dI5"]   = wide["dI5"] / wide["TotalCells"]
        wide["Frac_other"] = wide["other"] / wide["TotalCells"]

        # Sort for readability
        wide = wide.sort_values(["Frac_dI5", "Clone"], ascending=[False, True]).reset_index(drop=True)

        # ---------------------------
        # DIRECT stacked bar plot
        # ---------------------------
        x = np.arange(len(wide))
        frac_dI5 = wide["Frac_dI5"].to_numpy(dtype=float)
        frac_other = wide["Frac_other"].to_numpy(dtype=float)

        dI5_n = wide["dI5"].to_numpy(dtype=int)
        other_n = wide["other"].to_numpy(dtype=int)

        fig_w = max(10.0, 0.35 * len(wide))  # long clone IDs need more width
        fig_h = 5

        fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=150)

        ax.bar(x, frac_dI5, label="dI5")                         # bottom = 0
        ax.bar(x, frac_other, bottom=frac_dI5, label="other")    # stacked on top

        ax.set_ylim(0, 1.0)
        ax.set_ylabel("Fraction of cells (predicted identity)", fontweight="bold")
        ax.set_xlabel("Mixed clones (sorted by dI5 fraction)", fontweight="bold")
        ax.set_title("Purity of mixed clones: dI5 vs other")

        # Use clone IDs as tick labels
        ax.set_xticks(x)
        ax.set_xticklabels(wide["Clone"].tolist(), rotation=90, fontsize=7)

        ax.legend(frameon=False)

        # ---------------------------
        # Annotate counts inside bars
        # ---------------------------
        # Only draw text if the segment is tall enough (avoid unreadable clutter)
        min_height_for_text = 0.06  # ~6% of bar height

        for i in range(len(wide)):
            # dI5 segment annotation
            if frac_dI5[i] >= min_height_for_text and dI5_n[i] > 0:
                ax.text(
                    x[i],
                    frac_dI5[i] / 2.0,
                    str(dI5_n[i]),
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="white",
                    fontweight="bold"
                )
            elif dI5_n[i] > 0:
                # if too small, place just above the segment
                ax.text(
                    x[i],
                    frac_dI5[i] + 0.01,
                    str(dI5_n[i]),
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    color="black"
                )

            # other segment annotation
            if frac_other[i] >= min_height_for_text and other_n[i] > 0:
                ax.text(
                    x[i],
                    frac_dI5[i] + frac_other[i] / 2.0,
                    str(other_n[i]),
                    ha="center",
                    va="center",
                    fontsize=8,
                    color="black",
                    fontweight="bold"
                )
            elif other_n[i] > 0:
                # if too small, place near the very top
                ax.text(
                    x[i],
                    min(0.99, frac_dI5[i] + frac_other[i] + 0.01),
                    str(other_n[i]),
                    ha="center",
                    va="bottom",
                    fontsize=7,
                    color="black"
                )

        fig.tight_layout()
        fig.savefig(purity_out_base + ".png", dpi=300, bbox_inches="tight")
        fig.savefig(purity_out_base + ".svg", bbox_inches="tight")
        plt.close(fig)

        # ---------------------------
        # Save tables
        # ---------------------------
        wide.to_csv(purity_out_base + "_per_clone_purity_table.csv", index=False)

        plot_df_purity = pd.DataFrame({
            "Clone": np.repeat(wide["Clone"].values, 2),
            "Status": ["dI5", "other"] * len(wide),
            "CellCount": np.concatenate([wide["dI5"].values, wide["other"].values]),
            "Fraction": np.concatenate([wide["Frac_dI5"].values, wide["Frac_other"].values]),
        })
        plot_df_purity.to_csv(purity_out_base + "_plot_long.csv", index=False)
    
    ##############################################################################################################################################################################################################
    ### Histogram of AB combinations across clones (only plots clones with mixed cardinal classes) (thresholded presence)
    plot_lineage_histogram_from_groups(
    matrix=matrix,
    mix_threshold=MIX_THRESHOLD,
    lineage_groups=LINEAGE_COMP,                 # your list of AB-label lists
    labels=["dI1/dI2/dI3", "dI5", "dI6/V0/V1", "V2/MN", "MN/V3"],  # or None for Group 1..N
    out_dir=OUTPUT_DIR,
    base_name="lineage_subdivision_hist",
    mixed_only=True,              # only clones mixed across cardinal classes
    min_abs_mixed=2,              # require ≥2 ABs present
    min_classes_mixed=2           # require ≥2 cardinal classes present
    )

    plot_lineage_histogram_from_groups(
    matrix=matrix,
    mix_threshold=MIX_THRESHOLD,
    lineage_groups=DV_GROUPS,                 # your list of AB-label lists
    labels=["dI1/dI2/dI3/dI6/V0/V1/V2/MN/V3", "dI5"],  # or None for Group 1..N
    out_dir=OUTPUT_DIR,
    base_name="dv_groups_hist",
    mixed_only=True,              # only clones mixed across cardinal classes
    min_abs_mixed=2,              # require ≥2 ABs present
    min_classes_mixed=2           # require ≥2 cardinal classes present
    )
    ##############################################################################################################################################################################################################
    ### Scatter plot of clone size vs degree of cardinal class mixing
    # Scatter: clone size vs degree of mixing
    scatter_df = per_clone.reset_index()  # columns: Clone, DegreeMix, CloneSize
    scat_base = os.path.join(OUTPUT_DIR, "clone_size_vs_degree_mixing")

    # Use your snsScatterPlot for .svg
    cap.snsScatterPlot(
        data=scatter_df,
        x="CloneSize",
        y="DegreeMix",
        x_axis_title="Clone size (cells)",
        y_axis_title=f"Degree of mixing (# cardinal classes < {MIX_THRESHOLD}%)",
        output_path=scat_base,
        palette=["#4C72B0"],  # single-color palette
        regline=False,
        hue=None,
        transparent=False,
        figsize = (5,5)
    )

    # Save PNG version as well
    _fig, _ax = _plt.subplots(figsize=(10,5))
    sns.scatterplot(data=scatter_df, x="CloneSize", y="DegreeMix", ax=_ax)
    _ax.set_xlabel("Clone size (cells)", fontweight="bold")
    _ax.set_ylabel(f"Degree of mixing (# cardinal classes < {MIX_THRESHOLD}%)", fontweight="bold")
    _fig.tight_layout()
    _fig.savefig(scat_base + ".png", dpi=300, bbox_inches="tight")
    _plt.close(_fig)



if __name__ == "__main__":
    main()
