"""
Score query cells against reference AB distributions using KDE
percentile (HDR) contours, and write per-cell membership to CSV.

Overview
--------
This script compares the spatial distribution of individual "query" cells
(e.g. clone cells) to a set of "reference" AB distributions (e.g. marker-defined
cell types) in transverse spinal cord coordinates.

For each reference AB distribution, a 2D kernel density estimate (KDE) is
computed on a fixed rectilinear grid and converted into highest-density
regions (HDRs) at user-defined percentiles. Each query cell is then assigned
the *smallest* percentile region it falls into for each AB reference
(e.g. ">=10%", ">=20%", …, or "Outside"). Results are saved in a long-form
CSV table, one row per query cell.

In addition, for each query file, the script produces an overlay plot showing
the clone’s cells on top of the HDR contours of all reference distributions.

Inputs
------
1) `cord.mat` (set by CORD_MAT_PATH)
   - MATLAB file containing a struct `cord` with fields:
       - `cord.x` : 1D array of normalized X coordinates (0–1)
       - `cord.y` : 1D array of normalized Y coordinates (0–1)
   - These define the outline of a single hemicord in normalized space.
     The outline is scaled to the plotting limits (X_LIMS, Y_LIMS) and drawn
     as a pale background in the overlay plots. It does *not* affect the
     KDE computation or scoring; it is purely for visualization.

2) Reference CSV files (REF_DIR)
   - Directory of CSV files, each representing the spatial distribution of
     one AB-defined cell population.
   - Required columns:
       - X/Y coordinate columns that can be inferred by `pick_xy_columns()`.
         The function tries several name patterns (e.g. "X", "Y",
         "Coord_X", "Coord_Y", "Position X", "Position Y",
         "X coordinate (scaled to reference sc)", etc.), and if that fails
         falls back to the first two numeric columns.
       - Optionally a "Hemisphere" column when HEMISPHERES = 'LR', with
         values indicating left/right (e.g. "L", "R", "left", "right").
   - Optional: the last `DROP_LAST_N` rows can be dropped (e.g. landmarks).
   - All reference files are processed with the same hemisphere policy and
     coordinate scaling as the query files (see "Hemisphere handling" below).

3) Query CSV files (QUERY_DIR)
   - Directory of CSV files, each representing one "clone" (or other query
     population) to be scored.
   - Required:
       - Same coordinate conventions as reference files (X/Y columns).
   - Optional metadata columns (carried into the output CSV if present):
       - "Animal", "TM", "Segment", "Clone", "ID"
       - "Hemisphere" (used if HEMISPHERES = 'LR')
   - Each row is assumed to be a single cell. After hemisphere transforms
     and dropping NaNs, each remaining row becomes one scored "PointIndex".

4) Hemisphere handling (HEMISPHERES, SINGLE_HEMISPHERE, LEFT_FILES, RIGHT_FILES)
   - HEMISPHERES controls how left/right information is treated:
       - 'single'   : mirror all cells into a single hemicord.
                     - If SINGLE_HEMISPHERE == 'right': x ← |x|
                     - If SINGLE_HEMISPHERE == 'left' : x ← -|x|
       - 'LR'       : use a "Hemisphere" column in the CSV to flip sign:
                     - left rows:  x ← -x
                     - right rows: x ← +x
       - 'random'   : randomly mirror all cells left or right for each file.
       - 'assigned' : use filename patterns in LEFT_FILES / RIGHT_FILES to
                     force a file to the left or right side; otherwise leave
                     x as-is.
   - This policy is applied identically to both reference and query files so
     that all data live in a consistent coordinate system.

5) Thoracic / lumbar scaling (X_SCALE)
   - All X coordinates are multiplied by X_SCALE before any hemisphere
     handling. For lumbar, typically X_SCALE = 1.0; for thoracic, you can
     use the same factor as in the original MATLAB code (e.g. 1/1.49).

KDE / percentile computation
----------------------------
- All reference KDEs are computed on the same fixed grid:
    - X grid: X_GRID (derived from X_LIMS and HEMISPHERES)
    - Y grid: Y_GRID
    - Number of bins: N_BINS (per axis)
- `gaussian_kde` (SciPy) is used with bandwidth `BW_METHOD` ("scott" by
  default).
- Densities are normalized so that the integral over the grid is ~1, assuming
  uniform cell areas.
- For each reference KDE:
    - `PERCENTILES_QUERY` defines the HDR levels used for scoring.
      Typically something like [0.10, 0.20, ..., 0.90].
    - For each percentile p, the code finds the density threshold `thr(p)`
      such that the most dense p fraction of the probability mass is above
      that threshold (highest-density region).
    - `PERCENTILES_PLOT` similarly defines a (possibly smaller) set of
      percentiles for plotting contours.

Scoring logic
-------------
- For each query file:
    - Coordinates are loaded, optionally downsampled by row removal, and
      transformed by X_SCALE and the chosen hemisphere policy.
    - A simple overlay plot is generated: the spinal cord hemicord outline,
      all reference KDE HDRs (filled contours + contour lines), and the
      query cells as black dots.
    - For scoring, each cell (x, y) is evaluated against each reference
      KDE:
        - The KDE is evaluated via a `RegularGridInterpolator`, giving
          a density value `z` at that point (or NaN outside the grid).
        - `label_from_levels(z, levels_query, PERCENTILES_QUERY)` returns:
            - the *smallest* percentile for which z ≥ thr(p),
              as a string ">=10%", ">=20%", etc., or
            - "Outside" if z does not reach any HDR threshold or z is NaN.
      In other words, lower percentile labels correspond to denser core regions;
      a cell with label ">=10%" lies in the central 10% HDR of that AB
      distribution, ">=20%" in the 10–20% shell, etc., and "Outside" is
      beyond the largest HDR queried.

Outputs
-------
1) Per-cell membership CSV
   - Written to OUTPUT_DIR with a timestamped name, e.g.:
       query_cells_percentile_membership_YYYYMMDD_HHMMSS.csv
   - Columns:
       - "QueryFile"  : basename (without extension) of the query CSV.
       - "Animal"     : (optional) copied from query CSV if present.
       - "TM"         : (optional)
       - "Segment"    : (optional)
       - "Clone"      : (optional)
       - "ID"         : (optional)
       - "PointIndex" : 1-based index of the cell within the filtered query
                        file (after removing NaNs / dropped rows).
       - "X", "Y"     : transformed coordinates used for scoring.
       - One column per reference file:
           - Column name = reference file basename (without ".csv").
           - Entry = string label:
                ">=10%", ">=20%", … based on PERCENTILES_QUERY, or "Outside".
   - Each row is one cell from one query file, with a categorical label
     describing its membership in each AB reference KDE.

2) Overlay plots for each query file
   - For each query CSV, two figures are saved in OUTPUT_DIR:
       - PNG : clone_overlay_<QueryFile>.png
       - SVG : clone_overlay_<QueryFile>.svg
   - Each plot shows:
       - The hemicord outline from `cord.mat` scaled to X_LIMS / Y_LIMS.
       - All reference AB KDE HDRs as filled contours + contour lines:
           - Percentile levels defined by PERCENTILES_PLOT.
           - Colours assigned either via CUSTOM_COLORS (COLOR_MODE="assigned")
             or by cycling through a base colormap (e.g. "tab20").
       - Query cells as black points, plotted in the transformed coordinate
         system (single hemicord or bilateral, depending on HEMISPHERES).

Notes and assumptions
---------------------
- A reference CSV must have at least MIN_KDE_POINTS valid points after
  cleaning; otherwise it is skipped.
- All percentile lists (PERCENTILES_QUERY, PERCENTILES_PLOT) must be
  in ascending order (e.g. [0.10, 0.20, ..., 0.90]).
- The fixed KDE grid ensures that:
    - All references share the same spatial support.
    - Density values are comparable across reference AB distributions.
- Hemisphere handling and scaling are applied consistently to both
  references and queries, so scoring is always between distributions
  defined in the same coordinate frame.
"""

import os, glob, re, datetime
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.interpolate import RegularGridInterpolator
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap
from scipy.io import loadmat
import matplotlib.colors as mcolors  

# Find current date and time, which will be used in file naming
now = datetime.datetime.now()
curr_date = now.strftime("%Y%m%d_%H%M%S")
    
########################################################################################################################################################################################################
### User input
CORD_MAT_PATH = r"Z:\People\Sophie\3 Code\Matlab Plot Code_Mariano\cord.mat"  # Path of the cord.mat file
REF_DIR       = r"Z:\People\Sophie\10 AB Spatial Distribution Map\20251114 AB Spatial Map\20251114 Lumbar Matt and Giulia Counting\AB Files" # Path of the folder ("AB Files") containing the AB distributions (cell coordinates). The KDEs are built from these files
QUERY_DIR     = r"Z:\People\Sophie\10 AB Spatial Distribution Map\20260113 Clone Spatial Analysis\Sparsest dataset normalized clones - Neurons only 80% complete\Output_20260117_102658" # Path of the folder containing csv files of normalized individual clones (output of SCNormalization_BLT)
OUTPUT_DIR    = r"Z:\People\Sophie\10 AB Spatial Distribution Map\20260113 Clone Spatial Analysis\AB vs Clone distribution comparison\Output_Sparsest dataset neurons only 80% complete_" # Desired output folder 
os.makedirs(OUTPUT_DIR + curr_date, exist_ok=True)

# Thoracic vs Lumbar scaling (if needed): set to 1/1.49 for T (as in MATLAB), else 1.0
X_SCALE = 1.0  # e.g., 1.0 for L, or (1/1.49) for T

# If the last N rows are landmarks to drop:
DROP_LAST_N = 0  # set 0 if your CSVs have no trailing landmarks

# Hemisphere handling
HEMISPHERES = 'single'         # 'single' | 'LR' | 'random' | 'assigned'
SINGLE_HEMISPHERE = 'right'    # 'right' or 'left' (used only if HEMISPHERES='single')

# KDE / grid
N_BINS     = 256
X_LIMS     = (0.0, 650.0)
Y_LIMS     = (-400.0, 500.0)
# Set X_GRID based on chosen single hemisphere
if HEMISPHERES == 'single':
    X_GRID = (0.0, X_LIMS[1]) if SINGLE_HEMISPHERE == 'right' else (-X_LIMS[1], 0.0)
else:
    X_GRID = (-X_LIMS[1], X_LIMS[1])
Y_GRID     = (Y_LIMS[0],  Y_LIMS[1])
BW_METHOD  = "scott"
PERCENTILES_QUERY = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]  # ascending!

# Percentiles used for plotting contours (can be different)
PERCENTILES_PLOT  = [0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70]  # tweak as you like

MIN_KDE_POINTS = 5

CLONE_POINT_SIZE = 12
CLONE_POINT_ALPHA = 1.0

# Filled contour style for clone overlays
OVERLAY_FILL_ALPHA = 0.25   
OVERLAY_LINE_ALPHA = 0.8    

# Colour mode for AB distributions in overlays
COLOR_MODE = "assigned"   # options: "tab20" (repeat tab20), 'assigned' if you want to assign colour per filename in CUSTOM_COLORS
CMAP_NAME  = "tab20"      # base colormap (used in tab20/assigned modes)

PASTEL_S = 0.60           # pastel saturation for COLOR_MODE="unique"
PASTEL_V = 0.95           # pastel value/brightness

# Same custom colours as in PlotTransverseDistribution_20250411.py
CUSTOM_COLORS = {
    "V3 (Nkx2.2+).csv": "#E6F598",     # pale green
    "MN (Hb9+).csv": "#8079b7",        # purple
    "MN (ventral Isl+).csv": "#5e4fa2",# purple
    "V1 (En1+).csv": "#D53E4F",        # red
    "V1 (En1+FoxD3+).csv": "#db7785",  # red
    "V0 (Evx1+).csv": "#FDAE61",       # orange
    "V2 (Chx10+).csv": "#66C2A5",      # teal
    "dI1 (Barhl1+).csv": "#ABDDA4",    # green
    "dI2 (FoxD3+En1-).csv": "#FEE08B", # yellow
    "dI3 (Isl+Tlx3+Lbx1-).csv": "#9b9b9b",  # grey
    "dI3 (dorsal Isl+).csv": "#777777",     # grey
    "dI6 (Dmrt3+).csv": "#9B6981FF",        # mauve
    "dI5 (Tlx3+).csv": "#236687",          # blue
    "dI5 (Lmx1b+).csv": "#5baece",         # blue
    "dI5 (Lbx1+Tlx3+).csv": "#3288BD",     # blue
}

# --------------------------
# Helpers
# --------------------------
def get_color_for_file(path: str, idx: int, cmap):
    """
    Assigned colour mode: try to match the filename against CUSTOM_COLORS keys.
    Keys can match full basename, stem, or first 3 underscore-separated parts.
    Falls back to colormap if no match.
    """
    base = os.path.basename(path)
    stem, _ = os.path.splitext(base)
    pref3 = file_key_three_underscores(base)

    for key, col in CUSTOM_COLORS.items():
        if key in (base, stem, pref3):
            return col  # hex string is fine for matplotlib
    # fallback: cycle through colormap
    return cmap(idx % cmap.N)

def _extract_vector(maybe_array):
    """Flatten nested MATLAB struct arrays to a 1D float vector."""
    arr = np.array(maybe_array)
    # unwrap 1x1 object nests
    while arr.dtype == object and arr.size == 1:
        arr = np.array(arr.item())
    return np.ravel(arr).astype(float)

def load_cord_xy(mat_path: str):
    """Load 'cord.x' and 'cord.y' from cord.mat."""
    md = loadmat(mat_path, squeeze_me=True, struct_as_record=False)
    if "cord" not in md:
        raise KeyError(f"'cord' variable not found in {mat_path}")
    cord = md["cord"]

    # handle both struct-with-attributes and dict-like
    if hasattr(cord, "__dict__"):
        x_raw = getattr(cord, "x", None)
        y_raw = getattr(cord, "y", None)
    else:
        x_raw = cord["x"]
        y_raw = cord["y"]

    if x_raw is None or y_raw is None:
        raise KeyError("Could not access 'cord.x' and 'cord.y' in MAT file.")

    x = _extract_vector(x_raw)
    y = _extract_vector(y_raw)
    if x.size != y.size:
        raise ValueError("cord.x and cord.y lengths differ.")
    return x, y

def scale_and_translate_outline(x: np.ndarray, y: np.ndarray, x_lims, y_lims):
    """
    Scale normalized cord outline (0–1) to your X/Y limits and translate so
    the ventral-most and medial-most points align with (x_lims[0], y_lims[0]).
    """
    x0, x1 = x_lims
    y0, y1 = y_lims
    xs = x * (x1 - x0)
    ys = y * (y1 - y0)
    xs = xs - xs.min() + x0
    ys = ys - ys.min() + y0
    return xs, ys

def draw_hemisection_outline(ax, x_lims, y_lims,
                             facecolor=(0.96, 0.96, 0.96),
                             edgecolor="k",
                             lw=1.0):
    """
    Draw a single hemisection outline (no mirroring) on the given axes.
    Assumes cord.mat contains a normalized single hemicord outline.
    """
    x, y = load_cord_xy(CORD_MAT_PATH)
    xs, ys = scale_and_translate_outline(x, y, x_lims, y_lims)
    ax.fill(xs, ys,
            facecolor=facecolor,
            edgecolor=edgecolor,
            linewidth=lw,
            zorder=0)
    
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
    for a, b in (("x", "y"),
                 ("coord_x", "coord_y"),
                 ("pos_x", "pos_y"),
                 ("position x", "position y"),
                 ("x coordinate (scaled to reference sc)",
                  "y coordinate (scaled to reference sc)")):
        if a in lower and b in lower:
            return lower[a], lower[b]
    numeric = df.select_dtypes(include=["number"]).columns.tolist()
    if len(numeric) >= 2: return numeric[0], numeric[1]
    cols = df.columns.tolist()
    if len(cols) >= 2: return cols[0], cols[1]
    raise ValueError("Could not infer X/Y columns.")

def read_xy_points(csv_path: str, hemisphere_policy: str, return_df: bool = False): 
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
        # Find (case-insensitive) hemisphere column
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

    if return_df:
        return df, x, y
    return x, y

def kde2d_on_fixed_grid(x, y, x_grid, y_grid, n_bins=256, bw_method="scott"):
    gx = np.linspace(x_grid[0], x_grid[1], n_bins)
    gy = np.linspace(y_grid[0], y_grid[1], n_bins)
    GX, GY = np.meshgrid(gx, gy)
    kde = gaussian_kde(np.vstack([x, y]), bw_method=bw_method)
    Z = kde(np.vstack([GX.ravel(), GY.ravel()])).reshape(GX.shape)

    # Normalize to integrate to ~1 over the grid (uniform cell areas)
    dx = (gx[-1]-gx[0])/(n_bins-1)
    dy = (gy[-1]-gy[0])/(n_bins-1)
    Z = Z / (Z.sum() * dx * dy + 1e-16)
    return gx, gy, Z

def kde_percentile_levels(Z, percentiles):
    z = Z.ravel()
    z = z[np.isfinite(z)]
    z = z[z > 0]
    if z.size == 0:
        return None
    zs = np.sort(z)[::-1]
    cdf = np.cumsum(zs) / zs.sum()
    levels = []
    for p in percentiles:
        idx = np.searchsorted(cdf, p)
        levels.append(float(zs[min(idx, zs.size-1)]))
    return np.array(levels)

def label_from_levels(z_val, levels, percentiles):
    """
    Given a density value z_val, and descending thresholds for ascending percentiles,
    return the smallest percentile satisfied (e.g. '>=10%') or 'Outside'.
    """
    if z_val is None or not np.isfinite(z_val):
        return "Outside"
    # levels are for ascending percentiles; density threshold decreases as percentile increases
    # We want the FIRST percentile whose threshold is <= z_val
    for thr, p in zip(levels, percentiles):
        if z_val >= thr:
            return f">={int(round(p*100))}%"
    return "Outside"

def regular_grid_interpolator(gx, gy, Z):
    """Bilinear interpolator over (gx,gy) → Z."""
    return RegularGridInterpolator((gy, gx), Z, bounds_error=False, fill_value=np.nan)

def safe_var_name(s: str) -> str:
    # Similar to matlab.lang.makeValidName
    s = re.sub(r"[^\w]", "_", s)
    if re.match(r"^\d", s): s = "_" + s
    return s

def plot_clone_overlay(xq, yq, refs, qname, out_dir):
    """
    Plot the clone (xq,yq) in right hemicord on top of all reference KDEs,
    and save PNG + SVG in out_dir.

    - Uses filled contours (contourf) for KDE HPD regions
    - Contour levels come from the global PERCENTILES list
    - Colours per AB distribution follow COLOR_MODE:
        "assigned" -> CUSTOM_COLORS (by filename) + fallback to CMAP_NAME
        "unique"   -> pastel colours
        "tab20"    -> straight tab20 cycling
    """
    if xq.size == 0 or yq.size == 0:
        return

    fig, ax = plt.subplots(figsize=(5, 5), dpi=300)
    # Draw hemisection
    draw_hemisection_outline(ax, X_LIMS, Y_LIMS,
                             facecolor=(0.96, 0.96, 0.96),
                             edgecolor="k",
                             lw=1.0)

    # Colormap
    cmap = get_cmap(CMAP_NAME)

    # Decide colour sequence / strategy
    if COLOR_MODE == "assigned":
        # we'll call get_color_for_file per ref below
        color_seq = None
    else:  # "tab20" default
        color_seq = [cmap(i % cmap.N) for i in range(len(refs))]

    # --- Draw reference KDEs ---
    for idx, ref in enumerate(refs):
        gx = ref["gx"]          # 1D x grid
        gy = ref["gy"]          # 1D y grid
        Z  = ref["Z"]           # 2D density
        levels_raw = ref["levels_plot"]   # HPD thresholds for plotting (PERCENTILES_PLOT)

        finite_levels = levels_raw[np.isfinite(levels_raw)]
        if finite_levels.size == 0:
            continue

        levels_plot = np.sort(finite_levels)  # must be increasing

        # Pick colour for this ref
        if COLOR_MODE == "assigned":
            color = get_color_for_file(ref["path"], idx, cmap)
        else:
            color = color_seq[idx]

        # Filled contours (HPD regions)
        ax.contourf(
            gx,
            gy,
            Z,
            levels=np.append(levels_plot, Z.max()),
            colors=[color],
            alpha=OVERLAY_FILL_ALPHA,
            antialiased=True,
        )

        # Contour lines
        ax.contour(
            gx,
            gy,
            Z,
            levels=levels_plot,
            colors=[color],
            linewidths=0.7,
            alpha=OVERLAY_LINE_ALPHA,
        )

    # Clone cells as black filled circles
    ax.scatter(
        xq, yq,
        s=CLONE_POINT_SIZE,
        c="k",
        edgecolors="k",
        alpha=CLONE_POINT_ALPHA,
        zorder=5,
    )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlim(X_GRID)
    ax.set_ylim(Y_GRID)
    ax.set_xlabel("X (µm)")
    ax.set_ylabel("Y (µm)")

    hemi_label = "right" if HEMISPHERES == "single" and SINGLE_HEMISPHERE == "right" else "hemicord"
    ax.set_title(f"Clone {qname} vs AB KDEs ({hemi_label})", fontsize=9)

    fig.tight_layout()

    base = re.sub(r"[^\w\-]+", "_", qname)
    png_path = os.path.join(out_dir, f"clone_overlay_{base}.png")
    svg_path = os.path.join(out_dir, f"clone_overlay_{base}.svg")
    fig.savefig(png_path, dpi=300)
    fig.savefig(svg_path)
    plt.close(fig)
    print(f"[PLOT] {png_path}")
    
# --------------------------
# Main
# --------------------------
def main():
    # --- Collect reference files (build KDEs) ---
    ref_paths = sorted(glob.glob(os.path.join(REF_DIR, "*.csv")))
    if not ref_paths:
        raise FileNotFoundError(f"No reference CSVs found in {REF_DIR}")

    refs = []  # list of dicts: {name, gx, gy, Z, levels, interp}
    for rp in ref_paths:
        try:
            xr, yr = read_xy_points(rp, HEMISPHERES)
        except Exception as e:
            print(f"[REF SKIP READ] {os.path.basename(rp)}: {e}")
            continue
        valid = np.isfinite(xr) & np.isfinite(yr)
        xr = xr[valid]; yr = yr[valid]
        if xr.size < MIN_KDE_POINTS:
            print(f"[REF SKIP KDE] {os.path.basename(rp)}: too few points ({xr.size}).")
            continue

        gx, gy, Z = kde2d_on_fixed_grid(
            xr, yr, X_GRID, Y_GRID,
            n_bins=N_BINS,
            bw_method=BW_METHOD,
        )

        levels_query = kde_percentile_levels(Z, PERCENTILES_QUERY)
        levels_plot  = kde_percentile_levels(Z, PERCENTILES_PLOT)

        if levels_query is None or levels_plot is None:
            print(f"[REF SKIP KDE] {os.path.basename(rp)}: KDE levels undefined.")
            continue

        interp = regular_grid_interpolator(gx, gy, Z)
        name = os.path.splitext(os.path.basename(rp))[0]
        refs.append(dict(
            name=name,
            path=rp,           # keep full path for color assignment
            gx=gx,
            gy=gy,
            Z=Z,
            levels_query=levels_query,
            levels_plot=levels_plot,
            interp=interp,
        ))
        name = os.path.splitext(os.path.basename(rp))[0]

    if not refs:
        raise RuntimeError("No valid reference KDEs could be built.")

    # --- Prepare output table columns ---
    out_rows = []
    # --- Score each query file ---
    q_paths = sorted(glob.glob(os.path.join(QUERY_DIR, "*.csv")))
    if not q_paths:
        print(f"[WARN] No query CSVs found in {QUERY_DIR}")

    for qp in q_paths:
        try:
            qdf, xq, yq = read_xy_points(qp, HEMISPHERES, return_df=True)
        except Exception as e:
            print(f"[QUERY SKIP READ] {os.path.basename(qp)}: {e}")
            continue

        valid = np.isfinite(xq) & np.isfinite(yq)
        xq = xq[valid]
        yq = yq[valid]
        qdf = qdf.loc[valid].reset_index(drop=True)

        if xq.size == 0:
            continue

        qname = os.path.splitext(os.path.basename(qp))[0]

        plot_clone_overlay(xq, yq, refs, qname, OUTPUT_DIR)

        for i, (xv, yv) in enumerate(zip(xq, yq), start=1):
            meta = qdf.iloc[i-1]  # same row as this point

            row = {
                "QueryFile": qname,
                "Animal":   meta.get("Animal", ""),
                "TM":       meta.get("TM", ""),
                "Segment":  meta.get("Segment", ""),
                "Clone":    meta.get("Clone", ""),
                "ID":       meta.get("ID", ""),
                "PointIndex": i,
                "X": xv,
                "Y": yv,
            }
            for ref in refs:
                z = float(ref["interp"]((yv, xv)))  # note (y,x) order
                lbl = label_from_levels(z, ref["levels_query"], PERCENTILES_QUERY)
                col = ref["name"]          # <<< column name = ref file basename
                row[col] = lbl

            out_rows.append(row)

    # --- Write CSV ---
    df_out = pd.DataFrame(out_rows)
    ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_csv = os.path.join(OUTPUT_DIR, f"query_cells_percentile_membership_{ts}.csv")
    if not df_out.empty:
        df_out.to_csv(out_csv, index=False)
        print(f"[SAVED] {out_csv}")
    else:
        print("[WARN] No query rows scored; CSV not written.")

if __name__ == "__main__":
    main()
