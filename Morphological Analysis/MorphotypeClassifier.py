"""
================================================================================
UNIFIED NEURONAL ANGLE ANALYSIS — FINAL MORPHOTYPE CLASSIFICATION + COMPLEXITY
SWC (± CSV METADATA) PIPELINE • PER-CLASS PANELS • EXPORT TABLES
================================================================================
Purpose:
Classify reconstructed neurons into 5 geometry-based types using the angular
arrangement of primary projections measured at the FIRST BRANCH POINT or at 50 µm
from the soma (whichever occurs first), then add a Low/High complexity label using
T100 (Sholl total intersections to 100 µm computed directly from SWC) to create the Morphotype.

Core method:
- For each primary projection, find the analysis point:
    • first branch point along that branch, OR
    • the node closest to TARGET_R = 50 µm from soma
- Compute 2D polar angles (XY) of the resulting vectors and their gap structure:
    • L1 = largest gap, L2 = second largest gap, Ratio = L2/L1
- Assign FINAL geometry class (no hemidirectional):
    1) Unidirectional
    2) Bipolar
    3) Triangular
    4) Stellate
    5) Bidirectional
  with empirical decision boundaries (L1, Ratio, n_primary) and a distance-to-ideal
  fallback for ambiguous cases.

Complexity (SWC-only, no external merge):
- Sholl intersections are computed from SWC (3D or XY) with step SHOLL_STEP_UM.
- T100 = cumulative Sholl intersections from 0–100 µm.
- Complexity is assigned per base type:
    • Low  : T100 < mean(T100) for that base type
    • High : T100 ≥ mean(T100) for that base type
- Final category label:
    Morphotype = <Geometry>_<Low/High>   (stored as cell_type_with_complexity)

Input modes:
- CSV+SWC mode:
    • CSV_PATH provided → cells grouped from neurite rows and matched to SWCs
- SWC-only mode:
    • CSV_PATH missing/not found → cells inferred directly from folder hierarchy

Required SWC hierarchy (for CloneKey parsing):
  .../<Litter S###>/<Animal S###_a#>/<Section S1..S6>/<Clone clone*>/<cell>.swc

Outputs (saved to OUT_DIR):
1) cell_analysis_detailed_with_complexity.csv
   - IDs: CellKey, Label, CloneKey
   - Geometry (base type) + Morphotype (type+complexity)
   - n_primary, L1/L2/Ratio, T100, pairwise angles (wide format)
2) Figures (SVG):
   1_method_schematic_with_complexity.svg
   2_comprehensive_overview_with_complexity.svg
   3_cell_type_distribution_with_complexity.svg
   4_clone_composition_with_complexity.svg
   5_empirical_distributions.svg
   6_complexity_analysis.svg
3) Per-class overview folder:
   OUT_DIR/PER_CLASS_OVERVIEWS/class_<Morphotype>_pXX_of_YY.svg  (paged)

Key parameters (edit in CONFIG):
- TARGET_R (default 50 µm)
- SHOLL_STEP_UM, USE_3D_SHOLL, RADII (default [100] → T100)
- Empirical classification thresholds (L1/Ratio bounds)
================================================================================
"""

import os, re, glob, math, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import patheffects as pe
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy import stats
import seaborn as sns

# ==================== CONFIGURATION ====================
CSV_PATH = None #r"Z:\People\Francisco\Code_testing_folder\Input\imaris_cellinfo_with_total_length_with_map__len_gt400.csv" #csv file directory or None (to run swc only - the metadata is computed based on the hierarchy of folder)
SWC_ROOT = r"Z:\People\Francisco\Code_testing_folder\Input\Neuromorpho Inh_Exc"
OUT_DIR = r"Z:\People\Francisco\Code_testing_folder\Output"
PER_CLASS_DIR = os.path.join(OUT_DIR, "PER_CLASS_OVERVIEWS")

# Analysis parameters
TARGET_R = 50.0  # Analysis at 50μm or first branch point
# ==================== SHOLL (FROM SWC) ====================
USE_3D_SHOLL = True      # True=3D radii, False=XY only
SHOLL_STEP_UM = 2        # ring spacing
RADII = [100]  # will compute T100

# FINAL Cell type classification criteria (DATA-DRIVEN)
# Based on empirical data analysis from your distributions
UNIDIRECTIONAL_L1_THRESHOLD = 300.0  # L1 > 300° for multi-primary cells

BIPOLAR_MIN = 130.0
BIPOLAR_MAX = 230.0

# Bidirectional: True clean cluster from your data
BIDIRECTIONAL_L1_MIN = 160.0  # From your clean bidirectional cluster (was 160, but 165 better)
BIDIRECTIONAL_L1_MAX = 299.0
BIDIRECTIONAL_RATIO_MAX = 0.99  # Strict ratio for clean bidirectional

# Stellate: Tight cluster from your data
STELLATE_L1_MIN = 45.0
STELLATE_L1_MAX = 159.9  # From your stellate distribution (90-136°)
STELLATE_RATIO_MIN = 0.33  # High ratio for stellate

# Triangular
TRIANGULAR_MIN = 75.0
TRIANGULAR_MAX = 145.0

# Visualization settings
DPI = 150
LINE_WIDTH = 1.2  # Thicker lines for better visibility
VECTOR_WIDTH = 3.5  # Thicker vectors
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']
plt.rcParams['axes.linewidth'] = 1.5
warnings.filterwarnings("ignore")
os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PER_CLASS_DIR, exist_ok=True)

print("\n" + "="*80)
print("UNIFIED ANGLE ANALYSIS - FINAL CELL TYPE CLASSIFICATION WITH COMPLEXITY")
print("="*80)

# ==================== UTILITY FUNCTIONS ====================
def norm(s):
    return re.sub(r"[^a-z0-9]+", "_", str(s).strip().lower())

def find_col(df, must, exclude=None):
    """Flexible column finder"""
    if exclude is None:
        exclude = []
    nmap = {c: norm(c) for c in df.columns}
    cands = [
        (col, len(nc)) for col, nc in nmap.items()
        if all(t in nc for t in must) and not any(x in nc for x in exclude)
    ]
    if not cands:
        return None
    cands.sort(key=lambda x: x[1])
    return cands[0][0]

def canon_litter(s):
    if pd.isna(s):
        return ""
    m = re.search(r"[sS]\d+", str(s))
    return m.group(0).upper() if m else str(s).strip().upper()

def canon_animal(s):
    if pd.isna(s):
        return ""
    tok = str(s).strip().split("_")[-1]
    m = re.search(r"[aA]\d+", tok)
    if not m:
        m = re.search(r"\d+", tok)
        return ("a" + m.group(0).lstrip("0")) if m else tok.lower()
    return m.group(0).lower()

def canon_section(s):
    if pd.isna(s):
        return ""
    m = re.search(r"[sS]\d+", str(s))
    if m:
        return m.group(0).lower()
    m = re.search(r"\d+", str(s))
    return ("s" + m.group(0)) if m else str(s).lower()

def canon_clone(s):
    if s is None or (isinstance(s, float) and math.isnan(s)):
        return "nan"
    ss = str(s).strip()
    if ss == "" or ss.lower() == "nan":
        return "nan"

    s0 = ss.lower()
    s0 = re.sub(r"(?i)\bclone\b", "", s0).strip()

    m = re.fullmatch(r"c(\d+)", s0)
    if m:
        return f"clone_{m.group(1)}"

    if re.fullmatch(r"\d+", s0):
        return f"clone_{s0.lstrip('0') or '0'}"

    core = re.sub(r"[^a-z0-9]+", "_", s0).strip("_")
    return f"clone_{core}" if core else "clone"

def id_from_text(s):
    m = re.search(r"\d+", str(s))
    return m.group(0) if m else ""

# ==================== SWC PATH PARSING ====================
def parse_swc_hierarchy(path):
    parts = re.split(r"[\\/]+", os.path.normpath(path))

    litter = ""
    animal = ""
    section = ""
    clone = ""

    # First pass: detect combined litter_animal folder like "S127_a20"
    for p in parts:
        m = re.fullmatch(r"([sS]\d+)_([aA]\d+)", p)
        if m:
            litter = m.group(1).upper()
            animal = m.group(2).lower()
            break

    # Second pass: detect separate litter / animal / section
    seen_animal = False
    for p in parts:
        if not litter and re.fullmatch(r"[sS]\d+", p):
            litter = p.upper()
            continue

        if not animal and re.fullmatch(r"[aA]\d+", p):
            animal = p.lower()
            seen_animal = True
            continue

        # section = first s# after animal (not the litter)
        if (seen_animal or animal) and not section and re.fullmatch(r"[sS]\d+", p):
            if p.upper() != litter:
                section = p.lower()
                continue

        if not clone and "clone" in p.lower():
            clone = p
            continue

    return litter, animal, section, clone

def cell_id_from_filename(path):
    stem = os.path.splitext(os.path.basename(path))[0]   # removes .swc
    # also remove trailing extra suffix like ".CNG" if present
    stem = re.sub(r"\.(cng|CNG)$", "", stem)

    # Prefer explicit "cell-<n>" or "cell_<n>"
    m = re.search(r"(?:^|[-_])cell[-_]?(\d+)\b", stem, flags=re.IGNORECASE)
    if m:
        return m.group(1)

    # Fallback: last number in the string (more stable than first)
    nums = re.findall(r"\d+", stem)
    return nums[-1] if nums else stem

def create_cell_key(path):
        litter, animal, section, clone = parse_swc_hierarchy(path)

        lit = canon_litter(litter)

        # animal can be "S058_a6" OR "a6" depending on parser; normalize to "S058_a6"
        ani_token = ""
        if animal:
            # if already like S058_a6
            m = re.fullmatch(r"([sS]\d+)_([aA]\d+)", str(animal).strip())
            if m:
                lit = canon_litter(m.group(1))
                ani_token = canon_animal(m.group(2))
            else:
                ani_token = canon_animal(animal)

        ani_full = f"{lit}_{ani_token}" if (lit and ani_token) else (lit or ani_token or "unknown")

        sec = canon_section(section)
        cln = canon_clone(clone)

        # Cell ID: for your dataset filenames like "63.swc" this returns "63"
        cid = id_from_text(os.path.splitext(os.path.basename(path))[0])

        clone_key = f"{ani_full}-{sec}-{cln}"
        cell_label = f"{ani_full}-{sec}-{cln}-{cid}"
        cell_key = f"{ani_full}|{sec}|{cln}|{cid}"
        return cell_key, cell_label, clone_key

# ==================== CSV PROCESSING ====================
def load_cell_metadata(csv_path):
    print(f"\nLoading cell metadata from CSV (neurites per line)...")
    df = pd.read_csv(csv_path, dtype=str)

    # Find columns flexibly
    LIT = find_col(df, ["litter"]) or "Litter"
    ANI = find_col(df, ["animal"]) or "Animal"
    SEC = find_col(df, ["segment"]) or find_col(df, ["section"]) or "Segment"
    ID_COL = find_col(df, ["id"], exclude=["project"]) or "ID"
    CLONE_COL = find_col(df, ["clone"]) or "Clone"

    # Canonicalize to match SWC keys
    df["_lit"] = df[LIT].map(canon_litter)
    df["_ani"] = df[ANI].map(canon_animal)           # "S058_a6" -> "a6"
    df["_sec"] = df[SEC].map(canon_section)          # "s1" -> "s1"
    df["_cln"] = df[CLONE_COL].map(canon_clone)      # "c3" -> "clone_3" (with our improved canon_clone)
    df["_cid"] = df[ID_COL].map(id_from_text)        # "63" -> "63"

    # Build animal_full like "S058_a6"
    df["_ani_full"] = np.where(
        (df["_lit"] != "") & (df["_ani"] != ""),
        df["_lit"] + "_" + df["_ani"],
        np.where(df["_lit"] != "", df["_lit"], np.where(df["_ani"] != "", df["_ani"], "unknown"))
    )

    # Group rows into cells (each group = one cell)
    group_cols = ["_ani_full", "_sec", "_cln", "_cid"]
    g = df.groupby(group_cols, dropna=False)

    # One row per cell
    cells_df = g.first().reset_index()

    # How many neurites (rows) belonged to this cell
    cells_df["n_neurites_rows"] = g.size().values

    # Build keys used everywhere else
    cells_df["CellKey"] = (
        cells_df["_ani_full"].astype(str) + "|" +
        cells_df["_sec"].astype(str) + "|" +
        cells_df["_cln"].astype(str) + "|" +
        cells_df["_cid"].astype(str)
    )
    cells_df["Label"] = (
        cells_df["_ani_full"].astype(str) + "-" +
        cells_df["_sec"].astype(str) + "-" +
        cells_df["_cln"].astype(str) + "-" +
        cells_df["_cid"].astype(str)
    )
    cells_df["CloneKey"] = (
        cells_df["_ani_full"].astype(str) + "-" +
        cells_df["_sec"].astype(str) + "-" +
        cells_df["_cln"].astype(str)
    )

    print(f"Found {len(cells_df)} unique cells in CSV (grouped from {len(df)} neurite rows)")
    return cells_df

def build_cell_metadata_from_swc(swc_root):
    print("\nLoading cell metadata from SWC files (SWC-only mode)...")
    swc_files = glob.glob(os.path.join(swc_root, "**", "*.swc"), recursive=True)

    if not swc_files:
        print(f"ERROR: No SWC files found under: {swc_root}")
        return pd.DataFrame()

    rows = []
    key_counts = {}

    for swc_path in sorted(swc_files):
        cell_key, cell_label, clone_key = create_cell_key(swc_path)

        # If collision still happens, make key unique without losing hierarchy
        if cell_key in key_counts:
            key_counts[cell_key] += 1
            cell_key_unique = f"{cell_key}__dup{key_counts[cell_key]}"
            cell_label = f"{cell_label}__dup{key_counts[cell_key]}"
        else:
            key_counts[cell_key] = 0
            cell_key_unique = cell_key

        rows.append({
            "CellKey": cell_key_unique,
            "Label": cell_label,
            "CloneKey": clone_key,
            "SWC_Path": swc_path,
            "Source": "SWC_ONLY"
        })

    cells_df = pd.DataFrame(rows)
    print(f"Found {len(cells_df)} SWC files (kept all as cells).")
    print(f"Unique CellKey count: {cells_df['CellKey'].nunique()}")
    return cells_df

    print("SWC CellKey examples:", list(swc_lookup.keys())[:5])
    print("CSV CellKey examples:", cells_df["CellKey"].head(5).tolist())
    print("SWC CellKey examples:", list(swc_lookup.keys())[:5])

    csv_keys = set(cells_df["CellKey"])
    swc_keys = set(swc_lookup.keys())
    print("Exact matches:", len(csv_keys & swc_keys))

def load_cells_metadata():
    """
    Decide whether to load from CSV or SWC-only.
    """
    use_csv = isinstance(CSV_PATH, str) and CSV_PATH.strip() != "" and os.path.exists(CSV_PATH)

    if use_csv:
        return load_cell_metadata(CSV_PATH)

    print("\n[INFO] CSV not provided or not found -> running in SWC-only mode.")
    return build_cell_metadata_from_swc(SWC_ROOT)

# ==================== SWC ANALYSIS ====================
def load_swc(path):
    """Load SWC file"""
    nodes = {}
    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            p = s.split()
            if len(p) < 7:
                continue
            try:
                nid = int(float(p[0]))
                ntype = int(float(p[1]))
                x, y, z = float(p[2]), float(p[3]), float(p[4])
                r = float(p[5])
                parent = int(float(p[6]))
                nodes[nid] = dict(type=ntype, x=x, y=y, z=z, r=r, parent=parent)
            except:
                pass
    return nodes

def soma_center_and_ids(nodes):
    """Find soma center and IDs"""
    soma = [i for i, n in nodes.items() if n["type"] == 1]
    if soma:
        cx = np.mean([nodes[i]["x"] for i in soma])
        cy = np.mean([nodes[i]["y"] for i in soma])
        cz = np.mean([nodes[i]["z"] for i in soma])
        return (float(cx), float(cy), float(cz)), set(soma)

    roots = [i for i, n in nodes.items() if n["parent"] == -1]
    if roots:
        r0 = nodes[roots[0]]
        return (r0["x"], r0["y"], r0["z"]), {roots[0]}

    any_id = next(iter(nodes))
    n = nodes[any_id]
    return (n["x"], n["y"], n["z"]), {any_id}

def primary_children(nodes, soma_ids):
    """Find primary dendrites"""
    return [i for i, n in nodes.items() if n["parent"] in soma_ids]

def dist_3d(n, center):
    """3D distance"""
    dx = n["x"] - center[0]
    dy = n["y"] - center[1]
    dz = n["z"] - center[2]
    return math.sqrt(dx*dx + dy*dy + dz*dz)

def node_dist_from_soma(node, center, use_3d=True):
    dx = node["x"] - center[0]
    dy = node["y"] - center[1]
    if use_3d:
        dz = node["z"] - center[2]
        return math.sqrt(dx*dx + dy*dy + dz*dz)
    return math.sqrt(dx*dx + dy*dy)

def sholl_intersections_from_swc(nodes, center, rings, use_3d=True):
    """
    Count intersections at each radius ring.
    Segment crosses ring r if its endpoints straddle r.
    """
    dist = {nid: node_dist_from_soma(nd, center, use_3d=use_3d) for nid, nd in nodes.items()}
    counts = {r: 0 for r in rings}

    for nid, nd in nodes.items():
        pid = nd["parent"]
        if pid not in nodes:
            continue

        d1 = dist[nid]
        d2 = dist[pid]
        if d1 == d2:
            continue

        lo, hi = (d1, d2) if d1 < d2 else (d2, d1)
        for r in rings:
            if lo < r <= hi:
                counts[r] += 1

    return counts

def sholl_totals_Tcols(nodes, center, RADII, step_um=5, use_3d=True):
    """
    Returns dict: T50, T100, ... (sum of intersections from 0..R).
    """
    max_r = int(max(RADII))
    rings = list(range(step_um, max_r + step_um, step_um))  # 5,10,...,max

    I = sholl_intersections_from_swc(nodes, center, rings, use_3d=use_3d)

    # cumulative sum over rings
    cum = {}
    running = 0
    for r in rings:
        running += I[r]
        cum[r] = running

    out = {}
    for R in RADII:
        Rint = int(R)
        snapped = step_um * int(round(Rint / step_um))
        snapped = max(step_um, min(snapped, max_r))
        out[f"T{Rint}"] = cum.get(snapped, np.nan)

    return out

def find_first_branch_point(nodes, start_id):
    """Find the first branching point along a dendrite"""
    cur = start_id
    visited = set()

    for _ in range(10000):
        visited.add(cur)
        children = [i for i, n in nodes.items()
                    if n["parent"] == cur and i not in visited]

        if len(children) >= 2:
            # This is a branch point
            return cur
        elif len(children) == 1:
            cur = children[0]
        else:
            # Terminal point
            return cur

    return cur

def project_to_radius_or_branch(nodes, center, start_id, target_r=TARGET_R):
    """
    Project to first branch point OR 50μm, whichever comes first.
    Returns: (node_id, distance_from_soma, is_branch_point)
    """
    cur = start_id
    visited = set()

    for _ in range(10000):
        visited.add(cur)
        cur_dist = dist_3d(nodes[cur], center)

        # Check if we've reached 50μm
        if cur_dist >= target_r:
            # Find the closest point to exactly 50μm
            best = cur
            best_diff = abs(cur_dist - target_r)

            # Check if parent is closer to 50μm
            if nodes[cur]["parent"] in nodes:
                parent = nodes[cur]["parent"]
                parent_dist = dist_3d(nodes[parent], center)
                parent_diff = abs(parent_dist - target_r)
                if parent_diff < best_diff:
                    best = parent

            return best, dist_3d(nodes[best], center), False

        # Check for branch point
        children = [i for i, n in nodes.items()
                    if n["parent"] == cur and i not in visited]

        if len(children) >= 2:
            # Branch point found before 50μm
            return cur, cur_dist, True
        elif len(children) == 1:
            cur = children[0]
        else:
            # Terminal point before 50μm
            return cur, cur_dist, False

    return cur, dist_3d(nodes[cur], center), False

def compute_pairwise_angles_with_details(rays):
    """
    Compute all pairwise angles with full details.
    Returns: list of (angle, ray_i_idx, ray_j_idx)
    """
    angle_details = []
    m = len(rays)
    for i in range(m):
        for j in range(i + 1, m):
            dot = np.clip(np.dot(rays[i], rays[j]), -1.0, 1.0)
            ang = math.degrees(math.acos(dot))
            # Keep in 0-180 range
            if ang > 180:
                ang = 360 - ang
            angle_details.append((ang, i, j))
    return angle_details

def analyze_cell_angles(nodes):
    """Perform angle analysis at first branch point or 50μm"""
    center, soma_ids = soma_center_and_ids(nodes)
    prims = primary_children(nodes, soma_ids)

    if not prims:
        return None

    # Get direction vectors at branch point or 50μm
    rays = []
    analysis_points = []

    for p in prims:
        node_id, dist, is_branch = project_to_radius_or_branch(nodes, center, p, TARGET_R)
        nd = nodes[node_id]

        v = np.array([
            nd["x"] - center[0],
            nd["y"] - center[1],
            nd["z"] - center[2]
        ], dtype=float)
        norm_v = np.linalg.norm(v)
        if norm_v > 0:
            rays.append(v / norm_v)
            analysis_points.append({
                'node_id': node_id,
                'distance': dist,
                'is_branch': is_branch,
                'position': (nd["x"], nd["y"], nd["z"])
            })

    if len(rays) < 1:
        return None

    # Compute pairwise angles with details
    angle_details = compute_pairwise_angles_with_details(rays)
    angles = [a[0] for a in angle_details]

    # Calculate gaps (angles between consecutive rays when sorted)
    if len(rays) >= 2:
        # Convert to 2D angles for gap calculation
        angles_2d = []
        for ray in rays:
            angle_2d = math.degrees(math.atan2(ray[1], ray[0]))
            if angle_2d < 0:
                angle_2d += 360
            angles_2d.append(angle_2d)

        angles_2d_sorted = sorted(angles_2d)
        gaps = []
        for i in range(len(angles_2d_sorted)):
            next_i = (i + 1) % len(angles_2d_sorted)
            gap = angles_2d_sorted[next_i] - angles_2d_sorted[i]
            if gap < 0:
                gap += 360
            gaps.append(gap)

        gaps_sorted = sorted(gaps, reverse=True)
    else:
        gaps_sorted = []

    return {
        'rays': rays,
        'angles': angles,
        'angle_details': angle_details,
        'n_primary': len(prims),
        'n_pairs': len(angles),
        'center': center,
        'analysis_points': analysis_points,
        'gaps': gaps_sorted
    }

def classify_cell_type(n_primary, gaps):
    """
    IMPROVED Classification based on empirical data:

    1. UNIDIRECTIONAL:
       - 1 primary projection OR
       - >1 projections with L1 > 300°

    2. BIPOLAR: 2 primary projections

    3. TRIANGULAR: EXACTLY 3 primaries AND smallest gap > 50°

    4. STELLATE:
       - ≥4 primaries
       - AND L1 ≤ 135°
       - AND L2/L1 ratio ≥ 0.80

    5. BIDIRECTIONAL:
       - ≥3 primaries
       - AND L1 ≥ 165°
       - AND L2/L1 ratio ≤ 0.70

    For cells in the intermediate zone (L1 = 135-165°):
       - Use ratio to decide:
         • Ratio ≥ 0.85 → Stellate (high ratio = more even distribution)
         • Ratio ≤ 0.65 → Bidirectional (low ratio = two clear lobes)
         • Otherwise, assign to closest match based on distance to ideal criteria
    """
    # If no gaps (n_primary==1), then L1 is not defined
    # For n_primary>1, we compute L1
    if len(gaps) == 0:
        return "Unidirectional" if n_primary == 1 else "Stellate"

    L1 = gaps[0] if len(gaps) > 0 else 0

    # UNIDIRECTIONAL: 1 primary OR L1 > 300°
    if n_primary == 1 or (n_primary > 1 and L1 > UNIDIRECTIONAL_L1_THRESHOLD):
        return "Unidirectional"

    # Bipolar: ONLY 2 primaries
    if n_primary == 2:
        return "Bipolar"

    # Triangular: EXACTLY 3 primaries AND smallest gap > 50°
    if n_primary == 3:
        min_gap = min(gaps) if len(gaps) > 0 else 0
        if min_gap > 50.0:
            return "Triangular"
        # If it fails the triangular criterion, fall through

    # For ≥3 primaries
    if n_primary >= 3:
        # Compute L1, L2 and ratio
        if len(gaps) >= 2:
            L1 = gaps[0]
            L2 = gaps[1]
            R = L2 / L1 if L1 > 0 else 0
        else:
            L1 = gaps[0] if len(gaps) > 0 else 0
            L2 = 0
            R = 0

        # Clean Stellate: L1 ≤ 135°, R ≥ 0.80, and ≥4 primaries
        if n_primary >= 4 and L1 <= STELLATE_L1_MAX and R >= STELLATE_RATIO_MIN:
            return "Stellate"

        # Clean Bidirectional: L1 ≥ 165°, R ≤ 0.70, and ≥3 primaries
        if L1 >= BIDIRECTIONAL_L1_MIN and R <= BIDIRECTIONAL_RATIO_MAX:
            return "Bidirectional"

        # Intermediate zone: L1 between 135-165°
        # Use ratio to disambiguate
        if 135.0 <= L1 <= 165.0:
            if R >= 0.85:
                return "Stellate"
            elif R <= 0.65:
                return "Bidirectional"
            # If ratio is between 0.65-0.85, use distance calculation

        # For cells that don't fit clean criteria, use distance to ideal criteria
        # This is simpler and more transparent than point scoring

        # Calculate "distance" to each type's ideal
        # Smaller distance = better match

        # Distance to Stellate ideal: L1 ≤ 135°, R ≥ 0.80
        stellate_distance = 0
        if n_primary >= 4:
            stellate_distance += 0  # Good, meets primary count requirement
        else:
            stellate_distance += 20  # Penalize for not having ≥4 primaries

        if L1 <= STELLATE_L1_MAX:
            stellate_distance += (STELLATE_L1_MAX - L1) / 10.0  # Closer to 135° is better
        else:
            stellate_distance += (L1 - STELLATE_L1_MAX) * 2.0  # Penalize for being above 135°

        if R >= STELLATE_RATIO_MIN:
            stellate_distance += (R - STELLATE_RATIO_MIN) * 5.0  # Higher ratio is better
        else:
            stellate_distance += (STELLATE_RATIO_MIN - R) * 10.0  # Penalize for low ratio

        # Distance to Bidirectional ideal: L1 ≥ 165°, R ≤ 0.70
        bidirectional_distance = 0
        if L1 >= BIDIRECTIONAL_L1_MIN:
            bidirectional_distance += (L1 - BIDIRECTIONAL_L1_MIN) / 10.0  # Closer to 165° is better
        else:
            bidirectional_distance += (BIDIRECTIONAL_L1_MIN - L1) * 2.0  # Penalize for being below 165°

        if R <= BIDIRECTIONAL_RATIO_MAX:
            bidirectional_distance += (BIDIRECTIONAL_RATIO_MAX - R) * 5.0  # Lower ratio is better
        else:
            bidirectional_distance += (R - BIDIRECTIONAL_RATIO_MAX) * 10.0  # Penalize for high ratio

        # If 3 primaries, slightly favor bidirectional (since stellate requires ≥4)
        if n_primary == 3:
            stellate_distance += 10

        # Choose type with smaller distance (closer match)
        if stellate_distance <= bidirectional_distance:
            return "Stellate"
        else:
            return "Bidirectional"

    # Default fallback
    return "Stellate"

# ==================== COMPLEXITY (T100 FROM SWC) ====================
def add_complexity_classification(analysis_df):
    print("\n" + "=" * 80)
    print("ADDING COMPLEXITY CLASSIFICATION (T100 COMPUTED FROM SWC)")
    print("=" * 80)

    if "T100" not in analysis_df.columns:
        print("WARNING: T100 not found in analysis_df. Complexity will be 'Unknown'.")
        analysis_df["complexity"] = "Unknown"
        analysis_df["cell_type_with_complexity"] = analysis_df["cell_type"] + "_Unknown"
        return analysis_df

    analysis_df["T100"] = pd.to_numeric(analysis_df["T100"], errors="coerce")

    type_means = analysis_df.groupby("cell_type")["T100"].mean()

    analysis_df["complexity"] = "Unknown"
    for ct, mean_val in type_means.items():
        m = (analysis_df["cell_type"] == ct) & (analysis_df["T100"].notna())
        analysis_df.loc[m & (analysis_df["T100"] < mean_val), "complexity"] = "Low"
        analysis_df.loc[m & (analysis_df["T100"] >= mean_val), "complexity"] = "High"

    analysis_df["cell_type_with_complexity"] = analysis_df["cell_type"] + "_" + analysis_df["complexity"]

    print("\nT100 means by cell type:")
    for ct, mv in type_means.items():
        print(f"  {ct}: mean T100 = {mv:.2f}")

    print("\nComplexity distribution:")
    for ct in sorted(analysis_df["cell_type"].unique()):
        for comp in ["Low", "High", "Unknown"]:
            n = len(analysis_df[(analysis_df["cell_type"] == ct) & (analysis_df["complexity"] == comp)])
            if n > 0:
                print(f"  {ct}_{comp}: {n} cells")

    return analysis_df

# ==================== BUILD ANALYSIS TABLE ====================
def build_angle_analysis_table(cells_df):
    """Analyze SWCs and build the per-cell analysis table.
       SWC-only mode: use SWC_Path directly (no re-scanning / no matching).
       CSV mode: scan SWC_ROOT and match by CellKey.
    """
    print(f"\nScanning for SWC files...")

    # ---------------- SWC-ONLY MODE ----------------
    if "SWC_Path" in cells_df.columns and cells_df["SWC_Path"].notna().all():
        swc_lookup = {row["CellKey"]: row["SWC_Path"] for _, row in cells_df.iterrows()}
        print("CSV CellKey examples:", cells_df["CellKey"].head(5).tolist())
        print("SWC CellKey examples:", list(swc_lookup.keys())[:5])

        csv_keys = set(cells_df["CellKey"])
        swc_keys = set(swc_lookup.keys())
        print("Exact matches:", len(csv_keys & swc_keys))

        # Compare without the last field (cell id)
        csv_loose = set("|".join(k.split("|")[:3]) for k in csv_keys)
        swc_loose = set("|".join(k.split("|")[:3]) for k in swc_keys)
        print("Matches without CellID:", len(csv_loose & swc_loose))
        print(f"Using {len(swc_lookup)} SWC paths from metadata (SWC-only mode).")

        results = []
        skipped = []  # (CellKey, reason)

        for _, cell in cells_df.iterrows():
            cell_key = cell["CellKey"]
            swc_path = cell["SWC_Path"]

            try:
                nodes = load_swc(swc_path)
                if len(nodes) < 10:
                    skipped.append((cell_key, "too_few_nodes(<10)"))
                    continue

                analysis = analyze_cell_angles(nodes)
                if analysis is None:
                    skipped.append((cell_key, "angle_analysis_none(no_primaries?)"))
                    continue

                cell_type = classify_cell_type(analysis["n_primary"], analysis["gaps"])

                row = {
                    "CellKey": cell_key,
                    "Label": cell.get("Label", os.path.basename(swc_path)),
                    "CloneKey": cell.get("CloneKey", "ROOT"),
                    "SWC_Path": swc_path,
                    "n_primary": analysis["n_primary"],
                    "n_pairs": analysis["n_pairs"],
                    "cell_type": cell_type,
                    "largest_gap": analysis["gaps"][0] if len(analysis["gaps"]) > 0 else np.nan,
                    "second_gap": analysis["gaps"][1] if len(analysis["gaps"]) > 1 else np.nan,
                    "all_gaps": analysis["gaps"],
                    "pairwise_angles": analysis["angles"],
                    "angle_details": analysis["angle_details"],
                    "rays": analysis["rays"],
                    "center": analysis["center"],
                    "analysis_points": analysis["analysis_points"],
                }

                # ---- SHOLL TOTALS FROM SWC (adds T50/T100/...) ----
                tcols = sholl_totals_Tcols(
                    nodes,
                    analysis["center"],
                    RADII=RADII,
                    step_um=SHOLL_STEP_UM,
                    use_3d=USE_3D_SHOLL,
                )
                row.update(tcols)

                results.append(row)

            except Exception as e:
                skipped.append((cell_key, f"exception:{type(e).__name__}"))
                print(f"Error processing {swc_path}: {e}")

        df = pd.DataFrame(results)
        print(f"\nSuccessfully processed {len(df)} cells")

        if skipped:
            print(f"Skipped {len(skipped)} cells. Top reasons:")
            # summarize skip reasons
            from collections import Counter
            c = Counter([r for _, r in skipped])
            for reason, n in c.most_common(10):
                print(f"  {reason}: {n}")
            # optional: save skipped list
            pd.DataFrame(skipped, columns=["CellKey", "Reason"]).to_csv(
                os.path.join(OUT_DIR, "skipped_cells_debug.csv"), index=False
            )
            print(f"[Saved] Skip report: {os.path.join(OUT_DIR, 'skipped_cells_debug.csv')}")

        return df

    # ---------------- CSV MODE (your original matching) ----------------
    swc_files = glob.glob(os.path.join(SWC_ROOT, "**", "*.swc"), recursive=True)
    print(f"Found {len(swc_files)} SWC files")

    if len(swc_files) == 0:
        print(f"\nERROR: No SWC files found in {SWC_ROOT}")
        return pd.DataFrame()

    swc_lookup = {}
    for swc_path in swc_files:
        cell_key, _, _ = create_cell_key(swc_path)
        swc_lookup[cell_key] = swc_path

    print(f"Indexed {len(swc_lookup)} unique cells from SWC files")

    results = []
    skipped = []

    for _, cell in cells_df.iterrows():
        cell_key = cell["CellKey"]
        if cell_key not in swc_lookup:
            skipped.append((cell_key, "no_match_in_swc_lookup"))
            continue

        swc_path = swc_lookup[cell_key]

        try:
            nodes = load_swc(swc_path)
            if len(nodes) < 10:
                skipped.append((cell_key, "too_few_nodes(<10)"))
                continue

            analysis = analyze_cell_angles(nodes)
            if analysis is None:
                skipped.append((cell_key, "angle_analysis_none(no_primaries?)"))
                continue

            cell_type = classify_cell_type(analysis["n_primary"], analysis["gaps"])

            row = {
                "CellKey": cell_key,
                "Label": cell["Label"],
                "CloneKey": cell["CloneKey"],
                "SWC_Path": swc_path,
                "n_primary": analysis["n_primary"],
                "n_pairs": analysis["n_pairs"],
                "cell_type": cell_type,
                "largest_gap": analysis["gaps"][0] if len(analysis["gaps"]) > 0 else np.nan,
                "second_gap": analysis["gaps"][1] if len(analysis["gaps"]) > 1 else np.nan,
                "all_gaps": analysis["gaps"],
                "pairwise_angles": analysis["angles"],
                "angle_details": analysis["angle_details"],
                "rays": analysis["rays"],
                "center": analysis["center"],
                "analysis_points": analysis["analysis_points"],
            }

            tcols = sholl_totals_Tcols(
                nodes,
                analysis["center"],
                RADII=RADII,
                step_um=SHOLL_STEP_UM,
                use_3d=USE_3D_SHOLL,
            )
            row.update(tcols)

            results.append(row)

        except Exception as e:
            skipped.append((cell_key, f"exception:{type(e).__name__}"))
            print(f"Error processing {swc_path}: {e}")

    df = pd.DataFrame(results)
    print(f"\nSuccessfully processed {len(df)} cells")

    if skipped:
        print(f"Skipped {len(skipped)} cells. Top reasons:")
        from collections import Counter
        c = Counter([r for _, r in skipped])
        for reason, n in c.most_common(10):
            print(f"  {reason}: {n}")

    return df

# ==================== VISUALIZATIONS ====================

def get_cell_type_color(cell_type):
    """Get consistent color for each cell type"""
    # First, extract base type if it has complexity suffix
    base_type = cell_type
    if '_Low' in cell_type:
        base_type = cell_type.replace('_Low', '')
    elif '_High' in cell_type:
        base_type = cell_type.replace('_High', '')

    # Base colors for original types
    base_colors = {
        'Unidirectional': '#E63946',      # Red
        'Bipolar': '#F77F00',             # Orange
        'Bidirectional': '#06D6A0',       # Teal
        'Triangular': '#118AB2',          # Blue
        'Stellate': '#8338EC',            # Purple
    }

    # Get base color
    if base_type in base_colors:
        base_color = base_colors[base_type]
    else:
        base_color = '#999999'

    # Adjust brightness for complexity
    import colorsys
    # Convert hex to RGB
    r = int(base_color[1:3], 16) / 255.0
    g = int(base_color[3:5], 16) / 255.0
    b = int(base_color[5:7], 16) / 255.0

    # Convert to HLS
    h, l, s = colorsys.rgb_to_hls(r, g, b)

    # Adjust lightness based on complexity
    if '_Low' in cell_type:
        # Lighten for low complexity
        l = min(1.0, l * 1.4)
    elif '_High' in cell_type:
        # Darken for high complexity
        l = max(0.0, l * 0.6)

    # Convert back to RGB
    r, g, b = colorsys.hls_to_rgb(h, l, s)

    # Convert back to hex
    return f'#{int(r*255):02x}{int(g*255):02x}{int(b*255):02x}'

def plot_method_schematic(example_df, out_path, pick_label_contains=None, pick_cellkey=None, pick_index=0):
    """Create clean schematic showing the analysis method"""
    if len(example_df) == 0:
        return

    # Get first example
    # Choose which cell to use
    chosen_df = example_df

    if pick_cellkey is not None:
        chosen_df = chosen_df[chosen_df["CellKey"] == pick_cellkey]

    if pick_label_contains is not None:
        chosen_df = chosen_df[chosen_df["Label"].astype(str).str.contains(pick_label_contains, case=False, na=False)]

    if len(chosen_df) == 0:
        print("[WARN] Requested schematic cell not found; falling back to first cell.")
        chosen_df = example_df

    example = chosen_df.iloc[pick_index].copy()
    nodes = load_swc(example['SWC_Path'])
    center = example['center']
    soma_ids = set([i for i, n in nodes.items() if n["type"] == 1])
    prims = primary_children(nodes, soma_ids)

    fig = plt.figure(figsize=(22, 6), dpi=DPI)  # wider
    gs = GridSpec(
        1, 3,
        figure=fig,
        wspace=0.35,
        width_ratios=[1.25, 1.0, 1.75]  # give the text panel more room
    )

    # === Panel A: Full neuron with analysis points ===
    ax1 = fig.add_subplot(gs[0, 0])

    # Draw neuron with thicker lines
    for nid, nd in nodes.items():
        pid = nd["parent"]
        if pid in nodes:
            ax1.plot(
                [nodes[pid]["x"], nd["x"]],
                [nodes[pid]["y"], nd["y"]],
                'k-',
                linewidth=LINE_WIDTH,
                alpha=0.6
            )

    # Draw soma
    ax1.plot(
        center[0], center[1],
        'o',
        color='#E63946',
        markersize=14,
        zorder=10,
        label='Soma'
    )

    # Draw 50μm circle
    circle = plt.Circle(
        (center[0], center[1]),
        TARGET_R,
        fill=False,
        edgecolor='#118AB2',
        linewidth=2.5,
        linestyle='--',
        alpha=0.8,
        label=f'{int(TARGET_R)}μm radius'
    )
    ax1.add_patch(circle)

    # Highlight analysis points
    for i, point_info in enumerate(example['analysis_points']):
        pos = point_info['position']
        marker = 's' if point_info['is_branch'] else 'o'
        label = 'Branch point' if i == 0 and point_info['is_branch'] else (
            '50μm point' if i == 0 else None
        )
        ax1.plot(
            pos[0], pos[1],
            marker,
            color='#06D6A0',
            markersize=10,
            zorder=11,
            markeredgecolor='white',
            markeredgewidth=1.5,
            label=label
        )

        # Draw line from soma to analysis point with thicker line
        ax1.plot(
            [center[0], pos[0]],
            [center[1], pos[1]],
            'r-',
            linewidth=VECTOR_WIDTH,
            alpha=0.8,
            zorder=9
        )

    ax1.set_aspect('equal')
    ax1.set_title(
        'A. Analysis Points Detection',
        fontsize=14,
        fontweight='bold',
        pad=15
    )
    ax1.legend(loc='best', fontsize=10, framealpha=0.95)
    ax1.set_xlabel('X (μm)', fontsize=11)
    ax1.set_ylabel('Y (μm)', fontsize=11)
    ax1.grid(True, alpha=0.2, linestyle=':')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # === Panel B: Polar projection ===
    ax2 = fig.add_subplot(gs[0, 1], projection='polar')

    rays = example['rays']
    colors_ray = plt.cm.Set2(np.linspace(0, 1, len(rays)))

    for i, ray in enumerate(rays):
        theta = np.arctan2(ray[1], ray[0])
        ax2.plot(
            [0, theta],
            [0, 1],
            linewidth=4,
            color=colors_ray[i],
            label=f'Branch {i+1}',
            alpha=0.9
        )
        ax2.plot(
            theta, 1,
            'o',
            color=colors_ray[i],
            markersize=10,
            markeredgecolor='white',
            markeredgewidth=1.5
        )

    ax2.set_ylim(0, 1.2)
    ax2.set_title(
        'B. Polar Projection',
        fontsize=14,
        fontweight='bold',
        pad=20
    )
    ax2.grid(True, alpha=0.3, linestyle=':')
    ax2.set_theta_zero_location('E')
    ax2.legend(
        loc='upper left',
        bbox_to_anchor=(1.15, 1.0),
        fontsize=9,
        framealpha=0.95
    )

    # === Panel C: FINAL Classification criteria ===
    ax3 = fig.add_subplot(gs[0, 2])
    ax3.axis('off')

    criteria_text = """
IMPROVED CELL TYPE CLASSIFICATION WITH COMPLEXITY
(BASED ON EMPIRICAL DATA ANALYSIS)

Method: Analysis at first branch point
        OR 50μm (whichever comes first)

Types (with Complexity):
━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
I. UNIDIRECTIONAL (Low/High Complexity)
   • 1 primary projection OR
   • >1 projections with L1 > 300°

II. BIPOLAR (Low/High Complexity)
   • 2 primaries only

III. TRIANGULAR (Low/High Complexity)
   • EXACTLY 3 primaries
   • Smallest gap > 50°

IV. STELLATE (Low/High Complexity)
   • ≥4 primaries
   • L1 ≤ 135°
   • L2/L1 ratio ≥ 0.80

V. BIDIRECTIONAL (Low/High Complexity)
   • ≥3 primaries
   • L1 ≥ 165°
   • L2/L1 ratio ≤ 0.70

COMPLEXITY CLASSIFICATION:
   • Based on T100 (Total Sholl intersections)
   • Low: T100 < mean for that cell type
   • High: T100 ≥ mean for that cell type

INTERMEDIATE ZONE (L1=135-165°):
   • Ratio ≥ 0.85 → Stellate
   • Ratio ≤ 0.65 → Bidirectional
   • Otherwise → Closest match to ideal criteria
    """

    ax3.text(
        0.03, 0.98, criteria_text,
        transform=ax3.transAxes,
        fontsize=9,  # was 10
        va='top',
        family='monospace',
        bbox=dict(boxstyle='round,pad=1.0', facecolor='white', edgecolor='#118AB2', linewidth=2, alpha=0.95)
    )

    plt.suptitle(
        'Angle Analysis Method & Improved Empirical Classification WITH COMPLEXITY',
        fontsize=16,
        fontweight='bold',
        y=0.98
    )

    plt.savefig(out_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[Saved] Method schematic: {out_path}")

def plot_comprehensive_overview(df, out_path):
    """Create a comprehensive overview with multiple examples per type"""
    # Get all cell types with complexity
    cell_types = sorted(df['cell_type_with_complexity'].dropna().unique())

    base_types = ['Unidirectional', 'Bipolar', 'Triangular', 'Stellate', 'Bidirectional']
    comp_order = ["Low", "High", "Unknown"]

    ordered_types = []
    for base in base_types:
        for comp in comp_order:
            t = f"{base}_{comp}"
            if t in cell_types:
                ordered_types.append(t)

    # Fallback: plot whatever exists if Low/High/Unknown naming doesn't match
    if not ordered_types:
        ordered_types = cell_types

    # Hard safety
    if len(ordered_types) == 0:
        print("[Skip] Comprehensive overview: no categories found")
        return

    examples_per_type = 4  # Show 4 examples per type (reduced due to more types)

    fig = plt.figure(figsize=(20, len(ordered_types) * 3), dpi=DPI)
    gs = GridSpec(
        len(ordered_types),
        examples_per_type * 2,
        figure=fig,
        hspace=0.5,
        wspace=0.3
    )

    for type_idx, cell_type in enumerate(ordered_types):
        type_cells = df[df['cell_type_with_complexity'] == cell_type]
        if len(type_cells) == 0:
            continue

        # Take up to examples_per_type cells
        examples = type_cells.head(examples_per_type)

        for example_idx, (_, example) in enumerate(examples.iterrows()):
            nodes = load_swc(example['SWC_Path'])
            center = example['center']

            # Morphology subplot
            morph_col = example_idx * 2
            ax_morph = fig.add_subplot(gs[type_idx, morph_col])

            # Draw full morphology with thicker lines
            for nid, nd in nodes.items():
                pid = nd["parent"]
                if pid in nodes:
                    ax_morph.plot(
                        [nodes[pid]["x"], nd["x"]],
                        [nodes[pid]["y"], nd["y"]],
                        'k-',
                        linewidth=LINE_WIDTH,
                        alpha=0.4
                    )

            # Draw soma
            ax_morph.plot(
                center[0], center[1],
                'o',
                color='#E63946',
                markersize=8,
                zorder=10,
                markeredgecolor='white',
                markeredgewidth=1.5
            )

            # Draw vectors
            for point_info, ray in zip(example['analysis_points'], example['rays']):
                pos = point_info['position']
                color = get_cell_type_color(cell_type)

                # Draw vector line
                ax_morph.plot(
                    [center[0], pos[0]],
                    [center[1], pos[1]],
                    color=color,
                    linewidth=VECTOR_WIDTH,
                    alpha=0.8,
                    zorder=9
                )

                # Draw analysis point
                marker = 's' if point_info['is_branch'] else 'o'
                ax_morph.plot(
                    pos[0], pos[1],
                    marker,
                    color=color,
                    markersize=6,
                    zorder=11,
                    markeredgecolor='white',
                    markeredgewidth=1
                )

            ax_morph.set_aspect('equal')
            ax_morph.axis('off')

            # Add cell info
            gaps_text = f"L1: {example['largest_gap']:.0f}°"
            if not pd.isna(example['second_gap']):
                ratio = (
                    example['second_gap'] / example['largest_gap']
                    if example['largest_gap'] > 0 else 0
                )
                gaps_text += f"\nL2: {example['second_gap']:.0f}°\nRatio: {ratio:.2f}"

            ax_morph.set_title(
                f'{example["Label"]}\n{gaps_text}',
                fontsize=9,
                color=get_cell_type_color(cell_type),
                pad=5
            )

            # Polar subplot
            polar_col = example_idx * 2 + 1
            ax_polar = fig.add_subplot(gs[type_idx, polar_col], projection='polar')

            rays = example['rays']
            cell_color = get_cell_type_color(cell_type)

            # Convert rays to angles and plot
            for i, ray in enumerate(rays):
                angle_2d = math.degrees(math.atan2(ray[1], ray[0]))
                if angle_2d < 0:
                    angle_2d += 360

                theta = np.radians(angle_2d)

                # Plot vector
                ax_polar.plot(
                    [0, theta],
                    [0, 1],
                    linewidth=2.5,
                    color=cell_color,
                    alpha=0.8
                )
                ax_polar.plot(
                    theta, 1,
                    'o',
                    color=cell_color,
                    markersize=5,
                    markeredgecolor='white',
                    markeredgewidth=1
                )

                # Add angle label for first few vectors
                if i < 3:  # Label only first 3 to avoid clutter
                    ax_polar.text(
                        theta, 1.1,
                        f'{angle_2d:.0f}°',
                        ha='center',
                        va='bottom',
                        fontsize=7,
                        fontweight='bold',
                        bbox=dict(
                            boxstyle='round',
                            facecolor='white',
                            alpha=0.8
                        )
                    )

            ax_polar.set_ylim(0, 1.2)
            ax_polar.grid(True, alpha=0.3, linestyle=':')
            ax_polar.set_theta_zero_location('E')
            ax_polar.set_title(f'n={example["n_primary"]}', fontsize=9, pad=5)

            # Add type label only for first column
            if example_idx == 0:
                ax_morph.text(
                    -0.3, 0.5,
                    cell_type,
                    transform=ax_morph.transAxes,
                    fontsize=11,
                    fontweight='bold',
                    va='center',
                    ha='right',
                    color=get_cell_type_color(cell_type)
                )

    plt.suptitle(
        'Comprehensive Cell Type Overview - Multiple Examples per Type (WITH COMPLEXITY)\n'
        '(Improved Empirical Criteria with Low/High Complexity Classification)',
        fontsize=16,
        fontweight='bold',
        y=0.98
    )

    # Reserve top area for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    print(f"[Saved] Comprehensive overview: {out_path}")

def plot_cell_type_distribution(df, out_path):
    """Plot distribution of cell types with complexity"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6), dpi=DPI)

    # Count by type with complexity
    type_counts = df['cell_type_with_complexity'].value_counts()

    # Order types logically
    base_types = ['Unidirectional', 'Bipolar', 'Triangular', 'Stellate', 'Bidirectional']
    ordered_types = []
    for base in base_types:
        low_type = f"{base}_Low"
        high_type = f"{base}_High"
        if low_type in type_counts.index:
            ordered_types.append(low_type)
        if high_type in type_counts.index:
            ordered_types.append(high_type)

    type_counts = type_counts.reindex(ordered_types, fill_value=0)

    colors = [get_cell_type_color(ct) for ct in type_counts.index]

    # Bar plot
    bars = ax1.bar(
        range(len(type_counts)),
        type_counts.values,
        color=colors,
        edgecolor='black',
        linewidth=1.5,
        alpha=0.85
    )
    ax1.set_xticks(range(len(type_counts)))
    ax1.set_xticklabels(type_counts.index, rotation=45, ha='right', fontsize=10)
    ax1.set_ylabel('Number of Cells', fontsize=12, fontweight='bold')
    ax1.set_title(
        'Cell Type Distribution\n(With Complexity Classification)',
        fontsize=14,
        fontweight='bold',
        pad=15
    )
    ax1.grid(True, alpha=0.3, axis='y', linestyle=':')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Add count labels on bars
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax1.text(
                bar.get_x() + bar.get_width() / 2.,
                height,
                f'{int(height)}',
                ha='center',
                va='bottom',
                fontsize=9,
                fontweight='bold'
            )

    # Pie chart
    wedges, texts, autotexts = ax2.pie(
        type_counts.values,
        labels=type_counts.index,
        colors=colors,
        autopct='%1.1f%%',
        startangle=90,
        textprops={'fontsize': 9}
    )

    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(9)

    ax2.set_title(
        'Cell Type Proportions\n(With Complexity Classification)',
        fontsize=14,
        fontweight='bold',
        pad=15
    )

    # Reserve top area for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    print(f"[Saved] Cell type distribution: {out_path}")

def plot_clone_composition(df, out_path):
    """Plot cell type composition within clones (with complexity)"""
    # Only include clones with at least 2 cells
    clone_sizes = df.groupby('CloneKey').size()
    multi_cell_clones = clone_sizes[clone_sizes > 1].index
    multi_clone_df = df[df['CloneKey'].isin(multi_cell_clones)]

    if len(multi_clone_df) == 0:
        print("[Skip] No clones with multiple cells found")
        return

    clone_composition = multi_clone_df.groupby(['CloneKey', 'cell_type_with_complexity']).size().unstack(fill_value=0)

    # Order columns logically
    base_types = ['Unidirectional', 'Bipolar', 'Triangular', 'Stellate', 'Bidirectional']
    ordered_columns = []
    for base in base_types:
        low_type = f"{base}_Low"
        high_type = f"{base}_High"
        if low_type in clone_composition.columns:
            ordered_columns.append(low_type)
        if high_type in clone_composition.columns:
            ordered_columns.append(high_type)

    clone_composition = clone_composition.reindex(columns=ordered_columns, fill_value=0)

    # Sort by total cells
    clone_composition['Total'] = clone_composition.sum(axis=1)
    clone_composition = clone_composition.sort_values('Total', ascending=False)
    clone_composition = clone_composition.drop('Total', axis=1)

    # Only show top clones if there are too many
    max_clones_to_show = 25
    if len(clone_composition) > max_clones_to_show:
        clone_composition = clone_composition.head(max_clones_to_show)

    fig, ax = plt.subplots(figsize=(14, max(6, len(clone_composition) * 0.5)), dpi=DPI)

    colors = [get_cell_type_color(ct) for ct in ordered_columns]

    # Create stacked bar plot
    bottom = np.zeros(len(clone_composition))
    for i, cell_type in enumerate(ordered_columns):
        counts = clone_composition[cell_type].values
        ax.barh(
            range(len(clone_composition)),
            counts,
            left=bottom,
            color=colors[i],
            edgecolor='white',
            linewidth=0.5,
            label=cell_type,
            alpha=0.8
        )
        bottom += counts

    ax.set_yticks(range(len(clone_composition)))
    ax.set_yticklabels(clone_composition.index, fontsize=9)
    ax.set_xlabel('Number of Cells', fontsize=12, fontweight='bold')
    ax.set_ylabel('Clone', fontsize=12, fontweight='bold')
    ax.set_title(
        f'Cell Type Composition within Clones (With Complexity)\n(Showing {len(clone_composition)} clones with ≥2 cells)',
        fontsize=14,
        fontweight='bold',
        pad=15
    )
    ax.legend(
        bbox_to_anchor=(1.05, 1),
        loc='upper left',
        frameon=True,
        fontsize=8
    )
    ax.grid(True, alpha=0.3, axis='x', linestyle=':')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add count labels on bars
    for i, (idx, row) in enumerate(clone_composition.iterrows()):
        total = row.sum()
        ax.text(
            total + 0.1,
            i,
            f'{int(total)}',
            va='center',
            fontsize=9,
            fontweight='bold'
        )

    # Reserve top area for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    print(f"[Saved] Clone composition: {out_path}")

# ==================== ADDITIONAL ANALYSIS ====================
def plot_empirical_distributions(df, out_path):
    """Plot empirical distributions of L1 and ratios to show data-driven criteria"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=DPI)

    # Filter cells with at least 3 primaries for meaningful analysis
    multi_cells = df[df['n_primary'] >= 3]

    # Calculate ratios
    ratios = []
    for _, row in multi_cells.iterrows():
        if not pd.isna(row['largest_gap']) and not pd.isna(row['second_gap']):
            if row['largest_gap'] > 0:
                ratios.append(row['second_gap'] / row['largest_gap'])
            else:
                ratios.append(np.nan)
        else:
            ratios.append(np.nan)

    multi_cells = multi_cells.copy()
    multi_cells['ratio'] = ratios

    # Separate bidirectional and stellate for comparison (use base types for this plot)
    bidirectional_cells = multi_cells[multi_cells['cell_type'] == 'Bidirectional']
    stellate_cells = multi_cells[multi_cells['cell_type'] == 'Stellate']
    triangular_cells = multi_cells[multi_cells['cell_type'] == 'Triangular']

    # Panel A: L1 distribution comparison
    ax1 = axes[0, 0]

    # Plot histograms
    if len(bidirectional_cells) > 0:
        ax1.hist(bidirectional_cells['largest_gap'].dropna(),
                bins=15, alpha=0.6, color=get_cell_type_color('Bidirectional'),
                label=f'Bidirectional (n={len(bidirectional_cells)})', density=True)

    if len(stellate_cells) > 0:
        ax1.hist(stellate_cells['largest_gap'].dropna(),
                bins=15, alpha=0.6, color=get_cell_type_color('Stellate'),
                label=f'Stellate (n={len(stellate_cells)})', density=True)

    if len(triangular_cells) > 0:
        ax1.hist(triangular_cells['largest_gap'].dropna(),
                bins=10, alpha=0.6, color=get_cell_type_color('Triangular'),
                label=f'Triangular (n={len(triangular_cells)})', density=True)

    # Add threshold lines
    ax1.axvline(x=BIDIRECTIONAL_L1_MIN, color=get_cell_type_color('Bidirectional'),
                linestyle='--', linewidth=2, label=f'Bidirectional min ({BIDIRECTIONAL_L1_MIN}°)')
    ax1.axvline(x=STELLATE_L1_MAX, color=get_cell_type_color('Stellate'),
                linestyle='--', linewidth=2, label=f'Stellate max ({STELLATE_L1_MAX}°)')

    ax1.set_xlabel('Largest Gap (L1) [degrees]', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Density', fontsize=12, fontweight='bold')
    ax1.set_title('A. Largest Gap (L1) Distributions', fontsize=14, fontweight='bold', pad=15)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel B: Ratio distribution comparison
    ax2 = axes[0, 1]

    # Plot histograms
    if len(bidirectional_cells) > 0:
        ax2.hist(bidirectional_cells['ratio'].dropna(),
                bins=15, alpha=0.6, color=get_cell_type_color('Bidirectional'),
                label=f'Bidirectional', density=True)

    if len(stellate_cells) > 0:
        ax2.hist(stellate_cells['ratio'].dropna(),
                bins=15, alpha=0.6, color=get_cell_type_color('Stellate'),
                label=f'Stellate', density=True)

    # Add threshold lines
    ax2.axvline(x=BIDIRECTIONAL_RATIO_MAX, color=get_cell_type_color('Bidirectional'),
                linestyle='--', linewidth=2, label=f'Bidirectional max ({BIDIRECTIONAL_RATIO_MAX})')
    ax2.axvline(x=STELLATE_RATIO_MIN, color=get_cell_type_color('Stellate'),
                linestyle='--', linewidth=2, label=f'Stellate min ({STELLATE_RATIO_MIN})')

    ax2.set_xlabel('L2/L1 Ratio', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Density', fontsize=12, fontweight='bold')
    ax2.set_title('B. L2/L1 Ratio Distributions', fontsize=14, fontweight='bold', pad=15)
    ax2.legend(fontsize=10)
    ax2.grid(True, alpha=0.3)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    # Panel C: Scatter plot L1 vs Ratio
    ax3 = axes[1, 0]

    # Plot scatter points
    if len(bidirectional_cells) > 0:
        ax3.scatter(bidirectional_cells['largest_gap'], bidirectional_cells['ratio'],
                   color=get_cell_type_color('Bidirectional'), s=60, alpha=0.7,
                   label='Bidirectional', edgecolors='white', linewidth=1)

    if len(stellate_cells) > 0:
        ax3.scatter(stellate_cells['largest_gap'], stellate_cells['ratio'],
                   color=get_cell_type_color('Stellate'), s=60, alpha=0.7,
                   label='Stellate', edgecolors='white', linewidth=1)

    if len(triangular_cells) > 0:
        ax3.scatter(triangular_cells['largest_gap'], triangular_cells['ratio'],
                   color=get_cell_type_color('Triangular'), s=60, alpha=0.7,
                   label='Triangular', edgecolors='white', linewidth=1)

    # Add decision boundaries
    ax3.axvline(x=BIDIRECTIONAL_L1_MIN, color=get_cell_type_color('Bidirectional'),
                linestyle='--', linewidth=1.5, alpha=0.7)
    ax3.axvline(x=STELLATE_L1_MAX, color=get_cell_type_color('Stellate'),
                linestyle='--', linewidth=1.5, alpha=0.7)
    ax3.axhline(y=BIDIRECTIONAL_RATIO_MAX, color=get_cell_type_color('Bidirectional'),
                linestyle='--', linewidth=1.5, alpha=0.7)
    ax3.axhline(y=STELLATE_RATIO_MIN, color=get_cell_type_color('Stellate'),
                linestyle='--', linewidth=1.5, alpha=0.7)

    # Shade regions
    ax3.axvspan(BIDIRECTIONAL_L1_MIN, BIDIRECTIONAL_L1_MAX, alpha=0.1, color=get_cell_type_color('Bidirectional'))
    ax3.axvspan(STELLATE_L1_MIN, STELLATE_L1_MAX, alpha=0.1, color=get_cell_type_color('Stellate'))
    ax3.axhspan(0, BIDIRECTIONAL_RATIO_MAX, alpha=0.1, color=get_cell_type_color('Bidirectional'))
    ax3.axhspan(STELLATE_RATIO_MIN, 1.0, alpha=0.1, color=get_cell_type_color('Stellate'))

    ax3.set_xlabel('Largest Gap (L1) [degrees]', fontsize=12, fontweight='bold')
    ax3.set_ylabel('L2/L1 Ratio', fontsize=12, fontweight='bold')
    ax3.set_title('C. L1 vs Ratio: Classification Space', fontsize=14, fontweight='bold', pad=15)
    ax3.legend(fontsize=10)
    ax3.grid(True, alpha=0.3)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # Panel D: Summary statistics table
    ax4 = axes[1, 1]
    ax4.axis('off')

    # Calculate statistics
    stats_text = "IMPROVED EMPIRICAL CRITERIA SUMMARY (n≥3 primaries)\n"
    stats_text += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"

    for cell_type in ['Bidirectional', 'Stellate', 'Triangular']:
        cells = multi_cells[multi_cells['cell_type'] == cell_type]
        if len(cells) > 0:
            l1_mean = cells['largest_gap'].mean()
            l1_std = cells['largest_gap'].std()
            l1_min = cells['largest_gap'].min()
            l1_max = cells['largest_gap'].max()

            if 'ratio' in cells.columns:
                ratio_mean = cells['ratio'].mean()
                ratio_std = cells['ratio'].std()
                ratio_min = cells['ratio'].min()
                ratio_max = cells['ratio'].max()
            else:
                ratio_mean = ratio_std = ratio_min = ratio_max = np.nan

            stats_text += f"{cell_type.upper()} (n={len(cells)}):\n"
            stats_text += f"  L1: {l1_mean:.1f} ± {l1_std:.1f}° ({l1_min:.0f}-{l1_max:.0f}°)\n"
            if not np.isnan(ratio_mean):
                stats_text += f"  Ratio: {ratio_mean:.3f} ± {ratio_std:.3f} ({ratio_min:.3f}-{ratio_max:.3f})\n"
            stats_text += "\n"

    stats_text += "\nIMPROVED EMPIRICAL CRITERIA:\n"
    stats_text += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    stats_text += f"• Bidirectional: L1 ≥ {BIDIRECTIONAL_L1_MIN}°, Ratio ≤ {BIDIRECTIONAL_RATIO_MAX}\n"
    stats_text += f"• Stellate: L1 ≤ {STELLATE_L1_MAX}°, Ratio ≥ {STELLATE_RATIO_MIN}, n≥4\n"
    stats_text += "• Intermediate zone (135-165°):\n"
    stats_text += "  - Ratio ≥ 0.85 → Stellate\n"
    stats_text += "  - Ratio ≤ 0.65 → Bidirectional\n"
    stats_text += "  - Otherwise → Closest match to ideal criteria\n"

    ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes,
             fontsize=10, verticalalignment='top',
             family='monospace',
             bbox=dict(boxstyle='round,pad=0.8', facecolor='white',
                      edgecolor='#118AB2', linewidth=2, alpha=0.95))

    plt.suptitle('Improved Empirical Distributions & Classification Criteria',
                 fontsize=16, fontweight='bold', y=0.98)
    # Reserve top area for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    print(f"[Saved] Empirical distributions: {out_path}")

def plot_complexity_analysis(df, out_path):
    """Plot complexity analysis showing T100 distributions by cell type"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10), dpi=DPI)

    # Find T100 column
    t100_col = None
    for col in df.columns:
        if 'T100' in str(col):
            t100_col = col
            break

    if t100_col is None:
        print("WARNING: T100 column not found for complexity analysis")
        return

    # Panel A: T100 distribution by base cell type
    ax1 = axes[0, 0]
    base_types = ['Unidirectional', 'Bipolar', 'Triangular', 'Stellate', 'Bidirectional']

    for base_type in base_types:
        type_data = df[df['cell_type'] == base_type]
        if len(type_data) > 0:
            # Calculate mean
            mean_t100 = type_data[t100_col].mean()

            # Plot histogram
            ax1.hist(type_data[t100_col].dropna(),
                    bins=15, alpha=0.6,
                    color=get_cell_type_color(base_type),
                    label=f'{base_type} (n={len(type_data)}, mean={mean_t100:.1f})',
                    density=True)

    ax1.set_xlabel('T100 (Total Sholl Intersections at 100µm)', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Density', fontsize=12, fontweight='bold')
    ax1.set_title('A. T100 Distribution by Cell Type', fontsize=14, fontweight='bold', pad=15)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # Panel B: Box plot of T100 by cell type with complexity
    ax2 = axes[0, 1]

    # Prepare data for box plot
    box_data = []
    box_labels = []
    box_colors = []

    base_types = ['Unidirectional', 'Bipolar', 'Triangular', 'Stellate', 'Bidirectional']
    for base_type in base_types:
        for complexity in ['Low', 'High']:
            full_type = f"{base_type}_{complexity}"
            type_data = df[df['cell_type_with_complexity'] == full_type]
            if len(type_data) > 0:
                box_data.append(type_data[t100_col].dropna().values)
                box_labels.append(full_type)
                box_colors.append(get_cell_type_color(full_type))

    if box_data:
        box_plot = ax2.boxplot(box_data, labels=box_labels, patch_artist=True)

        # Color the boxes
        for patch, color in zip(box_plot['boxes'], box_colors):
            patch.set_facecolor(color)
            patch.set_alpha(0.6)

        ax2.set_xticklabels(box_labels, rotation=45, ha='right')
        ax2.set_ylabel('T100', fontsize=12, fontweight='bold')
        ax2.set_title('B. T100 Distribution by Cell Type with Complexity', fontsize=14, fontweight='bold', pad=15)
        ax2.grid(True, alpha=0.3, axis='y')
        ax2.spines['top'].set_visible(False)
        ax2.spines['right'].set_visible(False)

    # Panel C: Scatter plot of T100 vs L1
    ax3 = axes[1, 0]

    # Use color based on cell_type_with_complexity
    for cell_type in df['cell_type_with_complexity'].unique():
        type_data = df[df['cell_type_with_complexity'] == cell_type]
        if len(type_data) > 0:
            ax3.scatter(type_data['largest_gap'], type_data[t100_col],
                       color=get_cell_type_color(cell_type),
                       s=60, alpha=0.7, label=cell_type,
                       edgecolors='white', linewidth=1)

    ax3.set_xlabel('Largest Gap (L1) [degrees]', fontsize=12, fontweight='bold')
    ax3.set_ylabel('T100', fontsize=12, fontweight='bold')
    ax3.set_title('C. T100 vs Largest Gap', fontsize=14, fontweight='bold', pad=15)
    ax3.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    ax3.grid(True, alpha=0.3)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    # Panel D: Summary statistics
    ax4 = axes[1, 1]
    ax4.axis('off')

    # Calculate summary statistics
    stats_text = "COMPLEXITY ANALYSIS SUMMARY\n"
    stats_text += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n"

    for base_type in ['Unidirectional', 'Bipolar', 'Triangular', 'Stellate', 'Bidirectional']:
        type_data = df[df['cell_type'] == base_type]
        if len(type_data) > 0:
            mean_t100 = type_data[t100_col].mean()
            std_t100 = type_data[t100_col].std()
            min_t100 = type_data[t100_col].min()
            max_t100 = type_data[t100_col].max()

            low_count = len(df[(df['cell_type'] == base_type) & (df['complexity'] == 'Low')])
            high_count = len(df[(df['cell_type'] == base_type) & (df['complexity'] == 'High')])

            stats_text += f"{base_type.upper()}:\n"
            stats_text += f"  Mean T100: {mean_t100:.1f} ± {std_t100:.1f}\n"
            stats_text += f"  Range: {min_t100:.0f} - {max_t100:.0f}\n"
            stats_text += f"  Low complexity: {low_count} cells\n"
            stats_text += f"  High complexity: {high_count} cells\n\n"

    stats_text += "\nCOMPLEXITY CLASSIFICATION:\n"
    stats_text += "━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n"
    stats_text += "• Low complexity: T100 < mean for that cell type\n"
    stats_text += "• High complexity: T100 ≥ mean for that cell type\n"

    ax4.text(0.05, 0.95, stats_text, transform=ax4.transAxes,
             fontsize=9, verticalalignment='top',
             family='monospace',
             bbox=dict(boxstyle='round,pad=0.8', facecolor='white',
                      edgecolor='#06D6A0', linewidth=2, alpha=0.95))

    plt.suptitle('Complexity Analysis Based on T100 (Total Sholl Intersections)',
                 fontsize=16, fontweight='bold', y=0.98)
    # Reserve top area for suptitle
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(out_path, dpi=300, bbox_inches='tight', pad_inches=0.2)
    plt.close()
    print(f"[Saved] Complexity analysis: {out_path}")

# ==================== PER-CLASS OVERVIEW ====================
def plot_one_figure_per_class(df, out_dir, cells_per_page=20, ncols_cells=2):
    """
    UPDATED ONLY:
    Per-class overview is now PAGED + TWO-COLUMN (overview-like) layout,
    to prevent gigantic vertical SVGs that open empty/crash with many cells.

    Output:
      class_<Morphotype>_p01_of_YY.svg, class_<Morphotype>_p02_of_YY.svg, ...
    Call stays the same:
      plot_one_figure_per_class(analysis_df, PER_CLASS_DIR)
    """
    if df.empty:
        print("[Skip] Per-class overviews: empty dataframe")
        return

    os.makedirs(out_dir, exist_ok=True)

    # Get all cell types with complexity
    cell_types = sorted(df['cell_type_with_complexity'].dropna().unique())

    # Order them logically
    base_types = ['Unidirectional', 'Bipolar', 'Triangular', 'Stellate', 'Bidirectional']
    ordered_types = []
    for base in base_types:
        for comp in ["Low", "High", "Unknown"]:
            t = f"{base}_{comp}"
            if t in cell_types:
                ordered_types.append(t)

    if not ordered_types:
        ordered_types = cell_types

    def _safe_name(s):
        return str(s).replace(' ', '_').replace('/', '_').replace('\\', '_')

    for cell_type in ordered_types:
        subdf = df[df['cell_type_with_complexity'] == cell_type].copy()
        if len(subdf) == 0:
            print(f"[INFO] No cells found for type: {cell_type}")
            continue

        # stable ordering across runs
        if "Label" in subdf.columns:
            subdf = subdf.sort_values("Label", kind="mergesort")
        else:
            subdf = subdf.sort_values("CellKey", kind="mergesort")

        n_total = len(subdf)
        n_pages = int(np.ceil(n_total / float(cells_per_page)))
        print(f"[INFO] {cell_type}: {n_total} cells -> {n_pages} pages ({cells_per_page}/page)")

        for page_idx in range(n_pages):
            start = page_idx * cells_per_page
            end = min((page_idx + 1) * cells_per_page, n_total)
            page_df = subdf.iloc[start:end].copy()

            n = len(page_df)
            ncols = int(max(1, ncols_cells))
            nrows = int(np.ceil(n / float(ncols)))

            # Figure sizing: wide (fits 2 columns), height scales with rows
            fig_w = 22 if ncols <= 2 else 22 + (ncols - 2) * 10
            fig_h = max(6, nrows * 3.2)
            fig = plt.figure(figsize=(fig_w, fig_h), dpi=150)

            # Outer grid = cell blocks arranged in nrows x ncols
            outer = GridSpec(
                nrows, ncols,
                figure=fig,
                hspace=0.35,
                wspace=0.25
            )

            for k, (_, row) in enumerate(page_df.iterrows()):
                rr = k // ncols
                cc = k % ncols

                # Inner layout per cell block:
                #   [ morphology | polar ]
                #   [  data spans both  ]
                inner = GridSpecFromSubplotSpec(
                    2, 2,
                    subplot_spec=outer[rr, cc],
                    height_ratios=[1.0, 0.62],
                    width_ratios=[1.55, 1.0],
                    hspace=0.05,
                    wspace=0.08
                )

                ax_morph = fig.add_subplot(inner[0, 0])
                ax_polar = fig.add_subplot(inner[0, 1], projection='polar')
                ax_data  = fig.add_subplot(inner[1, :])
                ax_data.axis('off')

                nodes = load_swc(row['SWC_Path'])
                center = row['center']
                c = get_cell_type_color(cell_type)

                # ---- morphology ----
                for nid, nd in nodes.items():
                    pid = nd["parent"]
                    if pid in nodes:
                        ax_morph.plot([nodes[pid]["x"], nd["x"]],
                                      [nodes[pid]["y"], nd["y"]],
                                      'k-', linewidth=LINE_WIDTH, alpha=0.45)

                ax_morph.plot(center[0], center[1], 'o', color='#E63946',
                              markersize=7, markeredgecolor='white', markeredgewidth=1.2)

                for pinfo, ray in zip(row['analysis_points'], row['rays']):
                    pos = pinfo['position']
                    ax_morph.plot([center[0], pos[0]], [center[1], pos[1]],
                                  color=c, linewidth=VECTOR_WIDTH, alpha=0.85)
                    mkr = 's' if pinfo['is_branch'] else 'o'
                    ax_morph.plot(pos[0], pos[1], mkr, color=c,
                                  markersize=5.5, markeredgecolor='white', markeredgewidth=0.9)

                ax_morph.set_aspect('equal')
                ax_morph.axis('off')

                label = row.get("Label", "")
                nprim = row.get("n_primary", np.nan)
                ax_morph.set_title(f"{label}\n(n={nprim})", fontsize=8.5, color=c, pad=3)

                # ---- polar ----
                for ray in row['rays']:
                    ang = math.degrees(math.atan2(ray[1], ray[0]))
                    if ang < 0:
                        ang += 360
                    theta = np.radians(ang)
                    ax_polar.plot([0, theta], [0, 1], color=c, linewidth=2.0, alpha=0.9)
                    ax_polar.plot(theta, 1, 'o', color=c,
                                  markersize=5, markeredgecolor='white', markeredgewidth=0.8)

                ax_polar.set_ylim(0, 1.2)
                ax_polar.set_theta_zero_location('E')
                ax_polar.grid(True, alpha=0.25, linestyle=':')

                L1 = row.get('largest_gap', np.nan)
                L2 = row.get('second_gap', np.nan)
                ax_polar.set_title(f"L1={L1:.0f}°\nL2={L2:.0f}°", fontsize=8, pad=2)

                # ---- data ----
                R = (L2 / L1) if (pd.notna(L1) and pd.notna(L2) and L1 > 0) else np.nan
                gaps = row.get("all_gaps", [])
                gap_list_str = ", ".join(f"{g:.0f}" for g in gaps) if gaps else "N/A"

                t100 = row.get("T100", np.nan)
                t100_str = f"{t100:.1f}" if pd.notna(t100) else "N/A"

                info = (
                    f"L1={L1:.1f}°   L2={L2:.1f}°   R={R:.2f}   T100={t100_str}   "
                    f"Complexity={row.get('complexity','')}\n"
                    f"Gaps: {gap_list_str}"
                )
                ax_data.text(0.0, 0.95, info, va='top', ha='left', fontsize=7.5, family='monospace')

            safe = _safe_name(cell_type)
            out_path = os.path.join(out_dir, f"class_{safe}_p{page_idx+1:02d}_of_{n_pages:02d}.svg")

            fig.suptitle(
                f"{cell_type} — Cells {start+1}-{end} of {n_total}  (Page {page_idx+1}/{n_pages})",
                fontsize=16, fontweight='bold', y=0.995
            )

            fig.tight_layout(rect=[0, 0, 1, 0.98])
            fig.savefig(out_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"[Saved] Per-class page: {out_path}")

# ==================== CSV OUTPUT ====================
def save_detailed_csv(df, out_path):
    """Save detailed analysis results to CSV with standard nomenclature."""
    output_rows = []

    for _, row in df.iterrows():
        # Ratio for multi-primary cells
        ratio = np.nan
        if row['n_primary'] >= 3 and not pd.isna(row['largest_gap']) and not pd.isna(row['second_gap']):
            ratio = row['second_gap'] / row['largest_gap'] if row['largest_gap'] > 0 else np.nan

        # --- T100 column (you compute T100 as T100) ---
        t100_val = row.get("T100", np.nan)

        base_row = {
            # IDs
            "CellKey": row.get("CellKey", ""),
            "Label": row.get("Label", ""),
            "Clone": row.get("CloneKey", ""),   # <-- IMPORTANT: your histogram code expects "Clone"
            "CloneKey": row.get("CloneKey", ""),

            # Nomenclature requested
            "Geometry": row.get("cell_type", ""),                   # base type
            "Morphotype": row.get("cell_type_with_complexity", ""), # type+complexity

            # Keep old names too (optional but useful)
            "CellType": row.get("cell_type", ""),
            "CellType_With_Complexity": row.get("cell_type_with_complexity", ""),
            "Complexity": row.get("complexity", ""),

            # Core metrics
            "n_primary": row.get("n_primary", np.nan),
            "Largest_Gap_deg": row.get("largest_gap", np.nan),
            "Second_Gap_deg": row.get("second_gap", np.nan),
            "L2_L1_Ratio": ratio,
            "T100": t100_val,
        }

        # Add all pairwise angles
        for idx, angle in enumerate(row.get('pairwise_angles', []), 1):
            base_row[f'Pairwise_Angle_{idx}_deg'] = angle

        output_rows.append(base_row)

    output_df = pd.DataFrame(output_rows)
    output_df.to_csv(out_path, index=False)
    print(f"[Saved] Detailed CSV: {out_path}")

# ==================== MAIN PIPELINE ====================
def main():
    # Load data
    cells_df = load_cells_metadata()
    if cells_df.empty:
        print("ERROR: No cells loaded from CSV")
        return

    # Build angle analysis table
    analysis_df = build_angle_analysis_table(cells_df)
    if analysis_df.empty:
        print("ERROR: No cells successfully analyzed")
        return

    # ADD COMPLEXITY CLASSIFICATION
    analysis_df = add_complexity_classification(analysis_df)

    # Save detailed CSV with all pairwise angles and complexity
    save_detailed_csv(
        analysis_df,
        os.path.join(OUT_DIR, 'cell_analysis_detailed_with_complexity.csv')
    )

    # Generate all visualizations
    print("\n" + "="*80)
    print("GENERATING ALL VISUALIZATIONS (WITH COMPLEXITY)")
    print("="*80)

    # 1. Method schematic
    print("\n[1/7] Generating method schematic...")
    plot_method_schematic(
        analysis_df,
        os.path.join(OUT_DIR, '1_method_schematic_with_complexity.svg'),
        pick_cellkey="S058_a6|s4|clone_10|21"  # example
    )

    # 2. Comprehensive overview
    print("\n[2/7] Generating comprehensive overview...")
    plot_comprehensive_overview(
        analysis_df,
        os.path.join(OUT_DIR, '2_comprehensive_overview_with_complexity.svg')
    )

    # 3. Cell type distribution
    print("\n[3/7] Generating cell type distribution...")
    plot_cell_type_distribution(
        analysis_df,
        os.path.join(OUT_DIR, '3_cell_type_distribution_with_complexity.svg')
    )

    # 4. Clone composition
    print("\n[4/7] Generating clone composition...")
    plot_clone_composition(
        analysis_df,
        os.path.join(OUT_DIR, '4_clone_composition_with_complexity.svg')
    )

    # 5. Empirical distributions
    print("\n[5/7] Generating empirical distributions...")
    plot_empirical_distributions(
        analysis_df,
        os.path.join(OUT_DIR, '5_empirical_distributions.svg')
    )

    # 6. Complexity analysis
    print("\n[6/7] Generating complexity analysis...")
    plot_complexity_analysis(
        analysis_df,
        os.path.join(OUT_DIR, '6_complexity_analysis.svg')
    )

    # 7. PER-CLASS OVERVIEWS
    print("\n[7/7] Generating per-class overviews...")
    print(f"Output directory: {PER_CLASS_DIR}")
    plot_one_figure_per_class(analysis_df, PER_CLASS_DIR)

    # Print summary statistics
    print("\n" + "="*80)
    print("FINAL ANALYSIS SUMMARY (IMPROVED EMPIRICAL CRITERIA WITH COMPLEXITY)")
    print("="*80)
    print(f"\nTotal cells analyzed: {len(analysis_df)}")
    print(f"Total clones: {analysis_df['CloneKey'].nunique()}")
    print(
        "Clones with multiple cells: "
        f"{analysis_df.groupby('CloneKey').size()[analysis_df.groupby('CloneKey').size() > 1].count()}"
    )

    print(f"\nCell Type Distribution (with complexity):")
    base_types = ['Unidirectional', 'Bipolar', 'Triangular', 'Stellate', 'Bidirectional']
    for base_type in base_types:
        for complexity in ['Low', 'High']:
            full_type = f"{base_type}_{complexity}"
            count = len(analysis_df[analysis_df['cell_type_with_complexity'] == full_type])
            if count > 0:
                pct = 100 * count / len(analysis_df)
                print(f"  {full_type:25s}: {count:3d} ({pct:5.1f}%)")

    # Calculate T100 statistics
    t100_col = None
    for col in analysis_df.columns:
        if 'T100' in str(col):
            t100_col = col
            break

    if t100_col:
        print(f"\nT100 Statistics by Complexity:")
        for complexity in ['Low', 'High']:
            comp_data = analysis_df[analysis_df['complexity'] == complexity][t100_col].dropna()
            if len(comp_data) > 0:
                print(f"  {complexity:8s}: Mean = {comp_data.mean():.1f}, "
                      f"Std = {comp_data.std():.1f}, Range = {comp_data.min():.0f}-{comp_data.max():.0f}")

    # Show primary count statistics
    print(f"\nPrimary Branches Statistics:")
    print(f"  Mean: {analysis_df['n_primary'].mean():.2f} ± {analysis_df['n_primary'].std():.2f}")
    print(f"  Range: {analysis_df['n_primary'].min():.0f} - {analysis_df['n_primary'].max():.0f}")
    print(f"  Distribution:")
    for n in sorted(analysis_df['n_primary'].unique()):
        count = (analysis_df['n_primary'] == n).sum()
        pct = 100 * count / len(analysis_df)
        print(f"    {n} primaries: {count:3d} ({pct:5.1f}%)")

    print(f"\nAll outputs saved to: {OUT_DIR}")
    print(f"Per-class overviews saved to: {PER_CLASS_DIR}")
    print("="*80)

if __name__ == "__main__":
    main()
