"""
================================================================================
UNIFIED NEURONAL ANGLE ANALYSIS - FINAL CELL TYPE CLASSIFICATION WITH COMPLEXITY
================================================================================
Complete angle analysis pipeline with FINAL cell type classification:
- Analysis at first branching point OR 50μm (whichever comes first)
- FINAL Classification into 5 morphological types (no hemidirectional)
- All cells forced into closest matching type
- Enhanced visualizations with thicker lines and multiple examples
- INCLUDES per-class overview folder
- COMPLEXITY (Low/High) computed DIRECTLY from SWC via Sholl T100 (no Excel merge)
================================================================================
"""

import os, re, glob, math, warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import seaborn as sns

########################################################################################################################################################################################################
# User input
CSV_PATH = r"C:\Users\fdasilve\Morphological Analysis\imaris_cellinfo_with_total_length_with_map__len_gt400.csv" # Path of the Cell_Info csv file
SWC_ROOT = r"C:\Users\fdasilve\Morphological Analysis\Our_data_final_version (25112025)" # Path of folder containing clone folders
OUT_DIR = r"C:\Users\fdasilve\Morphological Analysis\Morphotype (Sholl+Angle)" # Desired output folder
PER_CLASS_DIR = os.path.join(OUT_DIR, "PER_CLASS_OVERVIEWS")

# Angle analysis parameters
TARGET_R = 50.0  # Analysis at 50μm or first branch point

# Sholl parameters (computed directly from SWC)
SHOLL_STEP = 2.0
SHOLL_R_MAX = 500.0
# Complexity defined from T100 = total Sholl intersections from 0..100µm (inclusive)
T_TOTAL_RADII = (100,)  # you can add more: (20,30,40,60,90,100,120,150)

# FINAL Cell type classification criteria (your data-driven values)
UNIDIRECTIONAL_L1_THRESHOLD = 300.0

BIPOLAR_MIN = 130.0
BIPOLAR_MAX = 230.0

BIDIRECTIONAL_L1_MIN = 160.0
BIDIRECTIONAL_L1_MAX = 299.0
BIDIRECTIONAL_RATIO_MAX = 0.99

STELLATE_L1_MIN = 45.0
STELLATE_L1_MAX = 159.9
STELLATE_RATIO_MIN = 0.33

TRIANGULAR_MIN = 75.0
TRIANGULAR_MAX = 145.0

# Visualization settings
DPI = 150
LINE_WIDTH = 1.2
VECTOR_WIDTH = 3.5

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial"]
plt.rcParams["axes.linewidth"] = 1.5
warnings.filterwarnings("ignore")

os.makedirs(OUT_DIR, exist_ok=True)
os.makedirs(PER_CLASS_DIR, exist_ok=True)

plt.style.use("default")
sns.set_context("paper")

print("\n" + "=" * 80)
print("UNIFIED ANGLE ANALYSIS - FINAL CELL TYPE CLASSIFICATION WITH COMPLEXITY (T100 FROM SWC)")
print("=" * 80)

########################################################################################################################################################################################################
# Utility functions
def norm(s):
    return re.sub(r"[^a-z0-9]+", "_", str(s).strip().lower())

def find_col(df, must, exclude=None):
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
    m = re.search(r"\d+", ss)
    return "c" + m.group(0) if m else ss.lower()

def id_from_text(s):
    m = re.search(r"\d+", str(s))
    return m.group(0) if m else ""

# SWC file parsing
def parse_swc_hierarchy(path):
    parts = re.split(r"[\\/]+", os.path.normpath(path))
    litter = animal = section = clone = None

    for p in parts:
        if re.match(r"^[sS]\d+$", p):
            litter = p.upper()
            break

    for p in parts:
        if re.match(r"^[sS]\d+_[aA]\d+", p, re.IGNORECASE):
            animal = p
            break

    for i, p in enumerate(parts):
        if animal and animal in parts and i > parts.index(animal):
            if re.match(r"^[sS]\d+$", p):
                section = p.lower()
                break

    for p in parts:
        if "clone" in p.lower():
            clone = p
            break

    return litter or "", animal or "", section or "", clone or ""

def create_cell_key(path):
    litter, animal, section, clone = parse_swc_hierarchy(path)
    ani_full = animal if animal else "unknown"
    sec = canon_section(section)
    cln = canon_clone(clone)
    cid = id_from_text(os.path.splitext(os.path.basename(path))[0])

    clone_key = f"{ani_full}-{sec}-{cln}"
    cell_label = f"{ani_full}-{sec}-{cln}-{cid}"
    cell_key = f"{ani_full}|{sec}|{cln}|{cid}"
    return cell_key, cell_label, clone_key

# CSV Processing
def load_cell_metadata(csv_path):
    print(f"\nLoading cell metadata from CSV (dendrites per line)...")
    df = pd.read_csv(csv_path, dtype=str)

    LIT = find_col(df, ["litter"]) or "Litter"
    ANI = find_col(df, ["animal"]) or "Animal"
    SEC = find_col(df, ["segment"]) or find_col(df, ["section"]) or "Segment"
    ID_COL = find_col(df, ["id"], exclude=["project"]) or "ID"
    CLONE_COL = find_col(df, ["clone"]) or "Clone"

    cells = {}
    for _, row in df.iterrows():
        lit = canon_litter(row.get(LIT, ""))
        ani = canon_animal(row.get(ANI, ""))
        sec = canon_section(row.get(SEC, ""))
        cln = canon_clone(row.get(CLONE_COL, ""))
        cid = id_from_text(row.get(ID_COL, ""))

        ani_full = f"{lit}_{ani}" if lit and ani else (lit or ani or "unknown")
        cell_key = f"{ani_full}|{sec}|{cln}|{cid}"

        if cell_key not in cells:
            cells[cell_key] = {
                "CellKey": cell_key,
                "Label": f"{ani_full}-{sec}-{cln}-{cid}",
                "CloneKey": f"{ani_full}-{sec}-{cln}",
                "Litter": row.get(LIT, ""),
                "Animal": row.get(ANI, ""),
                "Section": row.get(SEC, ""),
                "Clone": row.get(CLONE_COL, ""),
                "CellID": row.get(ID_COL, "")
            }

    cells_df = pd.DataFrame(list(cells.values()))
    print(f"Found {len(cells_df)} unique cells in CSV")
    return cells_df

# SWC IO
def load_swc(path):
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
    return [i for i, n in nodes.items() if n["parent"] in soma_ids]

def dist_3d(n, center):
    dx = n["x"] - center[0]
    dy = n["y"] - center[1]
    dz = n["z"] - center[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)

# Sholl analysis
def dist_to_center(nd, c):
    dx, dy, dz = nd["x"] - c[0], nd["y"] - c[1], nd["z"] - c[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)

def sholl_counts(nodes, center, step=2.0, r_max=500.0, exclude_types={2}):
    radii = np.arange(0.0, r_max + step, step, dtype=float)
    counts = np.zeros_like(radii, dtype=int)
    eps = 1e-9

    for nid, nd in nodes.items():
        pid = nd["parent"]
        if pid not in nodes:
            continue

        p = nodes[pid]

        if (p["type"] in exclude_types) or (nd["type"] in exclude_types):
            continue

        dp = dist_to_center(p, center)
        dc = dist_to_center(nd, center)
        lo, hi = (dp, dc) if dp <= dc else (dc, dp)

        i0 = int(max(0, math.floor(lo / step)))
        i1 = int(min(len(radii) - 1, math.ceil(hi / step)))

        for k in range(i0, i1 + 1):
            r = radii[k]
            dpp, dcc = dp - r, dc - r

            if abs(dpp) < eps:
                dpp = eps if dpp >= 0 else -eps
            if abs(dcc) < eps:
                dcc = eps if dcc >= 0 else -eps

            if (dpp < 0 and dcc > 0) or (dpp > 0 and dcc < 0):
                counts[k] += 1

    return radii, counts

def compute_sholl_metrics_from_swc(nodes, center, step=SHOLL_STEP, r_max=SHOLL_R_MAX, radii_totals=T_TOTAL_RADII):
    r, c = sholl_counts(nodes, center, step=step, r_max=r_max)

    totals = {}
    for xx in radii_totals:
        idx = int(round(float(xx) / float(step)))
        idx = max(0, min(idx, len(c) - 1))
        totals[f"T{int(xx)}"] = int(c[: idx + 1].sum())

    MAX = int(np.max(c)) if len(c) else 0
    R_AT_MAX = float(r[int(np.argmax(c))]) if len(c) else np.nan

    return {
        **totals,
        "SHOLL_MAX": MAX,
        "SHOLL_R_AT_MAX": R_AT_MAX,
    }

# Angle analysis
def project_to_radius_or_branch(nodes, center, start_id, target_r=TARGET_R):
    cur = start_id
    visited = set()

    for _ in range(10000):
        visited.add(cur)
        cur_dist = dist_3d(nodes[cur], center)

        # radius reached
        if cur_dist >= target_r:
            best = cur
            best_diff = abs(cur_dist - target_r)

            if nodes[cur]["parent"] in nodes:
                parent = nodes[cur]["parent"]
                parent_dist = dist_3d(nodes[parent], center)
                parent_diff = abs(parent_dist - target_r)
                if parent_diff < best_diff:
                    best = parent

            return best, dist_3d(nodes[best], center), False

        # branch check
        children = [i for i, n in nodes.items() if n["parent"] == cur and i not in visited]
        if len(children) >= 2:
            return cur, cur_dist, True
        elif len(children) == 1:
            cur = children[0]
        else:
            return cur, cur_dist, False

    return cur, dist_3d(nodes[cur], center), False

def compute_pairwise_angles_with_details(rays):
    angle_details = []
    m = len(rays)
    for i in range(m):
        for j in range(i + 1, m):
            dot = np.clip(np.dot(rays[i], rays[j]), -1.0, 1.0)
            ang = math.degrees(math.acos(dot))
            if ang > 180:
                ang = 360 - ang
            angle_details.append((ang, i, j))
    return angle_details

def analyze_cell_angles(nodes):
    center, soma_ids = soma_center_and_ids(nodes)
    prims = primary_children(nodes, soma_ids)
    if not prims:
        return None

    rays = []
    analysis_points = []

    for p in prims:
        node_id, dist, is_branch = project_to_radius_or_branch(nodes, center, p, TARGET_R)
        nd = nodes[node_id]

        v = np.array([nd["x"] - center[0], nd["y"] - center[1], nd["z"] - center[2]], dtype=float)
        nv = np.linalg.norm(v)
        if nv > 0:
            rays.append(v / nv)
            analysis_points.append({
                "node_id": node_id,
                "distance": dist,
                "is_branch": is_branch,
                "position": (nd["x"], nd["y"], nd["z"])
            })

    if len(rays) < 1:
        return None

    angle_details = compute_pairwise_angles_with_details(rays)
    angles = [a[0] for a in angle_details]

    # gaps from 2D angles
    if len(rays) >= 2:
        angles_2d = []
        for ray in rays:
            a2 = math.degrees(math.atan2(ray[1], ray[0]))
            if a2 < 0:
                a2 += 360
            angles_2d.append(a2)

        ang_sorted = sorted(angles_2d)
        gaps = []
        for i in range(len(ang_sorted)):
            ni = (i + 1) % len(ang_sorted)
            gap = ang_sorted[ni] - ang_sorted[i]
            if gap < 0:
                gap += 360
            gaps.append(gap)
        gaps_sorted = sorted(gaps, reverse=True)
    else:
        gaps_sorted = []

    return {
        "rays": rays,
        "angles": angles,
        "angle_details": angle_details,
        "n_primary": len(prims),
        "n_pairs": len(angles),
        "center": center,
        "analysis_points": analysis_points,
        "gaps": gaps_sorted
    }

def classify_cell_type(n_primary, gaps):
    if len(gaps) == 0:
        return "Unidirectional" if n_primary == 1 else "Stellate"

    L1 = gaps[0] if len(gaps) > 0 else 0

    # UNIDIRECTIONAL
    if n_primary == 1 or (n_primary > 1 and L1 > UNIDIRECTIONAL_L1_THRESHOLD):
        return "Unidirectional"

    # BIPOLAR
    if n_primary == 2:
        return "Bipolar"

    # TRIANGULAR
    if n_primary == 3:
        min_gap = min(gaps) if len(gaps) > 0 else 0
        if min_gap > 50.0:
            return "Triangular"

    # ≥3 primaries: Stellate vs Bidirectional
    if n_primary >= 3:
        if len(gaps) >= 2:
            L1 = gaps[0]
            L2 = gaps[1]
            R = L2 / L1 if L1 > 0 else 0
        else:
            L1 = gaps[0] if len(gaps) > 0 else 0
            L2 = 0
            R = 0

        # Stellate (clean)
        if n_primary >= 4 and L1 <= STELLATE_L1_MAX and R >= STELLATE_RATIO_MIN:
            return "Stellate"

        # Bidirectional (clean)
        if L1 >= BIDIRECTIONAL_L1_MIN and R <= BIDIRECTIONAL_RATIO_MAX:
            return "Bidirectional"

        # Intermediate zone
        if 135.0 <= L1 <= 165.0:
            if R >= 0.85:
                return "Stellate"
            elif R <= 0.65:
                return "Bidirectional"

        # distance fallback
        stellate_distance = 0.0
        if n_primary >= 4:
            stellate_distance += 0
        else:
            stellate_distance += 20

        if L1 <= STELLATE_L1_MAX:
            stellate_distance += (STELLATE_L1_MAX - L1) / 10.0
        else:
            stellate_distance += (L1 - STELLATE_L1_MAX) * 2.0

        if R >= STELLATE_RATIO_MIN:
            stellate_distance += (R - STELLATE_RATIO_MIN) * 5.0
        else:
            stellate_distance += (STELLATE_RATIO_MIN - R) * 10.0

        bidirectional_distance = 0.0
        if L1 >= BIDIRECTIONAL_L1_MIN:
            bidirectional_distance += (L1 - BIDIRECTIONAL_L1_MIN) / 10.0
        else:
            bidirectional_distance += (BIDIRECTIONAL_L1_MIN - L1) * 2.0

        if R <= BIDIRECTIONAL_RATIO_MAX:
            bidirectional_distance += (BIDIRECTIONAL_RATIO_MAX - R) * 5.0
        else:
            bidirectional_distance += (R - BIDIRECTIONAL_RATIO_MAX) * 10.0

        if n_primary == 3:
            stellate_distance += 10

        return "Stellate" if stellate_distance <= bidirectional_distance else "Bidirectional"

    return "Stellate"

# Complexity
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

# Analysis table
def build_angle_analysis_table(cells_df):
    print(f"\nScanning for SWC files...")
    swc_pattern = os.path.join(SWC_ROOT, "**", "*.swc")
    swc_files = glob.glob(swc_pattern, recursive=True)
    print(f"Found {len(swc_files)} SWC files")

    if len(swc_files) == 0:
        print(f"\nERROR: No SWC files found in {SWC_ROOT}")
        return pd.DataFrame()

    swc_lookup = {}
    for swc_path in swc_files:
        cell_key, cell_label, clone_key = create_cell_key(swc_path)
        swc_lookup[cell_key] = swc_path

    print(f"Indexed {len(swc_lookup)} unique cells from SWC files")

    results = []
    for _, cell in cells_df.iterrows():
        cell_key = cell["CellKey"]
        if cell_key not in swc_lookup:
            continue

        swc_path = swc_lookup[cell_key]

        try:
            nodes = load_swc(swc_path)
            if len(nodes) < 10:
                continue

            analysis = analyze_cell_angles(nodes)
            if analysis is None:
                continue

            cell_type = classify_cell_type(analysis["n_primary"], analysis["gaps"])

            # Sholl metrics (T100 etc.) directly from SWC
            center_for_sholl, _ = soma_center_and_ids(nodes)
            sholl_metrics = compute_sholl_metrics_from_swc(
                nodes,
                center_for_sholl,
                step=SHOLL_STEP,
                r_max=SHOLL_R_MAX,
                radii_totals=T_TOTAL_RADII
            )

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

                # Sholl outputs
                **sholl_metrics
            }
            results.append(row)

        except Exception as e:
            print(f"Error processing {swc_path}: {e}")

    if not results:
        print("\nNo cells successfully processed.")
        return pd.DataFrame()

    df = pd.DataFrame(results)
    print(f"\nSuccessfully processed {len(df)} cells")
    return df

# helpers
def get_cell_type_color(cell_type):
    base_type = cell_type
    if "_Low" in cell_type:
        base_type = cell_type.replace("_Low", "")
    elif "_High" in cell_type:
        base_type = cell_type.replace("_High", "")

    base_colors = {
        "Unidirectional": "#E63946",
        "Bipolar": "#F77F00",
        "Bidirectional": "#06D6A0",
        "Triangular": "#118AB2",
        "Stellate": "#8338EC",
    }

    base_color = base_colors.get(base_type, "#999999")

    import colorsys
    r = int(base_color[1:3], 16) / 255.0
    g = int(base_color[3:5], 16) / 255.0
    b = int(base_color[5:7], 16) / 255.0
    h, l, s = colorsys.rgb_to_hls(r, g, b)

    if "_Low" in cell_type:
        l = min(1.0, l * 1.4)
    elif "_High" in cell_type:
        l = max(0.0, l * 0.6)

    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return f"#{int(r * 255):02x}{int(g * 255):02x}{int(b * 255):02x}"

# Output CSVs
def save_detailed_csv(df, out_path):
    output_rows = []

    for _, row in df.iterrows():
        ratio = np.nan
        if row["n_primary"] >= 3 and not pd.isna(row["largest_gap"]) and not pd.isna(row["second_gap"]):
            ratio = row["second_gap"] / row["largest_gap"] if row["largest_gap"] > 0 else np.nan

        base_row = {
            "CellKey": row["CellKey"],
            "Label": row["Label"],
            "CloneKey": row["CloneKey"],
            "CellType": row["cell_type"],
            "CellType_With_Complexity": row.get("cell_type_with_complexity", row["cell_type"] + "_Unknown"),
            "Complexity": row.get("complexity", "Unknown"),
            "N_Primary_Branches": row["n_primary"],
            "Largest_Gap_deg": row["largest_gap"],
            "Second_Gap_deg": row["second_gap"],
            "L2_L1_Ratio": ratio,
            "T100": row.get("T100", np.nan),
            "SHOLL_MAX": row.get("SHOLL_MAX", np.nan),
            "SHOLL_R_AT_MAX": row.get("SHOLL_R_AT_MAX", np.nan),
        }

        for idx, angle in enumerate(row["pairwise_angles"], 1):
            base_row[f"Pairwise_Angle_{idx}_deg"] = angle

        output_rows.append(base_row)

    pd.DataFrame(output_rows).to_csv(out_path, index=False)
    print(f"[Saved] Detailed CSV: {out_path}")

# Simple plots
def plot_cell_type_distribution(df, out_path):
    fig, ax = plt.subplots(1, 1, figsize=(14, 6), dpi=DPI)

    type_counts = df["cell_type_with_complexity"].value_counts()

    base_types = ["Unidirectional", "Bipolar", "Triangular", "Stellate", "Bidirectional"]
    ordered = []
    for b in base_types:
        for comp in ["Low", "High", "Unknown"]:
            t = f"{b}_{comp}"
            if t in type_counts.index:
                ordered.append(t)

    type_counts = type_counts.reindex(ordered, fill_value=0)
    colors = [get_cell_type_color(ct) for ct in type_counts.index]

    bars = ax.bar(range(len(type_counts)), type_counts.values, color=colors, edgecolor="black", linewidth=1.2, alpha=0.85)
    ax.set_xticks(range(len(type_counts)))
    ax.set_xticklabels(type_counts.index, rotation=45, ha="right", fontsize=10)
    ax.set_ylabel("Number of Cells", fontsize=12, fontweight="bold")
    ax.set_title("Cell Type Distribution (With Complexity from T100 computed in SWC)", fontsize=14, fontweight="bold", pad=12)
    ax.grid(True, alpha=0.25, axis="y", linestyle=":")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    for bar in bars:
        h = bar.get_height()
        if h > 0:
            ax.text(bar.get_x() + bar.get_width() / 2.0, h, f"{int(h)}", ha="center", va="bottom", fontsize=9, fontweight="bold")

    plt.tight_layout()
    plt.savefig(out_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[Saved] Cell type distribution: {out_path}")

# ==================== MAIN PIPELINE ====================
def main():
    cells_df = load_cell_metadata(CSV_PATH)
    if cells_df.empty:
        print("ERROR: No cells loaded from CSV")
        return

    analysis_df = build_angle_analysis_table(cells_df)
    if analysis_df.empty:
        print("ERROR: No cells successfully analyzed")
        return

    analysis_df = add_complexity_classification(analysis_df)

    save_detailed_csv(
        analysis_df,
        os.path.join(OUT_DIR, "cell_analysis_detailed_with_complexity_T100_from_SWC.csv")
    )

    plot_cell_type_distribution(
        analysis_df,
        os.path.join(OUT_DIR, "cell_type_distribution_with_complexity_T100_from_SWC.svg")
    )

    print("\n" + "=" * 80)
    print("FINAL ANALYSIS SUMMARY (T100 FROM SWC)")
    print("=" * 80)
    print(f"\nTotal cells analyzed: {len(analysis_df)}")
    print(f"Total clones: {analysis_df['CloneKey'].nunique()}")
    print(f"\nT100 overall: mean={analysis_df['T100'].mean():.2f}  std={analysis_df['T100'].std():.2f}")

    print(f"\nCell Type Distribution (with complexity):")
    base_types = ["Unidirectional", "Bipolar", "Triangular", "Stellate", "Bidirectional"]
    for base_type in base_types:
        for complexity in ["Low", "High", "Unknown"]:
            full_type = f"{base_type}_{complexity}"
            count = len(analysis_df[analysis_df["cell_type_with_complexity"] == full_type])
            if count > 0:
                pct = 100 * count / len(analysis_df)
                print(f"  {full_type:28s}: {count:3d} ({pct:5.1f}%)")

    print(f"\nAll outputs saved to: {OUT_DIR}")
    print(f"Per-class overviews dir exists (not generated in this minimal build): {PER_CLASS_DIR}")
    print("=" * 80)

if __name__ == "__main__":
    main()
