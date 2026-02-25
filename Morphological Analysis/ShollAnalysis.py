# ==========================================================
#   SHOLL (paper-like) — per-cell & per-clone curves, Updated 26112025
#   totals (0–r), MAX intersections, bar plots,
#   and CSV report with full Sholl curves.
#
#   Folder hierarchy:
#     .../<Litter S###>/<Animal S###_a#>/<Section s1..s6>/<Clone clone*>/<cell>.swc
#   Example:
#     C:\...\our_swc_files_corrected\S127\S127_a20\S1\clone 1\3.swc
#   CloneKey = Litter-Animal-Section-Clone
#            = S127-S127_a20-S1-clone 1
# ==========================================================

import os, re, glob, math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Patch

plt.style.use("default")
sns.set_context("notebook")
DPI = 110

# ----------------------------
# CONFIG — edit these
# ----------------------------
SWC_ROOT = r"C:\Users\fdasilve\PycharmProjects\PythonProject\morphological_analysis_sc\our_swc_files_corrected"
SWC_GLOB = "**/*.swc"

STEP  = 2.0   # Sholl shell step (µm)
R_MAX = 500.0 # Max radius (µm)

# Radii used to compute per-radius intercepts (Ixx) and totals (Txx)
RADII = [20, 30, 40, 60, 90, 100, 120, 150]

# How many “big clones” to show in the example per-cell plots
N_TOP_CLONES_EXAMPLES = 3

# Output directory for SVG plots + CSV
OUT_DIR = r"C:\Users\fdasilve\Morphological Analysis\Sholl_analysis_per_cell (26112025)"
os.makedirs(OUT_DIR, exist_ok=True)

# ------------------ SWC parsing ------------------
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
                nid = int(p[0])
                ntype = int(float(p[1]))
                x, y, z = float(p[2]), float(p[3]), float(p[4])
                r = float(p[5])
                parent = int(p[6])
            except ValueError:
                continue
            nodes[nid] = dict(type=ntype, x=x, y=y, z=z, r=r, parent=parent)
    return nodes

def soma_center_and_ids(nodes):
    soma_ids = [i for i, n in nodes.items() if n["type"] == 1]
    if soma_ids:
        cx = float(np.mean([nodes[i]["x"] for i in soma_ids]))
        cy = float(np.mean([nodes[i]["y"] for i in soma_ids]))
        cz = float(np.mean([nodes[i]["z"] for i in soma_ids]))
        return (cx, cy, cz), set(soma_ids)
    roots = [i for i, n in nodes.items() if n["parent"] == -1]
    if roots:
        r0 = nodes[roots[0]]
        return (r0["x"], r0["y"], r0["z"]), set([roots[0]])
    any_id = next(iter(nodes))
    n = nodes[any_id]
    return (n["x"], n["y"], n["z"]), set([any_id])

def dist_to_center(n, c):
    dx, dy, dz = n["x"] - c[0], n["y"] - c[1], n["z"] - c[2]
    return math.sqrt(dx * dx + dy * dy + dz * dz)

def sholl_counts(nodes, center, step=2.0, r_max=500.0,
                 include_types={3, 4}, exclude_types={2}):
    radii = np.arange(0.0, r_max + step, step, dtype=float)
    counts = np.zeros_like(radii, dtype=int)
    eps = 1e-9
    for nid, nd in nodes.items():
        pid = nd["parent"]
        if pid in nodes:
            p = nodes[pid]
            if (p["type"] in exclude_types) or (nd["type"] in exclude_types):
                continue
            dp, dc = dist_to_center(p, center), dist_to_center(nd, center)
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

def count_primary_projections(nodes, soma_ids, include_types={3, 4}):
    soma_ids = set(soma_ids)
    return int(
        sum(
            1
            for nd in nodes.values()
            if nd["parent"] in soma_ids and nd["type"] in include_types
        )
    )

# ------------------ Folder HIERARCHY parsing (robust) ------------------
def _norm_parts(path):
    return re.split(r"[\\/]+", os.path.normpath(path))

def parse_hierarchy(path):
    """
    Parse by position from the file upward:
      .../<Litter S###>/<Animal S###_a#>/<Section s1..s6>/<Clone clone*>/<file>.swc

    Example:
      S127\S127_a20\S1\clone 1\3.swc
      -> Litter = S127
         Animal = S127_a20
         Section = S1
         Clone   = clone 1
    """
    parts = _norm_parts(path)
    dirs = parts[:-1]  # exclude filename

    # find last 'clone*' folder
    iclone = None
    for i in range(len(dirs) - 1, -1, -1):
        if dirs[i].lower().startswith("clone"):
            iclone = i
            break

    litter = animal = section = clone = None
    if iclone is not None:
        clone = dirs[iclone]

        # section: S1..S6 or s1..s6
        if iclone - 1 >= 0 and re.fullmatch(r"[sS][1-6]", dirs[iclone - 1]):
            section = dirs[iclone - 1]

        # animal: pattern like S###_a#, S127_a20, S058_a6 etc.
        if iclone - 2 >= 0:
            cand_animal = dirs[iclone - 2]
            if re.fullmatch(r"[sS]\d{2,4}_[aA]\d{1,3}", cand_animal):
                animal = cand_animal

        # litter: S### just before animal
        if iclone - 3 >= 0 and re.fullmatch(r"[sS]\d{2,4}", dirs[iclone - 3]):
            litter = dirs[iclone - 3].upper()

    litter  = litter  or "S?"
    animal  = animal  or "ANIMAL?"
    section = section or "S?"
    clone   = clone   or "clone ?"
    return litter, animal, section, clone

def make_clone_key(path):
    litter, animal, section, clone = parse_hierarchy(path)
    # keep original animal & section casing, litter forced to upper
    return f"{litter}-{animal}-{section}-{clone}"

def make_cell_id(path):
    return os.path.splitext(os.path.basename(path))[0]

# ------------------ Scan + compute ------------------
def scan_files(root):
    return sorted(set(glob.glob(os.path.join(root, SWC_GLOB), recursive=True)))

def compute_tables():
    files = scan_files(SWC_ROOT)
    if not files:
        raise RuntimeError("No SWC files found. Check SWC_ROOT.")
    all_r = sorted(set(RADII))
    common_r = np.arange(0.0, R_MAX + STEP, STEP)
    records = []
    curves_list = []   # list of dicts: {"ID", "CloneKey", "File", "curve"}

    for i, f in enumerate(files, 1):
        try:
            nodes = load_swc(f)
            if not nodes:
                continue
            center, soma_ids = soma_center_and_ids(nodes)
            r, c = sholl_counts(nodes, center, step=STEP, r_max=R_MAX)
            curve = np.interp(common_r, r, c, left=0, right=0)

            # per-radius & totals for all needed radii
            I = {f"I{int(x)}": int(curve[int(round(x / STEP))]) for x in all_r}
            T = {
                f"T{int(x)}": int(curve[: int(round(x / STEP)) + 1].sum())
                for x in all_r
            }

            # peak info
            MAX = int(curve.max())
            R_AT_MAX = float(common_r[int(np.argmax(curve))])

            n_primary = count_primary_projections(nodes, soma_ids)
            cid = make_cell_id(f)
            clone_key = make_clone_key(f)

            records.append(
                {
                    "ID": cid,
                    "File": f,
                    "CloneKey": clone_key,
                    **I,
                    **T,
                    "MAX": MAX,
                    "R_AT_MAX": R_AT_MAX,
                    "n_primary": n_primary,
                }
            )
            curves_list.append(
                {"ID": cid, "CloneKey": clone_key, "File": f, "curve": curve}
            )

            if i % 50 == 0:
                print(f"… {i}/{len(files)} processed")
        except Exception as e:
            print("Error:", f, e)

    return pd.DataFrame(records), curves_list, common_r

df, cell_curves, common_r = compute_tables()
if df.empty:
    raise RuntimeError("No SWC parsed into table.")

# ------------------ Aesthetic helpers ------------------
TAB20  = plt.get_cmap("tab20").colors
TAB20B = plt.get_cmap("tab20b").colors
TAB20C = plt.get_cmap("tab20c").colors
ALL_TAB_COLORS = list(TAB20) + list(TAB20B) + list(TAB20C)  # 60 distinct colors

def palette_from_names(names):
    names = list(dict.fromkeys(names))
    palette = {}
    for i, n in enumerate(names):
        palette[n] = ALL_TAB_COLORS[i % len(ALL_TAB_COLORS)]
    return palette

def tidy(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def safe_name(s):
    # for filenames
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", s)

clone_colors = palette_from_names(df["CloneKey"].unique())

# ==========================================================
# A) Sholl intersections per radius — per cell (all cells)
# ==========================================================
fig, ax = plt.subplots(figsize=(8.8, 5.1), dpi=DPI)
for entry in cell_curves:
    cid = entry["ID"]
    key = entry["CloneKey"]
    curve = entry["curve"]
    ax.plot(common_r, curve, lw=0.9, alpha=0.28, color=clone_colors[key])
ax.set_xlabel("Radius from soma (µm)")
ax.set_ylabel("Sholl intersections")
ax.set_title("Sholl intersections per radius — per cell")
tidy(ax)

# legend: ALL clones, outside axes
all_keys = list(pd.Series(df["CloneKey"]).value_counts().index)
handles = [
    Patch(facecolor=clone_colors[k], edgecolor="none",
          label=f"{k} (n={df['CloneKey'].eq(k).sum()})")
    for k in all_keys
]
ax.legend(handles=handles,
          frameon=True,
          fontsize=6,
          ncol=1 if False else 1,  # just to keep consistent structure
          loc="upper left",
          bbox_to_anchor=(1.02, 1.0),
          title="CloneKeys")
plt.tight_layout(rect=[0, 0, 0.8, 1])  # leave space on right
fig.savefig(os.path.join(OUT_DIR, "sholl_per_radius_all_cells.svg"),
            dpi=DPI, bbox_inches="tight")
plt.show()

# ==========================================================
# A2) Sholl intersections per radius — per cell
#     EXAMPLE: one figure per top clone (most cells)
# ==========================================================
clone_counts = df["CloneKey"].value_counts()
top_example_clones = list(clone_counts.head(N_TOP_CLONES_EXAMPLES).index)

for key in top_example_clones:
    fig, ax = plt.subplots(figsize=(8.8, 5.1), dpi=DPI)
    for entry in cell_curves:
        if entry["CloneKey"] != key:
            continue
        cid = entry["ID"]
        curve = entry["curve"]
        ax.plot(common_r, curve, lw=1.2, alpha=0.8, label=cid)
    ax.set_xlabel("Radius from soma (µm)")
    ax.set_ylabel("Sholl intersections")
    ax.set_title(
        f"Sholl intersections per radius — per cell\n"
        f"Clone: {key} (n={clone_counts[key]})"
    )
    tidy(ax)
    ax.legend(frameon=True, fontsize=7, ncol=1, loc="upper right", title="Cell ID")
    plt.tight_layout()
    fname = f"sholl_per_radius_cells_clone_{safe_name(key)}.svg"
    fig.savefig(os.path.join(OUT_DIR, fname), dpi=DPI, bbox_inches="tight")
    plt.show()

# ==========================================================
# B) Sholl intersections per radius — per CloneKey (mean ± SEM)
# ==========================================================
by_key = {}
for entry in cell_curves:
    key = entry["CloneKey"]
    curve = entry["curve"]
    by_key.setdefault(key, []).append(curve)

fig, ax = plt.subplots(figsize=(8.8, 5.1), dpi=DPI)
for key, curves in by_key.items():
    M = np.vstack(curves)
    mean = M.mean(axis=0)
    sem = M.std(axis=0, ddof=1) / np.sqrt(max(1, M.shape[0]))
    col = clone_colors[key]
    ax.plot(common_r, mean, lw=2.0, color=col, label=f"{key} (n={M.shape[0]})")
    ax.fill_between(common_r, mean - sem, mean + sem, color=col, alpha=0.18)
for x in RADII:
    ax.axvline(x, color="#e0a0a0", lw=0.6, alpha=0.35)
ax.set_xlabel("Radius from soma (µm)")
ax.set_ylabel("Sholl intersections")
ax.set_title("Sholl intersections per radius — per clone (mean ± SEM)")
tidy(ax)

handles = [
    Patch(facecolor=clone_colors[k], edgecolor="none",
          label=f"{k} (n={len(by_key[k])})")
    for k in all_keys
]
ax.legend(handles=handles,
          frameon=False,
          bbox_to_anchor=(1.02, 1.0),
          loc="upper left",
          fontsize=6,
          title="CloneKeys")
plt.tight_layout(rect=[0, 0, 0.8, 1])
fig.savefig(os.path.join(OUT_DIR, "sholl_per_radius_per_clone_mean_sem.svg"),
            dpi=DPI, bbox_inches="tight")
plt.show()

# ==========================================================
# C) Totals (0–r) — per CloneKey (group means)
# ==========================================================
T_cols = [f"T{int(r)}" for r in RADII]
g = df.groupby("CloneKey")[T_cols].mean().sort_index()
key_colors = palette_from_names(g.index)

fig, axes = plt.subplots(
    1,
    len(RADII),
    figsize=(max(12, 1.0 * len(g)), 5.2),
    dpi=DPI,
    sharey=True,
)
if len(RADII) == 1:
    axes = [axes]
for ax, r in zip(axes, RADII):
    col = f"T{int(r)}"
    cols = [key_colors[k] for k in g.index]
    ax.bar(
        np.arange(len(g)),
        g[col].values,
        color=cols,
        edgecolor="#222",
        linewidth=0.5,
    )
    ax.set_title(f"Total 0–{int(r)} µm", fontsize=11)
    ax.set_xticks(np.arange(len(g)))
    ax.set_xticklabels(g.index, rotation=90, fontsize=7)
    tidy(ax)
axes[0].set_ylabel("Total intersections (0–r)")
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "totals_0_r_per_clone.svg"),
            dpi=DPI, bbox_inches="tight")
plt.show()

# ==========================================================
# C2) Max intersections — per CloneKey (group means)
# ==========================================================
gmax = df.groupby("CloneKey")["MAX"].mean().sort_index()
key_colors = palette_from_names(gmax.index)

fig, ax = plt.subplots(
    figsize=(max(12, 1.0 * len(gmax)), 5.0),
    dpi=DPI,
)
ax.bar(
    np.arange(len(gmax)),
    gmax.values,
    color=[key_colors[k] for k in gmax.index],
    edgecolor="#222",
    linewidth=0.5,
)
ax.set_title("Max Sholl intersections — per CloneKey (group mean)", fontsize=11)
ax.set_xticks(np.arange(len(gmax)))
ax.set_xticklabels(gmax.index, rotation=90, fontsize=7)
ax.set_ylabel("Max intersections")
tidy(ax)
plt.tight_layout()
fig.savefig(os.path.join(OUT_DIR, "max_per_clone.svg"),
            dpi=DPI, bbox_inches="tight")
plt.show()

# ==========================================================
# D) Totals (0–r) — per cell (one figure per radius; bars colored by CloneKey)
#     X labels = CloneKey — ID
# ==========================================================
df_sorted = df.sort_values(["CloneKey", "ID"]).reset_index(drop=True)
vals = df_sorted["CloneKey"].values
cuts = [i for i in range(1, len(vals)) if vals[i] != vals[i - 1]]
df_sorted["_xlab"] = df_sorted["CloneKey"].astype(str) + " — " + df_sorted["ID"].astype(str)

for r in RADII:
    col = f"T{int(r)}"
    fig, ax = plt.subplots(
        figsize=(max(14, 0.24 * len(df_sorted)), 5.8),
        dpi=DPI,
    )
    bar_cols = [clone_colors[k] for k in df_sorted["CloneKey"]]
    ax.bar(
        np.arange(len(df_sorted)),
        df_sorted[col].values,
        color=bar_cols,
        edgecolor="#222222",
        linewidth=0.35,
        width=0.86,
        zorder=2,
    )
    for pos in cuts:
        ax.axvline(pos - 0.5, color="#d8d8d8", lw=0.8, alpha=0.9, zorder=0)
    ax.set_title(
        f"Total 0–{int(r)} µm — per cell (colored by CloneKey)",
        fontsize=11,
    )
    ax.set_xticks(np.arange(len(df_sorted)))
    ax.set_xticklabels(df_sorted["_xlab"].values, rotation=90, fontsize=6)
    ax.set_ylabel("Total intersections (0–r)")
    tidy(ax)

    handles = [
        Patch(facecolor=clone_colors[k], edgecolor="none",
              label=f"{k} (n={df_sorted['CloneKey'].eq(k).sum()})")
        for k in all_keys
    ]
    ax.legend(
        handles=handles,
        ncol=1,
        fontsize=6,
        frameon=True,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        title="CloneKeys",
    )
    plt.tight_layout(rect=[0, 0, 0.8, 1])
    fname = f"totals_0_{int(r)}_per_cell.svg"
    fig.savefig(os.path.join(OUT_DIR, fname),
                dpi=DPI, bbox_inches="tight")
    plt.show()

# ==========================================================
# D2) Max intersections — per cell (bars colored by CloneKey)
#     X labels = CloneKey — ID
# ==========================================================
fig, ax = plt.subplots(
    figsize=(max(14, 0.24 * len(df_sorted)), 5.6),
    dpi=DPI,
)
bar_cols = [clone_colors[k] for k in df_sorted["CloneKey"]]
ax.bar(
    np.arange(len(df_sorted)),
    df_sorted["MAX"].values,
    color=bar_cols,
    edgecolor="#222222",
    linewidth=0.35,
    width=0.86,
    zorder=2,
)
for pos in cuts:
    ax.axvline(pos - 0.5, color="#d8d8d8", lw=0.8, alpha=0.9, zorder=0)
ax.set_title(
    "Max Sholl intersections — per cell (colored by CloneKey)",
    fontsize=11,
)
ax.set_xticks(np.arange(len(df_sorted)))
ax.set_xticklabels(df_sorted["_xlab"].values, rotation=90, fontsize=6)
ax.set_ylabel("Max intersections")
tidy(ax)

handles = [
    Patch(facecolor=clone_colors[k], edgecolor="none",
          label=f"{k} (n={df_sorted['CloneKey'].eq(k).sum()})")
    for k in all_keys
]
ax.legend(
    handles=handles,
    ncol=1,
    fontsize=6,
    frameon=True,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    title="CloneKeys",
)
plt.tight_layout(rect=[0, 0, 0.8, 1])
fig.savefig(os.path.join(OUT_DIR, "max_per_cell.svg"),
            dpi=DPI, bbox_inches="tight")
plt.show()

# ==========================================================
# F) Radius at max intersections — per cell (bars colored by CloneKey)
#     X labels = CloneKey — ID
# ==========================================================
fig, ax = plt.subplots(
    figsize=(max(14, 0.24 * len(df_sorted)), 5.6),
    dpi=DPI,
)
bar_cols = [clone_colors[k] for k in df_sorted["CloneKey"]]

ax.bar(
    np.arange(len(df_sorted)),
    df_sorted["R_AT_MAX"].values.astype(float),
    color=bar_cols,
    edgecolor="#222222",
    linewidth=0.35,
    width=0.86,
    zorder=2,
)

for pos in cuts:
    ax.axvline(pos - 0.5, color="#d8d8d8", lw=0.8, alpha=0.9, zorder=0)

for xref in RADII:
    ax.axhline(y=float(xref), color="#e0e0e0", lw=0.6, zorder=0, alpha=0.7)

vals = df_sorted["R_AT_MAX"].values.astype(float)
ymin = float(np.nanmin(vals))
ymax = float(np.nanmax(vals))
pad = max(2.0, 0.05 * (ymax - ymin))
ax.set_ylim(ymin - pad, ymax + pad)

ax.set_title(
    "Radius at max Sholl intersections — per cell (colored by CloneKey)",
    fontsize=11,
)
ax.set_xticks(np.arange(len(df_sorted)))
ax.set_xticklabels(df_sorted["_xlab"].values, rotation=90, fontsize=6)
ax.set_ylabel("Radius of max intersections (µm)")
ax.yaxis.grid(True, linestyle=":", alpha=0.5)
tidy(ax)

handles = [
    Patch(facecolor=clone_colors[k], edgecolor="none",
          label=f"{k} (n={df_sorted['CloneKey'].eq(k).sum()})")
    for k in all_keys
]
ax.legend(
    handles=handles,
    ncol=1,
    fontsize=6,
    frameon=True,
    loc="upper left",
    bbox_to_anchor=(1.02, 1.0),
    title="CloneKeys",
)

plt.tight_layout(rect=[0, 0, 0.8, 1])
fig.savefig(os.path.join(OUT_DIR, "radius_at_max_per_cell.svg"),
            dpi=DPI, bbox_inches="tight")
plt.show()

# ==========================================================
# G) Export CSV report: df + full Sholl curve
# ==========================================================
curve_rows = []
for entry in cell_curves:
    row = {
        "File": entry["File"],
        "ID": entry["ID"],
        "CloneKey": entry["CloneKey"],
    }
    curve = entry["curve"]
    for r, val in zip(common_r, curve):
        col_name = f"Sholl_{int(r)}"
        row[col_name] = int(val)
    curve_rows.append(row)

df_curves = pd.DataFrame(curve_rows)
df_report = pd.merge(
    df,
    df_curves,
    on=["File", "ID", "CloneKey"],
    how="left",
)

csv_path = os.path.join(
    OUT_DIR,
    "sholl_per_cell_report_with_curves.csv"
)
df_report.to_csv(csv_path, index=False)

print(
    f"\n✅ Done. Cells: {len(df)} | CloneKeys: {df['CloneKey'].nunique()} "
    f"| Example clones: {top_example_clones}"
)
print(f"CSV report saved to: {csv_path}")
