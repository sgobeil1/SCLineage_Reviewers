# SCLineage_Review Clonal Analysis - Spatial Analysis

This file contains information about the clone spatial analysis pipeline used in **Gobeil et al. (2026)** and how to run the demo.

## Overview
This folder (`Clonal Analysis/Spatial Analysis/`) contains scripts for clone **spatial analyses**, focused on comparing the spatial distribution of clone neurons to **reference AB (cardinal class) distributions** in transverse spinal cord coordinates.

### Contents (scripts)
- **`CompareDistributions.py`**  
  Computes a 2D KDE for each AB reference distribution and scores each **query cell** against each AB distribution by the **smallest HDR percentile** it falls into.  
  **Output:** `query_cells_percentile_membership_*.csv` (+ overlay plots).

- **`CloneDistributionHeatmap.py`**  
  Takes the **per-cell** membership table from `CompareDistributions.py`, collapses to **per-clone × AB** (minimum percentile per clone), and generates heatmaps + mixing metrics.  
  **Output:** clone×AB matrix CSVs, mixing stats CSV, and heatmaps/plots.

- **`DistributionOverlapHeatmap.py`**  
  Quantifies overlap (e.g., Jaccard/Dice) between **any two sets** of distributions (often AB vs AB) based on a target HDR percentile and plots heatmaps.

- **`SCNormalization_BLT.py`**
  Normalizes raw cell coordinates into a standard reference spinal cord coordinate system using CC + bottom/left/top reference points.

- **`QuadColocSpotsWithTransformCoordinates.m`**
  MATLAB (R2022b) script that should be run in Imaris. Originally written by Mariano Gabitto and modified by Sophie Gobeil. Exports Spots objects (cell coordinates) from Imaris. Used to export the coordinates of cells expressing AB markers for the AB distribution spatial map. Used in conjunction with SpotNormScalingRotation.py. Not used in the demo.

- **`SpotNormScalingRotation.py `**
  Normalizes, rotates, and scales cell 2D coordinates to a spinal cord hemisection of fixed size. Used in conjunction with QuadColocSpotsWithTransformCoordinates.m. Not used in the demo.

---

## 1. System requirements

### Operating system
Tested on Windows 10.

### Software
Python **3.10+** (written on 3.12.8)

### Python dependencies
- `numpy`
- `pandas`
- `scipy`
- `matplotlib`
- `seaborn`
- `pypalettes`

> Note: Some scripts import local helper modules from the `Clonal Analysis/` folder
> (e.g. `ClonalAnalysisFunctions.py`, `ClonalAnalysisPlots.py`, `ColPals.py`).

---

## 2. Installation guide
From a terminal:

```bash
pip install numpy pandas scipy matplotlib seaborn pypalettes
```

**Typical install time:** ~5–10 minutes on a normal desktop (depends on internet speed).

---

## 3. Demo

This demo should be run after the Clonal Analysis demo.

The demo inputs live in:

- `Demo/AB Files/` — reference AB distributions (CSV per AB class)
- `Demo/Overall_rg_cells.csv` — example clone-annotated cell table (optional; for generating query clone CSVs)
- `Demo/Clone_Ref.csv` — example reference-point table (used by some workflows)
- `cord.mat` — hemicord outline used for overlay plots (in this folder)

### 3.1 Download the code + demo files
**Code → Download ZIP**, then unzip  

### 3.2 Run

#### Normalize clone cell coordinates
Open `SCNormalization_BLT.py` and edit the `### User input` section paths:
```python
input_cells_file_path = r'Overall_rg_cells.csv'
input_ref_folder_path = r'Ref Files'
output_folder_name = '' # Desired output folder name
overall_output_file_path = r'Output/' # Desired output path
```

Additional criteria for clone/cell exclusion can also be specified here.

**Expected outputs**
Individual csv files for each clone containing normalized coordinates.

#### Run per-cell KDE scoring (CompareDistributions)
Open `CompareDistributions.py` and edit the `### User input` section paths:

```python
CORD_MAT_PATH = r"cord.mat"
REF_DIR = r"Demo/AB Files"
QUERY_DIR = r"Demo/QueryClones"
OUTPUT_DIR = r"Demo/Output_CompareDistributions_"
```

**Expected outputs**
A timestamped folder is created under `OUTPUT_DIR + curr_date` containing:
- `query_cells_percentile_membership_*.csv`
- overlay plots per query clone (`clone_overlay_*.png/.svg`)

**Expected runtime (demo)**
Typically **minutes** on a normal desktop (depends on number of query clones and AB references).

---

#### Summarize per-clone overlap + plots (CloneDistributionHeatmap)

Open `CloneDistributionHeatmap.py` and edit the `### User input` section paths:

```python
INPUT_PATH = r"Demo/Output_CompareDistributions_<timestamp>/query_cells_percentile_membership_<timestamp>.csv"
OUTPUT_DIR = r"Demo/Output_CloneDistributionHeatmap"
```

**Expected outputs**
In `OUTPUT_DIR/Plotting_<timestamp>/` (and/or `OUTPUT_DIR` depending on settings), files include:
- `query_cells_percentile_membership_matrix.csv` (AB × clone minima)
- `query_cells_percentile_membership_clone_mixing.csv`
- heatmaps (`*.png`, `*.svg`)
- optional per-cell predicted-identity table (`per_cell_predicted_identity.csv`)

**Expected runtime (demo)**
Typically **minutes** on a normal desktop.

---

#### Optional: AB ↔ AB overlap heatmaps (DistributionOverlapHeatmap)

This script is useful to quantify overlap between AB distributions themselves.

Open `DistributionOverlapHeatmap.py` and edit the `### User input` section paths:

```python
REF_DIR = r"Demo/AB Files"
QUERY_DIR = r"Demo/AB Files"  # same as REF_DIR for AB vs AB
OUTPUT_DIR = r"Demo/Output_ABvsAB_Overlap_<timestamp>"
```

**Expected outputs**
- long-form overlap table (CSV)
- wide overlap matrix (CSV)
- heatmaps for selected metrics (`.png` + `.svg`)

---

## 4. Instructions for use (your own data)

### 4.1 Prepare inputs

To run the spatial pipeline you need two sets of 2D transverse coordinates in the **same coordinate system** (same scaling, same origin, same hemicord orientation):

1) **AB reference distributions (“AB Files”)**  
   One CSV per AB class (e.g., `dI1`, `dI5`, `V2`, `MN`), containing at minimum X/Y coordinates.  
   Recommended column names: `Position X`, `Position Y` (or `X coordinate (scaled to Reference SC)` and `Y coordinate (scaled to Reference SC)`).

   In the manuscript workflow, AB coordinates are exported from Imaris using:
   - `QuadColocSpotsWithTransformCoordinates.m` (MATLAB; run in Imaris)
   - then normalized/scaled with `SpotNormScalingRotation.py`

2) **Query clone distributions (“QueryClones”)**  
   One CSV per clone, containing at minimum X/Y coordinates in the same coordinate frame as the AB references.  
   If you already have a clone cell table in this frame, you can generate per-clone CSVs by splitting by clone:

```python
import os
import pandas as pd

df = pd.read_csv("Overall_rg_cells.csv")
out_dir = "QueryClones"
os.makedirs(out_dir, exist_ok=True)

for clone, g in df.dropna(subset=["Clone", "Position X", "Position Y"]).groupby("Clone"):
    g[["Position X", "Position Y"]].to_csv(os.path.join(out_dir, f"{clone}.csv"), index=False)
```

   If you are starting from 3D clonal outputs (e.g., `Overall_rg_cells_*.csv` from `ClonalAnalysis.py`), you must first generate the appropriate 2D transverse coordinates using your normalization workflow (e.g., `SCNormalization_BLT.py` if you have the required reference files).

3) **`cord.mat` (optional but recommended)**  
   Used for overlay plots of distributions on a hemicord outline.

### 4.2 Run the spatial analysis scripts

1) **Per-cell KDE scoring** (`CompareDistributions.py`)  
   Set in the `### User input` block:
   - `CORD_MAT_PATH` (optional)
   - `REF_DIR` (folder containing AB reference CSVs)
   - `QUERY_DIR` (folder containing per-clone CSVs)
   - `OUTPUT_DIR` (output base)

   Then run:
```bash
python CompareDistributions.py
```

2) **Per-clone summarization + heatmaps** (`CloneDistributionHeatmap.py`)  
   Set:
   - `INPUT_PATH` to the `query_cells_percentile_membership_*.csv` produced above
   - `OUTPUT_DIR` to your output base

   Then run:
```bash
python CloneDistributionHeatmap.py
```

3) *(Optional)* **Distribution overlap heatmaps** (`DistributionOverlapHeatmap.py`)  
   Use this to quantify AB↔AB overlap (or any distribution set ↔ distribution set) for a target HDR percentile.

### 4.3 Notes on coordinate conventions

- If you want to analyze left/right hemispheres separately, `CompareDistributions.py` can use a `Hemisphere` column (when `HEMISPHERES='LR'`).  
- Otherwise, you can choose to mirror or randomly assign hemispheres using the options in the script’s user‑input block.

---

## Code availability
Code repository: `https://github.com/sgobeil1/SCLineage_Review`

---
## Author
Sophie A. Gobeil

Institute of Science and Technology Austria
