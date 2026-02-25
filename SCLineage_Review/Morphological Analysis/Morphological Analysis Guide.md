# SCLineage_Review Morphological Analysis

This file contains information about the clone morphological analysis pipeline used in **Gobeil et al. (2026)** and how to run the demo.

## Overview
This folder (`Morphological Analysis/`) contains scripts used in **Gobeil et al. (2026)** to quantify morphometric features of reconstructed neurons and classify them into morphotypes. It also provides scripts to assess the morphological heterogeneity of clonally related neurons in MADM clones.

### Contents (scripts)
- **`MorphotypeClassifier.py`**  
   Main demo-ready pipeline: loads `Cell_Info.csv`, scans SWC files, performs angle analysis + final morphotype classification, and computes complexity from SWC (T100) internally.  
   **Outputs (examples):**
   - `cell_analysis_detailed_with_complexity_T100_from_SWC.csv`
   - `cell_type_distribution_with_complexity_T100_from_SWC.svg`
   - `PER_CLASS_OVERVIEWS/` (example cells per class)

- **`ShollAnalysis.py`** *(optional, Sholl curves)*  
   Computes per-cell and per-clone Sholl curves + totals, MAX intersections, and summary plots.  
   **Outputs (examples):**
   - `sholl_per_cell_report_with_curves.csv`
   - SVG plots in `OUT_DIR/`

- **`HeterogeneityIndex.py`**  
   Computes a clone heterogeneity index from a morphotype table (Excel).  
   **Outputs:**
   - `clone_heterogeneity_index_*.svg/.jpg`
   - `clone_heterogeneity_index_reconstructed_rows_with_HI.csv`

- **`bootstrap comparison.py`** *(optional)*  
   Permutation / bootstrap analysis of clone morphotype richness (observed vs expected), comparing real data to homogeneous/heterogeneous simulated datasets.

- **`HistPlotting.py`** *(paper-figure helper)*  
   Additional plotting helper; expects inputs produced by the scripts above (and/or manuscript tables). Not required for the demo.

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
- `openpyxl` (required by pandas to read `.xlsx`)

---

## 2. Installation guide
From a terminal:

```bash
pip install numpy pandas scipy matplotlib seaborn pypalettes openpyxl
```

**Typical install time:** ~5–10 minutes on a normal desktop (depends on internet speed).

---

## 3. Demo

The demo inputs live in `Morphological Analysis/Demo/`:

- `Cell_Info.csv` — cell metadata exported from Imaris (one row per dendrite; script deduplicates to one row per cell)
- `S107/`, `S140/`, `S141/` — example SWC reconstructions in the expected folder hierarchy
- `Clone_morphotypes.xlsx` — example morphotype table (used by `HeterogeneityIndex.py`)
- `Data Tables/` — simulated datasets used by `bootstrap comparison.py`

Folder hierarchy for SWCs (as in the demo):
```
Demo/<Litter S###>/<Animal S###_a#>/<Segment s1..s6>/<Clone clone #>/<cell_id>.swc
```

### 3.1 Download the code + demo files
GitHub UI: **Code → Download ZIP**, then unzip

### 3.2 Run MorphotypeClassifier

1) Open `MorphotypeClassifier.py` and update the `# User input` paths to:

```python
CSV_PATH = r"Demo/Cell_Info.csv"
SWC_ROOT = r"Demo"
OUT_DIR  = r"Demo/Output_MorphotypeClassifier"
```

2) Run from this folder:
```bash
cd "Morphological Analysis"
python MorphotypeClassifier.py
```

**Expected output**
In `Demo/Output_MorphotypeClassifier/`:
- `cell_analysis_detailed_with_complexity_T100_from_SWC.csv`
- `cell_type_distribution_with_complexity_T100_from_SWC.svg`
- `PER_CLASS_OVERVIEWS/`

**Expected runtime**
Typically **< 5 minutes** for the included demo SWCs.

---

### 3.3 (Optional) Run ShollAnalysis

1) Open `ShollAnalysis_20251217.py` and update:
```python
SWC_ROOT = r"Demo"
OUT_DIR  = r"Demo/Output_Sholl"
```

2) Run:
```bash
cd "Morphological Analysis"
python ShollAnalysis_20251217.py
```

**Expected output**
- `Demo/Output_Sholl/sholl_per_cell_report_with_curves.csv`
- SVG plots in `Demo/Output_Sholl/`

**Expected runtime**
Typically **minutes** for the demo.

---

### 3.4 Run HeterogeneityIndex (clone heterogeneity)

1) Open `HeterogeneityIndex.py` and update:
```python
EXCEL_PATH = r"Demo/Clone_morphotypes.xlsx"
OUT_DIR    = r"Demo/Output_HeterogeneityIndex"
# Demo file uses 'CellType' (not 'CellType_With_Complexity'):
MORPHO_COL  = "CellType"
```

2) Run:
```bash
cd "Morphological Analysis"
python HeterogeneityIndex.py
```

**Expected output**
In `Demo/Output_HeterogeneityIndex/`:
- `clone_heterogeneity_index_*.svg` and `.jpg`
- `clone_heterogeneity_index_reconstructed_rows_with_HI.csv`

**Expected runtime**
Typically **< 1 minute** for the demo.

---

### 3.5 (Optional) Run bootstrap / permutation richness analysis

`bootstrap comparison.py` compares a “real” dataset to two simulated datasets in `Demo/Data Tables/`.

1) Open `bootstrap comparison.py` and update:
```python
REAL_EXCEL_PATH = r"Demo/Clone_morphotypes.xlsx"
HOMO_EXCEL_PATH = r"Demo/Data Tables/Homogeneous Simulated Dataset.xlsx"
HET_EXCEL_PATH  = r"Demo/Data Tables/Heterogeneous Simulated Dataset.xlsx"
OUTPUT_DIR      = r"Demo/Output_Bootstrap"
```

2) Run:
```bash
cd "Morphological Analysis"
python "bootstrap comparison.py"
```

---

## 4) Instructions for use (your own data)

### 4.1 Prepare inputs

1) Export neuron reconstructions from Imaris (Filament Tracer) as **SWC** files.  
2) Organize SWCs in the same hierarchy used by the demo:
```
<SWC_ROOT>/<Litter>/<Animal>/<Segment>/<Clone>/<cell_id>.swc
```
- Litter folder should look like `S107`
- Animal folder should look like `S107_a7`
- Clone folder name must contain the word `clone` (e.g. `clone 2`)

3) Export a cell metadata table (CSV) with columns like:
- `Litter`, `Animal`, `Segment`, `Clone`, `ID`
The script is tolerant to capitalization and will try to infer columns by name.

### 4.2 Run
- Set `CSV_PATH`, `SWC_ROOT`, `OUT_DIR` in `MorphotypeClassifier.py`
- Run:
```bash
python MorphotypeClassifier.py
```

---

## Code availability
Code repository: `https://github.com/sgobeil1/SCLineage_Review`

---

## Author
Francisco Da Silveira Neto

Institute of Science and Technology Austria
