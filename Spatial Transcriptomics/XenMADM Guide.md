# Xenium Spatial Analysis

## Overview

This folder contains the Jupyter notebooks used to analyze Xenium spatial transcriptomics data from the developing mouse spinal cord and to relate the expression profiles of lineage-traced MADM cells to transcriptomically defined neuronal populations.

The workflow combines cell-level expression data from 16 Xenium slide regions, performs quality control and reference clustering, identifies and annotates neuronal populations using spinal cord marker genes, and projects manually segmented MADM cells into the resulting references. A separate notebook generates UMAPs, spatial transcript maps, marker-expression plots, and clone-specific figures.

**Data availability:** The notebooks and configuration files document the analysis, but the Xenium and MADM expression datasets are not included in this repository at present. 

## 1. System requirements

### Operating system
Tested on Windows 10.

## Software and Python dependencies

### Software

The analysis was performed using the following software:

| Software | Purpose |
|---|---|
| Python 3.11.8 | Programming language used for data processing and analysis |
| JupyterLab / Jupyter Notebook | Execution of the analysis notebooks |
| 10x Genomics Xenium Onboard Analysis | Generation of spatial transcriptomics data |
| 10x Genomics Xenium Explorer | Visualization and inspection of Xenium data (optional) |
| Custom Z-stack visualizer | Manual segmentation and transcript counting of MADM-labelled cells |

### Python dependencies

The following Python packages are imported by the analysis notebooks:

| Package | Purpose |
|---|---|
| `scanpy` | Single-cell analysis, normalization, PCA, clustering, UMAP, and reference mapping |
| `anndata` | Storage and manipulation of annotated expression matrices |
| `spatialdata` | Reading and manipulating spatial transcriptomics data |
| `spatialdata-io` | Import of Xenium output into SpatialData |
| `spatialdata-plot` | SpatialData plotting functionality |
| `squidpy` | Spatial omics analysis toolkit |
| `numpy` | Numerical calculations and array operations |
| `pandas` | Data manipulation and metadata handling |
| `scipy` | Sparse matrices and scientific computing |
| `scikit-learn` | Nearest-neighbor calculations |
| `h5py` | HDF5 data support |
| `matplotlib` | Data visualization and figure generation |
| `seaborn` | Statistical plotting and heatmaps |
| `adjustText` | Adjustment of text labels to prevent overlap in figures |

Additional dependencies required for clustering and visualization include:

| Package | Purpose |
|---|---|
| `umap-learn` | UMAP dimensionality reduction |
| `python-igraph` | Graph operations used in Leiden clustering |
| `leidenalg` | Leiden community detection |
| `ipykernel` | Python kernel for Jupyter |
| `jupyterlab` | JupyterLab notebook environment |

### Installation

The Python packages can be installed using pip:

```bash
python -m pip install scanpy anndata spatialdata spatialdata-io spatialdata-plot squidpy numpy pandas scipy scikit-learn h5py matplotlib seaborn adjustText umap-learn python-igraph leidenalg jupyterlab ipykernel
```

For reproducibility, installation in a dedicated Python environment is recommended. Exact versions of the Python packages used in the original analysis have not yet been documented.

## 2. Files and workflow organization

| File | Purpose |
| --- | --- |
| `Notebooks/XeniumAnalysis_Preprocessing.ipynb` | Load Xenium SpatialData objects; perform QC, normalization, scaling, and gene filtering; assemble manually segmented MADM cell counts. |
| `Notebooks/XeniumAnalysis_AllCells_Reference.ipynb` | Compute PCA, a nearest-neighbor graph, Leiden clusters, and UMAP for all Xenium cells. |
| `Notebooks/XeniumAnalysis_Neuron_Reference.ipynb` | Score and select neurons, exclude residual non-neuronal clusters, recluster neuronal populations, and assign marker-based cell-type identities. |
| `Notebooks/XeniumAnalysis_MADM.ipynb` | Map MADM cell expression profiles onto the existing all-cell and neuronal references using Scanpy `ingest`. |
| `Notebooks/XeniumAnalysis_Figures.ipynb` | Produce spatial expression plots, neuronal and MADM UMAPs, class-marker scores, heatmaps, and manuscript-oriented figures. |
| `Notebooks/xenium_config.py` | Define marker-group selections, cell-type colors, and gene sets used for neuronal analysis and plotting. |
| `Notebooks/cardinal_class_markers.py` | Define marker sets for scoring cardinal spinal cord neuronal classes in the figure notebook. |
| `Marker_Lists_v4.csv` | Additional marker lists loaded by the reference and figure notebooks; required input file to provide alongside the notebooks. |

The notebooks are intended to be run in the order shown. The figure notebook can also be used independently if the processed reference and mapped MADM objects already exist.

```text
Xenium spatial outputs (.zarr SpatialData objects)
             |                      Manually segmented MADM ROI-count CSVs
             |                                  |
             +---------------+------------------+
                             |
                   1. Preprocessing
                             |
               Xenium_preprocessed.h5ad
                  MADM_preprocessed.h5ad
                             |
                    2. All-cell reference
                             |
               Xenium_allcell_reference.h5ad
                             |
                    3. Neuronal reference
                             |
         Final neuronal reference + MADM neuron subset
                             |
                        4. MADM mapping
                             |
                   Mapped MADM .h5ad files
                             |
                       5. Figures
```

## 3) Required input data

### 3.1 Xenium data

The preprocessing notebook loads one `SpatialData` `.zarr` directory per region, for example `data_60_Region8.zarr`. In the original analysis, 16 regions are specified explicitly. Each SpatialData object must provide a `table` containing cell-level gene expression and metadata, and `cell_boundaries` shapes for extracting cell-centroid coordinates. Spatial transcript plots additionally require the `transcripts` points element.

The notebook includes commented examples of conversion from a 10x Xenium output directory to SpatialData using `spatialdata_io.xenium`; alternatively, already-converted `.zarr` objects can be supplied.

### 3.2 MADM cell transcript counts

MADM-labelled cells were manually delineated using fluorescence images and a custom Z-stack visualization/ROI-counting workflow. For each selected cell, transcripts were counted within a defined XY region and Z interval to limit contributions from overlapping cells.

The preprocessing notebook expects one CSV per manually segmented cell, with `Gene` and `Count` columns, in a local `MADM Cells/` directory. The file names, Z intervals, cell identifiers, and clone assignments are defined explicitly in the notebook and must be adapted for another dataset. The ROI-counting utility itself is upstream of, and is not included in, this set of notebooks.

### 3.3 Marker definitions

`Marker_Lists_v4.csv` supplies marker lists for broad cell groups and spinal cord neuronal classes. `xenium_config.py` selects which of these groups to use, defines neuronal gene panels and plotting colors, and should be placed on the Python import path. `cardinal_class_markers.py` supplies an additional marker dictionary for class-score visualizations. Gene names must match those in the Xenium expression matrix.

## Code availability
Code repository: `https://github.com/sgobeil1/SCLineage_Review`

---

## Author
Sophie A. Gobeil

Institute of Science and Technology Austria

