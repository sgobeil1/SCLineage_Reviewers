"""
================================================================================
HISTOGRAM PANEL GENERATOR — MORPHOLOGY + CLONE HETEROGENEITY
PRIMARY NEURITES • T100 SHOLL TOTALS • GEOMETRY • MORPHOTYPES • Hi DISTRIBUTION
================================================================================
Purpose:
Generate standardized histograms from the MorphotypeClassifier
cell-level output and a clone-level heterogeneity table.

Inputs:
- Cell table (CSV):
    • n_primary   (number of primary neurites)
    • T100        (total Sholl intersections 0–100 µm)
    • Geometry    (geometry class)
    • Morphotype  (geometry + complexity)

- Clone table (CSV):
    • HeterogeneityIndex  OR  observed_H  (auto-detected)

What this script produces (SVG, timestamped folder):
- Histogram: Number of primary neurites (per cell)
- Histogram: T100 (per cell)
- Histogram: Geometry categories (per cell)
- Histogram: Morphotypes (per cell; relabeled to Simple/Complex)
- Histogram: Heterogeneity index Hi (per clone; 0–1)

Output:
- Creates a timestamped directory and saves each plot as .svg with embedded text
  (fonts preserved for Illustrator/Inkscape workflows).

Notes:
- Hi column is automatically detected via `load_hi_values()` to support both:
    • clone_heterogeneity_index_per_clone_summary.csv  (HeterogeneityIndex)
    • clone_richness_H_summary_real_Morphotype.csv     (observed_H)
================================================================================
"""
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import distance_matrix
from scipy import stats
import numpy as np
import copy
import matplotlib.ticker as ticker
import sys
sys.path.append(r"Z:\People\Sophie\9 Mouse SC lineage paper\20251124 Data Analysis and Plots\Code Versions\ClonalAnalysis")
import ColPals
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import datetime
from statistics import mean
import random
from pypalettes import load_cmap
from itertools import combinations
import seaborn as sns

#########################################################################################################################################################
### Functions
#########################################################################################################################################################
def snsHistPlot(array, x_axis_title, y_axis_title, output_path, colour,
                binwidth=None, binrange=None, order=None, max_y=None,
                fig_size=(6.5,5)):

    fig, ax = plt.subplots(figsize=fig_size)

    array = pd.Series(array)

    if pd.api.types.is_numeric_dtype(array):
        sns.histplot(array, color=colour, ax=ax,
                     binwidth=binwidth, binrange=binrange)

        if order is not None:
            ax.xaxis.set_major_locator(MultipleLocator(order))
    else:
        # countplot handles ordering but we must avoid set_ticklabels warning
        sns.countplot(x=array, ax=ax, color=colour, order=order)

        # FIX: set ticks explicitly BEFORE setting ticklabels
        ticks = ax.get_xticks()
        ax.set_xticks(ticks)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')

    sns.despine()
    ax.tick_params(axis='x', labelsize=16)
    ax.tick_params(axis='y', labelsize=16)
    ax.set_xlabel(x_axis_title, fontsize=18, weight='bold')
    ax.set_ylabel(y_axis_title, fontsize=18, weight='bold')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.tick_params(width=1.5)

    plt.yticks(weight='bold')
    plt.xticks(weight='bold')

    if max_y is not None:
        plt.ylim(0, max_y)

    plt.tight_layout()
    fig.savefig(output_path + '.svg', format='svg', dpi=300, bbox_inches='tight')


def get_n(df):
    df = df.copy()
    df['Animal_ID'] = df['Clone'].str[:8]
    summary = (
        df.groupby('TM')
        .agg(
            Number_of_clones=('Clone', 'count'),
            Number_of_animals=('Animal_ID', lambda x: x.nunique())
        )
        .reset_index()
    )
    return summary


def extract_animal_id(clone_str: str) -> str:
    parts = str(clone_str).split("_")
    if len(parts) >= 3:
        return "_".join(parts[:2])


def read_csv(file_path, skip_rows):
    if file_path == 'None':
        points = None
    else:
        points = pd.read_csv(file_path, skiprows=skip_rows)
        points.dropna(axis=0, how='all', inplace=True)
    return points


def load_hi_values(Hi_df: pd.DataFrame) -> np.ndarray:
    """
    UPDATED: Supports both formats:
      A) clone_heterogeneity_index_per_clone_summary.csv
         -> uses column 'HeterogeneityIndex'
      B) clone_richness_H_summary_real_Morphotype.csv
         -> uses column 'observed_H'
    Also accepts a few common fallbacks.
    """
    candidate_cols = [
        "observed_H",
        "HeterogeneityIndex",
        "HeterogeneityIndex_in_clone",
        "HI_clonal",  # legacy, if ever present
        "HI",
        "H"
    ]
    for c in candidate_cols:
        if c in Hi_df.columns:
            return Hi_df[c].to_numpy()

    raise KeyError(
        "Could not find a heterogeneity column in Hi_df. "
        "Expected one of: " + ", ".join(candidate_cols) +
        f"\nAvailable columns: {Hi_df.columns.tolist()}"
    )

#########################################################################################################################################################
# Inputs
#########################################################################################################################################################
input_file_path_cells = r"Z:\People\Francisco\Code_testing_folder\Output\cell_analysis_detailed_with_complexity.csv"

# HI input can be EITHER:
# 1) clone_heterogeneity_index_per_clone_summary.csv  (HeterogeneityIndex)
# OR
# 2) clone_richness_H_summary_real_Morphotype.csv     (observed_H)
input_file_path_Hi = r"Z:\People\Francisco\Code_testing_folder\Output\Clone_Richness_Permutation\clone_richness_permutation_results_real_Morphotype.csv"
# input_file_path_Hi = r"Z:\People\Francisco\Code_testing_folder\Output\Clone_Richness_Permutation\clone_richness_H_summary_real_Morphotype.csv"

clones_df = read_csv(input_file_path_cells, skip_rows=0)
Hi_df = read_csv(input_file_path_Hi, skip_rows=0)

# Timestamp output folder
now = datetime.datetime.now()
curr_date = now.strftime("%Y%m%d_%H%M%S")

chosen_name = '_'
output_folder_name = curr_date + chosen_name
overall_output_file_path = r'Z:\People\Francisco\Morphological Analysis\Histograms Fig 3\Plots/' + output_folder_name
os.makedirs(overall_output_file_path, exist_ok=True)

# Keep text as text (not paths) in SVG/PDF/PS
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']

######################################################################################################################################################################
### Histogram of # primary neurites
neurites = clones_df['n_primary'].to_numpy()

output_path = overall_output_file_path + '/Primary_Neurites_Hist_' + curr_date
snsHistPlot(
    array=neurites,
    x_axis_title='Number of primary neurites',
    y_axis_title='Number of cells',
    output_path=output_path,
    colour=ColPals.lightgrey,
    binwidth=1,
    binrange=[0, 12],
    fig_size=(4, 5),
    order=2
)

######################################################################################################################################################################
### Histogram of total Sholl intersections by 100 um radius
t100 = clones_df['T100'].to_numpy()
print('\nt100\n', t100)

output_path = overall_output_file_path + '/T100_Hist_' + curr_date
snsHistPlot(
    array=t100,
    x_axis_title='Total neurite intercepts (0-100 \u03bcm)',
    y_axis_title='Number of cells',
    output_path=output_path,
    colour=ColPals.lightgrey,
    binwidth=50,
    binrange=[0, 1000],
    fig_size=(5, 5),
    order=200
)

######################################################################################################################################################################
### Histogram of geometries
geometry = clones_df['Geometry'].to_numpy()
print('\ngeometry\n', geometry)

output_path = overall_output_file_path + '/Geometry_Hist_' + curr_date
snsHistPlot(
    array=geometry,
    x_axis_title='Cell geometry',
    y_axis_title='Number of cells',
    output_path=output_path,
    colour=ColPals.lightgrey,
    fig_size=(3, 5),
    order=['Unidirectional', 'Bipolar', 'Bidirectional', 'Triangular', 'Stellate']
)

######################################################################################################################################################################
### Histogram of morphotypes
morphotypes = clones_df['Morphotype'].to_numpy()

morphotypes = pd.Series(morphotypes).replace({
    "Unidirectional_Low": "Simple Unidirectional",
    "Unidirectional_High": "Complex Unidirectional",
    "Bipolar_Low": "Simple Bipolar",
    "Bipolar_High": "Complex Bipolar",
    "Bidirectional_Low": "Simple Bidirectional",
    "Bidirectional_High": "Complex Bidirectional",
    "Triangular_Low": "Simple Triangular",
    "Triangular_High": "Complex Triangular",
    "Stellate_Low": "Simple Stellate",
    "Stellate_High": "Complex Stellate",
}).to_numpy()

output_path = overall_output_file_path + '/Morphotypes_Hist_' + curr_date
snsHistPlot(
    array=morphotypes,
    x_axis_title='Morphotypes',
    y_axis_title='Number of cells',
    output_path=output_path,
    colour=ColPals.lightgrey,
    fig_size=(4.5, 5),
    order=[
        'Simple Unidirectional', 'Complex Unidirectional',
        'Simple Bipolar', 'Complex Bipolar',
        'Simple Bidirectional', 'Complex Bidirectional',
        'Simple Triangular', 'Complex Triangular',
        'Simple Stellate', 'Complex Stellate'
    ]
)

######################################################################################################################################################################
### Histogram of Heterogeneity index (Hi)
Hi = load_hi_values(Hi_df)

output_path = overall_output_file_path + '/Hi_Hist_' + curr_date
snsHistPlot(
    array=Hi,
    x_axis_title='Hi',
    y_axis_title='Number of clones',
    output_path=output_path,
    colour=ColPals.lightgrey,
    fig_size=(4, 5),
    binwidth=0.1,
    binrange=[0, 1],
    order=0.2
)
