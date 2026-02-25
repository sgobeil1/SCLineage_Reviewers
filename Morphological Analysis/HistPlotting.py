

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
sys.path.append(r"Z:\People\Sophie\9 Mouse SC lineage paper\20251124 Data Analysis and Plots\Code Versions\ClonalAnalysis") # Folder containing ColPals
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
        sns.countplot(x=array, ax=ax, color=colour, order = order)
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
    # Make a copy to avoid modifying your original df
    df = df.copy()

    # Extract the first 8 characters from the Clone column as animal ID
    df['Animal_ID'] = df['Clone'].str[:8]

    # Group by TM and compute:
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
    # Gets Animal from Clone ID: e.g. S107_a7_c3 returns S107_a7
    parts = str(clone_str).split("_")
    if len(parts) >= 3:
        return "_".join(parts[:2])
    
    # Read CSV using pandas
def read_csv(file_path, skip_rows):
    if file_path == 'None':
        points = None
    else: 
        points = pd.read_csv(file_path, skiprows=skip_rows) # skip initial rows
        points.dropna(axis = 0, how = 'all', inplace = True) # drop empty rows
    
    return points

#########################################################################################################################################################
# Import the rg_clones output csv table from ClonalAnalysis
#########################################################################################################################################################
input_file_path_cells = r'Z:\People\Francisco\Morphological Analysis\Data Tables/Morphotype_Integrated_Full_FINAL_with_sholl.csv'
input_file_path_Hi = r'Z:\People\Francisco\Morphological Analysis\Histograms Fig 3/clone_heterogeneity_table.csv'
clones_df = read_csv(input_file_path_cells, skip_rows = 0)
Hi_df = read_csv(input_file_path_Hi, skip_rows = 0)

# Find current date and time, which will be used in file naming
now = datetime.datetime.now()
curr_date = now.strftime("%Y%m%d_%H%M%S")
    
# Create folder for the output
chosen_name = '_' # Choose the folder name
output_folder_name = curr_date + chosen_name
overall_output_file_path = r'Z:\People\Francisco\Morphological Analysis\Histograms Fig 3\Plots/' + output_folder_name
os.makedirs(overall_output_file_path) # Create folder

# Keep text as text (not paths) in SVG/PDF/PS
mpl.rcParams['svg.fonttype'] = 'none'  # <-- critical for Illustrator
mpl.rcParams['pdf.fonttype'] = 42      # TrueType; helps if you also export PDF
mpl.rcParams['ps.fonttype']  = 42

# (Optional) pick a common font that Illustrator will have
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']  # or 'Helvetica', etc.



# Save plotted data in an excel spreadsheet where first sheet contains the data and the second sheet contains some statistics. 
# df is the df inputted for plotting
# columns is the list of columns of df you want saved in the spreadsheet
# dep_var is the dependent variable (y column) on which stats should be performed
def saveDataTableGroupMeans(df, columns, dep_var, ind_var, output_path, order):
    df = df[columns]
    # Derive Animal_ID from first 8 characters of Clone
    df['Animal'] = df['Clone'].apply(extract_animal_id)
    # Move Animal_ID column to the front
    cols = ['Animal'] + [c for c in df.columns if c != 'Animal']
    df = df[cols]
    # Sort column in defined order
    df[ind_var] = pd.Categorical(df[ind_var], categories=order, ordered=True)
    df = df.sort_values(ind_var).reset_index(drop=True)
    # Per-TM summary stats
    stats_df = (
        df
        .groupby(ind_var, as_index=False)
        .agg(
            Mean=(dep_var, 'mean'),
            SEM=(dep_var, 'sem'),
            N_clones=('Clone', 'size'),
            N_animals=('Animal', lambda x: x.nunique())
        )
    )
    # Save to Excel
    excel_path = output_path + ".xlsx"  

    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        df.to_excel(writer, sheet_name="Data", index=False)
        stats_df.to_excel(writer, sheet_name="Stats", index=False)


######################################################################################################################################################################
### Histogram of # primary neurites


neurites = clones_df['n_primary'].to_numpy()

output_path = overall_output_file_path + '/Primary_Neurites_Hist_' + curr_date 
snsHistPlot(array = neurites, 
        x_axis_title = 'Number of primary neurites',
        y_axis_title = 'Number of cells',
        output_path = output_path,
        colour = ColPals.lightgrey,
        binwidth = 1,
        binrange = [0,12],
        fig_size = (4,5),
        order = 2)


######################################################################################################################################################################
### Histogram of total Sholl intersections by 100 um radius


t100 = clones_df['T100'].to_numpy()

print('\nt100\n', t100)

output_path = overall_output_file_path + '/T100_Hist_' + curr_date 
snsHistPlot(array = t100, 
        x_axis_title = 'Total neurite intercepts (0-100 \u03bcm)',
        y_axis_title = 'Number of cells',
        output_path = output_path,
        colour = ColPals.lightgrey,
        binwidth = 50,
        binrange = [0,1000],
        fig_size = (5,5),
        order = 200)


######################################################################################################################################################################
### Histogram of geometries


geometry = clones_df['Geometry'].to_numpy()

print('\ngeometry\n', geometry)

output_path = overall_output_file_path + '/Geometry_Hist_' + curr_date 
snsHistPlot(array = geometry, 
        x_axis_title = 'Cell geometry',
        y_axis_title = 'Number of cells',
        output_path = output_path,
        colour = ColPals.lightgrey,
        fig_size = (3,5),
        order = ['Unidirectional', 'Bipolar', 'Bidirectional', 'Triangular', 'Stellate']
        )

######################################################################################################################################################################
### Histogram of morphotypes


morphotypes = clones_df['CellType_With_Complexity'].to_numpy()
'Unidirectional_Low', 'Unidirectional_High', 'Bidirectional_Low', 'Bidirectional_High', 'Bipolar_Low', 'Bipolar_High,', 'Triangular_Low', 'Triangular_High', 'Stellate_Low', 'Stellate_High'

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


print('\ngeometry\n', geometry)

output_path = overall_output_file_path + '/Morphotypes_Hist_' + curr_date 
snsHistPlot(array = morphotypes, 
        x_axis_title = 'Morphotypes',
        y_axis_title = 'Number of cells',
        output_path = output_path,
        colour = ColPals.lightgrey,
        fig_size = (4.5,5),
        order = ['Simple Unidirectional', 'Complex Unidirectional', 'Simple Bipolar', 'Complex Bipolar', 'Simple Bidirectional', 'Complex Bidirectional', 'Simple Triangular', 'Complex Triangular', 'Simple Stellate', 'Complex Stellate']
        )


######################################################################################################################################################################
### Histogram of Heterogeneity index (Hi)

Hi = Hi_df['HI_clonal'].to_numpy()

output_path = overall_output_file_path + '/Hi_Hist_' + curr_date 
snsHistPlot(array = Hi, 
        x_axis_title = 'Hi',
        y_axis_title = 'Number of clones',
        output_path = output_path,
        colour = ColPals.lightgrey,
        fig_size = (4,5),
        binwidth = 0.1,
        binrange = [0,1],
        order = 0.2)
