### Script to make plots from the rg_clones clone summary file. Intended to be more straightforward than plotting
### in the main ClonalAnalysis script

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
import ClonalAnalysisFunctions as caf
import ClonalAnalysisPlots as cap
import ColPals as ColPals
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import datetime
from statistics import mean
import random
from pypalettes import load_cmap
from itertools import combinations
import seaborn as sns

### User input
input_file_path = r'Z:\People\Sophie\9 Mouse SC lineage paper\Initial Submission\Nature Submission\GitHub Repository\Output/_rg_clones_20260224_152812.csv' # Path of _rg_clones csv file
output_folder_name = '_simpleplotting_testing' # Desired output folder name
overall_output_file_path = r'Z:\People\Sophie\9 Mouse SC lineage paper\Initial Submission\Nature Submission\GitHub Repository\Output/' # Desired output folder

### Prep
# Import the rg_clones output csv table from ClonalAnalysis
rg_clones_df = caf.read_csv(input_file_path, skip_rows = 0)

# Find current date and time, which will be used in file naming
now = datetime.datetime.now()
curr_date = now.strftime("%Y%m%d_%H%M%S")
    
# Create folder for the output
chosen_name = '_simpleplotting_testing' # Choose the folder name
output_folder = curr_date + output_folder_name
os.makedirs(overall_output_file_path + output_folder) # Create folder

# Keep text as text (not paths) in SVG/PDF/PS
mpl.rcParams['svg.fonttype'] = 'none'  # <-- critical for Illustrator
mpl.rcParams['pdf.fonttype'] = 42      # TrueType; helps if you also export PDF
mpl.rcParams['ps.fonttype']  = 42

# (Optional) pick a common font that Illustrator will have
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial']  # or 'Helvetica', etc.

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

def plotSize(df, columns, output_name, y_label = 'Clone size', plot = 'bar'): 
    # Custom x-axis tick labels
    n = get_n(df)  
    x_axis_ticks = ('E9.5', 'E10.5', 'E11.5')
    #     f"E9.5\n$\it{{n={n.loc[n['TM'] == 'E9.5', 'Number_of_clones'].iloc[0]}}}$\n$\it{{a={n.loc[n['TM'] == 'E9.5', 'Number_of_animals'].iloc[0]}}}$",
    #     f"E10.5\n$\it{{n={n.loc[n['TM'] == 'E10.5', 'Number_of_clones'].iloc[0]}}}$\n$\it{{a={n.loc[n['TM'] == 'E10.5', 'Number_of_animals'].iloc[0]}}}$",
    #     f"E11.5\n$\it{{n={n.loc[n['TM'] == 'E11.5', 'Number_of_clones'].iloc[0]}}}$\n$\it{{a={n.loc[n['TM'] == 'E11.5', 'Number_of_animals'].iloc[0]}}}$",
    # )

    # Ensure output folder exists
    size_output_path = overall_output_file_path + '/Clone size'
    os.makedirs(size_output_path, exist_ok=True)

    # clone size per TM
    df['Sum'] = df[columns].sum(axis = 1)
    output_path = f"{size_output_path}/{output_name}_{curr_date}"

    if plot == 'bar':
        cap.snsBarPlot(
            x_column='TM',
            y_column='Sum',
            df=df[df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])],
            x_axis_title='TM',
            y_axis_title=y_label,
            x_axis_ticks=x_axis_ticks,
            output_path=output_path,
            colour=ColPals.lightgrey,
            fig_size=(3.5, 5),
            order = ['E9.5', 'E10.5', 'E11.5']
        )
    elif plot == 'violin':
        cap.snsViolinPlot(
        x_column='TM',
        y_column='Sum',
        df=df[df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])],
        x_axis_title='TM',
        y_axis_title=y_label,
        x_axis_ticks=x_axis_ticks,
        output_path=output_path,
        colour=ColPals.lightgrey,
        fig_size=(3.5, 5),
        order=['E9.5', 'E10.5', 'E11.5'],      
        max_y=None
    )
    saveDataTableGroupMeans(df = df,
              columns = ['Clone', 'TM', 'Sum'],
              dep_var = 'Sum', 
              ind_var = 'TM', 
              output_path = output_path,
              order = ['E9.5', 'E10.5', 'E11.5'])
    # Save data to CSV


def extract_animal_id(clone_str: str) -> str:
    # Gets Animal from Clone ID: e.g. S107_a7_c3 returns S107_a7
    parts = str(clone_str).split("_")
    if len(parts) >= 3:
        return "_".join(parts[:2])

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

# Save three-sheet Excel for proportions data:
#   Sheet 1: selected columns from df2 (clone-level), incl. Animal
#   Sheet 2: df1 summary (e.g. TM × Type × Proportions)
#   Sheet 3: counts (N_clones, N_animals) from df2 grouped by group_cols
#
# df1         : summary table (e.g. output of caf.cloneTypeProportions)
# df2         : clone-level table (e.g. rg_clones_df_glia_no_unicolor)
# df2_columns : list of columns from df2 to put on sheet 1
# group_cols  : columns to group on for counts, e.g. ['TM', 'Type']
# output_path : path *without* extension (".xlsx" will be appended)
# ind_var     : optional main independent variable for ordering (e.g. 'TM')
# order       : optional categorical order for ind_var (e.g. ['E9.5','E10.5','E11.5'])
def saveDataTableProps(df1, df2, df2_columns, group_cols, output_path,
                          ind_var=None, order=None):

    # Sheet 1: clone-level data from df2
    df2_clones = df2.copy()

    # Derive Animal from Clone
    if 'Clone' not in df2_clones.columns:
        raise KeyError("df2 must contain a 'Clone' column.")
    df2_clones['Animal'] = df2_clones['Clone'].apply(extract_animal_id)

    # Ensure required columns are included
    extra_needed = set(group_cols + ['Clone', 'Animal'])
    cols_sheet1 = list(extra_needed.union(df2_columns))

    # Re-order columns
    ordered_cols = []
    for c in ['Animal', 'Clone'] + group_cols + df2_columns:
        if c in cols_sheet1 and c not in ordered_cols:
            ordered_cols.append(c)
    for c in cols_sheet1:
        if c not in ordered_cols:
            ordered_cols.append(c)

    df2_clones = df2_clones[ordered_cols]

    # Sheet 2
    df1_summary = df1.copy()
    if ind_var is not None and order is not None and ind_var in df1_summary.columns:
        df1_summary[ind_var] = pd.Categorical(df1_summary[ind_var],
                                              categories=order,
                                              ordered=True)
        df1_summary = df1_summary.sort_values(ind_var).reset_index(drop=True)

    # Sheet 3 (number of clones and animals)
    # Reuse Animal column from df2_clones
    for gc in group_cols:
        if gc not in df2_clones.columns:
            raise KeyError(f"group column '{gc}' not found in df2.")

    # Write Excel with 3 sheets
    excel_path = output_path + ".xlsx"
    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        df2_clones.to_excel(writer, sheet_name="Clones", index=False)
        df1_summary.to_excel(writer, sheet_name="Summary", index=False)
        
######################################################################################################################################################################
### Plots related to clone size

### Total cells

plotSize(rg_clones_df, ['Total cells'], output_name = 'Total_cells_per_clone', plot = 'violin')

# Ensure output folder exists
size_output_path = overall_output_file_path + '/Clone size'
os.makedirs(size_output_path, exist_ok=True)

# Histogram of clone size
size_array = np.array(rg_clones_df['Total cells']) # Get an array from the 'Total cells' column
output_path = size_output_path + '/Clone_Size_Hist_' + curr_date 
cap.snsHistPlot(array = size_array, 
        x_axis_title = 'Clone size',
        y_axis_title = 'Counts',
        output_path = output_path,
        colour = 'k',
        binwidth = 5,
        binrange = [0,100],
        fig_size = (6,5))

######################################################################################################################################################################
### Plots related to clone spread

general_output_path = overall_output_file_path + '/General'
os.makedirs(general_output_path, exist_ok=True)

# Transverse intra-clone spread
tm_order = ['E9.5', 'E10.5', 'E11.5']
output_path = general_output_path + '/Intra-clone_spread_xy' + curr_date
cap.snsViolinPlot(
    x_column='TM',
    y_column='Spread XY',
    df = rg_clones_df[rg_clones_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
    x_axis_title='TM',
    y_axis_title='Transverse intra-clone spread (\u03bcm)',
    x_axis_ticks=['E9.5', 'E10.5', 'E11.5'],
    output_path=output_path,
    colour=ColPals.lightgrey,
    fig_size=(3, 5),
    order=['E9.5', 'E10.5', 'E11.5'],      
    max_y=None)

saveDataTableGroupMeans(df = rg_clones_df,
              columns = ['Clone', 'TM', 'Spread XY'],
              dep_var = 'Spread XY', 
              ind_var = 'TM', 
              output_path = output_path,
              order = ['E9.5', 'E10.5', 'E11.5'])


# Rostrocaudal intra-clone spread
tm_order = ['E9.5', 'E10.5', 'E11.5']
output_path = general_output_path + '/Intra-clone_spread_z' + curr_date
cap.snsViolinPlot(
    x_column='TM',
    y_column='Spread Z',
    df = rg_clones_df[rg_clones_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
    x_axis_title='TM',
    y_axis_title='Rostrocaudal intra-clone spread (\u03bcm)',
    x_axis_ticks=['E9.5', 'E10.5', 'E11.5'],
    output_path=output_path,
    colour=ColPals.lightgrey,
    fig_size=(3, 5),
    order=['E9.5', 'E10.5', 'E11.5'],      
    max_y=None)

saveDataTableGroupMeans(df = rg_clones_df,
              columns = ['Clone', 'TM', 'Spread Z'],
              dep_var = 'Spread Z', 
              ind_var = 'TM', 
              output_path = output_path,
              order = ['E9.5', 'E10.5', 'E11.5'])

### DV distribution clones
rg_clones_dv_prop_df = caf.cloneTypeProportions(rg_clones_df, ['TM'], ['DV'])
output_path = general_output_path + '/DV_distribution' + curr_date
cap.snsStackedBarPlot(data = rg_clones_dv_prop_df, 
                  x = 'TM', 
                  hue = 'DV',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = None,
                  output_path = output_path,
                  palette = ColPals.stckpal_dv,
                  fig_size = (3,5),
                  hue_order = ['Ventral', 'Mixed', 'Dorsal'])

saveDataTableProps(df1 = rg_clones_dv_prop_df,
                   df2 = rg_clones_df, 
                   df2_columns = ['Clone', 'TM', 'DV'], 
                   group_cols = ['TM', 'DV'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['Dorsal', 'Mixed', 'Ventral'])
######################################################################################################################################################################
### Plots related to neurons

# Ensure output folder exists
neuron_output_path = overall_output_file_path + '/Neurons'
os.makedirs(neuron_output_path, exist_ok=True)

### Neurons in neuron-containing clones
# Subset data
rg_clones_df_neurons = rg_clones_df[rg_clones_df['Neuron count (incl Comm,MN)'] > 0]
plotSize(rg_clones_df_neurons, ['Neuron count (incl Comm,MN)'], output_name = 'Neurons_per_clone', y_label = 'Neurons per clone', plot = 'violin')

### Proportion of unicolor clones
rg_clones_df_neurons_unicolor = rg_clones_df_neurons.copy()
rg_clones_df_neurons_unicolor['Type'] = rg_clones_df_neurons_unicolor['Type'].replace(['Asymmetric', 'Small neurogenic', 'Terminal neurogenic', 'Symmetric'], 'R/G')
print('\nrg_clones_df_neurons_unicolor\n', rg_clones_df_neurons_unicolor)
rg_clones_df_neurons_unicolor_df = caf.cloneTypeProportions(rg_clones_df_neurons_unicolor, ['TM'], ['Type'])
print('\nrg_clones_df_neurons_unicolor_df\n', rg_clones_df_neurons_unicolor_df)
output_path = neuron_output_path + '/Proportion of unicolor clones' + curr_date

cap.snsStackedBarPlot(data = rg_clones_df_neurons_unicolor_df[rg_clones_df_neurons_unicolor_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
                  x = 'TM', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['E9.5', 'E10.5', 'E11.5'],
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (3,5),
                  hue_order = ['Unicolor', 'R/G'])

saveDataTableProps(df1 = rg_clones_df_neurons_unicolor_df,
                   df2 = rg_clones_df_neurons_unicolor, 
                   df2_columns = ['Clone', 'TM', 'Type'], 
                   group_cols = ['TM', 'Type'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['E9.5', 'E10.5', 'E11.5'])

### Type of neuron-containing clones (excluding unicolor clones)
rg_clones_df_neurons_no_unicolor = rg_clones_df_neurons[(rg_clones_df_neurons['Type'].notna()) & (rg_clones_df_neurons['Type'] != 'Unicolor')] # Exclude unicolor clones
neuron_clone_type_df = caf.cloneTypeProportions(rg_clones_df_neurons_no_unicolor, ['TM'], ['Type'])
output_path = neuron_output_path + '/Type of neuron-containing clones' + curr_date

cap.snsStackedBarPlot(data = neuron_clone_type_df[neuron_clone_type_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
                  x = 'TM', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = 'TM',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['E9.5', 'E10.5', 'E11.5'],
                  output_path = output_path,
                  palette = ColPals.stckpal_4col,
                  fig_size = (3,5),
                  hue_order = ['Terminal neurogenic', 'Small neurogenic', 'Asymmetric', 'Symmetric'])

saveDataTableProps(df1 = neuron_clone_type_df,
                   df2 = rg_clones_df_neurons_no_unicolor, 
                   df2_columns = ['Clone', 'TM', 'Type'], 
                   group_cols = ['TM', 'Type'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['E9.5', 'E10.5', 'E11.5'])

# Histogram of neuron output size
neuron_size_array = np.array(rg_clones_df_neurons['Neuron count (incl Comm,MN)']) # Get an array from the 'Total cells' column
output_path = neuron_output_path + '/Neuron_output_Hist_' + curr_date 
cap.snsHistPlot(array = neuron_size_array, 
        x_axis_title = 'Neuron output',
        y_axis_title = 'Counts',
        output_path = output_path,
        colour = 'k',
        binwidth = 5,
        binrange = [0,80],
        fig_size = (5.3,5))

### Size of asymmetric neurogenic clones (asymmetric clones containing at least one neuron)
# Subset data
rg_clones_df_neurons_asymmetric = rg_clones_df_neurons[rg_clones_df_neurons['Type'] == 'Asymmetric']

plotSize(rg_clones_df_neurons_asymmetric, ['Total cells'], output_name = 'Neurons_per_asymmetric_neurogenic_clone', y_label = 'Size of asymmetric neurogenic clones', plot = 'violin')

### Size of symmetric neurogenic clones (symmetric clones containing at least one neuron)
rg_clones_df_neurons_symmetric = rg_clones_df_neurons[rg_clones_df_neurons['Type'] == 'Symmetric']
plotSize(rg_clones_df_neurons_symmetric, ['Total cells'], output_name = 'Neurons_per_symmetric_neurogenic_clone', y_label = 'Size of symmetric neurogenic clones', plot = 'violin')

### Size of small neurogenic clones (small neurogenic clones containing at least one neuron)
rg_clones_df_neurons_smallneurogenic = rg_clones_df_neurons[rg_clones_df_neurons['Type'] == 'Small neurogenic']
plotSize(rg_clones_df_neurons_smallneurogenic, ['Total cells'], output_name = 'Neurons_per_small_neurogenic_clone', y_label = 'Size of small neurogenic clones', plot = 'violin')

### Unicolor clone proportion
uni_comp_df = rg_clones_df.copy()
uni_comp_df['Unicolor'] = np.where(
    pd.to_numeric(uni_comp_df['Type'] == 'Unicolor'),
    'Unicolor',
    'R/G'
)

uni_comp_prop_df = caf.cloneTypeProportions(uni_comp_df, ['TM'], ['Unicolor'])

output_path = neuron_output_path + '/Unicolor_proportion' + curr_date

cap.snsStackedBarPlot(data = uni_comp_prop_df[uni_comp_prop_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
                  x = 'TM', 
                  hue = 'Unicolor',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['E9.5', 'E10.5', 'E11.5'],
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (3,5),
                  hue_order = ['R/G', 'Unicolor'])

saveDataTableProps(df1 = uni_comp_prop_df,
                   df2 = uni_comp_df, 
                   df2_columns = ['Clone', 'TM', 'Unicolor'], 
                   group_cols = ['TM', 'Unicolor'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['E9.5', 'E10.5', 'E11.5'])

######################################################################################################################################################################
### Plots related to glia

# Ensure output folder exists
glia_output_path = overall_output_file_path + '/Glia'
os.makedirs(glia_output_path, exist_ok=True)

### Proportion of clones containing glia (PA,FA,WF,GL,OG,OL)
glia_comp_df = rg_clones_df.copy()
glia_comp_df['Glia status'] = np.where(
    pd.to_numeric(glia_comp_df['Glia count (PA,FA,WF,GL,OG,OL)'], errors='coerce') > 0,
    'Glia-containing',
    'Non-glia-containing'
)

glia_comp_prop_df = caf.cloneTypeProportions(glia_comp_df, ['TM'], ['Glia status'])

output_path = glia_output_path + '/Glia-containing' + curr_date

cap.snsStackedBarPlot(data = glia_comp_prop_df[glia_comp_prop_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
                  x = 'TM', 
                  hue = 'Glia status',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['E9.5', 'E10.5', 'E11.5'],
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (3.3,5),
                  hue_order = ['Glia-containing', 'Non-glia-containing'])

saveDataTableProps(df1 = glia_comp_prop_df,
                   df2 = glia_comp_df, 
                   df2_columns = ['Clone', 'TM', 'Glia status'], 
                   group_cols = ['TM', 'Glia status'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['E9.5', 'E10.5', 'E11.5'])

### Number of glia per clone in glia-containing clones (PA,FA,WF,GL,OG,OL)
# Subset data
rg_clones_df_glia = rg_clones_df[rg_clones_df['Glia count (PA,FA,WF,GL,OG,OL)'] > 0] # glia-containing clones

# Bar graph
plotSize(rg_clones_df_glia, ['Glia count (PA,FA,WF,GL,OG,OL)'], output_name = 'Glia_per_clone', y_label = 'Glia per clone', plot = 'bar')

# Violin plot
plotSize(rg_clones_df_glia, ['Glia count (PA,FA,WF,GL,OG,OL)'], output_name = 'Glia_per_clone', y_label = 'Glia per clone', plot = 'violin')

### Type of glia-containing clones (excluding unicolor clones)
rg_clones_df_glia_no_unicolor = rg_clones_df_glia[(rg_clones_df_glia['Type'].notna()) & (rg_clones_df_glia['Type'] != 'Unicolor')] # Exclude unicolor clones
glia_clone_type_df = caf.cloneTypeProportions(rg_clones_df_glia_no_unicolor, ['TM'], ['Type'])
output_path = glia_output_path + '/Type of glia-containing clones' + curr_date

cap.snsStackedBarPlot(data = glia_clone_type_df[glia_clone_type_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
                  x = 'TM', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['E9.5', 'E10.5', 'E11.5'],
                  output_path = output_path,
                  palette = ColPals.stckpal_4col,
                  fig_size = (6,5),
                  hue_order = ['Terminal neurogenic', 'Small neurogenic', 'Asymmetric', 'Symmetric'])

saveDataTableProps(df1 = glia_clone_type_df,
                   df2 = rg_clones_df_glia_no_unicolor, 
                   df2_columns = ['Clone', 'TM', 'Type'], 
                   group_cols = ['TM', 'Type'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['E9.5', 'E10.5', 'E11.5'])

### Proportions of different types of glia in glia-containing clones
glia_cols = {
    "FA":    "FA count",
    "PA":    "PA count",
    "GL":    "GL count",
    "WF":    "WF count",
    "OG":    "OG count",
    "OL":    "OL count"
}

output_path = glia_output_path + '/Glia composition' + curr_date
cap.cellTypeHist(rg_clones_df_glia, glia_cols, output_path)

### Fraction of all asymmetric clones (regardless of time point) that contain glia (strictly glia, not progenitors or glia with endfeet)
rg_clones_df_asymmetric = rg_clones_df[rg_clones_df['Type'] == 'Asymmetric']
num_asym_glia = len(rg_clones_df_asymmetric[(rg_clones_df_asymmetric['Astrocyte count (PA,FA,GL)'] > 0) |
                                        (rg_clones_df_asymmetric['OG count'] > 0) |
                                        (rg_clones_df_asymmetric['OL count'] > 0)])
total_num_asymm = len(rg_clones_df_asymmetric)
prop_asym_glia = num_asym_glia / total_num_asymm
prop_asym_no_glia = 1 - prop_asym_glia
asym_prop_df = pd.DataFrame({
    "TM": ["combined", "combined"],
    "Type": ["with glia", "no glia"],
    "Proportions": [prop_asym_glia, prop_asym_no_glia],
    'Counts': [num_asym_glia, (total_num_asymm - num_asym_glia)],
    'SE': [caf.SECalculation(prop_asym_glia, total_num_asymm), caf.SECalculation(prop_asym_no_glia, total_num_asymm)]})

print('\nasym_prop_df\n', asym_prop_df)
output_path = glia_output_path + '/Proportion of all asymmetric clones that contain glia' + curr_date

cap.snsStackedBarPlot(data = asym_prop_df, 
                  x = 'TM', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion of asymmetric clones',
                  x_axis_ticks = None,
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (0.9,5),
                  hue_order = ['with glia', 'no glia'])

saveDataTableProps(df1 = asym_prop_df,
                   df2 = rg_clones_df_asymmetric, 
                   df2_columns = ['Clone', 'TM', 'Type'], 
                   group_cols = ['TM', 'Type'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = None)

### Fraction of neuron-containing asymmetric clones (regardless of time point) that contain glia (strictly glia, not progenitors or glia with endfeet)
rg_clones_df_asymmetric_neur = rg_clones_df_asymmetric[rg_clones_df_asymmetric['Neuron count (incl Comm,MN)'] > 0]
num_asym_neur_glia = len(rg_clones_df_asymmetric_neur[(rg_clones_df_asymmetric_neur['Astrocyte count (PA,FA,GL)'] > 0) |
                                        (rg_clones_df_asymmetric_neur['OG count'] > 0) |
                                        (rg_clones_df_asymmetric_neur['OL count'] > 0)])
total_num_neur_asymm = len(rg_clones_df_asymmetric_neur)
prop_asym_neur_glia = num_asym_neur_glia / len(rg_clones_df_asymmetric_neur)
prop_asym_neur_no_glia = 1 - prop_asym_neur_glia
asym_neur_prop_df = pd.DataFrame({
    "TM": ["combined", "combined"],
    "Type": ["with glia", "no glia"],
    "Proportions": [prop_asym_neur_glia, prop_asym_neur_no_glia],
    'Counts': [num_asym_neur_glia, (total_num_asymm - total_num_neur_asymm)],
    'SE': [caf.SECalculation(prop_asym_neur_glia, total_num_asymm), caf.SECalculation(prop_asym_neur_no_glia, total_num_asymm)]})

output_path = glia_output_path + '/Proportion of neuron-containing asymmetric clones that contain glia' + curr_date

cap.snsStackedBarPlot(data = asym_neur_prop_df, 
                  x = 'TM', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion of asymmetric clones',
                  x_axis_ticks = None,
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (0.9,5),
                  hue_order = ['with glia', 'no glia'])

saveDataTableProps(df1 = asym_neur_prop_df,
                   df2 = rg_clones_df_asymmetric_neur, 
                   df2_columns = ['Clone', 'TM', 'Type'], 
                   group_cols = ['TM', 'Type'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = None)

### DV distribution of glia-containing clones
glia_clone_dv_df = caf.cloneTypeProportions(rg_clones_df_glia, ['TM'], ['DV'])
output_path = glia_output_path + '/DV distribution of glia-containing clones' + curr_date

cap.snsStackedBarPlot(data = glia_clone_dv_df[glia_clone_dv_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
                  x = 'TM', 
                  hue = 'DV',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['E9.5', 'E10.5', 'E11.5'],
                  output_path = output_path,
                  palette = ColPals.stckpal_dv,
                  fig_size = (4.5,5),
                  hue_order = ['Ventral', 'Mixed', 'Dorsal'])

saveDataTableProps(df1 = glia_clone_dv_df,
                   df2 = rg_clones_df_glia, 
                   df2_columns = ['Clone', 'TM', 'DV'], 
                   group_cols = ['TM', 'DV'], 
                   output_path = output_path,
                   ind_var = 'DV', 
                   order = ['Dorsal', 'Mixed', 'Ventral'])

### Proportions of clones containing neurons-only, glia-only, or neurons and glia
# Output folder
ng_output_path = glia_output_path + '/Neuron_Glia_overlap'
os.makedirs(ng_output_path, exist_ok=True)

ng_df = rg_clones_df.copy()

# Define "has neuron" and "has glia"
ng_df['Has_neuron'] = pd.to_numeric(ng_df['Neuron count (incl Comm,MN)'], errors='coerce').fillna(0) > 0
ng_df['Has_glia']   = pd.to_numeric(ng_df['Glia count (PA,FA,WF,GL,OG,OL)'], errors='coerce').fillna(0) > 0

def ng_class(row):
    if row['Has_neuron'] and row['Has_glia']:
        return 'Neuron+Glia'
    elif row['Has_neuron'] and (not row['Has_glia']):
        return 'Neuron only'
    elif (not row['Has_neuron']) and row['Has_glia']:
        return 'Glia only'
    else:
        return 'Neither'  # optional; you can drop these

ng_df['Neuron/Glia class'] = ng_df.apply(ng_class, axis=1)

# If you want ONLY clones that have at least one of the two:
ng_df = ng_df[ng_df['Neuron/Glia class'] != 'Neither'].copy()

# Proportions per TM
ng_prop_df = caf.cloneTypeProportions(ng_df, ['TM'], ['Neuron/Glia class'])

# Enforce TM order (important for plotting + Excel ordering)
tm_order = ['E9.5', 'E10.5', 'E11.5']
ng_prop_df['TM'] = pd.Categorical(ng_prop_df['TM'], categories=tm_order, ordered=True)
ng_prop_df = ng_prop_df.sort_values('TM').reset_index(drop=True)

# Plot
output_path = ng_output_path + '/Neuron_Glia_proportions_' + curr_date

cap.snsStackedBarPlot(
    data=ng_prop_df[ng_prop_df['TM'].isin(tm_order)],
    x='TM',
    hue='Neuron/Glia class',
    weights='Proportions',
    x_axis_title='',
    y_axis_title='Proportion of clones',
    x_axis_ticks=tm_order,
    output_path=output_path,
    palette=ColPals.stckpal_4col,  # or another palette you prefer
    fig_size=(3, 5),
    hue_order=['Neither', 'Neuron+Glia', 'Glia only', 'Neuron only']
)

# Save data + counts (clones + animals) to Excel (3 sheets)
saveDataTableProps(
    df1=ng_prop_df,
    df2=ng_df,
    df2_columns=['Clone', 'TM', 'Neuron count (incl Comm,MN)', 'Glia count (PA,FA,WF,GL,OG,OL)', 'Neuron/Glia class'],
    group_cols=['TM', 'Neuron/Glia class'],
    output_path=output_path,
    ind_var='TM',
    order=tm_order
)

### Proportions of clones containing neurons-only, glia-only, or neurons and glia ONLY CLONES WITH COMPLETE MORPHOLOGICAL IDENTIFICATION
# Output folder
ng_output_path = glia_output_path + '/Neuron_Glia_overlap'
os.makedirs(ng_output_path, exist_ok=True)

ng_df = rg_clones_df.copy()
ng_df = ng_df[ng_df['Uncertain count'] == 0]

# Define "has neuron" and "has glia"
ng_df['Has_neuron'] = pd.to_numeric(ng_df['Neuron count (incl Comm,MN)'], errors='coerce').fillna(0) > 0
ng_df['Has_glia']   = pd.to_numeric(ng_df['Glia count (PA,FA,WF,GL,OG,OL)'], errors='coerce').fillna(0) > 0

def ng_class(row):
    if row['Has_neuron'] and row['Has_glia']:
        return 'Neuron+Glia'
    elif row['Has_neuron'] and (not row['Has_glia']):
        return 'Neuron only'
    elif (not row['Has_neuron']) and row['Has_glia']:
        return 'Glia only'
    else:
        return 'Neither'  # optional; you can drop these

ng_df['Neuron/Glia class'] = ng_df.apply(ng_class, axis=1)

# If you want ONLY clones that have at least one of the two:
ng_df = ng_df[ng_df['Neuron/Glia class'] != 'Neither'].copy()

# Proportions per TM
ng_prop_df = caf.cloneTypeProportions(ng_df, ['TM'], ['Neuron/Glia class'])

# Enforce TM order (important for plotting + Excel ordering)
tm_order = ['E9.5', 'E10.5', 'E11.5']
ng_prop_df['TM'] = pd.Categorical(ng_prop_df['TM'], categories=tm_order, ordered=True)
ng_prop_df = ng_prop_df.sort_values('TM').reset_index(drop=True)

# Plot
output_path = ng_output_path + '/Neuron_Glia_proportions_only_complete_clones' + curr_date

cap.snsStackedBarPlot(
    data=ng_prop_df[ng_prop_df['TM'].isin(tm_order)],
    x='TM',
    hue='Neuron/Glia class',
    weights='Proportions',
    x_axis_title='',
    y_axis_title='Proportion of clones',
    x_axis_ticks=tm_order,
    output_path=output_path,
    palette=ColPals.stckpal_4col,  # or another palette you prefer
    fig_size=(3, 5),
    hue_order=['Neither', 'Neuron+Glia', 'Glia only', 'Neuron only']
)

# Save data + counts (clones + animals) to Excel (3 sheets)
saveDataTableProps(
    df1=ng_prop_df,
    df2=ng_df,
    df2_columns=['Clone', 'TM', 'Neuron count (incl Comm,MN)', 'Glia count (PA,FA,WF,GL,OG,OL)', 'Neuron/Glia class'],
    group_cols=['TM', 'Neuron/Glia class'],
    output_path=output_path,
    ind_var='TM',
    order=tm_order
)
######################################################################################################################################################################
### Plots related to cell-type composition

# Ensure output folder exists
comp_output_path = overall_output_file_path + '/Cell-type composition'
os.makedirs(comp_output_path, exist_ok=True)

### General cell types including uncertain cells
general_cols = {
    "N":    "Neuron count (incl Comm,MN)",
    "A":    "Astrocyte count (PA,FA,GL)",
    "WF":   "WF count",
    "OG":   "OG count",
    "OL":   "OL count",
    "PR":   "PR count",
    'MC':   'Midline count (dv)',
    "U":    "Uncertain count"
}

output_path = comp_output_path + '/General cell types incl uncertain' + curr_date
cap.cellTypeHist(rg_clones_df, general_cols, output_path)

### General cell types including excluding cells
general_cols = {
    "N":    "Neuron count (incl Comm,MN)",
    "A":    "Astrocyte count (PA,FA,GL)",
    "WF":   "WF count",
    "OG":   "OG count",
    "OL":   "OL count",
    "PR":   "PR count",
    'MC':   'Midline count (dv)'
}

output_path = comp_output_path + '/General cell types' + curr_date
cap.cellTypeHist(rg_clones_df, general_cols, output_path)

### Neuron types in neuron-containing clones
neuron_cols = {
    "N":    "Neuron count (incl Comm,MN)",
    "DCN":  "Commissural d count",
    "VCN":  "Commissural v count",
    "MN":   "MN count"
}

output_path = comp_output_path + '/Neuron types' + curr_date
cap.cellTypeHist(rg_clones_df_neurons, neuron_cols, output_path)

### Proportion of clones containing 'Uncertain' cells
uncert_comp_df = rg_clones_df.copy()
uncert_comp_df['Uncertain status'] = np.where(
    pd.to_numeric(uncert_comp_df['Uncertain count'], errors='coerce') > 0,
    'Uncertain-containing',
    'Non-uncertain-containing'
)

uncert_comp_prop_df = caf.cloneTypeProportions(uncert_comp_df, ['TM'], ['Uncertain status'])

output_path = comp_output_path + '/Uncertain-containing' + curr_date
cap.snsStackedBarPlot(data = uncert_comp_prop_df[uncert_comp_prop_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
                  x = 'TM', 
                  hue = 'Uncertain status',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['E9.5', 'E10.5', 'E11.5'],
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (7,5),
                  hue_order = ['Uncertain-containing', 'Non-uncertain-containing'])

saveDataTableProps(df1 = uncert_comp_prop_df,
                   df2 = uncert_comp_df, 
                   df2_columns = ['Clone', 'TM', 'Uncertain status'], 
                   group_cols = ['TM', 'Uncertain status'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['Non-uncertain-containing', 'Uncertain-containing'])

### Proportion of total cells that are 'Uncertain' and average clone coverage
prop_uncert = rg_clones_df['Uncertain count'].sum() / rg_clones_df['Total cells'].sum() # Fraction of overall cells that are uncertain
per_clone_cov = rg_clones_df['Uncertain count'] / rg_clones_df['Total cells'] # Fraction of each clone that is uncertain
clone_cov = per_clone_cov.mean()
# Write to Excel file
summary_df = pd.DataFrame({
    'Metric': ['Fraction of total cells that are uncertain', 'Avg fraction of clones that is uncertain (avg [# uncertain cells/total cells])'],
    'Value': [prop_uncert, clone_cov]
})

excel_path = comp_output_path +'/' + curr_date + '_uncertain_cells_summary.xlsx'

with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
    summary_df.to_excel(writer, sheet_name="Summary", index=False)
    per_clone_cov.to_excel(writer, sheet_name="Uncertain fraction per clone", index=False)

### Proportion of clones containing commissural neurons per TM
rg_clones_comm_df = rg_clones_df.copy()
rg_clones_comm_df["Comm_status"] = np.where(rg_clones_comm_df["Commissural count_(dv)"] > 0, "comm", "no comm")
rg_clones_comm_prop_df = caf.cloneTypeProportions(rg_clones_comm_df, ['TM'], ['Comm_status'])
# Re-organize the order of the TM timepoints - important for snsStackedBarPlot
tm_order = ['E9.5', 'E10.5', 'E11.5']
rg_clones_comm_prop_df['TM'] = pd.Categorical(rg_clones_comm_prop_df['TM'], categories=tm_order, ordered=True) # Convert 'TM' column to a categorical type with the desired order
rg_clones_comm_prop_df = rg_clones_comm_prop_df.sort_values('TM').reset_index(drop=True) # Sort the DataFrame by the 'TM' column
output_path = comp_output_path + '/Commissural_neuron_proportion_TM' + curr_date
cap.snsStackedBarPlot(data = rg_clones_comm_prop_df, 
                  x = 'TM', 
                  hue = 'Comm_status',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['E9.5', 'E10.5', 'E11.5'],
                  output_path = output_path,
                  palette = ColPals.stckpal_lat,
                  fig_size = (3,5),
                  hue_order = ['comm', 'no comm'])

### Proportion of clones containing commissural neurons overall
num_comm = (rg_clones_comm_df['Comm_status'] == 'comm').sum()
total_num_clones = len(rg_clones_comm_df)
prop_comm = num_comm / total_num_clones
prop_no_comm = 1 - prop_comm
prop_comm_df = pd.DataFrame({
    "TM": ["combined", "combined"],
    "Type": ["with comm", "no comm"],
    "Proportions": [prop_comm, prop_no_comm],
    'Counts': [num_comm, (total_num_clones - num_comm)],
    'SE': [caf.SECalculation(prop_comm, total_num_clones), caf.SECalculation(prop_no_comm, total_num_clones)]})

output_path = comp_output_path + '/Commissural_neuron_proportion' + curr_date

cap.snsStackedBarPlot(data = prop_comm_df, 
                  x = 'TM', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion of clones',
                  x_axis_ticks = None,
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (0.6,5),
                  hue_order = ['with comm', 'no comm'])

saveDataTableProps(df1 = prop_comm_df,
                   df2 = rg_clones_comm_df, 
                   df2_columns = ['Clone', 'TM', 'Comm_status'], 
                   group_cols = ['TM', 'Comm_status'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['no comm', 'with comm'])



######################################################################################################################################################################
### Plots related to clone laterality

lat_output_path = overall_output_file_path + '/Laterality'
os.makedirs(lat_output_path, exist_ok=True)

### Proportion of bilateral clones. Clones not classified as uni or bi because they contain mainly midline cells (where Hemisphere = 'na') are not counted.
lat_rg_clones_df = rg_clones_df[rg_clones_df['Laterality'].fillna('') != '']
bilateral_prop_df = caf.cloneTypeProportions(lat_rg_clones_df, ['TM'], ['Laterality'])
# Re-organize the order of the TM timepoints - important for snsStackedBarPlot
tm_order = ['E9.5', 'E10.5', 'E11.5']
bilateral_prop_df['TM'] = pd.Categorical(bilateral_prop_df['TM'], categories=tm_order, ordered=True) # Convert 'TM' column to a categorical type with the desired order
bilateral_prop_df = bilateral_prop_df.sort_values('TM').reset_index(drop=True) # Sort the DataFrame by the 'TM' column
output_path = lat_output_path + '/Bilateral_proportion' + curr_date
cap.snsStackedBarPlot(data = bilateral_prop_df, 
                  x = 'TM', 
                  hue = 'Laterality',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['E9.5', 'E10.5', 'E11.5'],
                  output_path = output_path,
                  palette = ColPals.stckpal_lat,
                  fig_size = (3,5),
                  hue_order = ['bi', 'uni'])

saveDataTableProps(df1 = bilateral_prop_df,
                   df2 = lat_rg_clones_df, 
                   df2_columns = ['Clone', 'TM', 'Laterality'], 
                   group_cols = ['TM', 'Laterality'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['uni', 'bi'])

### DV distribution of bilateral clones
rg_clones_df_bilateral = rg_clones_df[rg_clones_df['Laterality'] == 'bi'] # Subset of bilateral clones
rg_clones_df_bilateral['TM'] = 'all' # Don't want to separate TM amounts

bilateral_dv_df = caf.cloneTypeProportions(rg_clones_df_bilateral, ['TM'], ['DV'])
bilateral_dv_df['TM'] = 'all'
output_path = lat_output_path + '/bilateral_DV_distribution' + curr_date
cap.snsStackedBarPlot(data = bilateral_dv_df, 
                  x = 'TM', 
                  hue = 'DV',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = None,
                  output_path = output_path,
                  palette = ColPals.stckpal_dv,
                  fig_size = (0.6,5),
                  hue_order = ['Ventral', 'Mixed', 'Dorsal'])

saveDataTableProps(df1 = bilateral_dv_df,
                   df2 = rg_clones_df_bilateral, 
                   df2_columns = ['Clone', 'TM', 'DV'], 
                   group_cols = ['TM', 'DV'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['Dorsal', 'Mixed', 'Ventral'])

### Plot the relative sizes of left and right hemisphere clone fragments (i.e. the number of cells on the side of the SC with the most cells vs the 
### number of cells on the side of the SC with the fewest cells)

# Organize the data for input into a seaborn line graph
le_list = rg_clones_df_bilateral['Total cells_le'].tolist()
ri_list = rg_clones_df_bilateral['Total cells_ri'].tolist()
frag_size_df = pd.DataFrame({
    "Left": le_list,
    "Right": ri_list,
})
frag_size_df = frag_size_df.T
output_path = lat_output_path + '/Bilateral_clone_fragment_size' + curr_date

cap.snsLinePlot(df = frag_size_df,
            x_axis_title = 'Hemisphere',
            y_axis_title = 'Number of cells',
            x_axis_ticks = ['Left', 'Right'],
            output_path = output_path,
            palette = [ColPals.purple],
            fig_size = (3.5,5))

### Size of bilateral clones vs unilateral clones
output_path = lat_output_path + '/Bilateral_vs_unilateral_clone_size' + curr_date
cap.snsViolinPlot(
    x_column='Laterality',
    y_column='Total cells',
    df = rg_clones_df[rg_clones_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
    x_axis_title='Laterality',
    y_axis_title='Clone size',
    x_axis_ticks=['Unilateral', 'Bilateral'],
    output_path=output_path,
    colour=ColPals.lightgrey,
    fig_size=(3.5, 5),
    order=['uni', 'bi'],      
    max_y=None)

saveDataTableGroupMeans(df = rg_clones_df,
              columns = ['Clone', 'Laterality', 'Total cells'],
              dep_var = 'Total cells', 
              ind_var = 'Laterality', 
              output_path = output_path,
              order = ['uni', 'bi'])

### Cell-type composition of left and right fragments
# Get a subset of the bilateral clones where neither hemisphere contains only 'uncertain' cells. We want to exclude bilateral clones with only uncertain cells on either side from all the cell-type related bilateral clone analyses.
df = rg_clones_df_bilateral.copy()

# Make sure NAs are treated as zero for the counts
count_cols = [
    'Neuron count (incl Comm,MN)_le','Neuron count (incl Comm,MN)_ri',
    'Glia count (PA,FA,WF,GL,OG,OL)_le','Glia count (PA,FA,WF,GL,OG,OL)_ri',
    'Commissural count_(dv)_le','Commissural count_(dv)_ri',
    'Uncertain count_le','Uncertain count_ri'
]
df[count_cols] = df[count_cols].fillna(0)

# Condition: side has ONLY Uncertain cells (and at least one Uncertain)
only_uncertain_le = (
    (df['Uncertain count_le'] > 0) &
    (df['Neuron count (incl Comm,MN)_le'] == 0) &
    (df['Glia count (PA,FA,WF,GL,OG,OL)_le'] == 0) &
    (df['Commissural count_(dv)_le'] == 0)
)

only_uncertain_ri = (
    (df['Uncertain count_ri'] > 0) &
    (df['Neuron count (incl Comm,MN)_ri'] == 0) &
    (df['Glia count (PA,FA,WF,GL,OG,OL)_ri'] == 0) &
    (df['Commissural count_(dv)_ri'] == 0)
)

# Remove rows where *either* side is only Uncertain
mask_bad = only_uncertain_le | only_uncertain_ri

rg_clones_df_bilateral_filtered = df[~mask_bad].copy()

# Left/right cell type composition
# Helper to safely fetch a column by trying multiple likely variants
def get_first_existing(row, candidates, default=0):
    for c in candidates:
        if c in row.index:
            val = row[c]
            return 0 if pd.isna(val) else val
    return default

# columns we will read (fill missing with 0 for safety)
count_cols = [
    'Total cells_le','Total cells_ri',
    'Neuron count (incl Comm,MN)_le','Neuron count (incl Comm,MN)_ri',
    'Glia count (PA,FA,WF,GL,OG,OL)_le','Glia count (PA,FA,WF,GL,OG,OL)_ri',
    'Commissural count_(dv)_le','Commissural count_(dv)_ri',
    'Uncertain count_le','Uncertain count_ri'
]
for c in count_cols:
    if c not in rg_clones_df_bilateral.columns:
        rg_clones_df_bilateral[c] = 0
rg_clones_df_bilateral[count_cols] = rg_clones_df_bilateral[count_cols].fillna(0)

def yesno(x: bool) -> str:
    return 'yes' if bool(x) else 'no'

records = []
for _, row in rg_clones_df_bilateral_filtered.iterrows():

    # Colour logic: check per-hemisphere red/green counts
    # Colour logic: check per-hemisphere red/green counts
    le_g = get_first_existing(row, ['Total cells_g_le','Total cells_le_g','Total_cells_le_g'], 0)
    le_r = get_first_existing(row, ['Total cells_r_le','Total cells_le_r','Total_cells_le_r'], 0)
    ri_g = get_first_existing(row, ['Total cells_g_ri','Total cells_ri_g','Total_cells_ri_g'], 0)
    ri_r = get_first_existing(row, ['Total cells_r_ri','Total cells_ri_r','Total_cells_ri_r'], 0)

    def colour_tags(g_cnt, r_cnt):
        g_pos = float(g_cnt) > 0
        r_pos = float(r_cnt) > 0

        # "simple" = uni/bi/none (backward compatible)
        if g_pos and r_pos:
            simple = 'bi'
        elif g_pos ^ r_pos:
            simple = 'uni'
        else:
            simple = 'none'

        # "detailed" = bi, g, r, none
        if g_pos and r_pos:
            detailed = 'bi'
        elif g_pos and not r_pos:
            detailed = 'g'      # uni-g
        elif r_pos and not g_pos:
            detailed = 'r'      # uni-r
        else:
            detailed = 'none'

        return simple, detailed

    le_colour, le_colour_type = colour_tags(le_g, le_r)
    ri_colour, ri_colour_type = colour_tags(ri_g, ri_r)

    # LEFT fragment row (based purely on *_le columns)
    records.append({
        'Clone': row['Clone'],
        'Fragment': 'Left',
        'Contains neuron': yesno(row['Neuron count (incl Comm,MN)_le'] > 0),
        'Contains glia': yesno(row['Glia count (PA,FA,WF,GL,OG,OL)_le'] > 0),
        'Contains commissural neuron': yesno(row['Commissural count_(dv)_le'] > 0),
        'Contains uncertain': yesno(row['Uncertain count_le'] > 0),
        'Colour': le_colour,              # uni/bi/none (old behaviour)
        'Colour_type': le_colour_type,    # bi / g / r / none (new)
    })

    # RIGHT fragment row (based purely on *_ri columns)
    records.append({
        'Clone': row['Clone'],
        'Fragment': 'Right',
        'Contains neuron': yesno(row['Neuron count (incl Comm,MN)_ri'] > 0),
        'Contains glia': yesno(row['Glia count (PA,FA,WF,GL,OG,OL)_ri'] > 0),
        'Contains commissural neuron': yesno(row['Commissural count_(dv)_ri'] > 0),
        'Contains uncertain': yesno(row['Uncertain count_ri'] > 0),
        'Colour': ri_colour,
        'Colour_type': ri_colour_type,
    })
# Final per-fragment dataframe
frag_df = pd.DataFrame.from_records(
    records,
    columns=[
        'Clone','Fragment',
        'Contains neuron','Contains glia',
        'Contains commissural neuron','Contains uncertain',
        'Colour','Colour_type'
    ]
)

# Order fragments Left, Right
frag_df['Fragment'] = pd.Categorical(
    frag_df['Fragment'],
    categories=['Left','Right'],
    ordered=True
)
frag_df = frag_df.sort_values(['Clone','Fragment']).reset_index(drop=True)

## Make a summary table to find the proportion of various bilateral cell type combinations. Hemispheres are just considered different sides, 'left' and 'right' don't matter
def fragment_type(row):
    flags = []
    if row['Contains neuron'] == 'yes':
        flags.append('N')   # neuron
    if row['Contains glia'] == 'yes':
        flags.append('G')   # glia
    if row['Contains commissural neuron'] == 'yes':
        flags.append('C')   # commissural
    if row['Contains uncertain'] == 'yes':
        flags.append('U')   # uncertain

    if not flags:
        return 'None'       # contains none of the 4 types
    return '+'.join(flags)  # e.g. "N", "G", "N+G", "N+G+U", ...
    
frag_df['FragType'] = frag_df.apply(fragment_type, axis=1)

# one row per clone: columns 'Left' and 'Right' holding FragType
per_clone_types = (
    frag_df
    .pivot(index='Clone', columns='Fragment', values='FragType')
    .reset_index()
)

def combo_key(row):
    # sort the two types so 'N & G' == 'G & N'
    types = sorted([row['Left'], row['Right']])
    return ' & '.join(types)

per_clone_types['Combo'] = per_clone_types.apply(combo_key, axis=1)

n_clones = len(per_clone_types)

combo_summary = (
    per_clone_types['Combo']
    .value_counts()
    .rename_axis('Left/Right combination')
    .reset_index(name='Count')
)

combo_summary['Proportion'] = combo_summary['Count'] / n_clones

### Make a summary table for cell-type composition regardless of hemisphere (all cells)
type_summary = rg_clones_df_bilateral_filtered.copy()

# Ensure all count columns exist and fill NAs
count_cols = [
    'Neuron count (incl Comm,MN)_le','Neuron count (incl Comm,MN)_ri',
    'Glia count (PA,FA,WF,GL,OG,OL)_le','Glia count (PA,FA,WF,GL,OG,OL)_ri',
    'Commissural count_(dv)_le','Commissural count_(dv)_ri',
    'Uncertain count_le','Uncertain count_ri'
]
type_summary[count_cols] = type_summary[count_cols].fillna(0)

# Combine L + R per clone
type_summary['Neuron_total']      = type_summary['Neuron count (incl Comm,MN)_le'] + type_summary['Neuron count (incl Comm,MN)_ri']
type_summary['Glia_total']        = type_summary['Glia count (PA,FA,WF,GL,OG,OL)_le'] + type_summary['Glia count (PA,FA,WF,GL,OG,OL)_ri']
type_summary['Comm_total']        = type_summary['Commissural count_(dv)_le']       + type_summary['Commissural count_(dv)_ri']
type_summary['Uncertain_total']   = type_summary['Uncertain count_le']              + type_summary['Uncertain count_ri']

def combined_comp(row):
    types = []
    if row['Neuron_total']     > 0: types.append('N')
    if row['Glia_total']       > 0: types.append('G')
    if row['Comm_total']       > 0: types.append('C')
    if row['Uncertain_total']  > 0: types.append('U')

    if not types:
        return 'None'
    return '+'.join(types)

type_summary['CombinedType'] = type_summary.apply(combined_comp, axis=1)

n_clones = len(type_summary)

combined_summary = (
    type_summary['CombinedType']
    .value_counts()
    .rename_axis('Cell-type composition (combined L+R)')
    .reset_index(name='Count')
)

combined_summary['Proportion'] = combined_summary['Count'] / n_clones

### Make a summary table for hemisphere r/g composition
colour_per_clone = (
    frag_df
    .pivot(index='Clone', columns='Fragment', values='Colour_type')
    .reset_index()
)

def colour_combo(row):
    # sort so, e.g., 'bi' + 'g' and 'g' + 'bi' collapse to the same combo
    types = sorted([row['Left'], row['Right']])
    return ' & '.join(types)

colour_per_clone['ColourCombo'] = colour_per_clone.apply(colour_combo, axis=1)

n_clones = len(colour_per_clone)

colour_summary = (
    colour_per_clone['ColourCombo']
    .value_counts()
    .rename_axis('Colour combination (order independent)')
    .reset_index(name='Count')
)

colour_summary['Proportion'] = colour_summary['Count'] / n_clones

# Save files
file_name = lat_output_path + '/_bilateral_clone_hemisphere_composition_' + curr_date + '.csv' # composition of each clone
frag_df.to_csv(file_name, mode = 'w', index = False) # Will overwrite existing file
file_name = lat_output_path + '/_bilateral_clone_hemisphere_composition_summary' + curr_date + '.csv' # hemisphere cell-type composition combinations
combo_summary.to_csv(file_name, mode = 'w', index = False) # Will overwrite existing file
file_name = lat_output_path + '/_bilateral_clone_overall_composition_summary' + curr_date + '.csv' # overall cell-type composition regardless of hemisphere
combined_summary.to_csv(file_name, mode = 'w', index = False) # Will overwrite existing file
file_name = lat_output_path + '/_bilateral_clone_hemisphere_colour_composition_summary' + curr_date + '.csv' # overall cell-type composition regardless of hemisphere
colour_summary.to_csv(file_name, mode = 'w', index = False) # Will overwrite existing file



### Proportion of bilateral clones with neurons or glia in both hemispheres
grouped_frag_df  = frag_df.groupby("Clone")
num_bi = len(grouped_frag_df)

# Number of clones where both majority and minority fragments have 'yes' in Contains neuron
bi_neuron = grouped_frag_df["Contains neuron"].apply(lambda s: (s == "yes").all()).sum()
bi_neuron_prop_df = pd.DataFrame({
    "Clones": ["bilateral", "bilateral"],
    "Type": ["Bilateral neurons", "Unilateral neurons"],
    "Proportions": [(bi_neuron/num_bi), (1 - bi_neuron/num_bi)],
    'Counts': [bi_neuron, (num_bi - bi_neuron)],
    'SE': [caf.SECalculation((bi_neuron/num_bi), num_bi), caf.SECalculation((1 - bi_neuron/num_bi), num_bi)]})
    
# Number of clones where both majority and minority fragments have 'yes' in Contains glia
bi_glia = grouped_frag_df["Contains glia"].apply(lambda s: (s == "yes").all()).sum()
bi_glia_prop_df = pd.DataFrame({
    "Clones": ["bilateral", "bilateral"],
    "Type": ["Bilateral glia", "Unilateral glia"],
    "Proportions": [(bi_glia/num_bi), (1 - bi_glia/num_bi)],
    'Counts': [bi_neuron, (num_bi - bi_glia)],
    'SE': [caf.SECalculation((bi_glia/num_bi), num_bi), caf.SECalculation((1 - bi_glia/num_bi), num_bi)]})

output_path = lat_output_path + '/Proportion of bilateral clones with bilateral neurons' + curr_date

cap.snsStackedBarPlot(data = bi_neuron_prop_df, 
                  x = 'Clones', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion of bilateral clones',
                  x_axis_ticks = None,
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (0.6,5),
                  hue_order = ["Bilateral neurons", "Unilateral neurons"])


output_path = lat_output_path + '/Proportion of bilateral clones with bilateral glia' + curr_date

cap.snsStackedBarPlot(data = bi_glia_prop_df, 
                  x = 'Clones', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion of bilateral clones',
                  x_axis_ticks = None,
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (0.6,5),
                  hue_order = ["Bilateral glia", "Unilateral glia"])

### Proportion of unilateral vs bilateral clones containing commissural neurons
# Proportion of unilateral clones containing commissural neurons
rg_clones_df_unilateral = rg_clones_df[rg_clones_df['Laterality'] == 'uni'] # Subset of unilateral clones
rg_clones_df_unilateral["Comm_status"] = np.where(rg_clones_df_unilateral["Commissural count_(dv)"] > 0, "comm", "no comm")
num_comm_uni = (rg_clones_df_unilateral['Comm_status'] == 'comm').sum()
total_num_uni = len(rg_clones_df_unilateral)
prop_comm_uni = num_comm_uni / total_num_uni
prop_no_comm_uni = 1 - prop_comm_uni
# Proportion of bilateral clones containing commissural neurons
rg_clones_df_bilateral["Comm_status"] = np.where(rg_clones_df_bilateral["Commissural count_(dv)"] > 0, "comm", "no comm")
num_comm_bi = (rg_clones_df_bilateral['Comm_status'] == 'comm').sum()
total_num_bi = len(rg_clones_df_bilateral)
prop_comm_bi = num_comm_bi / total_num_bi
prop_no_comm_bi = 1 - prop_comm_bi

prop_comm_unibi_df = pd.DataFrame({
    "Laterality": ["Unilateral", "Unilateral", "Bilateral", "Bilateral"],
    "Type": ["with comm", "no comm", "with comm", "no comm"],
    "Proportions": [prop_comm_uni, prop_no_comm_uni, prop_comm_bi, prop_no_comm_bi],
    'Counts': [num_comm_uni, (total_num_uni - num_comm_uni), num_comm_bi, (total_num_bi - num_comm_bi)],
    'SE': [caf.SECalculation(prop_comm_uni, total_num_uni), caf.SECalculation(prop_no_comm_uni, total_num_uni), caf.SECalculation(prop_comm_bi, total_num_bi), caf.SECalculation(prop_no_comm_bi, total_num_bi)]})

output_path = lat_output_path + '/Commissural_neuron_proportion_in_uni_and_bilateral' + curr_date

cap.snsStackedBarPlot(data = prop_comm_unibi_df, 
                  x = 'Laterality', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = '',
                  y_axis_title = 'Proportion of clones',
                  x_axis_ticks = ['Unilateral', 'Bilateral'],
                  output_path = output_path,
                  palette = ColPals.stckpal_2col,
                  fig_size = (1.2,5),
                  hue_order = ['with comm', 'no comm'])

saveDataTableProps(df1 = prop_comm_unibi_df,
                   df2 = pd.concat([rg_clones_df_unilateral, rg_clones_df_bilateral], axis = 0), 
                   df2_columns = ['Clone', 'Laterality', 'Comm_status'], 
                   group_cols = ['Laterality', 'Comm_status'], 
                   output_path = output_path,
                   ind_var = 'Laterality', 
                   order = ['no comm', 'with comm'])

######################################################################################################################################################################
### Plots related to comparisons between SC levels (RC)

rg_clones_levels_df = rg_clones_df.copy()
#rg_clones_levels_df['Level'] = rg_clones_levels_df['Level'].replace('sacral', 'limb')

levels_output_path = overall_output_file_path + '/Level Comparison'
os.makedirs(levels_output_path, exist_ok=True)

# Clone size per level (limb:s1,s2,s5,s6 vs thoracic: s3,s4) (all clones)
output_path = levels_output_path + '/Clone_size_per_level' + curr_date
cap.snsViolinPlot(
    x_column='Level',
    y_column='Total cells',
    df = rg_clones_levels_df[rg_clones_levels_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
    x_axis_title='Level',
    y_axis_title='Clone size',
    x_axis_ticks=['Limb', 'Thoracic', 'Sacral'],
    output_path=output_path,
    colour=ColPals.lightgrey,
    fig_size=(3.5, 5),
    order=['limb', 'thoracic', 'sacral'],      
    max_y=None)

# Symmetric clone size per level
rg_clones_levels_symm_df = rg_clones_levels_df[rg_clones_levels_df['Type'] == 'Symmetric']
output_path = levels_output_path + '/Symmetric_clone_size_per_level' + curr_date
cap.snsViolinPlot(
    x_column='Level',
    y_column='Total cells',
    df = rg_clones_levels_symm_df[rg_clones_levels_symm_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
    x_axis_title='Level',
    y_axis_title='Symmetric clone size',
    x_axis_ticks=['Limb', 'Thoracic', 'Sacral'],
    output_path=output_path,
    colour=ColPals.lightgrey,
    fig_size=(3.5, 5),
    order=['limb', 'thoracic', 'sacral'],      
    max_y=None)


# Asymmetric clone size per level
rg_clones_levels_asymm_df = rg_clones_levels_df[rg_clones_levels_df['Type'] == 'Asymmetric']
output_path = levels_output_path + '/Asymmetric_clone_size_per_level' + curr_date
cap.snsViolinPlot(
    x_column='Level',
    y_column='Total cells',
    df = rg_clones_levels_asymm_df[rg_clones_levels_asymm_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
    x_axis_title='Level',
    y_axis_title='Asymmetric clone size',
    x_axis_ticks=['Limb', 'Thoracic', 'Sacral'],
    output_path=output_path,
    colour=ColPals.lightgrey,
    fig_size=(3.5, 5),
    order=['limb', 'thoracic', 'sacral'],      
    max_y=None)

# Clone type per level
output_path = levels_output_path + '/Clone_type_per_level' + curr_date
rg_clones_levels_df_no_unicolor = rg_clones_levels_df[(rg_clones_levels_df['Type'].notna()) & (rg_clones_levels_df['Type'] != 'Unicolor')] # Exclude unicolor clones
rg_clones_levels_type_df = caf.cloneTypeProportions(rg_clones_levels_df_no_unicolor, ['Level'], ['Type'], sort_TM = False, sort_level = True)

cap.snsStackedBarPlot(data = rg_clones_levels_type_df, 
                  x = 'Level', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = 'Level',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['Limb', 'Thoracic', 'Sacral'],
                  output_path = output_path,
                  palette = ColPals.stckpal_4col,
                  fig_size = (3.5,5),
                  hue_order = ['Symmetric', 'Asymmetric', 'Small neurogenic', 'Terminal neurogenic'])

# Clone type per level per TM
# E9.5
output_path = levels_output_path + '/Clone_type_per_level_E9' + curr_date
rg_clones_levels_type_df_E9 = rg_clones_levels_df_no_unicolor[rg_clones_levels_df_no_unicolor['TM'] == 'E9.5']
rg_clones_levels_type_df_E9 = caf.cloneTypeProportions(rg_clones_levels_type_df_E9, ['Level'], ['Type'], sort_TM = False, sort_level = True)

cap.snsStackedBarPlot(data = rg_clones_levels_type_df_E9, 
                  x = 'Level', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = 'Level',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['Limb', 'Thoracic', 'Sacral'],
                  output_path = output_path,
                  palette = ColPals.stckpal_4col,
                  fig_size = (3.5,5),
                  hue_order = ['Symmetric', 'Asymmetric', 'Small neurogenic', 'Terminal neurogenic'])

# E10.5
output_path = levels_output_path + '/Clone_type_per_level_E10' + curr_date
rg_clones_levels_type_df_E10 = rg_clones_levels_df_no_unicolor[rg_clones_levels_df_no_unicolor['TM'] == 'E10.5']
rg_clones_levels_type_df_E10 = caf.cloneTypeProportions(rg_clones_levels_type_df_E10, ['Level'], ['Type'], sort_TM = False, sort_level = True)

cap.snsStackedBarPlot(data = rg_clones_levels_type_df_E10, 
                  x = 'Level', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = 'Level',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['Limb', 'Thoracic', 'Sacral'],
                  output_path = output_path,
                  palette = ColPals.stckpal_4col,
                  fig_size = (3.5,5),
                  hue_order = ['Symmetric', 'Asymmetric', 'Small neurogenic', 'Terminal neurogenic'])

# E11.5
output_path = levels_output_path + '/Clone_type_per_level_E11' + curr_date
rg_clones_levels_type_df_E11 = rg_clones_levels_df_no_unicolor[rg_clones_levels_df_no_unicolor['TM'] == 'E11.5']
rg_clones_levels_type_df_E11 = caf.cloneTypeProportions(rg_clones_levels_type_df_E11, ['Level'], ['Type'], sort_TM = False, sort_level = True)

cap.snsStackedBarPlot(data = rg_clones_levels_type_df_E11, 
                  x = 'Level', 
                  hue = 'Type',
                  weights = 'Proportions',
                  x_axis_title = 'Level',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['Limb', 'Thoracic', 'Sacral'],
                  output_path = output_path,
                  palette = ColPals.stckpal_4col,
                  fig_size = (3.5,5),
                  hue_order = ['Symmetric', 'Asymmetric', 'Small neurogenic', 'Terminal neurogenic'])

# Clone laterality per level
output_path = levels_output_path + '/Clone_laterality_per_level' + curr_date
rg_clones_levels_laterality_df = rg_clones_levels_df[
    rg_clones_levels_df['Laterality'].notna() &
    rg_clones_levels_df['Laterality'].str.strip().ne('')] # Only consider clones with laterality (exclude midline cells-only clones)
rg_clones_levels_laterality_prop_df = caf.cloneTypeProportions(rg_clones_levels_laterality_df, ['Level'], ['Laterality'], sort_TM = False, sort_level = True)

cap.snsStackedBarPlot(data = rg_clones_levels_laterality_prop_df, 
                  x = 'Level', 
                  hue = 'Laterality',
                  weights = 'Proportions',
                  x_axis_title = 'Level',
                  y_axis_title = 'Proportion',
                  x_axis_ticks = ['limb', 'thoracic', 'sacral'],
                  output_path = output_path,
                  palette = ColPals.stckpal_lat,
                  fig_size = (3.5,5),
                  hue_order = ['bi', 'uni'])

saveDataTableProps(df1 = rg_clones_levels_laterality_prop_df,
                   df2 = rg_clones_levels_laterality_df, 
                   df2_columns = ['Clone', 'Level', 'Laterality'], 
                   group_cols = ['Level' ,'Laterality'], 
                   output_path = output_path,
                   ind_var = 'TM', 
                   order = ['uni', 'bi'])
