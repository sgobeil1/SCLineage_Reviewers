import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
import numpy as np
import copy
import matplotlib.ticker as ticker
import seaborn as sns
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as mpatches  # For creating custom legend entries
from ClonalAnalysisFunctions_20260126 import *
import ColPals
import random
from statistics import mean
import os

# ### Set some global plot parameters
# plt.rcParams['font.family'] = 'Arial' # Set font family

# plt.rcParams['svg.fonttype'] = 'none' # Prevent conversion of text to paths (so it is editable in Illustator)

# Takes list of values and plots the distribution of those values
def plotDistribution(values, colour, bins):
    fig = plt.figure()
    ax = fig.add_subplot()
    ax.hist(values, color = colour, bins = bins)
    
    #plt.show()

# 3d plot of cells (points) where 3d coordinates are in dataframe columns labeled 'Position X', 'Position Y', and 'Position Z'.
# boundaries are the z positions that demarcate the borders between segments to be shown as vertical planes in the plot. 'colours' should be 'madm' or 'clone_colours'
def plotCells(cells, boundaries, output_path, endz, colours, image_format, red, green, yellow, cloneplots = False):
    if cloneplots == True:
        # Get unique clones
        unique_clones = cells['Clone'].unique() # WARNING: Currently counting excl_message as a unique clone too which is introducing extra spaces in the figure.
        n_clones = len(unique_clones)
    else:
        n_clones = 0
    
    # Create figures for individual clones (fig_clones) and the whole SC (fig_main). The main figure will always be created but only saved
    # if cloneplots == False
    fig_clones = plt.figure(figsize=(30,5))
    fig_main = plt.figure(figsize = (20,10))
    ax_main = fig_main.add_subplot(1, 1, 1, projection='3d')    
    
    # Scale the plot so that the z (horizontal) axis is longer
    x_scale = 1
    y_scale = 1
    z_scale = 5

    scale = np.diag([x_scale, z_scale, y_scale, 1.0])  # Swap z and y
    scale = scale * (1.0 / scale.max())
    scale[3, 3] = 1.0

    def short_proj():
        return np.dot(Axes3D.get_proj(ax_main), scale)

    ax_main.get_proj = short_proj
    
    if colours == 'madm':
        plot_colours = cells['Colour']
        plot_colours.replace({'r': red, 'g': green, 'y': yellow}) # change 'r,g,y' entries to red and green so you can choose the madmred and madmgreen colours
    elif colours == 'clone_colours':
        plot_colours = cells['Clone Colour']
    
    ax_main.scatter(cells['Position X'], cells['Position Z'], cells['Position Y'], color=plot_colours, alpha = 1)

    # Add a line through the central canal
    ax_main.plot([0, 0], [0, endz], [0, 0], linewidth=2, color='k')
    
    # Add labels
    ax_main.set_xlabel('Mediolateral (\u03bcm)', labelpad=10, fontsize = 18, weight = 'bold')
    ax_main.set_ylabel('Rostrocaudal (\u03bcm)', labelpad=50, fontsize = 18, weight = 'bold')  # Exchanged with z because of axis rotation
    ax_main.set_zlabel('Dorsoventral (\u03bcm)', labelpad=10, fontsize = 18, weight = 'bold')  # Exchanged with y because of axis rotation
    
    # Equal aspect ratio: set x and z axis limits to be the same. Not working properly
    min_limit = min(cells['Position X'].min(), cells['Position Y'].min())
    max_limit = max(cells['Position X'].max(), cells['Position Y'].max())

    ax_main.set_xlim([min_limit, max_limit])
    ax_main.set_zlim([min_limit, max_limit])
    
    # Axis parameters
    ax_main.xaxis.set_major_locator(MultipleLocator(200))
    ax_main.zaxis.set_major_locator(MultipleLocator(200))
    ax_main.yaxis.set_major_locator(MultipleLocator(2000))
    ax_main.spines['top'].set_linewidth(1) 
    ax_main.spines['right'].set_linewidth(1) 
    ax_main.spines['left'].set_linewidth(1)
    ax_main.spines['bottom'].set_linewidth(1)
    ax_main.tick_params(width = 2)
                
    # Add planes to demarcate the segments
    x = np.linspace(ax_main.get_xlim()[0], ax_main.get_xlim()[1], 10)
    z = np.linspace(ax_main.get_zlim()[0], ax_main.get_zlim()[1], 10)
    X, Z = np.meshgrid(x, z)
   
    for boundary in boundaries:
        curr_boundary = np.full_like(X, boundary)
        ax_main.plot_surface(X, curr_boundary, Z, color='gray', alpha=0.1, rstride=100, cstride=100)
    
    # Choose the viewing angle
    ax_main.view_init(elev=30, azim=350)
 
    # Plot for each clone
    if cloneplots == True: 
        clone_num = 1 # WARNING - THIS NUMBER DOES NOT CORRESPOND TO THE CLONE NUMBER IN THE SPREADSHEET. IF YOU WANT IT TO CORRESPOND, UNCOMMENT THE FLAGGED LINE BELOW
        for i, clone in enumerate(unique_clones):
            if clone.startswith('c') == False: # If this is not a clone but rather excluded cells, don't plot
                continue
            clone_cells = cells[cells['Clone'] == clone]
            ax_clone = fig_clones.add_subplot(1, n_clones, clone_num, projection='3d')
            
            ax_clone.scatter(clone_cells['Position X'], clone_cells['Position Z'], clone_cells['Position Y'], color = clone_cells['Colour'], alpha = 0.8)
            
            # Add labels
            ax_clone.set_xlabel('Mediolateral (\u03bcm)', labelpad=10, fontsize = 18, weight = 'bold')
            ax_clone.yaxis.set_ticklabels([]) # Hide z axis text
            ax_clone.set_zlabel('Dorsoventral (\u03bcm)', labelpad=10, fontsize = 18, weight = 'bold') 
            #ax_clone.set_title(f'Clone {clone}', color = clone_cells['Clone Colour'].iloc[0]) # UNCOMMENT IF YOU WANT THE CLONE PLOT TITLES TO CORRESPOND TO THE CLONE NUMBERS IN THE SPREADSHEET
            ax_clone.set_title(('Clone ' + str(clone_num)), color = clone_cells['Clone Colour'].iloc[0], fontsize = 15, weight = 'bold')
            
            # Set limits to match the main plot on X and Z, but truncate Y
            ax_clone.set_xlim(ax_main.get_xlim())
            ax_clone.set_ylim(clone_cells['Position Z'].min(), clone_cells['Position Z'].max())
            ax_clone.set_zlim(ax_main.get_zlim()) 
            
            # Axis parameters
            ax_clone.set_xlim([min_limit, max_limit])
            ax_clone.set_zlim([min_limit, max_limit])
            ax_clone.xaxis.set_major_locator(MultipleLocator(200))
            ax_clone.zaxis.set_major_locator(MultipleLocator(200))
            ax_clone.spines['top'].set_linewidth(2) 
            ax_clone.spines['right'].set_linewidth(2) 
            ax_clone.spines['left'].set_linewidth(2)
            ax_clone.spines['bottom'].set_linewidth(2)
            ax_clone.tick_params(width = 6)
            ax_clone.tick_params(axis = 'z', pad = 8)
                
            # Draw a line through the CC
            ax_clone.plot([0,0], [ax_clone.get_ylim()[0], ax_clone.get_ylim()[1]], [0,0], color = 'k', linewidth = 2)
            
            # Choose a transverse SC viewing angle
            ax_clone.view_init(elev=0, azim=269) # One degree off for azim to make sure z axis labels are on the left side
            
            clone_num += 1
    
    try: # Added because sometimes I get a 'filenotfound' error.        
        if cloneplots == True: 
            fig_clones.tight_layout()
            fig_clones.savefig(output_path + '.eps', format = 'eps', dpi=300)
            plt.close(fig_clones)
        else:
            fig_main.tight_layout()
            fig_main.savefig(output_path + '.eps', format = 'eps', dpi=300)
            plt.close(fig_main)
    except FileNotFoundError as e:
        print(f'FileNotFoundError encountered: {e}. Continuing execution.')    

# Input a df for a single clone to plot its cells in 3d
def plotClone(cells):
    fig_main = plt.figure(figsize = (3,3))
    ax_main = fig_main.add_subplot(1, 1, 1, projection='3d')  
    
    ax_main.scatter(cells['Position X'], cells['Position Z'], cells['Position Y'], color = cells['Colour'], alpha = 1) # Switched z and y position

    # Add a line through the central canal
    ax_main.plot([0, 0], [ax_main.get_ylim()[0], ax_main.get_ylim()[1]], [0, 0], linewidth=2, color='k')
    
    # Add labels
    ax_main.set_xlabel('Mediolateral (\u03bcm)', labelpad=5, fontsize = 18, weight = 'bold')
    ax_main.set_zlabel('Dorsoventral (\u03bcm)', labelpad=5, fontsize = 18, weight = 'bold')  # Exchanged with z because of axis rotation
    ax_main.set_ylabel('Rostrocaudal (\u03bcm)', labelpad=5, fontsize = 18, weight = 'bold')  # Exchanged with y because of axis rotation
    
    # Equal aspect ratio: set x and z axis limits to be the same. Not working properly
    min_limit = min(cells['Position X'].min(), cells['Position Y'].min())
    max_limit = max(cells['Position X'].max(), cells['Position Y'].max())

    ax_main.set_xlim([min_limit, max_limit])
    ax_main.set_zlim([min_limit, max_limit])
    
    # Axis parameters
    ax_main.xaxis.set_major_locator(MultipleLocator(200))
    ax_main.zaxis.set_major_locator(MultipleLocator(200))
    ax_main.yaxis.set_major_locator(MultipleLocator(50))
    ax_main.spines['top'].set_linewidth(1.5) 
    ax_main.spines['right'].set_linewidth(1.5) 
    ax_main.spines['left'].set_linewidth(1.5)
    ax_main.spines['bottom'].set_linewidth(1.5)
    ax_main.tick_params(width = 2)     
    
    # Choose the viewing angle
    ax_main.view_init(elev=10, azim=-110)
    
    fig_main.tight_layout()
    fig_main.show()
    
    input('Press enter key')
    
def plotClone2(cells):
    fig_main = plt.figure(figsize = (3,3))
    ax_main = fig_main.add_subplot(1, 1, 1, projection='3d')  
    
    ax_main.scatter(cells['Position X'], cells['Position Z'], cells['Position Y'], color = cells['Colour'], alpha = 1) # Switched z and y position

    # Add a line through the central canal
    ax_main.plot([0, 0], [6350, 6450], [0, 0], linewidth=2, color='k')
    
    # Add labels
    ax_main.set_xlabel('Mediolateral (\u03bcm)', labelpad=5, fontsize = 18, weight = 'bold')
    ax_main.set_zlabel('Dorsoventral (\u03bcm)', labelpad=5, fontsize = 18, weight = 'bold')  # Exchanged with z because of axis rotation
    ax_main.set_ylabel('Rostrocaudal (\u03bcm)', labelpad=5, fontsize = 18, weight = 'bold')  # Exchanged with y because of axis rotation
    
    ax_main.set_xlim([-25, 400])
    ax_main.set_zlim([-25, 400])
    ax_main.set_ylim([6350, 6450])
    
    # Axis parameters
    ax_main.xaxis.set_major_locator(MultipleLocator(200))
    ax_main.zaxis.set_major_locator(MultipleLocator(200))
    ax_main.yaxis.set_major_locator(MultipleLocator(50))
    ax_main.spines['top'].set_linewidth(1.5) 
    ax_main.spines['right'].set_linewidth(1.5) 
    ax_main.spines['left'].set_linewidth(1.5)
    ax_main.spines['bottom'].set_linewidth(1.5)
    ax_main.tick_params(width = 2)     
    
    # Choose the viewing angle
    ax_main.view_init(elev=10, azim=-110)
    
    fig_main.tight_layout()
    fig_main.show()
    
    input('Press enter key')
       
def boxPlot(df, xlabel, ylabel, output_path):
    fig, ax = plt.subplots()

    df.boxplot(ax = ax, grid = False, fontsize = 10, color = 'grey')
    ax.set_yticks(range(0, round(df.max(numeric_only=True).max()))) # Range of 0 to max value in the df as a whole
    ax.set_xlabel(xlabel, fontsize = 18)
    ax.set_ylabel(ylabel, fontsize = 18)
    
    fig.savefig(output_path + '.svg', format = 'svg', dpi=300)
    plt.close(fig)
    fig.tight_layout()
    
def barPlot(x, y, error, x_axis_title, y_axis_title, output_path, clone_counts, animal_counts, fig_size):
    fig = plt.figure(figsize = fig_size)
    ax = fig.add_subplot()
    ax.bar(x, y, yerr = error, color='lightgrey', capsize = 2, error_kw=dict(elinewidth=1))
    ax.set_xticklabels(('E8.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E8.5 Clones'], animal_counts['Num E8.5 Animals']),
                        'E9.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E9.5 Clones'], animal_counts['Num E9.5 Animals']),
                        'E10.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E10.5 Clones'], animal_counts['Num E10.5 Animals']),
                        'E11.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E11.5 Clones'], animal_counts['Num E11.5 Animals'])),
                       fontsize = 16, weight = 'bold')
    ax.set_xlabel(x_axis_title, fontsize = 18, weight = 'bold')
    ax.set_ylabel(y_axis_title, fontsize = 18, weight = 'bold')
    
    ax.spines['top'].set_visible(False) # Remove top of frame around plot
    ax.spines['right'].set_visible(False) # Remove right side of frame around plot
    ax.spines['left'].set_linewidth(1)
    ax.spines['bottom'].set_linewidth(1)
    ax.tick_params(width = 1)
    plt.yticks(weight = 'bold')
    plt.xticks(weight = 'bold')
    plt.tight_layout()
    
    fig.savefig(output_path + '.svg', format = 'svg', dpi=300, bbox_inches = 'tight')
    plt.close(fig)
    
def snsBarPlot(x_column, y_column, df, x_axis_title, y_axis_title, x_axis_ticks, output_path, colour, order = None, max_y = None, fig_size = (6.5,5)):
    fig = plt.figure(figsize = fig_size)
    ax = fig.add_subplot()
    sns.barplot(x = x_column, y = y_column, data = df, order = order, capsize = .05, color = colour, errorbar = 'se', errwidth = 2, errcolor = 'k')
    sns.swarmplot(x = x_column, y = y_column, data = df, order = order, color = 'black', alpha = 1, size = 3)
    sns.despine() # Remove top and right axis spines
    ax.set_xticklabels(x_axis_ticks,
                       fontsize = 16, 
                       weight = 'bold')
    ax.tick_params(axis = 'y', labelsize = 16)  # Adjust font size and weight for y-axis tick labels
    ax.set_xlabel(x_axis_title, fontsize = 18, weight = 'bold')
    ax.set_ylabel(y_axis_title, fontsize = 18, weight = 'bold')
    
    ax.spines['top'].set_visible(False) # Remove top of frame around plot
    ax.spines['right'].set_visible(False) # Remove right side of frame around plot
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.tick_params(width = 1.5)
    plt.yticks(weight = 'bold')
    plt.xticks(weight = 'bold')
    if max_y != None:
        plt.ylim(0, max_y)
    plt.tight_layout()
    fig.savefig(output_path + '.svg', format = 'svg', dpi=300, bbox_inches = 'tight')
    #plt.show()
    plt.close(fig)

# Make a stacked bar plot using seaborn. Input a dataframe and state what the index (x axis categories) should be (column in df)
def snsStackedBarPlot(
    data, x, hue, weights,
    x_axis_title, y_axis_title, x_axis_ticks,
    output_path, palette,
    fig_size=(6.5, 5),
    hue_order=None,
    se_col='SE',                
    err_color='black',
    capsize=3,
    bar_width=0.8
):
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot()

    # Order x categories
    if x_axis_ticks is not None:
        x_order = list(x_axis_ticks)
    else:
        x_order = list(pd.unique(data[x]))

    # Order hue categories
    if hue_order is None:
        hue_levels = list(pd.unique(data[hue]))
    else:
        hue_levels = list(hue_order)

    # Build pivot tables for heights and SEs
    h = (data.pivot_table(index=x, columns=hue, values=weights, aggfunc='sum', fill_value=0)
            .reindex(index=x_order, columns=hue_levels, fill_value=0))

    if se_col in data.columns:
        se = (data.pivot_table(index=x, columns=hue, values=se_col, aggfunc='first')
                .reindex(index=x_order, columns=hue_levels))
    else:
        se = None

    # Palette mapping
    if isinstance(palette, dict):
        color_map = palette
    else:
        color_map = {lvl: palette[i] for i, lvl in enumerate(hue_levels)}

    x_pos = np.arange(len(x_order))
    bottoms = np.zeros(len(x_order), dtype=float)

    handles = []
    labels = []

    # Draw stacked bars + one-sided lower error bars at top of each segment
    for lvl in hue_levels:
        heights = h[lvl].to_numpy(dtype=float)

        bars = ax.bar(
            x_pos, heights, bottom=bottoms,
            width=bar_width, color=color_map.get(lvl, None),
            edgecolor='k',
            alpha = 0.75
        )
        handles.append(bars[0])
        labels.append(lvl)

        if se is not None:
            se_vals = se[lvl].to_numpy(dtype=float)

            # One-sided (lower only) error, clamped so it doesn't go below the segment base
            yerr_lower = np.minimum(np.nan_to_num(se_vals, nan=0.0), heights)
            y_top = bottoms + heights

            ax.errorbar(
                x_pos, y_top,
                yerr=[yerr_lower, np.zeros_like(yerr_lower)],
                fmt='none',
                ecolor=err_color,
                elinewidth=1.5,
                capsize=capsize
            )

        bottoms += heights

    # Formatting
    ax.set_xticks(x_pos)
    ax.set_xticklabels(x_order if x_axis_ticks is None else x_axis_ticks,
                       fontsize=16, weight='bold')
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

    ax.legend(handles, labels, loc="upper left", bbox_to_anchor=(1, 1), frameon=False,
              title=hue)
    plt.setp(ax.get_legend().get_texts(), fontsize='16')
    plt.setp(ax.get_legend().get_title(), fontsize='16')

    fig.savefig(output_path + '.svg', format='svg', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
# Create a scatterplot of two continuous variables. If you want a lineage regression line, set regline = True. 
def snsScatterPlot(data, x, y, x_axis_title, y_axis_title, output_path, palette, regline = False, hue = None, transparent = False, figsize = (10,5)):
    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot()
    
    if transparent is True:
        alpha_val = 0.5
    else:
        alpha_val = 1
        
    sns.scatterplot(data = data, 
                    x = x, 
                    y = y, 
                    hue = hue,
                    palette = palette,
                    alpha=alpha_val)
    if regline == True:
        sns.regplot(data = data, x = x, y = y, scatter = False, color = 'black') # Add a regression line
    ax.set_xlabel(x_axis_title, fontsize = 18, weight = 'bold')
    ax.set_ylabel(y_axis_title, fontsize = 18, weight = 'bold') 
    ax.tick_params(axis = 'y', labelsize = 16)  # Adjust font size and weight for y-axis tick labels
    ax.tick_params(axis = 'x', labelsize = 16) 
    ax.spines['top'].set_visible(False) # Remove top of frame around plot
    ax.spines['right'].set_visible(False) # Remove right side of frame around plot
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.tick_params(width = 1.5)
    plt.yticks(weight = 'bold')
    plt.xticks(weight = 'bold')
    #sns.move_legend(ax, "upper left", bbox_to_anchor=(1, 1)) # Move legend outside of plot
    if hue is not None:
        plt.setp(ax.get_legend().get_texts(), fontsize='16') # Set font size for legend text
        plt.setp(ax.get_legend().get_title(), fontsize='16') # Set font size for legend title
    plt.tight_layout()
    fig.savefig(output_path + '.svg', format = 'svg', dpi=300, bbox_inches = 'tight')
    plt.close(fig)
    #plt.show()
    
# Create a grouped barplot
def snsGroupedBarPlot(data, x, y, hue, x_axis_title, y_axis_title, output_path, palette):
    g = sns.catplot(data = data,
                kind = 'bar', 
                x = x, 
                y = y, 
                hue = hue, 
                errorbar = 'se', 
                capsize = 0.2,
                err_kws = {'color': 'k', 'linewidth': 1},
                palette = palette,
                height = 5, 
                aspect = 1)
    
    ax = g.ax # This gets the Axe object from the FacetGrid which the catplot has created
    ax.set_xlabel(x_axis_title, fontsize = 18, weight = 'bold')
    ax.set_ylabel(y_axis_title, fontsize = 18, weight = 'bold') 
    ax.spines['top'].set_visible(False) # Remove top of frame around plot
    ax.spines['right'].set_visible(False) # Remove right side of frame around plot
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.tick_params(width = 1.5)
    plt.yticks(weight = 'bold')
    plt.xticks(weight = 'bold')
                
    plt.tight_layout()
    fig = g.fig
    fig.savefig(output_path + '.svg', format = 'svg', dpi=300, bbox_inches = 'tight')
    #plt.show()
    
# Create a grouped plot of stacked bars using customized seaborn plots thanks to ChatGPT
def snsGroupedStackedBarPlot(data, x, hue, weights, group, x_axis_title, y_axis_title, output_path):
    """
    Creates a grouped plot of stacked bars using seaborn.
    :param data: DataFrame containing the data
    :param x: Column for x-axis categories
    :param hue: Column to define stack categories
    :param weights: Column for bar heights (e.g., proportions or counts)
    :param group: Column for grouping (each group will have a set of stacked bars at each x-axis category)
    :param x_axis_title: Title for the x-axis
    :param y_axis_title: Title for the y-axis
    :param output_path: File path to save the plot
    """
    
    # Maintain the original order of x (TM) as in the data
    unique_x = pd.Categorical(data[x], categories=data[x].unique(), ordered=True)
    data[x] = unique_x
    unique_x_ordered = unique_x.categories
    
    # Sort the data by TM and group (Level) to ensure proper ordering
    data = data.sort_values(by=[x, group], key=lambda col: col.map({'limb': 0, 'thoracic': 1, 'sacral': 2}))
    
    unique_groups = ['limb', 'thoracic', 'sacral']  # Specify the fixed order of groups
    unique_hue = sorted(data[hue].unique())
    group_count = len(unique_groups)
    bar_width = 0.8 / group_count  # Divide space at each x axis tick for evenly for the number of groups
    colors = ('darksalmon', 'mediumpurple', 'springgreen', 'gold', 'grey')
    
    # Make sure the df contains all possible combinations of x, group, and hue, and fill empty ones with 0. This is necessary for the plotting
    # Create a complete grid of x, group, and hue. 
    complete_grid = pd.DataFrame(
        [(tm, lvl, typ) for tm in unique_x_ordered for lvl in unique_groups for typ in unique_hue],
        columns=[x, group, hue])
    data = complete_grid.merge(data, on=[x, group, hue], how="left").fillna({weights: 0}) # Merge with the original data to fill missing combinations with 0
        
    # Initialize the plot
    fig, ax = plt.subplots(figsize=(10, 5))
    
    # Iterate through each group to plot stacked bars
    # For each group, the function stacks bars based on the hue variable by incrementing the bottom variable.
    for i, grp in enumerate(unique_groups):
        # Filter data for the current group
        group_data = data[data[group] == grp]
        
        # Initialize bottom for stacking
        bottom = np.zeros(len(group_data[x].unique()))
        
        # Plot each stack category
        for j, hue_level in enumerate(data[hue].unique()):
            hue_data = group_data[group_data[hue] == hue_level]
            x_positions = np.arange(len(unique_x_ordered)) + (i - group_count / 2) * bar_width
            
            bars = ax.bar(
                x_positions,
                hue_data[weights],
                bar_width,
                bottom=bottom,
                label=f"{grp}-{hue_level}" if i == 0 else "",
                color=colors[j]
            )
            
            # Add labels to the bottom of the bars (group name: Ce, Th, Lu, Sa)
            for bar, label in zip(bars, hue_data[group].values):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,  # Center of the bar
                    0,  # Bottom of the bar
                    str(label),  # Label text
                    ha='center', va='bottom', fontsize=10, weight="bold"
                )
                
            bottom += hue_data[weights].values
        
    # Customize the plot
    ax.set_xlabel(x_axis_title, fontsize=18, weight="bold")
    ax.set_ylabel(y_axis_title, fontsize=18, weight="bold")
    ax.set_xticks(np.arange(len(data[x].unique())))
    ax.set_xticklabels(data[x].unique(), fontsize=16, weight="bold")
    ax.legend(title=f"{hue} by {group}", bbox_to_anchor=(1.05, 1), loc="upper left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.tick_params(width=1.5)
    plt.yticks(weight="bold")
    plt.xticks(weight="bold")
    plt.tight_layout()
    
    # Save the plot
    fig.savefig(output_path + ".svg", format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
    
# Make a raster plot with all clones in individual rows to show the cell distribution along the Z axis (R-C) axis.
# Mark the cell-cell distances that are greater than dist_threshold with a line.
# Thanks ChatGPT
def snsRasterPlot(df, dist_threshold, x_axis_title, xmin, xmax, interval, output_path, red, green, row):
    # Initialize lists for plotting
    group_labels = []
    positions = []
    colors = []
    avg_nn_distances = []
    data_counts = []
    gap_lines = []
    tm_values = []    
    
    TM_colour = {'E9.5': ColPals.light_orange, 'E10.5': ColPals.med_orange, 'E11.5': ColPals.dark_orange} # colours for rows based on TM time point
        
    if row == 'Clone': # 1 row per clone
        grouped = df.groupby(['Animal', 'Clone'])

        for (animal, clone), group in grouped:
            # Normalize Position Z
            pos_z = group['Position Z'].values
            normalized_z = pos_z - np.mean(pos_z)
                
            # Calculate average nearest neighbor distance
            distances = np.diff(np.sort(pos_z))
            avg_nn_distance = np.mean(distances) if len(distances) > 0 else 0
                
            # Count data points
            data_count = len(pos_z)
                
            # Identify gaps larger than dist_threshold and record lines for large gaps
            lines = [
                ((normalized_z[i] + 10), (normalized_z[i + 1] - 10))
                for i in range(len(normalized_z) - 1)
                if abs(pos_z[i] - pos_z[i + 1]) > dist_threshold
            ]
                
            # Collect data for plotting
            group_labels.append(f"{animal} - {clone}")
            positions.append(normalized_z)
            colors.append(group['Colour'].values)
            avg_nn_distances.append(avg_nn_distance)
            data_counts.append(data_count)
            gap_lines.append(lines)
                
            # Record TM value for the group
            tm_values.append(group['TM'].iloc[0])  # Assume TM is consistent within the group
            
    elif row == 'Animal': # 1 row per 'Animal'
        grouped = df.groupby(['Animal'])

        for animal, group in grouped:
            # Don't normalize Position Z
            pos_z = group['Position Z'].values
            
            # Calculate average nearest neighbor distance
            distances = np.diff(np.sort(pos_z))
            avg_nn_distance = np.mean(distances) if len(distances) > 0 else 0
                
            # Count data points
            data_count = len(pos_z)
                
            # Identify gaps larger than dist_threshold and record lines for large gaps
            lines = [
                ((pos_z[i] + 10), (pos_z[i + 1] - 10))
                for i in range(len(pos_z) - 1)
                if abs(pos_z[i] - pos_z[i + 1]) > dist_threshold
            ]
                
            # Collect data for plotting
            group_labels.append(f"{animal}")
            positions.append(pos_z)
            colors.append(group['Colour'].values)
            avg_nn_distances.append(avg_nn_distance)
            data_counts.append(data_count)
            gap_lines.append(lines)
                
            # Record TM value for the group
            tm_values.append(group['TM'].iloc[0])  # Assume TM is consistent within the group

    # Sort groups by number of data points (fewest to most)
    sorted_indices = np.argsort(data_counts)
    group_labels = [group_labels[i] for i in sorted_indices]
    positions = [positions[i] for i in sorted_indices]
    colors = [colors[i] for i in sorted_indices]
    colors = [[red if x == 'r' else green if x == 'g' else x for x in sublist] for sublist in colors]

    tm_values = [tm_values[i] for i in sorted_indices]

    # Determine dynamic x-axis limits
    all_positions = np.concatenate(positions)
    x_min, x_max = all_positions.min() - 10, all_positions.max() + 10

    ### Create the raster plot
    if row == 'Clone':
        row_height = 0.12
        y_label = 'Clones'
        point_size = 20
    elif row == 'Animal':
        row_height = 0.4
        y_label = 'Animal'
        point_size = 40
    fig, ax = plt.subplots(figsize=(16, len(group_labels) * row_height))

    for i, (label, pos, color, lines, tm_value) in enumerate(zip(group_labels, positions, colors, gap_lines, tm_values)):
        # Add background colour for the group based on TM value
        bg_color = TM_colour.get(tm_value, '#FFFFFF')  # Default to white if TM not in TM_colour
        ax.axhspan(i - 0.44, i + 0.44, color=bg_color, alpha=0.3, zorder=0)
        
        # Add group label outside the plot
        ax.text(x_min - 2, i, label, va='center', ha='right', fontsize=8)
        # Plot points
        ax.scatter(pos, [i] * len(pos), c=color, s=point_size, edgecolor='none', zorder=2, alpha = 0.5)

    # Add a legend for the background colours
    legend_patches = [
        mpatches.Patch(color=color, label=f'TM = {tm}') for tm, color in TM_colour.items()
    ]
    ax.legend(handles=legend_patches, 
              loc='right',
              title='TM', 
              fontsize=15, 
              title_fontsize=15, 
              frameon=False)
   
    # Adjust plot
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(False)
    ax.set_yticks([])
    ax.spines["bottom"].set_linewidth(1.5)
    ax.tick_params(width=1.5)
    plt.xticks(weight="bold")
    ax.set_xlabel(x_axis_title, fontsize=18, weight="bold")
    ax.set_ylabel(y_label, fontsize=18, weight="bold")
    ax.set_xlim(xmin, xmax)
    ax.set_xticks(np.linspace(xmin, xmax, interval))  # Adjust tick range and interval as needed

    # Set axis limits for proper spacing of labels
    ax.spines['left'].set_position(('outward', 100))

    # Tight layout for better appearance
    plt.tight_layout()
    
    # Save the plot
    fig.savefig(output_path + ".svg", format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
    
def snsHistPlot(array, x_axis_title, y_axis_title, output_path, colour, binwidth = None, binrange = None, order = None, max_y = None, fig_size = (6.5,5)):
    fig, ax = plt.subplots(figsize = fig_size)
    sns.histplot(array, color = colour, ax = ax, binwidth = binwidth, binrange = binrange) 
    
    sns.despine() # Remove top and right axis spines
    ax.tick_params(axis = 'x', labelsize = 16)
    ax.tick_params(axis = 'y', labelsize = 16) 
    ax.set_xlabel(x_axis_title, fontsize = 18, weight = 'bold')
    ax.set_ylabel(y_axis_title, fontsize = 18, weight = 'bold')
    
    ax.spines['top'].set_visible(False) # Remove top of frame around plot
    ax.spines['right'].set_visible(False) # Remove right side of frame around plot
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.tick_params(width = 1.5)
    plt.yticks(weight = 'bold')
    plt.xticks(weight = 'bold')
    if max_y != None:
        plt.ylim(0, max_y)
    plt.tight_layout()
    fig.savefig(output_path + '.svg', format = 'svg', dpi=300, bbox_inches = 'tight')
    #plt.show()
    #plt.close(fig)
    
def snsLinePlot(df, x_axis_title, y_axis_title, x_axis_ticks, output_path, palette, max_y = None, fig_size = (6.5,5)):
    fig = plt.figure(figsize = fig_size)
    ax = fig.add_subplot()
    sns.lineplot(data = df, 
                 palette = palette,
                 dashes = False,
                 legend = False)
    sns.despine() # Remove top and right axis spines
    ax.set_xticklabels(x_axis_ticks,
                       fontsize = 16, 
                       weight = 'bold')
    ax.tick_params(axis = 'y', labelsize = 16)  # Adjust font size and weight for y-axis tick labels
    ax.set_xlabel(x_axis_title, fontsize = 18, weight = 'bold')
    ax.set_ylabel(y_axis_title, fontsize = 18, weight = 'bold')
    
    ax.spines['top'].set_visible(False) # Remove top of frame around plot
    ax.spines['right'].set_visible(False) # Remove right side of frame around plot
    ax.spines['left'].set_linewidth(1.5)
    ax.spines['bottom'].set_linewidth(1.5)
    ax.tick_params(width = 1.5)
    plt.yticks(weight = 'bold')
    plt.xticks(weight = 'bold')
    if max_y != None:
        plt.ylim(0, max_y)
    plt.tight_layout()
    fig.savefig(output_path + '.svg', format = 'svg', dpi=300, bbox_inches = 'tight')
    # plt.show()
    plt.close(fig)

# Plot the distribution of nearest-neighbour cell XYZ distances, similar to how it's done in Gao et al., 2014
# Perform the nearest-neighbour analysis on the entire dataset without excluding cells. Look for nearest neighbour in xyz axis. Nearest neighbours are found within each animal, then all the distances are pooled to be plotted together.
# Get rg cells and y cells in separate dataframes
# df is the dataframe containing all cells (R/G/Y) 
def NNDPlot(df, average_length, average_x_range, average_y_range, average_num_cells, output_path, inter_cell_threshold, exclude_single_y = True):
    included_animals = []  # will store (Animal, TM)
    if exclude_single_y == True:
        grouped = df.groupby('Animal')
        
        excl_animals = []
        excl_ids = []
        
        for animal, group in grouped:
            y_z_positions = group[group['Colour'] == 'y']['Position Z'].values # Get the Z positions of all y cells in this animal
            y_cells = group[group['Colour'] == 'y'] # y cells in this animal
            
            for index, row in y_cells.iterrows(): # Iterate through each y cell in this animal
                z = row['Position Z']
                if not any(abs(z - x) < inter_cell_threshold and x != z for x in y_z_positions): # If there is no other y cell within the inter_cell_threshold (excluding itself), add it to lists for exclusion. 
                    excl_animals.append(animal)
                    excl_ids.append(row['ID'])
                    
        # Create a set of (animal, id) pairs to exclude
        exclude_pairs = set(zip(excl_animals, excl_ids))

        # Use a mask to filter
        df = df[~df[['Animal', 'ID']].apply(tuple, axis=1).isin(exclude_pairs)]
            
        
        # df = df.reset_index(drop=True)  # Reset index for easier sequential access
        # rows_to_keep = set()
        # i = 0

        # while i < len(df):
        #     current_row = df.iloc[i]

        #     if current_row['Colour'] == 'y':
        #         if i + 1 < len(df):
        #             next_row = df.iloc[i + 1]
        #             prev_row = df.iloc[i - 1]
        #             # If the distance between a y cell and a neighbouring y cell is less than the inter-cell threshold, keep it. Otherwise, exclude it.
        #             if next_row['Colour'] == 'y':
        #                 z_diff = abs(current_row['Position Z'] - next_row['Position Z'])
        #                 if z_diff < inter_cell_threshold:
        #                     rows_to_keep.add(i)
        #                     i += 1
        #                 elif prev_row['Colour'] == 'y':
        #                     z_diff = abs(current_row['Position Z'] - next_row['Position Z'])
        #                     if z_diff < inter_cell_threshold:
        #                         rows_to_keep.add(i)
        #                         i += 1
        #         # If the conditions fail, don't add current row
        #     else:
        #         rows_to_keep.add(i)  # Keep non-'y' Colour rows as-is

        #     i += 1

        # # Create filtered DataFrame
        # filtered_df = df.loc[list(rows_to_keep)].sort_index().reset_index(drop=True)
        
    animals = df.Animal.unique() # Get a list of the unique animals in the dataframe
    
    fig = plt.figure()
    ax = fig.add_subplot()
    
    # Initialize lists to store cumulative distribution functions (CDFs) for each distance type
    cdf_rg_rg_list = []
    cdf_y_y_list = []
    cdf_rg_y_list = []
    distances_rg_y_list = []
    
    # Define common bin edges for all histograms (adjust the range and bin size as needed)
    bin_edges = np.arange(0, 1000, 25)  # Example bin edges from 0 to 1000 with a step of 25

    # Loop over each animal to calculate distances and CDFs
    for animal in animals:
        subset = df[(df['Animal'] == animal)] # df subset for current animal
        subset_rg = subset[(subset['Colour'] == 'r') | (subset['Colour'] == 'g')]
        subset_y = subset[(subset['Colour'] == 'y')]
        
        if subset_rg.empty == True or subset_y.empty == True:
            continue
        
        else:
            
            # Record this animal as included
            tm = subset['TM'].iloc[0]  # TM is constant per animal
            included_animals.append({'Animal': animal, 'TM': tm})
            
            # Find the nearest neighbour of each rg cell among the rg cells
            distances_rg_rg = nearestNeighbourDistance(subset_rg, subset_rg, same = True)    
            
            # Find the nearest neighbour of each y cell among the y cells 
            distances_y_y = nearestNeighbourDistance(subset_y, subset_y, same = True)   
            
            # Find the nearest neighbour of each rg cell among the y cells
            distances_rg_y = nearestNeighbourDistance(subset_y, subset_rg, same = False)   
            distances_rg_y_list.append(distances_rg_y)
            
            # Calculate the CDFs
            cdf_rg_rg, _ = np.histogram(distances_rg_rg, bins = bin_edges, density=True)
            cdf_rg_rg = np.cumsum(cdf_rg_rg) / np.sum(cdf_rg_rg)
            
            cdf_y_y, _ = np.histogram(distances_y_y, bins = bin_edges, density=True)
            cdf_y_y = np.cumsum(cdf_y_y) / np.sum(cdf_y_y)
            
            cdf_rg_y, _ = np.histogram(distances_rg_y, bins = bin_edges, density=True)
            cdf_rg_y = np.cumsum(cdf_rg_y) / np.sum(cdf_rg_y)
            
            # Store the CDFs
            cdf_rg_rg_list.append(cdf_rg_rg)
            cdf_y_y_list.append(cdf_y_y)
            cdf_rg_y_list.append(cdf_rg_y)
    
    # Convert lists to numpy arrays for easier averaging and error calculation
    cdf_rg_rg_array = np.array(cdf_rg_rg_list)
    cdf_y_y_array = np.array(cdf_y_y_list)
    cdf_rg_y_array = np.array(cdf_rg_y_list)

    # Calculate mean and standard deviation (for error bars) of CDFs
    mean_cdf_rg_rg = np.mean(cdf_rg_rg_array, axis=0)
    sem_cdf_rg_rg = stats.sem(cdf_rg_rg_array, axis=0)

    mean_cdf_y_y = np.nanmean(cdf_y_y_array, axis=0) # Ignores nans when computing the mean
    sem_cdf_y_y = stats.sem(cdf_y_y_array, axis=0, nan_policy = 'omit') # Ignores nans when computing the sem

    mean_cdf_rg_y = np.nanmean(cdf_rg_y_array, axis=0)# Ignores nans when computing the mean
    sem_cdf_rg_y = stats.sem(cdf_rg_y_array, axis=0, nan_policy = 'omit')    # Ignores nans when computing the sem

      
    # As a comparison, plot the nearest-neighbour distances in a random dataset where the average number of cells (all cells, R/G/Y) per animal are randomly 
    # distributed in a 3d volume similar to the average spinal cord (average R-C length and range in transverse (xy) plane where MADM cells are found). 
    # This represents a distribution of labeled cells if individual cells were independent and randomly distributed throughout the SC as opposed to being clustered into clones of different progenitors.
    # Make 100 of these plot lines

    # The z range is 0 to the average spinal cord length
    z_range = (0, average_length)

    for j in range(0,100):
        random_x = []
        random_y = []
        random_z = []
        for i in range (1, average_num_cells):
            random_x.append(random.uniform(*average_x_range))
            random_y.append(random.uniform(*average_y_range))
            random_z.append(random.uniform(*z_range))
        random_cells = pd.DataFrame({'Position X': random_x, 'Position Y' : random_y, 'Position Z' : random_z})
        distances_random = nearestNeighbourDistance(random_cells, random_cells, same = True)
        min_dist = min(min(sublist) for sublist in distances_rg_y_list)
        max_dist = max(max(sublist) for sublist in distances_rg_y_list)
                
        ax.hist(distances_random, color = 'lightgrey', histtype='step', density = True, bins = np.arange(min_dist, max_dist + 25, 25), cumulative = True)
    
    # Plot the average CDFs with error bars (doing it down here so that they end up on top of the random data)
    bin_centers = bin_edges[:-1] + np.diff(bin_edges) / 2
    
    ax.errorbar(bin_centers, mean_cdf_y_y, yerr=sem_cdf_y_y, color='orange', capsize = 2, label='Y-Y')
    ax.errorbar(bin_centers, mean_cdf_rg_rg, yerr=sem_cdf_rg_rg, color='mediumblue', capsize = 2, label='RG-RG')
    ax.errorbar(bin_centers, mean_cdf_rg_y, yerr=sem_cdf_rg_y, color='deeppink', capsize = 2, label='RG-Y')
    
    ax.set_xlabel('Nearest Neighbour Distance (\u03bcm)', labelpad = 10, fontsize = 12)
    ax.set_ylabel('Cumulative Percentage', labelpad = 10, fontsize = 12)
    ax.set_title('Nearest-neighbour distance', fontsize = 12)
    ax.xaxis.set_minor_locator(MultipleLocator(100))
    ax.xaxis.set_major_locator(MultipleLocator(500))
    ax.set_xlim(0, 1000)
    fig.savefig(output_path + '.svg', format = 'svg', dpi = 1200)
    fig.savefig(output_path + '.png', format = 'png', dpi = 1200)
    plt.close(fig)
    
    included_animals_df = pd.DataFrame(included_animals)

    animals_per_tm = (
        included_animals_df
        .drop_duplicates()           # one row per animal
        .groupby('TM')
        .size()
        .reset_index(name='Num animals (NND)'))
        
    output_excel_path = output_path + "_NND_animal_number.xlsx"
    animals_per_tm.to_excel(output_excel_path, index=False)

    
def plotCloneSize(clone_type, data, x_axis_title, y_axis_title, output_path, colour = 'grey', fig_size = (3,5)):
    clones = data[(data['Type'] == clone_type)]
    
    # Get a df with the clone_counts and animal_counts
    expected_TMs = ['E9.5', 'E10.5', 'E11.5'] # List of all expected TM values
    counts_df = (
        clones.groupby("TM")  # Group by TM
        .agg(
            clone_num=("TM", "size"),  # Total entries with this TM value
            animal_num=("Animal", "nunique")  # Unique Animal entries with this TM value
        )
        .reset_index()
    )
    # Add 0s for missing TMs
    default_df = pd.DataFrame({'TM': expected_TMs})
    counts_df = default_df.merge(counts_df, on='TM', how='left').fillna(0)
    
    clone_counts = {'Num E9.5 Clones': counts_df.loc[counts_df['TM'] == 'E9.5', 'clone_num'].values[0],
                    'Num E10.5 Clones': counts_df.loc[counts_df['TM'] == 'E10.5', 'clone_num'].values[0],
                    'Num E11.5 Clones': counts_df.loc[counts_df['TM'] == 'E11.5', 'clone_num'].values[0]}
    animal_counts = {'Num E9.5 Animals': counts_df.loc[counts_df['TM'] == 'E9.5', 'animal_num'].values[0],
                    'Num E10.5 Animals': counts_df.loc[counts_df['TM'] == 'E10.5', 'animal_num'].values[0],
                    'Num E11.5 Animals': counts_df.loc[counts_df['TM'] == 'E11.5', 'animal_num'].values[0]}
    
    x_axis_ticks = ('E9.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E9.5 Clones'], animal_counts['Num E9.5 Animals']),
                    'E10.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E10.5 Clones'], animal_counts['Num E10.5 Animals']),
                    'E11.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E11.5 Clones'], animal_counts['Num E11.5 Animals']))
    
    snsBarPlot(x_column = 'TM', 
               y_column = 'Size', 
               df = clones[clones['TM'].isin(['E9.5', 'E10.5', 'E11.5'])], 
               x_axis_title = x_axis_title,
               y_axis_title = y_axis_title,
               x_axis_ticks = x_axis_ticks,
               output_path = output_path,
               colour = colour,
               fig_size = fig_size)
    
def plot_clone_size_metrics(overall_cells_df, animal_counts, general_output_path, curr_date, suffix=""):
    """
    Generate clone size summary plots and CSV outputs for a given dataset.
    
    Parameters
    ----------
    overall_cells_df : pd.DataFrame
        Dataframe containing all RG cells (possibly excluding single-cell clones).
    animal_counts : dict
        Dictionary with animal counts per TM stage.
    general_output_path : str
        Path where results will be saved.
    curr_date : str
        Date string used in output filenames.
    suffix : str
        Optional suffix for distinguishing datasets (e.g. "_not_excl_single-cell").
    """

    # Compute clone size dataframes per TM stage
    clone_size_df_E8 = cellsPerClone(overall_cells_df, 'E8.5')
    clone_size_df_E9 = cellsPerClone(overall_cells_df, 'E9.5')
    clone_size_df_E10 = cellsPerClone(overall_cells_df, 'E10.5')
    clone_size_df_E11 = cellsPerClone(overall_cells_df, 'E11.5')

    # Normalize lengths to handle NaN padding
    max_length = max(len(clone_size_df_E8), len(clone_size_df_E9), len(clone_size_df_E10), len(clone_size_df_E11))
    def pad(values): return values + [np.nan] * (max_length - len(values))
    
    cell_num_df = pd.DataFrame({
        'E8.5': pad(clone_size_df_E8['Total cells'].tolist()),
        'E9.5': pad(clone_size_df_E9['Total cells'].tolist()),
        'E10.5': pad(clone_size_df_E10['Total cells'].tolist()),
        'E11.5': pad(clone_size_df_E11['Total cells'].tolist())
    })

    # Melt into long format for seaborn
    cell_num_df_long = cell_num_df.melt(ignore_index=False, var_name='TM', value_name='Clone counts').reset_index()

    # Custom x-axis tick labels
    x_axis_ticks = (
        f"E9.5\n$\it{{n={cell_num_df['E9.5'].count()}}}$\n$\it{{a={animal_counts['Num E9.5 Animals']}}}$",
        f"E10.5\n$\it{{n={cell_num_df['E10.5'].count()}}}$\n$\it{{a={animal_counts['Num E10.5 Animals']}}}$",
        f"E11.5\n$\it{{n={cell_num_df['E11.5'].count()}}}$\n$\it{{a={animal_counts['Num E11.5 Animals']}}}$"
    )

    # Ensure output folder exists
    os.makedirs(general_output_path, exist_ok=True)

    # Bar plot: clone size per TM
    output_path = f"{general_output_path}/Cells_per_Clone{suffix}_{curr_date}"
    snsBarPlot(
        x_column='TM',
        y_column='Clone counts',
        df=cell_num_df_long[cell_num_df_long['TM'].isin(['E9.5', 'E10.5', 'E11.5'])],
        x_axis_title='TM',
        y_axis_title='Clone size',
        x_axis_ticks=x_axis_ticks,
        output_path=output_path,
        colour=ColPals.lightblue,
        fig_size=(3, 5)
    )

    # Save data to CSV
    cell_num_df.to_csv(f"{general_output_path}/Cells_per_Clone{suffix}_{curr_date}.csv", index=False)

    # Grouped bar plot: clone size per TM and level
    clone_size_df = pd.concat([clone_size_df_E8, clone_size_df_E9, clone_size_df_E10, clone_size_df_E11])
    n = clone_size_df.groupby(['TM', 'Level']).size().reset_index(name='count')
    print('\n n for clone size per level stacked bar plot: \n', n)

    snsGroupedBarPlot(
        data=clone_size_df,
        x='TM',
        y='Total cells',
        hue='Level',
        x_axis_title='TM',
        y_axis_title='Clone size',
        output_path=f"{general_output_path}/Cells_per_Clone_per_Level{suffix}_{curr_date}",
        palette=ColPals.TMpal
    )

    clone_size_df.to_csv(f"{general_output_path}/Cells_per_Clone_per_Level_and_TM{suffix}_{curr_date}.csv", index=False)

    # Bar plot: clone size per TM amount
    tm_amount_order = [0.25, 0.5, 1, 2]
    clone_size_df['TM amount'] = pd.Categorical(clone_size_df['TM amount'], categories=tm_amount_order, ordered=True)
    clone_size_df = clone_size_df.sort_values('TM amount').reset_index(drop=True)
    n = clone_size_df.groupby(['TM amount']).size().reset_index(name='count')
    print('\n n for clone size per TM amount bar plot: \n', n)

    snsBarPlot(
        x_column='TM amount',
        y_column='Total cells',
        df=clone_size_df[clone_size_df['TM'].isin(['E9.5', 'E10.5', 'E11.5'])],
        x_axis_title='TM amount',
        y_axis_title='Clone size',
        x_axis_ticks=['0.25', '0.5', '1', '2'],
        output_path=f"{general_output_path}/Cells_per_Clone_per_TM_amount{suffix}_{curr_date}",
        colour=ColPals.lightblue,
        fig_size=(3, 5)
    )


def cellTypeHist(df, cols, output_path):
    # Create and save:
    #   1) a vertical bar chart of clone composition frequencies, and
    #   2) a pie chart showing the share of each composition.

    # A "composition" is defined per row as the set of cell types whose corresponding
    # count column is > 0. Compositions are labeled in a fixed order based on the keys
    # of `cols` and rendered as strings like "Neurons + Astrocytes + OL". Plots are
    # sorted by frequency (descending). Bar x-labels are rotated for readability.

    # Parameters
    # ----------
    # df : pandas.DataFrame
    #     DataFrame where each row is a clone and columns contain counts for each cell type.
    # cols : dict[str, str]
    #     Mapping from human-readable type names (e.g., "Neurons") to column names in `df`.
    #     The key order determines the label order when building composition strings.
    # output_path : str or pathlib.Path
    #     Path prefix (without extension). The bar chart is saved to f"{output_path}.svg"
    #     and the pie chart to f"{output_path}_pie.svg".

    # Notes
    # -----
    # - A type is considered present if its count column > 0; otherwise it is omitted.
    # - Rows with no present types receive the label "None".
    # - To include every possible non-empty combination (even those with zero counts),
    #   uncomment the `all_labels` block inside the function.
    # - Requires: pandas, seaborn, matplotlib.
    
    type_order = list(cols.keys())  # fixes the order in the labels

    # --- Build a composition label per row (order is fixed; >0 means “present”) ---
    df = df.copy()  # your dataframe
    def composition_label(row):
        present = [t for t, col in cols.items() if row[col] > 0]
        return " + ".join(present) if present else "None"

    df["Composition"] = df.apply(composition_label, axis=1)

    # --- Count frequencies ---
    comp_counts = (df["Composition"]
                .value_counts()
                .rename_axis("Composition")
                .reset_index(name="Frequency"))

    # # If you really want every possible (non-empty) combination to appear (even 0-count ones):
    # all_labels = [" + ".join(c) for k in range(1, len(type_order)+1)
    #               for c in combinations(type_order, k)]
    # comp_counts = (comp_counts.set_index("Composition")
    #                .reindex(all_labels, fill_value=0)
    #                .reset_index())

    # --- Simple overall histogram (vertical bar chart) ---
    order = comp_counts.sort_values("Frequency", ascending=False)["Composition"].tolist()

    fig, ax = plt.subplots(figsize=(max(3, 0.6*len(order)), 7))
    sns.barplot(
        data=comp_counts,
        x="Composition", y="Frequency",
        order=order, edgecolor="black", orient="v", ax=ax
    )
    
    ax.set_xticks(np.arange(len(order)))
    ax.set_xticklabels(order, rotation=70, ha='right', rotation_mode='anchor')

    ax.set_xlabel("Clone composition (types with >0)")
    ax.set_ylabel("Number of clones")
    
    # Slightly reduce side padding so bars use the width evenly
    ax.margins(x=0.01)
    fig.tight_layout()
    fig.savefig(str(output_path) + ".svg", format="svg", dpi=300)
    plt.close(fig)
    
    # --- Pie chart (composition share) ---
    # Optionally group small slices into "Other" to keep the pie readable
    PIE_TOP_N = 12  # adjust if you want more/fewer labeled slices
    if len(comp_counts) > PIE_TOP_N:
        top = comp_counts.nlargest(PIE_TOP_N, "Frequency")
        other_sum = comp_counts.drop(top.index)["Frequency"].sum()
        pie_labels = top["Composition"].tolist() + (["Other (small slices)"] if other_sum > 0 else [])
        pie_sizes = top["Frequency"].tolist() + ([other_sum] if other_sum > 0 else [])
    else:
        pie_labels = comp_counts["Composition"].tolist()
        pie_sizes = comp_counts["Frequency"].tolist()

    plt.figure(figsize=(6.5, 6.5))
    # Show percentages >= 1% to avoid clutter
    def _pct_fmt(p):
        return f"{p:.1f}%" if p >= 1 else ""
    plt.pie(
        pie_sizes,
        labels=pie_labels,
        autopct=_pct_fmt,
        startangle=90,
        counterclock=False,
        wedgeprops=dict(linewidth=0.5, edgecolor="white"),
    )
    plt.title("Clone composition (share of clones)")
    plt.tight_layout()
    plt.savefig(str(output_path) + "_pie.svg", format="svg", dpi=300, bbox_inches="tight")
    plt.close()
    
def snsViolinPlot(
    x_column, y_column, df,
    x_axis_title, y_axis_title, x_axis_ticks,
    output_path, colour,
    order=None, max_y=None, fig_size=(6.5, 5),
    show_points=True,
    violin_width=0.8,
    mean_line=True,
    mean_line_kwargs=None,
    global_mean_line=False,
    global_mean_kwargs=None,
    sem_bars = True,
    sem_kwargs = None):
    """
    Violin plot of per-category distributions with optional points and mean lines.

    - Draws a Seaborn violin for each category in `x_column` (values from `y_column`).
    - Overlays individual data points for visibility (strip plot).
    - Marks the mean of each category with a short horizontal line centered on the violin.
    - (Optional) Adds a global horizontal mean line across all categories.

    Parameters
    ----------
    x_column : str
        Categorical column for grouping (x-axis).
    y_column : str
        Numeric column to plot (y-axis)
    df : pandas.DataFrame
        Data source.
    x_axis_title, y_axis_title : str
        Axis titles.
    x_axis_ticks : list[str]
        Labels to show under each category (same order as `order`).
    output_path : str or pathlib.Path
        Path prefix (SVG will be saved to f"{output_path}.svg").
    colour : str
        Violin color (e.g., "#4C72B0").
    order : list[str], optional
        Category order along x. If None, inferred from data.
    max_y : float, optional
        Upper y-limit.
    fig_size : tuple, optional
        Figure size in inches.
    show_points : bool, optional
        Overlay individual observations.
    violin_width : float, optional
        Width of each violin (default 0.8).
    mean_line : bool, optional
        Draw a short line at the per-category mean.
    mean_line_kwargs : dict, optional
        Matplotlib line kwargs for per-category mean lines.
        Defaults to {"color":"black","linewidth":2.2,"zorder":6}.
    global_mean_line : bool, optional
        Draw a dashed horizontal line at the overall mean across all data.
    global_mean_kwargs : dict, optional
        Matplotlib line kwargs for the global mean line.
        Defaults to {"color":"black","linestyle":"--","linewidth":1.5","alpha":0.8,"zorder":4}.
    """
    if order is None:
        order = list(pd.unique(df[x_column]))

    grp = df.groupby(x_column, observed=True)[y_column]

    means = grp.mean().reindex(order)
    # SEM = sample std / sqrt(n)
    sems = (grp.std(ddof=1) / np.sqrt(grp.count())).reindex(order)
    
    fig = plt.figure(figsize=fig_size)
    ax = fig.add_subplot()

    # Violin plot (distribution)
    sns.violinplot(
        x=x_column, y=y_column, data=df,
        order=order,
        color=colour, inner="quart",
        cut=0, linewidth=1, width=violin_width,
        fill = True, alpha = 0.8
    )

    # Optional: raw points
    if show_points:
        sns.stripplot(
            x=x_column, y=y_column, data=df,
            order=order,
            color="black", alpha=0.4, size=6, jitter=True, zorder=3
        )

    # Per-category mean line (short segment centered on each violin)
    if mean_line:
        if mean_line_kwargs is None:
            mean_line_kwargs = {"color": "black", "linewidth": 1, "zorder": 6}
        # Length of the mean segment as a fraction of the violin width
        seg_half = (violin_width * 0.45)
        for i, m in enumerate(means.values):
            if pd.isna(m):
                continue
            ax.plot([i - seg_half, i + seg_half], [m, m], **mean_line_kwargs)
      # SEM bars (centered at the mean)
    if sem_bars:
        if sem_kwargs is None:
            sem_kwargs = {
                "color": "black",
                "elinewidth": 1.5,
                "capsize": 3,
                "capthick": 1.5,
                "zorder": 7
            }
        for i, (m, s) in enumerate(zip(means.values, sems.values)):
            if pd.isna(m) or pd.isna(s) or s == 0:
                continue
            ax.errorbar(i, m, yerr=s, fmt="none", **sem_kwargs)
            
    # Optional global mean
    if global_mean_line:
        if global_mean_kwargs is None:
            global_mean_kwargs = {"color": "black", "linestyle": "--", "linewidth": 1.5, "alpha": 0.8, "zorder": 4}
        overall_mean = df[y_column].mean()
        ax.axhline(overall_mean, **global_mean_kwargs)

    # Style
    sns.despine()
    ax.set_xlabel(x_axis_title, fontsize=18, weight="bold")
    ax.set_ylabel(y_axis_title, fontsize=18, weight="bold")
    ax.set_xticklabels(x_axis_ticks, fontsize=16, weight="bold")
    ax.tick_params(axis="y", labelsize=16)

    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_linewidth(1.5)
    ax.spines["bottom"].set_linewidth(1.5)
    ax.tick_params(width=1.5)
    plt.yticks(weight="bold")
    plt.xticks(weight="bold")

    if max_y is not None:
        plt.ylim(0, max_y)

    #plt.tight_layout()
    fig.savefig(str(output_path) + ".svg", format="svg", dpi=300, bbox_inches="tight")
    plt.close(fig)
