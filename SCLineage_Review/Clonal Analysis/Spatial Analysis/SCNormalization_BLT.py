### Uses CC plus three reference points (bottom edge (b), left/side edge(l), and top edge (t)) to scale coordinates of cells within a transverse SC section to a reference SC with specified dimensions. 

import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import distance_matrix
from scipy import stats
import numpy as np
import copy
import matplotlib.ticker as ticker
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import datetime
from statistics import mean
import random
from pypalettes import load_cmap
import glob

########################################################################################################################################################################################################
### User input

# Size of graphing SC (from make_contours_mariano.m script)
width = 650 # CC to outermost edge of SC
vheight = -400 # CC to lowest edge of SC
dheight = 500 # CC to highest edge of SC. Only works if it's negative

input_cells_file_path = r'Z:\People\Sophie\9 Mouse SC lineage paper\20251124 Data Analysis and Plots\Code Versions\ClonalAnalysis\Output\20251126_114705_125-200um_exclusion_all_complete_animals/Overall_rg_cells_20251126_114705.csv' # Path of overall_rg_cells csv file
# Input folder path for reference coordinates. All individual reference files should be saved in this folder. 
input_ref_folder_path = r'Z:\People\Sophie\10 AB Spatial Distribution Map\20251128 Clone Spatial Analysis\Code Versions\SCNormalization\Input Files for script\Ref Files'
output_folder_name = '_all_complete_animals' # Desired output folder name
overall_output_file_path = r'Z:\People\Sophie\9 Mouse SC lineage paper\Initial Submission\Nature Submission\GitHub Repository\Output/' # Desired output path

# Find current date and time, which will be used in file naming
now = datetime.datetime.now()
curr_date = now.strftime("%Y%m%d_%H%M%S")

# Output file path 
output_path = overall_output_file_path + curr_date
os.makedirs(output_path)

# Do you want the output clone csvs to contain only neurons (no glia or uncertain cells) or only glia? Set one to True or leave them both False.
neurons_only = True
glia_only = False

# Do you want the output clone csvs to contain only clones with complete cell-type information (i.e. no uncertain cells)
# or clones with 20% or fewer uncertain cells?
complete_only = False
complete_20 = True

max_y = 2239 # Max y value in the image (see image properties)

keep_z = False # Set to true if you want a column with the z coordinates in the output csv files

########################################################################################################################################################################################################    
### Define functions

# Read CSV using pandas
def read_csv(file_path, skip_rows):
    if file_path == 'None':
        points = None
    else: 
        points = pd.read_csv(file_path, skiprows=skip_rows) # skip initial rows
        points.dropna(axis = 0, how = 'all', inplace = True) # drop empty rows
    return points

# Rotate x,y coordinates in columns of a data frame named 'Position X', 'Position Y' by specified degrees ("rotation" argument)
def rotatePoints(points, rotation):
    # Rotation
    angle_radians = np.deg2rad(rotation)
    cos_angle = np.cos(angle_radians)
    sin_angle = np.sin(angle_radians)
    
    # Define the rotation matrix for the given axis
    rotation_matrix = np.array([
        [cos_angle, -sin_angle],
        [sin_angle,  cos_angle]
    ])
    
    # Extract the coordinates as a numpy array
    coordinates = points[['Position X', 'Position Y']].values

    # Apply the rotation matrix
    rotated_coordinates = coordinates @ rotation_matrix.T

    # Update a new DataFrame with the new rotated coordinates
    rotated_points = copy.deepcopy(points)
    rotated_points['Position X'], rotated_points['Position Y'] = rotated_coordinates.T
    
    return rotated_points

# Find the middle point between two points in 1D
def findMid(rostral_point, caudal_point):
    if rostral_point > caudal_point:
        mid = caudal_point + ((rostral_point - caudal_point) / 2)
    else:
        mid = rostral_point + ((caudal_point - rostral_point) / 2)
    return mid

#######################################################################################################################################################
if __name__ == "__main__":
    ### Import an overall_rg_cells file that contains all cells with clone annotations into a DataFrame
    input_cells = read_csv(input_cells_file_path, skip_rows = 0)
    input_cells = input_cells[~input_cells['Clone'].str.contains('exclude', na=False)] # Remove clones that are excluded by ClonalAnalysis.py based on proximity, damage, etc.
    
    ### Import all the Ref files and store them in a big DataFrame
    csv_files = glob.glob(os.path.join(input_ref_folder_path, "*.csv"))
    input_ref_file_list = [pd.read_csv(f) for f in csv_files]
    input_ref = pd.concat(input_ref_file_list, ignore_index=True)
    
    ### Rotate the reference points so that the DV axis is oriented vertically. The cells have already been rotated in ClonalAnalysis.py
    grouped_ref = input_ref.groupby(['Animal', 'Segment'])
    for (animal, segment), group in grouped_ref:
        if group['Rotation'].isna().all(): # If there is no data for this segment, skip it
            continue
        
        mask = (input_ref['Animal'] == animal) & (input_ref['Segment'] == segment)
        subset = input_ref.loc[mask].copy()# make a copy of the subset to rotate  
        
        ### In Imaris, y axis values increase downward. To flip them so they increase upward like on a regular cartesian plane, we need to invert all the y values
        subset['Position Y'] = (subset['Position Y'] - max_y).abs()
        
        # Get the rotation angle
        rotation = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == segment) & (input_ref['Rotation'].notnull()), 'Rotation'].values[0]
        
        # Rotate the reference points by a specified amount to orient DV axis. 
        subset_rotated = rotatePoints(subset, rotation) # rotate the subset
        subset_rotated.index = subset.index # ensure the index matches
        input_ref.loc[mask] = subset_rotated # assign back
    

    ### Make sure all cells are in the right hemisphere (have an x greater than 0) by taking the absolute value of x.
    input_cells['Position X'] = abs(input_cells['Position X'])
        
    ### Scale each clone based on its reference coordinates and the size of the graphing SC.
    # Iterate through each clone in the overall_rg_cells DataFrame and find the corresponding reference coordinates in the Ref DataFrame.
    normalized_clones = []
    grouped_clones = input_cells.groupby(['Animal', 'Clone'])

    
    for (animal, clone), group in grouped_clones:
        current_segment = group.iloc[0]['Segment']
        ### WARNING: Mid vs RC ends 
        if complete_only:
            if (group['Uncertain'] == 'yes').any():
                continue
        if complete_20:
            uncert_mask  = group['Uncertain'].astype(str).str.strip().str.lower().eq('yes')
            if (uncert_mask.sum() / len(group)) > 0.2:   # if the number of uncertain cells in the clone is greater than 20%, skip the clone
                continue
        if neurons_only:
            group = group[group['Neuron'] == 'yes']
            if group.empty:
                continue
        if glia_only:
            group = group[group['Glial cell'] != 'no']
            if group.empty:
                continue
        # For each clone, first check if the clone has its own precise set of reference points (saved in Clone_Ref csv file for clones for whom the average between the rostral and caudal ends wasn't good enough)
        if ((input_ref["Animal"] == animal) & (input_ref["Clone"] == clone)).any():
            print("specific ref points for: ", animal, clone, current_segment)
            rostral_xOffset = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'CC') & (input_ref['Location'] == clone), 'Position X'].iloc[0] # rostral CC x coordinate
            caudal_xOffset = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'CC') & (input_ref['Location'] == clone), 'Position X'].iloc[0] # caudal CC x coordinate
            rostral_yOffset = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'CC') & (input_ref['Location'] == clone), 'Position Y'].iloc[0] # rostral CC y coordinate
            caudal_yOffset = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'CC') & (input_ref['Location'] == clone), 'Position Y'].iloc[0] # caudal CC y coordinate
            rostral_lPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'L') & (input_ref['Location'] == clone), 'Position X'].iloc[0] # rostral L (side) x coordinate
            caudal_lPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'L') & (input_ref['Location'] == clone), 'Position X'].iloc[0] # caudal L (side) x coordinate
            rostral_bPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'B') & (input_ref['Location'] == clone), 'Position Y'].iloc[0] # rostral B (bottom) y coordinate
            caudal_bPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'B') & (input_ref['Location'] == clone), 'Position Y'].iloc[0]# caudal B (bottom) y coordinate
            rostral_tPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'T') & (input_ref['Location'] == clone), 'Position Y'].iloc[0] # rostral T (top) y coordinate
            caudal_tPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Clone'] == clone) & (input_ref['Ref'] == 'T') & (input_ref['Location'] == clone), 'Position Y'].iloc[0]# caudal T (top) y coordinate

               
        # Else, find the XY positions of the rostral and caudal CC, L, and B reference points for the corresponding animal and SC segment.  
        else: 
            rostral_xOffset = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'CC') & (input_ref['Location'] == 'rostral'), 'Position X'].iloc[0] # rostral CC x coordinate
            caudal_xOffset = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'CC') & (input_ref['Location'] == 'caudal'), 'Position X'].iloc[0] # caudal CC x coordinate
            rostral_yOffset = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'CC') & (input_ref['Location'] == 'rostral'), 'Position Y'].iloc[0] # rostral CC y coordinate
            caudal_yOffset = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'CC') & (input_ref['Location'] == 'caudal'), 'Position Y'].iloc[0] # caudal CC y coordinate
            rostral_lPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'L') & (input_ref['Location'] == 'rostral'), 'Position X'].iloc[0] # rostral L (side) x coordinate
            caudal_lPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'L') & (input_ref['Location'] == 'caudal'), 'Position X'].iloc[0] # caudal L (side) x coordinate
            rostral_bPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'B') & (input_ref['Location'] == 'rostral'), 'Position Y'].iloc[0] # rostral B (bottom) y coordinate
            caudal_bPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'B') & (input_ref['Location'] == 'caudal'), 'Position Y'].iloc[0]# caudal B (bottom) y coordinate
            rostral_tPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'T') & (input_ref['Location'] == 'rostral'), 'Position Y'].iloc[0] # rostral T (top) y coordinate
            caudal_tPoint = input_ref.loc[(input_ref['Animal'] == animal) & (input_ref['Segment'] == current_segment) & (input_ref['Ref'] == 'T') & (input_ref['Location'] == 'caudal'), 'Position Y'].iloc[0]# caudal T (top) y coordinate

        ### WARNING: It would be more precise if you got the coordinates of the point between the rostral and caudal ends that's at the z level of the current clone, but this is challenging because after ClonalAnalysis.py the clone is no longer on the same z axis as the reference points 
        # Find the middle coordinates between the rostral and caudal points
        xOffset = findMid(rostral_xOffset, caudal_xOffset)
        yOffset = findMid(rostral_yOffset, caudal_yOffset)
        lPoint = findMid(rostral_lPoint, caudal_lPoint)
        bPoint = findMid(rostral_bPoint, caudal_bPoint)
        tPoint = findMid(rostral_tPoint, caudal_tPoint)
        
        # Calculate the scaling factors for x (l), ventral y (b) and dorsal y (t)
        lScale = width / abs(lPoint - xOffset)
        bScale = vheight / (bPoint - yOffset)
        tScale = dheight / (tPoint - yOffset)

        ### Scale the cell XY positions. The cells are already normalized to the CC midpoint as 0,0 by ClonalAnalysis.py
        group = group.copy()
        group['Position X'] = (group['Position X']) * lScale
        
        group['Position Y'] = group['Position Y'].apply( 
            lambda x: x * tScale if x > 0 else x * bScale # If the cell is above the CC (0), multiply by the dorsal y scaling factor. If it is below, multiply by the ventral factor
            )
        
        normalized_clones.append(group)
        
        ### For each clone, save a csv file with the normalized X and Y coordinates and the necessary placeholder values in the bottom 4 rows. Be careful that the format matches what is needed as input by the make_contours_mariano.m script
        if keep_z:
            output = group[['Position X', 'Position Y', 'Position Z']]
            output.columns = ['X coordinate (scaled to Reference SC)', 'Y coordinate (scaled to Reference SC)', 'Position Z']
            cols = ['X coordinate (scaled to Reference SC)', 'Y coordinate (scaled to Reference SC)', 'Position Z']
        else:
            output = group[['Position X', 'Position Y']]
            output.columns = ['X coordinate (scaled to Reference SC)', 'Y coordinate (scaled to Reference SC)']
            cols = ['X coordinate (scaled to Reference SC)', 'Y coordinate (scaled to Reference SC)']
            
        file_name = output_path + '/' + animal + '_' + clone + '_' + curr_date + '.csv'
        output.to_csv(file_name, mode = 'w', index = False)
   
    ### Plots for troubleshooting     
    # Concatenate all the normalized clones back together
    normalized_clones_df = pd.concat(normalized_clones, ignore_index=True)

    # Group original, rotated, and normalized clones
    grouped_original = input_cells.groupby(['Animal', 'Clone'])
    grouped_rotated  = rotatePoints(input_cells, rotation).groupby(['Animal', 'Clone'])
    grouped_normalized = normalized_clones_df.groupby(['Animal', 'Clone'])
    
    print('\normalized_clones_df\n', normalized_clones_df)

    # Plot per clone
    for (animal, clone), group_orig in grouped_original:
        fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=False, sharey=False)
        
        # Match rotated and normalized groups
        group_rot = grouped_rotated.get_group((animal, clone))
        group_norm = grouped_normalized.get_group((animal, clone))
        
        # Plot original
        axes[0].scatter(group_orig['Position X'], group_orig['Position Y'], s=10, alpha=0.6)
        axes[0].set_title("Original")
        
        # Plot rotated
        axes[1].scatter(group_rot['Position X'], group_rot['Position Y'], s=10, alpha=0.6)
        axes[1].set_title("Rotated")
        
        # Plot normalized
        axes[2].scatter(group_norm['Position X'], group_norm['Position Y'], s=10, alpha=0.6)
        axes[2].set_title("Normalized")
        
        # Add overall label
        fig.suptitle(f"Clone: {animal}, {clone}", fontsize=14)
        plt.tight_layout()
        plt.show()
        #plt.close()





