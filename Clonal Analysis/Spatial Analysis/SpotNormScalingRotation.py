# This script normalizes and scales cell coordinates outputed from Imaris, as well as rotating coordinates
# for sections without perfectly vertical DV orientation. It should be used in conjunction with the Imaris MATLAB script QuadColocSpotsWithTransformCoordinates
# The QuadColocSpots script does not normalize or scale coordinates. It only finds the coloc and outputs
# csv files in the oridinal format. Importantly, it also takes an additional reference point, 'to' which should be placed on
# the uppermost dorsal edge of the SC as the last Spots object. It also asks the user whether the section is rotated, and if so,
# by how many degrees.
#
# Cell coordinates are scaled so they fit into a reference spinal hemisection, like in CompareDistributions.py.
# The reference spinal cord has dimensions 650 (width: CC X to side X); -400 (vheight: CC Y to bottom Y); 500 
# (dheight: CC Y to top Y). The following scaling factors are calculated:
# 
# tScale = 500 / (y coordinate of 'to'  - y coordinate of 'ce') 
# bScale = -400 / (y coordinate of 'bo'  - y coordinate of 'ce')
# lScale = 650 / abs(x coordinate of 'si' - x coordinate of 'ce')
#
# tScale is applied to cells lying above the CC. vScale is applied to cells lying below the CC.
#
# Input: a parent folder that contains subfolders for individual sections. Each subfolder should contain only the Imaris csv 
# output files for one image. It can also contain the .ims image. 
#
# Output: a folder will be created in the location as the input folder with suffix 'NormScaleRotate'. 

import os
from pathlib import Path
import datetime
import pandas as pd
import numpy as np
import copy
import matplotlib.pyplot as plt

###################################################################################################################################
# USER INPUT    
input_folder = r'Z:\People\Sophie\10 AB Spatial Distribution Map\20251120 Iterative Clone Spatial Analysis\Input for SpotsNormScalingRotation'

# Size of graphing SC (from make_contours_mariano.m script)
width = 650 # CC to outermost edge of SC
vheight = -400 # CC to lowest edge of SC
dheight = 500 # CC to highest edge of SC

x_column = 'X coordinate (scaled to Reference SC)' # Name of column containing x coordinates
y_column = 'Y coordinate (scaled to Reference SC)' # Name of column containing y coordinates
###################################################################################################################################
# DEFINE FUNCTIONS

# Rotate x,y coordinates in columns of a data frame named x_column, y_column by specified degrees ("rotation" argument)
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
    coordinates = points[[x_column, y_column]].values

    # Apply the rotation matrix
    rotated_coordinates = coordinates @ rotation_matrix.T

    # Update a new DataFrame with the new rotated coordinates
    rotated_points = copy.deepcopy(points)
    rotated_points[x_column], rotated_points[y_column] = rotated_coordinates.T
    
    return rotated_points
###################################################################################################################################

now = datetime.datetime.now()
curr_date = now.strftime("%Y%m%d_%H%M%S")

# Create an output folder
output_folder = input_folder + ' NormScaleRotate ' + curr_date
os.makedirs(output_folder) 

# Iterate through the folders in the input folder
parent_folder = Path(input_folder)
for folder in parent_folder.iterdir():
    # Create output folders for each of them
    curr_output_folder = output_folder + '/' + folder.name
    os.makedirs(curr_output_folder) 

    # Find all csv files in the parent folder
    csv_files = list(folder.glob("*.csv"))
    print('csv_files', csv_files)
    if not csv_files:
        print(f"No CSV files found in {folder}")
        continue
    
    # Iterate through the other files in the folder
    for f in csv_files:
        print('f', f)
        try:
            df = pd.read_csv(f)
            df.columns = [c.strip() for c in df.columns]
            print(f"\nLoaded file: {f.name} ({len(df)} rows)")
            
            if df.shape[0] == 1:
                continue # skip the 1-row coloc_spots_data file
            
            # get last 5 rows, which contain the reference points and rotation angle.
            tail5 = df.tail(5).reset_index(drop=True)
            rotation = tail5.loc[4, x_column]
            
            # In Imaris, y axis values increase downward. To flip them so they increase upward like on a regular cartesian plane, we need to invert all the y values
            df[y_column] *= -1
            
            # fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=False, sharey=False)
            # df_subset = df.iloc[:-6]
            # axes[0].scatter(df_subset[x_column], df_subset[y_column], s=10, alpha=0.6)
            # axes[0].set_title("Before rotation")
            
            
            # Rotate all coordinates, including reference points, before calculating the scaling factors
            df = rotatePoints(df, rotation)
            
            # df_subset = df.iloc[:-6]
            # axes[1].scatter(df_subset[x_column], df_subset[y_column], s=10, alpha=0.6)
            # axes[1].set_title("After rotation")
            
            
            # get last 5 rows again after rotation, which contain the reference points and rotation angle.
            tail5 = df.tail(5).reset_index(drop=True)
            xOffset = tail5.loc[0, x_column]
            yOffset = tail5.loc[0, y_column]
            bpointx = tail5.loc[1, x_column]
            bpointy = tail5.loc[1, y_column]
            lpointx = tail5.loc[2, x_column]
            lpointy = tail5.loc[2, y_column]
            tpointx = tail5.loc[3, x_column]
            tpointy = tail5.loc[3, y_column]
            
            print(f"  tpointx : {tpointx}")
            print(f"  tpointy : {tpointy}")
            print(f"  xOffset : {xOffset}")
            print(f"  yOffset : {yOffset}")

            # Calculate the scaling factors for x (l), ventral y (b) and dorsal y (t)
            tScale = dheight / (tpointy - yOffset)
            bScale = vheight / (bpointy - yOffset)
            lScale = width / abs(lpointx - xOffset)
            
            # Normalization: subtract the x yOffset from all cells to set the CC to 0,0
            df[x_column] = df[x_column] - xOffset
            df[y_column] = df[y_column] - yOffset
            
            # Scaling: Multiply Y coordinates < 0 by bScale and > 0 by tScale and X coordinates by lScale
            df = df.copy()
            df[x_column] = abs(df[x_column] * lScale) # All cells need to be on the right side of the CC (positive) for plotting on the reference SC. 
                                                      # If the DV axis was not perfectly straightened and some cells fall on the left, they'll be reflected onto the right.
            
            df[y_column] = df[y_column].apply( 
                lambda x: x * tScale if x > 0 else x * bScale # If the cell is above the CC (0), multiply by the dorsal y scaling factor. If it is below, multiply by the ventral factor
                )
            
            # df_subset = df.iloc[:-6]
            # axes[2].scatter(df_subset[x_column], df_subset[y_column], s=10, alpha=0.6)
            # axes[2].set_title("After normalization")
            
            # plt.tight_layout()
            # plt.show()
        except Exception as e:
            print(f"Error reading {f}: {e}")
            continue

        # Save new csv file
        output_file = curr_output_folder + '/' + f.stem + '.csv'

        df.to_csv(output_file, mode = 'w', index = False)
