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
import ClonalAnalysisFunctions as caf
import ClonalAnalysisPlots as cap
import ClonalAnalysisOther as cao
import ColPals
import os
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import datetime
from statistics import mean
import random
from pypalettes import load_cmap

### User input
input_file_path = r'Z:\People\Sophie\9 Mouse SC lineage paper\Initial Submission\Nature Submission\GitHub Repository\demo\ClonalAnalysis_Input_File.csv' # Path of ClonalAnalysis_Input_File csv file
additional_refs_path = r'Z:\People\Sophie\9 Mouse SC lineage paper\Initial Submission\Nature Submission\GitHub Repository\demo\Clone_Ref.csv'# Path of Clone_Ref csv file
output_folder_name = '_all_complete_animals' # Desired output folder name
overall_output_file_path = r'Z:\People\Sophie\9 Mouse SC lineage paper\Initial Submission\Nature Submission\GitHub Repository\Output/' # Desired output path


### Set some global plot parameters
plt.rcParams['font.family'] = 'Arial'           # Set font family
# plt.rcParams['axes.titlesize'] = 14             # Set title font size
# plt.rcParams['xtick.labelsize'] = 12            # Set font size for x-axis tick labels
# plt.rcParams['ytick.labelsize'] = 12            # Set font size for y-axis tick labels
# plt.rcParams['font.weight'] = 'bold'            # Global font weight (affects all text)
plt.rcParams['svg.fonttype'] = 'none' # Prevent conversion of text to paths (so it is editable in Illustator)

if __name__ == "__main__":
    
    input_file = caf.read_csv(input_file_path, skip_rows = 0)
    additional_refs = caf.read_csv(additional_refs_path, skip_rows = 0)
    # In Imaris, the y axis is inverted such that 0 is at the top and the max value is at the bottom. So, we first want to express the Y Positions 
    # as the distance from the x axis. This should be done once for the additional_refs df and then within the for loop for each individual animal ref file/segment
    max_y = 2239 # Max y value in the image (see image properties)
    additional_refs['Position Y'] = (additional_refs['Position Y'] - max_y).abs()
        
    # Set thresholds for clone parsing
    inter_cell_threshold = 125 # Cells less than inter_cell_threshold apart are considered part of the same clone
    inter_clone_threshold = 200 # Inter-clone distance threshold (min). Clones < inter_clone_threshold apart are excluded
    cell_boundary_threshold = 100 # Cell-segment boundary distance threshold (min)
    exclusion_zone = True # Set to true if you want to include an inter_clone_threshold and have a clone exclusion zone
    normalize_z_length = True # Set to true if you want to normalize the cell z coordinates to a spinal cord of average length. This is useful if there is a lot of variability in SC size
    # WARNING: Only the z coordinates get normalized to max values. The xy coordinates only get normalized to the 0,0 of the central canal.
    
    # Do you want to exclude 1-cell clones from most plots?
    exclude_1_cell = True
       
    # Find current date and time, which will be used in file naming
    now = datetime.datetime.now()
    curr_date = now.strftime("%Y%m%d_%H%M%S")
    
    # Create folder for the output
    if exclusion_zone == True:
        threshold_descrip = '_' + str(inter_cell_threshold) + '-' + str(inter_clone_threshold) + 'um_exclusion'
    else:
        threshold_descrip = '_' + str(inter_cell_threshold) + 'um_threshold'
    output_folder = curr_date + threshold_descrip + output_folder_name
    os.makedirs(overall_output_file_path + output_folder) # Create folder
    
    first_pass = True # First pass of the for loop
    
    num_rows_input = input_file.shape[0] # Get number of rows in the dataframe
    
    length_list = [] # Initialize list to collect the lengths of each SC for finding the average length later
    num_cells_list_all = [] # Initialize list to collect the number of cells per SC for finding the average number later
    num_cells_list_rg = [] # Initialize list to collect the number of cells per SC for finding the average number later
    x_range_list = [] # Initialize list to collect the spread of cells in x axis to find the average later
    y_range_list = []# Initialize list to collect the spread of cells in y axis to find the average later
    no_clones = pd.DataFrame({'E8.5':[], 'E9.5':[], 'E10.5':[], 'E11.5':[], 'E12.5':[], 'E13.5':[]}) # Initialize a dataframe that will be used to record 0s in a dictionary if there's an animal with no clones. This will be important for counting the number of clones per animal later.
    all_exclud_clones_damage = []
    all_exclud_clones_boundary = []
    all_exclud_clones_proximity = []
    
    # For loop that will iterate through each row in the input file (i.e. each animal)
    for row in range(0, num_rows_input):
        print('\n current animal: \n', input_file.iloc[row])
        curr_animal = str(input_file.iloc[row]['Litter']) + '_' + str(input_file.iloc[row]['Animal'])
        
        curr_dict = input_file.iloc[row].to_dict() # Convert the current row to a dictionary
        ################################################################################################################################################
        ### Import the reference spinal cord coordinates and specify the output paths

        # Specify the file paths
        if isinstance(curr_dict['s1'], str): # Check that the variable contains a str (file path) and not a nan
            s1_cells_file_path = curr_dict['s1'].strip("'") # The file paths are getting stored with extra quotation marks 
        else:
            s1_cells_file_path = None
        if isinstance(curr_dict['s2'], str):
            s2_cells_file_path = curr_dict['s2'].strip("'") 
        else:
            s2_cells_file_path = None
        if isinstance(curr_dict['s3'], str):
            s3_cells_file_path = curr_dict['s3'].strip("'") 
        else:
            s3_cells_file_path = None
        if isinstance(curr_dict['s4'], str):
            s4_cells_file_path = curr_dict['s4'].strip("'")
        else:
            s4_cells_file_path = None
        if isinstance(curr_dict['s5'], str):
            s5_cells_file_path = curr_dict['s5'].strip("'")
        else:
            s5_cells_file_path = None 
        if isinstance(curr_dict['s6'], str):
            s6_cells_file_path = curr_dict['s6'].strip("'")
        else:
            s6_cells_file_path = None 
        
        if isinstance(curr_dict['ref'], str):
            a_ref_file_path = curr_dict['ref'].strip("'") 
               
        # Read the input csvs into panda dataframes
        if s1_cells_file_path is not None:
            s1_cells = caf.read_csv(s1_cells_file_path, skip_rows = 0)
        else :
            s1_cells = None
        if s2_cells_file_path is not None:
            s2_cells = caf.read_csv(s2_cells_file_path, skip_rows = 0)
        else :
            s2_cells = None
        if s3_cells_file_path is not None:
            s3_cells = caf.read_csv(s3_cells_file_path, skip_rows = 0)
        else :
            s3_cells = None
        if s4_cells_file_path is not None:
            s4_cells = caf.read_csv(s4_cells_file_path, skip_rows = 0)
        else :
            s4_cells = None
        if s5_cells_file_path is not None:
            s5_cells = caf.read_csv(s5_cells_file_path, skip_rows = 0)
        else :
            s5_cells = None
        if s6_cells_file_path is not None:
            s6_cells = caf.read_csv(s6_cells_file_path, skip_rows = 0)
        else :
            s6_cells = None
        
        if a_ref_file_path is not None:
            a_ref = caf.read_csv(a_ref_file_path, skip_rows = 0)

        ### Run a few checks on the dataframes to make sure there are no missing values or inconsistencies
        allowed_cols = ['Litter', 'TM', 'TM Amount', 'Animal', 'Segment', 'ID', 'Position X', 'Position Y', 'Position Z', 
                        'Clone', 'Colour', 'D/V', 'G/W', 'Hemisphere', 'Damage', 'Neuron', 'MN', 'Commissural neuron',
                        'Midline cell', 'Glial cell', 'Uncertain', 'Chx10', 'ChAT', 'Lhx5', 'Notes']

        
        segment_dfs = {'s1_cells': s1_cells, 's2_cells': s2_cells, 's3_cells': s3_cells,
               's4_cells': s4_cells, 's5_cells': s5_cells, 's6_cells': s6_cells} # Dictionary of dataframes so I can store them after modifying them in the for loop below
        segment_count = 1
        
        for name, df in segment_dfs.items():
            if df is None: # If there is no data for this segment, skip it
                segment_count += 1
                continue
                       
            if not df['Colour'].isin(['r', 'g']).any(): # If there are no r or g cells in this segment, skip it
                segment_count += 1
                continue
            
            # Do all columns have expected names?
            if not set(df.columns).issubset(allowed_cols):
                invalid = set(df.columns) - set(allowed_cols)
                raise ValueError(f"Unexpected columns are present in {curr_dict['Litter']}_{curr_dict['Animal']} segment {segment_count}")            
                    
            df = df.copy() 
            df.fillna('no', inplace=True) # When filling in the Spots Data spreadsheets I sometimes leave cells blank instead of writing in "no" to save time.
                        
            # Create a mask for only g/r cells
            colour_mask = df['Colour'].isin(['g', 'r'])
            
            # If a cell is a commissural neuron or a MN, it should be considered a neuron too
            df.loc[colour_mask & (
                (df['Commissural neuron'].isin(['d', 'v'])) |
                (df['MN'].isin(['yes', 'maybe']))
            ), 'Neuron'] = 'yes' 
            
            if (df.loc[colour_mask, 'Commissural neuron'] == "yes").all(axis = 0).any():
                raise ValueError(f"Commissural neuron contains a 'yes' in {curr_dict['Litter']}_{curr_dict['Animal']} segment {segment_count}.")
            
            if (df.loc[colour_mask, 'Midline cell'] == "yes").all(axis = 0).any():
                raise ValueError(f"Midline cell contains a 'yes' in {curr_dict['Litter']}_{curr_dict['Animal']} segment {segment_count}.")
            
            # If a cell is in the white matter but is currently classified as "Uncertain", it should be re-assigned as a "og" (other glia)
            og_mask = colour_mask & (df['G/W'] == 'w') & (df['Uncertain'] == 'yes')
            df.loc[og_mask, 'Glial cell'] = 'og'
            df.loc[og_mask, 'Uncertain'] = 'no'
            
            # Check for more than one positive cell type among g/r rows
            cols_to_check = ['Neuron', 'Midline cell', 'Glial cell', 'Uncertain']# Cell type columns (including only one instance of a neuron type)
            df[cols_to_check] = df[cols_to_check].astype(str)
            df['Glial cell'] = df['Glial cell'].astype(str)
            yes_count = (df.loc[colour_mask, cols_to_check] == 'yes').sum(axis=1) # Count number of yeses per row
            glial_count = df.loc[colour_mask, 'Glial cell'].isin(['pa', 'fa', 'wf', 'gl', 'og', 'ol', 'pr']).astype(int) # Count glial types per row
            if (yes_count > 1).any():
                raise ValueError(f"More than one 'yes' found in cell type columns in at least one row of {curr_dict['Litter']}_{curr_dict['Animal']} segment {segment_count}.")
            mask = (yes_count == 1) & (glial_count == 1) # Check if any row has exactly 1 'yes' AND exactly 1 glial type
            if mask.any():
                raise ValueError(f"More than one 'yes' found in cell type columns in at least one row of {curr_dict['Litter']}_{curr_dict['Animal']} segment {segment_count}.")
            
            # Check for rows with no positive assignments
            if (df.loc[colour_mask, cols_to_check] == "no").all(axis = 1).any():
                raise ValueError(f"No 'yes' or glial type found in any cell type column in at least one row of {curr_dict['Litter']}_{curr_dict['Animal']} s{segment_count}.")
            
            # Save back the modified dataframe
            segment_dfs[name] = df

            segment_count += 1
        
        # Get the modified dataframes
        def get_segment(df_dict, key): # If there is no dictionary entry, there is no data for that segment so it should be None
            val = df_dict.get(key)
            if val is None or (isinstance(val, pd.DataFrame) and val.empty):
                return None
            return val

        s1_cells = get_segment(segment_dfs, 's1_cells')
        s2_cells = get_segment(segment_dfs, 's2_cells')
        s3_cells = get_segment(segment_dfs, 's3_cells')
        s4_cells = get_segment(segment_dfs, 's4_cells')
        s5_cells = get_segment(segment_dfs, 's5_cells')
        s6_cells = get_segment(segment_dfs, 's6_cells')
            
        ################################################################################################################################################
        ### Normalize cell coordinates

        if s1_cells is not None:
            # Retrieve the rotation value where 'segment' is 's1' and 'rotation' is not null
            s1_rotate = a_ref.loc[(a_ref['Segment'] == 's1') & (a_ref['Rotation'].notnull()), 'Rotation'].values[0]

        if s2_cells is not None:
            s2_rotate = a_ref.loc[(a_ref['Segment'] == 's2') & (a_ref['Rotation'].notnull()), 'Rotation'].values[0]
        
        if s3_cells is not None:
            s3_rotate = a_ref.loc[(a_ref['Segment'] == 's3') & (a_ref['Rotation'].notnull()), 'Rotation'].values[0]
            
        if s4_cells is not None:
            s4_rotate = a_ref.loc[(a_ref['Segment'] == 's4') & (a_ref['Rotation'].notnull()), 'Rotation'].values[0]
            
        if s5_cells is not None:
            s5_rotate = a_ref.loc[(a_ref['Segment'] == 's5') & (a_ref['Rotation'].notnull()), 'Rotation'].values[0]
        
        if s6_cells is not None:
            s6_rotate = a_ref.loc[(a_ref['Segment'] == 's6') & (a_ref['Rotation'].notnull()), 'Rotation'].values[0]
        
        # In Imaris, the y axis is inverted such that 0 is at the top and the max value is at the bottom. So, we first want to express the Y Positions 
        # as the distance from the x axis.
        
        if s1_cells is not None:
            s1_cells['Position Y'] = (s1_cells['Position Y'] - max_y).abs()
            columns = s1_cells.columns # Need to get the dataframe column names from a non-None segment to reset dataframes below
        if s2_cells is not None:
            s2_cells['Position Y'] = (s2_cells['Position Y'] - max_y).abs()
            columns = s2_cells.columns
            
        if s3_cells is not None:
            s3_cells['Position Y'] = (s3_cells['Position Y'] - max_y).abs()
            columns = s3_cells.columns
            
        if s4_cells is not None:
            s4_cells['Position Y'] = (s4_cells['Position Y'] - max_y).abs()
            columns = s4_cells.columns
            
        if s5_cells is not None:
            s5_cells['Position Y'] = (s5_cells['Position Y'] - max_y).abs()
            columns = s5_cells.columns
            
        if s6_cells is not None:
            s6_cells['Position Y'] = (s6_cells['Position Y'] - max_y).abs() 
            columns = s6_cells.columns
            
        if all(segments is None for segments in [s1_cells, s2_cells, s3_cells, s4_cells, s5_cells, s6_cells]): # If all segment files are None this means the animal had no clones
            key = 'E' + str(curr_dict['TM Injection Timepoint E']) # Get the TM injection timepoint of the current animal
            no_clones.loc[len(no_clones), key] = 0 # Add this animal to a dataframe to keep track of animals with no clones by adding a new row with a 0 to the appropriate column
            continue # Skip to the next iteration of the for loop.
            
        a_ref['Position Y'] = (a_ref['Position Y'] - max_y).abs() 
         
            
        # Set the first central canal point to 0,0,0 and normalize all other points to this. Also rotate the points by a specified amount. 
        # Do this first individually for all segments.
        # The midpoint between the start and end CC points is used for the reference x and y coordinates. The start CC z position is used as the z coordinate.
        ### WARNING: Since you set the midpoint between the start and end CC points to 0,0,0, you should not rely on the position relative to 0 to judge D/V or L/R hemisphere. Best to rely on hand annotations for this.
        
        # Reset all dfs
        s1_cells_normalized = pd.DataFrame(columns = columns)
        s2_cells_normalized = pd.DataFrame(columns = columns)
        s3_cells_normalized = pd.DataFrame(columns = columns)
        s4_cells_normalized = pd.DataFrame(columns = columns)
        s5_cells_normalized = pd.DataFrame(columns = columns)
        s6_cells_normalized = pd.DataFrame(columns = columns)
        
        if s1_cells is not None:
            mask = (
                (additional_refs['Animal'] == curr_animal) &
                (additional_refs['Segment'] == 's1') &
                (additional_refs['Ref'] == 'CC')
            )
            if not additional_refs.loc[mask].empty:
                print('used additional refs for ', curr_animal, '_s1')
                x1_startCC = additional_refs.loc[mask, 'Position X'].iloc[0]
                y1_startCC = additional_refs.loc[mask, 'Position Y'].iloc[0]
                x1_endCC   = additional_refs.loc[mask, 'Position X'].iloc[0]
                y1_endCC   = additional_refs.loc[mask, 'Position Y'].iloc[0]
            else:
                x1_startCC = a_ref.loc[(a_ref['Segment']=='s1') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position X'].iloc[0]
                y1_startCC = a_ref.loc[(a_ref['Segment']=='s1') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Y'].iloc[0]
                x1_endCC   = a_ref.loc[(a_ref['Segment']=='s1') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position X'].iloc[0]
                y1_endCC   = a_ref.loc[(a_ref['Segment']=='s1') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Y'].iloc[0]
            
            z1_startCC = a_ref.loc[(a_ref['Segment']=='s1') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Z'].iloc[0]
            z1_endCC   = a_ref.loc[(a_ref['Segment']=='s1') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Z'].iloc[0]

            s1_midpoint = (np.array([x1_startCC, y1_startCC, z1_startCC]) + np.array([x1_endCC, y1_endCC, z1_endCC])) / 2
            s1_cells_normalized = caf.normalizePoints2(s1_cells, s1_midpoint[0], s1_midpoint[1], z1_startCC, s1_rotate)
            s1_cells_normalized['Segment'] = 's1'
            
        if s2_cells is not None:
            mask = (
                (additional_refs['Animal'] == curr_animal) &
                (additional_refs['Segment'] == 's2') &
                (additional_refs['Ref'] == 'CC')
            )
            if not additional_refs.loc[mask].empty:
                print('used additional refs for ', curr_animal, '_s2')
                x2_startCC = additional_refs.loc[mask, 'Position X'].iloc[0]
                y2_startCC = additional_refs.loc[mask, 'Position Y'].iloc[0]
                x2_endCC   = additional_refs.loc[mask, 'Position X'].iloc[0]
                y2_endCC   = additional_refs.loc[mask, 'Position Y'].iloc[0]
            else:
                x2_startCC = a_ref.loc[(a_ref['Segment']=='s2') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position X'].iloc[0]
                y2_startCC = a_ref.loc[(a_ref['Segment']=='s2') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Y'].iloc[0]
                x2_endCC   = a_ref.loc[(a_ref['Segment']=='s2') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position X'].iloc[0]
                y2_endCC   = a_ref.loc[(a_ref['Segment']=='s2') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Y'].iloc[0]
            
            z2_startCC = a_ref.loc[(a_ref['Segment']=='s2') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Z'].iloc[0]    
            z2_endCC   = a_ref.loc[(a_ref['Segment']=='s2') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Z'].iloc[0]

            s2_midpoint = (np.array([x2_startCC, y2_startCC, z2_startCC]) + np.array([x2_endCC, y2_endCC, z2_endCC])) / 2
            s2_cells_normalized = caf.normalizePoints2(s2_cells, s2_midpoint[0], s2_midpoint[1], z2_startCC, s2_rotate)
            s2_cells_normalized['Segment'] = 's2'

        # ---- s3 ----
        if s3_cells is not None:
            mask = (
                (additional_refs['Animal'] == curr_animal) &
                (additional_refs['Segment'] == 's3') &
                (additional_refs['Ref'] == 'CC')
            )
            if not additional_refs.loc[mask].empty:
                print('used additional refs for ', curr_animal, '_s3')
                x3_startCC = additional_refs.loc[mask, 'Position X'].iloc[0]
                y3_startCC = additional_refs.loc[mask, 'Position Y'].iloc[0]
                x3_endCC   = additional_refs.loc[mask, 'Position X'].iloc[0]
                y3_endCC   = additional_refs.loc[mask, 'Position Y'].iloc[0]
            else:
                x3_startCC = a_ref.loc[(a_ref['Segment']=='s3') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position X'].iloc[0]
                y3_startCC = a_ref.loc[(a_ref['Segment']=='s3') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Y'].iloc[0]
                x3_endCC   = a_ref.loc[(a_ref['Segment']=='s3') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position X'].iloc[0]
                y3_endCC   = a_ref.loc[(a_ref['Segment']=='s3') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Y'].iloc[0]
                
            z3_startCC = a_ref.loc[(a_ref['Segment']=='s3') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Z'].iloc[0]    
            z3_endCC   = a_ref.loc[(a_ref['Segment']=='s3') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Z'].iloc[0]

            s3_midpoint = (np.array([x3_startCC, y3_startCC, z3_startCC]) + np.array([x3_endCC, y3_endCC, z3_endCC])) / 2
            s3_cells_normalized = caf.normalizePoints2(s3_cells, s3_midpoint[0], s3_midpoint[1], z3_startCC, s3_rotate)
            s3_cells_normalized['Segment'] = 's3'

        # ---- s4 ----
        if s4_cells is not None:
            print('\n HELLOO\n')
            mask = (
                (additional_refs['Animal'] == curr_animal) &
                (additional_refs['Segment'] == 's4') &
                (additional_refs['Ref'] == 'CC')
            )
            print('\nadditional_refs.loc[mask\n', additional_refs.loc[mask])
            if not additional_refs.loc[mask].empty:
                print('used additional refs for ', curr_animal, '_s4')
                x4_startCC = additional_refs.loc[mask, 'Position X'].iloc[0]
                y4_startCC = additional_refs.loc[mask, 'Position Y'].iloc[0]
                print('\ny4_startCC\n', y4_startCC)
                x4_endCC   = additional_refs.loc[mask, 'Position X'].iloc[0]
                y4_endCC   = additional_refs.loc[mask, 'Position Y'].iloc[0]
            else:
                print('nope')
                x4_startCC = a_ref.loc[(a_ref['Segment']=='s4') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position X'].iloc[0]
                y4_startCC = a_ref.loc[(a_ref['Segment']=='s4') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Y'].iloc[0]
                x4_endCC   = a_ref.loc[(a_ref['Segment']=='s4') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position X'].iloc[0]
                y4_endCC   = a_ref.loc[(a_ref['Segment']=='s4') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Y'].iloc[0]
            
            z4_startCC = a_ref.loc[(a_ref['Segment']=='s4') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Z'].iloc[0]
            z4_endCC   = a_ref.loc[(a_ref['Segment']=='s4') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Z'].iloc[0]

            s4_midpoint = (np.array([x4_startCC, y4_startCC, z4_startCC]) + np.array([x4_endCC, y4_endCC, z4_endCC])) / 2
            s4_cells_normalized = caf.normalizePoints2(s4_cells, s4_midpoint[0], s4_midpoint[1], z4_startCC, s4_rotate)
            s4_cells_normalized['Segment'] = 's4'
            print('\nx4_startCC, y4_startCC, z4_startCC\n', x4_startCC, y4_startCC, z4_startCC)

        # ---- s5 ----
        if s5_cells is not None:
            mask = (
                (additional_refs['Animal'] == curr_animal) &
                (additional_refs['Segment'] == 's5') &
                (additional_refs['Ref'] == 'CC')
            )
            if not additional_refs.loc[mask].empty:
                print('used additional refs for ', curr_animal, '_s5')
                x5_startCC = additional_refs.loc[mask, 'Position X'].iloc[0]
                y5_startCC = additional_refs.loc[mask, 'Position Y'].iloc[0]
                x5_endCC   = additional_refs.loc[mask, 'Position X'].iloc[0]
                y5_endCC   = additional_refs.loc[mask, 'Position Y'].iloc[0]
            else:
                x5_startCC = a_ref.loc[(a_ref['Segment']=='s5') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position X'].iloc[0]
                y5_startCC = a_ref.loc[(a_ref['Segment']=='s5') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Y'].iloc[0]
                x5_endCC   = a_ref.loc[(a_ref['Segment']=='s5') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position X'].iloc[0]
                y5_endCC   = a_ref.loc[(a_ref['Segment']=='s5') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Y'].iloc[0]
            
            z5_startCC = a_ref.loc[(a_ref['Segment']=='s5') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Z'].iloc[0]
            z5_endCC   = a_ref.loc[(a_ref['Segment']=='s5') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Z'].iloc[0]

            s5_midpoint = (np.array([x5_startCC, y5_startCC, z5_startCC]) + np.array([x5_endCC, y5_endCC, z5_endCC])) / 2
            s5_cells_normalized = caf.normalizePoints2(s5_cells, s5_midpoint[0], s5_midpoint[1], z5_startCC, s5_rotate)
            s5_cells_normalized['Segment'] = 's5'

        # ---- s6 ----
        if s6_cells is not None:
            mask = (
                (additional_refs['Animal'] == curr_animal) &
                (additional_refs['Segment'] == 's6') &
                (additional_refs['Ref'] == 'CC')
            )
            if not additional_refs.loc[mask].empty:
                print('used additional refs for ', curr_animal, '_s6')
                x6_startCC = additional_refs.loc[mask, 'Position X'].iloc[0]
                y6_startCC = additional_refs.loc[mask, 'Position Y'].iloc[0]
                x6_endCC   = additional_refs.loc[mask, 'Position X'].iloc[0]
                y6_endCC   = additional_refs.loc[mask, 'Position Y'].iloc[0]
            else:
                x6_startCC = a_ref.loc[(a_ref['Segment']=='s6') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position X'].iloc[0]
                y6_startCC = a_ref.loc[(a_ref['Segment']=='s6') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Y'].iloc[0]
                x6_endCC   = a_ref.loc[(a_ref['Segment']=='s6') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position X'].iloc[0]
                y6_endCC   = a_ref.loc[(a_ref['Segment']=='s6') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Y'].iloc[0]
            
            z6_startCC = a_ref.loc[(a_ref['Segment']=='s6') & (a_ref['Ref']=='CC') & (a_ref['Location']=='rostral'), 'Position Z'].iloc[0]
            z6_endCC   = a_ref.loc[(a_ref['Segment']=='s6') & (a_ref['Ref']=='CC') & (a_ref['Location']=='caudal'),  'Position Z'].iloc[0]

            s6_midpoint = (np.array([x6_startCC, y6_startCC, z6_startCC]) + np.array([x6_endCC, y6_endCC, z6_endCC])) / 2
            s6_cells_normalized = caf.normalizePoints2(s6_cells, s6_midpoint[0], s6_midpoint[1], z6_startCC, s6_rotate)
            s6_cells_normalized['Segment'] = 's6'
            
        # Determine the length of each segment
        if s1_cells is not None:
            s1_length = abs(a_ref.loc[(a_ref['Segment'] == 's1') & (a_ref['Ref'] == 'CC') & (a_ref['Location'] == 'caudal'), 'Position Z'].iloc[0] - z1_startCC)
        else:
            s1_length = 1399 # The average segment lengths used here are calculated empirically from 3-5 samples of each segment (see Excel doc 'Avg Segment Length' in Code folder)
        s1_startz = 0 # s1 start Z Position in the space where all segments are plotted together
        s1_endz = copy.deepcopy(s1_length) # s1 end Z Position in the space where all segments are plotted together
        
        if s2_cells is not None:
            s2_length = abs(a_ref.loc[(a_ref['Segment'] == 's2') & (a_ref['Ref'] == 'CC') & (a_ref['Location'] == 'caudal'), 'Position Z'].iloc[0] - z2_startCC)    
        else:
            s2_length = 2069
        s2_startz = s1_endz
        s2_endz = s2_startz + s2_length
        
        if s3_cells is not None:
            s3_length = abs(a_ref.loc[(a_ref['Segment'] == 's3') & (a_ref['Ref'] == 'CC') & (a_ref['Location'] == 'caudal'), 'Position Z'].iloc[0] - z3_startCC)
        else:
            s3_length = 2720
        s3_startz = s2_endz
        s3_endz = s3_startz + s3_length
        
        if s4_cells is not None:
            s4_length = abs(a_ref.loc[(a_ref['Segment'] == 's4') & (a_ref['Ref'] == 'CC') & (a_ref['Location'] == 'caudal'), 'Position Z'].iloc[0] - z4_startCC)
        else:
            s4_length = 2608
        s4_startz = s3_endz
        s4_endz = s4_startz + s4_length
        
        if s5_cells is not None:
            s5_length = abs(a_ref.loc[(a_ref['Segment'] == 's5') & (a_ref['Ref'] == 'CC') & (a_ref['Location'] == 'caudal'), 'Position Z'].iloc[0] - z5_startCC)
        else:
            s5_length = 2945
        s5_startz = s4_endz
        s5_endz = s5_startz + s5_length
        
        if s6_cells is not None:
            s6_length = abs(a_ref.loc[(a_ref['Segment'] == 's6') & (a_ref['Ref'] == 'CC') & (a_ref['Location'] == 'caudal'), 'Position Z'].iloc[0] - z6_startCC)
        else:
            s6_length = 3009
        s6_startz = s5_endz
        s6_endz = s6_startz + s6_length
              
        boundary_list = [s1_startz, s1_endz, s2_startz, s2_endz, s3_startz, s3_endz, s4_startz, s4_endz, s5_startz, s5_endz, s6_startz, s6_endz] # List of all start and end points for segments
        in_bt_boundary_list = [s2_startz, s3_startz,s4_startz,s5_startz,s6_startz]# List of the boundaries between segments
        
        # Add the cumulative length of the cord to each of the Z position of cells in each segment. Skip s1, as those cells can start at CC Z Position = 0
        # Concatenate all the individual segment dataframes into one dataframe
        cumulative_length = s1_length
        a_cells_normalized = copy.deepcopy(s1_cells_normalized)
        if s2_cells is not None:
            s2_cells_normalized['Position Z'] += s1_length
            a_cells_normalized = pd.concat([a_cells_normalized, s2_cells_normalized]) 
        
        cumulative_length += s2_length
        if s3_cells is not None:
            s3_cells_normalized['Position Z'] += cumulative_length
            a_cells_normalized = pd.concat([a_cells_normalized, s3_cells_normalized]) 
        
        cumulative_length += s3_length
        if s4_cells is not None:
            s4_cells_normalized['Position Z'] += cumulative_length
            a_cells_normalized = pd.concat([a_cells_normalized, s4_cells_normalized]) 
        
        cumulative_length += s4_length 
        if s5_cells is not None:
            s5_cells_normalized['Position Z'] += cumulative_length
            a_cells_normalized = pd.concat([a_cells_normalized, s5_cells_normalized]) 
            
        cumulative_length += s5_length 
        if s6_cells is not None:
            s6_cells_normalized['Position Z'] += cumulative_length
            a_cells_normalized = pd.concat([a_cells_normalized, s6_cells_normalized]) 
        
        # Add the spinal cord length, max spread in x axis, and max spread in y axis to lists
        length_list.append(s6_endz)
        x_range_list.append(a_cells_normalized['Position X'].max() - a_cells_normalized['Position X'].min())
        y_range_list.append(a_cells_normalized['Position Y'].max() - a_cells_normalized['Position Y'].min())
        
        # If the option is selected, normalize the Z positions to an average spinal cord length (see Excel doc 'Avg Segment Length' in Code folder)
        avg_total_length = 14749
        if normalize_z_length is True:
            a_cells_normalized['Position Z'] = (a_cells_normalized['Position Z'] / s6_endz) * avg_total_length
            boundary_list = (boundary_list / max(boundary_list)) * avg_total_length
        

        ################################################################################################################################################
        ### Sort cells into clones 
        
        ## TO DO currently, clones that are beyond the boundary, but at a greater distance than the threshold from it, may be included. (happens when the end is uneven and the boundary is a bit inside the segment)
        # Sort cells in ascending Z Position order 
        a_cells_normalized = a_cells_normalized.sort_values('Position Z')
        # Subset the dataframe to extract only the red and green cells
        rg_cells_normalized = a_cells_normalized[(a_cells_normalized['Colour'] == 'r') | (a_cells_normalized['Colour'] == 'g') ]
        
        num_rows = rg_cells_normalized.shape[0] # Get number of rows in dataframe
        
        clone_list = ['c1'] # Initialize list of clone labels with the label for the first cell
        i = 0 # Initialize index for clone_list
        excluded_clones_boundary = [] # Initialize a list for clones that are excluded based on proximity to a boundary
        excluded_clones_proximity = [] # # Initialize a list for clones that are excluded specifically because of proximity to another clone (exclusion zone)
        curr_clone = 1
        excl_message_boundary = 'exclude - too close to boundary'
        excl_message_damage = 'exclude - in vicinity of tissue damage'
        excl_message_proximity = 'exclude - too close to neighbouring clone'
        excl_message_1_cell = 'exclude - single-cell clone'
        # Iterate through all the cells and assign them clone labels
        for row in range (1, num_rows):
            previous_row = row - 1
            diff = rg_cells_normalized.iloc[row].loc['Position Z'] - rg_cells_normalized.iloc[previous_row].loc['Position Z'] # Z Position of current cell minus z position of previous cell
            if diff < inter_cell_threshold:
                clone_list.append(str('c' + str(curr_clone))) # If the distance between this cell and the previous cell is less than the threshold, label it as part of the same clone.
                i += 1 
            else:
                curr_clone += 1
                clone_list.append(str('c' + str(curr_clone))) # If the distance between this cell and the previous cell is greater than the threshold, label it as a new clone for now.
                if exclusion_zone and (inter_cell_threshold <= diff <= inter_clone_threshold): # If exclusion_zone has been set to True and the cell falls within this zone, add its name and the name of the clone before it to the list of clones to exclude.
                    excluded_clones_proximity.append(str('c' + str(curr_clone)))
                    all_exclud_clones_proximity.append(str('c' + str(curr_clone)))
                    excluded_clones_proximity.append(str('c' + str(curr_clone - 1)))
                    all_exclud_clones_proximity.append(str('c' + str(curr_clone - 1)))
                i += 1
            # Calculate the z distance between the current cell and all the segment boundaries in the a_ref dataframe. If any of the distances is less than the threshold, write "exclude" in the
            # "Clone" list/future column and skip the cell-cell threshold step.
            for z in boundary_list:
                dist = abs(rg_cells_normalized.iloc[row].loc['Position Z'] - z)
                # If the cell-boundary distance is less than the threshold, add the current clone name to a list that you will use to exclude all the other cells
                # belonging to the same clone later.
                if dist < cell_boundary_threshold:
                    excluded_clones_boundary.append(str('c' + str(curr_clone)))
                    all_exclud_clones_boundary.append(str('c' + str(curr_clone)))
                    break

        rg_cells_normalized['Clone'] = clone_list # Add the list of clone labels as a new column
              
        # Add sample information to the dataframes
        desired_columns =  ['Litter', 'TM', 'TM Amount', 'Animal', 'Segment', 'ID', 'Position X', 'Position Y', 'Position Z', 
                   'Clone', 'Colour', 'D/V', 'G/W', 'Hemisphere', 'Damage', 'Neuron', 'MN', 'Commissural neuron',
                   'Midline cell', 'Glial cell', 'Uncertain', 'Chx10', 'ChAT', 'Lhx5', 'Notes']
        animal_name = curr_dict['Litter'] + '_' + curr_dict['Animal'] # Useful as a unique identifier for each animal
        TM_amount = curr_dict['TM Amount']
        a_cells_normalized['Animal'] = animal_name
        a_cells_normalized['TM Amount'] = TM_amount
        a_cells_normalized['Litter'] = curr_dict['Litter']
        a_cells_normalized['TM'] = 'E' + str(curr_dict['TM Injection Timepoint E'])
        for col in desired_columns:
            if col not in a_cells_normalized.columns:
                a_cells_normalized[col] = np.nan
        a_cells_normalized = a_cells_normalized[desired_columns]
        rg_cells_normalized['Animal'] = animal_name
        rg_cells_normalized['TM Amount'] = TM_amount
        rg_cells_normalized['Litter'] = curr_dict['Litter']
        rg_cells_normalized['TM'] = 'E' + str(curr_dict['TM Injection Timepoint E'])
        for col in desired_columns:
            if col not in rg_cells_normalized.columns:
                rg_cells_normalized[col] = np.nan
        rg_cells_normalized = rg_cells_normalized[desired_columns]
        
        rg_cells_normalized_noexcl = copy.deepcopy(rg_cells_normalized) # Make a copy of the df where clones have no exclusion messages
        
        ### Label clones with exclusion messages so that they can be excluded later.
        # Label all cells belonging to clones to be excluded based on damage in the transverse plane in the same or opposite hemisphere
        excluded_clones_damage = rg_cells_normalized.loc[(rg_cells_normalized['Damage'] == 'yes') | (rg_cells_normalized['Damage'] == 'con') | (rg_cells_normalized['Damage'] == 'ips'), 'Clone'].tolist()
        all_exclud_clones_damage.append(excluded_clones_damage)
        rg_cells_normalized['Clone'] = rg_cells_normalized['Clone'].apply(lambda x: excl_message_damage if x in excluded_clones_damage else x)
        
        # Label all cells belonging to clones to be excluded based on proximity to another clone with the exclusion message
        rg_cells_normalized['Clone'] = rg_cells_normalized['Clone'].apply(lambda x: excl_message_proximity if x in excluded_clones_proximity else x)
        
        # Label all cells belonging to clones to be excluded based on proximity to the boundary with the exclusion message
        rg_cells_normalized['Clone'] = rg_cells_normalized['Clone'].apply(lambda x: excl_message_boundary if x in excluded_clones_boundary else x)
        
        # Collect all the cells from all animals into a single file
        if first_pass == True:
            overall_rg_cells = copy.deepcopy(rg_cells_normalized)
            overall_rg_cells_noexcl = copy.deepcopy(rg_cells_normalized_noexcl)
            overall_a_cells = copy.deepcopy(a_cells_normalized)
            first_pass = False
        else: 
            overall_rg_cells = pd.concat([overall_rg_cells, rg_cells_normalized])
            overall_rg_cells_noexcl = pd.concat([overall_rg_cells_noexcl, rg_cells_normalized_noexcl])
            overall_a_cells = pd.concat([overall_a_cells, a_cells_normalized])
        
        num_cells_list_all.append(a_cells_normalized.shape[0])
        num_cells_list_rg.append(rg_cells_normalized.shape[0])
        
    ### Exclude cells and save files  
    # Label all 1-cell clones with the exclusion message
    overall_rg_cells_withexcl_no_1_cell_temp = copy.deepcopy(overall_rg_cells)
    clone_counts = overall_rg_cells_withexcl_no_1_cell_temp.groupby(["Animal", "Litter", "Segment", "Clone"]).size().reset_index(name = 'Clone Counts')
        
    overall_rg_cells_withexcl_no_1_cell_temp = overall_rg_cells_withexcl_no_1_cell_temp.merge(clone_counts, on = ["Animal", "Litter", "Segment", "Clone"], how = "left")
        
    single_cell_mask = overall_rg_cells_withexcl_no_1_cell_temp['Clone Counts'] == 1 # Find rows where the combination is unique (count == 1) (i.e. single-cell clones)
        
    num_exclud_cells_1_cell = single_cell_mask.sum() # Count the number of single-cell clones     
        
    overall_rg_cells_withexcl_no_1_cell_temp.loc[single_cell_mask, 'Clone'] = excl_message_1_cell # Add exclusion message
        
    overall_rg_cells_withexcl_no_1_cell_temp.drop(columns = ['Clone Counts'], inplace = True)
       
    # Write the file containing all rg cells from all animals to csv file
    file_name = overall_output_file_path + '/Overall_rg_cells_' + curr_date + '.csv'
    overall_rg_cells_withexcl_no_1_cell_temp.to_csv(file_name, mode = 'w', index = False) # Will overwrite existing file 
    file_name = overall_output_file_path + '/overall_rg_cells_noexcl' + curr_date + '.csv'
    overall_rg_cells_noexcl.to_csv(file_name, mode = 'w', index = False) # Will overwrite existing file
    file_name = overall_output_file_path +  '/Overall_all_cells_' + curr_date + '.csv'
    overall_a_cells.to_csv(file_name, mode = 'w', index = False) # Will overwrite existing file     
    
    # Get a dataframe that omits all excluded clones except single-cell clones
    overall_rg_cells_withexcl = overall_rg_cells[~overall_rg_cells['Clone'].isin([excl_message_boundary, excl_message_damage, excl_message_proximity])]
    
    # Get a dataframe that omits all excluded clones
    overall_rg_cells_withexcl_no_1_cell = overall_rg_cells_withexcl_no_1_cell_temp[~overall_rg_cells_withexcl_no_1_cell_temp['Clone'].isin([excl_message_boundary, excl_message_damage, excl_message_proximity, excl_message_1_cell])]
    
    # Initialize a df that will contain all clones and their quantified features
    clone_codes = overall_rg_cells_withexcl_no_1_cell.groupby(['Animal', 'Clone']).size().reset_index().drop(columns = 0) # Get all the unique clones as clone codes
    clone_codes = clone_codes['Animal'].astype(str) + '_' + clone_codes['Clone'].astype(str)
    rg_clones_df = pd.DataFrame(clone_codes, columns = ['Clone'])
    
    ################################################################################################################################################  
    ### Plot the distribution of nearest-neighbour cell XYZ distances, similar to how it's done in Gao et al., 2014
    # Perform the nearest-neighbour analysis on the entire dataset without excluding cells. Look for nearest neighbour in xyz axis. Nearest neighbours are found within each animal, then all the distances are pooled to be plotted together.
    # Get rg cells and y cells in separate dataframes
    
    # Average number of cells
    num_cells_list_all = [int(a) for a in num_cells_list_all]
    average_num_cells = int(mean(num_cells_list_all))
    # Average spinal cord length (z)
    length_list = [int(a) for a in length_list]
    average_length = mean(length_list)
    # The x range in the average max spread of cells in a spinal along the x axis
    average_x_range = (0, mean(x_range_list))
    # The y range in the average max spread of cells in a spinal along the y axis
    average_y_range = (0, mean(y_range_list))
    

    output_path = overall_output_file_path + '/Nearest-neighbour Distance_' + curr_date
    cap.NNDPlot(overall_a_cells, 
            average_length, 
            average_x_range, 
            average_y_range,
            average_num_cells,
            output_path = output_path,
            inter_cell_threshold = 125,
            exclude_single_y = True)
    
    


    ################################################################################################################################################  
    ### Clone size (cells per clone)
    # Store the number of cells of each type per clone in rg_clones_df
    clone_size_E9 = caf.cellsPerClone(overall_rg_cells_withexcl_no_1_cell, 'E9.5')
    clone_size_E10 = caf.cellsPerClone(overall_rg_cells_withexcl_no_1_cell, 'E10.5')
    clone_size_E11 = caf.cellsPerClone(overall_rg_cells_withexcl_no_1_cell, 'E11.5')
    clone_size_df = pd.concat([clone_size_E9, clone_size_E10, clone_size_E11])
    rg_clones_df = pd.merge(rg_clones_df, clone_size_df, on = 'Clone', how = 'left')

################################################################################################################################################  
    # ### Plot the intra-clone spread
    # # In xy plane
    intra_spread_E8_df = caf.intraCloneSpread(overall_rg_cells_withexcl_no_1_cell, 'E8.5', ['Position X', 'Position Y'])
    intra_spread_E9_df = caf.intraCloneSpread(overall_rg_cells_withexcl_no_1_cell, 'E9.5', ['Position X', 'Position Y'])
    intra_spread_E10_df = caf.intraCloneSpread(overall_rg_cells_withexcl_no_1_cell, 'E10.5', ['Position X', 'Position Y'])
    intra_spread_E11_df = caf.intraCloneSpread(overall_rg_cells_withexcl_no_1_cell, 'E11.5', ['Position X', 'Position Y'])

    # Get data into long form to make sns plots
    # Pad the lists with nans to make them all the same length
    max_length = max(len(intra_spread_E8_df['Spread'].tolist()), len(intra_spread_E9_df['Spread'].tolist()), len(intra_spread_E10_df['Spread'].tolist()), len(intra_spread_E11_df['Spread'].tolist()))
    intra_spread_E8_list = intra_spread_E8_df['Spread'].tolist() + [np.nan] * (max_length - len(intra_spread_E8_df['Spread'].tolist()))
    intra_spread_E9_list = intra_spread_E9_df['Spread'].tolist() + [np.nan] * (max_length - len(intra_spread_E9_df['Spread'].tolist()))
    intra_spread_E10_list = intra_spread_E10_df['Spread'].tolist() + [np.nan] * (max_length - len(intra_spread_E10_df['Spread'].tolist()))
    intra_spread_E11_list = intra_spread_E11_df['Spread'].tolist() + [np.nan] * (max_length - len(intra_spread_E11_df['Spread'].tolist()))
    
    intra_spread_df = pd.DataFrame({'E8.5': intra_spread_E9_list,
                                    'E9.5': intra_spread_E9_list,
                                    'E10.5': intra_spread_E10_list, 
                                    'E11.5': intra_spread_E11_list})
    intra_spread_full_df =  pd.concat([intra_spread_E8_df, intra_spread_E9_df, intra_spread_E10_df, intra_spread_E11_df], ignore_index = True)   
    intra_spread_df_long = intra_spread_df.melt(ignore_index = False, var_name = 'TM', value_name = 'Intra-spread').reset_index()

    
    intra_spread_subset = intra_spread_full_df[['Clone','Spread', 'Cell_dist_avg', 'Cell_dist_sd']]
    intra_spread_subset = intra_spread_subset.rename(columns = {'Spread': 'Spread XY', 'Cell_dist_avg': 'Cell_dist_avg_xy', 'Cell_dist_sd': 'Cell_dist_sd_xy'})
    rg_clones_df = pd.merge(rg_clones_df, intra_spread_subset, on = 'Clone', how = 'left')
    
    


    # Along z axis
    intra_spread_z_E8_df = caf.intraCloneSpread(overall_rg_cells_withexcl_no_1_cell, 'E8.5', ['Position Z'])
    intra_spread_z_E9_df = caf.intraCloneSpread(overall_rg_cells_withexcl_no_1_cell, 'E9.5', ['Position Z'])
    intra_spread_z_E10_df = caf.intraCloneSpread(overall_rg_cells_withexcl_no_1_cell, 'E10.5', ['Position Z'])
    intra_spread_z_E11_df = caf.intraCloneSpread(overall_rg_cells_withexcl_no_1_cell, 'E11.5', ['Position Z'])

    intra_spread_z_full_df =  pd.concat([intra_spread_z_E8_df, intra_spread_z_E9_df, intra_spread_z_E10_df, intra_spread_z_E11_df], ignore_index = True)
    intra_spread_z_subset = intra_spread_z_full_df[['Clone','Spread', 'Cell_dist_avg', 'Cell_dist_sd']]
    intra_spread_z_subset = intra_spread_z_subset.rename(columns = {'Spread': 'Spread Z', 'Cell_dist_avg': 'Cell_dist_avg_Z', 'Cell_dist_sd': 'Cell_dist_sd_Z'})
    rg_clones_df = pd.merge(rg_clones_df, intra_spread_z_subset, on = 'Clone', how = 'left')

################################################################################################################################################  
    ### D/V cell fractions
    num_clones_E8, num_animals_E8, dv_data_df_E8 = caf.DVRatio(overall_rg_cells_withexcl_no_1_cell, 'E8.5')
    num_clones_E9, num_animals_E9, dv_data_df_E9 = caf.DVRatio(overall_rg_cells_withexcl_no_1_cell, 'E9.5')
    num_clones_E10, num_animals_E10, dv_data_df_E10 = caf.DVRatio(overall_rg_cells_withexcl_no_1_cell, 'E10.5')
    num_clones_E11, num_animals_E11, dv_data_df_E11 = caf.DVRatio(overall_rg_cells_withexcl_no_1_cell, 'E11.5')
    
    dv_data_df = pd.concat([dv_data_df_E8, dv_data_df_E9, dv_data_df_E10, dv_data_df_E11])
    dv_data_df_subset = dv_data_df[['Clone','DV']]
    rg_clones_df = pd.merge(rg_clones_df, dv_data_df_subset, on = 'Clone', how = 'left')
    



    ##################################################################################################################################################  
    ### Lhx5 expression (adds per-clone expression classification columns to rg_clones_df)
    # Make folder for co-expression analysis output
    expression_output_path = overall_output_file_path + '/Clone Co-expression'
    os.makedirs(expression_output_path, exist_ok = True)

    Lhx5_proportions, x_axis_ticks, rg_clones_df, mixed_express_info = cao.plotExpression(
        overall_rg_cells_withexcl_no_1_cell,
        rg_clones_df,
        'Lhx5'
    )

    # Save the information about colour composition of mixed clones
    file_name = expression_output_path + '/Mixed_expression_clones_info_' + curr_date + '.csv'
    mixed_express_info.to_csv(file_name, mode = 'w', index = False)

    # ################################################################################################################################################  
    ### Proportions of clone types (asymmetric, symmetric, small neurogenic, unicolor)
    
    #### WARNING Check that the n reporting here is correct!
    num_clones_E8, num_animals_E8, clone_type_E8 = caf.cloneType(overall_rg_cells_withexcl_no_1_cell, 'E8.5', unicolor = True)
    num_clones_E9, num_animals_E9, clone_type_E9 = caf.cloneType(overall_rg_cells_withexcl_no_1_cell, 'E9.5', unicolor = True)
    num_clones_E10, num_animals_E10, clone_type_E10 = caf.cloneType(overall_rg_cells_withexcl_no_1_cell, 'E10.5', unicolor = True)
    num_clones_E11, num_animals_E11, clone_type_E11 = caf.cloneType(overall_rg_cells_withexcl_no_1_cell, 'E11.5', unicolor = True)
        
    type_data = pd.concat([clone_type_E8, clone_type_E9, clone_type_E10, clone_type_E11])
    type_data_subset = type_data[['Clone','Type']]
    rg_clones_df = pd.merge(rg_clones_df, type_data_subset, on = 'Clone', how = 'left')

    
    ################################################################################################################################################  
    ### Clone laterality, i.e. the proportions of clones that are uni- and bilateral (comprised of cells in only one or both hemispheres)
    # Excludes one-cell clones (that are necessarily unilateral)
    #num_clones_E8, num_animals_E8, clone_lat_E8, clone_lat_size_E8 = cloneLat(overall_rg_cells_withexcl_no_1_cell, 'E8.5')
    num_clones_E9, num_animals_E9, clone_lat_E9, clone_lat_size_E9 = caf.cloneLat(overall_rg_cells_withexcl_no_1_cell, 'E9.5')
    num_clones_E10, num_animals_E10, clone_lat_E10, clone_lat_size_E10 = caf.cloneLat(overall_rg_cells_withexcl_no_1_cell, 'E10.5')
    num_clones_E11, num_animals_E11, clone_lat_E11, clone_lat_size_E11 = caf.cloneLat(overall_rg_cells_withexcl_no_1_cell, 'E11.5')
    
    # Concatenate the individual TM dataframes
    clone_lat = pd.concat([clone_lat_E9, clone_lat_E10, clone_lat_E11])
    clone_lat_size = pd.concat([clone_lat_size_E9, clone_lat_size_E10, clone_lat_size_E11])
    clone_lat_size_2 = clone_lat_size[clone_lat_size['TM'] == 'E10.5']
    clone_lat_size_subset = clone_lat_size[['Clone','Laterality']]
    rg_clones_df = pd.merge(rg_clones_df, clone_lat_size_subset, on = 'Clone', how = 'left')


    ################################################################################################################################################
    ### Write a csv file containing the parameters used for the current analysis as well as descriptive statistics
    file_name = overall_output_file_path + '/Descriptive_Stats_' + curr_date + '.csv'
    parameters = pd.DataFrame()
    num_clones_after_exclusion = overall_rg_cells_withexcl_no_1_cell.groupby(['Litter', 'Animal', 'Clone']).ngroups
    num_clones_after_exclusion_names = list(overall_rg_cells_withexcl_no_1_cell.groupby(['Litter', 'Animal', 'Clone']).groups.keys())
    num_exclud_cells_boundary = overall_rg_cells['Clone'].value_counts()[excl_message_boundary] # Count the number of R/G cells excluded by the chosen cell-boundary threshold
    num_exclud_cells_damage = overall_rg_cells['Clone'].value_counts()[excl_message_damage] # Count the number of R/G cells excluded due to tissue damage
    if exclusion_zone:
        vc = overall_rg_cells["Clone"].value_counts()
        num_exclud_cells_proximity = int(vc.get(excl_message_proximity, 0))
    else:
        num_exclud_cells_proximity = 0
    num_exclud_clones_boundary = len(set(all_exclud_clones_boundary))
    num_exclud_clones_proximity = len(set(all_exclud_clones_proximity))
    all_exclud_clones_damage = [item for sublist in all_exclud_clones_damage for item in sublist] # Flatten list of lists
    num_exclud_clones_damage = len(set(all_exclud_clones_damage))
    
    parameters = pd.DataFrame([{'Date and Time': curr_date, 
                                'Inter-clone threshold' : inter_clone_threshold, 
                                'Cell-boundary threshold' : cell_boundary_threshold,
                                'Exclusion zone': exclusion_zone,
                                'Normalize Z Length': normalize_z_length, 
                                'Total number of R/G cells': len(overall_rg_cells),
                                'Number of R/G clones': num_clones_after_exclusion,
                                'Number of excluded cells due to boundary' : num_exclud_cells_boundary, 
                                'Number of excluded cells due to tissue damage': num_exclud_cells_damage, 
                                'Number of excluded cells due to proximity to another clone': num_exclud_cells_proximity,
                                'incorrect -Number of excluded R/G clones': num_exclud_clones_boundary + num_exclud_clones_proximity + num_exclud_clones_damage + num_exclud_cells_1_cell,
                                'incorrect -Number of excluded clones due to boundary': num_exclud_clones_boundary,
                                'incorrect -Number of excluded clones due to tissue damage': num_exclud_clones_damage,
                                'incorrect -Number of excluded clones due to proximity to another clone': num_exclud_clones_proximity,
                                'Number of excluded single-cell clones': num_exclud_cells_1_cell}])
    # #parameters.transpose()
    # for key, value in clone_counts.items(): # Add the clone counts to the parameters dataframe
    #     parameters[key] = value
    # for key, value in animal_counts.items(): # Add the clone counts to the parameters dataframe
    #     parameters[key] = value
    parameters.to_csv(file_name, mode = 'w', index = False)
    
    # Write the df containing clone descriptors (rg_clone_df) to a csv file
    file_name = overall_output_file_path + '/_rg_clones_' + curr_date + '.csv'
    rg_clones_df.to_csv(file_name, mode = 'w', index = False) # Will overwrite existing file
    print('\nrg_clones_df', rg_clones_df)
