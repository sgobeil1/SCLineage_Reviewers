import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.spatial import ConvexHull
import numpy as np
import copy
import matplotlib.ticker as ticker
from scipy.spatial import distance_matrix
from scipy.spatial.distance import pdist, squareform
from scipy import stats

# Read CSV using pandas
def read_csv(file_path, skip_rows):
    if file_path == 'None':
        points = None
    else: 
        points = pd.read_csv(file_path, skiprows=skip_rows) # skip initial rows
        points.dropna(axis = 0, how = 'all', inplace = True) # drop empty rows
    
    return points

# Find x, y, z coordinates in 'points' DataFrame where value 'ID' is found in column named 'column'
def findRefPoint(points, column, ID):
    x, y, z = points.loc[points[column] == ID, 'Position X'].iloc[0], points.loc[points[column] == ID, 'Position Y'].iloc[0], points.loc[points[column] == ID, 'Position Z'].iloc[0]
    
    return x, y, z

# Find the distance between the 2 z positions of the reference points in the two lists of lists where each item has the format [ID, x, y, z]
def findDistance(points1, points2, axis, ID):
    if axis == "x":
        index = 0
    if axis == "y":
        index = 1
    if axis == "z":
        index = 2
    ref_1 = findRefPoint(points1, ID)[index]
    ref_2 = findRefPoint(points2, ID)[index]
    dist = abs(ref_1 - ref_2)
    
    return dist

# Find the scaling factor relating the distance between two points ('ID1' and 'ID2') in 'points' to the distnce between two points in 'ref' DataFrames. Axis specifies wether we're interested in the distance along the
# x, y, or z axis.
def scalingFactor(points, ref, axis, ID1, ID2):
    if axis == "x":
        index = 0
    if axis == "y":
        index = 1
    if axis == "z":
        index = 2
    points_ID1 = findRefPoint(points, 'Ref', ID1)[index]
    points_ID2 = findRefPoint(points, 'Ref', ID2)[index]
    points_range = points_ID1 - points_ID2
    ref_ID1 = findRefPoint(ref, 'Ref', ID1)[index]
    ref_ID2 = findRefPoint(ref, 'Ref', ID2)[index]
    ref_range = ref_ID1 - ref_ID2
    scale_fact = points_range / ref_range
    
    return scale_fact

# Rotate x,y,z coordinates in columns of a data frame named 'Position X', 'Position Y', 'Position Z' by specified degrees ("rotation" argument)
def rotatePoints(points, rotation):
    # Rotation
    angle_radians = np.deg2rad(rotation)
    cos_angle = np.cos(angle_radians)
    sin_angle = np.sin(angle_radians)
    
    # Define the rotation matrix for the given axis
    rotation_matrix = np.array([
            [cos_angle, -sin_angle, 0],
            [sin_angle, cos_angle, 0],
            [0, 0, 1]
        ])
    
    # Extract the coordinates as a numpy array
    coordinates = points[['Position X', 'Position Y', 'Position Z']].values

    # Apply the rotation matrix
    rotated_coordinates = coordinates @ rotation_matrix.T

    # Update a new DataFrame with the new rotated coordinates
    points1 = copy.deepcopy(points)
    points1['Position X'], points1['Position Y'], points1['Position Z'] = rotated_coordinates.T
    
    return points1
    
# Set the central canal point to 0,0,0. Its ID must be 0. Rotate by specified degrees so that ventral is down.
def normalizePoints(points, ref_points, rotation, column, ID):
    # Find the reference point where ID = 0
    ref_x, ref_y, ref_z = findRefPoint(ref_points, column, ID)
    
    # Normalize all points relative to the reference point. Create a copy of the points dataframe and subtract the ref point from every value in the respective column.
    normalized_points = copy.deepcopy(points)
    normalized_points['Position X'] = points['Position X'] - ref_x
    normalized_points['Position Y'] = points['Position Y'] - ref_y
    normalized_points['Position Z'] = points['Position Z'] - ref_z
    
    # Rotate the points so by specified degrees
    rotated_points = rotatePoints(normalized_points, rotation)
    
    return rotated_points

# Set the central canal point to 0,0,0. Its ID must be 0. Rotate by specified degrees so that ventral is down.
def normalizePoints2(points, ref_x, ref_y, ref_z, rotation):
    
    # Normalize all points relative to the reference point. Create a copy of the points dataframe and subtract the ref point from every value in the respective column.
    normalized_points = copy.deepcopy(points)
    normalized_points['Position X'] = points['Position X'] - ref_x
    normalized_points['Position Y'] = points['Position Y'] - ref_y
    normalized_points['Position Z'] = (points['Position Z'] - ref_z).abs() # Take the absolute value because we don't want any negative z positions
    
    # Rotate the points so by specified degrees
    rotated_points = rotatePoints(normalized_points, rotation)
    
    return rotated_points

# points is a list of lists with 3 items [x, y, z]
def scatter_plot(points, color, s, label, ax):
    xs = [sublist[0] for sublist in points]
    ys = [sublist[1] for sublist in points]
    zs = [sublist[2] for sublist in points]
    ax.scatter(xs, zs, ys,  c=color, s = s, label=label) # The axes input order is such that the y and z axis are switched
    
    #For z-axis limit
    ax.set_ylim(0, 17000) 
    
# Inputs points1 and points2 should be lists of xyz coordinates, without IDs
def plot_points_and_hull(points, z_dist, color, label, ax):
    
    # Plot points from the first set
    scatter_plot(points, color = color, s = 1, label = label, ax = ax)
    
    # Swap the second and third items in each sublist so that you can change the orientation of the z and y axis display later
    for sublist in points:
        sublist[1], sublist[2] = sublist[2], sublist[1]
    
    hull = ConvexHull(points)
        
    # Create the convex hull
    for simplex in hull.simplices:
        poly3d = [[points[vertex] for vertex in simplex]]
        ax.add_collection3d(Poly3DCollection(poly3d, facecolors='grey', linewidths=1, alpha=0.3))

    # Add the central canal
    CC = points[0]
    plt.plot([CC[0], CC[0]], [CC[1], CC[1] + z_dist], [CC[2], CC[2]], color = 'k', linestyle="-")

# cells is the cell coordinate dataframe, cells_ref is the cell reference coordinate dataframe, segment_ref is the reference coordinate dataframe for that segment, rotation is the rotation required to make ventral down.
def scaleCells(cells, cells_ref, ref_segment, ref_segment_length, rotation):
                
    cells = copy.deepcopy(cells)
        
    # Scale in x.
    x_scale = scalingFactor(cells_ref, ref_segment, axis = 'x', ID1 = 'LR', ID2 = 'CC')
    cells['Position X'] = cells['Position X'] / x_scale
            
    # Scale in y. Scaling depends on whether the cell is positioned above or below the CC (0 on the y axis)
    y_scale_belowCC = scalingFactor(cells_ref, ref_segment, axis = 'y', ID1 = 'VR', ID2 = 'CC')
    y_scale_aboveCC = scalingFactor(cells_ref, ref_segment, axis = 'y', ID1 = 'DR', ID2 = 'CC')
        
    cells['Position Y'] = cells['Position Y'].apply(lambda x: x / y_scale_belowCC if x < 0 else x / y_scale_aboveCC)
        
    # Scale in z.
    cells_length = float(cells_ref.loc[cells_ref['Ref'] == 'endz', 'Position Z']) - float(cells_ref.loc[cells_ref['Ref'] == 'startz', 'Position Z']) # Length of cells segment
    z_scale = cells_length / ref_segment_length
    cells['Position Z'] = cells['Position Z'] / z_scale
        
    return cells

# Takes two dataframes with xyz coordinates in columns named 'Position X', 'Position Y', 'Position Z', finds the nearest neighbour of
# points in df1 among points in df2, and returns a list of the distances between those nearest neighbours. Set same = True if the dataframes
# are the same.
def nearestNeighbourDistance(df1, df2, same):
    dist_matrix = distance_matrix(df1[['Position X', 'Position Y', 'Position Z']], df2[['Position X', 'Position Y', 'Position Z']]) # Calculate the pairwise distance matrix for cells in the dataframe
    if same == True:
        np.fill_diagonal(dist_matrix, np.inf) # Replace diagonal with np.inf to avoid finding the cell itself as the nearest neighbor
    nearest_indices = np.argmin(dist_matrix, axis = 1) # Find the index of the nearest neighbour for each point
    differences = df1[['Position X', 'Position Y', 'Position Z']].values - df2[['Position X', 'Position Y', 'Position Z']].values[nearest_indices] # Calculate the differences between each cell and its nearest neighbour
    distances = np.sqrt(np.sum(differences**2, axis = 1)) # Calculate the Euclidean distances and store in a list
    distances_list = distances.tolist()
    
    return distances_list

# Returns a list of counts of the number of cells per clone, as well as the average and SEM. Input a dataframe that includes columns 'TM', 'Animal', and 'Clone'. Choose desired TM injection time point.
def cellsPerClone(df, TM):
    # Returns a dataframe summarizing the number of cells per clone for a given TM stage, including counts of specific types of cells
    
    # Subset for TM time point
    df_TM = df[df['TM'] == TM]
    grouped = df_TM.groupby(['Animal', 'Clone'])

    num_data = {
        'Clone': [],
        'TM': [],
        'TM amount': [],
        'Segment': [],
        'Level': [],
        'Total cells': [],
        'Total cells_g': [],
        'Total cells_r': [],
        'Total cells_le': [],
        'Total cells_ri': [],
        'Total cells_le_g': [],
        'Total cells_le_r': [],
        'Total cells_ri_g': [],
        'Total cells_ri_r': [],

        'Neuron count (incl Comm,MN)': [],
        'Neuron count (incl Comm,MN)_g': [],
        'Neuron count (incl Comm,MN)_r': [],
        'Neuron count (incl Comm,MN)_le': [],
        'Neuron count (incl Comm,MN)_ri': [],

        'Glia count (PA,FA,WF,GL,OG,OL)': [],
        'Glia count (PA,FA,WF,GL,OG,OL)_g': [],
        'Glia count (PA,FA,WF,GL,OG,OL)_r': [],
        'Glia count (PA,FA,WF,GL,OG,OL)_le': [],
        'Glia count (PA,FA,WF,GL,OG,OL)_ri': [],

        'Astrocyte count (PA,FA,GL)': [],
        'Astrocyte count (PA,FA,GL)_g': [],
        'Astrocyte count (PA,FA,GL)_r': [],
        'Astrocyte count (PA,FA,GL)_le': [],
        'Astrocyte count (PA,FA,GL)_ri': [],

        'PA count': [],
        'PA count_g': [],
        'PA count_r': [],
        'PA count_le': [],
        'PA count_ri': [],

        'FA count': [],
        'FA count_g': [],
        'FA count_r': [],
        'FA count_le': [],
        'FA count_ri': [],

        'WF count': [],
        'WF count_g': [],
        'WF count_r': [],
        'WF count_le': [],
        'WF count_ri': [],

        'GL count': [],
        'GL count_g': [],
        'GL count_r': [],
        'GL count_le': [],
        'GL count_ri': [],

        'OG count': [],
        'OG count_g': [],
        'OG count_r': [],
        'OG count_le': [],
        'OG count_ri': [],

        'OL count': [],
        'OL count_g': [],
        'OL count_r': [],
        'OL count_le': [],
        'OL count_ri': [],

        'PR count': [],
        'PR count_g': [],
        'PR count_r': [],
        'PR count_le': [],
        'PR count_ri': [],

        'Uncertain count': [],
        'Uncertain count_g': [],
        'Uncertain count_r': [],
        'Uncertain count_le': [],
        'Uncertain count_ri': [],

        'Midline count (dv)': [],
        'Midline count (dv)_g': [],
        'Midline count (dv)_r': [],
        'Midline count (dv)_le': [],
        'Midline count (dv)_ri': [],

        'Midline d count': [],
        'Midline d count_g': [],
        'Midline d count_r': [],
        'Midline d count_le': [],
        'Midline d count_ri': [],

        'Midline v count': [],
        'Midline v count_g': [],
        'Midline v count_r': [],
        'Midline v count_le': [],
        'Midline v count_ri': [],

        'Commissural count_(dv)': [],
        'Commissural count_(dv)_g': [],
        'Commissural count_(dv)_r': [],
        'Commissural count_(dv)_le': [],
        'Commissural count_(dv)_ri': [],

        'Commissural d count': [],
        'Commissural d count_g': [],
        'Commissural d count_r': [],
        'Commissural d count_le': [],
        'Commissural d count_ri': [],

        'Commissural v count': [],
        'Commissural v count_g': [],
        'Commissural v count_r': [],
        'Commissural v count_le': [],
        'Commissural v count_ri': [],

        'MN count': [],
        'MN count_g': [],
        'MN count_r': [],
        'MN count_le': [],
        'MN count_ri': []
    }

    for name, group in grouped:
        # Basic clone info
        clone_code = str(name).strip('()').strip("'").replace("', '", '_')
        num_data['Clone'].append(str(clone_code)) # Store clone code
        num_data['TM'].append(TM)
        num_data['Total cells'].append(len(group))
        num_data['Total cells_g'].append(group['Colour'].eq('g').sum())
        num_data['Total cells_r'].append(group['Colour'].eq('r').sum())
        num_data['Total cells_le'].append(group['Hemisphere'].eq('l').sum())
        num_data['Total cells_ri'].append(group['Hemisphere'].eq('r').sum())
        num_data['Total cells_le_g'].append(((group['Colour'] == 'g') & (group['Hemisphere'] == 'l')).sum())
        num_data['Total cells_le_r'].append(((group['Colour'] == 'r') & (group['Hemisphere'] == 'l')).sum())
        num_data['Total cells_ri_g'].append(((group['Colour'] == 'g') & (group['Hemisphere'] == 'r')).sum())
        num_data['Total cells_ri_r'].append(((group['Colour'] == 'r') & (group['Hemisphere'] == 'r')).sum())
        num_data['TM amount'].append(group.iloc[0]['TM Amount'])
        segment = group.iloc[0]['Segment']
        num_data['Segment'].append(segment)
        num_data['Level'].append(getLevel(segment))

        # Category counts
        # Neurons (incl Comm, MN)
        num_data['Neuron count (incl Comm,MN)'].append((group['Neuron'] == 'yes').sum())
        num_data['Neuron count (incl Comm,MN)_g'].append(((group['Neuron'] == 'yes') & (group['Colour'] == 'g')).sum())
        num_data['Neuron count (incl Comm,MN)_r'].append(((group['Neuron'] == 'yes') & (group['Colour'] == 'r')).sum())
        num_data['Neuron count (incl Comm,MN)_le'].append(((group['Neuron'] == 'yes') & (group['Hemisphere'] == 'l')).sum())
        num_data['Neuron count (incl Comm,MN)_ri'].append(((group['Neuron'] == 'yes') & (group['Hemisphere'] == 'r')).sum())

        # Glia (PA, FA, WF, GL, OG, OL)
        num_data['Glia count (PA,FA,WF,GL,OG,OL)'].append(group['Glial cell'].isin(['pa', 'fa', 'wf', 'gl', 'og', 'ol']).sum())
        num_data['Glia count (PA,FA,WF,GL,OG,OL)_g'].append((group['Glial cell'].isin(['pa', 'fa', 'wf', 'gl', 'og', 'ol']) & (group['Colour'] == 'g')).sum())
        num_data['Glia count (PA,FA,WF,GL,OG,OL)_r'].append((group['Glial cell'].isin(['pa', 'fa', 'wf', 'gl', 'og', 'ol']) & (group['Colour'] == 'r')).sum())
        num_data['Glia count (PA,FA,WF,GL,OG,OL)_le'].append((group['Glial cell'].isin(['pa', 'fa', 'wf', 'gl', 'og', 'ol']) & (group['Hemisphere'] == 'l')).sum())
        num_data['Glia count (PA,FA,WF,GL,OG,OL)_ri'].append((group['Glial cell'].isin(['pa', 'fa', 'wf', 'gl', 'og', 'ol']) & (group['Hemisphere'] == 'r')).sum())

        # Astrocytes (PA, FA, GL)
        num_data['Astrocyte count (PA,FA,GL)'].append(group['Glial cell'].isin(['pa', 'fa', 'gl']).sum())
        num_data['Astrocyte count (PA,FA,GL)_g'].append((group['Glial cell'].isin(['pa', 'fa', 'gl']) & (group['Colour'] == 'g')).sum())
        num_data['Astrocyte count (PA,FA,GL)_r'].append((group['Glial cell'].isin(['pa', 'fa', 'gl']) & (group['Colour'] == 'r')).sum())
        num_data['Astrocyte count (PA,FA,GL)_le'].append((group['Glial cell'].isin(['pa', 'fa', 'gl']) & (group['Hemisphere'] == 'l')).sum())
        num_data['Astrocyte count (PA,FA,GL)_ri'].append((group['Glial cell'].isin(['pa', 'fa', 'gl']) & (group['Hemisphere'] == 'r')).sum())

        # PA
        num_data['PA count'].append((group['Glial cell'] == 'pa').sum())
        num_data['PA count_g'].append(((group['Glial cell'] == 'pa') & (group['Colour'] == 'g')).sum())
        num_data['PA count_r'].append(((group['Glial cell'] == 'pa') & (group['Colour'] == 'r')).sum())
        num_data['PA count_le'].append(((group['Glial cell'] == 'pa') & (group['Hemisphere'] == 'l')).sum())
        num_data['PA count_ri'].append(((group['Glial cell'] == 'pa') & (group['Hemisphere'] == 'r')).sum())

        # FA
        num_data['FA count'].append((group['Glial cell'] == 'fa').sum())
        num_data['FA count_g'].append(((group['Glial cell'] == 'fa') & (group['Colour'] == 'g')).sum())
        num_data['FA count_r'].append(((group['Glial cell'] == 'fa') & (group['Colour'] == 'r')).sum())
        num_data['FA count_le'].append(((group['Glial cell'] == 'fa') & (group['Hemisphere'] == 'l')).sum())
        num_data['FA count_ri'].append(((group['Glial cell'] == 'fa') & (group['Hemisphere'] == 'r')).sum())

        # WF
        num_data['WF count'].append((group['Glial cell'] == 'wf').sum())
        num_data['WF count_g'].append(((group['Glial cell'] == 'wf') & (group['Colour'] == 'g')).sum())
        num_data['WF count_r'].append(((group['Glial cell'] == 'wf') & (group['Colour'] == 'r')).sum())
        num_data['WF count_le'].append(((group['Glial cell'] == 'wf') & (group['Hemisphere'] == 'l')).sum())
        num_data['WF count_ri'].append(((group['Glial cell'] == 'wf') & (group['Hemisphere'] == 'r')).sum())

        # GL
        num_data['GL count'].append((group['Glial cell'] == 'gl').sum())
        num_data['GL count_g'].append(((group['Glial cell'] == 'gl') & (group['Colour'] == 'g')).sum())
        num_data['GL count_r'].append(((group['Glial cell'] == 'gl') & (group['Colour'] == 'r')).sum())
        num_data['GL count_le'].append(((group['Glial cell'] == 'gl') & (group['Hemisphere'] == 'l')).sum())
        num_data['GL count_ri'].append(((group['Glial cell'] == 'gl') & (group['Hemisphere'] == 'r')).sum())

        # OG
        num_data['OG count'].append((group['Glial cell'] == 'og').sum())
        num_data['OG count_g'].append(((group['Glial cell'] == 'og') & (group['Colour'] == 'g')).sum())
        num_data['OG count_r'].append(((group['Glial cell'] == 'og') & (group['Colour'] == 'r')).sum())
        num_data['OG count_le'].append(((group['Glial cell'] == 'og') & (group['Hemisphere'] == 'l')).sum())
        num_data['OG count_ri'].append(((group['Glial cell'] == 'og') & (group['Hemisphere'] == 'r')).sum())

        # OL
        num_data['OL count'].append((group['Glial cell'] == 'ol').sum())
        num_data['OL count_g'].append(((group['Glial cell'] == 'ol') & (group['Colour'] == 'g')).sum())
        num_data['OL count_r'].append(((group['Glial cell'] == 'ol') & (group['Colour'] == 'r')).sum())
        num_data['OL count_le'].append(((group['Glial cell'] == 'ol') & (group['Hemisphere'] == 'l')).sum())
        num_data['OL count_ri'].append(((group['Glial cell'] == 'ol') & (group['Hemisphere'] == 'r')).sum())

        # PR
        num_data['PR count'].append((group['Glial cell'] == 'pr').sum())
        num_data['PR count_g'].append(((group['Glial cell'] == 'pr') & (group['Colour'] == 'g')).sum())
        num_data['PR count_r'].append(((group['Glial cell'] == 'pr') & (group['Colour'] == 'r')).sum())
        num_data['PR count_le'].append(((group['Glial cell'] == 'pr') & (group['Hemisphere'] == 'l')).sum())
        num_data['PR count_ri'].append(((group['Glial cell'] == 'pr') & (group['Hemisphere'] == 'r')).sum())

        # Uncertain
        num_data['Uncertain count'].append((group['Uncertain'] == 'yes').sum())
        num_data['Uncertain count_g'].append(((group['Uncertain'] == 'yes') & (group['Colour'] == 'g')).sum())
        num_data['Uncertain count_r'].append(((group['Uncertain'] == 'yes') & (group['Colour'] == 'r')).sum())
        num_data['Uncertain count_le'].append(((group['Uncertain'] == 'yes') & (group['Hemisphere'] == 'l')).sum())
        num_data['Uncertain count_ri'].append(((group['Uncertain'] == 'yes') & (group['Hemisphere'] == 'r')).sum())

        # Midline (d or v)
        num_data['Midline count (dv)'].append(group['Midline cell'].isin(['d', 'v']).sum())
        num_data['Midline count (dv)_g'].append((group['Midline cell'].isin(['d', 'v']) & (group['Colour'] == 'g')).sum())
        num_data['Midline count (dv)_r'].append((group['Midline cell'].isin(['d', 'v']) & (group['Colour'] == 'r')).sum())
        num_data['Midline count (dv)_le'].append((group['Midline cell'].isin(['d', 'v']) & (group['Hemisphere'] == 'l')).sum())
        num_data['Midline count (dv)_ri'].append((group['Midline cell'].isin(['d', 'v']) & (group['Hemisphere'] == 'r')).sum())

        # Midline d
        num_data['Midline d count'].append((group['Midline cell'] == 'd').sum())
        num_data['Midline d count_g'].append(((group['Midline cell'] == 'd') & (group['Colour'] == 'g')).sum())
        num_data['Midline d count_r'].append(((group['Midline cell'] == 'd') & (group['Colour'] == 'r')).sum())
        num_data['Midline d count_le'].append(((group['Midline cell'] == 'd') & (group['Hemisphere'] == 'l')).sum())
        num_data['Midline d count_ri'].append(((group['Midline cell'] == 'd') & (group['Hemisphere'] == 'r')).sum())

        # Midline v
        num_data['Midline v count'].append((group['Midline cell'] == 'v').sum())
        num_data['Midline v count_g'].append(((group['Midline cell'] == 'v') & (group['Colour'] == 'g')).sum())
        num_data['Midline v count_r'].append(((group['Midline cell'] == 'v') & (group['Colour'] == 'r')).sum())
        num_data['Midline v count_le'].append(((group['Midline cell'] == 'v') & (group['Hemisphere'] == 'l')).sum())
        num_data['Midline v count_ri'].append(((group['Midline cell'] == 'v') & (group['Hemisphere'] == 'r')).sum())

        # Commissural (d or v)
        num_data['Commissural count_(dv)'].append(group['Commissural neuron'].isin(['d', 'v']).sum())
        num_data['Commissural count_(dv)_g'].append((group['Commissural neuron'].isin(['d', 'v']) & (group['Colour'] == 'g')).sum())
        num_data['Commissural count_(dv)_r'].append((group['Commissural neuron'].isin(['d', 'v']) & (group['Colour'] == 'r')).sum())
        num_data['Commissural count_(dv)_le'].append((group['Commissural neuron'].isin(['d', 'v']) & (group['Hemisphere'] == 'l')).sum())
        num_data['Commissural count_(dv)_ri'].append((group['Commissural neuron'].isin(['d', 'v']) & (group['Hemisphere'] == 'r')).sum())

        # Commissural d
        num_data['Commissural d count'].append((group['Commissural neuron'] == 'd').sum())
        num_data['Commissural d count_g'].append(((group['Commissural neuron'] == 'd') & (group['Colour'] == 'g')).sum())
        num_data['Commissural d count_r'].append(((group['Commissural neuron'] == 'd') & (group['Colour'] == 'r')).sum())
        num_data['Commissural d count_le'].append(((group['Commissural neuron'] == 'd') & (group['Hemisphere'] == 'l')).sum())
        num_data['Commissural d count_ri'].append(((group['Commissural neuron'] == 'd') & (group['Hemisphere'] == 'r')).sum())

        # Commissural v
        num_data['Commissural v count'].append((group['Commissural neuron'] == 'v').sum())
        num_data['Commissural v count_g'].append(((group['Commissural neuron'] == 'v') & (group['Colour'] == 'g')).sum())
        num_data['Commissural v count_r'].append(((group['Commissural neuron'] == 'v') & (group['Colour'] == 'r')).sum())
        num_data['Commissural v count_le'].append(((group['Commissural neuron'] == 'v') & (group['Hemisphere'] == 'l')).sum())
        num_data['Commissural v count_ri'].append(((group['Commissural neuron'] == 'v') & (group['Hemisphere'] == 'r')).sum())

        # MN
        num_data['MN count'].append((group['MN'] == 'yes').sum())
        num_data['MN count_g'].append(((group['MN'] == 'yes') & (group['Colour'] == 'g')).sum())
        num_data['MN count_r'].append(((group['MN'] == 'yes') & (group['Colour'] == 'r')).sum())
        num_data['MN count_le'].append(((group['MN'] == 'yes') & (group['Hemisphere'] == 'l')).sum())
        num_data['MN count_ri'].append(((group['MN'] == 'yes') & (group['Hemisphere'] == 'r')).sum())

    clone_size_df = pd.DataFrame(num_data)
    return clone_size_df

# Get the SC level (Ce, Th, Lu, Sa) based on the segment (str value s1, s2, s3, s4, s5, s6)
def getLevel(segment):
    # if segment == 's1' or segment == 's2':
    #     level = 'Ce'
    # elif segment == 's3' or segment == 's4':
    #     level = 'Th'
    # elif segment == 's5':
    #     level = 'Lu'
    # elif segment == 's6':
    #     level = 'Sa'
        
    if segment == 's1' or segment == 's2':
        level = 'limb'
    elif segment == 's3' or segment == 's4':
        level = 'thoracic'
    elif segment == 's5':
        level = 'limb'
    elif segment == 's6':
        level = 'sacral'
        
    return level    

# Input a group object from a dataframe that corresponds to one clone (Animal + Clone) and return the max pair-wise 
# Euclidean distance within the clone. 'positions' is a list of the column names corresponding to the positions (i.e. coordinates)
# you want to find the distance over. e.g. positions = ['Position X', 'Position Y'] if you want to find Euclidian distance in xy
# Excludes 1-cell clones
def intraCloneSpread(df, TM, positions):
    # Filter the dataframe for the given TM
    df_TM = df[(df['TM'] == TM)]
    grouped = df_TM.groupby(['Animal', 'Clone'])
    
    # Dictionary to store data
    spread_data = {
        'Clone': [],
        'TM': [],
        'Spread': [],
        'Size': [],
        'List of NNDs': [],
        'Cell_dist_avg': [],
        'Cell_dist_sd': []} # standard deviation of average cell-cell distances within a clone
    
    for name, group in grouped:
        position_array = group[positions].values  # numpy array of positions
        
        num_cells = len(group)  # Number of cells in the clone
        
        # Calculate intra-clone spread if clone size is greater than 1
        if num_cells > 1:
            distances = pdist(position_array, metric='euclidean')  # Pairwise distances between all cells in the clone
            max_distance = np.max(distances)  # Maximum pairwise distance
            neighbour_distances = np.diff(position_array, axis=0) # Distances between neighbouring cells in the clone.
            cell_dist_avg = np.mean(neighbour_distances) # Average cell-cell nearest-neighbour distance
            cell_dist_sd = np.std(neighbour_distances) # SD of average cell-cell nearest-neighbour distance
            
            # Store results in the dictionary
            clone_code = str(name).strip('()').strip("'").replace("', '", '_')
            spread_data['Clone'].append(clone_code)
            spread_data['TM'].append(TM)
            spread_data['Spread'].append(max_distance)
            spread_data['Size'].append(num_cells)
            spread_data['List of NNDs'].append([item for sublist in neighbour_distances for item in sublist]) # flatten the list of lists
            spread_data['Cell_dist_avg'].append(cell_dist_avg)
            spread_data['Cell_dist_sd'].append(cell_dist_sd)
    
    # Create a DataFrame from the dictionary
    intra_spread_df = pd.DataFrame(spread_data)
    
    # # Calculate intra-clone spread statistics
    # intra_spread_avg = size_spread_df['Spread'].mean()
    # intra_spread_sem = stats.sem(size_spread_df['Spread']) if not size_spread_df.empty else np.nan
    
    return intra_spread_df

# Input the dataframe containing all clones and output one containing the average inter-cell distance and variation in inter-cell distance per clone
# 'positions' is a list of the column names corresponding to the positions (i.e. coordinates)
# you want to find the distance over. e.g. positions = ['Position X', 'Position Y', 'Position Z'] if you want to find Euclidian distance in xyz
# Excludes 1-cell clones
def interCellDist(df, TM, positions):
    grouped = df.groupby(['Animal', 'Clone'])
    
    # Dictionary to store data
    inter_dist_data = {
        'Clone': [],
        'TM': [],
        'Avg Dist': [],
        'Var Dist': []}
    
    for name, group in grouped:
        position_array = group[positions].values # numpy array of positions
        num_cells = np.shape(group)[0]
        # Find the intra-clone spread
        if num_cells != 1: # If it's a 1-cell clone, exclude it
            distances = pdist(position_array, metric = 'euclidean') # Calculate the pairwise Euclidean distances
            inter_dist_data['Clone'].append(name)
            inter_dist_data['TM'].append(TM)
            inter_dist_data['Avg Dist'].append(sum(distances) / len(distances)) # Find the average inter-cell distance
            inter_dist_data['Var Dist'].append(np.var(distances, ddof = 1)) # Find the sample variance        

    # Make a dataframe containing the clone identifier, clone inter-cell distance mean and variation.
    inter_dist_df = pd.DataFrame(inter_dist_data)
    
    return inter_dist_df

# ToDo change data collection to dictionary instead of lists
def interCloneSpread(df, TM):
    # Subset for each animal then group by clone and find z upper and lower bounds for each clone.
    df_TM = df[(df['TM'] == TM)]
    print('\ndf_TM\n', df_TM)
    inter_spread_list = []
    animals = df_TM.Animal.unique() # Get a list of the unique animals in the dataframe
    print('\nanimals\n', animals)
    for animal in animals:
        subset = df_TM[(df_TM['Animal'] == animal)]
        print('\nsubset\n', subset)
        grouped = subset.groupby('Clone', sort = False) # set sort = False to avoid sorting groups in alphabetical order (want to maintain ascending z)
        if grouped.ngroups == 1: # If there is only one clone in the animal, skip the measurement
            continue
        else:
            print('yes')
            first_group = True
            for name, group in grouped:
                if first_group == True:
                    upper_z = group['Position Z'].max() # Set upper z bound of "previous" clone for the next iteration
                    first_group = False
                else: 
                    lower_z = group['Position Z'].min()# Lower z bound of current clone
                    inter_spread_list.append(lower_z - upper_z) # Inter-clone distance
                    upper_z = group['Position Z'].max() # Set upper z bound of "previous" clone for the next iteration
                    
    inter_spread_avg = np.mean(inter_spread_list)
    inter_spread_sem = stats.sem(inter_spread_list)
    print('\ninter_spread_list\n', inter_spread_list)
    
    return inter_spread_avg, inter_spread_sem, inter_spread_list # Return average clone spread for the inputted df

# Convert a long form dataframe that has binary (1 or 0) counts in some columns to a dataframe with proportions
# counts is the name of the column containing the binary counts. groupby is a list containing the names of the columns that you want to group the data by to find the proportions
# order is the order the rows should be in based on values of the ordered column
def binaryToProportion(df, counts, groupby, ordered_column, order):
    df_filtered = df[df[counts] == 1] # Filter rows where Count == 1   
    df_grouped = df_filtered.groupby(groupby).size().reset_index(name=counts) # Calculate the total number of counts per groupby columns to normalize the data
    df_grouped['Proportion'] = df_grouped.groupby(groupby[0])[counts].transform(lambda x: x / x.sum()) # Calculate the proportions (relative counts) for each groupby combination  
    df_grouped[ordered_column] = pd.Categorical(df_grouped[ordered_column], categories=order, ordered=True) # Convert the first column to a categorical type with the specific order
    df_grouped = df_grouped.sort_values(ordered_column).reset_index(drop=True) # Sort the DataFrame by the ordered_column
    
    return df_grouped

# Return the counts of dorsal (D) and ventral (V) cells in a clone
def DVRatio(df, TM):
    df_TM = df[(df['TM'] == TM)]
    grouped = df_TM.groupby(['Animal', 'Clone'])
    
    # Dictionary to store data
    dv_data = {
        'Clone': [],
        'Animal': [],
        'TM': [],
        'Dorsal Fraction': [],
        'Ventral Fraction': [],
        'DV': [],
        'Segment': [],
        'Level': []}
    
    for name, group in grouped:
        if '?' in group['D/V'].values: # If the d/v position is unknown (marked by a question mark), skip this clone. This might result in a smaller n for this time point compared to other graphs.
            continue
        clone_code = str(name).strip('()').strip("'").replace("', '", '_')
        dv_data['Clone'].append(str(clone_code)) # Store clone code
        dv_data['TM'].append(TM)
        dv_data['Animal'].append(group.iloc[0]['Animal'])
        segment = group.iloc[0]['Segment']
        dv_data['Segment'].append(segment)
        dv_data['Level'].append(getLevel(segment))
        if 'd' in group['D/V'].values:
            num_dorsal_cells = group['D/V'].value_counts()['d']
        else:
            num_dorsal_cells = 0
        if 'v' in group['D/V'].values:
            num_ventral_cells = group['D/V'].value_counts()['v']
        else:
            num_ventral_cells = 0
                    
        total_cells = num_dorsal_cells + num_ventral_cells
        dorsal_fraction = num_dorsal_cells / total_cells
        dv_data['Dorsal Fraction'].append(dorsal_fraction)
        ventral_fraction = num_ventral_cells / total_cells
        dv_data['Ventral Fraction'].append(ventral_fraction)
        
        if dorsal_fraction == 1:
            dv_data['DV'].append('Dorsal')
        elif ventral_fraction == 1:
            dv_data['DV'].append('Ventral')
        else:
            dv_data['DV'].append('Mixed')
     
    dv_data_df = pd.DataFrame(dv_data)
    
    num_clones = len(dv_data_df)
    num_animals = dv_data_df['Animal'].nunique() # Get a list of the unique animals in the dataframe
    
    return num_clones, num_animals, dv_data_df

# Return the frequency of MN-containing clones.
def MNFrequency(df, TM):
    df_TM = df[(df['TM'] == TM)]
    grouped = df_TM.groupby(['Animal', 'Clone'])
    
    # Dictionary to store data
    mn_data = {
        'Clone': [],
        'Animal': [],
        'TM': [],
        'MN fraction': [],
        'Non-MN fraction': [],
        'Clone size': [],
        'Segment': [],
        'Level': []}
    
    for name, group in grouped:
        if ~group['MN'].isin(['yes', 'no']).any(): # If the MN column is not annotated (contains something other than 'yes' or 'no'), skip this clone. This might result in a smaller n for this time point compared to other graphs.
            continue
        # if 'yes' in group['MN'].values:
        #     print('\n\n', name)
        #     print(TM)
        #     print('MN-containing')
        MN_fraction = (group['MN'] == 'yes').mean() # Get the fraction of 'yes' in the column
        mn_data['MN fraction'].append(MN_fraction)
        mn_data['Non-MN fraction'].append(1 - MN_fraction)
        mn_data['Clone size'].append(len(group['MN'])) # Get the number of cells in the clone (length of the clone dataframe)
        mn_data['Clone'].append(name) 
        mn_data['TM'].append(TM)
        mn_data['Animal'].append(group.iloc[0]['Animal'])
        segment = group.iloc[0]['Segment']
        mn_data['Segment'].append(segment)
        mn_data['Level'].append(getLevel(segment))
     
    mn_data_df = pd.DataFrame(mn_data)
    
    num_clones = len(mn_data_df)
    num_animals = mn_data_df['Animal'].nunique() # Get the number of the unique animals in the dataframe
    
    return num_clones, num_animals, mn_data_df

# Return the counts of clones of different types (symmetric, asymmetric, terminal neurogenic, small neurogenic, and unicellular)
# Set unicolor to False if you want to exclude unicolor clones
def cloneType(df, TM, unicolor = True):
    df_TM = df[(df['TM'] == TM)]
    grouped = df_TM.groupby(['Animal', 'Clone'])
    
    # Dictionary to store data
    type_data = {
        'Clone': [],
        'Animal': [],
        'TM': [],
        'Type': [],
        'Num_r': [],
        'Num_g': [],
        'Size': [],
        'Segment': [],
        'Level': []}
    
    for name, group in grouped:
        clone_code = str(name).strip('()').strip("'").replace("', '", '_')
        type_data['Clone'].append(clone_code) # Store clone code
        type_data['TM'].append(TM)
        type_data['Animal'].append(group.iloc[0]['Animal'])
        segment = group.iloc[0]['Segment']
        type_data['Segment'].append(segment)
        type_data['Level'].append(getLevel(segment))
        if 'g' in group['Colour'].values:
            num_g_cells = group['Colour'].value_counts()['g']
        else:
            num_g_cells = 0
        if 'r' in group['Colour'].values:
            num_r_cells = group['Colour'].value_counts()['r']
        else:
            num_r_cells = 0
            
        # Assign clone type
        if num_g_cells == 0 or num_r_cells == 0:
            type_data['Type'].append('Unicolor') # Store clone type 
        elif (num_g_cells >= 4 and num_r_cells < 4) or (num_r_cells >= 4 and num_g_cells < 4):
            type_data['Type'].append('Asymmetric') 
        elif num_g_cells == 1 and num_r_cells == 1:
            type_data['Type'].append('Terminal neurogenic') 
        elif num_g_cells >= 4 and num_r_cells >= 4:
            type_data['Type'].append('Symmetric') 
        else:
            type_data['Type'].append('Small neurogenic')  

        type_data['Num_g'].append(num_g_cells)
        type_data['Num_r'].append(num_r_cells)
        clone_size = num_g_cells + num_r_cells
        type_data['Size'].append(clone_size)
        
    type_data_df = pd.DataFrame(type_data) # Convert dictionary to df
    
    if unicolor == False:
        type_data_df = type_data_df[type_data_df['Type'] != 'Unicolor']
    
    num_clones = len(type_data_df)
    num_animals = type_data_df['Animal'].nunique() # Get a list of the unique animals in the dataframe
    
    return num_clones, num_animals, type_data_df

# Return the counts of clones with unilateral or bilateral cells. Can be 'uni' or 'bi'. Midline cells are considered ipsilateral
def cloneLat(df, TM):
    df_TM = df[df['TM'] == TM]
    grouped = df_TM.groupby(['Animal', 'Clone'])
    
    data = []  # A list of dictionaries to hold the data for each clone

    for name, group in grouped:
        hem_values = group['Hemisphere'].values
        if '' in hem_values:
            continue  # Skip clones with unannotated hemispheres
        cell_num = group.shape[0]
        if (group['Hemisphere'] == 'na').all():
            continue # Skip clones with only midline cells
        if cell_num <= 1:
            continue  # Skip one-cell clones
        
        clone_code = str(name).strip('()').strip("'").replace("', '", '_')
        clone_data = {'Clone': clone_code, 'Size': cell_num}  # Store size and clone code info
        
        # Determine laterality
        if ('l' in hem_values) != ('r' in hem_values):  # True if l or r, but not both
            clone_data['Laterality'] = 'uni'
        elif ('l' in hem_values) == ('r' in hem_values):  # True if both l and r are present
            clone_data['Laterality'] = 'bi'
            # Record the size of the majority and minority clone fragments
            num_l = np.sum(hem_values == 'l') 
            num_r = np.sum(hem_values == 'r')
            # Determine the colour representation in the two clone fragments
            l_cells = group[(group['Hemisphere'] == 'l')]
            if ('r' in l_cells['Colour'].values) and ('g' in l_cells['Colour'].values):
                l_colour = 'mix'
            else:
                l_colour = 'uni'
            r_cells = group[(group['Hemisphere'] == 'r')]
            if ('r' in r_cells['Colour'].values) and ('g' in r_cells['Colour'].values):
                r_colour = 'mix'
            else:
                r_colour = 'uni'
            if num_l > num_r: # If majority fragment is on left side
                clone_data['Majority fragment'] = int(num_l)
                clone_data['Minority fragment'] = int(num_r)
                clone_data['Majority fragment colour'] = l_colour
                clone_data['Minority fragment colour'] = r_colour
            # Minority fragment
            else: # If majority fragment is on right side
                clone_data['Majority fragment'] = int(num_r)
                clone_data['Minority fragment'] = int(num_l)
                clone_data['Majority fragment colour'] = r_colour
                clone_data['Minority fragment colour'] = l_colour
                
        else:
            raise ValueError("Unrecognized entry in Hemisphere column.")
        
        # Add any additional information from the group
        segment = group.iloc[0]['Segment'] # Assuming 'Segment' is constant in the group
        clone_data['Segment'] = segment  
        clone_data['Level'] = getLevel(segment)
        clone_data['TM'] = TM
        
        # Append the data dictionary to the list
        data.append(clone_data)

    # Create DataFrame from the list of dictionaries
    df_lat_size = pd.DataFrame(data)
    
    # Proportions for the overall summary DataFrame
    num_clones = len(df_lat_size)
    num_animals = len(set(df_lat_size['Clone'].apply(lambda x: x[0])))  # Extract 'Animal' from 'Clone'
    if num_clones > 0:  # Avoid division by zero
        uni_prop = sum(d['Laterality'] == 'uni' for d in data) / num_clones
        bi_prop = sum(d['Laterality'] == 'bi' for d in data) / num_clones
    else:
        uni_prop, bi_prop, mid_prop = 0, 0, 0

    df_proportions = pd.DataFrame({
        'TM': [TM] * 2,
        'Type': ['Unilateral', 'Bilateral'],
        'Proportion': [uni_prop, bi_prop]
    })

    return num_clones, num_animals, df_proportions, df_lat_size

def extract_animal_id(clone_str: str) -> str:
    # Gets Animal from Clone ID: e.g. S107_a7_c3 returns S107_a7
    parts = str(clone_str).split("_")
    if len(parts) >= 3:
        return "_".join(parts[:2])

def SECalculation(sample_prop, n):
    # Standard error of a proportion: sqrt(p(1-p)/n)
    se = np.sqrt(sample_prop * (1.0 - sample_prop) / n)
    return se

# Get the proportions of different types of clone in a dataframe by grouping the clones by the criteria of choice. 
# df is the input dataframe. If it doesn't have a column named 'TM' or you don't want to sort the rows by ascending TM timepoint, set sort_TM to False
# "criteria" is a list of the names of the dataframe columns (one or several) that you want to group by, in addition to the column with the clone type
# prop_of is a list containing the name of the column with the feature you want to get the proportion of (e.g. clone 'Type')
def cloneTypeProportions(df, criteria, prop_of, sort_TM = True, sort_level = False):
    counts_df = df.groupby(criteria + prop_of).size().reset_index(name = 'Count') # Count the occurrences of each type of clone for each TM timepoint. Concatenate the criteria and prop_of lists to input them to groupby
    df['Animal'] = df['Clone'].apply(extract_animal_id)
    animals_df = df.groupby(criteria + prop_of)['Animal'].nunique().reset_index(name='Animals') # Get the number of unique animals per criteria + prop_by group
    counts_df = counts_df.merge(animals_df, on=criteria + prop_of, how='left') # Add a column with the number of unique animals to the counts_df
    
    # Get the proportions of the diff clone types
    total_counts = df.groupby(criteria).size().reset_index(name = 'Total')
    counts_df = counts_df.merge(total_counts, on = criteria)
    counts_df['Proportions'] = counts_df['Count'] / counts_df['Total']
    
    # Final df
    prop_df = counts_df[criteria + prop_of + ['Proportions']]
    print('\nprop_df\n', prop_df)
    
    # Standard error of a proportion: sqrt(p(1-p)/n)
    counts_df['SE'] = SECalculation(counts_df['Proportions'], counts_df['Total'])
    print('\ncount_df\n', counts_df)
    
    if sort_TM == True:
        # Sort the rows in the correct TM order to avoid issues during plotting where xticks don't match the bars
        custom_order = ['E8.5', 'E9.5', 'E10.5', 'E11.5']
        prop_df['TM'] = pd.Categorical(prop_df['TM'], categories = custom_order, ordered = True)
        prop_df = prop_df.sort_values(by = 'TM')
        counts_df['TM'] = pd.Categorical(counts_df['TM'], categories = custom_order, ordered = True) # also sort counts_df for neatness
        counts_df = counts_df.sort_values(by = 'TM')
        
    if sort_level == True:
        custom_order = ['limb', 'thoracic', 'sacral']
        prop_df['Level'] = pd.Categorical(prop_df['Level'], categories = custom_order, ordered = True)
        prop_df = prop_df.sort_values(by = 'Level')
    
    return counts_df

# Return the gene expression make-up of clones
# gene is the antibody name used as a column in Spots Data as a string e.g. 'Lhx5'
def coExpression(df, TM, gene):
    df_TM = df[(df['TM'] == TM)]
    gene_subset = df_TM[df_TM[gene].notna()]
    grouped = gene_subset.groupby(['Animal', 'Clone'])
    
    # Dictionary to store data
    expression_data = {
        'Clone': [],
        'Animal': [],
        'Segment': [],
        'Level': [],
        'TM': [],
        'Expression': [],
        'Pos cell colour comp': [],
        'Neg cell colour comp': []
        }
    
    for name, group in grouped:
        clone_code = str(name).strip('()').strip("'").replace("', '", '_')
        expression_data['Clone'].append(clone_code) # Store clone code
        expression_data['TM'].append(TM)
        expression_data['Animal'].append(group.iloc[0]['Animal'])
        segment = group.iloc[0]['Segment']
        expression_data['Segment'].append(segment)
        expression_data['Level'].append(getLevel(segment))
        
        # Find expression information. Convert the column to a list then convert the list to a set
        express_list = group[gene].values
        express_list = set(express_list)
        if (express_list == {'yes', 'no'}) or (express_list == {'no', 'yes'}):
            expression_data['Expression'].append('mixed express') # mixed gene expression
            # Get information about the colour-composition of the mixed expression clones
            yes_gene_col = group[group[gene] == 'yes']['Colour'].tolist() # List of colours of cells expressing the gene
            yes_gene_col = set(yes_gene_col) # Get the set
            no_gene_col = group[group[gene] == 'no']['Colour'].tolist() # List of colours of cells not expressing the gene
            no_gene_col = set(no_gene_col) # Get the set
            if (yes_gene_col == {'g', 'r'}) or (yes_gene_col == {'r', 'g'}):
                yes_gene_col_comp = 'mixed'
            elif yes_gene_col == {'g'} or yes_gene_col == {'r'}:
                yes_gene_col_comp = 'uni'
            if (no_gene_col == {'g', 'r'}) or (no_gene_col == {'r', 'g'}):
                no_gene_col_comp = 'mixed'
            elif no_gene_col == {'g'} or no_gene_col == {'r'}:
                no_gene_col_comp = 'uni'
            # Record if the gene-expressing cells (pos) are of a single colour or mixed. Do the same with the gene-non-expressing cells (neg)
            expression_data['Pos cell colour comp'].append(yes_gene_col_comp)
            expression_data['Neg cell colour comp'].append(no_gene_col_comp)
        elif express_list == {'yes'}:
            expression_data['Expression'].append('coexpress-only') # all cells co-express the gene
            # Record NA in the colour composition column 
            expression_data['Pos cell colour comp'].append('na')
            expression_data['Neg cell colour comp'].append('na')
        elif express_list == {'no'}:
            expression_data['Expression'].append('no coexpress') # no cells express the gene
            # Record NA in the colour composition column 
            expression_data['Pos cell colour comp'].append('na')
            expression_data['Neg cell colour comp'].append('na')
        else:
            print('Unrecognized input in co-expression column.')
    
    expression_data_df = pd.DataFrame(expression_data) # Convert dictionary to df
    
    num_clones = len(expression_data_df)
    num_animals = expression_data_df['Animal'].nunique() # Get a list of the unique animals in the dataframe
    
    return num_clones, num_animals, expression_data_df
