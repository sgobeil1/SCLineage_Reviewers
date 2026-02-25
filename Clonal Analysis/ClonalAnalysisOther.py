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

def plotExpression(df, rg_clones_df, gene):
    num_clones_E8, num_animals_E8, express_data_df_E8 = coExpression(df, 'E8.5', gene)
    num_clones_E9, num_animals_E9, express_df_E9 = coExpression(df, 'E9.5', gene)
    num_clones_E10, num_animals_E10, express_df_E10 = coExpression(df, 'E10.5', gene)
    num_clones_E11, num_animals_E11, express_df_E11 = coExpression(df, 'E11.5', gene)

    clone_counts = {'Num E8.5 Clones': num_clones_E8,
                    'Num E9.5 Clones': num_clones_E9,
                    'Num E10.5 Clones': num_clones_E10,
                    'Num E11.5 Clones': num_clones_E11}
    animal_counts = {'Num E8.5 Animals': num_animals_E8,
                     'Num E9.5 Animals': num_animals_E9,
                    'Num E10.5 Animals': num_animals_E10,
                    'Num E11.5 Animals': num_animals_E11}
    
    express_data_df = pd.concat([express_data_df_E8, express_df_E9, express_df_E10, express_df_E11])
    print('\nexpress_data_df', express_data_df)
    column_name = gene + '_Expression'
    express_data_df_subset = express_data_df[['Clone','Expression']]
    express_data_df_subset = express_data_df_subset.rename(columns = {'Expression': column_name})
    rg_clones_df = pd.merge(rg_clones_df, express_data_df_subset, on = 'Clone', how = 'left') # merge, keeping matching and non-matching clone codes

    # Make a dataframe with the relative proportions of different clone types at different TM timepoints.
    proportions = cloneTypeProportions(express_data_df, ['TM'], ['Expression'])
           
    x_axis_ticks = ('E9.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E9.5 Clones'], animal_counts['Num E9.5 Animals']),
                    'E10.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E10.5 Clones'], animal_counts['Num E10.5 Animals']),
                    'E11.5\n$\it{n=%s}$\n$\it{a=%d}$' % (clone_counts['Num E11.5 Clones'], animal_counts['Num E11.5 Animals']))
    
    # Re-organize the order of the TM timepoints - important for snsStackedBarPlot
    tm_order = ['E9.5', 'E10.5', 'E11.5']
    proportions['TM'] = pd.Categorical(proportions['TM'], categories=tm_order, ordered=True) # Convert 'TM' column to a categorical type with the desired order
    proportions = proportions.sort_values('TM').reset_index(drop=True) # Sort the DataFrame by the 'TM' column

    # Collect information about the co-expressing clones that will be saved in a csv file
    express_data_df_mixed = express_data_df[express_data_df['Expression'] == 'mixed express']
    num_mixed = express_data_df_mixed.shape[0] # Number of mixed expression clones
    both_mixed = (express_data_df[(express_data_df['Pos cell colour comp'] == 'mixed') & (express_data_df['Neg cell colour comp'] == 'mixed')].shape[0]) / num_mixed
    uni_mixed = (express_data_df[(express_data_df['Pos cell colour comp'] == 'uni') & (express_data_df['Neg cell colour comp'] == 'mixed')].shape[0]) / num_mixed
    mixed_uni = (express_data_df[(express_data_df['Pos cell colour comp'] == 'mixed') & (express_data_df['Neg cell colour comp'] == 'uni')].shape[0]) / num_mixed
    both_uni = (express_data_df[(express_data_df['Pos cell colour comp'] == 'uni') & (express_data_df['Neg cell colour comp'] == 'uni')].shape[0]) / num_mixed
    mixed_express_info = {'The proportion of mixed-expression clones where both expressing and non-expressing cells are two colours is': both_mixed,
                          'The proportion of mixed-expression clones where the expressing cells are unicolour and the non-expressing cells are two-colour is': uni_mixed,
                          'The proportion of mixed-expression clones where the expressing cells are two-colour and the non-expressing cells are unicolour is': mixed_uni,
                          'The proportion of mixed-expression clones that are unicolour clones is': both_uni}
    mixed_express_info = pd.DataFrame(list(mixed_express_info.items()), columns = ['Measure', 'Fraction'])
    return proportions, x_axis_ticks, rg_clones_df, mixed_express_info
