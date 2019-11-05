import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from functools import reduce

from loguru import logger
from GEN_Utils import FileHandling

logger.info('Import OK')

# Print all lone variables during execution
# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = 'all'
# # Set plotting backgrounds to white
# matplotlib.rcParams.update(_VSCode_defaultMatplotlib_Params)
# matplotlib.rcParams.update({'figure.facecolor': (1,1,1,1)})

input_folder = 'cellprofiler_results/quantification/'
output_folder = 'Python_results/FUS/spot_cleanup/'

if not os.path.exists(output_folder):
    os.mkdir(output_folder)


list_25Q = [str(x) for x in np.arange(2, 11)] # image 1 has large horrible spot in - exclude
list_97Q_i = [str(x) for x in np.arange(1, 11)] # excluded by AO during ROI stage
list_97Q_i.pop(3)
list_97Q_ni = [str(x) for x in np.arange(1, 11)]

mutants = {'25Q': list_25Q, '97Q_i': list_97Q_i, '97Q_ni': list_97Q_ni}

# Collect all raw filteredspots into list
raw_data_collected = []
for mutant, images in mutants.items():
    for image_number in images:
        image_number
        input_path = f'{input_folder}{mutant}/{image_number}/filtered_spots.csv'
        raw_data = pd.read_csv(input_path)
        raw_data['mutant'] = mutant
        raw_data['image_number'] = image_number

        raw_data_collected.append(raw_data)

# concat into single summary dataframe
compiled_raw_data = pd.concat(raw_data_collected)
# Save to csv
compiled_raw_data.to_csv(f'{output_folder}compiled_raw_spots.csv')

compiled_raw_data.columns.tolist()
# create copy, collect columns of interest
cols_of_interest = ['AreaShape_Area', 'Intensity_MeanIntensity_GFP', 'Intensity_MinIntensity_GFP',  'Location_Center_X', 'Location_Center_Y',  'Number_Object_Number', 'mutant', 'image_number']
new_col_names = ['area', 'mean_intensity', 'min_intensity',  'center_X',
                 'center_Y',  'object_number', 'mutant', 'image_number']
compiled_spots = compiled_raw_data[cols_of_interest]
compiled_spots.columns = new_col_names

# Save to csv
compiled_spots.to_csv(f'{output_folder}compiled_spots.csv')


# Calculate interesting bits
## mean of area and mean_intensity
spot_calculations = compiled_spots.groupby(['mutant', 'image_number']).mean()[
    ['area', 'mean_intensity']].reset_index()

## count per cell
spot_count = compiled_spots.groupby(['mutant', 'image_number']).count()['object_number'].reset_index()
spot_min = compiled_spots.groupby(['mutant', 'image_number']).min()['min_intensity'].reset_index()

# concat into final summary df
spot_summary = reduce(lambda df1, df2: pd.merge(df1, df2, on=['mutant', 'image_number'], how='outer'), [spot_calculations, spot_count, spot_min])
spot_summary.rename(columns={'object_number': 'count'}, inplace=True)

FileHandling.df_to_excel(f'{output_folder}spot_summary.xlsx',
                         sheetnames=['spot_summary'], data_frames=[spot_summary])
