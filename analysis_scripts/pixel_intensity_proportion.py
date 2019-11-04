import os
import pandas as pd
import numpy as np
from loguru import logger

import matplotlib.pyplot as plt
import seaborn as sns

from GEN_Utils.PlotUtils import FileHandling
from GEN_Utils.CalcUtils import sorted_nicely

# Print all lone variables during execution
# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = 'all'
# # Set plotting backgrounds to white
# matplotlib.rcParams.update(_VSCode_defaultMatplotlib_Params)
# matplotlib.rcParams.update({'figure.facecolor': (1, 1, 1, 1)})

logger.info("Import OK")

input_path = 'python_results/EIF4A3/pixel_calculations/GFP_pixel_summary.xlsx'
output_path = 'python_results/EIF4A3/pixel_intensity_proportion/'

if not os.path.exists(output_path):
    os.mkdir(output_path)

# Read in compiled pixel summary
raw_data = pd.read_csv(input_path)
pixels = raw_data.copy()
pixels.shape
pixels.sample(15)

# Threshold values
threshold = 3000

# Collect only nuclear pixels, then generate dfs for above and below threshold
# thresholded = pixels[pixels['pixel_type'] == 'nucleus'].copy()
thresholded = pixels.copy()
above_thresh = thresholded[thresholded['GFP_intensity'] > threshold]
below_thresh = thresholded[thresholded['GFP_intensity'] <= threshold]
# Check number pixels in each - expected majority in below!
above_thresh.shape
below_thresh.shape

# Group according to info columns, collect count of pixels in each category
grouping_columns = ['mutant', 'replicate', 'inclusion', 'pixel_type']
above_count = above_thresh.groupby(grouping_columns).count()['GFP_intensity'].reset_index()
above_count.rename(columns={'GFP_intensity': 'above_count'}, inplace=True)
below_count = below_thresh.groupby(grouping_columns).count()['GFP_intensity'].reset_index()
below_count.rename(columns={'GFP_intensity': 'below_count'}, inplace=True)

proportion = pd.merge(above_count, below_count, on=grouping_columns, how='outer')
proportion['proportion_above'] = proportion['above_count'] / proportion['below_count']

# Save summary to excel
FileHandling.df_to_excel(output_path+f'thresholded_summary.xlsx',
                         sheetnames=[f'thresholded_{threshold}'], data_frames=[proportion])

# quick look at what the distribution looks like!
sns.swarmplot(x="mutant", y="proportion_above", hue="pixel_type", data=proportion, palette="Set2", dodge=True)




