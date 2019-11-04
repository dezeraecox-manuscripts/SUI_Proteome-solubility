import os
import numpy as np
import pandas as pd
import math

from statistics import mode #, geometric_mean
from loguru import logger
from GEN_Utils import FileHandling
from GEN_Utils.CalcUtils import sorted_nicely
from analysis_scripts.utils import image_processor

# Print all lone variables during execution
# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = 'all'
# # Set plotting backgrounds to white
# matplotlib.rcParams.update(_VSCode_defaultMatplotlib_Params)
# matplotlib.rcParams.update({'figure.facecolor': (1,1,1,1)})

def geomean(val_list):
    return math.exp(math.fsum(math.log(x) for x in val_list) / len(val_list))

logger.info("Import OK")

# Set input and output folders
input_folder = 'imageJ_results/EIF4A3/Results/'
output_path = 'python_results/EIF4A3/pixel_calculations/'
# Set threshold for excluding pixels that are part of the inclusions
mCherry_threshold = 3000

# Check output folder exists - if not, create it
if not os.path.exists(output_path):
    os.mkdir(output_path)

# collect list of files in the input_folder that have .csv extension with pixels
file_list = [filename for filename in os.listdir(input_folder) if 'pixels.csv' in filename]
logger.info(f'The following files were detected from {input_folder}:\n {file_list}')

# Generate list of images with results in that folder
image_names = sorted_nicely(list(set([("_".join(name.split("_")[0:2])) for name in file_list])))
logger.info(f'The following image names were detected:\n {image_names}')

channel_names = ['GFP', 'mCherry']

# Clean ROI pixels to simple dataframe for each image
pixel_dict = {}
summary_dict = {}
for image_name in image_names:
    image = f'{image_name}_GFP'
    cleaned_pixels = image_processor(input_folder, image, output_path,  save_file=False)
    pixel_dict[image_name] = cleaned_pixels

logger.info('All images processed successfully.')

# Add mCherry information by mapping pixel locations
mCherry_dict = {}
for image, dictionary_1 in pixel_dict.items():
    cell_input = f'{input_folder}{image}_inc_mCherry_cells_pixels.csv'
    cell_raw = pd.read_csv(cell_input)
    mCherry_cells = cell_raw.copy()
    # Rename columns and drop previous column names
    new_col_list = list(cell_raw.T.reset_index()['index'].str.split('_').str[1:3])
    mCherry_cells.columns = ['_'.join(x) for x in new_col_list]
    # Add x, y position for mapping
    mCherry_cells['X,Y'] = list(zip(mCherry_cells.X_pos, mCherry_cells.Y_pos))
    mCherry_dict[image] = mCherry_cells

    # Generate mapper to map location to mCherry intensity
    mCherry_mapper = dict(zip(mCherry_cells['X,Y'], mCherry_cells['inc_mCherry']))

    # Map these pixels onto the dfs stored in pixel_dict
    for pixel_type, df in dictionary_1.items():
        df.columns = ['X_pos',	'Y_pos', 'GFP_intensity', 'ROI_name', 'image_name', 'X,Y']
        df['mCherry_intensity'] = df['X,Y'].map(mCherry_mapper)

# Generate compiled df with all pixels from nuc/cyto types for all images
types = ['nucleus', 'cyto']
compiled = []
for image, dictionary_1 in pixel_dict.items():
    for pixel_type, df in dictionary_1.items():
        if pixel_type in types:
            df['pixel_type'] = pixel_type
            compiled.append(df)
            df.shape
compiled = pd.concat(compiled)
compiled.head()
compiled.shape
# Fix some of the sample-specific details
compiled['mutant'], compiled['replicate'], _ = compiled['image_name'].str.split(
    '_').str
inclusion_list_25Q = [0 for pixel in compiled[compiled['mutant'] == '25Q']['mutant']]
inclusion_list_97Q = [0 if 'n' in x else 1 for x in compiled[compiled['mutant'] == '97Q']['replicate']]
inclusion_list_UT = [0 for pixel in compiled[compiled['mutant'] == 'UT']['mutant']]
compiled['inclusion'] = inclusion_list_25Q + inclusion_list_97Q + inclusion_list_UT
compiled['replicate'] = compiled['replicate'].str.extract('(\d+)', expand=False)

# save to csv
compiled.to_csv(output_path+f'GFP_pixel_summary.xlsx')


# Filter pixels from inclusions by thresholding any pixels above threshold in inclusion channel, then calculate stats of remaining pixels
filtered_pixels = pixel_dict.copy()
types = ['nucleus', 'cyto']
cell_summary = pd.DataFrame(index=['image_name', 'pixel_type', 'mean', 'median', 'std', 'mode', 'geomean'])
for image, dictionary_1 in filtered_pixels.items():
    for pixel_type, df in dictionary_1.items():
        if pixel_type in types:
            logger.info(f'Processing {pixel_type} for {image}.')
            logger.info(f'Original shape: {df.shape}')
            df = df[df['mCherry_intensity'] < mCherry_threshold]
            logger.info(f'New shape: {df.shape}')
            # calculate relevant info
            stats = [image, pixel_type, df['GFP_intensity'].mean(), df['GFP_intensity'].median(), df['GFP_intensity'].std(), df['GFP_intensity'].mode()[0], geomean(list(df['GFP_intensity']))]
            cell_summary[f'{image}_{pixel_type}'] = stats


# Clean up stats df
cell_summary = cell_summary.T.reset_index(drop=True)
# Add column for the mutant, replicate
cell_summary['mutant'], cell_summary['replicate'] = cell_summary['image_name'].str.split('_').str

# Save summary to excel
FileHandling.df_to_excel(output_path+f'GFP_stats_summary.xlsx',
                         sheetnames=['GFP_summary'], data_frames=[cell_summary])
