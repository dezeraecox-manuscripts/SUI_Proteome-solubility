import os
import pandas as pd
import numpy as np
from loguru import logger

import matplotlib.pyplot as plt
import seaborn as sns

from GEN_Utils.PlotUtils import FileHandling
from matplotlib.colors import ListedColormap

from analysis_scripts.utils import image_maker

logger.info("Import OK")

# Print all lone variables during execution
# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = 'all'
# # Set plotting backgrounds to white
# matplotlib.rcParams.update(_VSCode_defaultMatplotlib_Params)
# matplotlib.rcParams.update({'figure.facecolor': (1, 1, 1, 1)})

input_path = 'Python_results/FUS/pixel_calculations/GFP_pixel_summary.xlsx'
output_path = 'Python_results/FUS/ROI_masks/'

if not os.path.exists(output_path):
    os.mkdir(output_path)

# Read in compiled pixel summary
raw_data = pd.read_csv(input_path)
pixels = raw_data.copy()
pixels.shape
pixels.head()

# Generate binary masked arrays for all images cytoplasm
for image_name, df in pixels.groupby('image_name'):
    cyto_pixels = df[df['pixel_type'] == 'cyto']
    array_image = image_maker(cyto_pixels)
    plt.imshow(array_image, cmap='gray')
    plt.axis('off')
    np.save(f'{output_path}{image_name}_cyto_array.npy', array_image)

# Generate binary masked arrays for all images nucleus
for image_name, df in pixels.groupby('image_name'):
    nuc_pixels = df[df['pixel_type'] == 'nucleus']
    array_image = image_maker(nuc_pixels)
    plt.imshow(array_image, cmap='gray')
    plt.axis('off')
    np.save(f'{output_path}{image_name}_nucl_array.npy', array_image)

# Generate binary masked arrays for all images all pixels
for image_name, df in pixels.groupby('image_name'):
    array_image = image_maker(df)
    plt.imshow(array_image, cmap='gray')
    plt.axis('off')
    np.save(f'{output_path}{image_name}_all_array.npy', array_image)
