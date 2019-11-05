import os
import pandas as pd
import numpy as np
from loguru import logger

import matplotlib.pyplot as plt
import seaborn as sns

from GEN_Utils.PlotUtils import FileHandling
from matplotlib.colors import ListedColormap

logger.info("Import OK")

# Print all lone variables during execution
# from IPython.core.interactiveshell import InteractiveShell
# InteractiveShell.ast_node_interactivity = 'all'
# # Set plotting backgrounds to white
# matplotlib.rcParams.update(_VSCode_defaultMatplotlib_Params)
# matplotlib.rcParams.update({'figure.facecolor': (1, 1, 1, 1)})


def image_maker(df, channel='GFP_intensity'):
    array_im = np.zeros((512, 512))
    pixel_map = list(zip(df['X_pos'], df['Y_pos'], df[channel]))
    for pixel in pixel_map:
        x, y, val = pixel
        array_im[x, y] = val
    # flip to maintain correct orientation
    new_arr = np.rot90(array_im, 3)
    final_arr = np.fliplr(new_arr)
    return final_arr

input_path = 'Python_results/FUS/pixel_calculations/GFP_pixel_summary.xlsx'
output_path = 'Python_results/FUS/plot_single_cells/'

if not os.path.exists(output_path):
    os.mkdir(output_path)

# Read in compiled pixel summary
raw_data = pd.read_csv(input_path)
pixels = raw_data.copy()
pixels.shape

pixels.head()

# read in spot_arrays from spot_to_ROI
# Generate default dict to store spot-based arrays
input_spot = 'Python_results/FUS/spot_to_ROI/'
file_list = [filename for filename in os.listdir(
    input_spot)]

spot_dict = {}
for filename in file_list:
    filename
    mutant = filename.split('_')[0]
    if mutant == '97Q':
        inc, image_number = filename.split('_')[1:3]
        inc = [1 if inc == 'i' else 0][0]
        spot_dict[(mutant, int(image_number), inc)] = np.load(
            f'{input_spot}{filename}')
    elif mutant == '25Q':
        image_number = filename.split('_')[1]
        spot_dict[(mutant, int(image_number), 0)] = np.load(
            f'{input_spot}{filename}')
# convert to normal dict to only be able to access keys that have spots associated


# To plot representative images
image_name = '97Q_7_GFP'
df = pixels[pixels['image_name'] == image_name]
image_array = image_maker(df)
spot_array = spot_dict[('97Q', 7, 1)]
threshold = np.where(spot_array == 0, np.nan, 1)
# Plot whole cell
plt.imshow(image_array, cmap='gray')
plt.imshow(threshold, cmap=ListedColormap('#e6000b'))
plt.axis('off')
plt.savefig(f'{output_path}{image_name}_thresholded.png')
plt.show()


# to plot all images

for image_info, df in pixels.groupby(['mutant', 'replicate', 'inclusion']):
    mutant, image_number, inclusion = image_info
    image_info
    set(df['image_name'])
    image_array = image_maker(df)
    try:
        spot_array = spot_dict[image_info]
        threshold = np.where(spot_array == 0, np.nan, 1)
        # Plot whole cell
        plt.imshow(image_array, cmap='gray')
        plt.imshow(threshold, cmap=ListedColormap('#e6000b'))
        plt.axis('off')
        plt.savefig(f'{output_path}{mutant}_{image_number}_{inclusion}_thresholded.png')
        plt.show()
    except:
        logger.info(f'No spot array found for {image_info}')

# choosing colour map
# sns.palplot(sns.color_palette('hot')[2:3])
