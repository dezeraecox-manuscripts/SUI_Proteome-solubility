import os
import collections
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

input_path = 'cellprofiler_results/spots/'
output_path = 'Python_results/FUS/spot_to_ROI/'

if not os.path.exists(output_path):
    os.mkdir(output_path)

# Collect list of images/folders needed
# image 1 has large horrible spot in - exclude
list_25Q = [str(x) for x in np.arange(2, 11)]
# excluded by AO during ROI stage
list_97Q_i = [str(x) for x in np.arange(1, 11)]
list_97Q_i.pop(3)
list_97Q_ni = [str(x) for x in np.arange(1, 11)]
mutants = {'25Q': list_25Q, '97Q_i': list_97Q_i, '97Q_ni': list_97Q_ni}

# Generate default dict to store spot-based arrays
def tree():
    return collections.defaultdict(tree)

spot_dict = tree()

# Collect all raw filteredspots into list

for mutant, images in mutants.items():
    for image_number in images:

        # construct path to each image folder
        image_folder = f'{input_path}{mutant}/{image_number}/'
        # Generate list of the spot files in that folder
        file_list = [filename for filename in os.listdir(image_folder)]
        spots = []
        for filename in file_list:
            spot_number = (os.path.splitext(filename)[0]).split('_')[-1]
            # For every file in that folder, read in the tiff to nparray and store in list
            spots.append(plt.imread(f'{image_folder}{filename}'))
        spot_image = sum(spots)
        spot_dict[mutant][image_number]['spot_array'] = spot_image
        plt.imshow(spot_image)

# Generate compiled array for each image, save as png
for mutant, images in mutants.items():
    for image_number in images:
        spot_image = spot_dict[mutant][image_number]['spot_array']
        fig, ax = plt.subplots()
        plt.imshow(spot_image)
        plt.axis('off')
        # plt.savefig(f'{output_path}{mutant}_{image_number}_rawspots.png')

        # Save image spot array to np format to read into single cell plotting
        np.save(f'{output_path}{mutant}_{image_number}_spot_array.npy', spot_image)
        
