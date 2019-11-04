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
# matplotlib.rcParams.update({'figure.facecolor': (1,1,1,1)})

input_path = 'python_results/EIF4A3/pixel_calculations/GFP_pixel_summary.xlsx'
output_path = 'python_results/EIF4A3/single_cell_plotting/'

if not os.path.exists(output_path):
    os.mkdir(output_path)


def single_cell_plotter(image_name, output_path, threshold=3000, cmap='Reds_r'):
    # plot single cell nuclear and cyto pixels for Eif4a3 according to intensity above/below
    data = thresholded[thresholded['image_name'] == image_name]
    nuc = data[data['pixel_type'] == 'nucleus']
    cyto = data[data['pixel_type'] == 'cyto']

    im_nuc = image_maker(nuc)
    im_cyto = image_maker(cyto)

    # Plot nuc, cytoplasm separate
    cyto_plot = plt.imshow(
        np.where(im_cyto == 0, np.nan, im_cyto), cmap='Blues')
    nuc_plot = plt.imshow(np.where(im_nuc == 0, np.nan, im_nuc), cmap='Reds')
    plt.savefig(f'{output_path}{image_name}_NvC.png')
    # plt.show()

    # Plot whole cell
    plt.imshow(image_maker(data), cmap='gray')
    plt.savefig(f'{output_path}{image_name}_wholecell.png')
    # plt.show()

    # Plot whole cell with threshold
    threshold = np.where(im_nuc < threshold, np.nan, 1)
    plt.imshow(image_maker(data), cmap='gray')
    plt.imshow(threshold, cmap=cmap)
    plt.axis('off')
    plt.savefig(f'{output_path}{image_name}_thresholded.png')
    plt.show()


# Read in compiled pixel summary
raw_data = pd.read_csv(input_path)
pixels = raw_data.copy()
pixels.shape

pixels.head()

# Generate copy, add colour column according to thresholding
threshold = 3000
thresholded = pixels.copy()

for image_name in list(set(thresholded['image_name'])):
    single_cell_plotter(image_name, output_path,
                        cmap=ListedColormap('#e6000b'), threshold=threshold)
