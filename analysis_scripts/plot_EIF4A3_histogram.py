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
# matplotlib.rcParams.update({'figure.facecolor': (1,1,1,1)})

logger.info("Import OK")

input_path = 'python_results/EIF4A3/pixel_calculations/GFP_pixel_summary.xlsx'
output_path = 'python_results/EIF4A3/histograms/'

if not os.path.exists(output_path):
    os.mkdir(output_path)

# Read in compiled pixel summary
raw_data = pd.read_csv(input_path)
pixels = raw_data.copy()
pixels.shape
pixels.head()

bins = np.arange(0, 6000, 100)

# plot histogram of 25Q versus 97Q nuclear pixels
nuclear_pixels = pixels[pixels['pixel_type'] == 'nucleus'].copy()
mutant_pos = {'25Q': 0, '97Q': 1}
col_pal = {'25Q': '#c9003c', '97Q': '#0600a6'}

# plot all pixel distributions together
fix, ax = plt.subplots(figsize=(8, 5))
for mutant, data in nuclear_pixels.groupby('mutant'):
    if not mutant == 'UT':
        for replicate, df in data.groupby('replicate'):
            sns.distplot(df['GFP_intensity'], hist=False, color=col_pal[mutant], label=mutant)
ax.set_ylim(0, 0.0015)
ax.set_xlim(0, 6000)
plt.ylabel('Density')
plt.xlabel('GFP Fluorescence Intensity (A.U.)')
plt.tight_layout()
plt.savefig(f'{output_path}nuclear_compiled_histogram.png')


# Plot on separate axes
fix, axes = plt.subplots(2, 1, figsize=(8, 8))
for mutant, data in nuclear_pixels.groupby('mutant'):
    if not mutant == 'UT':
        for replicate, df in data.groupby('replicate'):
            sns.distplot(df['GFP_intensity'], hist=False,
                         ax=axes[mutant_pos[mutant]], color=col_pal[mutant])
axes[0].set_title('25Q')
axes[1].set_title('97Q')
for ax in axes:
    ax.set_ylim(0, 0.0015)
    ax.set_xlim(0, 6000)
plt.ylabel('Density')
plt.xlabel('GFP Fluorescence Intensity (A.U.)')
plt.tight_layout()
plt.savefig(f'{output_path}nuclear_panel_histogram.png')


# plot compiled pixels from all cells
fix = plt.subplots(figsize=(8, 5))
for mutant, data in nuclear_pixels.groupby('mutant'):
    if not mutant == 'UT':
        sns.distplot(data['GFP_intensity'], hist=False, color=col_pal[mutant], label=mutant)
plt.ylim(0, 0.0025)
plt.xlim(0, 6000)
plt.ylabel('Density')
plt.xlabel('GFP Fluorescence Intensity (A.U.)')
plt.tight_layout()
plt.savefig(f'{output_path}nuclear_average_histogram.png')


# plot histogram of 25Q versus 97Q CYTOPLASM pixels
cyto_pixels = pixels[pixels['pixel_type'] == 'cyto'].copy()
mutant_pos = {'25Q': 0, '97Q': 1}
col_pal = {'25Q': '#c9003c', '97Q': '#0600a6'}

# plot all pixel distributions together
fix, ax = plt.subplots(figsize=(8, 5))
for mutant, data in cyto_pixels.groupby('mutant'):
    if not mutant == 'UT':
        for replicate, df in data.groupby('replicate'):
            sns.distplot(df['GFP_intensity'], hist=False,
                         color=col_pal[mutant], label=mutant)
ax.set_ylim(0, 0.005)
ax.set_xlim(0, 6000)
plt.ylabel('Density')
plt.xlabel('GFP Fluorescence Intensity (A.U.)')
plt.tight_layout()
plt.savefig(f'{output_path}cyto_compiled_histogram.png')


# Plot on separate axes
fix, axes = plt.subplots(2, 1, figsize=(8, 8))
for mutant, data in cyto_pixels.groupby('mutant'):
    if not mutant == 'UT':
        for replicate, df in data.groupby('replicate'):
            sns.distplot(df['GFP_intensity'], hist=False,
                         ax=axes[mutant_pos[mutant]], color=col_pal[mutant])
axes[0].set_title('25Q')
axes[1].set_title('97Q')
for ax in axes:
    ax.set_ylim(0, 0.005)
    ax.set_xlim(0, 6000)
plt.ylabel('Density')
plt.xlabel('GFP Fluorescence Intensity (A.U.)')
plt.tight_layout()
plt.savefig(f'{output_path}cyto_panel_histogram.png')


# plot compiled pixels from all cells
fix = plt.subplots(figsize=(8, 5))
for mutant, data in cyto_pixels.groupby('mutant'):
    if not mutant == 'UT':
        sns.distplot(data['GFP_intensity'], hist=False,
                     color=col_pal[mutant], label=mutant)
plt.ylim(0, 0.005)
plt.xlim(0, 6000)
plt.ylabel('Density')
plt.xlabel('GFP Fluorescence Intensity (A.U.)')
plt.tight_layout()
plt.savefig(f'{output_path}cyto_average_histogram.png')


# Plot compilation plots as svg, with inclusion separated
# plot compiled pixels from all cells
fix = plt.subplots(figsize=(8, 5))
for mutant_info, data in nuclear_pixels.groupby(['mutant', 'inclusion']):
    if not 'UT' in mutant_info:
        mutant, inclusion = mutant_info
        sns.distplot(data['GFP_intensity'], hist=False, label=mutant_info)
plt.ylim(0, 0.002)
plt.xlim(0, 6000)
plt.ylabel('Density')
plt.xlabel('GFP Fluorescence Intensity (A.U.)')
plt.tight_layout()
plt.savefig(f'{output_path}nuclear_average_histogram.svg')

# plot compiled pixels from all cells
fix = plt.subplots(figsize=(8, 5))
for mutant_info, data in cyto_pixels.groupby(['mutant', 'inclusion']):
    if not 'UT' in mutant_info:
        mutant, inclusion = mutant_info
        sns.distplot(data['GFP_intensity'], hist=False, label=mutant_info)
plt.ylim(0, 0.002)
plt.xlim(0, 6000)
plt.ylabel('Density')
plt.xlabel('GFP Fluorescence Intensity (A.U.)')
plt.tight_layout()
plt.savefig(f'{output_path}cyto_average_histogram.svg')

# plot ALL compiled pixels from ALL cells
fix = plt.subplots(figsize=(8, 5))
for mutant_info, data in pixels.groupby(['mutant', 'inclusion', 'pixel_type']):
    if not 'UT' in mutant_info:
        mutant, inclusion, pixel_type = mutant_info
        sns.distplot(data['GFP_intensity'], hist=False, label=mutant_info)
plt.ylim(0, 0.002)
plt.xlim(0, 6000)
plt.ylabel('Density')
plt.xlabel('GFP Fluorescence Intensity (A.U.)')
plt.tight_layout()
plt.savefig(f'{output_path}average_histogram.svg')
