import os
import numpy as np

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from functools import wraps
import time

from GEN_Utils import FileHandling
from loguru import logger

logger.info('Hi new user')

def timed(func):
    """This decorator prints the execution time for the decorated function."""
    @wraps(func)
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        end = time.time()
        logger.debug("{} ran in {}s".format(func.__name__, round(end - start, 2)))
        return result
    return wrapper

@timed
def pixel_cleaner(results_path, image_name):
    # read into pandas df
    raw_pixels = pd.read_csv(results_path)
    raw_pixels

    # Clean pixels info to generate single list with ROI, coords and intensity
    cleaning_pixels = raw_pixels.T.reset_index()
    cleaning_pixels

    cleaning_pixels['ROI_name'] = cleaning_pixels['index'].str.split('_').str[0]
    grouped_ROI = cleaning_pixels.set_index('ROI_name', drop=True).groupby('ROI_name')

    ROI_data_list = []
    for ROI, group in grouped_ROI:
        new_col_list = list(group['index'].str.split('_').str[1:3])
        ROI_data = group.T
        # Rename columns and drop previous column names
        ROI_data.columns = ['_'.join(x) for x in new_col_list]
        #ROI_data.columns = ['pos_X', 'pos_Y', 'Intensity_C1', 'Intensity_C2', 'Intensity_C3', 'Intensity_C4']
        ROI_data.drop('index', axis=0, inplace=True)
        # Remove leftover pixels from where imageJ assigns 0, 0 to the row - if all values are zero, we assume this is what happened
        ROI_data = ROI_data.replace(0, np.nan)
        ROI_data = ROI_data.dropna(how='all')
        ROI_data = ROI_data.replace(np.nan, 0)
        # add description columns for ROI name and image name
        ROI_data['ROI_name'] = ROI
        ROI_data['image_name'] = image_name
        ROI_data_list.append(ROI_data)

    ROI_pixels = pd.concat(ROI_data_list)
    ROI_pixels.reset_index(inplace=True, drop=True)
    ROI_pixels['X,Y'] = list(zip(ROI_pixels.X_pos, ROI_pixels.Y_pos))

    return ROI_pixels

@timed
def image_processor(input_folder, image_name, output_path, save_file=False):
    logger.info(f'Processing ROI pixels from {image_name}')
    results_path = input_folder+image_name+'_cells_pixels.csv'
    ROI_pixels = pixel_cleaner(results_path, image_name)
    ROI_pixels.rename(columns={'cells' : 'Intensity'}, inplace=True)
    logger.debug(ROI_pixels.columns.tolist())

    # Repeat for nuclear pixel info
    nuclear_path = input_folder+image_name+'_nuclei_pixels.csv'
    logger.info(f'Processing nuclear pixels from {image_name}')
    nuclei_pixels = pixel_cleaner(nuclear_path, image_name)
    nuclei_pixels.rename(columns={'nuclei' : 'Intensity'}, inplace=True)
    logger.debug(nuclei_pixels.columns.tolist())


    # Filter ROI_pixels to remove any pixels that are found in the nuclear pixels list
    logger.info(f'Filtering nuclear pixels from ROI pixels for {image_name}')
    cyto_pixels = ROI_pixels[~ROI_pixels['X,Y'].isin(nuclei_pixels['X,Y'])]
    logger.debug(cyto_pixels.columns.tolist())


    # Collect nuclear pixels that were excluded from whole cell pixels as nuclear pixels for each cell
    nuclear_pixels = ROI_pixels[~ROI_pixels['X,Y'].isin(cyto_pixels['X,Y'])]

    if save_file:
        # Save to excel file:
        logger.info(f'Saving individual excel files for {image_name}')
        FileHandling.df_to_excel(output_path+f'{image_name}_pixel_filtered.xlsx', ['cytoplasm', 'nucleus', 'whole_cell'], data_frames=[cyto_pixels, nuclear_pixels, ROI_pixels])

    logger.info(f'Pixels processed successfully for {image_name}! Results saved to {output_path}')

    return {'cyto':cyto_pixels, 'nuclei_pixels':nuclei_pixels, 'nucleus':nuclear_pixels, 'whole_cell':ROI_pixels}


def image_maker(df, channel='GFP_intensity'):
    array_im = np.zeros((512, 512))
    pixel_map = list(zip(df['X_pos'], df['Y_pos'], df[channel]))
    for pixel in pixel_map:
        x, y, val = pixel
        # array_im[x, y] = val
        # to make binary
        array_im[x, y] = 1
        # flip to maintain correct orientation
        new_arr = np.rot90(array_im, 3)
        final_arr = np.fliplr(new_arr)

    return final_arr


def image_maker(df, channel='GFP_intensity'):
    array_im = np.zeros((512, 512))
    pixel_map = list(zip(df['X_pos'], df['Y_pos'], df[channel]))
    for pixel in pixel_map:
        x, y, val = pixel
        array_im[x, y] = val
    return array_im


