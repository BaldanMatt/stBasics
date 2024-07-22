"""
author: matteo baldan
date: 2024-07-09
"""

# Importing libraries
import numpy as np
import pandas as pd
from math import cos, sin
import h5py 
import cv2
import tifffile
import matplotlib.pyplot as plt
import seaborn as sns

# GLOBAL VARIABLES
DATA_PATH = "/home/dati/merfish_mouse_brain/BrainReceptorShowcase"
DO_LOAD_HEAVY_DAPI = False 
DO_LOAD_HEAVY_TRANSCRIPTS = False
FRAC_LIGHT_TRANSCRIPTS = .25
SCALE_PERCENT = 5
RANDOM_STATE = 14

def main():
    slice_number = 1
    replicate_number = 1
    z_index_number = 2
    fov = 75
    data_path = DATA_PATH + f"/Slice{slice_number}" + f"/Replicate{replicate_number}"
    data_suffix = '_S' + str(slice_number) + 'R' + str(replicate_number)

    # Load transformation matrix
    filename = data_path + "/images/micron_to_mosaic_pixel_transform.csv"
    transformation_matrix = pd.read_csv(filename, header=None, sep = " ").values
    print("\033[35mThis is the transformation matrix from micro to mosaic pixel...\n\033[32m",transformation_matrix)
    # Load metadata for the cells
    meta_cell = pd.read_csv(data_path + "/cell_metadata" + data_suffix + \
                            ".csv", index_col = 0)
    meta_cell = meta_cell[meta_cell.fov == fov]
    print(meta_cell)

    # Load cell boundaries
    print("Reading cell boundaries...")
    filename = data_path + "/cell_boundaries/feature_data_" + str(fov) + '.hdf5'
    cellBoundaries = h5py.File(filename)
    
    # Collect boundaries in specific fov
    z_index = 'zIndex_' + str(z_index_number)
    currentCells = []
    for inst_cell in meta_cell.index.tolist():
        temp = cellBoundaries['featuredata'][inst_cell][z_index]['p_0']['coordinates'][0]
        boundaryPolygon = np.ones((temp.shape[0], temp.shape[1]+1))
        boundaryPolygon[:, :-1] = temp
        transformedBoundary = np.matmul(transformation_matrix,
                                        np.transpose(boundaryPolygon))[:-1]
        currentCells.append(transformedBoundary)
    minCoord = np.min([np.min(x, axis=1) for x in currentCells],
                      axis=0).astype(int)
    maxCoord = np.max([np.max(x, axis=1) for x in currentCells],
                      axis=0).astype(int)
    # Load DAPI Image and downsample image
    if DO_LOAD_HEAVY_DAPI:
        image = tifffile.imread(data_path + "/images/mosaic_DAPI_z" + \
            str(z_index_number) + '.tif')
        print('Loaded image in memory')
        width = int(image.shape[1] * SCALE_PERCENT / 100)
        height = int(image.shape[0] * SCALE_PERCENT / 100)
        dim = (width, height)

        # Resize
        resizedImage = cv2.resize(image, dim, interpolation = cv2.INTER_AREA)
        tifffile.imsave('figures/resized_image.tif', resizedImage)
        tifffile.imsave('figures/fov_image.tif', image[minCoord[1]:maxCoord[1],
                                                       minCoord[0]:maxCoord[0]])
        print("Saved resized image in /figures/fov_image.tif ...")

        del image, resizedImage
    
    img = tifffile.imread("figures/resized_image.tif")
    plt.figure(figsize=(15,15))
    plt.imshow(img, vmax=30000)
    plt.show()
    
    # Load transcript data
    print("Loading transcripts in memory...")
    if DO_LOAD_HEAVY_TRANSCRIPTS:
        transcripts = pd.read_csv(data_path + "/detected_transcripts" + data_suffix
                              + ".csv", index_col = 0)
        fractions = [0.25, 0.5, 0.75]
        # save reduced sizes of transcript data
        save_heavy_transcript = input("Do you want to save the heavy transcripts in lighter ones...? (Y/N)")
        if save_heavy_transcript == "Y":
            for inst_frac in fractions:
                transcripts_frac = transcripts.sample(frac=inst_frac,
                                                replace=False,
                                                random_state=RANDOM_STATE)
                transcripts_frac.to_csv("data/detected_transcripts" +
                                        data_suffix + "_" + str(inst_frac))
    else:
        transcripts = pd.read_csv("data/detected_transcripts" + data_suffix +
                                  "_" + str(FRAC_LIGHT_TRANSCRIPTS))
        print("Transcripts loaded in memory...")
    
    # Set gene colors
    gene_colors = {
            'Ntsr2':'teal',
            'Grin2b':'orangered',
            'olig1':'indigo',
            }
    # process data for data viz
    keep_genes = gene_colors.keys()
    transcripts_keep = transcripts[transcripts.gene.isin(keep_genes)]
    # rotate the mouse brain to the upright position
    transcripts_keep['-global_y'] = -transcripts_keep['global_y']
    # data viz
    fig, axs = plt.subplots(nrows=1,
                            ncols=1,
                            figsize=(15,15))

    sns.scatterplot(transcripts_keep,
                    x='global_x',
                    y='-global_y',
                    hue='gene',
                    s=0.1,
                    ax=axs)
    plt.show()

if __name__ == "__main__":
    main()

