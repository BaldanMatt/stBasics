"""
author: matteo baldan
date: 04/07/2024
"""
#import squidpy as sq
import h5py
import pandas as pd
import numpy as np
import tifffile
import cv2
import matplotlib.pyplot as plt

DATA_PATH = '/home/dati/merfish_mouse_brain/BrainReceptorShowcase'

def h5_tree(val, pre=''):
    items = len(val)
    for key, val in val.items():
        items -= 1
        if items == 0:
            # the last item
            if type(val) == h5py._hl.group.Group:
                print(pre + '└── ' + key)
                h5_tree(val, pre+'    ')
            else:
                try:
                    print(pre + '└── ' + key + ' (%d)' % len(val))
                except TypeError:
                    print(pre + '└── ' + key + ' (scalar)')
        else:
            if type(val) == h5py._hl.group.Group:
                print(pre + '├── ' + key)
                h5_tree(val, pre+'│   ')
            else:
                try:
                    print(pre + '├── ' + key + ' (%d)' % len(val))
                except TypeError:
                    print(pre + '├── ' + key + ' (scalar)')
def main():
    slice_num = 1
    replicate_num = 1
    data_path = DATA_PATH + '/Slice' + str(slice_num) + '/Replicate' + str(replicate_num)
    data_suffix = '_S' + str(slice_num) + 'R' + str(replicate_num)

    """
    adata = sq.read.vizgen(path=data_path,
                           counts_file="cell_by_gene"+data_suffix+".csv",
                           meta_file="cell_metadata"+data_suffix+".csv",
                           transformation_file="micron_to_mosaic_pixel_transform.csv")
    print(adata)

    # Evaluate the average position of all cells
    adata.obs["avg_x"] = (adata.obs["min_x"] + adata.obs["max_x"])/2
    adata.obs["avg_y"] = (adata.obs["min_y"] + adata.obs["max_y"])/2
    """

    # Load transformation matrix
    transformation_matrix = pd.read_csv(data_path + "/images/micron_to_mosaic_pixel_transform.csv",
                                        header=None,
                                        sep=' ').values

    # Load cell boundaries for a single fov
    fov = 75
    cellBoundaries = h5py.File(data_path + "/cell_boundaries" + "/feature_data_" + str(fov) + ".hdf5")
    print(f"\033[37mThis is the structure tree of the hdf5 file\033[32m")
    h5_tree(cellBoundaries)
    meta_cell = pd.read_csv(data_path + "/cell_metadata" + data_suffix + ".csv",
                            index_col = 0)
    meta_cell = meta_cell[meta_cell.fov == fov]
    print(f"\033[37mThis is the metadata dataframe for fov: {fov}\033[32m\n",meta_cell)
    print(f"\033[37mThis is the metadata dataframe for fov: {fov}\033[32m\n",meta_cell.columns)
    print(f"\033[37mThis is type of data structure that can read the cell boundaries of fov: {fov}\033[32m\n",cellBoundaries)

    z_index_number = 2
    z_index = "zIndex_" + str(z_index_number)
    # collect boundaries in fov
    currentCells = []
    for instance_cell in meta_cell.index.tolist(): # for each cell
        foo = cellBoundaries['featuredata'][instance_cell][z_index]['p_0']['coordinates'][0]
        boundaryPolygon = np.ones((foo.shape[0], foo.shape[1]+1))
        boundaryPolygon[:, :-1] = foo
        transformedBoundary = np.matmul(transformation_matrix, np.transpose(boundaryPolygon))[:-1]
        currentCells.append(transformedBoundary)
    minCoord = np.min([np.min(x, axis=1) for x in currentCells], axis=0).astype(int)
    maxCoord = np.max([np.max(x, axis=1) for x in currentCells], axis=0).astype(int)

    # Load DAPI image and downsample image
    image = tifffile.imread(data_path + "/images" + "/mosaic_DAPI_z" + str(z_index_number) + ".tif")
    print('Loaded image into memory')
    scale_percent = 5
    width = int(image.shape[1] * scale_percent / 100)
    height = int(image.shape[0] * scale_percent / 100)
    dim = (width, height)

    # reduce size
    resized = cv2.resize(image, dim, interpolation = cv2.INTER_AREA)
    print('Resized image..')
    tifffile.imsave(data_path + '/images' + f'/fov_{fov}_zIndex_{z_index_number}_image_size_{width}_{height}.tif',
                    image[minCoord[1]:maxCoord[1], minCoord[0]:maxCoord[0]])
    print('Saved image subset for single FOV section')

    # delete image
    del image
    fig = plt.figure(figsize=(15,15))
    plt.imshow(resized, vmax=30000)
    fig.savefig(f"figures/" + f"fov_{fov}_zIndex_{z_index_number}_image_size_{width}_{height}.png")

if __name__ == "__main__":
    main()