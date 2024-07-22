"""
author: matteo baldan
date: 04/07/2024
"""

# import squidpy as sq
import h5py
import pandas as pd
import numpy as np
import tifffile
import cv2
import matplotlib.pyplot as plt
import seaborn as sns

DATA_PATH = "/home/dati/merfish_mouse_brain/BrainReceptorShowcase"
DO_LOAD_IMAGE = True
DO_PLOT_IMAGE = False


def h5_tree(val, pre=""):
    items = len(val)
    for key, val in val.items():
        items -= 1
        if items == 0:
            # the last item
            if type(val) == h5py._hl.group.Group:
                print(pre + "└── " + key)
                h5_tree(val, pre + "    ")
            else:
                try:
                    print(pre + "└── " + key + " (%d)" % len(val))
                except TypeError:
                    print(pre + "└── " + key + " (scalar)")
        else:
            if type(val) == h5py._hl.group.Group:
                print(pre + "├── " + key)
                h5_tree(val, pre + "│   ")
            else:
                try:
                    print(pre + "├── " + key + " (%d)" % len(val))
                except TypeError:
                    print(pre + "├── " + key + " (scalar)")


def main():
    slice_num = 1
    replicate_num = 1
    data_path = (
            DATA_PATH + "/Slice" + str(slice_num) + "/Replicate" + str(replicate_num)
    )
    data_suffix = "_S" + str(slice_num) + "R" + str(replicate_num)

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
    transformation_matrix = pd.read_csv(
        data_path + "/images/micron_to_mosaic_pixel_transform.csv", header=None, sep=" "
    ).values

    # Load cell boundaries for a single fov
    fov = 75
    cellBoundaries = h5py.File(
        data_path + "/cell_boundaries" + "/feature_data_" + str(fov) + ".hdf5"
    )
    print(f"\033[37mThis is the structure tree of the hdf5 file\033[32m")
    h5_tree(cellBoundaries)
    meta_cell = pd.read_csv(
        data_path + "/cell_metadata" + data_suffix + ".csv", index_col=0
    )
    meta_cell = meta_cell[meta_cell.fov == fov]
    print(f"\033[37mThis is the metadata dataframe for fov: {fov}\033[32m\n", meta_cell)
    print(
        f"\033[37mThis is the metadata dataframe for fov: {fov}\033[32m\n",
        meta_cell.columns,
    )
    print(
        f"\033[37mThis is type of data structure that can read the cell boundaries of fov: {fov}\033[32m\n",
        cellBoundaries,
    )

    z_index_number = 2
    z_index = "zIndex_" + str(z_index_number)
    # collect boundaries in fov
    currentCells = []
    for instance_cell in meta_cell.index.tolist():  # for each cell
        foo = cellBoundaries["featuredata"][instance_cell][z_index]["p_0"][
            "coordinates"
        ][0]
        boundaryPolygon = np.ones((foo.shape[0], foo.shape[1] + 1))
        boundaryPolygon[:, :-1] = foo
        transformedBoundary = np.matmul(
            transformation_matrix, np.transpose(boundaryPolygon)
        )[:-1]
        currentCells.append(transformedBoundary)
    minCoord = np.min([np.min(x, axis=1) for x in currentCells], axis=0).astype(int)
    maxCoord = np.max([np.max(x, axis=1) for x in currentCells], axis=0).astype(int)

    if DO_LOAD_IMAGE:
        # Load DAPI image and downsample image
        image = tifffile.imread(
            data_path + "/images" + "/mosaic_DAPI_z" + str(z_index_number) + ".tif"
        )
        print("Loaded image into memory")
        scale_percent = 5
        width = int(image.shape[1] * scale_percent / 100)
        height = int(image.shape[0] * scale_percent / 100)
        dim = (width, height)

        # reduce size
        resized = cv2.resize(image, dim, interpolation=cv2.INTER_AREA)
        print("Resized image..")
        tifffile.imsave(
            data_path
            + "/images"
            + f"/fov_{fov}_zIndex_{z_index_number}_image_size_{width}_{height}.tif",
            image[minCoord[1]: maxCoord[1], minCoord[0]: maxCoord[0]],
        )
        print("Saved image subset for single FOV section")

        # delete image
        del image
        fig = plt.figure(figsize=(15, 15))
        plt.imshow(resized, vmax=30000)
        fig.savefig(
            f"figures/"
            + f"fov_{fov}_zIndex_{z_index_number}_image_size_{width}_{height}.png"
        )

    """
    Load transcripts 
    the data are contianed in a .csv file of all detected transcripts in the sample where each row is a detected transcirpt.
    The column names and descriptions of what the column represents are:
    1. barcode_id: internally used id of the gene. Each gene has a unique barcode
    2. global_x and global_y: micron x and micron y coordinates
    3. x,y: pixel coorinate within field of view
    4. global_z: index of z slice that this transcript was detected in. Each z slice is separated by 1.5 microns
    5. fov: the index of the field of view where this transcript was detected
    6. gene: the gene name of this detected transcript
    """
    print("Loading transcripts...")
    transcripts = pd.read_csv(
        data_path + "/detected_transcripts" + data_suffix + ".csv", index_col=0
    )
    print(transcripts.describe)
    if DO_PLOT_IMAGE:
        # Select genes, in parituclar the most highly expressed ones
        transcripts["gene"].value_counts().head(20)
        gene_colors = {"Ntsr2": "teal", "Grin2b": "orangered", "Olig1": "indigo"}
        # code for compressing data for visuals
        from math import sin, cos

        # process data for interactive visual
        keep_genes = gene_colors.keys()
        transcripts_keep = transcripts[transcripts.gene.isin(keep_genes)]
        transcripts_keep["-global_y"] = -transcripts_keep["global_y"]
        # transcripts_keep = transcripts_keep.sample(frac=0.5, replace=False, random_state=44)
        # rotate the mouse brain to the upright position
        theta = np.deg2rad(-195)  # -15
        rot = np.array([[cos(theta), -sin(theta)], [sin(theta), cos(theta)]])
        transcripts_keep[["global_x", "-global_y"]] = transcripts_keep[
            ["global_x", "-global_y"]
        ].dot(rot)

        transcripts_keep["name"] = transcripts_keep.loc[:, "gene"]
        print(type(transcripts_keep))
        print(transcripts_keep.describe)
        fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(15, 15))
        sns.scatterplot(
            transcripts_keep,
            x="global_x",
            y="-global_y",
            hue="name",
            palette=gene_colors,
            ax=axs,
            s=0.3,
        )
        plt.show()
    # Let's keep in memory only the transcripts in our field of view of interest
    transcripts = transcripts[transcripts["fov"] == fov].copy()

    """
    Load single Field of view data
    """
    from PIL import Image
    import base64
    from io import BytesIO
    from copy import deepcopy

    # transcripts
    tmp = transcripts[["global_x", "global_y"]].values
    transcript_positions = np.ones((tmp.shape[0], tmp.shape[1] + 1))
    transcript_positions[:, :-1] = tmp
    # transform coordinates to mosaic pixel coordinates
    transformed_positions = np.matmul(
        transformation_matrix, np.transpose(transcript_positions)
    )[:-1]
    transcripts.loc[:, "local_x"] = transformed_positions[0, :]
    transcripts.loc[:, "local_y"] = transformed_positions[1, :]
    print("Loading small tiff in memory...")
    image = tifffile.imread(
        data_path
        + "/images"
        + f"/fov_{fov}_zIndex_{z_index_number}_image_size_{width}_{height}.tif"
    )

    # segmentation data
    polygon_data = []
    print("Looping across cells....")
    for inst_index in range(len(currentCells)):
        inst_cell = currentCells[inst_index]
        df_poly_z = pd.DataFrame(inst_cell).transpose()
        print(df_poly_z.columns.tolist())
        inst_name = meta_cell.iloc[inst_index].name
        inst_poly = {"coordinates": df_poly_z.values.tolist(), "name": inst_name}
        polygon_data.append(inst_poly)

    df_obs = transcripts[["local_x", "local_y", "global_z", "gene"]]
    df_obs.columns = ["x", "y", "z", "name"]
    scatter_data = df_obs.to_dict("records")

    def make_image_url(image_data):
        buffered = BytesIO()
        noise_img = Image.fromarray(image_data, "RGB")
        noise_img.save(buffered, format="PNG")
        noise_img_str = base64.b64encode(buffered.getvalue())
        noise_url = '"data:image/png;base64,' + str(noise_img_str)[2:-1] + '"'
        return noise_url

    # image background
    vmax = 30000
    tmpImage = deepcopy(np.flip(image, axis=0))
    tmpImage[tmpImage > vmax] = vmax
    tmpImage = tmpImage * (255 / vmax)
    tmpImage[tmpImage > 255] = 255
    tmpImage = tmpImage.astype(np.int64)

    image_data = np.zeros((image.shape[0], image.shape[1], 3), dtype=np.uint8)
    image_data[:, :, 0] = tmpImage
    image_data[:, :, 1] = tmpImage
    image_data[:, :, 2] = tmpImage
    image_url = make_image_url(image_data)

    # gene colors
    cats = df_obs.name.unique().tolist()
    colors = []
    for inst_index in range(len(cats)):
        if "Blank" in cats[inst_index]:
            colors.append("#d3d3d3")
        else:
            mod_index = inst_index % 60
            if mod_index < 20:
                color_rgba = plt.cm.tab20(mod_index % 20)
            elif mod_index < 40:
                color_rgba = plt.cm.tab20b(mod_index % 20)
            elif mod_index < 60:
                color_rgba = plt.cm.tab20b(mod_index % 20)
            inst_color = "#{:02x}{:02x}{:02x}".format(
                int(color_rgba[0] * 255),
                int(color_rgba[1] * 255),
                int(color_rgba[2] * 255),
            )
            colors.append(inst_color)
    cat_colors = dict(zip(cats, colors))

    x_offset = np.float64(minCoord[0])
    y_offset = np.float64(minCoord[1])
    image_bounds = [
        [x_offset, y_offset],
        [x_offset, y_offset + image.shape[0]],
        [x_offset + image.shape[1], y_offset + image.shape[0]],
        [x_offset + image.shape[1], y_offset],
    ]

    # Preparing dictionary for fov_input
    fov_inputs = {
        "polygon_data": polygon_data,
        "scatter_data": scatter_data,
        "center_x": np.mean([minCoord[0], maxCoord[0]]),
        "center_y": np.mean([minCoord[1], maxCoord[1]]),
        "image": image_url,
        "image_bounds": image_bounds,
        "zoom": -1.5,
        "min_zoom": -3,
        "height": 800,
        "ini_view_type": "2D",
        "cat_colors": cat_colors,
    }

    print("Finish preprocessing data to visualize...")
    print(scatter_data)
    print(polygon_data)

    # fig, axs = plt.subplots(nrows = 1, ncols = 1, figsize=(15,15))
    # sns.scatterplot(data=scatter_data,
    #            x = '',
    #            y = '',
    #            )


if __name__ == "__main__":
    main()
