"""
author: matteo baldan
date: 28-06-2024
description: SPARC tutorial - the one on their documentation. We are using MERFISH data and not their VISIUM dataset
"""

import squidpy as sq
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import SPARC
import numpy as np
import graphtools

VIZGEN_PATH = "/home/dati/merfish_mouse_brain/BrainReceptorShowcase/Slice1/Replicate1/"


def main():
    adata = sq.read.vizgen(path=VIZGEN_PATH,
                           counts_file="cell_by_gene_S1R1.csv",
                           meta_file="cell_metadata_S1R1.csv",
                           transformation_file="micron_to_mosaic_pixel_transform.csv", )
    print(adata)

    """
    Remove not unique genes and collapse them in some others with string "MT-"
    """
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # calculate quality control metrics
    sc.pp.calculate_qc_metrics(adata,
                               percent_top=(50, 100, 200, 300),
                               inplace=True)

    # visualize the quality of the sample
    fig, axs = plt.subplots(2, 2, figsize=(15, 10))
    axs = axs.ravel()
    sns.distplot(adata.obs["total_counts"], kde=False, bins=100, ax=axs[0])
    axs[0].title.set_text("Total counts within cells\nNcells:78329")
    sns.distplot(adata.obs["total_counts"][adata.obs["total_counts"] < 1e2], kde=False, bins=100, ax=axs[1])
    axs[1].title.set_text("Zoom to cells with low counts")
    sns.distplot(adata.obs["n_genes_by_counts"], kde=False, bins=100, ax=axs[2])
    axs[2].title.set_text("Num of unique genes per cell")
    sns.distplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4e1], kde=False, bins=100, ax=axs[3])
    axs[3].title.set_text("Zoom to cells with low num of unique genes")
    plt.show()

    # filter cells
    sc.pp.filter_cells(adata, min_counts=100)
    sc.pp.filter_cells(adata, max_counts=150)
    print(f"#cells after filter: {adata.n_obs}")

    """
    preprocessing
    """
    print("Normalize...")
    sc.pp.normalize_total(adata, inplace=True)
    print("log transformation...")
    sc.pp.log1p(adata)
    print("evaluate Variable genes...")
    sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=400)
    print("PCA...")
    sc.pp.pca(adata)
    print("Find neighbors...")
    sc.pp.neighbors(adata)
    print("Umap reduction...")
    sc.tl.umap(adata)
    print("Leiden clustering...")
    sc.tl.leiden(adata,
                 flavor="igraph",
                 n_iterations=2,
                 key_added="clusters")

    # visualize reduction results
    plt.rcParams["figure.figsize"] = (4, 4)
    sc.pl.umap(adata,
               color=["total_counts",
                      "n_genes_by_counts",
                      "clusters"],
               wspace=0.4)
    plt.rcParams["figure.figsize"] = (8, 8)
    sq.pl.spatial_scatter(adata,
                          color=["total_counts",
                                 "n_genes_by_counts",
                                 "clusters"],
                          size=0.5,
                          library_id=["spatial"],
                          shape=None)
    plt.show()

    """
    testing spARC
    """
    sparc_op = SPARC.spARC()
    adata.obs["avg_x"] = adata.obs["min_x"] + adata.obs["max_x"] / 2
    adata.obs["avg_y"] = adata.obs["min_y"] + adata.obs["max_y"] / 2
    data_sparc = sparc_op.fit_transform(expression_X=adata.to_df(),
                                        spatial_X=adata.obs[['avg_x', 'avg_y']], )

    adata.layers["rna_raw"] = adata.to_df()
    adata.layers["rna_sparc"] = data_sparc
    # 14819x483 uses about 12 GB of memory

    plt.rcParams["figure.figsize"] = (2, 2)
    sc.set_figure_params(facecolor="white", transparent=True, frameon=True, fontsize=20)
    sq.pl.spatial_scatter(adata,
                          size=1.5,
                          layer="rna_raw",
                          color=["Npy2r",
                                 "Slc17a7"],
                          library_id=["spatial"],
                          shape=None)
    plt.show()
    sq.pl.spatial_scatter(adata,
                          size=1.5,
                          layer="rna_sparc",
                          color=["Npy2r",
                                 "Slc17a7",
                                 ],
                          library_id=["spatial"],
                          shape=None)
    plt.show()

    """
    Check soluble factors
    """
    spatial_data = np.array(adata.obs[['avg_x', 'avg_y']])
    G_spatial = graphtools.Graph(spatial_data, knn=5)
    soluable_sparc = sparc_op.diffuse_soluble_factors(soluble_spatial_graph=G_spatial)
    adata.layers["soluble_sparc"] = soluable_sparc
    plt.rcParams["figure.figsize"] = (2, 2)
    sc.set_figure_params(facecolor="white", transparent=True, frameon=True, fontsize=20)
    sq.pl.spatial_scatter(adata,
                          size=1.5,
                          layer="soluble_sparc",
                          color=["Npy2r",
                                 "Slc17a7",
                                 "Gabbr2"],
                          library_id=["spatial"],
                          shape=None)
    plt.show()
if __name__ == "__main__":
    main()
