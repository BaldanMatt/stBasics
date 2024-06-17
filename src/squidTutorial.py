"""
author: matteo
date: 2024-06-11
description:
    This tutorial shows how to apply Squidpy for the analysis of Merfish data
"""
import scanpy as sc
import squidpy as sq
import seaborn as sns
import matplotlib.pyplot as plt

VIZGEN_PATH = "/home/dati/merfish_mouse_brain/BrainReceptorShowcase/Slice1/Replicate1/"


def tutorialMerfishData():
    sc.logging.print_header()
    print(f"squidpy=={sq.__version__}")

    # Load the data
    adata = sq.datasets.merfish()
    # Let's familiarize with the AnnData object
    """
    The anndata object:
        is specifically designed for matrix-like data. 
        We usually have:
            n observations, each of which can be represented as d-dimensional vectors
            d dimensions, corresponding to variable or features
        both rows and columns of this n x d matrix are indexed
        
    """
    print(adata)
    print(adata.X)
    print("\033[36mThese are the observation column names: \033[37m\n", \
          adata.obs.columns)
    print("\033[36mThese are the observation names, i.e., the Cell_ID in this case: \033[37m\n", \
          adata.obs_names)
    print("\033[36mThese are the observation shape: \033[37m\n", \
          adata.obs.shape, end="\n" + "%" * 10 + "\n")
    print("\033[36mThese are the var column names: \033[37m\n", \
          adata.var.columns)
    print("\033[36mThese are the var shape, i.e., the Cell_ID in this case: \033[37m\n", \
          adata.var.shape, end="\n" + "%" * 10 + "\n")
    # Apparently there are no var columns... Let's inspect the values
    print("\033[36mThese are the var values, i.e., the 161x0 values: \033[37m\n", \
          adata.var.values, end="\n" + "%" * 10 + "\n")  # empty...
    # Let's try to index it
    print("\033[36mLet's inspect the whole first 3 rows of the anndata: \033[37m\n", \
          adata[0:3, ].X, end="\n" + "%" * 10 + "\n")

    """
    We see that our anndata has observation/variable-level matrices.
        such data can describe metadata at either level and can have many dimensions to it. 
        you could think of a UMAP embedding of the data. 
        The restriction is that .obsm must length equal to the number of observations (n_obs). 
    i.e., the obsm or varsm. In this case we only have the obsm.
      obsm: 'spatial', 'spatial3d'
      
    Our data has also unstructured metadata, which allows for any unstructured metadata. This
    can be anything, like a list or a dictionary with general information.
    """
    print("\033[36mThese are the observation-level metadata: \033[37m\n", \
          adata.obsm, end="\n" + "%" * 10 + "\n")
    print("\033[36mThese are the observation-level 'spatial' metadata: \033[37m\n", \
          adata.obsm["spatial"], end="\n" + "%" * 10 + "\n")  # should be (x,y) for each obs
    print("\033[36mThese are the unstructured data: \033[37m\n", \
          adata.uns, end="\n" + "%" * 10 + "\n")  # cell class colors

    """
    In the anndata there may be multiple layers. Layers might be different forms of our original
    core data, perhaps one that is normalized and one that is not.
    """
    print("\033[36mThese are the layers: \033[37m\n", \
          adata.layers, end="\n" + "%" * 10 + "\n")  # no layers are present yet

    # This dataset consists of consecutive slices from the mouse hypothalamic preoptic region.
    # It represents an interesting example of how to work with 3D spatial data in squidpy.

    # let's visualize the 3D stack of slides
    sc.pl.embedding(adata,
                    basis="spatial3d",
                    projection="3d",
                    color="Cell_class",
                    show=True, )

    # Neighborhood enrichment analysis in 3D
    """
    It is important to consider whether the analysis should be performed on the 3D spatial coordinates or the 2D
    coordinates for a single slice.
    Squidpy functions that make use of the spatial graph can already support 3D coordinates.
    """
    # Building a spatial neighbors graph
    sq.gr.spatial_neighbors(adata,
                            coord_type="generic",
                            spatial_key="spatial3d", )
    print("\033[36mThese are the new unstructured metadata: \033[37m\n", \
          adata.uns["spatial_neighbors"], end="\n" + "%" * 10 + "\n")
    print("\033[36mThese are the evaluated spatial connectivities: \033[37m\n", \
          adata.obsp["spatial_connectivities"], end="\n" + "%" * 10 + "\n")  # it appears to be an adjacent matrix
    # adata.obsp appears to be saving only those tuples that have a 1 value, therefore and edge between obs
    print("\033[36mThese are the evaluated spatial_distances: \033[37m\n", \
          adata.obsp["spatial_distances"], end="\n" + "%" * 10 + "\n")  # it contains all the actual absolute distances
    # the distances are between the same tuples of obs that are saved in "spatial_connectivities"
    # Compute the nhood enhancement
    sq.gr.nhood_enrichment(adata,
                           cluster_key="Cell_class")
    sq.pl.nhood_enrichment(adata,
                           cluster_key="Cell_class",
                           method="single",
                           cmap="inferno",
                           vmin=-50,
                           vmax=100,
                           figsize=(8, 6))
    plt.show()

    # We can visualize some of the co-enriched clusters with the embedding plot
    sc.pl.embedding(adata,
                    basis="spatial3d",
                    groups=["Microglia", "Inhibitory", "Excitatory"],
                    na_color=(1, 1, 1, 0),  # (R,G,B,a)
                    projection="3d",
                    color="Cell_class",
                    show=True, )

    # We can also visualize the gene expression in 3D coordinates
    """
    To do so, usually we mean to perform differnetial expression testing.
    Which can be done using scanpy.tl.rank_genes_groups()
    """
    sc.tl.rank_genes_groups(adata,
                            groupby="Cell_class", )
    print("\033[36mThis is what rank genes do: \033[37m\n", \
          adata.uns["rank_genes_groups"], "\n",
          type(adata.uns["rank_genes_groups"]), end="\n" + "%" * 10 + "\n")
    print("\033[36mThese are the keys f such ranking: \033[37m\n", \
          adata.uns["rank_genes_groups"].keys(), end="\n" + "%" * 10 + "\n")

    sc.pl.rank_genes_groups(adata,
                            groupby="Cell_class",
                            show=True)
    # And the respective expression in 3D
    sc.pl.embedding(adata,
                    basis="spatial3d",
                    projection="3d",
                    color=["Aldh1l1", "Sgk1"])
    # If the same needs to be performed on a single slice
    adata_slice = adata[adata.obs.Bregma == -9].copy()
    sq.pl.spatial_scatter(adata_slice,
                          color="Cell_class",
                          shape=None,
                          groups=["Ependymal", "Microglia", "OD Mature 2"],
                          size=10)
    plt.show()

    # Spatially variable genes with spatial autocorrelation statistics
    """"
    With squidpy we can use autocorrelation statistics (Moran's I and Geary's C scores) to understand if there
    are any SVGs. These two statistics provide a score on the degree of spatial variability of gene expression.
    The statistic as well as the significativity are computed for each gene and FDR correction is performed.
    """
    # Let's first compute the Moran's I index
    sq.gr.spatial_autocorr(adata_slice,
                           mode="moran", )
    print("\033[36mThese are the genes that have the highest moran I: \033[37m\n", \
          adata_slice.uns["moranI"].sort_values(by="I", ascending=False).head(10), end="\n" + "%" * 10 + "\n")
    sq.pl.spatial_scatter(adata_slice,
                          shape=None,
                          color=["Nnat", "Sln", "Cd24a"],
                          size=3)
    plt.show()
    sq.gr.spatial_autocorr(adata_slice,
                           mode="geary", )
    print("\033[36mThese are the genes that have the lowest geary C: \033[37m\n", \
          adata_slice.uns["gearyC"].sort_values(by="C", ascending=True).head(10), end="\n" + "%" * 10 + "\n")
    sq.pl.spatial_scatter(adata_slice,
                          shape=None,
                          color=["Nnat", "Sln", "Cd24a"],  # the 3 lowest geary C are the 3 highest moran I
                          size=3)
    plt.show()


def tutorialVizgenData(path: str):
    """
    In this second tutorial, we will learn how to read VIZGEN DATA. First by using the squidpy utilities.
    The results is an AnnData that we will thoroughly investigate.
    """
    # You must have the structure of the vizgen data:
    # Replicate/
    #   counts_file.csv
    #   meta_file.csv
    #   ...
    #   images/
    #       transformation_file.csv
    #       ...
    adata = sq.read.vizgen(path=VIZGEN_PATH,
                           counts_file="cell_by_gene_S1R1.csv",
                           meta_file="cell_metadata_S1R1.csv",
                           transformation_file="micron_to_mosaic_pixel_transform.csv", )
    print("\033[36mThese are the observation column names: \033[37m\n", \
          adata.obs.columns)
    print("\033[36mThese are the observation names (cells?): \033[37m\n", \
          adata.obs_names)  # 78329 unique cells?
    print("\033[36mThese are the var names (genes?): \033[37m\n", \
          adata.var_names)  # 483 unique genes
    print("\033[36mWhat is fov? Field of view? Capiamo...: \033[37m\n", \
          adata.obs["fov"], "\n",
          type(adata.obs["fov"]))  # seems to be the ids of where in the slice the cell is

    # Calculate the quality control metrics
    """
    We can calculate the quality by using scanpy.pp.calculate_qc_metrics
        What does it do?
            How does it do it?
            ...
            Adds the following obs fields:
                n_genes_by_counts: ...
                log1p_n_genes_by_counts: ...
                total_counts: ...
                log1p_total_counts: ...
                pct_counts_in_top_50_genes: ...
                pct_counts_in_top_100_genes: ...
                pct_counts_in_top_200_genes: ...
                pct_counts_in_top_300_genes: ...
            Adds the following var fields:
                n_cells_by_counts: ...
                mean_counts: ...
                log1p_mean_counts: ...
                pct_dropout_by_counts: ...
                total_counts: ...
                log1p_total_counts: ...
    """
    sc.pp.calculate_qc_metrics(adata,
                               percent_top=(50, 100, 200, 300),
                               inplace=True)
    sns.jointplot(data=adata.obs,
                  x="log1p_total_counts",
                  y="log1p_n_genes_by_counts",
                  kind="hex")
    print("\033[36madata after having calculated qc metrics (new field in var or uns?): \033[37m\n", \
          adata)
    print("\033[36mPercentage of blank genes: \033[37m\n", \
          adata.obsm["blank_genes"].to_numpy().sum() / adata.var["total_counts"].sum() * 100)
    # this percentage can be used to estimate FDR
    """
    What are the blank genes? they are unassigned transcripts """  # Let's plot distribution of transcripts per cell, etc..
    fig, axs = plt.subplots(4, 1, figsize=(16, 9))
    axs[0].set_title("Total transcripts per cell")
    sns.histplot(data=adata.obs["total_counts"],
                 ax=axs[0],
                 kde=False, )
    axs[1].set_title("Unique transcripts per cell")
    sns.histplot(data=adata.obs["n_genes_by_counts"],
                 ax=axs[1],
                 kde=False, )
    axs[2].set_title("Transcripts per FOV")
    sns.histplot(data=adata.obs.groupby("fov").sum()["total_counts"],
                 ax=axs[2],
                 kde=False, )
    axs[3].set_title("Volume of segmented cells")  # metti unità di misura
    sns.histplot(data=adata.obs["volume"],
                 ax=axs[3],
                 kde=False)
    plt.subplots_adjust(hspace=0.5)
    plt.show()

    """
    All the cells that do not contain at least 10 transcripts are filtered out with sc.pp.filter_cells
    Also the genes can be similarly filtered with sc.pp.filter_genes.
    """
    sc.pp.filter_cells(adata, min_counts=10)

    """
    After the filtering, one can proceed with the preprocessing and try to understand the following. 
    We would like to reduce the dimensionality, however, the PCA needs to be applied over normalized values. 
    -> Normalization:  sc.pp.normalize_total() 
        normalize counts per cell. Normalize each cell by total counts over all genes, so that every cell has the same
        total count after normalization. If choosing `target_sum=1e6`, this is CPM normalization.
        If `exclude_highly_expressed=True`, very highly expressed genes are excluded from the computation of the
        normalization factor (size factor) for each cell. This is meaningful as these can strongly influence the
        resulting normalized values for all other genes.
        similar functions...
            Seurat
            Cell Ranger
            SPRING
    """
    # We copy the X layer and call it "counts"
    adata.layers["counts"] = adata.X.copy()
    # We evaluate the highly-variable-genes
    sc.pp.highly_variable_genes(adata,
                                flavor="seurat_v3",
                                n_top_genes=4000)
    sc.pp.normalize_total(adata, inplace=True)  # normalize
    sc.pp.log1p(adata)  # transform
    sc.pp.pca(adata)  # principal coordinate analysis
    sc.pp.neighbors(adata)  # evaluate the neighborhood graph. K-NN?
    sc.tl.umap(adata)  # visualize the reduced dimensionality
    sc.tl.leiden(adata,
                 flavor="igraph",
                 n_iterations=2,
                 directed=False)  # cluster the reduced dimensionality data.

    # Let's visualize the umap using the leidne clustering to display the different cell types
    sc.pl.umap(adata,
               color=["total_counts",
                      "n_genes_by_counts",
                      "leiden",
                      ],
               wspace=0.4)
    # Let's visualize in a scatter map
    sq.pl.spatial_scatter(adata,
                          shape=None,
                          color=[
                              "leiden",
                          ], )
    """
    From this point clusters can be annotated by differentially expressed genes, e.g., `sc.tl.rank_genes_groups` or
    by integrating scRNA and transferring labels e.g. Tangram or Harmony
    """
    # Annotating clusters
    ...

    """
    Computation of spatial statistics - Building the spatial neighbors graphs
        this example shows how to compute centrality scores, given a spatial graph and cell type annotation
        scores:
            - closeness centrality : measure of how close the group is to other nodes
            - clustering coefficient : measure of the degree to which nodes cluster together
            - degree centrality : fraction of non-group members connected to group members
            all scores are descriptive statistics of the spatial graph.
        In this example we have used leiden algorithm to cluster groups of cells. We will use such clustering to
        compute the centrality scores.
    """
    # sq.gr.spatial_neighbors() evaluate the connectivity matrix from spatial coordinates
    sq.gr.spatial_neighbors(adata,
                            coord_type="generic",
                            delaunay=True,)

    # Compute centrality scores
    sq.gr.centrality_scores(adata,
                            cluster_key="leiden")
    sq.pl.centrality_scores(adata,
                            cluster_key="leiden",
                            figsize=(16,5))

    """
    Computation of spatial statistics - Compute co-occurrence probability
    The co-occurrence score is defined as:
        p(exp|cond)/p(exp)
    where the numerator is the conditional probability of observing a cluster exp conditioned on the presence
    of a cluster cond, whereas p(exp) is the probability of observing exp in the radius size of interest.
    """
    # let's take a subsample (computation costs??) -> co_occurrence creates a (obs x obs) matrix.
    adata_subsample = sc.pp.subsample(adata,
                                      fraction=0.5,
                                      copy=True)
    # compute co-occurrence
    sq.gr.co_occurrence(adata_subsample,
                        cluster_key="leiden",)
    # show the co-occurrence for a number of clusters that we want
    sq.pl.co_occurrence(adata_subsample,
                        cluster_key="leiden",
                        clusters="12",
                        figsize=(10,10))
    sq.pl.spatial_scatter(adata_subsample,
                          color="leiden",
                          shape=None,
                          size=2)

    """
    Computation of spatial statistics - Neighbors enrichment analysis
    This example show neighbors enrichment routine. It calculates an enrichment score based on proximity on the
    connectivity graph of cell clusters. The number of observed events is compared against N permutations and a 
    z-score is computed.
    """
    # compute neighborhood enrichment
    sq.gr.nhood_enrichment(adata,
                           cluster_key="leiden",)
    # visualize
    fig, ax = plt.subplots(1,2, figsize=(13,7))
    sq.pl.nhood_enrichment(adata,
                           cluster_key="leiden",
                           figsize=(8,8),
                           title="Neighborhood enrichment adata",
                           ax=ax[0],)
    sq.pl.spatial_scatter(adata_subsample,
                          color="leiden",
                          shape=None,
                          size=2,
                          ax=ax[1])

    """
    Computation of spatial statistics - Ripley's statistics
    The Ripley's L function is a descriptive statistics function generally used to determine whether points have a 
    random, dispersed or clustered distribution pattern at certain scale. 
    The Ripley's L is a variance-normalized version of the Ripley's K statistic.
    there are two statistics:
        - 'G':  monitors the portion of points for which the nearest neighbor is within a given threshold
                and plots the cumulative percentage against the increasing distance radii.
        - 'F':  assembles the percentage of points which can be found in the aforemntioned range from an arbitrary point
                pattern spawned in the expanse of the noticed pattern.
    """
    # compute ripley L function with mode, either G or F
    fig, ax = plt.subplots(1,2,figsize=(15,7))
    mode="L"
    sq.gr.ripley(adata,
                 cluster_key="leiden",
                 mode=mode)
    sq.pl.ripley(adata,
                 cluster_key="leiden",
                 mode=mode,
                 ax=ax[0],)
    sq.pl.spatial_scatter(
        adata_subsample,
        color="leiden",
        groups=["0", "1", "26"],
        shape=None,
        size=2,
        ax=ax[1],
    )

if __name__ == "__main__":
    tutorialVizgenData(path=VIZGEN_PATH)
