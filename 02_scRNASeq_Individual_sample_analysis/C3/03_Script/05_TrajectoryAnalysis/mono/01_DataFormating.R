################################                                                   #######################################################
                                 # Data loading and CellDataSet object creation #
################################                                                   #######################################################

# ...............sc10x from GlobalHeterogeneity step ........................

# Loading the SeuratObject from the GlobalHeterogeneity analysis
sc10x <- readRDS(PATH_SEURAT_RDS_FILE)


# Creating a "monocle3_clusters" and "monocle3_partitions" metadata column for clustering information to be transfered from the Seurat object to CellDataSet object
sc10x$monocle3_clusters <- sc10x$seurat_clusters
sc10x$monocle3_partitions <- sc10x$seurat_clusters

# Converting the seurat object to a CellDataSet object
library(monocle3)
data_CDS <- as.cell_data_set(sc10x, import_all = TRUE)

# Creating a gene_short_name metadata column for better subsequent handling
rowData(data_CDS)$gene_short_name <- rownames(data_CDS)


# ...............sc10x_SciGeneX from the SciGeneX analysis ..................


# Loading the seurat object from the SciGeneX analysis containing expression data
sc10x_SciGeneX <- readRDS(PATH_SEURAT_RDS_SGX_FILE)


# Creating a "monocle3_clusters" and "monocle3_partitions" metadata column for clustering information to be transfered from the Seurat object to CellDataSet object
sc10x_SciGeneX$monocle3_clusters <- sc10x_SciGeneX$seurat_clusters
sc10x_SciGeneX$monocle3_partitions <- sc10x_SciGeneX$seurat_clusters

# Converting the seurat object to a CellDataSet object
data_CDS_SciGeneX <- as.cell_data_set(sc10x_SciGeneX, import_all = TRUE)

# Creating a gene_short_name metadata column for better subsequent handling
rowData(data_CDS_SciGeneX)$gene_short_name <- rownames(data_CDS_SciGeneX)

# .............. Loading the SciGeneX ClusterSet object containing modules ...............
# NOT USED
res_SciGeneX <- readRDS(PATH_SGX_CLUSTER_SET)
