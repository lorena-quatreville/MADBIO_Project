# #########################################
# This script reads and filters sc10x  data
# #########################################


## .................................................................
## CREATE THE OUTPUT DIRECTORIES
## .................................................................

## @knitr createOutputDir

### Create output directory
dir.create( PATH_ANALYSIS_OUTPUT, showWarnings = FALSE, recursive = TRUE)
PATH_ANALYSIS_EXTRA_OUTPUT = file.path( PATH_ANALYSIS_OUTPUT, "Extra")
dir.create( PATH_ANALYSIS_EXTRA_OUTPUT, showWarnings = FALSE, recursive = TRUE)

## .................................................................
## LOAD AND PREPARE THE DATA
## .................................................................

## @knitr loadData
sc10x_orig = sc10x
SELECTED_CLUSTERS = c("1","2", "3", "4", "5", "6", "7")

# Read the seurat object from previous analysis
PATH_SEURAT_OBJECT = file.path(PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "_seuratObject_unsupervised_clustering.RDS"))
sc10x = readRDS( PATH_SEURAT_OBJECT)

# Plot the data before and after selection

cat("\n \n")
cat("### Data before cluster selection")
Idents( sc10x) = "heatmap_clusters"
DimPlot( sc10x)

cat("\n \n")
cat("### Selected clusters", paste( SELECTED_CLUSTERS, collapse = "; "))
Idents( sc10x) = "heatmap_clusters"
sc10x = subset( sc10x,idents = SELECTED_CLUSTERS)
sc10x = AddMetaData( sc10x, metadata = Idents( sc10x), col.name = "orig.seurat_clusters")
DimPlot( sc10x)

cat("\n \n")
cat("### SciGeneX analysis")

# Select informative genes
set_verbosity(3)

# Selecting genes and clustering
scigx <- select_genes(
    data = sc10x,
    row_sum=5,
    k = SCIGENEX_K_NEIGHBORHOOD
)

scigx <- gene_clustering(
    scigx,
    method = SCIGENEX_METHOD,
    inflation = SCIGENEX_INFLATION
)


# Filter the modules
scigx <- filter_cluster_size(scigx, min_cluster_size = 10)
scigx <- filter_cluster_sd(scigx, min_sd = 0.15)


# SAVING GENE MODULE LISTS
module_gene_list <- list()

for (module_nb in 1:length(scigx@gene_clusters)) {
    gene_vector <- scigx@gene_clusters[[module_nb]]
    gene_list <- list(paste(gene_vector, collapse = ","))
    module_gene_list <- list.append(module_gene_list, gene_list)
   }
# Convert list elements to a single string where elements are separated by commas
    module_gene_string <- sapply(module_gene_list, function(x) paste(x, collapse = ","))

write.csv(module_gene_string, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_ModuleGenes_list_module_activity_scoring.csv")) , na='')



# Select the top genes per module
scigx_res <- top_genes(scigx)

# Plot the heatmap of cells and gene modules
Idents(sc10x) = "heatmap_clusters"
plot_heatmap(scigx_res, cell_clusters = Seurat::Idents(sc10x), interactive = FALSE, row_labels = FALSE, show_legend = FALSE)

# Compute the UCell score for each cell and each module
sgx_module_ordered_features_set = c()
sgx_module_feature_list = list()
sgx_module_name_set = c()
for( cluster in names( scigx_res@top_genes)){
  module_name = paste0( "Module", cluster)
  sgx_module_name_set = c( sgx_module_name_set, module_name)
  feature_list = list( unlist( scigx_res@top_genes[ cluster], use.names = FALSE))
  names( feature_list) = module_name
  sgx_module_feature_list[[ module_name]] = unlist( scigx_res@top_genes[ cluster], use.names = FALSE)
  sc10x = UCell::AddModuleScore_UCell( sc10x, features = feature_list, name = "")
  sgx_module_ordered_features_set = append( sgx_module_ordered_features_set, unlist( scigx_res@top_genes[ cluster], use.names = FALSE))
}

pheatmap( as.matrix(sc10x@meta.data[ sgx_module_name_set]),
          scale = "column",
          show_rownames = FALSE)

# ###########################################
# Compute activity of gene modules
# ###########################################

# Get the UCell score for each module
ordered_data = scigx_res@data[ sgx_module_ordered_features_set, ]

gene_means = rowMeans( ordered_data)
gene_sd = matrixStats::rowSds( ordered_data)

module_activity_ordered_data = data.frame()
for( col_index in 1:ncol( ordered_data)){
  module_activity_ordered_data = rbind( module_activity_ordered_data, as.numeric( ordered_data[ ,col_index] > gene_means))
}

module_activity_ordered_data = t( module_activity_ordered_data)
colnames( module_activity_ordered_data) = colnames( ordered_data)
rownames( module_activity_ordered_data) = rownames( ordered_data)

# Compute the module activity
# The activity is based on a score which is the mean of the 
module_activity_score_names = c()
module_activity_names = c()
for( cluster in names( scigx_res@top_genes)){
  module_name = paste0( "Module", cluster)
  feature_set = unlist( scigx_res@top_genes[ cluster], use.names = FALSE)
  cat( "\n", module_name, "=", length( feature_set))
  score = colSums( module_activity_ordered_data[ feature_set, ])/length( feature_set)
  score_bin = as.factor( as.numeric( score >= 0.5))
  names( score_bin) = names( score)
  score_name =  paste0( "Activity.Score.SGX", module_name)
  score_bin_name =  paste0( "Activity.SGX", module_name)
  module_activity_score_names = c( module_activity_score_names, score_name)
  module_activity_names = c( module_activity_names, score_bin_name)
  sc10x = AddMetaData( sc10x, metadata = score[ Cells( sc10x)], col.name = score_name)
  sc10x = AddMetaData( sc10x, metadata = score_bin[ Cells( sc10x)], col.name = score_bin_name)
}

# module_activity_names_df = sc10x@meta.data[ module_activity_names]
# pheatmap( as.numeric( t(module_activity_names_df[ do.call(order, as.list( module_activity_names_df)), ])),  cluster_rows = FALSE, cluster_cols = FALSE, show_colnames = FALSE, scale = "none")


# Visualization

# Annotation column for clusters
cell_clusters <- sc10x@meta.data$heatmap_clusters  # change if needed
annotation_col <- data.frame(
  Cluster = factor(cell_clusters)
)
rownames(annotation_col) <- rownames(sc10x@meta.data)


module_cols <- grep("^Activity\\.SGXModule", colnames(sc10x@meta.data), value = TRUE)


mat_bin <-
as.matrix(
  sc10x@meta.data[
      ,
      module_cols,
      drop=FALSE
  ]
)

mode(mat_bin) <- "numeric"

rownames(mat_bin) <-
rownames(sc10x@meta.data)

mat_bin <- t(mat_bin)


png(file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "module_activity_heatmap.png")),
    width = 2500,
    height = 2000,
    res = 300)

pheatmap(
  mat_bin,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = c("white", "red"),
  breaks = c(-0.01, 0.5, 1.01)
)

dev.off()











# Same but cells ordered by cluster

# -----------------------------
# 1. Select module columns
# -----------------------------
module_cols <- grep("^Activity\\.SGXModule", colnames(sc10x@meta.data), value = TRUE)

# -----------------------------
# 2. Build binary matrix (modules x cells)
# -----------------------------
mat_bin <-
as.matrix(
  sc10x@meta.data[
      ,
      module_cols,
      drop=FALSE
  ]
)

mode(mat_bin) <- "numeric"

rownames(mat_bin) <-
rownames(sc10x@meta.data)

mat_bin <- t(mat_bin)

# -----------------------------
# 3. Add cluster annotation
# -----------------------------

cell_order <- hc_cols$order

ann_col <- data.frame(
    Cluster = factor(
        sc10x$heatmap_clusters[cell_order]
    )
)

rownames(ann_col) <- colnames(mat_bin)

# -----------------------------
# ADD COLORS (IMPORTANT FIX) (=annotation_colors)
# -----------------------------
cluster_ordered <- sc10x$heatmap_clusters[cell_order]
cluster_factor <- factor(cluster_ordered, levels = unique(cluster_ordered))

n_clusters <- length(unique(cluster_factor))

cluster_colors <- setNames(
    colorRampPalette(brewer.pal(8, "Set2"))(n_clusters),
    levels(cluster_factor)
)

ann_colors <- list(Cluster = cluster_colors)
# -----------------------------
# 4. Plot heatmap
# -----------------------------
png(file.path(PATH_ANALYSIS_OUTPUT,
              paste0(outputFilesPrefix, "_module_activity_heatmap_ordered_by_cluster.png")),
    width = 2500,
    height = 2000,
    res = 300)

pheatmap(
  mat_bin,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = ann_col,
  show_colnames = FALSE,
  color = c("white", "red"),
  breaks = c(-0.01, 0.5, 1.01)
)

dev.off()

