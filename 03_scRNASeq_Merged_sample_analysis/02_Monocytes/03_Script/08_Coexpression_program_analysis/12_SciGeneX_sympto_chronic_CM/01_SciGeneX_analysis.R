
###############################################################################
#                           Subclustering + SciGeneX                         #
###############################################################################

## If needed to load the merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x) <- "cluster_renames"

# Subset to only include cells from Individual A and Individual B
sc10x_subset <- subset(
  sc10x,
  subset = multiplet_type == "Classical monocyte" &
           orig.ident %in% c(sample_id_A, sample_id_B)
)




# Select informative genes
res_subset <- select_genes(sc10x_subset,
                             distance = "pearson",
                             row_sum=5,
                             k = SCIGENEX_K_NEIGHBORHOOD)

# Cluster informative features
res_subset
# Create the output directory
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME))
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME))
saveRDS(res_subset, file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_clusterset_object.RDS")))
res_cluster_subtype <- readRDS(file = file.path(PATH_ANALYSIS_OUTPUT,paste0(outputFilesPrefix, "_clusterset_object.RDS")))

## Construct and partition the graph

# testing parameter 
# reciprocal_neighborhood with inflation 1.8 -> module 9 cDC1 XCR1+ 
# closest_neighborhood with inflation 2 -> no good

res_cluster_subtype <- gene_clustering(res_cluster_subtype, method= SCIGENEX_METHOD, 
                                        inflation = SCIGENEX_INFLATION, 
                                        threads = 4,
                                        output_path = NULL,
                                        s = 5)

# Filter gene modules
stat_gene_modules <- cluster_stats(res_cluster_subtype)
stat_gene_modules

# stat_gene_modules[order(stat_gene_modules$sd_total, decreasing = T),]

res_cluster_subtype <- filter_cluster_size(res_cluster_subtype, min_cluster_size = 4)
res_cluster_subtype <- filter_cluster_sd(res_cluster_subtype, min_sd = 0) #Set minimal sd

# Display the heatmap of gene clusters
res_cluster_subtype <- top_genes(res_cluster_subtype)

# Display top gene of each module
for( cluster_id in 1:length(res_cluster_subtype@gene_clusters) ){
  cat("Module SciGeneX:", as.numeric(cluster_id),"\n")
  cat(get_genes(res_cluster_subtype , as.numeric(cluster_id) ,top = T),"\n \n")
}


# Compute the a list of gene in each module and the name of each module
cluster_feature_combine <- NULL
Module_name <- NULL
for( cluster_id in 1:length(res_cluster_subtype@gene_clusters) ){
  cluster_feature <- list(get_genes(res_cluster_subtype , cluster_id))
  cluster_feature_combine <- c(cluster_feature_combine,cluster_feature)
  Module_name <- c(Module_name,paste("Module",cluster_id, sep="_"))
}

# Create au module score of each module and compute score for each cell
sc10x_subset <- AddModuleScore(sc10x_subset, features = cluster_feature_combine ,name = "Module_")

# Plot of each module and cell in the module 
p <-FeaturePlot_scCustom(sc10x_subset, features = Module_name, reduction  =  "umap", colors_use = viridis_light_high, pt.size = 0.05) & 
  theme(text = element_text(size = 10))

for(i in 1:length(Module_name)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

module_umap <- p+plot_layout(ncol = 7)
module_umap


#save plot 
ggsave(plot = module_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_modules_umap_subclusters.png"))
       , width = 23, height = 13, units ="cm")


# Plot top gene in module 

gene_in_module <- cluster_feature_combine[[2]]
selection_gene_in_module <- gene_in_module[1:24]
p_gene_module <-FeaturePlot_scCustom(sc10x_subset, features =selection_gene_in_module, reduction  =  "umap", colors_use = viridis_light_high)& 
  theme(text = element_text(size = 10))

for(i in 1:length(selection_gene_in_module)) {
  p_gene_module[[i]] <- p_gene_module[[i]] + NoLegend() + NoAxes()
}

p_gene_module

ggsave(plot = p_gene_module, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_modules_original_umap.png"))
       , width = 23, height = 13, units ="cm")



#### Module heatmap
##Display heatmap for a specific module and top genes
heatmap_one_module <- function(module_nb){

  p_heatmap <- plot_ggheatmap(res_cluster_subtype[module_nb,], ident = Seurat::Idents(sc10x_subset), use_top_genes = T) + ggtitle("All clusters (top genes)") +  theme(strip.text.y = element_text(size=4))
  # save plot
  ggsave(plot = p_heatmap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0("meaningful_modules_", outputFilesPrefix,SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_heatmap_chronic.png"))
       , width = 23, height = 13, units ="cm")

  plot_heatmap(res_cluster_subtype[module_nb,],use_top_genes = T, cell_clusters =Seurat::Idents(sc10x_subset)) #becarful heatmap interactive 
}
# Choice of meaningful modules
meaningful_modules <- c(1:37) # module numbers must be in order (same order as it will be displayed interactively on Rstudio
heatmap_one_module(meaningful_modules) # module numbers must be in order (same order as it will be displayed interactively on Rstudio




################################## Save gene list with one row per module, list of genes separated by a comma (for CellMarker2.0 handling)
library(rlist)


module_gene_list <- list()

for (module_nb in 1:length(res_cluster_subtype@gene_clusters)) {
    gene_vector <- res_cluster_subtype@gene_clusters[[module_nb]]
    gene_list <- list(paste(gene_vector, collapse = ","))
    module_gene_list <- list.append(module_gene_list, gene_list)
   }
# Convert list elements to a single string where elements are separated by commas
    module_gene_string <- sapply(module_gene_list, function(x) paste(x, collapse = ","))

write.csv(module_gene_string, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_ModuleGenes_list.csv")) , na='')

# --- Flatten all genes from all modules into a single vector ---
all_SciGenX_gene <- unlist(res_cluster_subtype@gene_clusters)  # vector of all genes

# Save CSV with one gene per row
write.csv(
  data.frame(Gene = all_SciGenX_gene),
  file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_all_SciGeneX_genes.csv")),
  row.names = FALSE,
  na = ''
)

# To save a csv for STRING handling, one column per module, one gene per row
# --- Extract module genes into a list of vectors ---
gene_module_vectors <- res_cluster_subtype@gene_clusters  # already a list of vectors

# --- Find the max module size to pad shorter columns ---
max_len <- max(sapply(gene_module_vectors, length))

# --- Pad each module vector with NA to make all equal length ---
gene_module_padded <- lapply(gene_module_vectors, function(x) {
  c(x, rep(NA, max_len - length(x)))
})

# --- Combine into a data.frame (columns = modules) ---
gene_module_df <- as.data.frame(gene_module_padded)

# Optional: rename columns
colnames(gene_module_df) <- paste0("Module_", seq_along(gene_module_vectors))

# --- Save CSV ---
write.csv(
  gene_module_df,
  file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_gene_module_in_column.csv")),
  row.names = FALSE,
  na = ""
)





















