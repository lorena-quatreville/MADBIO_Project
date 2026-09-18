###############################################################################
#                           Subclustering                            #
###############################################################################

## If needed to load the merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x) <- "cluster_renames"

# Subset to only include cells from Individual A and Individual B
sc10x_subset_sympto <- subset(
  sc10x,
  subset = multiplet_type == "Classical monocyte" &
           orig.ident %in% c(sample_id_A)
)


# Set cluster identity for this subset
Idents(sc10x_subset_sympto) <- sc10x_subset_sympto$orig.ident

sc10x_subset_sympto <- FindNeighbors(object    = sc10x_subset_sympto,
                                k.param   = FINDNEIGHBORS_K, 
                                reduction = "pca",
                                dims      = 1:nbPC_findclusters,
                                verbose   = .VERBOSE);

sc10x_subset_sympto <- FindClusters(object       = sc10x_subset_sympto,
                               resolution         = 0.55,
                               algorithm          = FINDCLUSTERS_ALGORITHM,
                               temp.file.location = "/tmp/",
                               verbose            = .VERBOSE);

#Parameter of SciGeneX
SCIGENEX_K_NEIGHBORHOOD = 80
SCIGENEX_METHOD         = "closest_neighborhood"  # 2 methodes : reciprocal_neighborhood or closest_neighborhood
SCIGENEX_INFLATION      = 1.8                        # number of module


# Select informative genes
res_subset_sympto <- select_genes(sc10x_subset_sympto,
                             distance = "pearson",
                             row_sum=5,
                             k = SCIGENEX_K_NEIGHBORHOOD)

# Cluster informative features
res_subset_sympto
# Create the output directory
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME))
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME))
saveRDS(res_subset_sympto, file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_clusterset_object_sympto.RDS")))
res_cluster_subtype_sympto <- readRDS(file = file.path(PATH_ANALYSIS_OUTPUT,paste0(outputFilesPrefix, "_clusterset_object_sympto.RDS")))

## Construct and partition the graph

# testing parameter 
# reciprocal_neighborhood with inflation 1.8 -> module 9 cDC1 XCR1+ 
# closest_neighborhood with inflation 2 -> no good

res_cluster_subtype_sympto <- gene_clustering(res_cluster_subtype_sympto, method= SCIGENEX_METHOD, 
                                        inflation = SCIGENEX_INFLATION, 
                                        threads = 4,
                                        output_path = NULL,
                                        s = 5)

# Filter gene modules
stat_gene_modules <- cluster_stats(res_cluster_subtype_sympto)
stat_gene_modules

# stat_gene_modules[order(stat_gene_modules$sd_total, decreasing = T),]

res_cluster_subtype_sympto <- filter_cluster_size(res_cluster_subtype_sympto, min_cluster_size = 4)
res_cluster_subtype_sympto <- filter_cluster_sd(res_cluster_subtype_sympto, min_sd = 0) #Set minimal sd

# Display the heatmap of gene clusters
res_cluster_subtype_sympto <- top_genes(res_cluster_subtype_sympto)

# Display top gene of each module
for( cluster_id in 1:length(res_cluster_subtype_sympto@gene_clusters) ){
  cat("Module SciGeneX:", as.numeric(cluster_id),"\n")
  cat(get_genes(res_cluster_subtype_sympto , as.numeric(cluster_id) ,top = T),"\n \n")
}


# Compute the a list of gene in each module and the name of each module
cluster_feature_combine <- NULL
Module_name <- NULL
for( cluster_id in 1:length(res_cluster_subtype_sympto@gene_clusters) ){
  cluster_feature <- list(get_genes(res_cluster_subtype_sympto , cluster_id))
  cluster_feature_combine <- c(cluster_feature_combine,cluster_feature)
  Module_name <- c(Module_name,paste("Module",cluster_id, sep="_"))
}

# Create au module score of each module and compute score for each cell
sc10x_subset_sympto <- AddModuleScore(sc10x_subset_sympto, features = cluster_feature_combine ,name = "Module_")

# Plot of each module and cell in the module 
p <-FeaturePlot_scCustom(sc10x_subset_sympto, features = Module_name, reduction  =  "umap", colors_use = viridis_light_high, pt.size = 0.05) & 
  theme(text = element_text(size = 10))

for(i in 1:length(Module_name)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

module_umap <- p+plot_layout(ncol = 7)
module_umap


#save plot 
ggsave(plot = module_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_modules_umap_subclusters_sympto.png"))
       , width = 23, height = 13, units ="cm")


# Plot top gene in module 

gene_in_module <- cluster_feature_combine[[2]]
selection_gene_in_module <- gene_in_module[1:24]
p_gene_module <-FeaturePlot_scCustom(sc10x_subset_sympto, features =selection_gene_in_module, reduction  =  "umap", colors_use = viridis_light_high)& 
  theme(text = element_text(size = 10))

for(i in 1:length(selection_gene_in_module)) {
  p_gene_module[[i]] <- p_gene_module[[i]] + NoLegend() + NoAxes()
}

p_gene_module

ggsave(plot = p_gene_module, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_modules_original_umap_sympto.png"))
       , width = 23, height = 13, units ="cm")

################################## Save gene list with one row per module, list of genes separated by a comma (for CellMarker2.0 handling)
library(rlist)


module_gene_list <- list()

for (module_nb in 1:length(res_cluster_subtype_sympto@gene_clusters)) {
    gene_vector <- res_cluster_subtype_sympto@gene_clusters[[module_nb]]
    gene_list <- list(paste(gene_vector, collapse = ","))
    module_gene_list <- list.append(module_gene_list, gene_list)
   }
# Convert list elements to a single string where elements are separated by commas
    module_gene_string <- sapply(module_gene_list, function(x) paste(x, collapse = ","))

write.csv(module_gene_string, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_ModuleGenes_list_sympto.csv")) , na='')



################################ Plotting cell groups isolated by the subclustering on UMAP

umap_subclusters <-  DimPlot( sc10x_subset_sympto, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "seurat_clusters",
           label = TRUE) + 
    ggtitle( paste0("Map of cells with subclustering "))
print(umap_subclusters)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "_DimPlot_sub_clusters_sympto.png")), plot = umap_subclusters, device ="png", width = 23, height = 13, units ="cm")



# Plot the UMAP grouped by sample of origin

# Plot UMAP with clusters grouped by sample of origin
plot_sample_orig <- DimPlot(sc10x_subset_sympto, group.by = "orig.ident") +
  ggtitle("UMAP with Cluster Names") +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank())
  print(plot_sample_orig)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "DimPlot_grouped_sample_orig_sympto.png")), plot = plot_sample_orig, device ="png",  width = 12, height =8, dpi = 300)



##Display heatmap for a specific module and top genes
heatmap_one_module <- function(module_nb){

  p_heatmap <- plot_ggheatmap(res_cluster_subtype_sympto[module_nb,], ident = Seurat::Idents(sc10x_subset_sympto), use_top_genes = T) + ggtitle("All clusters (top genes)") +  theme(strip.text.y = element_text(size=4))
  # save plot
  ggsave(plot = p_heatmap, filename = file.path(PATH_ANALYSIS_OUTPUT,  paste0("meaningful_modules_", outputFilesPrefix,SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_heatmap_sympto.png"))
       , width = 23, height = 13, units ="cm")

  plot_heatmap(res_cluster_subtype_sympto[module_nb,],use_top_genes = T, cell_clusters =Seurat::Idents(sc10x_subset_sympto)) #becarful heatmap interactive 
}
# Choice of meaningful modules
meaningful_modules <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) # module numbers must be in order (same order as it will be displayed interactively on Rstudio
heatmap_one_module(meaningful_modules) # module numbers must be in order (same order as it will be displayed interactively on Rstudio


## Get the percentage of cells belonging to each individual in each cluster

# Extract metadata
metadata <- sc10x_subset_sympto@meta.data

# Count cells per cluster and individual
cluster_composition <- metadata %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  ungroup()

# Calculate percentage of each individual in each cluster
cluster_composition <- cluster_composition %>%
  group_by(seurat_clusters) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  select(seurat_clusters, orig.ident, percent) %>%
  pivot_wider(names_from = orig.ident, values_from = percent, values_fill = 0)  # Wide format for better readability

# Write to a TSV file
write.table(cluster_composition, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cluster_compo_indiv_sympto.csv")), sep = "\t", row.names = FALSE, quote = FALSE)














                                ####### Chronic individuals ################
                                
                                
                                
                                

# Subset to only include cells from Individual A and Individual B
sc10x_subset_chronic <- subset(
  sc10x,
  subset = multiplet_type == "Classical monocyte" &
           orig.ident %in% c(sample_id_B)
)


# Set cluster identity for this subset
Idents(sc10x_subset_chronic) <- sc10x_subset_chronic$orig.ident

sc10x_subset_chronic <- FindNeighbors(object    = sc10x_subset_chronic,
                                k.param   = FINDNEIGHBORS_K, 
                                reduction = "pca",
                                dims      = 1:nbPC_findclusters,
                                verbose   = .VERBOSE);

sc10x_subset_chronic <- FindClusters(object       = sc10x_subset_chronic,
                               resolution         = 0.55,
                               algorithm          = FINDCLUSTERS_ALGORITHM,
                               temp.file.location = "/tmp/",
                               verbose            = .VERBOSE);

#Parameter of SciGeneX
SCIGENEX_K_NEIGHBORHOOD = 80
SCIGENEX_METHOD         = "closest_neighborhood"  # 2 methodes : reciprocal_neighborhood or closest_neighborhood
SCIGENEX_INFLATION      = 1.8                        # number of module


# Select informative genes
res_subset_chronic <- select_genes(sc10x_subset_chronic,
                             distance = "pearson",
                             row_sum=5,
                             k = SCIGENEX_K_NEIGHBORHOOD)

# Cluster informative features
res_subset_chronic
# Create the output directory
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME))
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME))
saveRDS(res_subset_chronic, file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_clusterset_object_chronic.RDS")))
res_cluster_subtype_chronic <- readRDS(file = file.path(PATH_ANALYSIS_OUTPUT,paste0(outputFilesPrefix, "_clusterset_object_chronic.RDS")))

## Construct and partition the graph

# testing parameter 
# reciprocal_neighborhood with inflation 1.8 -> module 9 cDC1 XCR1+ 
# closest_neighborhood with inflation 2 -> no good

res_cluster_subtype_chronic <- gene_clustering(res_cluster_subtype_chronic, method= SCIGENEX_METHOD, 
                                        inflation = SCIGENEX_INFLATION, 
                                        threads = 4,
                                        output_path = NULL,
                                        s = 5)

# Filter gene modules
stat_gene_modules <- cluster_stats(res_cluster_subtype_chronic)
stat_gene_modules

# stat_gene_modules[order(stat_gene_modules$sd_total, decreasing = T),]

res_cluster_subtype_chronic <- filter_cluster_size(res_cluster_subtype_chronic, min_cluster_size = 4)
res_cluster_subtype_chronic <- filter_cluster_sd(res_cluster_subtype_chronic, min_sd = 0) #Set minimal sd

# Display the heatmap of gene clusters
res_cluster_subtype_chronic <- top_genes(res_cluster_subtype_chronic)

# Display top gene of each module
for( cluster_id in 1:length(res_cluster_subtype_chronic@gene_clusters) ){
  cat("Module SciGeneX:", as.numeric(cluster_id),"\n")
  cat(get_genes(res_cluster_subtype_chronic , as.numeric(cluster_id) ,top = T),"\n \n")
}


# Compute the a list of gene in each module and the name of each module
cluster_feature_combine <- NULL
Module_name <- NULL
for( cluster_id in 1:length(res_cluster_subtype_chronic@gene_clusters) ){
  cluster_feature <- list(get_genes(res_cluster_subtype_chronic , cluster_id))
  cluster_feature_combine <- c(cluster_feature_combine,cluster_feature)
  Module_name <- c(Module_name,paste("Module",cluster_id, sep="_"))
}

# Create au module score of each module and compute score for each cell
sc10x_subset_chronic <- AddModuleScore(sc10x_subset_chronic, features = cluster_feature_combine ,name = "Module_")

# Plot of each module and cell in the module 
p <-FeaturePlot_scCustom(sc10x_subset_chronic, features = Module_name, reduction  =  "umap", colors_use = viridis_light_high, pt.size = 0.05) & 
  theme(text = element_text(size = 10))

for(i in 1:length(Module_name)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

module_umap <- p+plot_layout(ncol = 7)
module_umap


#save plot 
ggsave(plot = module_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_modules_umap_subclusters_chronic.png"))
       , width = 23, height = 13, units ="cm")


# Plot top gene in module 

gene_in_module <- cluster_feature_combine[[2]]
selection_gene_in_module <- gene_in_module[1:24]
p_gene_module <-FeaturePlot_scCustom(sc10x_subset_chronic, features =selection_gene_in_module, reduction  =  "umap", colors_use = viridis_light_high)& 
  theme(text = element_text(size = 10))

for(i in 1:length(selection_gene_in_module)) {
  p_gene_module[[i]] <- p_gene_module[[i]] + NoLegend() + NoAxes()
}

p_gene_module

ggsave(plot = p_gene_module, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_modules_original_umap_chronic.png"))
       , width = 23, height = 13, units ="cm")

################################## Save gene list with one row per module, list of genes separated by a comma (for CellMarker2.0 handling)
library(rlist)


module_gene_list <- list()

for (module_nb in 1:length(res_cluster_subtype_chronic@gene_clusters)) {
    gene_vector <- res_cluster_subtype_chronic@gene_clusters[[module_nb]]
    gene_list <- list(paste(gene_vector, collapse = ","))
    module_gene_list <- list.append(module_gene_list, gene_list)
   }
# Convert list elements to a single string where elements are separated by commas
    module_gene_string <- sapply(module_gene_list, function(x) paste(x, collapse = ","))

write.csv(module_gene_string, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_ModuleGenes_list_chronic.csv")) , na='')



################################ Plotting cell groups isolated by the subclustering on UMAP

umap_subclusters <-  DimPlot( sc10x_subset_chronic, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "seurat_clusters",
           label = TRUE) + 
    ggtitle( paste0("Map of cells with subclustering "))
print(umap_subclusters)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "_DimPlot_sub_clusters_chronic.png")), plot = umap_subclusters, device ="png", width = 23, height = 13, units ="cm")



# Plot UMAP with clusters grouped by sample of origin
plot_sample_orig <- DimPlot(sc10x_subset_chronic, group.by = "orig.ident") +
  ggtitle("UMAP with Cluster Names") +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank())
  print(plot_sample_orig)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "DimPlot_grouped_sample_orig_chronic.png")), plot = plot_sample_orig, device ="png",  width = 12, height =8, dpi = 300)



##Display heatmap for a specific module and top genes
heatmap_one_module <- function(module_nb){

  p_heatmap <- plot_ggheatmap(res_cluster_subtype_chronic[module_nb,], ident = Seurat::Idents(sc10x_subset_chronic), use_top_genes = T) + ggtitle("All clusters (top genes)") +  theme(strip.text.y = element_text(size=4))
  # save plot
  ggsave(plot = p_heatmap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0("meaningful_modules_", outputFilesPrefix,SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_heatmap_chronic.png"))
       , width = 23, height = 13, units ="cm")

  plot_heatmap(res_cluster_subtype_chronic[module_nb,],use_top_genes = T, cell_clusters =Seurat::Idents(sc10x_subset_chronic)) #becarful heatmap interactive 
}
# Choice of meaningful modules
meaningful_modules <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10) # module numbers must be in order (same order as it will be displayed interactively on Rstudio
heatmap_one_module(meaningful_modules) # module numbers must be in order (same order as it will be displayed interactively on Rstudio



## Get the percentage of cells belonging to each individual in each cluster

# Extract metadata
metadata <- sc10x_subset_chronic@meta.data

# Count cells per cluster and individual
cluster_composition <- metadata %>%
  group_by(seurat_clusters, orig.ident) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  ungroup()

# Calculate percentage of each individual in each cluster
cluster_composition <- cluster_composition %>%
  group_by(seurat_clusters) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  select(seurat_clusters, orig.ident, percent) %>%
  pivot_wider(names_from = orig.ident, values_from = percent, values_fill = 0)  # Wide format for better readability

# Write to a TSV file
write.table(cluster_composition, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cluster_compo_indiv_chronic.csv")), sep = "\t", row.names = FALSE, quote = FALSE)







# Building a global heatmap of all modules accross all cells

###############################################################################
# Extract symptomatic modules
###############################################################################

modules_sympto <- list()

for(cluster_id in seq_len(length(res_cluster_subtype_sympto@gene_clusters))){

  modules_sympto[[paste0("Sympto_M", cluster_id)]] <-
    get_genes(
      res_cluster_subtype_sympto,
      cluster_id,
      top = TRUE
    )

}

###############################################################################
# Extract chronic modules
###############################################################################

modules_chronic <- list()

for(cluster_id in seq_len(length(res_cluster_subtype_chronic@gene_clusters))){

  modules_chronic[[paste0("Chronic_M", cluster_id)]] <-
    get_genes(
      res_cluster_subtype_chronic,
      cluster_id,
      top = TRUE
    )

}

###############################################################################
# Merge all modules
###############################################################################

all_modules <- c(
  modules_sympto,
  modules_chronic
)

# Reloading the original seurat object with both groups

sc10x <- subset(
  sc10x,
  subset = multiplet_type == "Classical monocyte"
)

# Computing the module score on all cells
sc10x <- AddModuleScore(
  object = sc10x,
  features = all_modules,
  name = "SciGeneX_"
)

# Renaiming module metadatacolumn with identifyers
score_cols <- grep(
  "^SciGeneX_",
  colnames(sc10x@meta.data),
  value = TRUE
)

new_names <- names(all_modules)

colnames(sc10x@meta.data)[
  match(score_cols, colnames(sc10x@meta.data))
] <- new_names

# Creating condition labels
sc10x$condition <- ifelse(
  sc10x$orig.ident %in% sample_id_A,
  "Symptomatic",
  "Chronic"
)

# Building a module score matriw for the plotting
module_names <- names(all_modules)

score_matrix <- t(
  as.matrix(
    sc10x@meta.data[, module_names]
  )
)

# Ordering cells by condition/individuals
cells_sympto <- rownames(sc10x@meta.data)[
  sc10x$condition == "Symptomatic"
]

cells_chronic <- rownames(sc10x@meta.data)[
  sc10x$condition == "Chronic"
]

cell_order <- c(
  cells_sympto,
  cells_chronic
)

score_matrix <- score_matrix[, cell_order]

# Plotting the heatmap
library(ComplexHeatmap)
library(circlize)

global_heatmap = Heatmap(
  score_matrix,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  show_column_names = FALSE,

  row_split = c(
    rep("Symptomatic modules", length(modules_sympto)),
    rep("Chronic modules", length(modules_chronic))
  ),

  column_split = c(
    rep("Symptomatic", length(cells_sympto)),
    rep("Chronic", length(cells_chronic))
  ),

  col = colorRamp2(
    c(-2,0,2),
    c("blue","white","red")
  ),

  name = "Score"
)

  # Saving the heatmap
  png(file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix,
"sympto_chronic_heatmap.png")),
  width = 3000,
  height = 2500,
  res = 300
)

draw(global_heatmap)

dev.off()






























#### To plot one module
umap_one_module <- function(module_nb){
  module_name_list = Module_name
  module_name = module_name_list[module_nb]
  p <-FeaturePlot_scCustom(sc10x, features = module_name, reduction  =  "umap", colors_use = viridis_light_high,pt.size = 0.5) & 
  theme(text = element_text(size = 10))

  p[module_name] <- p[module_name] + NoLegend() + NoAxes()

  module_nb_umap <- p+plot_layout(ncol = 1)
  module_nb_umap
  ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, module_nb, "FeaturePlot.png")), plot = p, device ="png", width = 23, height = 13, units ="cm")

}









