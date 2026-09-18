

cellType = "effector_CD4_T_cell"

###############################################################################
#                           Subclustering CD4_T_cells                           #
###############################################################################

# open the seurat object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)    

Idents(sc10x) <- "seurat_clusters"

sc10x_subtype <- subset(sc10x, idents = c(2))


                                        ## Clustering on the subcluster extracted ##
                                        
                                        
sc10x_subtype <- FindNeighbors(object    = sc10x_subtype,
                                k.param   = FINDNEIGHBORS_K, 
                                reduction = "pca",
                                dims      = 1:nbPC_findclusters,
                                verbose   = .VERBOSE);

sc10x_subtype <- FindClusters(object       = sc10x_subtype,
                               resolution         = 0.55,
                               algorithm          = FINDCLUSTERS_ALGORITHM,
                               temp.file.location = "/tmp/",
                               verbose            = .VERBOSE);

################################ Plotting cell groups isolated by the subclustering on UMAP

umap_subclusters <-  DimPlot( sc10x_subtype, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "seurat_clusters",
           label = TRUE) + 
    ggtitle( paste0("Map of cells with subclustering ", cellType))
print(umap_subclusters)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, cellType, "_DimPlot_sub_clusters.png")), plot = umap_subclusters, device ="png", width = 23, height = 13, units ="cm")


################################ Plotting cell groups showing to which original cluster they belonged to
umap_original_clusters <-  DimPlot( sc10x_subtype, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "cluster_renames",
           label = TRUE) + 
    ggtitle( paste0("Map of cells with their original cluster ", cellType))
print(umap_original_clusters)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, cellType, "_DimPlot_cluster_of_origin.png")), plot = umap_original_clusters, device ="png", width = 23, height = 13, units ="cm")


################################ Identifying doublet clusters

# Plot the presence of doublets (by DoubletFinder)

doublets <- DimPlot( sc10x_subtype, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "multiplet_prediction",
           label = TRUE) + 
    ggtitle( "Map of doublet cells")
print(doublets)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, "DimPlot_doublets.png")), plot = doublets, device ="png", width = 23, height = 13, units ="cm")
    
# Plot the cells percentage of mitochondrial genes in UMAP  
perc_mitok <- FeaturePlot( sc10x_subtype, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "percent.mito") + 
  scale_colour_gradientn(colors = c(low = "#D0E1F9", high = "#053061")) +
  ggtitle( "Map of cells with level of percentage of\n mitochondrial genes")
print(perc_mitok)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, "FeaturePlot_perc_mitok.png")), plot = perc_mitok, device ="png", width = 23, height = 13, units ="cm")

# Plot the cells RNA counts in UMAP
rna_count <- FeaturePlot( sc10x_subtype, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nCount_RNA") + 
  scale_colour_gradientn(colors = c(low = "#D0E1F9", high = "#053061")) +
  ggtitle( "Map of cells with RNA counts")
 print(rna_count)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, "FeaturePlot_rna_count.png")), plot = rna_count, device ="png", width = 23, height = 13, units ="cm")

# Plot the cells genes counts in UMAP
gene_count <- FeaturePlot( sc10x_subtype, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nFeature_RNA") + 
  scale_colour_gradientn(colors = c(low = "#D0E1F9", high = "#053061")) +
  ggtitle( "Map of cells with \nnumber of detected genes")
print(gene_count)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, "FeaturePlot_gene_count.png")), plot = gene_count, device ="png", width = 23, height = 13, units ="cm")



                                                           ## Running SciGeneX on the subcluster ##
                                                           
                                                           
# Select informative genes
res_subtype <- select_genes(sc10x_subtype,
                             distance = "pearson",
                             row_sum=5,
                             k = SCIGENEX_K_NEIGHBORHOOD)

# Cluster informative features
res_subtype
dir.create(file.path(PATH_ANALYSIS_OUTPUT))
dir.create(file.path(PATH_ANALYSIS_OUTPUT, cellType))
# save the ClusterSet object generated by "select_gene"
saveRDS(res_subtype, file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "_clusterset_object.RDS")))
res_cluster_subtype <- readRDS(file = file.path(PATH_ANALYSIS_OUTPUT, cellType ,paste0(outputFilesPrefix, cellType, "_clusterset_object.RDS")))
# save the SeuratObject with subclustering
saveRDS(sc10x_subtype, file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "_seuratObject_final.RDS")))




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


# Save gene list with one row per module, list of genes separated by a comma (for CellMarker2.0 handling)
library(rlist)


module_gene_list <- list()

for (module_nb in 1:length(res_cluster_subtype@gene_clusters)) {
    gene_vector <- res_cluster_subtype@gene_clusters[[module_nb]]
    gene_list <- list(paste(gene_vector, collapse = ","))
    module_gene_list <- list.append(module_gene_list, gene_list)
   }
# Convert list elements to a single string where elements are separated by commas
    module_gene_string <- sapply(module_gene_list, function(x) paste(x, collapse = ","))

write.csv(module_gene_string, file= file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, "_ModuleGenes_list.csv")) , na='')


# Creating a module score for each module for visualization

res_cluster_subtype <- top_genes(res_cluster_subtype) # creates a metadata column with 20 top genes for each module

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
sc10x_subtype <- AddModuleScore(sc10x_subtype, features = cluster_feature_combine ,name = "Module_")



# Plot of each module and cell in the module (UMAP)

    # Loop over all modules in steps of 35
    for (mod in seq(1, length(Module_name), by = 35)) {
    
    # Define the start and end index for this chunk
    start_idx <- mod
    end_idx <- min(mod + 34, length(Module_name))  # Ensure it doesn't exceed max
    
    # Subset module names for this group
    Module_name_sub <- paste0("Module_", start_idx:end_idx)
    
    # Generate FeaturePlot list
    feature_plot_list <- FeaturePlot_scCustom(sc10x_subtype, features = Module_name_sub, 
                            reduction = "umap", 
                            colors_use = viridis_light_high, 
                            pt.size = 0.05) & 
       theme(text = element_text(size = 10))
    
    for(i in 1:length(Module_name_sub)) {
    feature_plot_list[[i]] <- feature_plot_list[[i]] + NoLegend() + NoAxes()
  }
    
    # Combine into layout
    module_umap <- feature_plot_list + plot_layout(ncol = 7)
    module_umap
    
    # Save to file
    filename <- paste0(outputFilesPrefix, cellType, "_", SCIGENEX_METHOD, 
                     "_inflation_", SCIGENEX_INFLATION, 
                     "_modules_", start_idx, "_to_", end_idx, "_umap_subclusters.png")
    
    ggsave(plot = module_umap, 
         filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, filename),
         width = 23, height = 13, units = "cm")
    
 }



# Plot top gene in the first module (UMAP)

gene_in_module <- cluster_feature_combine[[1]] # it extracts genes from module 1
selection_gene_in_module <- gene_in_module[1:24] # it selects the first top 24 genes in this module

p_gene_module <-FeaturePlot_scCustom(sc10x_subtype, features =selection_gene_in_module, reduction  =  "umap", colors_use = viridis_light_high)& 
  theme(text = element_text(size = 10))

for(i in 1:length(selection_gene_in_module)) {
  p_gene_module[[i]] <- p_gene_module[[i]] + NoLegend() + NoAxes()
}

p_gene_module

ggsave(plot = p_gene_module, filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_top_genes_module_1_umap.png"))
       , width = 23, height = 13, units ="cm")



##Display heatmap for a specific module and top genes
heatmap_one_module <- function(module_nb){

  p_heatmap <- plot_ggheatmap(res_cluster_subtype[module_nb,], ident = Seurat::Idents(sc10x_subtype), use_top_genes = T) + ggtitle("All clusters (top genes)") +  theme(strip.text.y = element_text(size=4))
  # save plot
  ggsave(plot = p_heatmap, filename = file.path(PATH_ANALYSIS_OUTPUT,  cellType, paste0("meaningful_modules_", outputFilesPrefix,SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_heatmap.png"))
       , width = 23, height = 13, units ="cm")

  plot_heatmap(res_cluster_subtype[module_nb,],use_top_genes = T, cell_clusters =Seurat::Idents(sc10x_subtype)) #becarful heatmap interactive 
}

# Choice of meaningful modules
meaningful_modules <- c(1,6, 10, 41, 44, 94) # module numbers must be in order (same order as it will be displayed interactively on Rstudio
heatmap_one_module(meaningful_modules) # module numbers must be in order (same order as it will be displayed interactively on Rstudio







## Module of interest characterization

#Biological process enrichment analysis
res_cluster_subtype_enriched <- enrich_go(res_cluster_subtype[meaningful_modules], specie = "Hsapiens", ontology = "BP") #enrich_go overwrites the 


plot_BP_enrich <- plot_clust_enrichments(res_cluster_subtype_enriched, gradient_palette=colors_for_gradient("Je1"), 
                       floor=50,
                       nb_go_term = 2) + 
  theme(axis.text.y = element_text(size=8))
print(plot_BP_enrich)


# Save plot as PNG
ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_plot_BP_enrich_analysis_selected_modules.png")), plot = plot_BP_enrich, width = 10, height = 6, units = "in")

# Mapping the cell population markers onto the gene modules # to modify
library(tidyr)
sctype <- "https://zenodo.org/record/8269433/files/sctype.app.tsv"
marker_table <- read.table(sctype, head=TRUE, sep="\t")
marker_table <- marker_table %>% 
                filter(Tissue == "Immune system") %>% 
                filter(!grepl('Pro-|Pre-|HSC|precursor|Progenitor', Cell.type)) %>%
                separate_rows(Marker_genes, convert = TRUE)

markers <- split(marker_table$Marker_genes, 
                 marker_table$Cell.type)
plot_module_cellType <- plot_markers_to_clusters(res_cluster_subtype_enriched, 
                         markers=markers, background = rownames(sc10x_subtype ))

ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_plot_modules_cellType_identification.png")), plot = plot_module_cellType, width = 10, height = 6, units = "in")





                                                ## Refining the clustering using SciGeneX selected genes ##

# Drop previous scaled data
sc10x_subtype_SciGeneX <- sc10x_subtype # create a new seurat object before removing the scaled expression values
sc10x_subtype_SciGeneX$seurat_original_clusters <- sc10x_subtype_SciGeneX$seurat_clusters # saving the original seurat_clusers into a new metadata column to be able to compare with new clusters afterwards
DefaultAssay(sc10x_subtype_SciGeneX) <- "RNA"  # just to be sure
sc10x_subtype_SciGeneX[["RNA"]]@scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)


# Scaling of the SciGeneX genes
sc10x_subtype_SciGeneX <- ScaleData(sc10x_subtype_SciGeneX, features = row_names(res_cluster_subtype), verbose = FALSE)


# Performing PCA on the SciGeneX genes that have been rescaled
sc10x_subtype_SciGeneX <- RunPCA(sc10x_subtype_SciGeneX, features = row_names(res_cluster_subtype), verbose = FALSE)

# Visualizing PCA results # A TESTER
print(sc10x_subtype_SciGeneX[["pca"]], dims = 1:5, nfeatures = 10)
DimHeatmap(sc10x_subtype_SciGeneX, dims = 1, cells = 500, balanced = TRUE)
# Data dimensionality
ElbowPlot(sc10x_subtype_SciGeneX)


# Cell clustering
sc10x_subtype_SciGeneX <- FindNeighbors(sc10x_subtype_SciGeneX, dims = 1:10, verbose = FALSE)
sc10x_subtype_SciGeneX <- FindClusters(sc10x_subtype_SciGeneX, resolution = 0.5, verbose = FALSE)

# Visualizing clustering in PCs space # A TESTER
DimPlot(sc10x_subtype_SciGeneX, group.by = "seurat_clusters", label = TRUE)

# UMAP visualization of cells which changed clusters
sc10x_subtype_SciGeneX$cluster_change <- as.integer(as.character(sc10x_subtype_SciGeneX$seurat_clusters)) != 
                                         as.integer(as.character(sc10x_subtype_SciGeneX$seurat_original_clusters))



umap_cluster_change <- DimPlot(sc10x_subtype_SciGeneX, 
                        reduction = "umap", 
                        group.by = "cluster_change", 
                        cols = c("gray80", "red")) +
  ggtitle(paste0("Cells that changed clusters (", cellType, ")")) +
  theme_minimal()

print(umap_cluster_change )


# UMAP visualization

# renaming clusters starting from "1"
current_clusters <- as.numeric(as.character(Idents(sc10x_subtype_SciGeneX)))
new_clusters <- factor(current_clusters + 1)

# Assign back the new cluster identities
Idents(sc10x_subtype_SciGeneX) <- new_clusters
sc10x_subtype_SciGeneX$seurat_clusters <- new_clusters

sc10x_subtype_SciGeneX <- suppressWarnings(RunUMAP(sc10x_subtype_SciGeneX, dim = 1:10, verbose = FALSE))

umap_SciGeneX <-  DimPlot( sc10x_subtype_SciGeneX, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "seurat_clusters",
           label = TRUE) + 
    ggtitle( paste0("Map of cells with subclustering from SciGeneX ", cellType))
print(umap_SciGeneX)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, cellType, "_DimPlot_SciGeneX_subclusters.png")), plot = umap_SciGeneX, device ="png", width = 23, height = 13, units ="cm")

################################ Plotting cell groups showing to which original cluster they belonged to
umap_SGX_original_clusters <-  DimPlot( sc10x_subtype_SciGeneX, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "cluster_renames",
           label = TRUE) + 
    ggtitle( paste0("Map of cells with their original cluster ", cellType))
print(umap_SGX_original_clusters)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, cellType, "_DimPlot_cluster_of_origin_after_SciGeneX.png")), plot = umap_SGX_original_clusters, device ="png", width = 23, height = 13, units ="cm")

# Plot the presence of doublets (by DoubletFinder)

doublets_sgx <- DimPlot( sc10x_subtype_SciGeneX, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "multiplet_prediction",
           label = TRUE) + 
    ggtitle( "Map of doublet cells")
print(doublets_sgx)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, "DimPlot_doublets_after_SciGeneX.png")), plot = doublets_sgx, device ="png", width = 23, height = 13, units ="cm")
    
# Plot the cells percentage of mitochondrial genes in UMAP  
perc_mitok_sgx <- FeaturePlot( sc10x_subtype_SciGeneX, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "percent.mito") + 
  scale_colour_gradientn(colors = c(low = "#D0E1F9", high = "#053061")) +
  ggtitle( "Map of cells with level of percentage of\n mitochondrial genes")
print(perc_mitok_sgx)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, "FeaturePlot_perc_mitok_after_SciGeneX.png")), plot = perc_mitok_sgx, device ="png", width = 23, height = 13, units ="cm")

# Plot the cells RNA counts in UMAP
rna_count_sgx <- FeaturePlot( sc10x_subtype_SciGeneX, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nCount_RNA") + 
  scale_colour_gradientn(colors = c(low = "#D0E1F9", high = "#053061")) +
  ggtitle( "Map of cells with RNA counts")
 print(rna_count_sgx)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, "FeaturePlot_rna_count_after_SciGeneX.png")), plot = rna_count_sgx, device ="png", width = 23, height = 13, units ="cm")

# Plot the cells genes counts in UMAP
gene_count_sgx <- FeaturePlot(sc10x_subtype_SciGeneX, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nFeature_RNA") + 
  scale_colour_gradientn(colors = c(low = "#D0E1F9", high = "#053061")) +
  ggtitle( "Map of cells with \nnumber of detected genes")
print(gene_count_sgx)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, "FeaturePlot_gene_count_after_SciGeneX.png")), plot = gene_count_sgx, device ="png", width = 23, height = 13, units ="cm")


                                               ## Redoing SciGeneX analysis on the new seurat object to check the improvement##




# Select informative genes
res_subtype_SciGeneX <- select_genes(sc10x_subtype_SciGeneX,
                             distance = "pearson",
                             row_sum=5,
                             k = SCIGENEX_K_NEIGHBORHOOD)

# Cluster informative features
res_subtype_SciGeneX



## Construct and partition the graph

# testing parameter 
# reciprocal_neighborhood with inflation 1.8 -> module 9 cDC1 XCR1+ 
# closest_neighborhood with inflation 2 -> no good

res_subtype_SciGeneX <- gene_clustering(res_subtype_SciGeneX, method= SCIGENEX_METHOD, 
                                        inflation = SCIGENEX_INFLATION, 
                                        threads = 4,
                                        output_path = NULL,
                                        s = 5)

# Filter gene modules
stat_gene_modules_SciGeneX <- cluster_stats(res_subtype_SciGeneX)
stat_gene_modules_SciGeneX

# stat_gene_modules[order(stat_gene_modules$sd_total, decreasing = T),]

res_subtype_SciGeneX <- filter_cluster_size(res_subtype_SciGeneX, min_cluster_size = 4)
res_subtype_SciGeneX <- filter_cluster_sd(res_subtype_SciGeneX, min_sd = 0) #Set minimal sd


# Save gene list with one row per module, list of genes separated by a comma (for CellMarker2.0 handling)
library(rlist)


module_gene_list <- list()

for (module_nb in 1:length(res_subtype_SciGeneX@gene_clusters)) {
    gene_vector <- res_cluster_subtype@gene_clusters[[module_nb]]
    gene_list <- list(paste(gene_vector, collapse = ","))
    module_gene_list <- list.append(module_gene_list, gene_list)
    }
# Convert list elements to a single string where elements are separated by commas
    module_gene_string <- sapply(module_gene_list, function(x) paste(x, collapse = ","))

write.csv(module_gene_string, file= file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, "_ModuleGenes_list_after_SciGeneX.csv")) , na='')


# Creating a module score for each module for visualization

res_subtype_SciGeneX <- top_genes(res_subtype_SciGeneX) # creates a metadata column with 20 top genes for each module

# Display top gene of each module
for( cluster_id in 1:length(res_subtype_SciGeneX@gene_clusters) ){
  cat("Module SciGeneX:", as.numeric(cluster_id),"\n")
  cat(get_genes(res_subtype_SciGeneX , as.numeric(cluster_id) ,top = T),"\n \n")
}


# Compute the a list of gene in each module and the name of each module
cluster_feature_combine <- NULL
Module_name <- NULL
for( cluster_id in 1:length(res_subtype_SciGeneX@gene_clusters) ){
  cluster_feature <- list(get_genes(res_subtype_SciGeneX , cluster_id))
  cluster_feature_combine <- c(cluster_feature_combine,cluster_feature)
  Module_name <- c(Module_name,paste("Module",cluster_id, sep="_"))
}

# Create au module score of each module and compute score for each cell
sc10x_subtype_SciGeneX <- AddModuleScore(sc10x_subtype_SciGeneX, features = cluster_feature_combine ,name = "Module_")

# Plot of each module and cell in the module (UMAP)
p <-FeaturePlot_scCustom(sc10x_subtype_SciGeneX, features = Module_name, reduction  =  "umap", colors_use = viridis_light_high, pt.size = 0.05) & 
  theme(text = element_text(size = 10))

for(i in 1:length(Module_name)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

module_umap <- p+plot_layout(ncol = 7)
module_umap


#save plot 
#ggsave(plot = module_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_modules_umap_subclusters.png"))
      # , width = 23, height = 13, units ="cm")


# Plot of each module and cell in the module only for the 30 first modules (UMAP)

if (length(Module_name) >= 30) {
Module_name_sub <- paste0("Module_", 1:30)

p <-FeaturePlot_scCustom(sc10x_subtype_SciGeneX, features = Module_name_sub, reduction  =  "umap", colors_use = viridis_light_high, pt.size = 0.05) & 
  theme(text = element_text(size = 10))

for(i in 1:length(Module_name_sub)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

module_umap <- p+plot_layout(ncol = 7)
module_umap

#save plot 
#ggsave(plot = module_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_30_modules_umap_subclusters.png"))
     #  , width = 23, height = 13, units ="cm")
}

# Plot top gene in the first module (UMAP)

gene_in_module <- cluster_feature_combine[[1]] # it extracts genes from module 2
selection_gene_in_module <- gene_in_module[1:24] # it selects the first top 24 genes in this module

p_gene_module <-FeaturePlot_scCustom(sc10x_subtype_SciGeneX, features =selection_gene_in_module, reduction  =  "umap", colors_use = viridis_light_high)& 
  theme(text = element_text(size = 10))

for(i in 1:length(selection_gene_in_module)) {
  p_gene_module[[i]] <- p_gene_module[[i]] + NoLegend() + NoAxes()
}

p_gene_module

#ggsave(plot = p_gene_module, filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_top_genes_module_1_umap.png"))
       #, width = 23, height = 13, units ="cm")



##Display heatmap for a specific module and top genes
heatmap_one_module <- function(module_nb){

  p_heatmap <- plot_ggheatmap(res_subtype_SciGeneX[module_nb,], ident = Seurat::Idents(sc10x_subtype_SciGeneX), use_top_genes = T) + ggtitle("All clusters (top genes)") +  theme(strip.text.y = element_text(size=4))
  # save plot
  ggsave(plot = p_heatmap, filename = file.path(PATH_ANALYSIS_OUTPUT,  cellType, paste0("meaningful_modules_", outputFilesPrefix,SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_heatmap_after_SciGeneX.png"))
       , width = 23, height = 13, units ="cm")

  plot_heatmap(res_subtype_SciGeneX[module_nb,],use_top_genes = T, cell_clusters =Seurat::Idents(sc10x_subtype_SciGeneX)) #becarful heatmap interactive 
}
heatmap_one_module(meaningful_modules) # module numbers must be in order (same order as it will be displayed interactively on Rstudio

# save the SeuratObject with the new subclustering based on SciGeneX
saveRDS(sc10x_subtype_SciGeneX, file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "_seuratObject_after_SciGeneX.RDS")))


# save the ClusterSet object after re-running SciGeneX on "sc10x_subtype_SciGeneX"
saveRDS(res_subtype_SciGeneX, file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "_clusterset_object_after_SciGeneX.RDS")))

# Plotting each meaningful module using a violin plot
for (mod_nb in meaningful_modules) {
    mod_col_name <- paste0("Module_", mod_nb)
    # Plotting the module score across seurat_clusters
    vln_plot <- VlnPlot(
    object   = sc10x_subtype_SciGeneX,
    features = mod_col_name,
    group.by = "seurat_clusters",
    pt.size  = 0
)
    # Creating the file name
    filename <- paste0(outputFilesPrefix, cellType, "_", SCIGENEX_METHOD, 
                     "_inflation_", SCIGENEX_INFLATION, 
                     "_", mod_col_name, "_vln_plot.png")
    # Saving the plot
    ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, filename),
  plot = vln_plot,
  width = 20,
  height = 12,
  bg = "white"
)
}

#plotting violin plots for CD4 effector T cell markers
CD4_gene_list <- c("ICOS", "CXCR5", "PDCD1", "CXCR3", "CCR6", "IL21", "IL4", "IFNG", "IL10", "TBX21", "KLRB1", "TNF", "GATA3", "AHR", "CCR6", "RORA", "RORC", "CCL20", "DAB1", "FURIN", "TIGIT", "CCR4", "FOXP1", "FOXP3", "IL2RA", "CTLA4", "TGFB1", "ICOS", "CXCR5", "PDCD1", "IL21", "CD4", "CD8A")

for (gene in CD4_gene_list) {
    # Plotting the module score across seurat_clusters
    CD4_vln_plot <- VlnPlot(
    object   = sc10x_subtype_SciGeneX,
    features = gene,
    group.by = "seurat_clusters",
    pt.size  = 0
)
    # Creating the file name
    filename <- paste0(outputFilesPrefix, cellType, "_", SCIGENEX_METHOD, 
                     "_inflation_", SCIGENEX_INFLATION, 
                     "_", gene, "_vln_plot.png")
    # Saving the plot
    ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, filename),
  plot = CD4_vln_plot,
  width = 20,
  height = 12,
  bg = "white"
)
}




