
###############################################################################
#                           Subclustering                            #
###############################################################################

## If needed to load the merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x) <- "cluster_renames"

# Subset to only include cells from Individual A (Cluster 1) and Individual B (Cluster 2)
sc10x_subset <- subset(
  sc10x,
  subset = multiplet_type == "Non-classical monocyte" &
           orig.ident %in% c(sample_id_A, sample_id_B)
)


# Set cluster identity for this subset
Idents(sc10x_subset) <- sc10x_subset$orig.ident

sc10x_subset <- FindNeighbors(object    = sc10x_subset,
                                k.param   = FINDNEIGHBORS_K, 
                                reduction = "pca",
                                dims      = 1:nbPC_findclusters,
                                verbose   = .VERBOSE);

sc10x_subset <- FindClusters(object       = sc10x_subset,
                               resolution         = 0.55,
                               algorithm          = FINDCLUSTERS_ALGORITHM,
                               temp.file.location = "/tmp/",
                               verbose            = .VERBOSE);

#Parameter of SciGeneX
SCIGENEX_K_NEIGHBORHOOD = 80
SCIGENEX_METHOD         = "closest_neighborhood"  # 2 methodes : reciprocal_neighborhood or closest_neighborhood
SCIGENEX_INFLATION      = 1.8                        # number of module


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



################################ Plotting cell groups isolated by the subclustering on UMAP

umap_subclusters <-  DimPlot( sc10x_subset, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "seurat_clusters",
           label = TRUE) + 
    ggtitle( paste0("Map of cells with subclustering "))
print(umap_subclusters)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "_DimPlot_sub_clusters.png")), plot = umap_subclusters, device ="png", width = 23, height = 13, units ="cm")


################################ Identifying doublet clusters

# Plot the presence of doublets (by DoubletFinder)

doublets <- DimPlot( sc10x_subset, 
           reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
           group.by = "multiplet_prediction",
           label = TRUE) + 
    ggtitle( "Map of doublet cells")
print(doublets)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "DimPlot_doublets.png")), plot = doublets, device ="png", width = 23, height = 13, units ="cm")
    
# Plot the cells percentage of mitochondrial genes in UMAP  
perc_mitok <- FeaturePlot( sc10x_subset, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "percent.mito") + 
  scale_colour_gradientn(colors = c(low = "#D0E1F9", high = "#053061")) +
  ggtitle( "Map of cells with level of percentage of\n mitochondrial genes")
print(perc_mitok)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "FeaturePlot_perc_mitok.png")), plot = perc_mitok, device ="png", width = 23, height = 13, units ="cm")

# Plot the cells RNA counts in UMAP
rna_count <- FeaturePlot( sc10x_subset, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nCount_RNA") + 
  scale_colour_gradientn(colors = c(low = "#D0E1F9", high = "#053061")) +
  ggtitle( "Map of cells with RNA counts")
 print(rna_count)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "FeaturePlot_rna_count.png")), plot = rna_count, device ="png", width = 23, height = 13, units ="cm")

# Plot the cells genes counts in UMAP
gene_count <- FeaturePlot( sc10x_subset, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "nFeature_RNA") + 
  scale_colour_gradientn(colors = c(low = "#D0E1F9", high = "#053061")) +
  ggtitle( "Map of cells with \nnumber of detected genes")
print(gene_count)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "FeaturePlot_gene_count.png")), plot = gene_count, device ="png", width = 23, height = 13, units ="cm")


# Plot the UMAP grouped by sample of origin

# Plot UMAP with clusters grouped by sample of origin
plot_sample_orig <- DimPlot(sc10x_subset, group.by = "orig.ident") +
  ggtitle("UMAP with Cluster Names") +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank())
  print(plot_sample_orig)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "DimPlot_grouped_sample_orig.png")), plot = plot_sample_orig, device ="png",  width = 23, height = 13, units ="cm")





############################### GO over-representation test

module_nb = 1

dir.create(file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb)))


# Creating a gene_list extracted from the module of your choice
gene_list <- res_cluster_subtype@gene_clusters[[module_nb]]


# Mapping genes from gene_list to gene symbols from the human database
mapped_gene_list <- mapIds(org.Hs.eg.db, keys = gene_list, keytype = "SYMBOL", column = "ENTREZID", multiVals = "first")

# Remove NAs (unmapped genes)
mapped_gene_list <- na.omit(mapped_gene_list)

# BP
go_results_BP <- enrichGO(
  gene          = mapped_gene_list,  # List of differentially expressed genes (Entrez IDs)
  OrgDb         = org.Hs.eg.db,  # Organism annotation database
  keyType       = "ENTREZID",  # The type of input gene identifiers
  ont           = "BP",  # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
  pAdjustMethod = "BH",  # Multiple testing correction method
  pvalueCutoff  = 0.05,  # Adjusted p-value threshold
  qvalueCutoff  = 0.05,  # FDR threshold
  readable      = TRUE  # Converts Entrez IDs back to gene symbols
)

# MF
go_results_MF <- enrichGO(
  gene          = mapped_gene_list,  # List of differentially expressed genes (Entrez IDs)
  OrgDb         = org.Hs.eg.db,  # Organism annotation database
  keyType       = "ENTREZID",  # The type of input gene identifiers
  ont           = "MF",  # Ontology: "BP" (Biological Process), "MF" (Molecular Function), or "CC" (Cellular Component)
  pAdjustMethod = "BH",  # Multiple testing correction method
  pvalueCutoff  = 0.05,  # Adjusted p-value threshold
  qvalueCutoff  = 0.05,  # FDR threshold
  readable      = TRUE  # Converts Entrez IDs back to gene symbols
)

# View results
head(go_results_BP)
head(go_results_MF)

write.csv(as.data.frame(go_results_BP), file = file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "GO_analysis_results_deg_BP.csv")), row.names = FALSE)
write.csv(as.data.frame(go_results_MF), file = file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "GO_analysis_results_deg_MF.csv")), row.names = FALSE)


            

