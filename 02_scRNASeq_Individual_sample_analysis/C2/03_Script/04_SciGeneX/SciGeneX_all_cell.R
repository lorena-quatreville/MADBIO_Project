library(Seurat)
library(scigenex)
library(dplyr)
library(dbplyr)
library(patchwork)
library(ggplot2)
library(viridis)
library(scCustomize)
library(DT)
library( scales) 

#install.packages("dbplyr")
#devtools::install_github("RobinHankin/Brobdingnag")
set_verbosity(1)

cellType = "all_cell"

###############################################################################
#                                 All data analysis                           #
###############################################################################


#Parameter of SciGeneX
SCIGENEX_K_NEIGHBORHOOD = 80
SCIGENEX_METHOD         = "reciprocal_neighborhood"  # 2 methodes : reciprocal_neighborhood or closest_neighborhood
SCIGENEX_INFLATION      = 2                          # number of module

# Select informative genes
res <- select_genes(sc10x,
                    distance   = "pearson",
                    row_sum    = 5,
                    k          = SCIGENEX_K_NEIGHBORHOOD)

# Cluster informative features
res
dir.create(file.path(PATH_ANALYSIS_OUTPUT))
dir.create(file.path(PATH_ANALYSIS_OUTPUT, cellType))
saveRDS(res, file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, "clusterset_object.RDS")))
res <- readRDS(file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, "clusterset_object.RDS")))

## Construct and partition the graph
res_cluster <- gene_clustering(res, method= SCIGENEX_METHOD, 
                               inflation = SCIGENEX_INFLATION, 
                               threads = 4,
                               output_path = NULL,
                               s = 5)

# Filter gene modules
stat_gene_modules <- cluster_stats(res_cluster)
stat_gene_modules

# stat_gene_modules[order(stat_gene_modules$sd_total, decreasing = T),]

res_cluster <- filter_cluster_size(res_cluster, min_cluster_size = 4)
res_cluster <- filter_cluster_sd(res_cluster, min_sd = 0) #Set minimal sd

# Display the heatmap of gene clusters
res_cluster <- top_genes(res_cluster)

# Display top gene of each module
for( cluster_id in 1:length(res_cluster@gene_clusters) ){
  cat("Module SciGeneX:", as.numeric(cluster_id),"\n")
  cat(get_genes(res_cluster , as.numeric(cluster_id) ,top = T),"\n \n")
}


# Create a dataframe with 2 columns, one for gene names (that is not rownames) and one for gene clusters
Gene_in_cluster <- as.data.frame(gene_cluster(res_cluster))
Gene_in_cluster <- cbind(as.data.frame(rownames(Gene_in_cluster)),as.data.frame(gene_cluster(res_cluster)))

# Define a set of colors for clusters (based on ggplot default)
Idents( sc10x) = "seurat_clusters"
clustersColor = hue_pal()( nlevels( Idents( sc10x)));
names( clustersColor) = levels( Idents( sc10x));


# intercative table of gene in each module
datatable( Gene_in_cluster,
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Module"),
           caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP_TABLE), "All", paste("Top", FINDMARKERS_SHOWTOP_TABLE)), "marker genes for each cluster"),
           extensions = c('Buttons', 'Select'),
           options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT, # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT,
                          columnDefs = list(
                            list( # Center all columns except first one
                              targets = 1:(ncol( Gene_in_cluster)-1),
                              className = 'dt-center')),
                          fixedColumns = TRUE, # Does not work well with filter on this column
                          fixedHeader = TRUE, # Does not work well with 'scrollX'
                          lengthMenu = list(c( 10, 50, 100, -1),
                                            c( 10, 50, 100, "All")),
                          orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                          processing = TRUE,
                          #search.regex= TRUE, # Does not work well with 'search.smart'
                          search.smart = TRUE,
                          select = TRUE, # Enable ability to select rows
                          scrollCollapse = TRUE,
                          scroller = TRUE,  # Only load visible data
                          scrollX = TRUE,
                          scrollY = "525px",
                          stateSave = TRUE)) %>%
  # Add color from cluster
  formatStyle( columns = "gene_cluster(res_cluster)",
               backgroundColor = styleEqual( names(clustersColor),
                                             scales::alpha(clustersColor, 0.3)))



# Compute the a list of gene in each module and the name of each module
cluster_feature_combine <- NULL
Module_name <- NULL
for( cluster_id in 1:length(res_cluster@gene_clusters) ){
  cluster_feature <- list(get_genes(res_cluster , cluster_id))
  cluster_feature_combine <- c(cluster_feature_combine,cluster_feature)
  Module_name <- c(Module_name,paste("Module",cluster_id, sep="_"))
}

# Create a module score of each module and compute score for each cell
sc10x <- AddModuleScore(sc10x, features = cluster_feature_combine ,name = "Module_")

# Plot of each module and cell in the module 
p <-FeaturePlot_scCustom(sc10x, features = Module_name, reduction  =  "umap", colors_use = viridis_light_high,pt.size = 0.05) & 
  theme(text = element_text(size = 10))

for(i in 1:length(Module_name)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}

module_umap <- p+plot_layout(ncol = 7)
  
module_umap


# save plot
ggsave(plot = module_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(cellType, "_", outputFilesPrefix,SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_module.png"))
       , width = 23, height = 13, units ="cm")


# Heatmap of module gene

p_heatmap <- plot_ggheatmap(res_cluster, ident = Seurat::Idents(sc10x),use_top_genes = T)
# save plot
ggsave(plot = p_heatmap, filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(cellType, "_", outputFilesPrefix,SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_heatmap.png"))
       , width = 23, height = 13, units ="cm")

plot_heatmap(res_cluster[4:5,],use_top_genes = T) #becarful heatmap interactive 


##### Display heatmap of modules and top genes per cell cluster

Idents(sc10x) <- "seurat_clusters"
p_heatmap <- plot_ggheatmap(res_cluster, ident = Seurat::Idents(sc10x),use_top_genes = T)


# save plot
ggsave(plot = p_heatmap, filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(cellType, "_", outputFilesPrefix,SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_heatmap.png"))
       , width = 23, height = 13, units ="cm")

plot_heatmap(res_cluster[4:5,],use_top_genes = T) #becarful heatmap interactive 

################################## Save gene list with one row per module, list of genes separated by a comma (for CellMarker2.0 handling)
library(rlist)


module_gene_list <- list()

for (module_nb in 1:length(res_cluster@gene_clusters)) {
    gene_vector <- res_cluster@gene_clusters[[module_nb]]
    gene_list <- list(paste(gene_vector, collapse = ","))
    module_gene_list <- list.append(module_gene_list, gene_list)
   }
# Convert list elements to a single string where elements are separated by commas
    module_gene_string <- sapply(module_gene_list, function(x) paste(x, collapse = ","))

write.csv(module_gene_string, file= file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, "_ModuleGenes_list.csv")) , na='')




