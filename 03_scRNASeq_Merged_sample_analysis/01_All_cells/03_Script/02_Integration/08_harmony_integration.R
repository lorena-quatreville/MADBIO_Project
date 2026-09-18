# ########################################################
# The script integrates the data from the two datasets
# ########################################################

# ..........................................................................................................
## @knitr integration_harmony
# ..........................................................................................................

# saving seurat object without integration
sc10x_noHarmony <- sc10x

# creating a "batch" variable in the metadata of the seurat object
sc10x$batch <- factor(sc10x$orig.ident, levels = names(batch_info), labels = batch_info)

# Run Harmony
sc10x <- harmony::RunHarmony(object = sc10x, group.by.vars = "batch", reduction.use = "pca", dims.use = 1:nbPC_dimreduc, reduction.save = "harmony")

sc10x <- RunUMAP( sc10x, dims = 1:nbPC_dimreduc, reduction = "harmony", seed.use = SEED)

# Plot UMAP with clusters grouped by sample of origin
plot_sample_orig <- DimPlot(sc10x, group.by = "orig.ident", reduction = "umap", pt.size = 1.5) +
  ggtitle("UMAP with Cluster Names") +
  theme(legend.position = "bottom") +
  theme(legend.title = element_blank())
  print(plot_sample_orig)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "DimPlot_grouped_sample_orig_harmony.png")), plot = plot_sample_orig, width = 15, height = 10, dpi = 300, device ="png")


## Plot UMAP with cluster names in the legend and split by sample of origin
umap_split <- DimPlot(sc10x, group.by = "cell_type_simplified", split.by = "orig.ident", ncol= 3, label = FALSE) +
  ggtitle("UMAP with Cluster Names") +
  theme(legend.title = element_blank())
print(umap_split)
# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "umap_plot_split_orig_harmony.png")), plot = umap_split, width = 15, height = 10, dpi = 300)


## Same but highlighting possible T cell/B cell doublets

# Generate UMAP plots split by sample origin and highlighting gene expression of a T marker

gene_of_interest = "CD3D"

# Get unique sample identities
idents <- unique(sc10x$orig.ident)

# Generate individual plots for each sample of origin
umap_plots <- lapply(idents, function(ident) {
  FeaturePlot(sc10x, features = gene_of_interest, cells = rownames(sc10x@meta.data)[sc10x$orig.ident == ident]) + 
    ggtitle(paste(ident))
})

# Arrange the plots in a grid with 3 columns
umap_gene_split <- wrap_plots(umap_plots, ncol = 3)

# Display the arranged plots
print(umap_gene_split)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, gene_of_interest, "FeaturePlot_split_orig_harmony.png")), plot = umap_gene_split, width = 15, height = 10, dpi = 300)

# # Generate UMAP plots split by sample origin and highlighting gene expression of a B marker

gene_of_interest = "MS4A1"
# Get unique sample identities
idents <- unique(sc10x$orig.ident)

# Generate individual plots for each sample of origin
umap_plots <- lapply(idents, function(ident) {
  FeaturePlot(sc10x, features = gene_of_interest, cells = rownames(sc10x@meta.data)[sc10x$orig.ident == ident]) + 
    ggtitle(paste(ident))
})

# Arrange the plots in a grid with 3 columns
umap_gene_split <- wrap_plots(umap_plots, ncol = 3)

# Display the arranged plots
print(umap_gene_split)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, gene_of_interest, "FeaturePlot_split_orig_harmony.png")), plot = umap_gene_split, width = 15, height = 10, dpi = 300)




## Plot using harmony axes with individual names in the legend and group by sample of origin

# Plot the cells sample of origin in UMAP
umap_group_orig <- DimPlot( sc10x, reduction = "harmony", group.by = "orig.ident") + 
  ggtitle( "Map of cells by sample of origin")
print(umap_group_orig)
# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "harmony_plot_group_orig.png")), plot = umap_group_orig, width = 12, height =8, dpi = 300)


## Plot a UMAP grouped by condition (chronic vs symptomatic)

sc10x$condition <- factor(sc10x$orig.ident, levels = names(condition_info), labels = condition_info)

# Plot UMAP grouped by condition
umap_condition <- DimPlot(sc10x, reduction = "umap", group.by = "condition") + 
  ggtitle("UMAP Grouped by Condition") + 
  theme(legend.title = element_blank())

# Display the plot
print(umap_condition)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "umap_plot_group_condition_harmony.png")), plot = umap_condition, width = 12, height =8, dpi = 300)





# IDENTIFY CLUSTERS
###################

# ..........................................................................................................
## @knitr integration_identifyClusters
# ..........................................................................................................


Idents(sc10x)<-"seurat_clusters"

# Identify clusters of cells by graph approach
nbPC_findclusters=FINDCLUSTERS_USE_PCA_NBDIMS
if(FINDCLUSTERS_USE_PCA_NBDIMS>nbPC){
  warning( paste0( "Number of computed PCs  (", nbPC, ") smaller than requested PCs for 'findclusters' (", FINDCLUSTERS_USE_PCA_NBDIMS,"), setting lower PC number (", nbPC, ")..." ))
  nbPC_findclusters = nbPC
}  
sc10x <- FindNeighbors(object    = sc10x,
                       k.param   = FINDNEIGHBORS_K, 
                       reduction = "harmony",
                       dims      = 1:nbPC_findclusters,
                       verbose   = .VERBOSE);

sc10x <- FindClusters(object             = sc10x,
                      resolution         = FINDCLUSTERS_RESOLUTION,
                      algorithm          = FINDCLUSTERS_ALGORITHM,
                      method = "igraph",
                      temp.file.location = "/tmp/",
                      verbose            = .VERBOSE);


# Rename the clusters if required
if( exists( "RENAME_CLUSTERS") && SAMPLE_NAME %in% names( RENAME_CLUSTERS)){
  Idents( sc10x) = "seurat_clusters"
  sc10x = AddMetaData( sc10x, metadata = Idents( sc10x), col.name = "original_seurat_clusters")
  new_names = Idents( sc10x)
  for( cluster_id in levels( Idents( sc10x))){
    new_names[ which( Idents( sc10x) == cluster_id)] = RENAME_CLUSTERS[[ SAMPLE_NAME]][[ cluster_id]]
  }
  sc10x = AddMetaData( sc10x, metadata = new_names, col.name = "seurat_clusters")
  
}

Idents( sc10x) = "seurat_clusters"



# Plot the map with clusters found from harmony with ggplot

cluster_plot_harmony <- DimPlot( sc10x, 
         reduction = ifelse(exists("useReduction"), useReduction, "umap"), 
         group.by = "seurat_clusters",
         label = TRUE) + 
  ggtitle( "Map of cells with clusters")
print(cluster_plot_harmony)
# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "umap_plot_cluster_harmony.png")), plot = cluster_plot_harmony, width = 12, height =8, dpi = 300)



## Get the percentage of cells belonging to each individual in each cluster

# Extract metadata
metadata <- sc10x@meta.data

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
write.table(cluster_composition, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cluster_compo_indiv_harmony.csv")), sep = "\t", row.names = FALSE, quote = FALSE)

## Same but with the cell percentage belonging to each condition (chronic vs symptomatic)

# Extract metadata
metadata <- sc10x@meta.data

# Count cells per cluster and condition
cluster_condition_composition <- metadata %>%
  group_by(seurat_clusters, condition) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  ungroup()

# Calculate percentage of cells in each condition per cluster
cluster_condition_composition <- cluster_condition_composition %>%
  group_by(seurat_clusters) %>%
  mutate(percent = (cell_count / sum(cell_count)) * 100) %>%
  select(seurat_clusters, condition, percent) %>%
  pivot_wider(names_from = condition, values_from = percent, values_fill = 0)  # Wide format for better readability

# Write to a TSV file
write.table(cluster_condition_composition, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cluster_compo_condition_harmony.csv")), sep = "\t", row.names = FALSE, quote = FALSE)






# MARKER GENES
##############

# ..........................................................................................................
## @knitr integration_markerGenes
# ..........................................................................................................

# Identify marker genes
markers = FindAllMarkers( object          = sc10x,
                          test.use        = FINDMARKERS_METHOD,
                          only.pos        = FINDMARKERS_ONLYPOS,
                          min.pct         = FINDMARKERS_MINPCT,
                          logfc.threshold = FINDMARKERS_LOGFC_THR,
                          return.thresh   = FINDMARKERS_PVAL_THR,
                          random.seed     = SEED,
                          verbose         = .VERBOSE);

# Save markers list as 'tsv' table
write.table( markers,
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "integration_MarkerGenes.tsv")),
             quote = FALSE,
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers_table = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP_TABLE)) x else head( x, n = FINDMARKERS_SHOWTOP_TABLE));
});

# Filter markers by cluster (TODO: check if downstream code works when no markers found)
topMarkers_heatmap = by( markers, markers[["cluster"]], function(x)
{
  # Filter markers based on adjusted PValue
  x = x[ x[["p_val_adj"]] < FINDMARKERS_PVAL_THR, , drop = FALSE];
  # Sort by decreasing logFC
  x = x[ order(abs(x[["avg_log2FC"]]), decreasing = TRUE), , drop = FALSE ]
  # Return top ones
  return( if(is.null( FINDMARKERS_SHOWTOP_HEATMAP)) x else head( x, n = FINDMARKERS_SHOWTOP_HEATMAP));
});

# Merge marker genes in a single data.frame and render it as datatable
topMarkers_tableDF = do.call( rbind, topMarkers_table);
topMarkers_heatmapDF = do.call( rbind, topMarkers_heatmap);

# Select and order columns to be shown in datatable
topMarkers_tableDT = topMarkers_tableDF[c("gene", "cluster", "avg_log2FC", "p_val_adj")]



# ..........................................................................................................
## @knitr integration_markerGenes_table
# ..........................................................................................................

# Create datatable
datatable( topMarkers_tableDT,
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Cluster", "Avg. LogFC", "Adj. Pvalue"),
           caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP_TABLE), "All", paste("Top", FINDMARKERS_SHOWTOP_TABLE)), "marker genes for each cluster"),
           extensions = c('Buttons', 'Select'),
           options = list(dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                          autoWidth = FALSE,
                          buttons = exportButtonsListDT,
                          columnDefs = list(
                            list( # Center all columns except first one
                              targets = 1:(ncol( topMarkers_tableDT)-1),
                              className = 'dt-center'),
                            list( # Set renderer function for 'float' type columns (LogFC)
                              targets = ncol( topMarkers_tableDT)-2,
                              render = htmlwidgets::JS("function ( data, type, row ) {return type === 'export' ? data : data.toFixed(4);}")),
                            list( # Set renderer function for 'scientific' type columns (PValue)
                              targets = ncol( topMarkers_tableDT)-1,
                              render = htmlwidgets::JS( "function ( data, type, row ) {return type === 'export' ? data : data.toExponential(4);}"))),
                          #fixedColumns = TRUE, # Does not work well with filter on this column
                          #fixedHeader = TRUE, # Does not work well with 'scrollX'
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
  # Add bar relative to logFC
  formatStyle( columns = "avg_log2FC",
               background = styleColorBar( data = range( topMarkers_tableDT[["avg_log2FC"]]), 'lightblue', angle = -90),
               backgroundSize = '95% 50%',      # Set horizontal and vertical span in cell
               backgroundRepeat = 'no-repeat',
               backgroundPosition = 'center') %>%
  # Add color from cluster
  formatStyle( columns = "cluster",
               backgroundColor = styleEqual( names(clustersColor),
                                             scales::alpha(clustersColor, 0.3)));



# ..........................................................................................................
## @knitr integration_markerGenes_expression_projection
# ..........................................................................................................
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of marker genes on dimreduc figures for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers_heatmap), function(clusterName)
{
  cat("##### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");
  
  # Highlight cells of current cluster on a dimreduc plot
  highlightClusterPlot(clusterName, seuratObject = sc10x, reduction = ifelse( exists("useReduction"), useReduction, "umap"));
  
  # Plots expression on projected cells
  invisible( lapply( topMarkers_heatmap[[clusterName]][["gene"]], function(featureName)
  {
    print( FeaturePlot( sc10x, features = featureName, reduction = ifelse( exists("useReduction"), useReduction, "umap"), order = TRUE) +
             theme( axis.title.x = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = "none"));
  }));
  
  cat(" \n \n"); # Required for '.tabset'
}));



# ..........................................................................................................
## @knitr integration_markerGenes_expression_violin
# ..........................................................................................................

# Plot expression values of marker genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( topMarkers_heatmap), function(clusterName)
{
  cat("##### Cl. <span style='border-radius: 3px; border: 3px solid ", clustersColor[clusterName], "; padding:0px 2px'>", clusterName, "</span>\n");
  
  # Remind cluster name in an empty figure to keep consistent alignment of panels between tabs
  plot( c( 0, 1), c( 0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n');
  text( x = 0.5, y = 0.5, paste( "Cluster", clusterName), cex = 2, col = clustersColor[clusterName]);
  
  # Violinplot for expression value of marker genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( topMarkers_heatmap[[clusterName]][["gene"]], violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor));
  
  cat(" \n \n"); # Required for '.tabset'
}))


# Save markers list with one row per cluster, marker genes separated by a comma (for CellMarker2.0 handling)


markerGenes_clusters_list <- list()

for (df_nb in seq_along(topMarkers_table)) {
    df_oneCluster <- topMarkers_table[[df_nb]]
    markerGenes_col <- df_oneCluster$gene
    markerGenes_list <- strsplit(markerGenes_col, ",")
    markerGenes_clusters_list <- list.append(markerGenes_clusters_list, markerGenes_list)
    # Convert list elements to a single string where elements are separated by commas
    markerGenes_string <- sapply(markerGenes_clusters_list, function(x) paste(x, collapse = ","))
}
    
    
write.csv(markerGenes_string, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "MarkerGenes_list_harmony.csv")) , na='')






# ..........................................................................................................
## @knitr clean_RDATA
# ..........................................................................................................




















