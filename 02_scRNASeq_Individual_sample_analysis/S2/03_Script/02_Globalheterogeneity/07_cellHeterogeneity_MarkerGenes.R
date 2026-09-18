# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################

# MARKER GENES
##############

## @knitr heterogeneity_markerGenes

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
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "MarkerGenes.tsv")),
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




## @knitr heterogeneity_markerGenes_table

# Create datatable
datatable( topMarkers_tableDT,
           class = "compact",
           filter="top",
           rownames = FALSE,
           colnames = c("Gene", "Cluster", "Avg. LogFC", "Adj. Pvalue"),
           caption = paste(ifelse( is.null( FINDMARKERS_SHOWTOP_TABLE), "All", paste("Top", FINDMARKERS_SHOWTOP_TABLE)), "marker genes for each cluster"),
           extensions = c('Buttons', 'Select'),
           options = list(selection = 'none', dom = "<'row'<'col-sm-8'B><'col-sm-4'f>> <'row'<'col-sm-12'l>> <'row'<'col-sm-12'rt>> <'row'<'col-sm-12'ip>>", # Set elements for CSS formatting ('<Blf><rt><ip>')
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
                          select = TRUE, # Enable ability to select rows #################################################### à enlever car conflit?
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



## @knitr heterogeneity_markerGenes_heatmap_mean

# Get the matrix of expression and associated clusters from Seurat object
# ........................................................................
expMat = as.matrix( GetAssayData( sc10x));
Idents( sc10x) = "seurat_clusters"
clusterID = Idents( sc10x);

# Select marker genes and reorder cells to group clusters together
topMarkersGenes = topMarkers_heatmapDF[["gene"]];
clusterOrdering = order( clusterID);

expMat = expMat[topMarkersGenes, clusterOrdering];
clusterID = clusterID[clusterOrdering];

# Compute the mean of top markers in clusters and produce matrix with the result
# ........................................................................
all_mean_expression = vector()
cluster_set = levels( clusterID)
clusters_marker_genes = topMarkers_heatmapDF[ which( topMarkers_heatmapDF$cluster %in% cluster_set), "gene"]
for( cluster_id in cluster_set){
  clusters_cells = names( clusterID)[ which( clusterID == cluster_id)]
  all_mean_expression = c( all_mean_expression, BiocGenerics::rowMeans( expMat[ clusters_marker_genes, clusters_cells]))
}
meanExpMat = t( matrix( all_mean_expression, byrow = TRUE, ncol = FINDMARKERS_SHOWTOP_HEATMAP*length( cluster_set)))
colnames( meanExpMat) = cluster_set
rownames( meanExpMat) = paste0( topMarkers_heatmapDF[ which( topMarkers_heatmapDF$cluster %in% cluster_set), "cluster"], '.', clusters_marker_genes)

# Plot a heatmap of the matrix of mean expression of top marker genes in clusters
# ............................................................................

cat("\n \n")
pheatmap( expMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(80),
          breaks = seq( -4.2, 4.2, 0.1),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "row",
          annotation_row = data.frame(Markers = topMarkers_heatmapDF[ rownames( meanExpMat), "cluster"], stringsAsFactors = FALSE, row.names = rownames( meanExpMat) ),
          annotation_col = data.frame(Cluster = clusterID, stringsAsFactors = FALSE, row.names = colnames( expMat)),
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = FALSE,
          fontsize_row = 5,
          main = "Scaled normalized expression by cell\nof top markers");

cat("\n \n")
DotPlot( sc10x, 
         features = unique( clusters_marker_genes),
         cols = "RdBu") +
  theme( axis.text.x = element_text( angle = 45, hjust = 1)) +
  ggtitle( "Scaled normalized mean expression of topmarker genes in clusters")

cat("\n \n")
pheatmap( meanExpMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          scale = "row",
          annotation_row = data.frame(Markers = topMarkers_heatmapDF[ rownames( meanExpMat), "cluster"], stringsAsFactors = FALSE, row.names = rownames( meanExpMat) ),
          annotation_col = data.frame(Cluster = cluster_set, stringsAsFactors = FALSE, row.names = cluster_set),
          annotation_colors = list( Markers = clustersColor,
                                    Cluster = clustersColor),
          show_colnames = TRUE,
          fontsize_row = 5,
          main = "Scaled mean normalized expression\nby cluster of top markers");


## @knitr heterogeneity_markerGenes_expression_projection
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




## @knitr heterogeneity_markerGenes_expression_violin

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
library(rlist)

markerGenes_clusters_list <- list()

for (df_nb in seq_along(topMarkers_table)) {
    df_oneCluster <- topMarkers_table[[df_nb]]
    markerGenes_col <- df_oneCluster$gene
    markerGenes_list <- strsplit(markerGenes_col, ",")
    markerGenes_clusters_list <- list.append(markerGenes_clusters_list, markerGenes_list)
    # Convert list elements to a single string where elements are separated by commas
    markerGenes_string <- sapply(markerGenes_clusters_list, function(x) paste(x, collapse = ","))
}
    
    
write.csv(markerGenes_string, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "MarkerGenes_list.csv")) , na='')

