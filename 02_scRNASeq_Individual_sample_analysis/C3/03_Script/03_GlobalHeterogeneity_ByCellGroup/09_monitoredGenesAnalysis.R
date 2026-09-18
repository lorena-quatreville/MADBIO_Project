# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# MONITORED GENES
#################

## @knitr heterogeneity_monitoredGenes

# Just remind the warning for genes names not in object
if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( monitoredGenesNotFound, collapse=" - ")));




## @knitr heterogeneity_monitoredGenes_heatmap

# DoHeatmap replaced by use of iHeatmapr after testing several options
#htmltools::tagList( ggplotly( DoHeatmap(object = sc10x, features = unlist(MONITORED_GENES))));

# Get the matrix of expression and associated clusters from Seurat object
expMat = as.matrix( GetAssayData( sc10x));
clusterID = Idents( sc10x);

# Select monitored genes and reorder cells to group clusters together
monitoredGenes = unlist( MONITORED_GENES);
clusterOrdering = order( clusterID);

expMat = expMat[monitoredGenes, clusterOrdering];
clusterID = clusterID[ clusterOrdering];

# Prepare rows and columns annotation bars (monitored group and cluster respectively)
rowsAnnot = data.frame( Monitored = fct_rev(fct_inorder(rep( names( MONITORED_GENES), sapply( MONITORED_GENES, length)))));
colsAnnot = data.frame( Cluster = clusterID);

# Prepare unique rows and cols names for pheatmap (annotation rows) and match with rowAnnots and colAnnots row names
originalRowNames = rownames( expMat);
originalColNames = colnames( expMat);
rownames( expMat) = make.unique( originalRowNames);
colnames( expMat) = make.unique( originalColNames);
rownames( rowsAnnot) = rownames( expMat);
rownames( colsAnnot) = colnames( expMat);

# Prepare colors of monitored genes groups for pheatmap (requires named vector matching data factor levels)
monitoredColors = rainbow( nlevels( rowsAnnot[["Monitored"]]), s = 0.8);
names( monitoredColors) = levels( rowsAnnot[["Monitored"]]);

pheatmap( expMat,
          color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
          cluster_rows = FALSE,
          cluster_cols = FALSE,
          annotation_row = rowsAnnot,
          annotation_col = colsAnnot,
          labels_row = originalRowNames,
          annotation_colors = list( Monitored = monitoredColors,
                                    Cluster = clustersColor),
          show_colnames = FALSE);


## @knitr heterogeneity_monitoredGenes_expression_projection
# Can define useReduction='tsne' or useReduction='umap' before (defaults to 'umap' otherwise)

# Plot expression values of monitored genes (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("#####", monitoredGroup, "\n");

  # Plots expression on projected cells (or error message if feature not found)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], function(featureName)
  {
    print(
      tryCatch( suppressWarnings( # Raise a warning # A warning is raised when all values are 0 (e.g. contamination genes after filtering)
                FeaturePlot( sc10x, features = featureName, reduction = ifelse(exists("useReduction"), useReduction, "umap"), order = TRUE)  +
                  theme( axis.title.x = element_blank(),
                         axis.title.y = element_blank(),
                         legend.position = "none")),
                error = function(e){return(conditionMessage(e))})); # If gene not found in object (should have been removed earlier anyway)
  }));

  cat(" \n \n"); # Required for '.tabset'
}));




## @knitr heterogeneity_monitoredGenes_expression_violin

# Plot expression values of monitored genes as violinplot for each cluster (TODO: message if list empty)
invisible( lapply( names( MONITORED_GENES), function(monitoredGroup)
{
  cat("#####", monitoredGroup, "\n");

  # Violinplot for expression value of monitored genes by cluster (+ number of 'zero' and 'not zero' cells)
  invisible( lapply( MONITORED_GENES[[monitoredGroup]], violinFeatureByCluster, seuratObject = sc10x, clustersColor = clustersColor));

  cat(" \n \n"); # Required for '.tabset'
}));


