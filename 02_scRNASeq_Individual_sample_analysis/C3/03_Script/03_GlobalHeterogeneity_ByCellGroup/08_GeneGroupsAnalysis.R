# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# GROUP GENES
#################

## @knitr heterogeneity_chosenGroupGenes_heatmap_mean

# Compute the mean of chosen group of genes in clusters and produce matrix with the result
# ........................................................................

for( groupid in names( GROUP_GENES)){ 
  
  cat("\n \n")
  cat("#### ", groupid, "\n");
  
  # Get the genes of the current group
  group_genes = GROUP_GENES[[ groupid]]
  
  # Get the matrix of counts
  expMat = as.matrix( GetAssayData( sc10x));
  group_genes = intersect( group_genes, rownames( expMat))
  
  # Get the cluster ID of the cells
  Idents( sc10x) = "seurat_clusters"
  clusterID = Idents( sc10x);
  
  # Reorder cells to group clusters together
  clusterOrdering = order( clusterID);
  clusterID = clusterID[clusterOrdering];
  
  # Get the count matrix with the right genes and the right cell order
  expMat = expMat[ group_genes, clusterOrdering];
  
  # Compute the mean of cells by gene and cluster and build a matrix with the results
  all_mean_expression = vector()
  cluster_set = levels( clusterID)
  for( cluster_id in cluster_set){
    clusters_cells = names( clusterID)[ which( clusterID == cluster_id)]
    all_mean_expression = c( all_mean_expression, BiocGenerics::rowMeans( expMat[ group_genes, clusters_cells]))
  }
  meanExpMat = t( matrix( all_mean_expression, byrow = TRUE, ncol = length( group_genes)))
  colnames( meanExpMat) = cluster_set
  rownames( meanExpMat) = group_genes
  
  # Plot the heatmap of mean counts
  # ................................
  
  DotPlot( sc10x, 
           features = sort( group_genes),
           cols = "RdBu") +
    theme( axis.text.x = element_text( angle = 45, hjust = 1)) +
    ggtitle( paste( "Scaled mean expression of genes in clusters for\ngroup of genes:", groupid))
  
  pheatmap( meanExpMat,
            color = colorRampPalette( rev( RColorBrewer::brewer.pal(n = 7, name = "RdBu")))(100),
            cluster_rows = TRUE,
            cluster_cols = FALSE,
            scale = "row",
            annotation_col = data.frame(Cluster = cluster_set, stringsAsFactors = FALSE, row.names = cluster_set),
            annotation_colors = list( Cluster = clustersColor),
            show_colnames = TRUE,
            fontsize_row = 10,
            border_color = NA, 
            main = paste( "Scaled mean expression of genes in clusters for\ngroup of genes:", groupid));

  
  # Compute the Module Score of group
  # ..................................
  
  sc10x = AddModuleScore( sc10x, features = list( group_genes), name = groupid)
  module_name = paste0( groupid, "1")
  
  
  # Plot the modules scores by cluster and condition
  module_score_df = sc10x@meta.data[ , c( "seurat_clusters", module_name)]
  print(
    ggplot( module_score_df,
            aes_string( x="seurat_clusters", y=module_name, fill = "seurat_clusters", color = "seurat_clusters")) +
      geom_violin( position = "dodge") +
      stat_summary( fun = mean, geom = "point", color = "black", position = position_dodge( width = 0.9)) +
      stat_summary(fun.data = mean_se, geom = "errorbar", color = "black", size = 0.5, width = 0.3, position = position_dodge( width = 0.9)) +
      theme_minimal() +
      theme( axis.text.x = element_text( size = "12"), axis.text.y = element_text( size = "12")) +
      ggtitle( paste( "Module scores of", groupid, "by cluster", "\nwith mean and standard deviation"))
  )
  
  # # Compute the comparison statistics of module score by cluster between conditions with wilcoxon test
  # module_stats_df = data.frame()
  # for( clusterid in levels( sc10x@meta.data$seurat_clusters)){
  #   data_cluster = sc10x@meta.data[ which( sc10x@meta.data$seurat_clusters == clusterid), c( "seurat_clusters", "HTO_classification", module_name)]
  #   stats_df = wilcox_test( data = data_cluster, formula = as.formula( paste( module_name, "~ HTO_classification")))
  #   stats_df$seurat_clusters = clusterid
  #   module_stats_df = rbind( module_stats_df,
  #                            stats_df
  #                           )
  # }
  
  # # Adjust the p-value for multitesting with BH method to get FDR
  # module_stats_df$FDR = p.adjust( module_stats_df$p, method = "BH")
  # module_stats_df$signif = sapply( module_stats_df$FDR, function( fdr){ if( fdr < 0.05) "*" else "" })
  # module_stats_df = module_stats_df[ , c( "seurat_clusters", "group1", "group2", "n1", "n2", "p", "FDR", "signif")]
  # 
  # # Display the statistics result table
  # print(
  #   module_stats_df %>% 
  #     mutate_if(is.double, format, digits=4,nsmall = 0) %>% 
  #     kable( caption = paste( groupid, " : Comparison between conditions by clusters with Wilcoxon test")) %>% 
  #     kable_styling( full_width = FALSE)
  # )
  
}
