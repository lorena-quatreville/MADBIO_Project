# ########################################################
# This script aims to analyse the heterogeneity of cells
#  of cells using dimension reduction techniques
# ########################################################


# CLUSTER FUNCTIONAL ENRICHMENT
################################

## @knitr heterogeneity_markerGenes_functional_enrichment

for( clusterName in levels( markers$cluster)){

  # Add a sub-tab with the Cluster name
  # ...............................................................................
  cat( "\n \n")
  cat( "\n##### Cluster", clusterName, "{.tabset}")
  cat( "\n \n")
  
  
  # Add a sub-tab with the enrichment of the positive markers in the GO BP terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n###### Markers genes GO Biological Process")
  cat( "\n \n")
  
  positive_ego <- enrichGO(gene          = tolower(markers[ which( markers$cluster == clusterName), "gene"]),
                           universe      = tolower(rownames( sc10x@assays$RNA)),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = "SYMBOL",
                           ont           = "BP",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                           qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
  
  selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
  
  if( selected_term_number> 0){  
    
    print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of marker genes in BP GO terms"))
    positive_ego_df = data.frame( positive_ego)
    
    positive_ego_simMatrix <- calculateSimMatrix( positive_ego_df$ID,
                                                  orgdb="org.Mm.eg.db",
                                                  ont="BP",
                                                  method="Rel")
    
    positive_ego_scores <- setNames(-log10( positive_ego_df$qvalue), positive_ego_df$ID)
    positive_ego_reducedTerms <- reduceSimMatrix(positive_ego_simMatrix,
                                                 positive_ego_scores,
                                                 threshold=0.7,
                                                 orgdb="org.Mm.eg.db")
    
    print( scatterPlot( positive_ego_simMatrix, positive_ego_reducedTerms))
    
    treemapPlot( positive_ego_reducedTerms)
  }
  else{
    cat("<BR>No result<BR>")
  }
  
  
  # Add a sub-tab with the enrichment of the positive markers in the GO MF terms
  # ...............................................................................
  cat( "\n \n")
  cat( "\n###### Markers genes GO Molecular Function")
  cat( "\n \n")
  
  positive_ego <- enrichGO(gene          = markers[ which( markers$cluster == clusterName), "gene"],
                           universe      = rownames( sc10x@assays$RNA),
                           OrgDb         = org.Mm.eg.db,
                           keyType       = "SYMBOL",
                           ont           = "MF",
                           pAdjustMethod = "BH",
                           pvalueCutoff  = ENRICHMENT_GO_PVALUECUTOFF,
                           qvalueCutoff  = ENRICHMENT_GO_QVALUECUTOFF)
  
  selected_term_number = length( which( positive_ego@result$p.adjust <= positive_ego@pvalueCutoff & positive_ego@result$qvalue <= positive_ego@qvalueCutoff))
  
  if( selected_term_number> 0){  
    
    print( dotplot( positive_ego, showCategory=10) + ggtitle("Enrichment of positive markers in MF GO terms"))
    positive_ego_df = data.frame( positive_ego)
    
    positive_ego_simMatrix <- calculateSimMatrix( positive_ego_df$ID,
                                                  orgdb="org.Mm.eg.db",
                                                  ont="MF",
                                                  method="Rel")
    
    positive_ego_scores <- setNames(-log10( positive_ego_df$qvalue), positive_ego_df$ID)
    positive_ego_reducedTerms <- reduceSimMatrix(positive_ego_simMatrix,
                                                 positive_ego_scores,
                                                 threshold=0.7,
                                                 orgdb="org.Mm.eg.db")
    
    print( scatterPlot( positive_ego_simMatrix, positive_ego_reducedTerms))
    
    treemapPlot( positive_ego_reducedTerms)
  }
  else{
    cat("<BR>No result<BR>")
  }
}

