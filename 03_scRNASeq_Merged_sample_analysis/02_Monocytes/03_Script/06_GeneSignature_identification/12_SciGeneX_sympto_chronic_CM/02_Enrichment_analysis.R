
############################### GO over-representation test

# Total nuber of modules
n_modules <- length(res_cluster_subtype@gene_clusters)
gene_signatures <- list()
relevant_ontologies <- c("heat", "fold", "cytokine", "MT-")

for (module_nb in seq_len(n_modules)) {
    message("Processing Module ", module_nb)
    
    # Creating a directory
    module_dir <- file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb))
  if (!dir.exists(module_dir)) {
    dir.create(module_dir, recursive = TRUE)
  }

    # Extracting a gene_list from a specific module
    gene_list <- res_cluster_subtype@gene_clusters[[module_nb]]

    if (length(gene_list) < 5) {
    message("Skipping Module ", module_nb, " (too few genes)")
    next
  }
  
   # Skip module if mitochondrial genes are present
    if (any(grepl("^MT-", gene_list, ignore.case = TRUE))) {
  message("Skipping Module ", module_nb, " (contains mitochondrial genes)")
  next
    }  
    
  # Map SYMBOL -> ENTREZID
  mapped_gene_list <- mapIds(
    org.Hs.eg.db,
    keys = gene_list,
    keytype = "SYMBOL",
    column = "ENTREZID",
    multiVals = "first"
  )

    # Remove NAs (unmapped genes)
    mapped_gene_list <- na.omit(mapped_gene_list)

    if (length(mapped_gene_list) < 5) {
    message("Skipping Module ", module_nb, " (too few mapped genes)")
    next
  }
  
 ########################
  #### GO ENRICHMENT #####
  ########################
  
    
# BP
run_enrichment <- function(ontology) {
    enrichGO(
      gene          = mapped_gene_list,
      OrgDb         = org.Hs.eg.db,
      keyType       = "ENTREZID",
      ont           = ontology,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.05,
      readable      = TRUE
    )
  }
   go_results_BP <- run_enrichment("BP")
  go_results_MF <- run_enrichment("MF")


  ########################
  #### SAVE RESULTS ######
  ########################

  save_results <- function(go_results, ontology) {

    if (is.null(go_results) || nrow(as.data.frame(go_results)) == 0) {
      message("No significant enrichment for Module ", module_nb, " (", ontology, ")")
      return(NULL)
    }
    
# Save table
    write.csv(
      as.data.frame(go_results),
      file = file.path(
        module_dir,
        paste0(outputFilesPrefix, "GO_analysis_", ontology, ".csv")
      ),
      row.names = FALSE
    )
    
 ####### Plots #######

    p_bar  <- barplot(go_results, showCategory = 10,
                      title = paste("Module", module_nb, "-", ontology))

    p_dot  <- dotplot(go_results, showCategory = 10,
                      title = paste("Module", module_nb, "-", ontology))

    p_net  <- cnetplot(go_results, showCategory = 5)

    ggsave(
      file.path(module_dir, paste0(outputFilesPrefix, "barplot_", ontology, ".png")),
      plot = p_bar, width = 12, height = 8, dpi = 300
    )

    ggsave(
      file.path(module_dir, paste0(outputFilesPrefix, "dotplot_", ontology, ".png")),
      plot = p_dot, width = 12, height = 8, dpi = 300
    )

    ggsave(
      file.path(module_dir, paste0(outputFilesPrefix, "cnetplot_", ontology, ".png")),
      plot = p_net, width = 12, height = 8, dpi = 300
    )
    
 } 
  save_results(go_results_BP, "BP")
  save_results(go_results_MF, "MF")              
                
  
    # Filtering BP terms containing migration/adhesion
    if (!is.null(go_results_BP) && nrow(as.data.frame(go_results_BP)) > 0) {

    bp_df <- as.data.frame(go_results_BP)
    
    for (ont in relevant_ontologies) {
  
  # Identify relevant GO terms

  bp_terms <- bp_df[
    grepl(paste(ont, collapse = "|"),
           bp_df$Description,
           ignore.case = TRUE),
  ]
   if (nrow(bp_terms) > 0) {
  # Extract gene lists from selected GO terms
    relevant_genes <- strsplit(bp_terms$geneID, "/")

    # Flatten + remove duplicates
    relevant_genes <- unique(unlist(relevant_genes))

    # Check not empty
    if (length(relevant_genes) > 0) {
    # Append to the correct ontology
    gene_signatures[[ont]] <- c(gene_signatures[[ont]], relevant_genes)
     } else {
    message("No ", ont, " terms in Module ", module_nb)
  }
}
}
}
    # Remove duplicates
    gene_signatures <- lapply(gene_signatures, unique)

    # Save all signatures
    for (ont in names(gene_signatures)) {
  
    write.csv(
    data.frame(gene = gene_signatures[[ont]]),
    file = file.path(
      PATH_ANALYSIS_OUTPUT,
      paste0(outputFilesPrefix, ont, "_Gene_sign.csv")
    ),
    row.names = FALSE
  )
}
}


          
          
