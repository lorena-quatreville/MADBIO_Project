
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

                
                

              #######  Visualise results with a plot
                
                

                                    # all genes

# BP
if (!is.null(go_results_BP) && nrow(as.data.frame(go_results_BP)) > 0) {
barplot <- barplot(go_results_BP, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in BP")
print(barplot)
dotplot <- dotplot(go_results_BP, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in BP") 
print(dotplot)
networkplot <- cnetplot(go_results_BP, showCategory = 5)
print(networkplot)

ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "barplot_GO_analysis_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "dotplot_GO_analysis_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "cnetplot_GO_analysis_BP.png")), width = 12, height = 8, dpi=300)
}

# MF
if (!is.null(go_results_MF) && nrow(as.data.frame(go_results_MF)) > 0) {
barplot <- barplot(go_results_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF")
print(barplot)
dotplot <- dotplot(go_results_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF")
print(dotplot)
networkplot <- cnetplot(go_results_MF, showCategory = 5)
print(networkplot)


ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "barplot_GO_analysis_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "dotplot_GO_analysis_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "cnetplot_GO_analysis_MF.png")), width = 12, height = 8, dpi=300)
}

