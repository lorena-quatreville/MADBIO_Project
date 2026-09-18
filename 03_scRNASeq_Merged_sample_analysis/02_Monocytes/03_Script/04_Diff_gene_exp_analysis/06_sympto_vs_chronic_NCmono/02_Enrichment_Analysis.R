########################                                                                   ##################################
                         #  Enrichment analysis (GO enrichment test and GSEA)  #
########################                                                                   ##################################





                            ########## GO overrepresentation test ##############
                            
# By taking all significant up AND down regulated genes

# Creating a gene_list extracted from the deg_filtered matrix
gene_list <- deg_filtered$gene

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

write.csv(as.data.frame(go_results_BP), file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "GO_analysis_results_deg_BP.csv")), row.names = FALSE)
write.csv(as.data.frame(go_results_MF), file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "GO_analysis_results_deg_MF.csv")), row.names = FALSE)



# By taking only significant up OR down regulated genes


# Separate upregulated and downregulated genes
upregulated_genes <- deg_filtered$gene[deg_filtered$avg_log2FC > 0]
# Save the sorted results to a CSV file
upregulated_genes_df <- deg_filtered[deg_filtered$avg_log2FC > 0,]
write.csv(upregulated_genes_df , file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "deg_res_up.csv")), row.names = TRUE)


downregulated_genes <- deg_filtered$gene[deg_filtered$avg_log2FC < 0]
# Save the sorted results to a CSV file
downregulated_genes_df <- deg_filtered[deg_filtered$avg_log2FC < 0,]
write.csv(downregulated_genes_df  , file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "deg_res_down.csv")), row.names = TRUE)


# Convert to Entrez IDs
mapped_upregulated <- na.omit(mapIds(org.Hs.eg.db, keys = upregulated_genes, keytype = "SYMBOL", column = "ENTREZID"))
mapped_downregulated <- na.omit(mapIds(org.Hs.eg.db, keys = downregulated_genes, keytype = "SYMBOL", column = "ENTREZID"))

# GO enrichment for upregulated genes

# BP
go_results_up_BP <- enrichGO(
  gene          = mapped_upregulated,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# MF
go_results_up_MF <- enrichGO(
  gene          = mapped_upregulated,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# GO enrichment for downregulated genes
# BP
go_results_down_BP <- enrichGO(
  gene          = mapped_downregulated,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# MF
go_results_down_MF <- enrichGO(
  gene          = mapped_downregulated,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# View results
head(go_results_up_BP)
head(go_results_up_MF)
head(go_results_down_BP)
head(go_results_down_MF)

write.csv(as.data.frame(go_results_up_BP), file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "GO_analysis_results_upReg_BP.csv")), row.names = FALSE)
write.csv(as.data.frame(go_results_up_MF), file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "GO_analysis_results_upReg_MF.csv")), row.names = FALSE)
write.csv(as.data.frame(go_results_down_BP), file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "GO_analysis_results_downReg_BP.csv")), row.names = FALSE)
write.csv(as.data.frame(go_results_down_MF), file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "GO_analysis_results_downReg_MF.csv")), row.names = FALSE)



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

ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "barplot_GO_analysis_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GO_analysis_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cnetplot_GO_analysis_BP.png")), width = 12, height = 8, dpi=300)
}

# MF
if (!is.null(go_results_MF) && nrow(as.data.frame(go_results_MF)) > 0) {
barplot <- barplot(go_results_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF")
print(barplot)
dotplot <- dotplot(go_results_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF")
print(dotplot)
networkplot <- cnetplot(go_results_MF, showCategory = 5)
print(networkplot)


ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "barplot_GO_analysis_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GO_analysis_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cnetplot_GO_analysis_MF.png")), width = 12, height = 8, dpi=300)
}


                                # up-regulated genes


# BP
if (!is.null(go_results_up_BP) && nrow(as.data.frame(go_results_up_BP)) > 0) {
barplot <- barplot(go_results_up_BP, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in BP for up-regulated genes")
print(barplot)
dotplot <- dotplot(go_results_up_BP, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in BP for up-regulated genes")
print(dotplot)
networkplot <- cnetplot(go_results_up_BP, showCategory = 5)
print(networkplot)

ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "barplot_GO_analysis_upReg_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GO_analysis_upReg_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cnetplot_GO_analysis_upReg_BP.png")), width = 12, height = 8, dpi=300)
}

# MF
if (!is.null(go_results_up_MF) && nrow(as.data.frame(go_results_up_MF)) > 0) {
barplot <- barplot(go_results_up_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF for up-regulated genes")
print(barplot)
dotplot <- dotplot(go_results_up_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF for up-regulated genes")
print(dotplot)
networkplot <- cnetplot(go_results_up_MF, showCategory = 5)
print(networkplot)

ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "barplot_GO_analysis_upReg_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GO_analysis_upReg_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cnetplot_GO_analysis_upReg_MF.png")), width = 12, height = 8, dpi=300)
}


                                        # down-regulated genes

# BP
if (!is.null(go_results_down_BP) && nrow(as.data.frame(go_results_down_BP)) > 0) {
barplot <- barplot(go_results_down_BP, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in BP for down-regulated genes")
print(barplot)
dotplot <- dotplot(go_results_down_BP, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in BP for down-regulated genes")
print(dotplot)
networkplot <- cnetplot(go_results_down_BP, showCategory = 5)
print(networkplot)

ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "barplot_GO_analysis_downReg_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GO_analysis_downReg_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cnetplot_GO_analysis_downReg_BP.png")), width = 12, height = 8, dpi=300)
}


# MF
if (!is.null(go_results_down_MF) && nrow(as.data.frame(go_results_down_MF)) > 0) {
barplot <- barplot(go_results_down_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF for down-regulated genes")
print(barplot)
dotplot <- dotplot(go_results_down_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF for down-regulated genes")
print(dotplot)
networkplot <- cnetplot(go_results_down_MF, showCategory = 5)
print(networkplot)

ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "barplot_GO_analysis_downReg_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GO_analysis_downReg_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cnetplot_GO_analysis_downReg_MF.png")), width = 12, height = 8, dpi=300)
}





                                                    ########## GSEA analysis ##############




# Extract log2 fold change and set gene names
gene_fc <- deg_results$avg_log2FC  # Extract log2 fold change
names(gene_fc) <- rownames(deg_results)  # Use gene symbols as names

# Convert gene symbols to Entrez IDs
mapped_gene_fc <- mapIds(org.Hs.eg.db, 
                         keys = names(gene_fc), 
                         keytype = "SYMBOL", 
                         column = "ENTREZID", 
                         multiVals = "first")

# Remove NAs and create a named vector
gene_fc <- gene_fc[!is.na(mapped_gene_fc)]
names(gene_fc) <- na.omit(mapped_gene_fc)

# Sort genes by log fold change (highly expressed at the top)
gene_fc <- sort(gene_fc, decreasing = TRUE)

# Print the top genes
head(gene_fc)

gse_results <- gseGO(
  geneList      = gene_fc,      # Ranked gene list
  OrgDb         = org.Hs.eg.db, # Human annotation database
  keyType       = "ENTREZID",   # Gene ID type
  ont           = "BP",         # GO category (BP, MF, or CC)
  minGSSize     = 10,           # Minimum genes in a gene set
  maxGSSize     = 500,          # Maximum genes in a gene set
  pvalueCutoff  = 0.05,         # P-value threshold
  pAdjustMethod = "BH",        # Multiple testing correction
  verbose       = FALSE
)

# View top enriched GO terms
head(gse_results)

# Convert gse_results to a data frame
gse_df <- as.data.frame(gse_results)

write.csv(gse_df, file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "GSEA_results.csv")), row.names = FALSE)


                                 ######## VISUALISATION #######



# Visualisation: (only top 10 most significant BP)
if (!is.null(gse_results) && nrow(as.data.frame(gse_results)) > 0) {
ridgeplot <- ridgeplot(gse_results, showCategory = 10) + 
  labs(title = "GSEA Ridge Plot")
print(ridgeplot)
ggsave(plot = ridgeplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "ridgeplot_GSEA_results.png")), width = 12, height = 10, dpi=300)


# doplot to visiualise p-values of enrichment
dotplot <- dotplot(gse_results, showCategory = 10, title = "GSEA Dot Plot")
print(dotplot)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GSEA_results_p_values.png")), width = 12, height = 8, dpi=300)


# same but colored by enrichment score
dotplot <- dotplot(gse_results, showCategory = 10, title = "GSEA Dot Plot", color = "NES") # red/orange positive enrichment score and blue/purple negative enrichment score (downregulated in the first indivdual compared to the second)
print(dotplot)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GSEA_results_enrich_score.png")), width = 12, height = 8, dpi=300)

# network plot
# Converting back ENTREZID to gene symbols
# Set readable gene names (convert Entrez IDs to gene symbols)
gse_results <- setReadable(gse_results, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
# Create the CnetPlot with readable gene names
networkplot <- cnetplot(gse_results, showCategory = 5, node_label = "all")
print(networkplot)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cnetplot_GSEA_results.png")), width = 12, height = 8, dpi=300)
}


# or doing the same by directly filtering the genes based on the enrichment score:
# Show only UPregulated GO terms

# Filter for Upregulated GO Terms (NES > 0)
gse_results_up <- gse_results
gse_results_up@result <- gse_results_up@result %>%
  filter(NES > 0) %>%
  arrange(desc(NES))  # Sort by NES
  
# Plot the Top 10 Upregulated GO Terms
if (!is.null(gse_results_up) && nrow(as.data.frame(gse_results_up)) > 0) {
dotplot <- dotplot(gse_results_up, showCategory = 10, 
                      title = "Upregulated GO Terms in Individual 1")
print(dotplot)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GSEA_results_upReg.png")), width = 12, height = 8, dpi=300)

# upregulated genes only in newtwork plot
# Set readable gene names (convert Entrez IDs to gene symbols)
gse_results_up <- setReadable(gse_results_up, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
# Create the CnetPlot with readable gene names
networkplot <- cnetplot(gse_results_up, showCategory = 5, node_label = "all")
print(networkplot)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cnetplot_GSEA_results_upReg.png")), width = 12, height = 8, dpi=300)
}

# Show only DOWNregulated GO terms

# Filter for Downregulated GO Terms (NES > 0)
gse_results_down <- gse_results
gse_results_down@result <- gse_results_down@result %>%
  filter(NES < 0) %>%
  arrange(desc(NES))  # Sort by NES
  
# Plot the Top 10 Upregulated GO Terms
if (!is.null(gse_results_down) && nrow(as.data.frame(gse_results_down)) > 0) {
dotplot <- dotplot(gse_results_down, showCategory = 10, 
                      title = "Downregulated GO Terms in Individual 1")
print(dotplot)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "dotplot_GSEA_results_downReg.png")), width = 12, height = 8, dpi=300)

# with downregulated genes only in network plot
# Set readable gene names (convert Entrez IDs to gene symbols)
gse_results_down <- setReadable(gse_results_down, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
# Create the CnetPlot with readable gene names
networkplot <- cnetplot(gse_results_down, showCategory = 5, node_label = "all")
print(networkplot)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cnetplot_GSEA_results_downReg.png")), width = 12, height = 8, dpi=300)
}






# A voir si utile: Perform pathway activity scoring (e.g., AUCell, GSVA) to measure activation levels of the process in each sample.








