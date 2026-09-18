#####################################################
                # Visualisation of the doublets
####################################################


# Create a consensus column
sc10x$doublet_consensus <- "Singlet"
sc10x$doublet_consensus[
  sc10x$DF_prediction == "Doublet" & sc10x$DD_prediction == "Doublet"
] <- "Consensus Doublet"

# Optional: mark cells predicted as doublet by only one method
sc10x$doublet_consensus[
  sc10x$DF_prediction == "Doublet" & sc10x$DD_prediction != "Doublet"
] <- "DF-only Doublet"

sc10x$doublet_consensus[
  sc10x$DD_prediction == "Doublet" & sc10x$DF_prediction != "Doublet"
] <- "DD-only Doublet"

# Plotting a UMAP for DoubletFinder prediction done on the subset of cells

DF_umap <- DimPlot(sc10x, group.by = "DF_prediction", reduction = "umap") +
  scale_color_manual(values = c("gray", "red")) +
  ggtitle("DoubletFinder on a cell subset Predictions")
  
#save plot 
ggsave(plot = DF_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, cellType, "_DF_prediction_umap.png"))
       , width = 23, height = 13, units ="cm")

# Plotting a UMAP for DoubletDecon prediction

DD_umap <- DimPlot(sc10x, group.by = "DD_prediction", reduction = "umap") +
  scale_color_manual(values = c("gray", "blue")) +
  ggtitle("DoubletDecon Predictions")

#save plot 
ggsave(plot = DD_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, cellType, "_DD_prediction_umap.png"))
       , width = 23, height = 13, units ="cm")
  
# Plotting a UMAP for consensus (comparison)

consensus_umap <- DimPlot(sc10x, group.by = "doublet_consensus", reduction = "umap") +
  scale_color_manual(values = c(
    "Singlet" = "gray",
    "DF-only Doublet" = "red",
    "DD-only Doublet" = "blue",
    "Consensus Doublet" = "purple"
  )) +
  ggtitle("Doublet Comparison: DF vs DD")
  
#save plot 
ggsave(plot = consensus_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, cellType, "_consensus_prediction_umap.png"))
       , width = 23, height = 13, units ="cm")  

# Plotting a UMAP for DoubletFinder prediction done on all the dataset (not with the cluster identified)

DF_umap <- DimPlot(sc10x, group.by = "multiplet_prediction", reduction = "umap") +
  scale_color_manual(values = c("gray", "red")) +
  ggtitle("DoubletFinder on all cells Predictions")
  
#save plot 
ggsave(plot = DF_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, cellType, "_DF_multiplet_prediction_all_umap.png"))
       , width = 23, height = 13, units ="cm")




  
# Summary

table_DF_prediction = table(sc10x@meta.data[["DF_prediction"]])
table_DD_prediction = table(sc10x[["DD_prediction"]])
table_doublet_consensus = table(sc10x[["doublet_consensus"]])
table_multiplet_prediction = table(sc10x[["multiplet_prediction"]])

write.table(table_DF_prediction, file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, cellType,"_table_DF_prediction.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(table_DD_prediction, file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, cellType,"_table_DD_prediction.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(table_doublet_consensus, file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, cellType, "_table_doublet_consensus.tsv")), #(DF_prediction and DD_prediction)
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(table_multiplet_prediction, file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, cellType, "_table_multiplet_prediction.tsv")),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Saving the updated seurat object
saveRDS(sc10x, file = PATH_SEURAT_RDS_FILE)

