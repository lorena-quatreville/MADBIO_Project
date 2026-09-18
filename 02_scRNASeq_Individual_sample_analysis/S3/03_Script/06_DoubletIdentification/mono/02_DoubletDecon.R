
# ##############################################################
       #     Doublet Detection Script: DoubletDecon   #
# ##############################################################



# Perform a FindAllMarkers() analysis to retrieve marker genes that will be used by DoubletDecon

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
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME))
dir.create(file.path(PATH_ANALYSIS_OUTPUT))
dir.create(file.path(PATH_ANALYSIS_EXTRA_OUTPUT))
write.table( markers,
             file= file.path( PATH_ANALYSIS_EXTRA_OUTPUT, paste0( outputFilesPrefix, "MarkerGenes.tsv")),
             quote = FALSE,
             row.names = TRUE, 
             col.names = NA, # Add a blank column name for row names (CSV convention)
             sep="\t");


## Formatting data to pas it to DoubletDecon

# Number of top markers per cluster to use
num_genes <- 50

# Keep all columns, but only top N genes per cluster
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n = num_genes, with_ties = FALSE) %>%
  ungroup()
  
top_markers$gene <- as.character(top_markers$gene)
top_markers$cluster <- as.character(top_markers$cluster)

# Rename explicitly
names(top_markers)[names(top_markers) == "avg_log2FC"] <- "avg_logFC"

# top_markers has to a df
top_markers_df <- as.data.frame(top_markers)

# Make sure rownames are numeric (like original FindAllMarkers output)
rownames(top_markers_df) <- 1:nrow(top_markers_df)

# Add an artificial unique first column
top_markers_df <- cbind(ID = 1:nrow(top_markers_df), top_markers_df)


# Top markers file (from your FindAllMarkers): # Write to file
genes_file <- file.path(PATH_ANALYSIS_EXTRA_OUTPUT, "Top50Genes_file_DoubletDecon.txt")
write.table(top_markers_df,
            file = genes_file,
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)


# DoubletDecon expects an expression matrix, cluster labels and markers per clusters:



# Creating the 3 files to feed Seurat pre-processed:

# Expression matrix file (log-normalized)
expr_mat <- as.matrix(GetAssayData(sc10x, slot = "data"))  # log-normalized
expr_file <- file.path(PATH_ANALYSIS_EXTRA_OUTPUT, "expression_file_DoubletDecon.txt")
write.table(expr_mat,
            file = expr_file,
            sep = "\t",
            quote = FALSE,
            col.names = NA)

# Cluster file
clusters <- sc10x$seurat_clusters           # cluster IDs
clusters_file <- file.path(PATH_ANALYSIS_EXTRA_OUTPUT, "clusters_file_DoubletDecon.txt")
clusters_df <- data.frame(Cell = names(clusters),
                          x = as.character(clusters))
write.table(clusters_df,
            file = clusters_file,
            sep = "\t",
            quote = TRUE,
            row.names = FALSE,
            col.names = c("", "x")) # match the expected header


# test
read.table(expr_file, header = TRUE, sep = "\t")
read.table(genes_file, header = TRUE, sep = "\t")
read.table(clusters_file, header = TRUE, sep = "\t")



# Running Seurat pre-process to have 3 formated inputs for Main_Doublet_Decon
DD_input <- Seurat_Pre_Process(
  expressionFile = expr_file,
  genesFile = genes_file,
  clustersFile = clusters_file
)


# Remove first row if it’s the label
if (rownames(DD_input$newExpressionFile)[1] == "column_clusters-flat") {
  DD_input$newExpressionFile <- DD_input$newExpressionFile[-1, ]
}
if (rownames(DD_input$newFullExpressionFile)[1] == "column_clusters-flat") {
  DD_input$newFullExpressionFile <- DD_input$newFullExpressionFile[-1, ]
}
# Remove first column if it’s the label
if (colnames(DD_input$newExpressionFile)[1] == "row_clusters-flat") {
  DD_input$newExpressionFile <- DD_input$newExpressionFile[, -1]
}
if (colnames(DD_input$newFullExpressionFile)[1] == "row_clusters-flat") {
  DD_input$newFullExpressionFile <- DD_input$newFullExpressionFile[, -1]
}

# Saving genrated files to pass on a file path to doubletDecon
newExpressionFile <- file.path(PATH_ANALYSIS_EXTRA_OUTPUT, "_DD_input_ExpressionFile.txt")
newFullExpressionFile <- file.path(PATH_ANALYSIS_EXTRA_OUTPUT, "_DD_input_FullExpressionFile.txt")
newGroupsFile <- file.path(PATH_ANALYSIS_EXTRA_OUTPUT, "_DD_input_GroupsFile.txt")

write.table(DD_input$newExpressionFile, file = newExpressionFile, sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(DD_input$newFullExpressionFile, file = newFullExpressionFile, sep="\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
write.table(DD_input$newGroupsFile, file = newGroupsFile , sep="\t", row.names = TRUE, col.names = FALSE)

# for DoubletDecon to save files
DD_path = file.path(PATH_ANALYSIS_OUTPUT, cellType) # tricks Main_Doublet_Decon function so that it stores output files in the last folder of the path (cellType folder)

# Running DoubletDecon
DD_results <- Main_Doublet_Decon(
  rawDataFile = newExpressionFile,        # supply the matrix directly
  groupsFile  = newGroupsFile,         # cluster assignments
  filename = "DD_results",
  location = DD_path,
  fullDataFile = newFullExpressionFile,   # can reuse the same matrix
  removeCC = FALSE,
  species = "hsa",
  rhop = 1,
  write = TRUE,
  PMF = TRUE, # Pseudo Marker Finder: t checks whether the putative doublets actually express unique marker genes that aren’t just a mixture of two clusters (limit false positive doublets)
  useFull = FALSE, # FALSE: only use the marker genes list (faster, but limited)
  heatmap = FALSE, #FALSE: skip heatmap plotting (common in large datasets)
  centroids = TRUE # Centroids are usually more stable for large clusters, while medoids can be more robust to outliers in small clusters
)

# Generating two files for Doublet Decon results from raw outputs

# Load results
final_groups <- read.table(
  file.path(PATH_ANALYSIS_OUTPUT, paste0(cellType, "Final_doublets_groups_DD_results.txt")),
  header = FALSE,
  sep = "\t",
  quote = "\"",
  stringsAsFactors = FALSE,
  skip = 1
)
                                                     
                           
colnames(final_groups) <- c("CellBarcode", "Cluster", "Group")

# Replacing ".1" by "-1" in DD output file cell barcode so that they are comparable to seurat cell barcode
final_groups$CellBarcode <- sub("\\.1$", "-1", final_groups$CellBarcode)


# Retrieving cell barcodes classified as doublet 
doublets <- final_groups$CellBarcode

# All barcodes
all_cells <- colnames(sc10x) 



# Build classification table
classification <- data.frame(
  Barcode = all_cells,
  DoubletDecon_DropletType = ifelse(all_cells %in% doublets, "Doublet", "Singlet")
)

# Save as DoubletDecon_doublets_singlets.tsv
write.table(classification, file.path(PATH_ANALYSIS_OUTPUT, "DoubletDecon_doublets_singlets.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

# Build summary table
summary_table <- as.data.frame(table(classification$DoubletDecon_DropletType))
colnames(summary_table) <- c("Classification", "Droplet_N")

write.table(summary_table, file.path(PATH_ANALYSIS_OUTPUT, "DoubletDecon_doublet_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)


# Saving DoubletDecon results in the Seurat object
# Create a named vector
DD.prediction = setNames(classification$DoubletDecon_DropletType, classification$Barcode)
sc10x <- AddMetaData(sc10x, metadata = DD.prediction, col.name = "DD_prediction")




