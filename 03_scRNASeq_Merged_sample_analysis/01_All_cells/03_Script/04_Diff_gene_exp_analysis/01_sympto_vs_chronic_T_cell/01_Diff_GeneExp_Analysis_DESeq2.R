########################                                                                   ##################################
                         #  Differential gene expression analysis using FindMarkers()  #
########################                                                                   ##################################


library(zinbwave)
library(DESeq2)
library(SingleCellExperiment)



## If needed to load the merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x) <- sc10x$sgx_clusters

# Subsetting T cell
sc10x_subset <- subset(
  sc10x,
  idents = grep("T cell", levels(Idents(sc10x)), value = TRUE)
)

# Removing doublets
sc10x_subset <- subset(
  sc10x_subset,
  subset = !grepl("doublet|triplet", Idents(sc10x_subset), ignore.case = TRUE) &
           orig.ident %in% c(sample_id_A, sample_id_B)
)

# Set cluster identity for this subset
Idents(sc10x_subset) <- sc10x_subset$orig.ident

# Running FindMarkers to find differentially expressed genes:
deg_results <- FindMarkers(sc10x_subset, ident.1 = sample_id_A , ident.2 = sample_id_B, test.use = "MAST", latent.vars = c("nCount_RNA", "batch"))

# View top differentially expressed genes
head(deg_results, 10)

# Sort the results by p-value (ascending order)
deg_results_sorted <- deg_results[order(deg_results$p_val), ]

# Create the output directory
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME))
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME))


# Save the sorted results to a CSV file
write.csv(deg_results_sorted, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "deg_res.csv")), row.names = TRUE)





# Filtering of only significant differentially expressed genes

deg_filtered <- deg_results[deg_results$p_val_adj < 0.05 & abs(deg_results$avg_log2FC) > 0.5, ]
head(deg_filtered, 10)

# Add a small constant to the adjusted p-values and create a new column 'p_val_adj_min'
deg_filtered$p_val_adj_min <- deg_filtered$p_val_adj + 1e-303

deg_filtered$gene <- rownames(deg_filtered)

# Select top 10 genes with the smallest adjusted p-values
top_10_genes <- deg_filtered[order(deg_filtered$p_val_adj_min), ][1:10, ]


## Visualisation

# Volcano plot



volcano_plot <- ggplot(deg_filtered, aes(x = avg_log2FC, y = -log10(p_val_adj_min))) +
    geom_point(aes(color = p_val_adj_min < 0.05), alpha = 0.5) +
    scale_color_manual(values = c("gray", "red")) +
    theme_minimal() +
    labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value") +
    geom_text_repel(data = top_10_genes, aes(label = gene),
                    size = 3, color = "blue", max.overlaps = 100) +
    scale_x_continuous(limits = c(-max(abs(deg_filtered$avg_log2FC)), max(abs(deg_filtered$avg_log2FC))))


ggsave(plot = volcano_plot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "VolcanoPlot_deg_res.png")), width = 12, height = 8, dpi = 300)

# heatmap
# Take the top 20 DEGs
deg_filtered_top_20 <- deg_filtered[order(deg_filtered$p_val_adj_min),][1:20,]
top_20_genes <- rownames(deg_filtered_top_20[order(deg_filtered_top_20$avg_log2FC), ])


# Ordering samples
sample_order <- c(
  "08_scRNASeq_A82_S_tp2",
  "12_scRNASeq_A86_S_tp2",
  "15_scRNASeq_A76_S_tp2",
  "04_scRNASeq_A81_AS_tp2",
  "16_scRNASeq_A83_AS_tp2",
  "13_scRNASeq_A62_AS_tp2",
  "14_scRNASeq_A75_AS_tp2"
)

# Setting factor levels in the seurat object
sc10x_subset$orig.ident <- factor(
  sc10x_subset$orig.ident,
  levels = sample_order
)

# Plot heatmap directly from Seurat
heatmap_plot <- DoHeatmap(sc10x_subset, features = top_20_genes, group.by = "orig.ident") +
  ggtitle("Heatmap of Top Differentially Expressed Genes") +
  theme(axis.text.y = element_text(size = 8))  # smaller text for readability

print(heatmap_plot)

# Save to file
ggsave(
  plot = heatmap_plot,
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "Heatmap_top_DEGs.png")),
  width = 10, height = 8, dpi = 300
)


# Plot heatmap with genes ordered by expression and cells by their similarity in expression across the selected genes

# Ordering genes by expression accross all cells in subsets
# Average expression across all cells
avg_exp <- Matrix::rowMeans(GetAssayData(sc10x_subset, slot = "scale.data")[top_20_genes, , drop = FALSE])

# Reorder genes by descending expression
top_20_genes_ordered <- names(sort(avg_exp, decreasing = TRUE))

# Ranking cells by pseudo-bulk expression
# Get scaled expression of top 20 genes
expr_mat <- GetAssayData(sc10x_subset, slot = "scale.data")[top_20_genes_ordered, ]

# Hierarchical clustering of cells
cell_order <- hclust(dist(t(expr_mat)))$order
cell_names_ordered <- colnames(expr_mat)[cell_order]

# Expression matrix for the top 20 ordered genes
expr_mat_sub <- expr_mat[, cell_names_ordered]

# Annotation for sample IDs
cell_annot <- data.frame(Sample = sc10x_subset$orig.ident[cell_names_ordered])
rownames(cell_annot) <- cell_names_ordered

# Plot
heatmap_ordered <- pheatmap(expr_mat_sub,
         cluster_rows = FALSE,   # keep genes in given order
         cluster_cols = FALSE,   # keep cells in given order
         annotation_col = cell_annot,
         show_colnames = FALSE,
         fontsize_row = 8,
         main = "Top 20 DEGs Heatmap (ranked)")
print(heatmap_ordered)
# Save to file
ggsave(
  plot = heatmap_ordered,
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "Heatmap_ordered_top_DEGs.png")),
  width = 10, height = 8, dpi = 300
)












