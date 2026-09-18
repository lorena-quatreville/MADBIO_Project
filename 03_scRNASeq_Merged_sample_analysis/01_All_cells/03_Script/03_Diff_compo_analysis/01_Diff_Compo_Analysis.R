##### Using DCATS #####


## If needed to load the integrated merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

# Removing spaces after underscores
sc10x$sgx_clusters <- gsub("\\s*_\\s*", "_", sc10x$sgx_clusters)

# Remove accidental line breaks in cluster names
sc10x$sgx_clusters <- gsub("\n", " ", sc10x$sgx_clusters)
sc10x$sgx_clusters <- trimws(sc10x$sgx_clusters)  # Remove leading/trailing spaces
sc10x$sgx_clusters <- gsub("\\s+", " ", sc10x$sgx_clusters)  # Replace multiple spaces with single space
sc10x$sgx_clusters <- as.factor(sc10x$sgx_clusters)  # Ensure it's a factor

# compute KNN similaruty matrix
knn_mat <- knn_simMat(sc10x@graphs$RNA_snn, sc10x$sgx_clusters)  # Adjust if using another graph
print(knn_mat)


# Construct a count matrix
count_mat <- table(sc10x$orig.ident, sc10x$sgx_clusters)
print(count_mat)

# to take into account the condition
# Reorder condition_info so it matches rownames(count_mat)
condition_vec <- condition_info[rownames(count_mat)]

# Build the design matrix taking into account the batch effect
design_mat <- data.frame(
   condition = condition_info[rownames(count_mat)],
   batch = batch_info[rownames(count_mat)]
)
print(design_mat)


# Run dacats_GLM()
glm_results <- dcats_GLM(count_mat, design_mat, similarity_mat = knn_mat)
print(glm_results)

# Convert the list of lists into a single data frame
glm_results_df <- data.frame(
  Cell_Type = rownames(glm_results$coeffs),  # Keep only one column for cell types
  Coefficients = glm_results$ceoffs[, 1],
  Coeff_Error = glm_results$coeffs_err[, 1],
  LR_Values = glm_results$LR_vals[, 1],
  LRT_pvals = glm_results$LRT_pvals[, 1],
  FDR = glm_results$fdr[, 1]
)

# Creat the new direcotory
dir.create(file.path(PATH_ANALYSIS_OUTPUT))

# Save to CSV file
write.csv(glm_results_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glm_results_batch.csv")), row.names = FALSE)



## Doing the same but simplifying clusters

# Simplify cluster labels
sc10x$cluster_simplified <- case_when(
  str_detect(sc10x$sgx_clusters, "Classical monocyte") ~ "Classical monocyte",
  str_detect(sc10x$sgx_clusters, "Non-classical monocyte") ~ "Non-classical monocyte",
  str_detect(sc10x$sgx_clusters, "Intermediate monocyte") ~ "Intermediate monocyte",
  TRUE ~ sc10x$sgx_clusters
)

# Ensure the simplified cluster names are factors
sc10x$cluster_simplified <- as.factor(sc10x$cluster_simplified)

# Compute KNN similarity matrix
knn_mat <- knn_simMat(sc10x@graphs$RNA_snn, sc10x$cluster_simplified)
print(knn_mat)

# Construct a count matrix with the simplified clusters
count_mat <- table(sc10x$orig.ident, sc10x$cluster_simplified)
print(count_mat)

# Run DCATS differential composition analysis
glm_results <- dcats_GLM(count_mat, design_mat, similarity_mat = knn_mat)
print(glm_results)

# Convert the list of lists into a single data frame
glm_results_df <- data.frame(
  Cell_Type = rownames(glm_results$coeffs),  # Keep only one column for cell types
  Coefficients = glm_results$coeffs[, 1],
  Coeff_Error = glm_results$coeffs_err[, 1],
  LR_Values = glm_results$LR_vals[, 1],
  LRT_pvals = glm_results$LRT_pvals[, 1],
  FDR = glm_results$fdr[, 1]
)

# Create the output directory if it doesn't exist
dir.create(file.path(PATH_ANALYSIS_OUTPUT), showWarnings = FALSE)

# Save results to a CSV file
write.csv(glm_results_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glm_results_cluster_simplified_batch.csv")), row.names = FALSE)




# to visualise the data with a boxplot (and significance bar)
# Extracting the number of monocytes for each individual
classical_counts <- count_mat[, "Classical monocyte"]

# Retrieving the condition column in the design_mat
conditions <- design_mat$condition[match(rownames(count_mat), rownames(design_mat))]

# Build a data frame for plotting
glm_results_plot <- data.frame(
  Sample = rownames(count_mat),
  Condition = conditions,
  Count = classical_counts
)

# Optional: compute proportion if you want
glm_results_plot$Proportion <- glm_results_plot$Count / rowSums(count_mat)

# Get DCATS FDR for annotation
fdr_val <- glm_results_df$FDR[glm_results_df$Cell_Type == "Classical monocyte"]

# Custom colors
cols <- c("Chronic" = "orange", "Symptomatic" = "red")

# Dot plot per sample with group mean + FDR annotation
glm_results_dotplot <- ggplot(glm_results_plot, aes(x = Condition, y = Proportion, color = Condition)) +
  geom_jitter(width = 0.1, size = 5, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.8,
               color = "black", fatten = 0, alpha = 0.8) +
     # --- Significance bar with DCATS p-value ---
  geom_signif(
    comparisons = list(c("Chronic", "Symptomatic")),
    annotations = paste0("p = ", formatC(fdr_val, format = "e", digits = 1)),
    y_position = max(glm_results_plot$Proportion) * 1.08,
    tip_length = 0.01,
    textsize = 7,
    color = "black"     # <- makes the bar black
  ) +
  
  # Colors
  scale_color_manual(values = cols) +
  
  theme_classic(base_size = 16) +
  labs(
    title = "Classical Monocyte Proportions per Condition",
    y = "Proportion",
    x = ""
  ) +
  theme(
    legend.position = "none",
    # Bigger axis titles
     # Bigger axis tick labels
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    plot.title = element_blank()

  )+ 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
print(glm_results_dotplot)            
              

# Save figure
ggsave(plot = glm_results_dotplot,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "DCATS_classical_monocytes_dotplot.png")),
       width = 8, height = 6)








## Doing the same but simplifying clusters of monocytes, CD4 effector, CD8 effector-NK and B cells

# Duplicating original label first
sc10x$cluster_simplified <- as.character(sc10x$sgx_clusters)


# Collapse all monocytes variants automatically
sc10x$cluster_simplified <- ifelse(
  str_detect(sc10x$cluster_simplified, regex("doublet|triplet", ignore_case = TRUE)),
  "Multiplet",
  sc10x$cluster_simplified
)

# Collapse all monocytes variants automatically
sc10x$cluster_simplified <- ifelse(
  str_detect(sc10x$cluster_simplified, regex("monocyte", ignore_case = TRUE)),
  "Monocyte",
  sc10x$cluster_simplified
)

# Collapse all unconventional T cells cell variants automatically
sc10x$cluster_simplified <- ifelse(
  str_detect(sc10x$cluster_simplified, regex("γδ T|NKT|DN T cell|MAIT", ignore_case = TRUE)),
  "Unconventional T cell",
  sc10x$cluster_simplified
)

# Collapse all B cell variants automatically
sc10x$cluster_simplified <- ifelse(
  str_detect(sc10x$cluster_simplified, regex("B cell", ignore_case = TRUE)),
  "B cell",
  sc10x$cluster_simplified
)

# Collapse all NK cell variants automatically
sc10x$cluster_simplified <- ifelse(
  str_detect(sc10x$cluster_simplified, regex("NK", ignore_case = TRUE)),
  "NK",
  sc10x$cluster_simplified
)

# Collapse all effector CD4 T cells cell variants automatically
sc10x$cluster_simplified <- ifelse(
  str_detect(sc10x$cluster_simplified, regex("Effector CD4", ignore_case = TRUE)),
  "Effector CD4 T cell",
  sc10x$cluster_simplified
)

# Collapse all memory CD4 T cells cell variants automatically
sc10x$cluster_simplified <- ifelse(
  str_detect(sc10x$cluster_simplified, regex("Memory CD4 T cell", ignore_case = TRUE)),
  "Effector CD4 T cell",
  sc10x$cluster_simplified
)

# Collapse all memory CD4 T cells cell variants automatically
sc10x$cluster_simplified <- ifelse(
  str_detect(sc10x$cluster_simplified, regex("CD4 T cell", ignore_case = TRUE))
  & 
  !str_detect(sc10x$cluster_simplified, regex("Naive CD4 T cell", ignore_case = TRUE)),
  "Effector CD4 T cell",
  sc10x$cluster_simplified
)

# Collapse all effector CD8 T cells cell variants automatically
sc10x$cluster_simplified <- ifelse(
  (str_detect(sc10x$cluster_simplified, regex("CD8 T cell", ignore_case = TRUE)) |
  str_detect(sc10x$cluster_simplified, regex("CD8 cytotoxic T cell", ignore_case = TRUE))
  )& 
    !str_detect(sc10x$cluster_simplified, regex("Naive CD8 T cell", ignore_case = TRUE)),
  "Effector CD8 T cell",
  sc10x$cluster_simplified
)



# Convert back to factor if needed
sc10x$cluster_simplified <- factor(sc10x$cluster_simplified)



# ============================================================
# Create cell-type counts and proportions for each individual
# ============================================================

# Count cells for each individual and simplified cell type
celltype_counts <- as.data.frame(
  table(
    orig.ident = sc10x$orig.ident,
    cell_type = sc10x$cluster_simplified
  )
)

# Rename count column
colnames(celltype_counts)[3] <- "Cell_Count"

# Calculate total number of cells per individual
total_cells <- aggregate(
  Cell_Count ~ orig.ident,
  data = celltype_counts,
  FUN = sum
)

colnames(total_cells)[2] <- "Total_Cells"

# Merge total cell numbers back into the cell-type table
celltype_summary <- merge(
  celltype_counts,
  total_cells,
  by = "orig.ident"
)

# Calculate proportion of each cell type within each individual
celltype_summary$Proportion <- (
  celltype_summary$Cell_Count /
    celltype_summary$Total_Cells
)

# Optional: proportion as percentage
celltype_summary$Percentage <- (
  celltype_summary$Proportion * 100
)

# Rename columns for clarity
colnames(celltype_summary) <- c(
  "Individual",
  "Cell_Type",
  "Cell_Count",
  "Total_Cells",
  "Proportion",
  "Percentage"
)

# Sort by individual and cell type
celltype_summary <- celltype_summary[
  order(celltype_summary$Individual,
        celltype_summary$Cell_Type),
]

# Inspect
print(celltype_summary)

# ============================================================
# Save CSV
# ============================================================

write.csv(
  celltype_summary,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "celltype_counts_proportions.csv"
    )
  ),
  row.names = FALSE
)













# Compute KNN similarity matrix
knn_mat <- knn_simMat(sc10x@graphs$RNA_snn, sc10x$cluster_simplified)
print(knn_mat)

# Construct a count matrix with the simplified clusters
count_mat <- table(sc10x$orig.ident, sc10x$cluster_simplified)
print(count_mat)

# Colnames must be in the same order
knn_mat <- knn_mat[
    colnames(count_mat),
    colnames(count_mat)
]
identical(colnames(count_mat), rownames(knn_mat))
identical(colnames(count_mat), colnames(knn_mat))

# Run DCATS differential composition analysis
glm_results <- dcats_GLM(count_mat, design_mat, similarity_mat = knn_mat)
print(glm_results)

# Convert the list of lists into a single data frame
glm_results_df <- data.frame(
  Cell_Type = rownames(glm_results$coeffs),  # Keep only one column for cell types
  Coefficients = glm_results$coeffs[, 1],
  Coeff_Error = glm_results$coeffs_err[, 1],
  LR_Values = glm_results$LR_vals[, 1],
  LRT_pvals = glm_results$LRT_pvals[, 1],
  FDR = glm_results$fdr[, 1]
)

# Create the output directory if it doesn't exist
dir.create(file.path(PATH_ANALYSIS_OUTPUT), showWarnings = FALSE)

# Save results to a CSV file
write.csv(glm_results_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glm_results_cluster_simplified_batch.csv")), row.names = FALSE)




# to visualise the data with a boxplot (and significance bar)
# Extracting the number of monocytes for each individual
mono_counts <- count_mat[, "Monocyte"]

# Retrieving the condition column in the design_mat
conditions <- design_mat$condition[match(rownames(count_mat), rownames(design_mat))]

# Build a data frame for plotting
glm_results_plot <- data.frame(
  Sample = rownames(count_mat),
  Condition = conditions,
  Count = mono_counts
)

# Optional: compute proportion if you want
glm_results_plot$Proportion <- glm_results_plot$Count / rowSums(count_mat)

# Get DCATS FDR for annotation
fdr_val <- glm_results_df$FDR[glm_results_df$Cell_Type == "Monocyte"]

# Custom colors
cols <- c("Chronic" = "orange", "Symptomatic" = "red")

library(ggbeeswarm)
# Dot plot per sample with group mean + FDR annotation
glm_results_dotplot <- ggplot(glm_results_plot, aes(x = Condition, y = Proportion, color = Condition)) +
  geom_quasirandom(width = 0.1, size = 5, alpha = 0.9) + # size to increase the dot size
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.8, # size for the thickness of the mean bar
               color = "black", fatten = 0, alpha = 0.8) +
     # --- Significance bar with DCATS p-value ---
  geom_signif(
    comparisons = list(c("Chronic", "Symptomatic")),
    annotations = paste0("p = ", formatC(fdr_val, format = "e", digits = 1)),
    y_position = max(glm_results_plot$Proportion) * 1.08,
    tip_length = 0.01,
    textsize = 7, # size of the p-value
    color = "black"     # <- makes the bar black
  ) +
  
  # Colors
  scale_color_manual(values = cols) +
  
  theme_classic(base_family = "Arial", base_size = 16) + # increase overall base text
  labs(
    title = "Monocyte Proportions per Condition",
    y = "Proportion",
    x = ""
  ) +
  theme(
    legend.position = "none",
    # Bigger axis titles
     # Bigger axis tick labels
    axis.text.y = element_text(size = 18, family = "Arial"),
    axis.text.x = element_text(size = 18, family = "Arial"),
    axis.title.y = element_text(size = 20, family = "Arial"),
    axis.title.x = element_text(size = 20, family = "Arial"),
    plot.title = element_blank()

  )+ 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

print(glm_results_dotplot)            
              

# Save figure
ggsave(plot = glm_results_dotplot,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "DCATS_monocytes_dotplot.pdf")),
       width = 8, height = 6, device = cairo_pdf , dpi = 600)


# Making histogram for each of the major cell types

# Total cells per individual (ALL cells)
total_cells <- sc10x@meta.data %>%
  group_by(orig.ident) %>%
  summarise(total_n = n(), .groups = "drop")

# Counts of selected cell types
celltype_counts <- sc10x@meta.data %>%
  group_by(orig.ident, sgx_clusters) %>%
  summarise(n = n(), .groups = "drop")

# Merge and calculate proportion relative to TOTAL cells
celltype_prop <- celltype_counts %>%
  left_join(total_cells, by = "orig.ident") %>%
  mutate(proportion = n / total_n)

# Saving proportions
# Save to CSV file
write.csv(as.data.frame(celltype_prop), file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "cell_type_proportions.csv")), row.names = FALSE)



# Monocyte
# to change sample order
# to order individuals for the plotting  
  celltype_prop$orig.ident <- factor(
  celltype_prop$orig.ident,
  levels = sample_order
)

# Changing identifyers
id_map <- c(
  "16_scRNASeq_A83_AS_tp2" = "C1",
  "04_scRNASeq_A81_AS_tp2" = "C2",
  "13_scRNASeq_A62_AS_tp2" = "C3",
  "14_scRNASeq_A75_AS_tp2" = "C4",
  "15_scRNASeq_A76_S_tp2" = "S1",
  "12_scRNASeq_A86_S_tp2" = "S2",
    "08_scRNASeq_A82_S_tp2" = "S3") 
    


hist_mono <- ggplot(filter(celltype_prop, cluster_simplified == "Monocyte"),
                 aes(x = orig.ident, y = proportion)) +
  geom_bar(stat = "identity", fill = "darkgreen", width = 0.5) +
  scale_x_discrete(labels = id_map) +
  labs(title = NULL,
       y = "Proportion of total cells",
       x = "Individual") +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
  axis.title = element_text(size = 16, family = "Arial"),
  axis.text.y = element_text(size = 16, family = "Arial"))  # <- increase y-axis numbers size)
# Save figure
ggsave(plot = hist_mono,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "histogram_monocyte.pdf")),
       width = 6, height = 4,  device = cairo_pdf, dpi = 600)
       
# CD4
hist_cd4 <- ggplot(filter(celltype_prop, cluster_simplified == "Effector CD4 T cell"),
                aes(x = orig.ident, y = proportion)) +
  geom_bar(stat = "identity", fill = "firebrick", width = 0.5) +
  scale_x_discrete(labels = id_map) +
  labs(title = NULL,
       y = "Proportion of total cells",
       x = "Individual") +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
  axis.title = element_text(size = 16, family = "Arial"),
  axis.text.y = element_text(size = 16, family = "Arial"))  # <- increase y-axis numbers size)
# Save figure
ggsave(plot = hist_cd4,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "histogram_effector_CD4.pdf")),
       width = 6, height = 4, device = cairo_pdf, dpi = 600)
       
# B cell
hist_b <- ggplot(filter(celltype_prop, cluster_simplified == "B cell"),
              aes(x = orig.ident, y = proportion)) +
  geom_bar(stat = "identity", fill = "steelblue", width = 0.5) +
scale_x_discrete(labels = id_map) +
  labs(title = NULL,
       y = "Proportion of total cells",
       x = "Individual") +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
  axis.title = element_text(size = 16, family = "Arial"),
  axis.text.y = element_text(size = 16, family = "Arial"))  # <- increase y-axis numbers size)
# Save figure
ggsave(plot = hist_b,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "histogram_B_cell.pdf")),
       width = 6, height = 4, device = cairo_pdf, dpi = 600)
       

# For CD8 effector and NK cells

cd8_nk_prop <- celltype_prop %>%
  filter(cluster_simplified %in% c(
    "Effector CD8 T cell",
    "NK"
  ))

hist_cd8_nk = ggplot(cd8_nk_prop,
       aes(x = orig.ident,
           y = proportion,
           fill = cluster_simplified)) +
  geom_bar(stat = "identity", width =0.5) +
scale_x_discrete(labels = id_map) +
  labs(title = NULL,
       y = "Proportion of total cells",
       x = "Individual") +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
  legend.position="none",
  axis.title = element_text(size = 16, family = "Arial"),
  axis.text.y = element_text(size = 16, family = "Arial"))  # <- increase y-axis numbers size)
  

# Save figure
ggsave(plot = hist_cd8_nk,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "histogram_effector_CD8_NK_no_legend.pdf")),
       width = 6, height = 4, device = cairo_pdf, dpi = 600)

# Same with legend

cd8_nk_prop <- celltype_prop %>%
  filter(cluster_simplified %in% c(
    "Effector CD8 T cell",
    "NK"
  ))

hist_cd8_nk = ggplot(cd8_nk_prop,
       aes(x = orig.ident,
           y = proportion,
           fill = cluster_simplified)) +
  geom_bar(stat = "identity", width =0.5) +
  scale_x_discrete(labels = id_map) +
  labs(title = NULL,
       y = "Proportion of total cells",
       x = "Individual") +
  theme_classic(base_family = "Arial") +
  theme(axis.text.x = element_text(size = 16, family = "Arial"),
  axis.title = element_text(size = 16, family = "Arial"),
  axis.text.y = element_text(size = 16, family = "Arial"))  # <- increase y-axis numbers size)
# Save figure
ggsave(plot = hist_cd8_nk,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "histogram_effector_CD8_NK.pdf")),
       width = 6, height = 4, device = cairo_pdf, dpi = 600)
# saveRDS(PATH_SEURAT_RDS_FILE)




                                    ###                                     ###
                                        # Focusing on monocytes doublets #
                                    ###                                     ###
                                    
                                    
# 1) Multiplets



## Plotting doublets based on SciGeneX clustering and identification

sc10x$sgx_mono_multiplet <- as.character(sc10x$sgx_clusters)

# Create sgx_mono_multiplet column
sc10x$sgx_mono_multiplet <- ifelse(
  grepl("monocyte", sc10x$sgx_mono_multiplet, ignore.case = TRUE) &
    grepl("doublet|triplet", sc10x$sgx_mono_multiplet, ignore.case = TRUE),
  "multiplet",
  sc10x$sgx_mono_multiplet
)

## Plotting doublets based on SciGeneX clustering and identification


# Convert to factor for consistency


sgx_mono_multiplet_umap <- DimPlot(sc10x, group.by = "sgx_mono_multiplet", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "orange", "multiplet" = "purple")) +
  ggtitle("Multiplet prediction based on SciGeneX signatures")
 

#save plot 
ggsave(plot = sgx_mono_multiplet_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_sgx_mono_multiplet_prediction_umap.png"))
       , width = 23, height = 13, units ="cm")




## same but split by condition (chronic vs sympto)
sgx_mono_multiplet_umap_split <- DimPlot(sc10x, group.by = "sgx_mono_multiplet", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "orange", "multiplet" = "purple")) +
  ggtitle("Multiplet prediction based on SciGeneX signatures")
  
#save plot 
ggsave(plot = sgx_mono_multiplet_umap_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_sgx_mono_multiplet_prediction_umap_split.png"))
       , width = 30, height = 13, units ="cm")


# Same without legend and title

sgx_mono_multiplet_umap_split <- DimPlot(sc10x, group.by = "sgx_mono_multiplet", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "orange", "multiplet" = "purple")) +
  theme( plot.title = element_blank(), legend.position = "none")
  
#save plot 
ggsave(plot = sgx_mono_multiplet_umap_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_sgx_mono_multiplet_prediction_umap_split_no_legend.png"))
       , width = 30, height = 13, units ="cm")



######### To compare multiplet type proportion between conditions (DCATS tool)


# Construct a count matrix
count_mat <- table(sc10x$orig.ident, sc10x$sgx_mono_multiplet)
print(count_mat)


# Create design matrix for DCATS
# Reorder condition_info so it matches rownames(count_mat)
condition_vec <- condition_info[rownames(count_mat)]

# Build the design matrix taking into account the batch effect
design_mat <- data.frame(
   condition = condition_info[rownames(count_mat)],
   batch = batch_info[rownames(count_mat)]
)
print(design_mat)


# Creating the similarity matrix
knn_mat <- knn_simMat(sc10x@graphs$RNA_snn, sc10x$sgx_mono_multiplet)

# Run DCATS GLM
glm_results <- dcats_GLM(
  count_mat,
  design_mat,
  similarity_mat = knn_mat
)
print(glm_results)

# Convert the list of lists into a single data frame
glm_results_df <- data.frame(
  Cell_Type = rownames(glm_results$coeffs),
  Coefficient = glm_results$coeffs[, 1],
  SE = glm_results$coeffs_err[, 1],
  LRT = glm_results$LR_vals[, 1],
  pval = glm_results$LRT_pvals[, 1],
  FDR = glm_results$fdr[, 1]
)

print(glm_results_df)


# Save to CSV file
write.csv(glm_results_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glm_res_multiplet_singlet_batch.csv")), row.names = FALSE)





# Visualization

# Extracting the number of monocytes for each individual
sig_doublet_counts <- count_mat[, "multiplet"]

# Retrieving the condition column in the design_mat
conditions <- design_mat$condition[match(rownames(count_mat), rownames(design_mat))]

# Build a data frame for plotting
glm_results_plot <- data.frame(
  Sample = rownames(count_mat),
  Condition = conditions,
  Count = sig_doublet_counts
)

# Optional: compute proportion if you want
glm_results_plot$Proportion <- glm_results_plot$Count / rowSums(count_mat)

# Get DCATS FDR for annotation
fdr_val <- glm_results_df$FDR[glm_results_df$Cell_Type == "multiplet"]

# Custom colors
cols <- c("Chronic" = "orange", "Symptomatic" = "red")

# Dot plot per sample with group mean + FDR annotation
glm_results_dotplot <- ggplot(glm_results_plot, aes(x = Condition, y = Proportion, color = Condition)) +
  geom_quasirandom(width = 0.1, size = 5, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, linewidth = 0.8, # size for the thickness of the mean bar
               color = "black", fatten = 0, alpha = 0.8) +
     # --- Significance bar with DCATS p-value ---
  geom_signif(
    comparisons = list(c("Chronic", "Symptomatic")),
    annotations = paste0("p = ", formatC(fdr_val, format = "e", digits = 1)),
    y_position = max(glm_results_plot$Proportion) * 1.08,
    tip_length = 0.01,
    textsize = 7,
    color = "black"     # <- makes the bar black
  ) +
  
  # Colors
  scale_color_manual(values = cols) +
  
  theme_classic(base_size = 16) +
  labs(
    title = "Multiplet proportions per condition",
    y = "Proportion",
    x = ""
  ) +
  theme(
    legend.position = "none",
     # Bigger axis titles
     # Bigger axis tick labels
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    plot.title = element_blank()

  )+ 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
print(glm_results_dotplot)            
              

# Save figure
ggsave(plot = glm_results_dotplot,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "DCATS_multiplet_dotplot.png")),
       width = 8, height = 6)



# 2) Multiplet types: singlet, doublet, triplet






















