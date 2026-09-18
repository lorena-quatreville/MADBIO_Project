
BP_name = "mitochondrion"
score_name = paste0(BP_name, "_Score")
features_name = paste0(score_name, "1")

# Create the output directory
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME, BP_name))

# Loading seuratObject
## If needed to load the merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x) <- "cluster_renames"

# Subset to only include cells from Individual A and Individual B
sc10x_subset <- subset(
  sc10x,
  subset = multiplet_type == "Classical monocyte" &
           orig.ident %in% c(sample_id_A, sample_id_B)
)


sc10x_subset_down <- sc10x_subset
sc10x_subset_up <- sc10x_subset

# Loading the csv file containing DEG
DEG_file_path <- file.path(PATH_EXPERIMENT_OUTPUT, "04_Diff_gene_exp_analysis", "07_sympto_vs_chronic_Cmono", "01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_07_sympto_vs_chronic_Cmono_deg_res.csv")


DEG_file <- read.csv(DEG_file_path, header = TRUE)
downregulated_in_symptomatic <- DEG_file[
  DEG_file$p_val_adj < 0.05 &
  DEG_file$avg_log2FC < 0,
]
downregulated_gene_list <- downregulated_in_symptomatic[[1]]

# Retrieving the list of gene from SciGeneX

scigenex_module_list <- read.csv(file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_ModuleGenes_list.csv")))

# Intersection
# Create a list of genes per module
module_gene_lists <- lapply(
  scigenex_module_list[[2]],   # second column contains comma-separated genes
  function(x) {
    trimws(unlist(strsplit(x, ",")))
  }
)

# Name the list using module IDs from first column
names(module_gene_lists) <- scigenex_module_list[[1]]

# Intersect each module with downregulated genes
module_intersections <- lapply(
  module_gene_lists,
  function(gene_list) {
    intersect(gene_list, downregulated_gene_list)
  }
)

# View overlap for each module
module_intersections





                                                        ####### Same for up-regulated genes ######



upregulated_in_symptomatic <- DEG_file[
  DEG_file$p_val_adj < 0.05 &
  DEG_file$avg_log2FC > 0,
]
upregulated_gene_list <- upregulated_in_symptomatic[[1]]


# Retrieving the list of gene from SciGeneX

scigenex_module_list <- read.csv(file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_ModuleGenes_list.csv")))

# Intersection
# Create a list of genes per module
module_gene_lists <- lapply(
  scigenex_module_list[[2]],   # second column contains comma-separated genes
  function(x) {
    trimws(unlist(strsplit(x, ",")))
  }
)

# Name the list using module IDs from first column
names(module_gene_lists) <- scigenex_module_list[[1]]

# Intersect each module with downregulated genes
module_intersections <- lapply(
  module_gene_lists,
  function(gene_list) {
    intersect(gene_list, upregulated_gene_list)
  }
)

# View overlap for each module
module_intersections

# Extract mitochondrial module genes
mt_genes <- trimws(
    unlist(strsplit(scigenex_module_list[[2]][12], ","))
)

# Keep only genes present in upregulated DEGs
final_filtered_genes <- intersect(
    mt_genes,
    upregulated_gene_list
)

final_filtered_genes


# Saving the list in a csv file
write.csv(data.frame(gene = final_filtered_genes),
          file = file.path(
        PATH_ANALYSIS_OUTPUT, BP_name,
        paste0(outputFilesPrefix, BP_name, "gene_signature_up.csv")), row.names = FALSE)

# Compute module scores
sc10x_subset_up <- AddModuleScore(sc10x_subset_up, features = list(final_filtered_genes), name = score_name)

Seurat::DoHeatmap(
  sc10x_subset_down,
  features = final_filtered_genes,
  group.by = "cluster_renames"
)

# to change sample order 
sample_order <- c("08_scRNASeq_A82_S_tp2", "12_scRNASeq_A86_S_tp2", "15_scRNASeq_A76_S_tp2", "04_scRNASeq_A81_AS_tp2", "16_scRNASeq_A83_AS_tp2", "13_scRNASeq_A62_AS_tp2", "14_scRNASeq_A75_AS_tp2" ) 
# to order individuals for the plotting 
sc10x_subset_up$orig.ident <- factor(
  sc10x_subset_up$orig.ident,
  levels = sample_order
)


# Global heatmap of the genes contained in the signature

heatmap_plot <- DoHeatmap(sc10x_subset_up, features = final_filtered_genes, group.by = "orig.ident") +
  ggtitle(paste0("Heatmap of Genes contained in ", BP_name, " signature")) +
  theme(axis.text.y = element_text(size = 8))  # smaller text for readability

print(heatmap_plot)

# Save to file
ggsave(
  plot = heatmap_plot,
  filename = file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, "_", BP_name, "_signature_Heatmap_up_reg_sympto.png")),
  width = 10, height = 8, dpi = 300
)



## Displaying the UMAP with the scores

# Plot UMAP FeaturePlot for the scores
umap_score <- FeaturePlot(
  sc10x_subset_up,
  features = features_name
) + 
  scale_color_gradientn(
    colors = c("green", "blue", "purple", "red", "yellow")
  ) +
  ggtitle(score_name)


# Display plots
print(umap_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,   paste0( outputFilesPrefix, "umap_", score_name, "_up.png")), plot = umap_score, width = 12, height =8, dpi = 300)



# Visualisation using violin plots and violin boxplot


# Visualisation of the score accross clusters:

# violin plot:
violin_plot <- VlnPlot(sc10x_subset_up, features = features_name, group.by = "condition") +
  ggtitle(paste0(score_name, " by Condition")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,  paste0( outputFilesPrefix, "vlnPlot_", score_name, "_up.png")), plot = violin_plot, width = 12, height =8, dpi = 300)
print(violin_plot)


# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x_subset_up, features = features_name, group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle(paste0(score_name, " by Condition"))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,   paste0( outputFilesPrefix, "vlnBoxPlot_", score_name, "_up.png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)
print(violin_box_plot)



################## statistical test to compare means (Welsh Anova which assumes normal distribution but unequal variances) ############################

# Test
features_scores <- sc10x_subset_up@meta.data %>%
  filter(condition == "Symptomatic") %>%
  pull(.data[[features_name]])
  # Shapiro–Wilk is very sensitive; With large numbers of cells, it will almost always reject normality; A significant p-value does not automatically invalidate Welch ANOVA
#shapiro.test(features_scores) Not doable because observations shuold not go up to 5000

# Histogram
hist(features_scores, breaks = 50,
     main = "Histogram of the score accross cells",
     xlab = score_name)

# QQplot
qqnorm(features_scores)
qqline(features_scores, col = "red")
# Mild deviations in the tails → usually acceptable; Strong S-shape or heavy skew → non-normal

###### Test for equal variance
#leveneTest(
#  .data[[features_name]] ~ condition,
#  data = sc10x_subset@meta.data,
# center = median   # more robust than mean
#)
# p < 0.05 → variances are significantly different; p ≥ 0.05 → no evidence of unequal variance
# If not significant, we have to use: Wesh anova (correction for unequal variance) and Games-Howell post-hoc test (for unequal variance and sample size


# Visualisation of the score accross  specific clusters:

# Wilcoxon rank sum test instead

wilcox_test <- wilcox.test(
  formula = reformulate("condition", response = features_name),
  data = sc10x_subset_up@meta.data,
  exact = FALSE
)

print(wilcox_test)

p_val <- wilcox_test$p.value
p_label <- paste0("p = ", signif(p_val, 3))

wilcox_df <- data.frame(
  group1 = "Chronic",
  group2 = "Symptomatic",
  p.adj = p_val,
  p_label = p_label
)

wilcox_df <- wilcox_df %>%
  mutate(label = paste0("p = ", scales::pvalue(p.adj)))
  
max_y <- max(sc10x_subset_up@meta.data[[features_name]], na.rm = TRUE)
wilcox_df$y.position <- max_y + 0.05

# Define custom colors
condition_colors <- c(
  "Chronic" = "#E69F00",      # orange
  "Symptomatic" = "#D55E00"   # red
)

# Plot (PAPER)
violin_boxplot_specific_clusters <- ggplot(sc10x_subset_up@meta.data, 
                                           aes(x = condition, y = .data[[features_name]])) +
    geom_violin(aes(fill = condition, color = condition), trim = FALSE, alpha = 0.5, width = 0.8) +
    geom_boxplot(aes(fill = condition, color = condition), 
                 width = 0.2, position = position_dodge(1.5), outlier.shape = NA) +  
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    stat_pvalue_manual(wilcox_df, label = "p_label", y.position = "y.position", size = 12, step.increase = 0.1) +  
    theme_classic(base_size = 22) +   # removes grey background & grids
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 25),
    axis.text = element_text(size = 25), legend.position = "none") +
    labs(y = score_name, x = "Condition")

print(violin_boxplot_specific_clusters)

# Print and save the plot
print(violin_boxplot_specific_clusters)
ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT, BP_name, paste0(outputFilesPrefix, "vlnBoxPlot_", score_name, "_chronic_sympto_wilcox_up_reg_sympto.png")),
       plot = violin_boxplot_specific_clusters, width = 12, height = 8, dpi = 300)


   
       
