
BP_name = "migration"
score_name = paste0(BP_name, "_Score")
features_name = paste0(score_name, "1")

# Loading seuratObject
## If needed to load the merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x) <- "cluster_renames"

# Subset to only include cells from Individual A and Individual B
sc10x_subset <- subset(
  sc10x,
  subset = multiplet_type %in% c("Classical monocyte", "Non-classical monocyte") &
           orig.ident %in% c(sample_id_A, sample_id_B)
)

# Loading the csv file containing DEG
BP_down_file_path <- file.path(PATH_EXPERIMENT_OUTPUT, "04_Diff_gene_exp_analysis", "12_sympto_vs_chronic_all_mono", "01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_12_sympto_vs_chronic_all mono_deg_res.csv")


BP_down_file <- read.csv(BP_down_file_path, header = TRUE)
downregulated_in_symptomatic <- BP_down_file[
  BP_down_file$p_val_adj < 0.05 &
  BP_down_file$avg_log2FC < 0,
]
downregulated_gene_list <- downregulated_in_symptomatic[[1]]

# Retrieving the list of gene from SciGeneX
proliferation_df <- read.csv(
  file.path(PATH_ANALYSIS_OUTPUT,
            paste0(outputFilesPrefix, paste0("Gene_sign_", BP_name ,".csv"))),
  stringsAsFactors = FALSE
)

scigenex_genes <- proliferation_df$gene
assign(paste0("all_", BP_name, "_genes"), scigenex_genes)

# Intersection between down-regulated genes and genes found in modules
final_filtered_genes <- intersect(
  get(paste0("all_", BP_name, "_genes")),
  downregulated_gene_list
)


# Compute module scores
sc10x_subset <- AddModuleScore(sc10x_subset, features = list(final_filtered_genes), name = score_name)


## Displaying the UMAP with the scores

# Plot UMAP FeaturePlot for the scores
umap_score <- FeaturePlot(
  sc10x_subset,
  features = features_name
) + 
  scale_color_gradientn(
    colors = c("green", "blue", "purple", "red", "yellow")
  ) +
  ggtitle(score_name)


# Display plots
print(umap_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT,  paste0( outputFilesPrefix, "umap_", score_name, ".png")), plot = umap_score, width = 12, height =8, dpi = 300)



# Visualisation using violin plots and violin boxplot


# Visualisation of the score accross clusters:

# violin plot:
violin_plot <- VlnPlot(sc10x_subset, features = features_name, group.by = "condition") +
  ggtitle(paste0(score_name, " by Condition")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "vlnPlot_", score_name, ".png")), plot = violin_plot, width = 12, height =8, dpi = 300)
print(violin_plot)


# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x_subset, features = features_name, group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle(paste0(score_name, " by Condition"))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT,  paste0( outputFilesPrefix, "vlnBoxPlot_", score_name, ".png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)
print(violin_box_plot)



################## statistical test to compare means (Welsh Anova which assumes normal distribution but unequal variances) ############################

# Test
features_scores <- sc10x_subset@meta.data %>%
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


# Welch's ANOVA
welch_anova <- oneway.test(
  formula = reformulate("condition", response = features_name),
  data = sc10x_subset@meta.data,
  var.equal = FALSE
)
print(welch_anova)


# Ensure custom_group exists
sc10x_subset@meta.data$condition <- as.factor(sc10x_subset@meta.data$condition)
    
    
# Perform pairwise Games-Howell test # Handles unequal sample size and variances)
games_howell_results <- gamesHowellTest(
  formula = reformulate("condition", response = features_name),
  data = sc10x_subset@meta.data
)


# Convert p-value matrix to tidy data frame
pvals <- games_howell_results$p.value

games_howell_df <- as.data.frame(as.table(pvals))
colnames(games_howell_df) <- c("group1", "group2", "p.adj")


# Convert to character before filtering NA
games_howell_df <- games_howell_df %>%
  mutate(group1 = as.character(group1),
         group2 = as.character(group2))%>%
  filter(!is.na(p.adj))


# Add significance label
games_howell_df <- games_howell_df %>%
  mutate(p.adj.signif = p_format(p.adj))

# Remove comparisons that failed or produce NA
games_howell_df <- games_howell_df %>% 
    filter(!is.na(p.adj), !is.na(group1), !is.na(group2))

# Add the y.position column
max_y <- max(sc10x_subset@meta.data[[features_name]], na.rm = TRUE)


games_howell_df <- games_howell_df %>%
  mutate(p.adj.signif = case_when(
    p.adj <= 0.001 ~ "***",   # p-value <= 0.001
    p.adj <= 0.01  ~ "**",    # p-value <= 0.01
    p.adj <= 0.05  ~ "*",     # p-value <= 0.05
    p.adj > 0.05   ~ "ns"     # p-value > 0.05 (not significant)
  ))
  
games_howell_df$y.position <- max_y + seq(0.05, by = 0.05, length.out = nrow(games_howell_df))

# Define custom colors
condition_colors <- c(
  "Chronic" = "#E69F00",      # orange
  "Symptomatic" = "#D55E00"   # red
)

# Reorder condition factor
sc10x_subset@meta.data$condition <- factor(
  sc10x_subset@meta.data$condition,
  levels = c("Chronic","Symptomatic")
)

# Plot
violin_boxplot_specific_clusters <- ggplot(sc10x_subset@meta.data, 
                                           aes(x = condition, y = .data[[features_name]])) +
    geom_violin(aes(fill = condition, color = condition), trim = FALSE, alpha = 0.5, width = 0.8) +
    geom_boxplot(aes(fill = condition, color = condition), 
                 width = 0.2, position = position_dodge(1.5), outlier.shape = NA) +  
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    stat_pvalue_manual(games_howell_df, label = "p.adj.signif", y.position = "y.position", size = 12, step.increase = 0.1) +  
    theme_classic(base_size = 22) +   # removes grey background & grids
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 25),
    axis.text = element_text(size = 25), legend.position = "none") +
    labs(y = score_name, x = "Condition")

print(violin_boxplot_specific_clusters)

# Print and save the plot
print(violin_boxplot_specific_clusters)
ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT,  paste0(outputFilesPrefix, "vlnBoxPlot_", score_name, "_chronic_sympto_all_mono.png")),
       plot = violin_boxplot_specific_clusters, width = 12, height = 8, dpi = 300)
       




# Separating classical monocytes from non-classical monocytes



sc10x_subtype = subset(sc10x_subset, (multiplet_type == "Classical monocyte"))
subtype = "CM"
sc10x_subtype = subset(sc10x_subset, (multiplet_type == "Non-classical monocyte"))
subtype = "NCM"


# Test
features_scores <- sc10x_subtype@meta.data %>%
  filter(condition == "Chronic") %>%
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
library(car)
#leveneTest(
#  .data[[features_name]] ~ condition,
#  data = sc10x_subtype@meta.data,
#  center = median   # more robust than mean
#)
# p < 0.05 → variances are significantly different; p ≥ 0.05 → no evidence of unequal variance
# If not significant, we have to use: Wesh anova (correction for unequal variance) and Games-Howell post-hoc test (for unequal variance and sample size


# Visualisation of the score accross  specific clusters:


# Welch's ANOVA
welch_anova <- oneway.test(
  formula = reformulate("condition", response = features_name),
  data = sc10x_subtype@meta.data,
  var.equal = FALSE
)
print(welch_anova)


# Ensure custom_group exists
sc10x_subtype@meta.data$condition <- as.factor(sc10x_subtype@meta.data$condition)
    
    
# Perform pairwise Games-Howell test # Handles unequal sample size and variances)
games_howell_results <- gamesHowellTest(
  formula = reformulate("condition", response = features_name),
  data = sc10x_subtype@meta.data
)


# Convert p-value matrix to tidy data frame
pvals <- games_howell_results$p.value

games_howell_df <- as.data.frame(as.table(pvals))
colnames(games_howell_df) <- c("group1", "group2", "p.adj")


# Convert to character before filtering NA
games_howell_df <- games_howell_df %>%
  mutate(group1 = as.character(group1),
         group2 = as.character(group2))%>%
  filter(!is.na(p.adj))


# Add significance label
games_howell_df <- games_howell_df %>%
  mutate(p.adj.signif = p_format(p.adj))

# Remove comparisons that failed or produce NA
games_howell_df <- games_howell_df %>% 
    filter(!is.na(p.adj), !is.na(group1), !is.na(group2))

# Add the y.position column
max_y <- max(sc10x_subtype@meta.data[[features_name]], na.rm = TRUE)


games_howell_df <- games_howell_df %>%
  mutate(p.adj.signif = case_when(
    p.adj <= 0.001 ~ "***",   # p-value <= 0.001
    p.adj <= 0.01  ~ "**",    # p-value <= 0.01
    p.adj <= 0.05  ~ "*",     # p-value <= 0.05
    p.adj > 0.05   ~ "ns"     # p-value > 0.05 (not significant)
  ))
  
games_howell_df$y.position <- max_y + seq(0.05, by = 0.05, length.out = nrow(games_howell_df))

# Define custom colors
condition_colors <- c(
  "Chronic" = "#E69F00",      # orange
  "Symptomatic" = "#D55E00"   # red
)

# Plot
violin_boxplot_specific_clusters <- ggplot(sc10x_subtype@meta.data, 
                                           aes(x = condition, y = .data[[features_name]])) +
    geom_violin(aes(fill = condition, color = condition), trim = FALSE, alpha = 0.5, width = 0.8) +
    geom_boxplot(aes(fill = condition, color = condition), 
                 width = 0.2, position = position_dodge(1.5), outlier.shape = NA) +  
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    stat_pvalue_manual(games_howell_df, label = "p.adj.signif", y.position = "y.position", size = 12, step.increase = 0.1) +  
    theme_classic(base_size = 22) +   # removes grey background & grids
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 25),
    axis.text = element_text(size = 25), legend.position = "none") +
    labs(y = score_name, x = "Condition")

print(violin_boxplot_specific_clusters)

# Print and save the plot
print(violin_boxplot_specific_clusters)
ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT,  paste0(outputFilesPrefix, "vlnBoxPlot_", score_name, "_", subtype, "_chronic_sympto.png")),
       plot = violin_boxplot_specific_clusters, width = 12, height = 8, dpi = 300)
       
       
       
       
# SciGeneX heatmap with only genes involved in the signature and their corresponding module

modules_list <- res_cluster_subtype@gene_clusters

# Identifying which modules contain genes from the signature
modules_with_signature <- modules_list[
  sapply(modules_list, function(m){
    any(m %in% final_filtered_genes)
  })
]

# Filtering genes inside each module (to remove the ones that are not part of the signature)
modules_filtered <- lapply(modules_with_signature, function(m){
  intersect(m, final_filtered_genes)
})

#  Creating a df of genes and modules
gene_module_df <- data.frame(
  module = names(modules_filtered),
  genes = sapply(modules_filtered, function(g) paste(g, collapse = "; "))
)

# Saving the gene signature list with the module correspondance
write.csv(gene_module_df,
          file = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, BP_name, subtype, "Gene_signature_modules.csv")), row.names = FALSE)

# Function definition for plotting a heatmap only for genes contained in the signature # NOT working
heatmap_signature_module <- function(module_list){

  p_heatmap <- plot_ggheatmap(
    res_cluster_subtype[module_list,
    ident = Seurat::Idents(sc10x_subtype)
  ) +
  theme(strip.text.y = element_text(size = 4))

  ggsave(
    plot = p_heatmap,
    filename = file.path(
      PATH_ANALYSIS_OUTPUT,
      paste0(BP_name, "_", subtype, "_signature_heatmap_modules.png")
    ),
    width = 23, height = 13, units = "cm"
  )

  plot_heatmap(
    module_list,
    cell_clusters = Seurat::Idents(sc10x_subtype)
  )
}

heatmap_signature_module(modules_with_signature)

Seurat::DoHeatmap(
  sc10x_subset,
  features = final_filtered_genes,
  group.by = "cluster_renames"
)

# to change sample order 
sample_order <- c( "08_scRNASeq_A82_S_tp2", "12_scRNASeq_A86_S_tp2", "15_scRNASeq_A76_S_tp2", "04_scRNASeq_A81_AS_tp2", "16_scRNASeq_A83_AS_tp2", "13_scRNASeq_A62_AS_tp2", "14_scRNASeq_A75_AS_tp2" ) 
# to order individuals for the plotting 
sc10x_subset$orig.ident <- factor(
  sc10x_subset$orig.ident,
  levels = sample_order
)


# Global heatmap of the genes contained in the signature

heatmap_plot <- DoHeatmap(sc10x_subset, features = final_filtered_genes, group.by = "orig.ident") +
  ggtitle(paste0("Heatmap of Genes contained in ", BP_name, " signature")) +
  theme(axis.text.y = element_text(size = 8))  # smaller text for readability

print(heatmap_plot)

# Save to file
ggsave(
  plot = heatmap_plot,
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_", BP_name, subtype, "_signature_Heatmap.png")),
  width = 10, height = 8, dpi = 300
)
