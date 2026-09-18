# Migration signature score

# Parameter for this step:
BP_name = "migration"
score_name = paste0(BP_name, "_Score")
features_name = paste0(score_name, "1")

# Loading data
library(readr)  
go_results <- read_csv(file = file.path(PATH_ANALYSIS_OUTPUT,  paste0("Module_", module_nb), paste0(outputFilesPrefix, "GO_analysis_results_deg_BP.csv")))

# Extract relevant columns
description_col <- "Description"  # Adjust if the column name differs
genes_col <- names(go_results)[8] # Assuming the 8th column contains gene lists

# Filter terms with a specific biological process in the description
BP_gene <- go_results %>%
  filter(grepl(BP_name, .[[description_col]], ignore.case = TRUE)) %>%
  pull(genes_col)

# Flatten lists (in case genes are separated by semicolons or commas)
BP_gene_list <- unique(unlist(strsplit(paste(BP_gene , collapse = ";"), ";")))


# Display results
print(BP_name)
print(BP_gene_list)


# FUNCTION to clean and extract unique genes from a list of GO term genes
extract_unique_genes <- function(gene_list) {
  unique_genes <- unique(trimws(unlist(strsplit(paste(gene_list, collapse = "/"), "/"))))
  return(unique_genes[unique_genes != ""]) # Remove empty entries
}

# Extract unique heat response genes
BP_gene_list <- extract_unique_genes(BP_gene_list)

# Display results
print(BP_name)
print(BP_gene_list)


## Computing a heat response score and a unfolded response score:


# Ensure genes are present in Seurat object
available_genes <- rownames(sc10x@assays$RNA@data)  # Extract gene names from the RNA assay
BP_gene_list <- BP_gene_list[BP_gene_list %in% available_genes]


# Compute module scores
sc10x <- AddModuleScore(sc10x, features = list(BP_gene_list), name = score_name)


## Displaying the UMAP with the scores

# Plot UMAP FeaturePlot for the scores
umap_score <- FeaturePlot(
  sc10x,
  features = features_name
) + 
  scale_color_gradientn(
    colors = c("green", "blue", "purple", "red", "yellow")
  ) +
  ggtitle(score_name)


# Display plots
print(umap_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb),  paste0( outputFilesPrefix, "umap_", score_name, ".png")), plot = umap_score, width = 12, height =8, dpi = 300)





# Visualisation using violin plots and violin boxplot


# Visualisation of the score accross clusters:

# violin plot:
violin_plot <- VlnPlot(sc10x, features = features_name, group.by = "condition") +
  ggtitle(paste0(score_name, " by Cluster")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb),  paste0( outputFilesPrefix, "vlnPlot_", score_name, ".png")), plot = violin_plot, width = 12, height =8, dpi = 300)
print(violin_plot)


# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x, features = features_name, group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle(paste0(score_name, " by Cluster"))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb),  paste0( outputFilesPrefix, "vlnBoxPlot_", score_name, ".png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)
print(violin_box_plot)

# Visualisation of the score accross  specific clusters:


################## statistical test to compare means (Welsh Anova which assumes normal distribution but unequal variances) ############################

# Visualisation of the score accross  specific clusters:


# Subset data and create a temporary column for grouping
sc10x_subset <- subset(
  sc10x,
  multiplet_type == "Classical monocyte" & orig.ident %in% c(
    "08_scRNASeq_A82_S_tp2",
    "12_scRNASeq_A86_S_tp2",
    "15_scRNASeq_A76_S_tp2",
    "04_scRNASeq_A81_AS_tp2",
    "16_scRNASeq_A83_AS_tp2",
    "13_scRNASeq_A62_AS_tp2",
    "14_scRNASeq_A75_AS_tp2"
  )
)

# Create custom ordered grouping
sc10x_subset@meta.data <- sc10x_subset@meta.data %>%
  mutate(custom_group = case_when(
  orig.ident %in% c("04_scRNASeq_A81_AS_tp2",
                      "16_scRNASeq_A83_AS_tp2",
                      "13_scRNASeq_A62_AS_tp2",
                      "14_scRNASeq_A75_AS_tp2") ~ "CM – Chronic",
    orig.ident == "15_scRNASeq_A76_S_tp2" ~ "CM – A76",
    orig.ident == "12_scRNASeq_A86_S_tp2" ~ "CM – A86",
     orig.ident == "08_scRNASeq_A82_S_tp2" ~ "CM – A82",
    TRUE ~ NA_character_
  ))

# Welch's ANOVA
welch_anova <- oneway.test(
  formula = reformulate("custom_group", response = features_name),
  data = sc10x_subset@meta.data,
  var.equal = FALSE
)
print(welch_anova)


# Ensure custom_group exists
sc10x_subset@meta.data$custom_group <- as.factor(sc10x_subset@meta.data$custom_group)
    
    
# Perform pairwise Games-Howell test # Handles unequal sample size and variances)
games_howell_results <- gamesHowellTest(
  formula = reformulate("custom_group", response = features_name),
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


sc10x_subset@meta.data$custom_group <- factor(sc10x_subset@meta.data$custom_group, levels = c("CM – Chronic", "CM – A76", "CM – A86", "CM – A82"))


# Plot
violin_boxplot_specific_clusters <- ggplot(sc10x_subset@meta.data, 
                                           aes(x = custom_group, y = .data[[features_name]])) +
    geom_violin(aes(fill = custom_group, color = custom_group), trim = FALSE, alpha = 0.5, width = 0.8) +
    geom_boxplot(aes(fill = custom_group, color = custom_group), 
                 width = 0.2, position = position_dodge(1.5), outlier.shape = NA) +  
    stat_pvalue_manual(games_howell_df, label = "p.adj.signif", y.position = "y.position", step.increase = 0.1) +  
    ggtitle(paste0(score_name, " in Selected Clusters")) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "none") +
    labs(y = score_name, x = "Cluster Group")

print(violin_boxplot_specific_clusters)

# Print and save the plot
print(violin_boxplot_specific_clusters)
ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb),  paste0(outputFilesPrefix, "vlnBoxPlot_", score_name, "_chronic_sympto.png")),
       plot = violin_boxplot_specific_clusters, width = 12, height = 8, dpi = 300)
       

