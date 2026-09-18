
# Variables
features_name <- "ICAM1"
score_name <- "ICAM1 expression"

## If needed to load the merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x) <- "cluster_renames"

# Subset to only include cells from Individual A and Individual B
sc10x_subset <- subset(
  sc10x,
  subset = multiplet_type == "Classical monocyte" &
           orig.ident %in% c(sample_id_A, sample_id_B)
)


plot = VlnPlot(sc10x_subset, features_name, group.by="condition")

ggsave(
filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "VlnPlot_ICAM1.png")),
  plot = plot,
  width = 20,
  height = 12
)
sc10x_subset@meta.data$ICAM1 <- FetchData(sc10x_subset, vars = features_name)[,1]
summary(sc10x_subset@meta.data$ICAM1)

# Statistical test
# Wilcoxon rank sum test instead

wilcox_test <- wilcox.test(
  formula = reformulate("condition", response = features_name),
  data = sc10x_subset@meta.data,
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
  
max_y <- max(sc10x_subset@meta.data[[features_name]], na.rm = TRUE)
wilcox_df$y.position <- max_y + 0.05




# Visualization
# Define custom colors
condition_colors <- c(
  "Chronic" = "#E69F00",      # orange
  "Symptomatic" = "#D55E00"   # red
)



# Use existing t.test result
wilcox_df <- data.frame(
  group1 = "Chronic",
  group2 = "Symptomatic",
  p = t_test_result$p.value
)

wilcox_df$p_label <- paste0("p_value= ", signif(wilcox_df$p, 3))
wilcox_df$y.position <- max(sc10x_subset@meta.data[[features_name]], na.rm = TRUE) * 1.08

# Plot (PAPER)
violin_boxplot_specific_clusters <- ggplot(sc10x_subset@meta.data, 
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
ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "vlnBoxPlot_", score_name, "_chronic_sympto_wilcox.png")),
       plot = violin_boxplot_specific_clusters, width = 12, height = 8, dpi = 300)
       
       
       
       
       


# T test if needed
t_test_result <- t.test(
    ICAM1 ~ condition,
    data = sc10x_subset@meta.data)

t_test_result

