

## To compare multiplet proportion between the two conditions

# Per-sample multiplet counts
multiplet_sample_df <- sc10x@meta.data %>%
  group_by(orig.ident, condition, sgx_multiplet) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = sgx_multiplet,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(
    total = singlet + multiplet,
    prop_multiplet = multiplet / total
  )

print(multiplet_sample_df)

# Binomial glm for comparing

glmm_multiplet <- glmer(
  cbind(multiplet, singlet) ~ condition + (1 | orig.ident),
  data = multiplet_sample_df,
  family = binomial()
)

summary(glmm_multiplet)

# Extract fixed effects summary
glmm_summary <- summary(glmm_multiplet)$coefficients

# Convert to dataframe
glmm_summary_df <- as.data.frame(glmm_summary)
glmm_summary_df$Term <- rownames(glmm_summary_df)
rownames(glmm_summary_df) <- NULL

print(glmm_summary_df)

# Save to CSV file
write.csv(glmm_summary_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glmm_multiplet_fixed_effects.csv")), row.names = FALSE)

# Saving random effects
ranef_df <- as.data.frame(ranef(glmm_multiplet)$orig.ident)
print(ranef_df)

# Save to CSV file
write.csv(ranef_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glmm_random_effects_results.csv")), row.names = FALSE)



###################### Visualization

# extracting glmm effect for annotation
glmm_pvalue <- summary(glmm_multiplet)$coefficients["conditionSymptomatic", "Pr(>|z|)"]
glmm_pvalue_label <- paste0("p = ", signif(glmm_pvalue, 3))

cols <- c("Chronic" = "orange", "Symptomatic" = "red")

# Create plotting dataframe (directly from multiplet_sample_df)
glm_results_plot <- multiplet_sample_df %>%
  dplyr::select(orig.ident, condition, prop_multiplet) %>%
  dplyr::rename(
    Sample = orig.ident,
    Condition = condition,
    Proportion = prop_multiplet
  )

# Build the dotplot with observed proportions
multiplet_dotplot <- ggplot(glm_results_plot,
                            aes(x = Condition, y = Proportion, color = Condition)) +
  
  geom_jitter(width = 0.1, size = 3, alpha = 0.9) +
  
  # Mean bar per group
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               color = "black", fatten = 0, alpha = 0.8) +
  
  # Significance bar with GLMM p-value
  geom_signif(
    comparisons = list(c("Chronic", "Symptomatic")),
    annotations = glmm_pvalue_label,
    y_position = max(glm_results_plot$Proportion) * 1.08,
    tip_length = 0.01,
    textsize = 4,
    color = "black"
  ) +
  
  scale_color_manual(values = cols) +
  
  theme_minimal(base_size = 12) +
  labs(
    title = "Multiplet Proportion per Condition",
    y = "Proportion of multiplets",
    x = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16)
  )

print(multiplet_dotplot)

# Save figure
ggsave(
  plot = multiplet_dotplot,
  filename = file.path(PATH_ANALYSIS_OUTPUT,
                       paste0(outputFilesPrefix, "GLMM_multiplet_proportions_dotplot.png")),
  width = 8, height = 6, dpi = 300
)









######### To compare multiplet type proportion between conditions (DCATS tool)


# Construct a count matrix
count_mat <- table(sc10x$orig.ident, sc10x$multiplet_type_simplified)
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
knn_mat <- knn_simMat(sc10x@graphs$RNA_snn, sc10x$multiplet_type_simplified)

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
write.csv(glm_results_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glm_results_batch.csv")), row.names = FALSE)







# Visualisation of Non-classical monocyte-T cell doublet proportion comparison between conditions

# Extracting the number of monocytes for each individual
sig_doublet_counts <- count_mat[, "Non-classical monocyte_T cell"]

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
fdr_val <- glm_results_df$FDR[glm_results_df$Cell_Type == "Non-classical monocyte_T cell"]

# Custom colors
cols <- c("Chronic" = "orange", "Symptomatic" = "red")

# Dot plot per sample with group mean + FDR annotation
glm_results_dotplot <- ggplot(glm_results_plot, aes(x = Condition, y = Proportion, color = Condition)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               color = "black", fatten = 0, alpha = 0.8) +
     # --- Significance bar with DCATS p-value ---
  geom_signif(
    comparisons = list(c("Chronic", "Symptomatic")),
    annotations = paste0("p = ", signif(fdr_val, 3)),
    y_position = max(glm_results_plot$Proportion) * 1.08,
    tip_length = 0.01,
    textsize = 4,
    color = "black"     # <- makes the bar black
  ) +
  
  # Colors
  scale_color_manual(values = cols) +
  
  theme_minimal(base_size = 12) +
  labs(
    title = "Non-classical monocyte-T cell complexe proportions per condition",
    y = "Proportion",
    x = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16)
  )

print(glm_results_dotplot)            
              

# Save figure
ggsave(plot = glm_results_dotplot,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "DCATS_NonClassical_monocyte_Tcell_dotplot.png")),
       width = 8, height = 6)






########################## To compare multiplet_number (singlet/doublet/triplet) proportions between conditions (DCATS tool)



# Construct a count matrix
count_mat <- table(sc10x$orig.ident, sc10x$multiplet_number)
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
knn_mat <- knn_simMat(sc10x@graphs$RNA_snn, sc10x$multiplet_number)

# Run DCATS GLM
glm_results <- dcats_GLM(
  count_mat,
  design_mat,
  similarity_mat = knn_mat
)
print(glm_results)

# Convert the list of lists into a single data frame
glm_results_multiplet_nbr_df <- data.frame(
  Cell_Type = rownames(glm_results$coeffs),
  Coefficient = glm_results$coeffs[, 1],
  SE = glm_results$coeffs_err[, 1],
  LRT = glm_results$LR_vals[, 1],
  pval = glm_results$LRT_pvals[, 1],
  FDR = glm_results$fdr[, 1]
)

print(glm_results_multiplet_nbr_df)


# Save to CSV file
write.csv(glm_results_multiplet_nbr_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glm_results_multiplet_nbr_batch.csv")), row.names = FALSE)







# Visualisation of triplet proportion comparison between conditions

# Extracting the number of monocytes for each individual
triplet_counts <- count_mat[, "triplet"]

# Retrieving the condition column in the design_mat
conditions <- design_mat$condition[match(rownames(count_mat), rownames(design_mat))]

# Build a data frame for plotting
glm_results_plot <- data.frame(
  Sample = rownames(count_mat),
  Condition = conditions,
  Count = triplet_counts
)

# Optional: compute proportion if you want
glm_results_plot$Proportion <- glm_results_plot$Count / rowSums(count_mat)

# Get DCATS FDR for annotation
fdr_val <- glm_results_multiplet_nbr_df$FDR[glm_results_multiplet_nbr_df$Cell_Type == "triplet"]

# Custom colors
cols <- c("Chronic" = "orange", "Symptomatic" = "red")

# Dot plot per sample with group mean + FDR annotation
glm_results_dotplot <- ggplot(glm_results_plot, aes(x = Condition, y = Proportion, color = Condition)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               color = "black", fatten = 0, alpha = 0.8) +
     # --- Significance bar with DCATS p-value ---
  geom_signif(
    comparisons = list(c("Chronic", "Symptomatic")),
    annotations = paste0("p = ", signif(fdr_val, 3)),
    y_position = max(glm_results_plot$Proportion) * 1.08,
    tip_length = 0.01,
    textsize = 4,
    color = "black"     # <- makes the bar black
  ) +
  
  # Colors
  scale_color_manual(values = cols) +
  
  theme_minimal(base_size = 12) +
  labs(
    title = "Triplet proportions per condition",
    y = "Proportion",
    x = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16)
  )

print(glm_results_dotplot)            
              

# Save figure
ggsave(plot = glm_results_dotplot,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "DCATS_triplet_dotplot.png")),
       width = 8, height = 6)





# Visualisation of doublet proportion comparison between conditions

# Extracting the number of monocytes for each individual
doublet_counts <- count_mat[, "doublet"]

# Retrieving the condition column in the design_mat
conditions <- design_mat$condition[match(rownames(count_mat), rownames(design_mat))]

# Build a data frame for plotting
glm_results_plot <- data.frame(
  Sample = rownames(count_mat),
  Condition = conditions,
  Count = doublet_counts
)

# Optional: compute proportion if you want
glm_results_plot$Proportion <- glm_results_plot$Count / rowSums(count_mat)

# Get DCATS FDR for annotation
fdr_val <- glm_results_multiplet_nbr_df$FDR[glm_results_multiplet_nbr_df$Cell_Type == "doublet"]

# Custom colors
cols <- c("Chronic" = "orange", "Symptomatic" = "red")

# Dot plot per sample with group mean + FDR annotation
glm_results_dotplot <- ggplot(glm_results_plot, aes(x = Condition, y = Proportion, color = Condition)) +
  geom_jitter(width = 0.1, size = 3, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.4,
               color = "black", fatten = 0, alpha = 0.8) +
     # --- Significance bar with DCATS p-value ---
  geom_signif(
    comparisons = list(c("Chronic", "Symptomatic")),
    annotations = paste0("p = ", signif(fdr_val, 3)),
    y_position = max(glm_results_plot$Proportion) * 1.08,
    tip_length = 0.01,
    textsize = 4,
    color = "black"     # <- makes the bar black
  ) +
  
  # Colors
  scale_color_manual(values = cols) +
  
  theme_minimal(base_size = 12) +
  labs(
    title = "Doublet proportions per condition",
    y = "Proportion",
    x = ""
  ) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, size = 16)
  )

print(glm_results_dotplot)            
              

# Save figure
ggsave(plot = glm_results_dotplot,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "DCATS_doublet_dotplot.png")),
       width = 8, height = 6)







################## To compare the interaction capacity of classical monocytes and non-classical monocytes between conditions (DCATS tool)



# Construct a count matrix
count_mat <- table(sc10x$orig.ident, sc10x$mono_interaction_capacity)
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
knn_mat <- knn_simMat(sc10x@graphs$RNA_snn, sc10x$mono_interaction_capacity)

# Run DCATS GLM
glm_results <- dcats_GLM(
  count_mat,
  design_mat,
  similarity_mat = knn_mat
)
print(glm_results)

# Convert the list of lists into a single data frame
glm_results_inter_mono_df <- data.frame(
  Cell_Type = rownames(glm_results$coeffs),
  Coefficient = glm_results$coeffs[, 1],
  SE = glm_results$coeffs_err[, 1],
  LRT = glm_results$LR_vals[, 1],
  pval = glm_results$LRT_pvals[, 1],
  FDR = glm_results$fdr[, 1]
)

print(glm_results_inter_mono_df)


# Save to CSV file
write.csv(glm_results_inter_mono_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glm_results_inter_mono_batch.csv")), row.names = FALSE)







# Visualisation of proportion comparison between conditions


## For classical monocytes

# Build the classical monocyte dataframe
df_classical <- data.frame(
  Sample = rownames(count_mat),
  Condition = design_mat$condition[match(rownames(count_mat), rownames(design_mat))],
  Singlet = count_mat[, "Classical monocyte_singlet"],
  Multiplet = count_mat[, "Classical monocyte_multiplet"]
)

# Compute proportions
df_classical$Singlet_prop   <- df_classical$Singlet   / rowSums(count_mat)
df_classical$Multiplet_prop <- df_classical$Multiplet / rowSums(count_mat)

# Reshape to 4 rows per sample
df_classical_long <- df_classical %>%
  tidyr::pivot_longer(
    cols = c(Singlet_prop, Multiplet_prop),
    names_to = "Interaction",
    values_to = "Proportion"
  ) %>%
  mutate(
    Interaction = ifelse(Interaction == "Singlet_prop", "Singlet", "Multiplet"),
    Group = paste0(Condition, "_", Interaction)   # final 4 groups on x-axis
  )

# Retrieving DCATS results
p_singlet <- glm_results_inter_mono_df$FDR[
  glm_results_inter_mono_df$Cell_Type == "Classical monocyte_singlet"
]

p_multiplet <- glm_results_inter_mono_df$FDR[
  glm_results_inter_mono_df$Cell_Type == "Classical monocyte_multiplet"
]


# Ploting classical monocytes with 4 categories
cols <- c(
  "Chronic_Singlet"     = "#88CCEE",
  "Chronic_Multiplet"   = "#332288",
  "Symptomatic_Singlet" = "#E69F00",
  "Symptomatic_Multiplet" = "red"
)

dotplot_classical_mono <- ggplot(df_classical_long, aes(x = Group, y = Proportion, color = Group)) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", color = "black", width = 0.5) +
  scale_color_manual(values = cols) +
  
  # significance bars: singlet vs multiplet inside each condition
  geom_signif(
    comparisons = list(c("Chronic_Singlet", "Symptomatic_Singlet")),
    annotations = paste0("p=", signif(p_singlet, 3)),
    y_position = max(df_classical_long$Proportion) * 1.05,
    tip_length = 0.01,
    color = "black"
  ) +
  geom_signif(
    comparisons = list(c("Chronic_Multiplet", "Symptomatic_Multiplet")),
    annotations = paste0("p=", signif(p_multiplet, 3)),
    y_position = max(df_classical_long$Proportion) * 1.15,
    tip_length = 0.01,
    color = "black"
  ) +

  theme_minimal(base_size = 12) +
  labs(
    title = "Classical monocyte singlet vs multiplet by condition",
    x = "",
    y = "Proportion"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dotplot_classical_mono

# Save figure
ggsave(plot = dotplot_classical_mono,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "DCATS_classicalMono_singlet_multiplet_dotplot.png")),
       width = 8, height = 6)




## For non-classical monocytes

# Build the classical monocyte dataframe
df_non_classical <- data.frame(
  Sample = rownames(count_mat),
  Condition = design_mat$condition[match(rownames(count_mat), rownames(design_mat))],
  Singlet = count_mat[, "Non-classical monocyte_singlet"],
  Multiplet = count_mat[, "Non-classical monocyte_multiplet"]
)

# Compute proportions
df_non_classical$Singlet_prop   <- df_non_classical$Singlet   / rowSums(count_mat)
df_non_classical$Multiplet_prop <- df_non_classical$Multiplet / rowSums(count_mat)

# Reshape to 4 rows per sample
df_non_classical_long <- df_non_classical %>%
  tidyr::pivot_longer(
    cols = c(Singlet_prop, Multiplet_prop),
    names_to = "Interaction",
    values_to = "Proportion"
  ) %>%
  mutate(
    Interaction = ifelse(Interaction == "Singlet_prop", "Singlet", "Multiplet"),
    Group = paste0(Condition, "_", Interaction)   # final 4 groups on x-axis
  )

# Retrieving DCATS results
p_singlet_NC <- glm_results_inter_mono_df$FDR[
  glm_results_inter_mono_df$Cell_Type == "Non-classical monocyte_singlet"
]

p_multiplet_NC <- glm_results_inter_mono_df$FDR[
  glm_results_inter_mono_df$Cell_Type == "Non-classical monocyte_multiplet"
]


# Ploting classical monocytes with 4 categories
cols <- c(
  "Chronic_Singlet"     = "#88CCEE",
  "Chronic_Multiplet"   = "#332288",
  "Symptomatic_Singlet" = "#E69F00",
  "Symptomatic_Multiplet" = "red"
)

dotplot_non_classical_mono <- ggplot(df_non_classical_long, aes(x = Group, y = Proportion, color = Group)) +
  geom_jitter(width = 0.15, size = 3, alpha = 0.9) +
  stat_summary(fun = mean, geom = "crossbar", color = "black", width = 0.5) +
  scale_color_manual(values = cols) +
  
  # significance bars: singlet vs multiplet inside each condition
  geom_signif(
    comparisons = list(c("Chronic_Singlet", "Symptomatic_Singlet")),
    annotations = paste0("p=", signif(p_singlet_NC, 3)),
    y_position = max(df_non_classical_long$Proportion) * 1.05,
    tip_length = 0.01,
    color = "black"
  ) +
  geom_signif(
    comparisons = list(c("Chronic_Multiplet", "Symptomatic_Multiplet")),
    annotations = paste0("p=", signif(p_multiplet_NC, 3)),
    y_position = max(df_non_classical_long$Proportion) * 1.15,
    tip_length = 0.01,
    color = "black"
  ) +

  theme_minimal(base_size = 12) +
  labs(
    title = "Non-classical monocyte singlet vs multiplet by condition",
    x = "",
    y = "Proportion"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

dotplot_non_classical_mono

# Save figure
ggsave(plot = dotplot_non_classical_mono,
       filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "DCATS_NonClassicalMono_singlet_multiplet_dotplot.png")),
       width = 8, height = 6)








