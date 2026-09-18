## To compare multiplet proportion between the two conditions

# Per-sample multiplet counts
multiplet_sample_df <- sc10x_Bcell@meta.data %>%
  group_by(orig.ident, condition, sgx_multiplet) %>%
  summarise(count = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = sgx_multiplet,
    values_from = count,
    values_fill = 0
  ) %>%
  mutate(
    total = singlet + doublet,
    prop_multiplet = doublet / total
  )

print(multiplet_sample_df)

# Binomial glm for comparing

glmm_multiplet <- glmer(
  cbind(doublet, singlet) ~ condition + (1 | orig.ident),
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
write.csv(glmm_summary_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glmm_multiplet_fixed_effects_Bcell.csv")), row.names = FALSE)

# Saving random effects
ranef_df <- as.data.frame(ranef(glmm_multiplet)$orig.ident)
print(ranef_df)

# Save to CSV file
write.csv(ranef_df, file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "glmm_random_effects_results_Bcell.csv")), row.names = FALSE)


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
  
  geom_jitter(width = 0.1, size = 6, alpha = 0.9) +
  
  # Mean bar per group
  stat_summary(fun = mean, geom = "crossbar", width = 0.4, size = 1,
               color = "black", fatten = 0, alpha = 0.8) +
  
  # Significance bar with GLMM p-value
  geom_signif(
    comparisons = list(c("Chronic", "Symptomatic")),
    annotations = glmm_pvalue_label,
    y_position = max(glm_results_plot$Proportion) * 1.08,
    tip_length = 0.01,
    textsize = 7,
    color = "black"
  ) +
  
  scale_color_manual(values = cols) +
  
  theme_classic(base_size = 20) +
  labs(
    title = "Multiplet Proportion per Condition",
    y = "Proportion of multiplets",
    x = ""
  ) +
  theme(
    legend.position = "none",
    
    # Bigger axis tick labels
    axis.text.y = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.title.y = element_text(size = 20),
    axis.title.x = element_text(size = 20),
    plot.title = element_blank()
  ) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))

print(multiplet_dotplot)

# Save figure
ggsave(
  plot = multiplet_dotplot,
  filename = file.path(PATH_ANALYSIS_OUTPUT,
                       paste0(outputFilesPrefix, "GLMM_multiplet_proportions_dotplot.png")),
  width = 8, height = 6, dpi = 300
)


