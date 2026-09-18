

# Subsetting objects by lineage 
sc10x_NCmono <- subset(sc10x, subset = cell_origin == "NCmono")
sc10x_NC_Tcell <- subset(sc10x, subset = cell_origin == "NCmono_Tcell")

# Choosing modules to plot
modules_to_plot_nc <- c(
  "NCmono_Tcell_M13", # FCGR3A
  "NCmono_M2" # S100A8
)

modules_to_plot_nc_tcell <- c(
  "NCmono_Tcell_M8", # LILRB1
  "NCmono_Tcell_M6", # S100A8
  "NCmono_Tcell_M4", #Tcell
  "NCmono_Tcell_M5" # Tcell activation
  
)


# Plotting data generation for NC_mono
df_NCmono <- sc10x_NCmono@meta.data %>%
  dplyr::select(orig.ident, all_of(modules_to_plot_nc)) %>%
  pivot_longer(
    cols = all_of(modules_to_plot_nc),
    names_to = "module",
    values_to = "score"
  )

# Visualization
NCmono_violin_indiv = ggplot(df_NCmono, aes(x = orig.ident, y = score, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2, alpha = 0.5) +
  facet_wrap(~module, scales = "free_y") +
  theme_classic() +
  labs(title = "NCmono module score by individual")
ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_NCmono_gene_signatures_by_indiv.png")
    ),
    plot = NCmono_violin_indiv,
    width = 8,
    height = 6,
    dpi = 300
)



# Plotting data generation for NC_Tcell
df_NC_Tcell <- sc10x_NC_Tcell@meta.data %>%
  dplyr::select(orig.ident, all_of(modules_to_plot_nc_tcell)) %>%
  pivot_longer(
    cols = all_of(modules_to_plot_nc_tcell),
    names_to = "module",
    values_to = "score"
  )

# Visualization
NC_Tcell_violin_indiv = ggplot(df_NC_Tcell, aes(x = orig.ident, y = score, fill = orig.ident)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2, alpha = 0.5) +
  facet_wrap(~module, scales = "free_y") +
  theme_classic() +
  labs(title = "NCmono module activity by orig.ident")
ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_NCmono_Tcell_gene_signatures_by_indiv.png")
    ),
    plot = NC_Tcell_violin_indiv,
    width = 8,
    height = 6,
    dpi = 300
)


# By condition/symptomatic individuals





# NC

# Creating a new metadata column for the grouping
sc10x_NCmono$signature_group = ifelse(
    sc10x$condition == "Chronic",
    "Chronic",
    sc10x$orig.ident
)

# Plotting data generation for NC_Tcell
df_NCmono <- sc10x_NCmono@meta.data %>%
  dplyr::select(signature_group, all_of(modules_to_plot_nc)) %>%
  pivot_longer(
    cols = all_of(modules_to_plot_nc),
    names_to = "module",
    values_to = "score"
  )

# Visualization
NCmono_violin_grouped = ggplot(df_NCmono, aes(x = signature_group, y = score, fill = signature_group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2, alpha = 0.5) +
  facet_wrap(~module, scales = "free_y") +
  theme_classic() +
  labs(title = "NCmono module activity by orig.ident")
ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_NCmono_gene_signatures_grouped.png")
    ),
    plot = NCmono_violin_grouped,
    width = 8,
    height = 6,
    dpi = 300
)



# NC_T

# Creating a new metadata column for the grouping
sc10x_NC_Tcell$signature_group = ifelse(
    sc10x$condition == "Chronic",
    "Chronic",
    sc10x$orig.ident
)

# Plotting data generation for NC_Tcell
df_NC_Tcell <- sc10x_NC_Tcell@meta.data %>%
  dplyr::select(signature_group, all_of(modules_to_plot_nc_tcell)) %>%
  pivot_longer(
    cols = all_of(modules_to_plot_nc_tcell),
    names_to = "module",
    values_to = "score"
  )

# Visualization
NC_Tcell_violin_grouped = ggplot(df_NC_Tcell, aes(x = signature_group, y = score, fill = signature_group)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.size = 0.2, alpha = 0.5) +
  facet_wrap(~module, scales = "free_y") +
  theme_classic() +
  labs(title = "NCmono-T module activity by orig.ident")
ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_NCmono_Tcell_gene_signatures_grouped.png")
    ),
    plot = NC_Tcell_violin_grouped,
    width = 8,
    height = 6,
    dpi = 300
)







### Visualization of all signatures on the same figure


##############################################################################
# NCmono-Tcell modules
##############################################################################

df_NC_Tcell <- sc10x_NC_Tcell@meta.data %>%
  dplyr::select(signature_group, all_of(modules_to_plot_nc_tcell)) %>%
  pivot_longer(
    cols = all_of(modules_to_plot_nc_tcell),
    names_to = "module",
    values_to = "score"
  ) %>%
  mutate(panel = "NCmono-Tcell")

##############################################################################
# NCmono modules
##############################################################################

df_NCmono <- sc10x_NCmono@meta.data %>%
  dplyr::select(signature_group, all_of(modules_to_plot_nc)) %>%
  pivot_longer(
    cols = all_of(modules_to_plot_nc),
    names_to = "module",
    values_to = "score"
  ) %>%
  mutate(panel = "NCmono")

##############################################################################
# Merge
##############################################################################

df_all <- bind_rows(df_NC_Tcell, df_NCmono)

##############################################################################
# Order individuals
##############################################################################

patient_order <- c(
  # chronic first
  unique(df_all$signature_group[df_all$signature_group == "Chronic"]),
  # symptomatic second
  unique(df_all$signature_group[df_all$signature_group == "15_scRNASeq_A76_S_tp2"]),
  unique(df_all$signature_group[df_all$signature_group == "12_scRNASeq_A86_S_tp2"]),
  unique(df_all$signature_group[df_all$signature_group == "08_scRNASeq_A82_S_tp2"])
)

df_all$signature_group <- factor(
  df_all$signature_group,
  levels = patient_order
)

# Define custom colors
sig_group_colors <- c(
  "Chronic" = "#E69F00",      # orange
  "15_scRNASeq_A76_S_tp2" = "#CC79A7", # purple/pink
  "12_scRNASeq_A86_S_tp2" = "#AA4499", # violet
  "08_scRNASeq_A82_S_tp2" = "#882255" # dark purple
)


##############################################################################
# Order modules
##############################################################################

module_order <- c(
  modules_to_plot_nc_tcell,
  modules_to_plot_nc
)

df_all$module <- factor(
  df_all$module,
  levels = module_order
)

##############################################################################
# Plot
##############################################################################

vlnplot_all = ggplot(
  df_all,
  aes(
    x = signature_group,
    y = score,
    fill = signature_group
  )
) +
  geom_violin(
    scale = "width",
    trim = TRUE
  ) +
  geom_boxplot(
    width = 0.10,
    alpha = 0.5,
    outlier.size = 0.2
  ) +
   scale_fill_manual(values = sig_group_colors) +
   scale_color_manual(values = sig_group_colors) +
    facet_grid(
  module ~ .,
  scales = "free_y"
) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    strip.text.y = element_text(
      face = "bold"
    ),
    strip.text.x = element_text(
      size = 14,
      face = "bold"
    )
  ) +
  labs(
    x = "Individual",
    y = "Module score"
  )

ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_all_signatures_grouped.png")
    ),
    plot = vlnplot_all,
    width = 10,
    height = 15,
    dpi = 300
)

















##### Plotting a signature in one module with the details (number of cells activated the module from each condition, score of these cells between conditions)

plot_module_cluster_condition(
    seurat_obj= sc10x,
    module= "NCmono_Tcell_M8",
    cluster = "NCmono_Tcell_5")

plot_module_cluster_condition <- function(
    seurat_obj,
    module,
    cluster,
    cluster_col = "cluster_origin",
    condition_col = "condition",
    score_cutoff = 0.5,
    min_cells = 10
) {

  library(dplyr)
  library(ggplot2)
  library(patchwork)

  # ------------------------------------------------------------
  # 1. Check inputs
  # ------------------------------------------------------------

  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster column not found: ", cluster_col)
  }

  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop("Condition column not found: ", condition_col)
  }

  if (!module %in% colnames(seurat_obj@meta.data)) {
    stop(
      "Module score not found: ", module,
      "\nAvailable module scores include:\n",
      paste(
        grep(
          "^NCmono_|^NCmono_Tcell_",
          colnames(seurat_obj@meta.data),
          value = TRUE
        ),
        collapse = ", "
      )
    )
  }

  # ------------------------------------------------------------
  # 2. Extract metadata
  # ------------------------------------------------------------

  metadata <- seurat_obj@meta.data

  # Restrict to selected cluster
  cluster_df <- metadata %>%
    filter(.data[[cluster_col]] == cluster)

  if (nrow(cluster_df) < min_cells) {
    warning(
      "Only ", nrow(df),
      " cells found in cluster ", cluster
    )
  }

  # Keep only the two conditions of interest
  cluster_df <- cluster_df %>%
    filter(.data[[condition_col]] %in% c(
      "Chronic",
      "Symptomatic"
    ))

  if (nrow(cluster_df) == 0) {
    stop(
      "No Chronic/Symptomatic cells found in cluster ",
      cluster
    )
  }

  # ------------------------------------------------------------
  # 3. Plot 1:
  #    Origin of cells ACTIVATING the module
  # ------------------------------------------------------------
binary_col <- paste0("Activity.Binary.", module)

  composition_df <- cluster_df %>%
  group_by(.data[[condition_col]]) %>%
  summarise(
    n_total = n(),
    n_active = sum(.data[[binary_col]] == 1, na.rm = TRUE),
    proportion = n_active / n_total,
    percentage = 100 * proportion,
    .groups = "drop"
  )

  composition_df[[condition_col]] <- factor(
    composition_df[[condition_col]],
    levels = c("Chronic", "Symptomatic")
  )

  p_composition <- ggplot(
    composition_df,
    aes(
      x = .data[[condition_col]],
      y = percentage,
      fill = .data[[condition_col]]
    )
  ) +
    geom_col(width = 0.65) +
    geom_text(
      aes(
        label = paste0(
          round(percentage, 1),
          "%\n(n = ", n_active, ")"
        )
      ),
      vjust = -0.2,
      size = 5
    ) +
    scale_fill_manual(
      values = c(
        "Chronic" = "#E69F00",
        "Symptomatic" = "#D55E00"
      )
    ) +
    scale_y_continuous(
      limits = c(0, 100),
      expand = expansion(mult = c(0, 0.08))
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none"
    ) +
    labs(
      title = paste(
        module,
        "active cells in",
        cluster
      ),
      subtitle = "Condition composition among module-active cells",
      x = "Condition",
      y = "% of module-active cells"
    )
ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(module, "_", cluster, "_percentage_active_cell.png")),
  plot = p_composition,
  width = 12,
  height = 6
)

  # ------------------------------------------------------------
  # 4. Plot 2:
  #    Continuous activity score by condition
  # ------------------------------------------------------------
score_col <- paste0("Activity.Score.", module)
  score_df <- cluster_df %>%
    select(
      condition = all_of(condition_col),
      score = all_of(score_col)
    ) %>%
    mutate(
      condition = factor(
        condition,
        levels = c("Chronic", "Symptomatic")
      )
    )

  p_score <- ggplot(
    score_df,
    aes(
      x = condition,
      y = score,
      fill = condition
    )
  ) +
    geom_violin(
      scale = "width",
      trim = TRUE,
      alpha = 0.7
    ) +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA,
      alpha = 0.8
    ) +
    geom_jitter(
      width = 0.08,
      size = 0.6,
      alpha = 0.25
    ) +
    scale_fill_manual(
      values = c(
        "Chronic" = "#E69F00",
        "Symptomatic" = "#D55E00"
      )
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none"
    ) +
    labs(
      title = paste(
        module,
        "score in",
        cluster
      ),
      subtitle = "Continuous activity score",
      x = "Condition",
      y = "Activity score"
    )
ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(module, "_", cluster, "_continuous_activity_score.png")),
  plot = p_score,
  width = 12,
  height = 6
)

}

plot_module_cluster_condition(
    seurat_obj= sc10x,
    module= "NCmono_Tcell_M8",
    cluster = "NCmono_Tcell_1")
















# If we want the percentage of each condition compared to the total number of cells
binary_col <- paste0("Activity.Binary.", module)
  active_df <- cluster_df %>%
    filter(.data[[binary_col]] == 1)

  if (nrow(active_df) == 0) {
    warning(
      "No cells activating ", module,
      " were found in cluster ", cluster
    )
  }

  composition_df <- active_df %>%
    count(.data[[condition_col]], name = "n_active") %>%
    mutate(
      proportion = n_active / nrow(cluster_df),
      percentage = 100 * proportion
    )

  composition_df[[condition_col]] <- factor(
    composition_df[[condition_col]],
    levels = c("Chronic", "Symptomatic")
  )


