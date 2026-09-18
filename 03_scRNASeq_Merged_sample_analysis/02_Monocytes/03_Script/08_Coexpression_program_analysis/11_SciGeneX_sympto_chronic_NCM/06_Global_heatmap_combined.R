

# Building a global heatmap of all modules accross all cells


###############################################################################
# Extract NC-Tcell modules
###############################################################################
scigx_res_NC_Tcell <- readRDS(file.path("/home/quatreville/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2/05_Output/08_Coexpression_program_analysis/13_SciGeneX_sympto_chronic_NC_Tcell/01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_13_SciGeneX_sympto_chronic_NC_Tcell__clusterset_object_module_activity.RDS"))

modules_NC_Tcell <- scigx_res_NC_Tcell@top_genes
names(modules_NC_Tcell) <- paste0("NCmono_Tcell_M", seq_along(modules_NC_Tcell))

###############################################################################
# Extract NC modules
###############################################################################
scigx_res_NCmono <- readRDS(file.path("/home/quatreville/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2/05_Output/08_Coexpression_program_analysis/11_SciGeneX_sympto_chronic_NCM/01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_11_SciGeneX_sympto_chronic_NCM__clusterset_object_module_activity.RDS")) 

modules_NCmono <- scigx_res_NCmono@top_genes
names(modules_NCmono) <- paste0("NCmono_M", seq_along(modules_NCmono))

###############################################################################
# Merge all modules
###############################################################################

all_modules <- c(
    modules_NCmono,
    modules_NC_Tcell
)
# Reloading the seurat objects for each group
sc10x_NC_Tcell = readRDS(file.path("/home/quatreville/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2/05_Output/08_Coexpression_program_analysis/13_SciGeneX_sympto_chronic_NC_Tcell/01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_13_SciGeneX_sympto_chronic_NC_Tcell__seuratObject_module_activity.RDS"))
sc10x_NCmono = readRDS(file.path("/home/quatreville/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2/05_Output/08_Coexpression_program_analysis/11_SciGeneX_sympto_chronic_NCM/01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_11_SciGeneX_sympto_chronic_NCM__seuratObject_module_activity.RDS"))

# Merging of the two seurat objects
sc10x = merge(sc10x_NC_Tcell, sc10x_NCmono)

# Keeping the information of which cell type each cluster comes from
sc10x$cell_origin <- ifelse(
    sc10x$multiplet_type == "Non-classical monocyte",
    "NCmono",
    "NCmono_Tcell"
)

sc10x$cluster_origin <- paste0(
    sc10x$cell_origin,
    "_",
    sc10x$heatmap_clusters
)

###############################################################################
# Computing the module score on all cells
###############################################################################


# Compute the UCell score for each cell and each module
  sc10x <- UCell::AddModuleScore_UCell(
    sc10x,
    features = all_modules,
    name = ""
  )


# ###########################################
# Compute activity of gene modules per cell
# ###########################################

# Defining expression matrix from sc10x object and not SciGeneX
expr <- sc10x@assays$RNA@data

# Mean computing
gene_means <- Matrix::rowMeans(expr)

# Binary expression assessment
activity_bin <- expr > gene_means

# Module activity scorign
module_activity_score_names <- c()
module_activity_names <- c()

for (module_name in names(all_modules)) {

  genes <- all_modules[[module_name]]

  genes <- intersect(genes, rownames(activity_bin))

  score <- Matrix::colMeans(activity_bin[genes, , drop = FALSE])

  score_bin <- as.numeric(score >= 0.5)

  sc10x <- AddMetaData(
    sc10x,
    metadata = score,
    col.name = paste0("Activity.Score.", module_name)
  )

  sc10x <- AddMetaData(
    sc10x,
    metadata = score_bin,
    col.name = paste0("Activity.Binary.", module_name)
  )
}

# ###########################################
# Compute activity of gene modules per cluster as the percentage of cells activating this module
# ###########################################


# Loop over the cluster to study them one by one
Idents( sc10x) = "cluster_origin"

# Creating module list
module_names <- names(all_modules)

# Computing cluster × module activity (% cells active)
clusters_modules_activity <- data.frame()

for (cluster in levels(Idents(sc10x))) {

  cluster_cells <- WhichCells(sc10x, idents = cluster)

  cat("\nCluster:", cluster, "n =", length(cluster_cells))

  for (module_name in module_names) {

    bin_col <- paste0("Activity.Binary.", module_name)

    values <- sc10x@meta.data[cluster_cells, bin_col]
    # mean(binary) = fraction of active cells
    percentage_activity <- mean(values, na.rm = TRUE)

    clusters_modules_activity <- rbind(
      clusters_modules_activity,
      data.frame(
        cluster = cluster,
        module = module_name,
        activity.pct = percentage_activity
      )
    )
  }
}

# Factor ordering
clusters_modules_activity$module <- factor(
  clusters_modules_activity$module,
  levels = module_names
)
clusters_modules_activity$cluster <- factor(
  clusters_modules_activity$cluster,
  levels = levels(Idents(sc10x))
)

# Formatting for heatmap plotting
library(tidyr)

clusters_modules_activity_wide <- pivot_wider(
  clusters_modules_activity,
  names_from = cluster,
  values_from = activity.pct
)

# Removing clusters that contain too few cells
cols_to_keep <- setdiff(
  colnames(clusters_modules_activity_wide),
  c("NCmono_Tcell_3", "NCmono_6")
)

clusters_modules_activity_wide_df <- as.data.frame(
  clusters_modules_activity_wide[, cols_to_keep, drop = FALSE]
)

# so that module column become row names and so that all column are numeric and pheatmap doesn't crash
rownames(clusters_modules_activity_wide_df) <- clusters_modules_activity_wide_df$module
clusters_modules_activity_wide_df$module <- NULL





###############################################################################
# Plotting the gloabl heatmap
###############################################################################
# Heatmap expects a matrix with a module column, a cluster_rigin column and an activity score
library(tidyr)

# Heatmap
png(
    file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_NCmono_vs_NCmonoTcell_module_activity_heatmap.png")
    ),
    width = 3000,
    height = 2500,
    res = 300
)

pheatmap(
    as.matrix(clusters_modules_activity_wide_df),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(
        brewer.pal(7, "YlOrRd")
    )(100)
)


dev.off()



# Ordering clusters and modules
# Extracting cluster names
cluster_levels <- colnames(clusters_modules_activity_wide_df)
#Ordering
cluster_levels <- c(
  grep("^NCmono_", cluster_levels, value = TRUE),
  grep("^NCmono_Tcell_", cluster_levels, value = TRUE)
)

# Ordering module names
module_levels <- names(all_modules)

module_levels <- c(
  grep("^NCmono_M", module_levels, value = TRUE),
  grep("^NCmono_Tcell_M", module_levels, value = TRUE)
)

# numeric sort within each group
module_levels <- c(
  module_levels[order(as.numeric(sub(".*_M", "", module_levels[grep("^NCmono_M", module_levels)])))],
  module_levels[order(as.numeric(sub(".*_M", "", module_levels[grep("^NCmono_Tcell_M", module_levels)])))]
)



png(
  file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "_NCmono_vs_NCmonoTcell_module_activity_heatmap_ordered.png")
  ),
  width = 3000,
  height = 2500,
  res = 300
)

pheatmap(
  as.matrix(clusters_modules_activity_wide_df),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(
    RColorBrewer::brewer.pal(7, "YlOrRd")
  )(100)
)

dev.off()

# In PDF
pdf(
    file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_NCmono_vs_NCmonoTcell_module_activity_heatmap_ordered.pdf")),
    width = 8,
  height = 6
)

# Draw heatmap
pheatmap(
  as.matrix(clusters_modules_activity_wide_df),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 35,
  cellheight = 35,
  color = colorRampPalette(
    RColorBrewer::brewer.pal(7, "YlOrRd")
  )(100)
)

# Close device
dev.off()







# ============================================================
# Module correlation matrix at the CELL level
# ============================================================

# Continuous custom activity scores
module_score_cols <- paste0(
  "Activity.Score.",
  names(all_modules)
)

# Keep only columns that actually exist
module_score_cols <- module_score_cols[
  module_score_cols %in% colnames(sc10x@meta.data)
]

# Extract cell × module matrix
mat <- as.matrix(
  sc10x@meta.data[, module_score_cols, drop = FALSE]
)

# Remove modules containing only NA
mat <- mat[, colSums(is.na(mat)) < nrow(mat), drop = FALSE]

# Spearman correlation across CELLS
cor_mat <- cor(
  mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

# Clean module names for display
clean_names <- sub(
  "^Activity\\.Score\\.",
  "",
  colnames(cor_mat)
)

rownames(cor_mat) <- clean_names
colnames(cor_mat) <- clean_names

# Save correlation matrix
write.csv(
  cor_mat,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "correlation_matrix_cell_level_module_activity_scoring.csv"
    )
  ),
  row.names = TRUE
)



library(ComplexHeatmap)
library(circlize)

global_corr_plot <- Heatmap(
  cor_mat,
  name = "Spearman",
  col = colorRamp2(
    c(-1, 0, 1),
    c("blue", "white", "red")
  )
)

png(
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "correlation_heatmap_cell_level_module_activity_scoring.png"
    )
  ),
  width = 3000,
  height = 2500,
  res = 300
)

draw(global_corr_plot)

dev.off()


### Scatterplot of specific correlation pairs


scatter_df <- sc10x@meta.data %>%
  select(
    M8 = Activity.Score.NCmono_Tcell_M8,
    M4 = Activity.Score.NCmono_Tcell_M4,
    multiplet_type
  ) %>%
  filter(
  multiplet_type == "Non-classical monocyte_T cell",
    !is.na(M8),
    !is.na(M4)
  )

p_scatter <- ggplot(
  scatter_df,
  aes(
    x = M8,
    y = M4,
    color = multiplet_type
  )
) +
  geom_point(
    alpha = 0.5,
    size = 1.5
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "NCmono_Tcell_M8 activity score",
    y = "NCmono_Tcell_M4 activity score",
    color = "Doublet status"
  ) + theme(
    legend.position = "none"
  )
  


ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, "NCmono_Tcell_M8_M4_scatterplot.png"),
  plot = p_scatter,
  width = 8,
  height = 6
)



### Scatterplot of specific correlation pairs Ucellscore


# Calculating the UCell Score
feature_list <- all_modules

sc10x <- UCell::AddModuleScore_UCell(
    sc10x,
    features = feature_list,
    name = ""
)

# Scatter plot
scatter_df <- sc10x@meta.data %>%
  select(
    M8 = Activity.Score.NCmono_Tcell_M8,
    M4 = Activity.Score.NCmono_Tcell_M4,
    multiplet_type
  ) %>%
  filter(
  multiplet_type == "Non-classical monocyte_T cell",
    !is.na(M8),
    !is.na(M4)
  )

p_scatter <- ggplot(
  scatter_df,
  aes(
    x = M8,
    y = M4,
    color = multiplet_type
  )
) +
  geom_point(
    alpha = 0.5,
    size = 1.5
  ) +
  theme_classic(base_size = 16) +
  labs(
    x = "NCmono_Tcell_M8 activity score",
    y = "NCmono_Tcell_M4 activity score",
    color = "Doublet status"
  ) + theme(
    legend.position = "none"
  )
  


ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, "NCmono_Tcell_M8_M4_scatterplot.png"),
  plot = p_scatter,
  width = 8,
  height = 6
)






###
        ## plot_sgx_modules_combined_scores function
####
plot_module_overlap <- function(
    seurat_obj,
    module1,
    module2,
    cluster_col = "cluster_origin"
){

    m1 <- as.character(seurat_obj[[module1]][,1])
    m2 <- as.character(seurat_obj[[module2]][,1])
    print(table(m1, m2, useNA = "ifany"))
    state <- rep(NA, length(m1))

    state[m1=="0" & m2=="0"] <- "Neither"
    state[m1=="1" & m2=="0"] <- paste0(module1, "_only")
    state[m1=="0" & m2=="1"] <- paste0(module2, "_only")
    state[m1=="1" & m2=="1"] <- "Both"
    print(table(state, useNA = "ifany"))
 
    df <- data.frame(
    cluster = seurat_obj[[cluster_col]][,1],
    state = state)
    df = df %>%
    filter(
        !cluster %in% c("NCmono_Tcell_3", "NCmono_6")
    )

    # Check which clusters remain
    print(sort(unique(df$cluster)))

    df_sum <- df %>%
        count(cluster, state) %>%
        group_by(cluster) %>%
        mutate(freq = n/sum(n))
    print(sort(unique(df_sum$cluster)))
    
    mat <- df_sum %>%
        select(cluster, state, freq) %>%
        pivot_wider(
            names_from = cluster,
            values_from = freq,
            values_fill = 0
        )
    mat2 <- as.matrix(mat[,-1])
    rownames(mat2) <- mat$state

    pheatmap(
        mat2,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        cellwidth = 25,
        cellheight = 25,
        fontsize_row = 11,
        fontsize_col = 10,
        main = paste(module1, "vs", module2)
    )
}
pdf(
    file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "NCM2_NCTcell_M4_combined_scores.pdf")),
    width = 10,
    height = 5
)

plot_module_overlap(
    sc10x,
    module1 = "Activity.Binary.NCmono_M2",
    module2 = "Activity.Binary.NCmono_Tcell_M4"
)

dev.off()





##### 
        #Plotting a signature in one module with the details (number of cells activated the module from each condition, score of these cells between  conditions)
#####

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
        label = paste0( n_active, " / ", n_total, "\n(", round(percentage, 1), "%)" )
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
      subtitle = "Module activation by condition",
      x = "Condition",
      y = "% of module-active cells"
    )
ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(module, "_", cluster, "_percentage_active_cell_globalHeatmap.pdf")),
  plot = p_composition,
  width = 12,
  height = 6,
  dpi = 600,
  device = cairo_pdf
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
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(module, "_", cluster, "_continuous_activity_score_globalHeatmap.pdf")),
  plot = p_score,
  width = 12,
  height = 6,
  dpi = 600,
  device = cairo_pdf
)

}

plot_module_cluster_condition(
    seurat_obj= sc10x,
    module= "NCmono_Tcell_M4",
    cluster = "NCmono_1")





##### 
        #Plotting the proportion of cells activating two modules in one cluster with the details (number of cells activated the module from each condition, score of these cells between  conditions)
#####




# ------------------------------------------------------------
# Proportion of cells activating BOTH modules,
# split by condition
# ------------------------------------------------------------

plot_both_modules_by_condition <- function(
    seurat_obj,
    module1,
    module2,
    cluster,
    condition_col = "condition",
    cluster_col = "cluster_origin"
){

  # ------------------------------------------------------------
  # 1. Check inputs
  # ------------------------------------------------------------

  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster column not found: ", cluster_col)
  }

  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop("Condition column not found: ", condition_col)
  }

  binary_col1 <- paste0("Activity.Binary.", module1)
  binary_col2 <- paste0("Activity.Binary.", module2)

  if (!binary_col1 %in% colnames(seurat_obj@meta.data)) {
    stop("Module binary score not found: ", binary_col1)
  }

  if (!binary_col2 %in% colnames(seurat_obj@meta.data)) {
    stop("Module binary score not found: ", binary_col2)
  }

  # ------------------------------------------------------------
  # 2. Extract metadata and restrict to selected cluster
  # ------------------------------------------------------------

  metadata <- seurat_obj@meta.data

  cluster_df <- metadata %>%
    filter(.data[[cluster_col]] == cluster) %>%
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
  # 3. Define "Both"
  # ------------------------------------------------------------

  cluster_df <- cluster_df %>%
    mutate(
      both = .data[[binary_col1]] == 1 &
             .data[[binary_col2]] == 1
    )

  # ------------------------------------------------------------
  # 4. Calculate proportion of "Both" cells
  #    WITHIN the selected cluster and condition
  # ------------------------------------------------------------

  composition_df <- cluster_df %>%
    group_by(.data[[condition_col]]) %>%
    summarise(
      n_total = n(),
      n_both = sum(both, na.rm = TRUE),
      proportion = n_both / n_total,
      percentage = 100 * proportion,
      .groups = "drop"
    )

  composition_df[[condition_col]] <- factor(
    composition_df[[condition_col]],
    levels = c("Chronic", "Symptomatic")
  )

  print(composition_df)

  # ------------------------------------------------------------
  # 5. Plot
  # ------------------------------------------------------------

  p_both <- ggplot(
    composition_df,
    aes(
      x = .data[[condition_col]],
      y = percentage,
      fill = .data[[condition_col]]
    )
  ) +
  geom_col(
    width = 0.65
  ) +
    geom_text(
      aes(
        label = paste0(
          n_both,
          " / ",
          n_total,
          "\n(",
          round(percentage, 1),
          "%)"
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
  coord_cartesian(
    ylim = c(0, 110),
    clip = "off"
    )+
    theme_classic(base_size = 16) +
    theme(
      legend.position = "none"
    ) +
    labs(
      title = paste(
        module1,
        "and",
        module2,
        "co-activation in",
        cluster
      ),
      subtitle = "Both-module activation by condition",
      x = "Condition",
      y = "% of cells in cluster"
    )

  # ------------------------------------------------------------
  # 6. Save
  # ------------------------------------------------------------

  ggsave(
    filename = file.path(
      PATH_ANALYSIS_OUTPUT,
      paste0(
        module1, "_", module2, "_",
        cluster,
        "_percentage_both_by_condition.pdf"
      )
    ),
    plot = p_both,
    width = 12,
    height = 6,
    dpi = 600,
    device = cairo_pdf
  )

  return(
    list(
      summary = composition_df,
      plot = p_both
    )
  )
}

plot_both_modules_by_condition(seurat_obj= sc10x,module1 = "NCmono_M2", module2 = "NCmono_Tcell_M4",  cluster = "NCmono_7")





# Same but grouping cluster 1, 2, 3 and 4, 5, 7

# ------------------------------------------------------------
# Proportion of cells activating BOTH modules,
# split by condition and cluster group
# ------------------------------------------------------------

plot_both_modules_by_condition_cluster_group <- function(
    seurat_obj,
    module1,
    module2,
    cluster_col = "cluster_origin",
    condition_col = "condition"
){

  # ------------------------------------------------------------
  # 1. Check inputs
  # ------------------------------------------------------------

  if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
    stop("Cluster column not found: ", cluster_col)
  }

  if (!condition_col %in% colnames(seurat_obj@meta.data)) {
    stop("Condition column not found: ", condition_col)
  }

  binary_col1 <- paste0("Activity.Binary.", module1)
  binary_col2 <- paste0("Activity.Binary.", module2)

  if (!binary_col1 %in% colnames(seurat_obj@meta.data)) {
    stop("Module binary score not found: ", binary_col1)
  }

  if (!binary_col2 %in% colnames(seurat_obj@meta.data)) {
    stop("Module binary score not found: ", binary_col2)
  }

  # ------------------------------------------------------------
  # 2. Extract metadata and define cluster groups
  # ------------------------------------------------------------

  metadata <- seurat_obj@meta.data

  metadata <- metadata %>%
    mutate(
      cluster_group = case_when(
        .data[[cluster_col]] %in% c("NCmono_1", "NCmono_2", "NCmono_3") ~ "Chronic clusters",
        .data[[cluster_col]] %in% c("NCmono_4", "NCmono_5", "NCmono_7") ~ "Symptomatic clusters",
        TRUE ~ NA_character_
      )
    ) %>%
    filter(
      !is.na(cluster_group),
      .data[[condition_col]] %in% c(
        "Chronic",
        "Symptomatic"
      )
    )

  if (nrow(metadata) == 0) {
    stop("No cells found in the selected cluster groups.")
  }

  # ------------------------------------------------------------
  # 3. Define "Both"
  # ------------------------------------------------------------

  metadata <- metadata %>%
    mutate(
      both = .data[[binary_col1]] == 1 &
             .data[[binary_col2]] == 1
    )

  # ------------------------------------------------------------
  # 4. Calculate proportion of "Both"
  #
  # Denominator:
  # total cells from each condition
  # within each cluster group
  # ------------------------------------------------------------

  composition_df <- metadata %>%
    group_by(
      cluster_group,
      .data[[condition_col]]
    ) %>%
    summarise(
      n_total = n(),
      n_both = sum(both, na.rm = TRUE),
      proportion = n_both / n_total,
      percentage = 100 * proportion,
      .groups = "drop"
    ) %>%
    mutate(
      cluster_group = factor(
        cluster_group,
        levels = c(
          "Chronic clusters",
          "Symptomatic clusters"
        )
      ),
      condition = factor(
        .data[[condition_col]],
        levels = c(
          "Chronic",
          "Symptomatic"
        )
      )
    )

  print(composition_df)

  # ------------------------------------------------------------
  # 5. Plot
  # ------------------------------------------------------------

  p_both <- ggplot(
    composition_df,
    aes(
      x = cluster_group,
      y = percentage,
      fill = condition
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.75),
      width = 0.65
    ) +
    geom_text(
      aes(
        label = paste0(
          n_both,
          " / ",
          n_total,
          "\n(",
          round(percentage, 1),
          "%)"
        )
      ),
      position = position_dodge(width = 0.75),
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
    coord_cartesian(
      ylim = c(0, 110),
      clip = "off"
    ) +
    theme_classic(base_size = 16) +
    theme(
      legend.position = "top"
    ) +
    labs(
      title = paste(
        module1,
        "and",
        module2,
        "co-activation"
      ),
      subtitle = "Both-module activation by condition and cluster group",
      x = "Cluster group",
      y = "% of cells activating both modules",
      fill = "Condition"
    )

  # ------------------------------------------------------------
  # 6. Save
  # ------------------------------------------------------------

  ggsave(
    filename = file.path(
      PATH_ANALYSIS_OUTPUT,
      paste0(
        module1,
        "_",
        module2,
        "_cluster_groups_percentage_both.pdf"
      )
    ),
    plot = p_both,
    width = 10,
    height = 6,
    dpi = 600,
    device = cairo_pdf
  )

  # ------------------------------------------------------------
  # 7. Return results
  # ------------------------------------------------------------

  return(
    list(
      summary = composition_df,
      plot = p_both
    )
  )
}


both_result <- plot_both_modules_by_condition_cluster_group(
  seurat_obj = sc10x,
  module1 = "NCmono_M2",
  module2 = "NCmono_Tcell_M4"
)


