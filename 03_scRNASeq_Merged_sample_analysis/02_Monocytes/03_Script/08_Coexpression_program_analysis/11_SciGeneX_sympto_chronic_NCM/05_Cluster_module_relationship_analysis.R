# ##############################################
# This script analyse the SciGeneX gene modules
# ##############################################


## .................................................................
## Analyze the content of genes in each SciGeneX modules
## and the activity of each SciGeneX Modules
## .................................................................

## @knitr heterogeneity_SGXGeneModules
options(warn = 1) # initially on 2
# Creating list of module name
sgx_module_name_set <- paste0(
    "Module_",
    seq_along(scigx_res@gene_clusters)
)

sgx_module_feature_list <- setNames(
    scigx_res@gene_clusters,
    sgx_module_name_set
)

# Creating a module activity name list
module_activity_names <- paste0(
    "Activity.SGXModule",
    seq_along(scigx_res@gene_clusters)
)

# Creating the list of genes for each module
sgx_module_name_set <- paste0(
    "Module_",
    seq_along(scigx_res@gene_clusters)
)

sgx_module_feature_list <- setNames(
    scigx_res@gene_clusters,
    sgx_module_name_set
)





cat(" \n \n")
cat("### Individual module analysis {.tabset .tabset-fade}")
cat(" \n \n")

for( module_index in 1: length( sgx_module_name_set)){
  

  # Get the name of the module
  module_name = sgx_module_name_set[ module_index]
  
    # Create a directory
  module_dir <- file.path(PATH_ANALYSIS_OUTPUT, module_name)
    dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)
    
    
  cat(" \n \n")  
  cat("#### ", module_name, " {.tabset .tabset-fade}")
  cat(" \n \n")
  
  cat(" \n \n")
  cat("##### Module top genes list")
  cat(" \n \n")
  
  # Get the list of module genes
  gene_ids = sgx_module_feature_list[[ module_name]]
  
  # Display the genes as a list

  cat("<div style='width:48%; display:inline-block; vertical-align:top;'>")
  print( table( sc10x@meta.data[[ module_activity_names[ module_index]]]) %>% kable( caption = paste( "Cell count of", module_name)) %>% kable_styling( full_width = FALSE))
  
  # Saving the count table 
  cell_count_df <- as.data.frame(
  table(sc10x@meta.data[[module_activity_names[module_index]]]))
    write.csv(
  cell_count_df,
  file.path(module_dir, "cell_counts.csv"),
  row.names = FALSE
    )
  cat("</div>")
    
  cat("<div style='width:48%; display:inline-block; vertical-align:top;'>")
  # df <- matrix( sort( gene_ids), ncol = 4, byrow = TRUE) |> as.data.frame()
  genes_sorted <- sort(gene_ids)

remainder <- length(genes_sorted) %% 4

if (remainder != 0) {
  genes_sorted <- c(
    genes_sorted,
    rep(NA, 4 - remainder)
  )
}

df <- as.data.frame(
  matrix(
    genes_sorted,
    ncol = 4,
    byrow = TRUE
  )
)
    write.csv(
  data.frame(gene = sort(gene_ids)),
  file.path(module_dir, "top_genes_long.csv"),
  row.names = FALSE
)


  print( kable(df, col.names = NULL, caption = paste( "Top genes of", module_name)) %>%  kable_styling(full_width = FALSE, position = "center", font_size = 12))
  cat("</div>")

  cat(" \n \n")  
  cat("##### UMAP of Module Activity")
  cat(" \n \n")
  
  print( DimPlot( sc10x, group.by = module_activity_names[ module_index]) + ggtitle( paste( "Activity of", module_name)))  
  # Saving
  p1 <- DimPlot(
  sc10x,
  group.by = module_activity_names[module_index]
) +
ggtitle(paste("Activity of", module_name))

ggsave(
  filename = file.path(module_dir, "UMAP_activity.png"),
  plot = p1,
  width = 8,
  height = 6
)

print(p1)


  print( DimPlot( sc10x, group.by = module_activity_names[ module_index], split.by = module_activity_names[ module_index]) + ggtitle( paste( "Activity of", module_name))) 
 # Saving
  p2 <- DimPlot(
  sc10x,
  group.by = module_activity_names[module_index],
  split.by = module_activity_names[module_index]
) +
ggtitle(paste("Activity of", module_name))

ggsave(
  filename = file.path(module_dir, "UMAP_activity_split.png"),
  plot = p2,
  width = 12,
  height = 6
)

print(p2)
  
  cat(" \n \n")  
  cat("##### UMAP of Module Score")
  cat(" \n \n")
  
  print( FeaturePlot( sc10x, feature = module_name) + 
           ggtitle( paste( "Score of", module_name)))  
  # Saving
  p3 <- FeaturePlot(
  sc10x,
  feature = module_name
) +
ggtitle(paste("Score of", module_name))

ggsave(
  filename = file.path(module_dir, "UMAP_score.png"),
  plot = p3,
  width = 8,
  height = 6
)

print(p3)
  
  cat(" \n \n")
  cat("##### Functional enrichment analysis")
  cat(" \n \n")
  
  # Compute GO BP functionnal enrichment
  # ...................................
  ego_bp <- enrichGO(gene          = gene_ids,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "SYMBOL",
                     ont           = "BP",        # Biological Process
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = FALSE)

  # Convert result to dataframe
  ego_df <- as.data.frame( ego_bp)
  
  # Check if some result is found
  if( nrow( ego_df) > 1){
  
    # Display result dotplot
    print( dotplot(ego_bp, showCategory = 10, title = "GO Biological Process Enrichment"))
    
    # Saving
    p_go <- dotplot(
  ego_bp,
  showCategory = 10,
  title = "GO Biological Process Enrichment"
    )

    ggsave(
  filename = file.path(module_dir, "GO_BP_dotplot.png"),
  plot = p_go,
  width = 10,
  height = 7
    )

    print(p_go)

    
    # Use REVIGO to reduce the GOTerms
    # ...................................
    # 0. Prepare GO IDs and scores (-log10 p.adjust)
    go_ids <- ego_df$ID
    scores <- setNames( -log10(ego_df$p.adjust), go_ids)
    
    # 1. Compute the semantic similarity matrix
    simMatrix <- calculateSimMatrix( go_ids,
                                     orgdb    = "org.Hs.eg.db",
                                     ont      = "BP",
                                     method   = "Rel")  # or "Resnik", "Lin", etc.
    
    # 2. Reduce redundant GO terms
    reduced <- reduceSimMatrix(simMatrix,
                               scores    = scores,
                               threshold = 0.7,               # similarity threshold (0–1)
                               orgdb     = "org.Hs.eg.db")
    
    reduced = reduced[ which( !duplicated( reduced$go)), ]
    
    # 3. Display the results as Treemap
    treemapPlot(reduced)
    # Saving
    png(
  file.path(module_dir, "REVIGO_treemap.png"),
  width = 2000,
  height = 1600,
  res = 300
    )

    treemapPlot(reduced)

    dev.off()

    cat(" \n \n")  
  }else{
    cat("<BR>No functional enrichment for this module.")
  }
}










###
        ## plot_sgx_modules_combined_scores function
####
plot_module_overlap <- function(
    seurat_obj,
    module1,
    module2,
    cluster_col = "heatmap_clusters"
){

    m1 <- as.character(seurat_obj[[module1]][,1])
    m2 <- as.character(seurat_obj[[module2]][,1])

    state <- rep(NA, length(m1))

    state[m1=="0" & m2=="0"] <- "Neither"
    state[m1=="1" & m2=="0"] <- paste0(module1, "_only")
    state[m1=="0" & m2=="1"] <- paste0(module2, "_only")
    state[m1=="1" & m2=="1"] <- "Both"

    df <- data.frame(
        cluster = seurat_obj[[cluster_col]][,1],
        state = state
    )
    df <- df %>%
    filter(cluster != 6)

    df_sum <- df %>%
        count(cluster, state) %>%
        group_by(cluster) %>%
        mutate(freq = n/sum(n))

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
        main = paste(module1, "vs", module2)
    )
}
pdf(
    file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix,"_mod2_mod3_combined_scores.pdf")),
   width = 6,
  height = 4
)

plot_module_overlap(
    sc10x,
    module1 = "Activity.SGXModule3",
    module2 = "Activity.SGXModule2"
)

dev.off()





cat("\n \n")


## ..........................................................................
## Analyze the relationship between clusters and sciegenex module activity
## ..........................................................................

## @knitr heterogeneity_SGXGeneModules_vs_clusters

cat(" \n \n")
cat("### Percentage of cells activating module in clusters")
cat(" \n \n")

# Loop over the cluster to study them one by one
Idents( sc10x) = "heatmap_clusters"
clusters_modules_activity = data.frame()
for( cluster in levels( Idents( sc10x))){

  # Get the list of cells in the cluster
  cluster_cells = names( Idents(sc10x))[ which( Idents(sc10x) == cluster)]
  cat("There are this number of cells in this cluster", length(cluster_cells))
  # For each SGX module, compute the number of cells in activity in the cluster
  for( module_index in 1:length( sgx_module_name_set)){
    
    module_activity_name = module_activity_names[ module_index]
    cluster_module_activity_table = table( sc10x@meta.data[ module_activity_name][ cluster_cells,])
    percentage_activity = cluster_module_activity_table[ "1"]/( cluster_module_activity_table[ "0"] + cluster_module_activity_table[ "1"])
    clusters_modules_activity = rbind( clusters_modules_activity, data.frame( cluster = cluster,
                                                                              module = sgx_module_name_set[ module_index],
                                                                              activity.pct = percentage_activity
                                                                              ))
  }
}

clusters_modules_activity$module = factor( clusters_modules_activity$module, levels = sgx_module_name_set)
clusters_modules_activity$cluster = factor( clusters_modules_activity$cluster, levels = levels( Idents( sc10x)))
clusters_modules_activity_wide = spread(clusters_modules_activity, cluster, activity.pct)
row.names( clusters_modules_activity_wide) = clusters_modules_activity_wide$module
clusters_modules_activity_wide$module = NULL

# Removing uninformative modules
modules_to_remove <- c(
    "Module_7",
    "Module_10",
    "Module_14",
    "Module_15"
)

clusters_modules_activity_wide <-
    clusters_modules_activity_wide[
        !rownames(clusters_modules_activity_wide) %in% modules_to_remove,
        ,
        drop = FALSE
    ]
    
    
# heatmap of module activity by cluster
print( pheatmap( clusters_modules_activity_wide, 
          cluster_rows = FALSE, 
          cluster_cols = FALSE, 
          cellwidth = 20, 
          cellheight = 20,
          colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))
       )

# Saving
# Open device
png(
    file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "module_activity_heatmap_by_cluster.png")),
    width = 2500,
    height = 1800,
    res = 300
)

# Draw heatmap
pheatmap(
    clusters_modules_activity_wide,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    cellwidth = 20,
    cellheight = 20,
    color = colorRampPalette(
        brewer.pal(n = 7, name = "YlOrRd")
    )(100)
)

# Close device
dev.off()


# Removing cluster 6:
clusters_modules_activity_wide <- clusters_modules_activity_wide[
    ,
    !colnames(clusters_modules_activity_wide) %in% c("6")
]

# Saving
# Open device
pdf(
    file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "module_activity_heatmap_by_cluster_without_6.pdf")),
    width = 8,
  height = 6
)

# Draw heatmap
pheatmap(
    clusters_modules_activity_wide,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    cellwidth = 20,
    cellheight = 20,
    color = colorRampPalette(
        brewer.pal(n = 7, name = "YlOrRd")
    )(100)
)

# Close device
dev.off()






# If I want to know for a specific cluster and 2 spcific modules if the same cells are expressing both modules


Idents(sc10x) <- "heatmap_clusters"

cluster1_cells <- WhichCells(sc10x, idents = "1")

m1 <- sc10x$Activity.SGXModule1[cluster1_cells]
m2 <- sc10x$Activity.SGXModule2[cluster1_cells]

table(m1, m2)



# Checking for cell composition between clusters (individuals and condition)

# Getting the percentage of cells belongign to each individual in each cluster
# Extract metadata
metadata <- sc10x@meta.data

# Count cells per cluster and individual
cluster_composition_indiv <- metadata %>%
  group_by(heatmap_clusters, orig.ident) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  ungroup()

# Calculate percentage of each individual within each cluster
cluster_composition_indiv <- cluster_composition_indiv %>%
  group_by(heatmap_clusters) %>%
  mutate(percent = 100 * cell_count / sum(cell_count)) %>%
  select(heatmap_clusters, orig.ident, percent) %>%
  pivot_wider(
    names_from = orig.ident,
    values_from = percent,
    values_fill = 0
  )






# Correlation between modules (only the ones from module_activity scoring)
module_scores <- grep("^Module[0-9]+$", colnames(sc10x@meta.data), value = TRUE)

# Remove Module11
module_scores <- setdiff(module_scores, "Module11")


patient_modules_activity_scoring <- sc10x@meta.data %>%
  group_by(orig.ident, condition) %>%
  summarise(across(all_of(module_scores), median))
  
matrix <- as.matrix(patient_modules_activity_scoring[, module_scores])

matrix_filtered <- matrix[, colSums(is.na(matrix)) == 0]

 cor_mat <- cor(matrix_filtered, use="pairwise.complete.obs", method="spearman") # option to handle NAs
 
 # Save correlation matrix as CSV
write.csv(
  cor_mat,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "correlation_matrix_module_activity_scoring.csv")
  ),
  row.names = TRUE
)


library(pheatmap)
cor_mat[is.na(cor_mat)] <- 0

correlation_plot <- Heatmap(
    cor_mat,
    name = "Spearman",
    col = colorRamp2(
        c(-1, 0, 1),
        c("blue", "white", "red")
    )
)

png(
    file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "correlation_heatmap_module_activity_scoring.png")),
    width = 3000,
    height = 2500,
    res = 300
)

draw(correlation_plot)

dev.off()




# Correlation between modules (only the ones from module_activity binary)
module_scores <- grep("^Activity.Score.SGXModule", colnames(sc10x@meta.data), value = TRUE)
patient_modules_activity_scoring <- sc10x@meta.data %>%
  group_by(orig.ident, condition) %>%
  summarise(across(all_of(module_scores), median))
  
matrix <- as.matrix(patient_modules_activity_scoring[, module_scores])

matrix_filtered <- matrix[, colSums(is.na(matrix)) == 0]

 cor_mat <- cor(matrix_filtered, use="pairwise.complete.obs", method="spearman") # option to handle NAs
  
library(pheatmap)
cor_mat[is.na(cor_mat)] <- 0

correlation_plot <- Heatmap(
    cor_mat,
    name = "Spearman",
    col = colorRamp2(
        c(-1, 0, 1),
        c("blue", "white", "red")
    )
)

png(
    file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "correlation_heatmap_module_activity_binary.png")),
    width = 3000,
    height = 2500,
    res = 300
)

draw(correlation_plot)

dev.off()

# Saving clusterset object andseurat object
saveRDS(scigx_res, file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_clusterset_object_module_activity.RDS")))
saveRDS(sc10x, file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_seuratObject_module_activity.RDS")))







# ============================================================
# Module correlation matrix at the CELL level
# ============================================================


# Correlation between modules (only the ones from module_activity scoring)
module_scores <- grep("^Module[0-9]+$", colnames(sc10x@meta.data), value = TRUE)
module_scores <- setdiff(
    module_scores,
    c("Module2", "Module9", "Module12","Module14", "Module16" )
)
# Extract cell × module matrix
mat <- as.matrix(
  sc10x@meta.data[, module_scores, drop = FALSE]
)
# Remove modules containing only NA
mat <- mat[
  ,
  colSums(is.na(mat)) < nrow(mat),
  drop = FALSE
]

# Spearman correlation across CELLS
cor_mat <- cor(
  mat,
  method = "spearman",
  use = "pairwise.complete.obs"
)

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


library(pheatmap)
cor_mat[is.na(cor_mat)] <- 0

correlation_plot <- Heatmap(
    cor_mat,
    name = "Spearman",
    col = colorRamp2(
        c(-1, 0, 1),
        c("blue", "white", "red")
    )
)

pdf(
    file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "correlation_heatmap_module_activity_scoring.pdf")),
  width = 10,
  height = 8
)

draw(correlation_plot)

dev.off()






























##### Plotting a signature in one module with the details (number of cells activated the module from each condition, score of these cells between conditions)


plot_module_cluster_condition <- function(
    seurat_obj,
    module,
    cluster,
    cluster_col = "heatmap_clusters",
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
binary_col <- paste0("Activity.SGX", module)

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
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(module, "_", cluster, "_percentage_active_cell.pdf")),
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
score_col <- paste0("Activity.Score.SGX", module)
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
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(module, "_", cluster, "_continuous_activity_score.pdf")),
  plot = p_score,
  width = 12,
  height = 6,
  dpi = 600,
  device = cairo_pdf
)

}

plot_module_cluster_condition(
    seurat_obj= sc10x,
    module= "Module3",
    cluster = "1")










# Same but grouping cluster 1, 2, 3 and 4, 5, 7


# ============================================================
# Define cluster groups
# ============================================================

metadata <- sc10x@meta.data

metadata <- metadata %>%
  mutate(
    cluster_group = case_when(
      heatmap_clusters %in% c(1, 2, 3) ~ "Chronic clusters",
      heatmap_clusters %in% c(4, 5, 7) ~ "Symptomatic clusters",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(cluster_group),
    condition %in% c("Chronic", "Symptomatic")
  )


# ============================================================
# Module to analyse
# ============================================================

module <- "Module2"

binary_col <- paste0("Activity.SGX", module)
score_col  <- paste0("Activity.Score.SGX", module)


# ============================================================
# Plot 1:
# Module activation by condition and cluster group
# ============================================================

composition_df <- metadata %>%
  group_by(cluster_group, condition) %>%
  summarise(
    n_total = n(),
    n_active = sum(.data[[binary_col]] == 1, na.rm = TRUE),
    proportion = n_active / n_total,
    percentage = 100 * proportion,
    .groups = "drop"
  ) %>%
  mutate(
    cluster_group = factor(
      cluster_group,
      levels = c("Chronic clusters", "Symptomatic clusters")
    ),
    condition = factor(
      condition,
      levels = c("Chronic", "Symptomatic")
    )
  )


p_composition <- ggplot(
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
        n_active,
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
    )+
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top"
  ) +
  labs(
    title = paste(module, "activity"),
    x = "Cluster group",
    y = "% of cells activating module",
    fill = "Condition"
  )


ggsave(
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      module,
      "_cluster_groups_percentage_active_cell.pdf"
    )
  ),
  plot = p_composition,
  width = 10,
  height = 6,
  dpi = 600,
  device = cairo_pdf
)


# ============================================================
# Plot 2:
# Continuous activity score by condition and cluster group
# ============================================================

score_df <- metadata %>%
  select(
    cluster_group,
    condition,
    score = all_of(score_col)
  ) %>%
  mutate(
    cluster_group = factor(
      cluster_group,
      levels = c("Chronic clusters", "Symptomatic clusters")
    ),
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
  facet_wrap(
    ~ cluster_group
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
    title = paste(module, "score"),
    subtitle = "Continuous activity score",
    x = "Condition",
    y = "Activity score"
  )


ggsave(
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      module,
      "_cluster_groups_continuous_activity_score.pdf"
    )
  ),
  plot = p_score,
  width = 10,
  height = 6,
  dpi = 600,
  device = cairo_pdf
)









