##' Plot a Complex Heatmap for Single-Cell Data
##'
##' This function generates a clustered heatmap from a ClusterSet object,
##' with support for top genes, transformations (log2, scaling), and cell identity annotations.
##'
##' @param object A ClusterSet object containing expression data and gene clusters.
##' @param use_top_genes Logical. If TRUE, only top genes (from object@top_genes) are plotted.
##' @param to_log2 Logical. If TRUE, applies log2 transformation to the matrix.
##' @param standardizing Logical. If TRUE, standardizes the matrix (z-score per gene).
##' @param ceil Numeric. Maximum value for color scaling (values above are capped).
##' @param floor Numeric. Minimum value for color scaling (values below are capped).
##' @param pseudocount Numeric. Small value added to avoid log(0) errors.
##' @param cell_identity Vector. Cell type/identity labels for column annotation.
##' @param centering Logical. If TRUE, centers the matrix (subtract row means).
##'
##' @note
##' - The 'magick' package is recommended for better rasterization (install with `install.packages("magick")`). #v2.7.4
##' - Set `ht_opt$message = FALSE` to suppress ComplexHeatmap messages.
##'
##' @examples
##' plot_ggComplexHeatmap(my_cluster_set, use_top_genes = TRUE, to_log2 = TRUE)


meaningful_modules <- c(1:34)

# Function for colors
gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}


plot_ggComplexHeatmap <- function(object,
                                  use_top_genes = FALSE,
                                  to_log2 = FALSE,
                                  standardizing = FALSE,
                                  ceil = 1,
                                  floor = -1,
                                  pseudocount = 0.1,
                                  cell_identity = NULL,
                                  centering = TRUE) {
  
  # --- Input Validation ---
  # Check if object is a valid ClusterSet
  check_format_cluster_set(object)
  
  # Check if cell_identity length matches matrix columns (if provided)
  if (!is.null(cell_identity)) {
    if (length(cell_identity) != ncol(object@data)) {
      stop("Error: 'cell_identity' length must match the number of columns in object@data.")
    }
  } else {
    warning("No cell_identity provided. Column annotations will be omitted.")
  }
  
  # --- Helper Function: Generate ggplot-like discrete colors ---
  # Creates a vector of n distinct colors using HCL space (hue-chroma-luminance)
  gg_color_hue <- function(n) {
    hues <- seq(15, 375, length = n + 1)  # Avoid red (0°) for better contrast
    hcl(h = hues, l = 65, c = 100)[1:n]     # Lightness=65, Chroma=100 for vibrant colors
  }
  
  # --- Data Selection ---
  # Subset matrix to top genes if requested
  if (use_top_genes) {
    if (length(object@top_genes) == 0) {
      stop("Error: No top genes found in object@top_genes. Use top_gene() method first.")
    }
    m <- object@data[unlist(object@top_genes), , drop = FALSE]
  } else {
    m <- object@data
  }
  
  # --- Data Transformation ---
  # Apply log2, scaling, centering, etc. (handled by transform_mat)
  m <- transform_mat(
    m,
    to_log2 = to_log2,
    pseudocount = pseudocount,
    centering = centering,
    standardizing = standardizing,
    ceil = ceil,
    floor = floor
  )
  
  # --- Column Clustering (Cells) ---
  # Compute Euclidean distance between cells (columns) and perform hierarchical clustering
  dist_cols <- dist(t(m), method = "euclidean")
  hc_cols <- hclust(dist_cols, method = "complete")
  
  # --- Row Clustering (Genes) ---
  # Get gene clusters from the object and ensure they match matrix rows
  gclust <- rep(NA, nrow(m))
    names(gclust) <- rownames(m)

    for (i in seq_along(object@gene_clusters)) {
  genes <- intersect(
        object@gene_clusters[[i]],
        rownames(m)
    )
  gclust[genes] <- as.character(i)
    }

    gclust <- factor(gclust, levels = as.character(seq_along(object@gene_clusters)))
  names(gclust) <- rownames(m)
  
  



  # Sort cluster levels numerically (avoid factor ordering issues)
  lev_clust <- suppressWarnings(sort(unique(as.numeric(gclust))))
  # Keeping only relevant modules
  keep_genes <- names(gclust)[gclust %in% meaningful_modules]
    cat("Genes in matrix:", nrow(m), "\n")
cat("Genes to keep:", length(keep_genes), "\n")
cat("Present:", sum(keep_genes %in% rownames(m)), "\n")

missing_genes <- setdiff(keep_genes, rownames(m))

cat("Missing:", length(missing_genes), "\n")

head(missing_genes)
cat("keep_genes lenght", length(keep_genes))
cat("unique keep_genes lenght", length(unique(keep_genes)))

    m <- m[keep_genes, , drop = FALSE]
    gclust <- gclust[keep_genes]
    
  gclust <- factor(
    gclust[rownames(m)],
    levels = lev_clust,
    ordered = TRUE
  )
  
  # --- Cell Identity Colors ---
  # Generate colors for cell identities (if provided)
  if (!is.null(cell_identity)) {
    ids <- unique(cell_identity)
    colors <- gg_color_hue(length(ids))
    names(colors) <- ids
  } else {
    colors <- NULL
    cell_identity <- rep("Unknown", ncol(m))  # Fallback if no identity provided
  }
  
  # --- Reorder Genes by Module and Clustering ---
  # For each gene module, cluster genes and reorder them
  modules <- unique(gclust)
  ordered_genes <- character(0)  # Initialize empty vector
  
  for (mod in modules) {
    genes_mod <- names(gclust)[gclust == mod]
    genes_mod <- intersect(genes_mod, rownames(m))  # Ensure genes exist in matrix
    
    if (length(genes_mod) > 1) {
      # Subset matrix for current module and cluster rows
      sub_m <- m[genes_mod, , drop = FALSE]
      dist_rows <- dist(sub_m, method = "euclidean")
      hc_rows <- hclust(dist_rows, method = "complete")
      ordered_genes <- c(ordered_genes, genes_mod[hc_rows$order])
    } else {
      # Single gene: no clustering needed
      ordered_genes <- c(ordered_genes, genes_mod)
    }
  }
  
  # Reorder the full matrix
  m <- m[ordered_genes, , drop = FALSE]
  
  # --- Column Annotation (Cell Identity) ---
  # Create HeatmapAnnotation only if cell_identity is provided
  if (!is.null(cell_identity)) {
    ha_col <- HeatmapAnnotation(
      Identity = cell_identity,
      col = list(Identity = colors),
      annotation_name_side = "right"  # Improve readability
    )
  } else {
    ha_col <- NULL
  }
  
  # --- Build Heatmap ---
  ht <- Heatmap(
    m,
    name = "Expression",
    col = colors_for_gradient("Ju1"),  # Predefined gradient (blue-white-red)
    
    # Clustering
    cluster_columns = hc_cols,
    cluster_rows = FALSE,  # Rows are pre-ordered by module
    
    # Performance
    use_raster = TRUE,     # Faster rendering for large matrices
    raster_by_magick = ifelse(requireNamespace("magick", quietly = TRUE), TRUE, FALSE),
    
    # Row splitting by gene clusters
    row_split = gclust,
    row_gap = unit(2, "mm"),  # Gap between modules
    
    # Annotations
    top_annotation = ha_col,
    
    # Aesthetics
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,  # Column dendrogram is redundant (pre-clustered)
    heatmap_legend_param = list(title = "Expression", title_position = "leftcenter-rot")
  )
  
  # Draw heatmap (suppress ComplexHeatmap messages if desired)
  ht_opt$message = FALSE  # Disable default messages
  draw(ht)
  
  # Saving the heatmap
  png(file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix,
"heatmap_no_cluster.png")),
  width = 3000,
  height = 2500,
  res = 300
)

draw(ht)

dev.off()

list(
        heatmap = ht,
        hc_cols = hc_cols,
        gclust = gclust,
        m = m
    )
}






transform_mat <- function(
    m,
    to_log2 = FALSE,
    pseudocount = 1,
    centering = FALSE,
    standardizing = FALSE,
    ceil = NULL,
    floor = NULL
) {

    # Log transform
    if (to_log2) {
        m <- log2(m + pseudocount)
    }

    # Center rows
    if (centering) {
        m <- t(scale(t(m), center = TRUE, scale = FALSE))
    }

    # Z-score rows
    if (standardizing) {
        m <- t(scale(t(m), center = TRUE, scale = TRUE))
    }

    # Apply limits
    if (!is.null(ceil)) {
        m[m > ceil] <- ceil
    }

    if (!is.null(floor)) {
        m[m < floor] <- floor
    }

    m
}
object = res_cluster_subtype
sample_labels <- sc10x_subset$orig.ident
res_heatmap = plot_ggComplexHeatmap(res_cluster_subtype, cell_identity=sample_labels, use_top_genes = TRUE, to_log2 = TRUE)
hc_cols <- res_heatmap$hc_cols
gclust <- res_heatmap$gclust
m <- res_heatmap$m

# To cluster cells and visualize the dendrogram
plot(hc_cols, cex = 0.3)

# Clusterizing cells
cell_clusters <- cutree(hc_cols, k = 7)
png(file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_cell_dendrogram.png")
    ),
    width = 3000,
    height = 2000,
    res = 300)
plot(hc_cols, cex = 0.3)
rect.hclust(hc_cols, k = 7)

dev.off()


# Visualize clusters on umap
table(cell_clusters)
cluster_sizes <- as.data.frame(table(cell_clusters))
write.csv(
    cluster_sizes,
    file = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_cell_cluster_sizes.csv")
    ),
    row.names = FALSE
)


sc10x_subset$cell_clusters <- factor(cell_clusters)
umap_new_clusters = DimPlot(sc10x_subset, group.by = "cell_clusters")
ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_cell_clusters_umap.png")
    ),
    plot = umap_new_clusters,
    width = 8,
    height = 6,
    dpi = 300
)

# Visualize through the heatmap


# Creating colors for clusters:
cluster_colors <- c(
    "1" = "#E69F00", # orange
    "2" = "#88CCEE", # sky blue
    "3" = "#009E73", # green
    "4" = "#F0E442", # yellow
    "5" = "#332288", # deep blue-violet
    "6" =  "#CC79A7", # purple/pink
    "7" =   "#AA4499", # violet
    "8" = "#8DB600", # apple green
    "9" = "#882255" # dark purple
    
)

# Defining cluster order:
cell_clusters <- factor(
    cell_clusters,
    levels = c("1", "2", "3", "4", "5", "6", "7")
)


ha_col <- HeatmapAnnotation(
  Cluster = factor(cell_clusters),
  col = list(Cluster = cluster_colors)
)

  # --- Build Heatmap ---
  ht <- Heatmap(
    m,
    name = "Expression",
    col = colors_for_gradient("Ju1"),  # Predefined gradient (blue-white-red)
    
    # Clustering
    cluster_columns = TRUE,    #hc_cols
    column_split = cell_clusters, # to split clusters
    cluster_column_slices = FALSE,
    cluster_rows = FALSE,  # Rows are pre-ordered by module
    
    # 🔥 REMOVE cluster labels on top
  column_title = NULL,
    
    # Performance
    use_raster = TRUE,     # Faster rendering for large matrices
    raster_by_magick = ifelse(requireNamespace("magick", quietly = TRUE), TRUE, FALSE),
    
    # Row splitting by gene clusters
    row_split = gclust,
    row_gap = unit(2, "mm"),  # Gap between modules
    
    # Annotations
    top_annotation = ha_col,
    
    # Aesthetics
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,  # Column dendrogram is redundant (pre-clustered)
    heatmap_legend_param = list(title = "Expression", title_position = "leftcenter-rot")
  )
  
  # Draw heatmap (suppress ComplexHeatmap messages if desired)
  ht_opt$message = FALSE  # Disable default messages
  draw(ht)
  
  # Saving the heatmap
  png(file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix,
  "heatmap_clusterized.png")),
  width = 3000,
  height = 2500,
  res = 300
)

draw(ht)

dev.off()

# Saving in PDF
pdf(
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "heatmap_clusterized.pdf")
  ),
  width = 12,
  height = 10
)

draw(ht)

dev.off()



library(dplyr)
library(tidyr)

# Getting the percentage of cells belongign to each individual in each cluster: Within each cluster → how cells are distributed across individuals
# Extract metadata
metadata <- sc10x_subset@meta.data

# Count cells per cluster and individual
cluster_composition <- metadata %>%
  group_by(cell_clusters, orig.ident) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  ungroup()

# Calculate percentage of each individual within each cluster
cluster_composition <- cluster_composition %>%
  group_by(cell_clusters) %>%
  mutate(percent = 100 * cell_count / sum(cell_count)) %>%
  select(cell_clusters, orig.ident, percent) %>%
  pivot_wider(
    names_from = orig.ident,
    values_from = percent,
    values_fill = 0
  )

# Save
write.table(
  cluster_composition,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "_cell_cluster_composition.csv")
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Visualzing using a stacked barplot
individual_cluster <- metadata %>%
    count(orig.ident, cell_clusters) %>%
    group_by(orig.ident) %>%
    mutate(percent = n / sum(n))

# Creating colors for clusters:
cluster_colors <- c(
    "1" = "#E69F00", # orange
    "2" = "#88CCEE", # sky blue
    "3" = "#009E73", # green
    "4" = "#F0E442", # yellow
    "5" = "#332288", # deep blue-violet
    "6" =  "#CC79A7", # purple/pink
    "7" =   "#AA4499", # violet
    "8" = "#8DB600", # apple green
    "9" = "#882255" # dark purple
    
)

# Order for patients
patient_order<- c(
    "16_scRNASeq_A83_AS_tp2",
  "04_scRNASeq_A81_AS_tp2",
  "13_scRNASeq_A62_AS_tp2",
  "14_scRNASeq_A75_AS_tp2",
  "15_scRNASeq_A76_S_tp2",
  "12_scRNASeq_A86_S_tp2",
    "08_scRNASeq_A82_S_tp2"  
)

individual_cluster$orig.ident <- factor(
    individual_cluster$orig.ident,
    levels = patient_order
)


    
barplot_cluster = ggplot(
    individual_cluster,
    aes(
        x = orig.ident,
        y = percent,
        fill = factor(cell_clusters)
    )
) +
    geom_col() +
    scale_fill_manual(
    values = cluster_colors) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw(base_size = 18) +
    labs(
        x = "Patient",
        y = "% cells",
        fill = "Cluster"
    )    
    
    ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_stacked_barplot_cluster_compo.png")
    ),
    plot = barplot_cluster,
    width = 8,
    height = 6,
    dpi = 300
)
    
# Same no legend    
    
  barplot_cluster = ggplot(
    individual_cluster,
    aes(
        x = orig.ident,
        y = percent,
        fill = factor(cell_clusters)
    )
) +
    geom_col() +
    scale_fill_manual(
    values = cluster_colors) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw(base_size = 18) +
    labs(
        x = "Patient",
        y = "% cells",
        fill = "Cluster"
    ) +
        theme(legend.position="none")   
    
    ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_stacked_barplot_cluster_compo_no_legend.png")
    ),
    plot = barplot_cluster,
    width = 8,
    height = 6,
    dpi = 300
)







  
    
    
    
    
    
    
    
    
    
    

# Now: Within each individual → how cells are distributed across clusters
# Extract metadata
metadata <- sc10x_subset@meta.data

# Count cells per individual per cluster
cluster_composition <- metadata %>%
  group_by(orig.ident, cell_clusters) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(percent = 100 * cell_count / sum(cell_count)) %>%
  ungroup() %>%
  select(orig.ident, cell_clusters, percent) %>%
  pivot_wider(
    names_from = cell_clusters,
    values_from = percent,
    values_fill = 0
  )

# Save
write.table(
  cluster_composition,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "_cell_individual_composition.csv")
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Visualizing using a stacked barplot
cluster_individual <- metadata %>%
    count(cell_clusters, orig.ident) %>%
    group_by(cell_clusters) %>%
    mutate(percent = n / sum(n))

# Creating colors for patients:
patient_colors <- c(
    "08_scRNASeq_A82_S_tp2" = "#882255", # dark purple
    "12_scRNASeq_A86_S_tp2" = "#AA4499", # violet
    "15_scRNASeq_A76_S_tp2" = "#CC79A7", # purple/pink
    "16_scRNASeq_A83_AS_tp2" = "#8DB600", # apple green
    "04_scRNASeq_A81_AS_tp2" = "#999933", # 14. Olive Green
    "13_scRNASeq_A62_AS_tp2" = "#009E73", # green
    "14_scRNASeq_A75_AS_tp2" = "#117733" # deep green
    
)

cluster_individual$orig.ident <- factor(
    cluster_individual$orig.ident,
    levels = patient_order
)

barplot_individual = ggplot(
    cluster_individual,
    aes(
        x = factor(cell_clusters),
        y = percent,
        fill = orig.ident
    )
) +
    geom_col() +
    scale_fill_manual(values = patient_colors, breaks = patient_order) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw(base_size = 18) +
    labs(
        x = "Cluster",
        y = "% cells",
        fill = "Patient"
    )
ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_stacked_barplot_indvidual_compo.png")
    ),
    plot = barplot_individual,
    width = 8,
    height = 6,
    dpi = 300
)

# same without legend

barplot_individual = ggplot(
    cluster_individual,
    aes(
        x = factor(cell_clusters),
        y = percent,
        fill = orig.ident
    )
) +
    geom_col() +
    scale_fill_manual(values = patient_colors, breaks = patient_order) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw(base_size = 18) +
    labs(
        x = "Cluster",
        y = "% cells",
        fill = "Patient"
    ) +
            theme(legend.position="none")   
ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_stacked_barplot_indvidual_compo_no_legend.png")
    ),
    plot = barplot_individual,
    width = 8,
    height = 6,
    dpi = 300
)




# Combining both individual and cluster information in a bubble plot

# Dataframe

bubble_df <- metadata %>%
    count(orig.ident, cell_clusters, name = "n") %>%
    
    group_by(cell_clusters) %>%
    mutate(
        pct_in_cluster = 100 * n / sum(n)
    ) %>%
    ungroup() %>%
    
    group_by(orig.ident) %>%
    mutate(
        pct_in_patient = 100 * n / sum(n)
    ) %>%
    ungroup()


# Renaming individuals for the plotting
id_map <- c(
  "16_scRNASeq_A83_AS_tp2" = "C1",
  "04_scRNASeq_A81_AS_tp2" = "C2",
  "13_scRNASeq_A62_AS_tp2" = "C3",
  "14_scRNASeq_A75_AS_tp2" = "C4",
  "15_scRNASeq_A76_S_tp2"  = "S1",
  "12_scRNASeq_A86_S_tp2"  = "S2",
  "08_scRNASeq_A82_S_tp2"  = "S3"
)

bubble_df$sample_id <- id_map[
  as.character(bubble_df$orig.ident)
]

bubble_df$sample_id <- factor(
  bubble_df$sample_id,
  levels = c("C1", "C2", "C3", "C4", "S1", "S2", "S3")
)


# Visualization
bubble_plot <- ggplot(
    bubble_df,
    aes(
        x = factor(cell_clusters),
        y = sample_id
    )
) +
    geom_point(
        aes(
            size = pct_in_cluster,   # <-- patient contribution to cluster
            fill = pct_in_patient    # <-- cluster abundance within patient
        ),
        shape = 21,
        colour = "black",
        stroke = 0.2
    ) +
    scale_x_discrete(
        limits = setdiff(
            levels(factor(bubble_df$cell_clusters)),
            "6"
        )
    ) +
    scale_size_continuous(
        name = "% of cluster",
        range = c(2, 15)
    ) +
    scale_fill_viridis_c(
        name = "% of patient"
    ) +
    theme_bw(base_size = 18) +
    labs(
        x = "Cluster",
        y = "Sample"
    ) +
    theme(
        axis.text.x = element_text(
            angle = 45,
            hjust = 1
        )
    )

ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_bubbleplot_patient_cluster_combined_info.png")
    ),
    plot = bubble_plot,
    width = 8,
    height = 6,
    dpi = 300
)
# in PDF
ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(
            outputFilesPrefix,
            "_bubbleplot_patient_cluster_combined_info.pdf"
        )
    ),
    plot = bubble_plot,
    width = 10,
    height = 7,
    dpi = 600
)




















# Getting the percentage of cells belongign to each condition in each cluster
# Extract metadata
metadata <- sc10x_subset@meta.data

# Count cells per cluster and individual
cluster_composition <- metadata %>%
  group_by(cell_clusters, condition) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  ungroup()

# Calculate percentage of each individual within each cluster
cluster_composition <- cluster_composition %>%
  group_by(cell_clusters) %>%
  mutate(percent = 100 * cell_count / sum(cell_count)) %>%
  select(cell_clusters, condition, percent) %>%
  pivot_wider(
    names_from = condition,
    values_from = percent,
    values_fill = 0
  )

# Save
write.table(
  cluster_composition,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "_cell_cluster_composition_condition.csv")
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Visualizing using a stacked barplot
cluster_condition <- metadata %>%
    count(cell_clusters, condition) %>%
    group_by(cell_clusters) %>%
    mutate(percent = n / sum(n))
    
    # condition colors
    condition_colors <- c(
    "Symptomatic" = "red",
    "Chronic" = "orange"
    
)

condition_order<- c(
    "Symptomatic",
    "Chronic"    
)
cluster_condition$condition <- factor(
    cluster_condition$condition,
    levels = condition_order
)

    
barplot_condition = ggplot(
    cluster_condition,
    aes(
        x = factor(cell_clusters),
        y = percent,
        fill = condition
    )
) +
    geom_col() +
     scale_fill_manual(values = condition_colors) +
    scale_y_continuous(labels = scales::percent) +
    theme_bw(base_size = 18) +
    labs(
        x = "Cluster",
        y = "% cells",
        fill = "Condition"
    ) +
    theme(legend.position="none")
    

ggsave(
    filename = file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_stacked_barplot_condition_compo.png")
    ),
    plot = barplot_condition,
    width = 8,
    height = 6,
    dpi = 300
)    
    
    
    
    
        




## Putting back clusters identity from the unsupervised clustrering in the SciGeneX object

sc10x_subset <- AddMetaData(
  sc10x_subset,
  metadata = as.factor(cell_clusters),
  col.name = "heatmap_clusters"
)

saveRDS(sc10x_subset, file = file.path(PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "_seuratObject_unsupervised_clustering.RDS"))
  )
