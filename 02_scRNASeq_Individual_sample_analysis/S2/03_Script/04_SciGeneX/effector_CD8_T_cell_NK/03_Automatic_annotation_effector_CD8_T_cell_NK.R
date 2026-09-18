library(SeuratData)
library(Seurat)
library(Azimuth)
library(ggplot2)
library(dplyr)
library(scales)
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(pheatmap)
library(SingleR)
library(tibble)
library(patchwork)
library(purrr) 

# Loading the SeuratObject
sc10x <- readRDS(PATH_SCIGENEX_RDS_FILE)

# Loading the reference data from Azimuth
PBMC_ref <- LoadReference(PBMC_REF_PATH)


# Setting RNA assay
DefaultAssay(sc10x) <- "RNA"

# Creating a temporary seurat object or SCTransform normalization
sc10x_SCT <- sc10x

# Applying SCTransform on the Seurat object
sc10x_SCT <- SCTransform(
  sc10x_SCT,
  verbose = TRUE
)

DefaultAssay(sc10x_SCT) <- "SCT"

# Extracting feature names from Azimuth reference dataset
ref_sct_features <- rownames(
  PBMC_ref$map@reductions$refDR@feature.loadings
)


# Finding anchors between reference and query datasets

anchors <- FindTransferAnchors(
  reference = PBMC_ref$map,
  query = sc10x_SCT,
  normalization.method = "SCT",
  reference.reduction = "refDR",
  dims = 1:30,
  features = ref_sct_features
)

# Label transfer

sc10x_SCT <- MapQuery(
  anchorset = anchors,
  query = sc10x_SCT,
  reference = PBMC_ref$map,
  refdata = list(
    celltype.l1 = "celltype.l1",
    celltype.l2 = "celltype.l2",
    celltype.l3 = "celltype.l3"
  ),
  reference.reduction = "refDR",
  reduction.model = "refUMAP"
)

# Compute Azimuth confidence score
sc10x_SCT$predicted.celltype.best <- ifelse(sc10x_SCT$predicted.celltype.l2.score >= 0.5 & sc10x_SCT$predicted.celltype.l3.score >= 0.6,
sc10x_SCT$predicted.celltype.l3,sc10x_SCT$predicted.celltype.l2)



# Creating cluster purity tables
annotation_levels <- c(
  "predicted.celltype.l1",
  "predicted.celltype.l2",
  "predicted.celltype.l3",
  "predicted.celltype.best")

cluster_annotation_tables <- lapply(annotation_levels, function(lvl) {
  prop.table(
    table(sc10x_SCT$cluster_renames,
          sc10x_SCT@meta.data[[lvl]]),
    margin = 1
  )
})

names(cluster_annotation_tables) <- annotation_levels

# To print them?
cluster_annotation_tables$predicted.celltype.l1
cluster_annotation_tables$predicted.celltype.l2
cluster_annotation_tables$predicted.celltype.l3
cluster_annotation_tables$predicted.celltype.best

# heatmap
heatmap_results <- data.frame()

clusters <- sort(unique(sc10x_SCT$seurat_clusters))

for (cl in clusters) {
  
  cells_cl <- sc10x_SCT@meta.data %>%
    filter(seurat_clusters == cl)
  
  for (lvl in annotation_levels) {
    
    if (!lvl %in% colnames(cells_cl)) next
    
    counts <- prop.table(table(cells_cl[[lvl]]))
    major_type <- names(counts)[which.max(counts)]
    major_score <- max(counts)
    
    heatmap_results <- rbind(
      heatmap_results,
      data.frame(
        cluster = cl,
        level = lvl,
        cell_type = major_type,
        score = major_score,
        stringsAsFactors = FALSE
      )
    )
  }
}

heatmap_results$cluster <- factor(
  heatmap_results$cluster,
  levels = clusters
)

plot_heatmap <- ggplot(
  heatmap_results,
  aes(x = level, y = cluster, fill = score)
) +
  geom_tile(color = "white") +
  geom_text(aes(label = cell_type), size = 3) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_bw(base_size = 12) +
  labs(
    title = "Azimuth annotation per Seurat cluster",
    x = "Annotation level",
    y = "Cluster",
    fill = "Proportion"
  )

ggsave(
  file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_Azimuth_heatmap.png")),
  plot = plot_heatmap,
  width = 12,
  height = 6,
  dpi = 300
)

# Plotting best levels
n_colors <- length(unique(sc10x_SCT$predicted.celltype.best))
colors <- hue_pal()(n_colors)

plot_best <- DimPlot(
  sc10x_SCT,
  reduction = "umap",
  group.by = "predicted.celltype.best",
  label = TRUE,
  repel = TRUE,
  label.size = 4
) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(title = "UMAP – Azimuth best") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_UMAP_best.png")),
  plot = plot_best,
  width = 10,
  height = 8,
  dpi = 300
)


# Plotting L1 level
n_colors <- length(unique(sc10x_SCT$predicted.celltype.l1))
colors <- hue_pal()(n_colors)

plot_l1 <- DimPlot(
  sc10x_SCT,
  reduction = "umap",
  group.by = "predicted.celltype.l1",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(title = "UMAP – Azimuth Level 1") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave(
  file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_UMAP_L1.png")),
  plot = plot_l1,
  width = 10,
  height = 8,
  dpi = 300
)


# Plotting L2 level
n_colors <- length(unique(sc10x_SCT$predicted.celltype.l2))
colors <- hue_pal()(n_colors)

plot_l2 <- DimPlot(
  sc10x_SCT,
  reduction = "umap",
  group.by = "predicted.celltype.l2",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(title = "UMAP – Azimuth Level 2") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave(
  file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_UMAP_L2.png")),
  plot = plot_l2,
  width = 10,
  height = 8,
  dpi = 300
)

# Plotting L3 level
n_colors <- length(unique(sc10x_SCT$predicted.celltype.l3))
colors <- hue_pal()(n_colors)

plot_l3 <- DimPlot(
  sc10x_SCT,
  reduction = "umap",
  group.by = "predicted.celltype.l3",
  label = TRUE,
  repel = TRUE,
  label.size = 5
) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(title = "UMAP – Azimuth Level 3") +
  theme(legend.position = "none", plot.title = element_text(hjust = 0.5))

ggsave(
  file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_UMAP_L3.png")),
  plot = plot_l3,
  width = 10,
  height = 8,
  dpi = 300
)


# UMAP with the most abundant cell group label (from predicted.celltype.best)
# Find dominant predicted.celltype.best per cluster
cluster_dominant_labels <- sc10x_SCT@meta.data %>%
  group_by(seurat_clusters, predicted.celltype.best) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(seurat_clusters) %>%
  slice_max(order_by = n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(seurat_clusters, dominant_celltype = predicted.celltype.best)
  
# Adding the cluster-level label back to all cells
sc10x_SCT@meta.data <- sc10x_SCT@meta.data %>%
  left_join(cluster_dominant_labels, by = "seurat_clusters")


# Ensure Idents match Seurat clusters
Idents(sc10x_SCT) <- sc10x_SCT$seurat_clusters

# Extracting UMAP embeddings and metadata
umap_df <- as.data.frame(sc10x_SCT@reductions$umap@cell.embeddings)
umap_df$cluster <- sc10x_SCT$seurat_clusters
umap_df$dominant <- sc10x_SCT$dominant_celltype

# Computing cluster centroids for labeling
centroids <- umap_df %>%
  group_by(cluster, dominant) %>%
  summarise(
    UMAP_1 = mean(UMAP_1),
    UMAP_2 = mean(UMAP_2),
    .groups = "drop"
  )

# Defining colors for clusters
clusters <- sort(unique(umap_df$cluster))
colors <- hue_pal()(length(clusters))
names(colors) <- clusters

# Plotting
plot_dominant_label <- ggplot(umap_df, aes(x = UMAP_1, y = UMAP_2, color = cluster)) +
  geom_point(size = 0.5) +
  geom_text_repel(
    data = centroids,
    aes(x = UMAP_1, y = UMAP_2, label = dominant),  # label = dominant cell type
    size = 3,
    color = "black",
    box.padding = 0.5,
    point.padding = 0.5
  ) +
  scale_color_manual(values = colors) +
  theme_minimal(base_size = 14) +
  labs(title = "UMAP – Seurat clusters colored, labeled by dominant cell type") +
  theme(legend.position = "right",
  panel.background = element_rect(fill = "white", color = NA), # white background
    plot.background = element_rect(fill = "white", color = NA))

plot_dominant_label

ggsave(
  file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_UMAP_best_dominant_cluster.png")),
  plot = plot_dominant_label,
  width = 10,
  height = 8,
  dpi = 300
)




# Creating a summary table
# Build the table
cluster_celltype_table <- sc10x_SCT@meta.data %>%
  select(
    seurat_clusters,
    cluster_renames,
    predicted.celltype.best
  ) %>%
  group_by(
    seurat_clusters,
    cluster_renames,
    predicted.celltype.best
  ) %>%
  summarise(
    n_cells = n(),
    .groups = "drop"
  ) %>%
  arrange(seurat_clusters, desc(n_cells))

# View table
print(cluster_celltype_table)

write.csv(
  cluster_celltype_table,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "_cluster_predicted_celltype_best.csv")
  ),
  row.names = FALSE
)

# Save a binary file of final Seurat object only (as RDS)
seuratObjectPath = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, cellType,  "_seuratObject_SCT_Azimuth.RDS"));
saveRDS( object = sc10x_SCT, file = seuratObjectPath)
