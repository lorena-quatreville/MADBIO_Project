


     
########################## To rename clusters on the Seurat object created after the SciGeneX analysis ##########################


# to make sure the cell type is defined
cellType = "effector_CD8_T_cell_NK"

# Define a named vector for clusters ("standardized")
cluster_renames <- c("1" = "CD8 T cell", "2" = "Unconventional T cell_γδ T", "3"="CD56 dim NK", "4" = "DN T cell_memory", "5" = "CD8 T cell", "6" = "CD56 dim NK", "7" = "DN T cell", "8" = "Uncharacterised")


# loading the Seurat object with the new clustering based on the SciGeneX analysis
sc10x <- readRDS(PATH_SCIGENEX_RDS_FILE)

# To save active.ident to "seurat_clusters" metadat
sc10x[["seurat_clusters"]] <- Idents(sc10x)

#Step 2: Add cluster names to metadata
sc10x$cluster_renames <- factor(sc10x$seurat_clusters, levels = names(cluster_renames), labels = cluster_renames)


# Convert UMAP embeddings to a data frame
umap_emb <- as.data.frame(sc10x@reductions$umap@cell.embeddings)

# Add the cluster labels
umap_emb$cluster <- as.character(sc10x$cluster_renames)  # Ensure it's character

# Compute centroids using aggregate
centroids <- aggregate(
  umap_emb[, c("UMAP_1", "UMAP_2")],     # Columns to average
  by = list(cluster = umap_emb$cluster), # Group by cluster
  FUN = mean,                           # Compute mean
  na.rm = TRUE                           # Ignore NAs if any
)

# Check the result
print(centroids)


# Step 3: Plot UMAP with cluster names in the legend
umap_clusters_rename <- DimPlot(sc10x, group.by = "cluster_renames", label = FALSE) +
  ggtitle("UMAP with Cluster Names") +
  geom_text_repel(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = str_wrap(cluster, width = 20)), box.padding = 0.5, point.padding = 0.5, size = 4, max.overlaps = 100)+
  theme(legend.title = element_blank())
 ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT,  paste0(outputFilesPrefix,  "DimPlot_clusters_rename.png")), plot = umap_clusters_rename, device ="png", width = 16, height = 10, dpi = 300)

# save the seurat object 
saveRDS(sc10x, file = file.path(PATH_SCIGENEX_RDS_FILE))
  

######################### UMAP for the paper

# Remove accidental line breaks in cluster names
sc10x$cluster_renames <- gsub("\n", " ", sc10x$cluster_renames)
sc10x$cluster_renames <- trimws(sc10x$cluster_renames)  # Remove leading/trailing spaces
sc10x$cluster_renames <- gsub("\\s+", " ", sc10x$cluster_renames)  # Replace multiple spaces with single space

# Compute a umap with cluster numbers on the graph and in the lgend + cluster labels in the legend

palette_10 <- c(
  "#E69F00", # orange
  "#88CCEE", # sky blue
  "#F0E442", # yellow
  "#332288", # deep blue-violet
  "#CC79A7", # purple/pink
  "#8DB600", # apple green
  "#882255", # dark purple
  "#117733", # deep green
  "#666666", # dark gray
  "#BBCCDA", # 16. Light Blue-Gray
  "#D55E00" # vermillion
)



# Ensure the combined column for the legend is created
sc10x$cluster_number_label <- paste0(sc10x$seurat_clusters, " - ", sc10x$cluster_renames)

# Sort cluster labels+number by number:
# 1. Create the combined descriptive label (still a character vector)
sc10x$cluster_number_label <- paste0(sc10x$seurat_clusters, " - ", sc10x$cluster_renames)

# 2. Get the unique combined labels in their current alphabetical order
unsorted_levels <- unique(sc10x$cluster_number_label)

# 3. Create a small data frame to match and sort
# This matches the combined label with the original numeric ID for sorting
level_df <- data.frame(
    original_id = as.numeric(gsub(" - .*", "", unsorted_levels)), # Extracts the cluster number
    combined_label = unsorted_levels
)

# 4. CRITICAL: Sort the combined labels based on the numeric cluster ID
sorted_combined_levels <- level_df[order(level_df$original_id), "combined_label"]

# 5. Set the levels of the combined factor using the sorted vector
sc10x$cluster_number_label <- factor(sc10x$cluster_number_label, levels = sorted_combined_levels)



# Plot the UMAP
umap_cluster_label <- DimPlot(
    sc10x,
    reduction = ifelse(exists("useReduction"), useReduction, "umap"),
    group.by = "cluster_number_label",
    cols = palette_10,
    label = FALSE
) + 
    ggtitle("UMAP with cluster identification")

# Merge the raw cluster numbers into the plot data frame
umap_cluster_label$data$seurat_clusters <- sc10x$seurat_clusters[
    match(rownames(umap_cluster_label$data), rownames(sc10x@meta.data))
]

# Calculate the median UMAP coordinates for each seurat_cluster
cluster_centers <- umap_cluster_label$data %>%
  dplyr::group_by(seurat_clusters) %>%
  dplyr::summarise(
    UMAP_1 = median(UMAP_1), # Column name confirmed here
    UMAP_2 = median(UMAP_2)  # Column name confirmed here
  )


# Add the labels using the new centers data frame
umap_cluster_label <- umap_cluster_label +
  # We use ggrepel::geom_text_repel for clarity and robustness
  ggrepel::geom_text_repel(
    data = cluster_centers,
    mapping = aes(
      x = UMAP_1,               # Confirmed column name
      y = UMAP_2,               # Confirmed column name
      label = seurat_clusters
    ), 
    inherit.aes = FALSE,        # Stops inheriting default aesthetics from DimPlot
    size = 4,
    color = "black",            # Add an explicit color for visibility
    bg.color = "white",
    bg.r = 0.15,
    min.segment.length = 0
  )
  
# Saving the UMAP
ggsave(plot = umap_cluster_label, filename = file.path(PATH_ANALYSIS_OUTPUT,  paste0(outputFilesPrefix, "_umap_cluster_label_nbr.png"))
       , width = 30, height = 16, units ="cm")


################################ Heatmap for the paper

# Get back the cluster set from SciGeneX
res_subtype_SciGeneX <- readRDS(file = file.path(PATH_ANALYSIS_OUTPUT,  paste0(outputFilesPrefix, cellType,  "_clusterset_object_after_SciGeneX.RDS")))
# Reconstruct sorted cluster labels (same as UMAP)
sorted_combined_levels <- levels(sc10x$cluster_number_label)

# Build a named color vector
cluster_colors <- palette_10
names(cluster_colors) <- sorted_combined_levels

# Forces the same factor order
Idents(sc10x) <- sc10x$seurat_clusters


heatmap_one_module <- function(module_nb){

    # Ensure consistent label set between Seurat object and color mapping
  cell_clusters <- Idents(sc10x)
  
  # Heatmap with top genes (non-interactive)
  p_heatmap <- plot_ggheatmap(
  res_subtype_SciGeneX[module_nb,], 
  ident = cell_clusters, 
  use_top_genes = T,
  color_ident = cluster_colors
  ) +
  ggtitle("All clusters (top genes)") +  theme(strip.text.y = element_text(size=4))
  
  # save plot
  ggsave(plot = p_heatmap, filename = file.path(PATH_ANALYSIS_OUTPUT,   paste0("meaningful_modules_", outputFilesPrefix,SCIGENEX_METHOD,"_inflation_",SCIGENEX_INFLATION,"_heatmap_after_SciGeneX_paper.png"))
       , width = 23, height = 13, units ="cm")
}






######################################### Plot cell proportion per cluster with an historgram (ranked from the biggest cluster to the lowest)


# Create a table of counts by cluster and sample
cluster_counts <- table(sc10x$cluster_renames)

# Calculate proportions within each sample
cluster_proportions <- prop.table(cluster_counts)

# Convert to a data frame for ggplot
cluster_df <- as.data.frame(cluster_proportions)
colnames(cluster_df) <- c("Cluster", "Proportion")

library(ggplot2)

# Reorder the factors based on proportions (from highest to lowest)
cluster_df$Cluster <- factor(cluster_df$Cluster, levels = cluster_df$Cluster[order(cluster_df$Proportion, decreasing = TRUE)])


# Create the grouped bar plot
plot_cell_proportion <- ggplot(cluster_df, aes(x = Cluster, y = Proportion, fill= Cluster)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Proportion of cells by cell type",
       x = "Clusters",
       y = "Proportion") +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::percent_format()) + # Show proportions as percentages
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot_cell_proportion)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT,  paste0(outputFilesPrefix,  "barplot_cell_proportion.png")), plot = plot_cell_proportion, device ="png",width=16, height = 8, dpi = 300)


# same with a stacked bar plot

#creation du tableau pour barplot et supprimé les lignes Negative et Doublet
cell_by_cluster <- data.frame(table(sc10x$cluster_renames))

# Calculer le total des fréquences pour chaque groupe 
cell_by_cluster <- cell_by_cluster %>%
  mutate(Total_Freq = sum(Freq)) %>%
  ungroup()

# Calculer le pourcentage de chaque Freq
cell_by_cluster <- cell_by_cluster %>%
  mutate(Percentage = (Freq / Total_Freq) * 100)

# Rename columns for clarity
colnames(cell_by_cluster) <- c("Cluster", "Frequency", "Total_Frequency", "Percentage")


stacked_barplot <- ggplot(cell_by_cluster, aes(fill = Cluster, y = Frequency, x = SAMPLE_NAME, 
                                 text = paste0("Cluster: ", Cluster, 
                                               "<br>Number of cells: ", Frequency, 
                                               "<br>Percentage: ", round(Percentage, 1), "%")))+ 
  geom_bar(width = 0.5, position = "fill", stat = "identity", color = "white") +
  scale_x_discrete(labels = SAMPLE_NAME) +
  xlab("") +
  ylab("Frequency") +
  labs(fill = "Cluster") +
  theme_classic(base_size = 15) +
  theme(legend.text = element_text(size = 18),  # Increase legend text size
        legend.title = element_text(size = 16)) + # Increase legend title size
  theme(axis.text.x = element_text(color = "black", face = "bold"))

stacked_barplot <- stacked_barplot + theme(plot.margin = margin(10, 10, 10, 10)) # Adjust margins



# Save the graph
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT,  paste0(outputFilesPrefix,  "stacked_barplot_cell_proportion.png")), plot = stacked_barplot, device ="png", width=12, height = 10, dpi = 300)


# Rendre le graphique interactif avec plotly
p_interactive <- ggplotly(stacked_barplot, tooltip = "text")

# Afficher le graphique interactif
p_interactive <- config(p_interactive, displaylogo = FALSE,
                 toImageButtonOptions = list( format='svg'),
                 modeBarButtons = list( list('toImage'),
                                        list( 'zoom2d', 'pan2d', 'resetScale2d')))

browsable(
  div( style = paste("display: flex; align-items: center; justify-content: center;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: green;"),
       div( style = paste("display: flex; flex-flow: row nowrap; overflow-x: auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: red;"),
            div(p_interactive, style = paste("flex : 0 0 auto;", if(.SHOWFLEXBORDERS) "border-style: solid; border-color: blue;"))))
)



### To save the cell proportions in a CSV file


# Create a dataframe combining cluster name and number
cluster_info <- data.frame(
  Cluster_Name = sc10x$cluster_renames,
  Cluster_Number = sc10x$seurat_clusters
)

cluster_summary <- cluster_info %>%
  dplyr::group_by(Cluster_Name, Cluster_Number) %>%
  dplyr::summarise(Cell_Count = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(
    Total_Cells = sum(Cell_Count),
    Proportion = Cell_Count / Total_Cells
  ) %>%
  dplyr::select(Cluster_Name, Cluster_Number, Cell_Count, Proportion, Total_Cells)



# Add a "Total" row
total_row <- data.frame(
  Cluster_Name = "Total",
  Cluster_Number = NA,
  Cell_Count = sum(cluster_summary$Cell_Count),
  Proportion = sum(cluster_summary$Proportion),
  Total_Cells = unique(cluster_summary$Total_Cells)
)

# Combine the data
cluster_summary <- bind_rows(cluster_summary, total_row)

# Save as CSV
write.csv(cluster_summary, file = file.path(PATH_ANALYSIS_OUTPUT,  paste0(outputFilesPrefix, "cluster_summary.csv")), row.names = FALSE)

cluster_detailed <- c("1" = "Effector CD8", "2" = "γδ T cell", "3"="CD56 dim NK1_cytotoxic", "4" = "Effector memory CD8_DN T cell", "5" = "Effector CD8_immediate ativation", "6" = "CD56 dim NK_recently activated", "7" = "DN T cell", "8" = "No mod")
