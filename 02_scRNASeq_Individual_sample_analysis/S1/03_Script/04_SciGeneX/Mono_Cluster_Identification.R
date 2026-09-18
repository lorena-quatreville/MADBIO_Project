


     
########################## To rename clusters on the Seurat object created after the SciGeneX analysis ##########################


# to make sure the cell type is defined
cellType = "mono"

# Define a named vector for clusters ("standardized")
cluster_renames <- c("1" = "Classical monocyte", "2" = "Non-classical monocyte", "3"="Classical monocyte", "4" = "Classical monocyte", "5" = "Intermediate monocyte_T cell_cytotoxic doublet", "6" = "Non-classical monocyte", "7" = "Non-classical monocyte_T cell_B cell_triplet", "8" = "Non-classical monocyte_T cell_cytotoxic doublet", "9" = "Acute inflammatory cell", "10"="T cell", "11" = "Intermediate monocyte_platelet_doublet")

# loading the Seurat object with the new clustering based on the SciGeneX analysis
sc10x <- readRDS(file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "_seuratObject_after_SciGeneX.RDS")))

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
 ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "DimPlot_clusters_rename.png")), plot = umap_clusters_rename, device ="png", width = 16, height = 10, dpi = 300)

# save the seurat object 
saveRDS(sc10x, file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "_seuratObject_after_SciGeneX.RDS")))
  
  
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
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "barplot_cell_proportion.png")), plot = plot_cell_proportion, device ="png",width=16, height = 8, dpi = 300)


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
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "stacked_barplot_cell_proportion.png")), plot = stacked_barplot, device ="png", width=12, height = 10, dpi = 300)


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
write.csv(cluster_summary, file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType,"cluster_summary.csv")), row.names = FALSE)

cluster_detailed <- c("1" = "Classical monocyte_structural regulation signature", "2" = "Non-classical monocyte_homeostatic signature", "3"="Classical monocyte_inflammatory acute activation", "4" = "Classical monocyte_structural regulation signature_2", "5" = "Intermediate monocyte_T cell_cytotoxic doublet", "6" = "Non-classical monocyte_homeostatic signature_2", "7" = "Non-classical monocyte_T cell_B cell_triplet", "8" = "Non-classical monocyte_T cell_cytotoxic doublet", "9" = "Acutely inflammatory cell", "10"="T cell", "11" = "Intermediate monocyte_platelet_doublet")


