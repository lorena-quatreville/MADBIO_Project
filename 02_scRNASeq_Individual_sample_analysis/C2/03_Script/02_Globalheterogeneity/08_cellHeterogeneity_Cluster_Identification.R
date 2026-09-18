

### Cluster labeling

# Load the seurat object
seuratObjectPath = file.path( PATH_ANALYSIS_OUTPUT, "Resolution_0.6", paste0( outputFilesPrefix, "seuratObject_final.RDS"))
sc10x <- readRDS(seuratObjectPath)

# Step 1: Define a named vector for clusters ("standardized")
cluster_renames <- c("1" = "Naive CD4 T cell", "2" = "Effector CD4 T cell", "3" = "Class-switched B cell", "4" = "Memory CD8 T cell", "5" = "Naive CD8 T cell", "6" = "Naive B cell", "7" = "Classical monocyte", "8" = "Effector CD8 T cell", "9" = "NK", "10" = "Transitional B cell", "11" = "Doublet", "12" = "Unconventional T cell_γδ T", "13" = "High metabolic cell", "14" = "Non-classical monocyte", "15" = "Non-classical monocyte", "16"= "Non-classical monocyte", "17" = "DN T cell", "18" = "Dendritic cell", "19" = "Plasmacitoid dendritic cell")


#Step 2: Add cluster names to metadata
sc10x$cluster_renames <- factor(sc10x$seurat_clusters, levels = names(cluster_renames), labels = cluster_renames)


# Function to break label into multiple lines if it's too long
break_label_lines <- function(label, max_chars = 20) {
  words <- strsplit(label, " ")[[1]]  # Split the label by spaces into words
  line <- ""
  lines <- c()  # Initialize empty vector for lines
  
  for (word in words) {
    # Check if adding the word exceeds max_chars
    if (nchar(paste(line, word, sep = " ")) <= max_chars) {
      line <- paste(line, word)
    } else {
      lines <- c(lines, line)  # Add current line to lines
      line <- word  # Start a new line with the current word
    }
  }
  lines <- c(lines, line)  # Add the last line
  return(paste(lines, collapse = "\n"))  # Return the lines as a single string
}

# Step 2: Convert the cluster names to a character vector and break the lines
sc10x$cluster_renames <- as.character(sc10x$cluster_renames)  # Ensure it's a character vector
sc10x$cluster_renames <- sapply(sc10x$cluster_renames, break_label_lines)

# Compute centroids for each cluster
centroids <- sc10x@reductions$umap@cell.embeddings %>%
  as.data.frame() %>%
  mutate(cluster = sc10x$cluster_renames) %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2))



# Step 3: Plot UMAP with cluster names in the legend
umap_clusters_rename <- DimPlot(sc10x, group.by = "cluster_renames", label = FALSE) +
  ggtitle("UMAP with Cluster Names") +
  geom_text_repel(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = cluster), box.padding = 0.5, point.padding = 0.5, size = 4, max.overlaps = 100)+
  theme(legend.title = element_blank())
  print(umap_clusters_rename)
 ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, "Resolution_0.6", paste0( outputFilesPrefix, "DimPlot_clusters_rename.png")), plot = umap_clusters_rename, device ="png", width = 23, height = 13, units ="cm")

# save the new seurat object 
# Save a binary file of final Seurat object only (as RDS)
saveRDS( object = sc10x, file = seuratObjectPath);
  
  
######################################### Plot cell proportion per cluster with an historgram (ranked from the biggest cluster to the lowest)

## remove line breaks

# Remove accidental line breaks in cluster names
sc10x$cluster_renames <- gsub("\n", " ", sc10x$cluster_renames)
sc10x$cluster_renames <- trimws(sc10x$cluster_renames)  # Remove leading/trailing spaces
sc10x$cluster_renames <- gsub("\\s+", " ", sc10x$cluster_renames)  # Replace multiple spaces with single space
sc10x$cluster_renames <- as.factor(sc10x$cluster_renames)  # Ensure it's a factor

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
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, "Resolution_0.6", paste0( outputFilesPrefix, "barplot_cell_proportion.png")), plot = plot_cell_proportion, device ="png",width=16, height = 8, dpi = 300)


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
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, "Resolution_0.6", paste0( outputFilesPrefix, "stacked_barplot_cell_proportion.png")), plot = stacked_barplot, device ="png", width=12, height = 10, dpi = 300)


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


###### UMAP and figures for the paper


# Loading the seurat object
sc10x = readRDS(seuratObjectPath)

## remove line breaks

# Remove accidental line breaks in cluster names
sc10x$cluster_renames <- gsub("\n", " ", sc10x$cluster_renames)
sc10x$cluster_renames <- trimws(sc10x$cluster_renames)  # Remove leading/trailing spaces
sc10x$cluster_renames <- gsub("\\s+", " ", sc10x$cluster_renames)  # Replace multiple spaces with single space

# Compute a umap with cluster numbers on the graph and in the lgend + cluster labels in the legend

palette_19 <- c(
  "#E69F00", # orange
  "#88CCEE", # sky blue
  "#009E73", # green
  "#F0E442", # yellow
  "#332288", # deep blue-violet
  "#D55E00", # vermillion
  "#CC79A7", # purple/pink
  "#AA4499", # violet
  "#8DB600", # apple green
  "#882255", # dark purple
  "#117733", # deep green
  "#666666" , # dark gray
  "#44AA99", # 13. Teal
  "#999933", # 14. Olive Green
  "#DDCC77", # 15. Pale Yellow/Gold
  "#BBCCDA", # 16. Light Blue-Gray
  "#E1BE6A", # 17. Sepia/Tan
  "#4477AA", # 18. Medium Blue
  "#77AADD"  # 19. Light Lavender/Periwinkle
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
    cols = palette_19,
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
    size = 5,
    family = "Arial",
    fontface = "bold",   
    color = "black",            # Add an explicit color for visibility
    bg.color = "white",
    bg.r = 0.15,
    min.segment.length = 0
  ) +
  theme(
  text = element_text(family = "Arial"),
    # Axis titles ("UMAP_1", "UMAP_2")
  axis.title = element_text(
    size = 20
  ),

  # Axis tick labels (the numbers)
  axis.text = element_text(
    size = 20)
  ) +
  labs(title = NULL)
  
# Saving the UMAP
ggsave(plot = umap_cluster_label, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_umap_cluster_label_nbr.pdf"))
       , width = 30, height = 16, units ="cm", dpi = 600, device = cairo_pdf)
ggsave(plot = umap_cluster_label, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_umap_cluster_label_nbr.png"))
       , width = 30, height = 16, units ="cm", dpi = 600)
       
# Same without the legend for the paper

umap_cluster_label_no_legend <- DimPlot(
    sc10x,
    reduction = ifelse(exists("useReduction"), useReduction, "umap"),
    group.by = "cluster_number_label",
    cols = palette_19,
    label = FALSE
) +
  ggrepel::geom_text_repel(
    data = cluster_centers,
    aes(
      x = UMAP_1,
      y = UMAP_2,
      label = seurat_clusters
    ),
    inherit.aes = FALSE,
    size = 8,
    family = "Arial",
    fontface = "bold",
    color = "black",
    bg.color = "white",
    bg.r = 0.15,
    min.segment.length = 0
  ) +
  theme(
    legend.position = "none",
    text = element_text(family = "Arial"),
    # Axis titles ("UMAP_1", "UMAP_2")
  axis.title = element_text(
    size = 20
  ),

  # Axis tick labels (the numbers)
  axis.text = element_text(
    size = 20
  )
  ) +
  labs(title = NULL)
  
# Saving the UMAP
ggsave(plot = umap_cluster_label_no_legend, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_umap_cluster_label_nbr_no_legend.pdf")),width = 10, height = 8, dpi = 600, device = cairo_pdf)
ggsave(plot = umap_cluster_label_no_legend, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_umap_cluster_label_nbr_no_legend.png")),width = 10, height = 8, dpi = 600)
       
       
       
## remove line breaks

# Remove accidental line breaks in cluster names
sc10x$cluster_renames <- gsub("\n", " ", sc10x$cluster_renames)
sc10x$cluster_renames <- trimws(sc10x$cluster_renames)  # Remove leading/trailing spaces
sc10x$cluster_renames <- gsub("\\s+", " ", sc10x$cluster_renames)  # Replace multiple spaces with single space
sc10x$cluster_renames <- as.factor(sc10x$cluster_renames)  # Ensure it's a factor

# Create a dataframe combining cluster name and number
cluster_info <- data.frame(
  Cluster_Name = sc10x$cluster_renames,
  Cluster_Number = sc10x$seurat_clusters
)

# Calculate absolute counts and proportions
cluster_summary <- cluster_info %>%
  group_by(Cluster_Name, Cluster_Number) %>%
  summarise(Cell_Count = n(), .groups = "drop") %>%
  mutate(
    Total_Cells = sum(Cell_Count),
    Proportion = Cell_Count / Total_Cells
  ) %>%
  select(Cluster_Name, Cluster_Number, Cell_Count, Proportion, Total_Cells)


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
write.csv(cluster_summary, file = file.path(PATH_ANALYSIS_OUTPUT, "Resolution_0.6", paste0(outputFilesPrefix, "cluster_summary.csv")), row.names = FALSE)




# Violin plot for a couple of markers per cluster
library(RColorBrewer)

# 1. Define immune markers with cell type annotation
marker_list <- list(
  "MS4A1" = "B cells", "CD79A" = "B cells", "CD79B" = "B cells", "CD19" = "B cells",
  "CD3D" = "T cells", "CD3E" = "T cells", "CD3G" = "T cells", "CD4" = "T cells", "CD8A" = "T cells",
  "NKG7" = "NK cells", "GNLY" = "NK cells", "KLRD1" = "NK cells", "FCGR3A" = "NK cells",
  "LYZ" = "Monocytes", "CD14" = "Monocytes", "S100A9" = "Monocytes"
)

marker_genes <- names(marker_list)
marker_groups <- unlist(marker_list)

# Breaking lines to increase the size
sc10x$cluster_renames <- as.character(sc10x$cluster_renames)  # Ensure it's a character vector
sc10x$cluster_renames <- sapply(sc10x$cluster_renames, break_label_lines)

# 2. Violin plot (stacked + flipped) for all marker genes by cluster
vln_plot <- VlnPlot(
  sc10x,
  features = marker_genes,
  group.by = "cluster_number_label",
  pt.size = 0.1,
  stack = TRUE,
  flip = TRUE
) +
  theme_minimal(base_size = 18) +
  labs(
    title = NULL,
    x = NULL
  ) +
  theme(
    legend.position = "none",
    text = element_text(family = "Arial"),

  # Keep expression scale numbers (left side)
  axis.text.y = element_text(
    family = "Arial",
    size = 46
  ),

  # Remove cluster labels at the bottom
  axis.text.x = element_blank(),

  # Remove gene names
  strip.text.y = element_blank(),

  # Optional: remove tick marks if you don't want them
  axis.ticks.y = element_blank(),

  axis.title = element_blank(),

  panel.background = element_rect(fill = "white"),
  plot.background = element_rect(fill = "white")
  )
    
 ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, "Resolution_0.6", paste0(outputFilesPrefix, "VlnPlot_ImmuneMarkers.pdf")),
  plot = vln_plot,
  width = 30,
  height = 20,
  bg = "white",
  device = cairo_pdf,
  dpi = 600
)   
    
 ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, "Resolution_0.6", paste0(outputFilesPrefix, "VlnPlot_ImmuneMarkers.png")),
  plot = vln_plot,
  width = 30,
  height = 20,
  bg = "white",
  dpi = 600
)   
        


cluster_detailed <- c("1" = "Naive CD4 T cell", "2" = "Effector CD4 T cell", "3" = "Class-switched memory B cell", "4" = "Memory CD8 T cell_central", "5" = "Naive CD8 T cell", "6" = "Naive B cell", "7" = "Classical monocyte", "8" = "Effector CD8 T cell", "9" = "NK", "10" = "Transitional B cell", "11" = "Doublet", "12" = "γδ T cell", "13" = "High metabolic cell", "14" = "Non-classical monocyte", "15" = "Non-classical monocyte_2", "16"= "Non-classical monocyte_3", "17" = "DN T cell_exhausted", "18" = "Dendritic cell", "19" = "Plasmacitoid dendritic cell")
