######### Data merging and normalization


# READ DATA
# ---------

## @knitr loadData

# Increase power
options(future.globals.maxSize = 16 * 1024^3)  # 16GB if needed
getOption("future.globals.maxSize")

# Retreive monocytes subset and merge them in a new Seurat object

sc10x_monocytes <- NULL  # Initialize before the loop

for (n in 1:nb_datasets) {
    
    ## Retrieve Seurat object from the first analysis
    PATH_SEURAT_RDS_FILE_n <- PATH_SEURAT_RDS_FILE[[n]]
    sc10x = readRDS(PATH_SEURAT_RDS_FILE_n)
    
    ##  Extracting monocytes subsets

    # Retrieve cluster names
    # Check cluster names
    Idents(sc10x) <- "cluster_renames"
    cluster_names <- levels(Idents(sc10x))
    
    # Retrieve the name of all clusters having "monocytes" in their label
    monocyte_clusters <- grep("monocyte", cluster_names, value = TRUE, ignore.case = TRUE)
    
    # Subset Seurat object with selected clusters
    sc10x_subset <- subset(sc10x, idents = monocyte_clusters)
    
    # Merge the new dataset subset with others
    if (is.null(sc10x_monocytes)) { 
        sc10x_monocytes <- sc10x_subset  # Initialize with first dataset
    } else {
        sc10x_monocytes <- merge(sc10x_monocytes, y = sc10x_subset)
    }


} 

sc10x <- sc10x_monocytes
gc()



############################# Plot cell proportions between samples
# Remove accidental line breaks in cluster names
sc10x$cluster_renames <- gsub("\n", " ", sc10x$cluster_renames)
sc10x$cluster_renames <- trimws(sc10x$cluster_renames)  # Remove leading/trailing spaces
sc10x$cluster_renames <- gsub("\\s+", " ", sc10x$cluster_renames)  # Replace multiple spaces with single space
sc10x$cluster_renames <- as.factor(sc10x$cluster_renames)  # Ensure it's a factor

# Create a table of counts by cluster and sample
cluster_sample_counts <- table(sc10x$cluster_renames, sc10x$orig.ident)

# Calculate proportions within each sample
cluster_sample_proportions <- prop.table(cluster_sample_counts, margin = 2)

# Convert to a data frame for ggplot
cluster_sample_df <- as.data.frame(cluster_sample_proportions)
colnames(cluster_sample_df) <- c("Cluster", "Sample", "Proportion")

library(ggplot2)

# Create the grouped bar plot
plot_cell_proportion <- ggplot(cluster_sample_df, aes(x = Cluster, y = Proportion, fill = Sample)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Proportion of cells by cell type between samples",
       x = "Clusters",
       y = "Proportion") +
  scale_fill_viridis_d() +
  scale_y_continuous(labels = scales::percent_format()) + # Show proportions as percentages
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(plot_cell_proportion)
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "barplot_cell_proportion.png")), plot = plot_cell_proportion, device ="png", width = 12, height = 6, dpi = 300)


############################## Same but using a stacked interactive plot

#creation du tableau pour barplot et supprimé les lignes Negative et Doublet
cell_by_cluster <- data.frame(table(sc10x$cluster_renames, sc10x$orig.ident))

# Calculer le total des fréquences pour chaque groupe 
cell_by_cluster <- cell_by_cluster %>%
  group_by(Var2) %>%
  mutate(Total_Freq = sum(Freq)) %>%
  ungroup()

# Calculer le pourcentage de chaque Freq
cell_by_cluster <- cell_by_cluster %>%
  mutate(Percentage = (Freq / Total_Freq) * 100)


p <- ggplot(cell_by_cluster, aes(fill = Var1, y = Freq, x = Var2, 
                                 text = paste0("Cluster: ", Var1, 
                                               "<br>Nombre de cellules: ", Freq, 
                                               "<br>Pourcentage: ", round(Percentage, 1), "%")))+ 
  geom_bar(position = "fill", stat = "identity", color = "white") +
  scale_x_discrete(labels = c(SAMPLE_NAME_1, SAMPLE_NAME_2, SAMPLE_NAME_3, SAMPLE_NAME_4, SAMPLE_NAME_5, SAMPLE_NAME_6, SAMPLE_NAME_7, SAMPLE_NAME_8, SAMPLE_NAME_9)) +
  xlab("") +
  ylab("Frequence") +
  labs(fill = "Cluster") +
  theme_classic(base_size = 15) +
  theme(axis.text.x = element_text(color = "black", face = "bold", angle = 45, hjust = 1))

# Save the graph
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_cell_proportion.png")), plot = p , device ="png", width = 12, height = 10, dpi = 300)


# Rendre le graphique interactif avec plotly
p_interactive <- ggplotly(p, tooltip = "text")

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




# NORMALIZE DATA
# --------------

## @knitr normalizeData

sc10x = NormalizeData( object = sc10x,
                       normalization.method = DATA_NORM_METHOD,
                       scale.factor = DATA_NORM_SCALEFACTOR,
                       verbose = .VERBOSE)


gc()

sc10x = ScaleData( object    = sc10x,
                   do.center = DATA_CENTER,
                   do.scale  = DATA_SCALE,
                   vars.to.regress = DATA_VARS_REGRESS,
                   verbose = .VERBOSE)
gc()








