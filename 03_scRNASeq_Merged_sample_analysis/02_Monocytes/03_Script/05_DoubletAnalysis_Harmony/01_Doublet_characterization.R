####################  Data loading and Doublet Visualization ###############################

## Loading the sc10x merged object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

## Plotting doublets based on SciGeneX clustering and identification

# Create sgx_multiplet column
sc10x$sgx_multiplet <- ifelse(
  grepl("doublet|triplet", sc10x$cluster_renames, ignore.case = TRUE),
  "multiplet",
  "singlet"
)

# Convert to factor for consistency
sc10x$sgx_multiplet <- factor(sc10x$sgx_multiplet, levels = c("singlet", "multiplet"))

# Check distribution
table(sc10x$sgx_multiplet)

sgx_multiplet_umap <- DimPlot(sc10x, group.by = "sgx_multiplet", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "#BDBDBD", "multiplet" = "orange")) +
  ggtitle("Multiplet prediction based on SciGeneX signatures")
  
# Creat the new direcotory
dir.create(file.path(PATH_ANALYSIS_OUTPUT))

#save plot 
ggsave(plot = sgx_multiplet_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_SGX_multiplet_prediction_umap.png"))
       , width = 23, height = 13, units ="cm")



# Multiplet prediction UMAP split by condition
sgx_multiplet_umap_split <- DimPlot(
  sc10x,
  group.by = "sgx_multiplet",
  split.by = "condition",
  reduction = "umap",
  pt.size = 0.2
) +
  scale_color_manual(
    values = c("singlet" = "#BDBDBD", "multiplet" = "orange")
  ) +
  guides(
  color = guide_legend(
    override.aes = list(size = 5)
  )
  ) +
  theme(
    text = element_text(family = "Arial"),
    
    # Axis labels and numbers
    axis.title = element_text(
      family = "Arial",
      size = 22
    ),
    axis.text = element_text(
      family = "Arial",
      size = 22
    ),
    
    # Remove Chronic / Symptomatic labels above panels
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(
      family = "Arial",
      size = 25
    ),
    legend.key.size = unit(1.2, "cm"),
)
# Save plot
ggsave(
  plot = sgx_multiplet_umap_split,
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "_SGX_multiplet_prediction_umap_split.pdf"
    )
  ),
  width = 30,
  height = 13,
  units = "cm",
  dpi = 600, device = cairo_pdf
)

# Same without the legend:

# Multiplet prediction UMAP split by condition
sgx_multiplet_umap_split <- DimPlot(
  sc10x,
  group.by = "sgx_multiplet",
  split.by = "condition",
  reduction = "umap",
  pt.size = 0.2
) +
  scale_color_manual(
    values = c("singlet" = "#BDBDBD", "multiplet" = "orange")
  ) +
  guides(
  color = guide_legend(
    override.aes = list(size = 5)
  )
  ) +
  theme(
    text = element_text(family = "Arial"),
    
    # Axis labels and numbers
    axis.title = element_text(
      family = "Arial",
      size = 22
    ),
    axis.text = element_text(
      family = "Arial",
      size = 22
    ),
    
    # Remove Chronic / Symptomatic labels above panels
    strip.text = element_blank(),
    plot.title = element_blank(),
 # Remove legend
    legend.position = "none"
)
# Save plot
ggsave(
  plot = sgx_multiplet_umap_split,
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "_SGX_multiplet_prediction_umap_split_no_legend.pdf"
    )
  ),
  width = 30,
  height = 13,
  units = "cm",
  dpi = 600, device = cairo_pdf
)















## same but split by condition (chronic vs sympto)

sgx_multiplet_umap_split <- DimPlot(sc10x, group.by = "sgx_multiplet", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "#BDBDBD", "multiplet" = "orange")) +
  ggtitle("Multiplet prediction based on SciGeneX signatures")
  
#save plot 
ggsave(plot = sgx_multiplet_umap_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_SGX_multiplet_prediction_umap_split.png"))
       , width = 30, height = 13, units ="cm")


# Same without legend and title

sgx_multiplet_umap_split <- DimPlot(sc10x, group.by = "sgx_multiplet", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "#BDBDBD", "multiplet" = "orange")) +
  theme( plot.title = element_blank(), legend.position = "none")
  
#save plot 
ggsave(plot = sgx_multiplet_umap_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_SGX_multiplet_prediction_umap_split_no_legend.png"))
       , width = 30, height = 13, units ="cm")
       
      
       
       

#####################################" To plot singlet/multiplet proportion ####################
# Creating a multiplet summary dataframe
singlet_multi_prop <- sc10x@meta.data %>%
  group_by(orig.ident, sgx_multiplet) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = count / sum(count))

# to change sample order
sample_order <- c(
  "16_scRNASeq_A83_AS_tp2",
  "04_scRNASeq_A81_AS_tp2",
  "13_scRNASeq_A62_AS_tp2",
  "14_scRNASeq_A75_AS_tp2",
  "15_scRNASeq_A76_S_tp2",
  "12_scRNASeq_A86_S_tp2",
    "08_scRNASeq_A82_S_tp2"
)

# to order individuals for the plotting  
  singlet_multi_prop$orig.ident <- factor(
  singlet_multi_prop$orig.ident,
  levels = sample_order
)
# Changing identifyers
id_map <- c(
  "16_scRNASeq_A83_AS_tp2" = "C1",
  "04_scRNASeq_A81_AS_tp2" = "C2",
  "13_scRNASeq_A62_AS_tp2" = "C3",
  "14_scRNASeq_A75_AS_tp2" = "C4",
  "15_scRNASeq_A76_S_tp2" = "S1",
  "12_scRNASeq_A86_S_tp2" = "S2",
    "08_scRNASeq_A82_S_tp2" = "S3") 
    
 singlet_multi_prop$sample_id <- id_map[as.character(singlet_multi_prop$orig.ident)]
 
    
# histogram stacked bar plot
doublet_type_proportion <- ggplot(singlet_multi_prop, aes(x = sample_id, y = prop, fill = sgx_multiplet)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("singlet" = "#BDBDBD", "multiplet" = "orange")) +
  scale_y_continuous(limits = c(0, 1),
  breaks = seq(0, 1, 0.2)) +
  labs(
    title = NULL,
    x = "Sample",
    y = "Proportion of cells",
    fill = "Singlet Multiplet proportion"
  ) +
  theme_bw(base_size = 30) +
  theme(
  text = element_text(family = "Arial"),
  axis.title = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text.x = element_text(angle = 0,
      hjust = 0.5,
      family = "Arial",
      size = 30),
      legend.title = element_blank(),
    legend.text = element_text(
      family = "Arial",
      size = 30
    ),
    legend.key.size = unit(1.2, "cm"),
    panel.grid.minor = element_blank(),
        plot.title = element_blank()
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_singlet_multiplet_proportion.pdf")), plot = doublet_type_proportion , width = 12, height = 10, dpi = 600, device = cairo_pdf)


# Same no legend

# histogram stacked bar plot
doublet_type_proportion <- ggplot(singlet_multi_prop, aes(x = sample_id, y = prop, fill = sgx_multiplet)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width=0.8) +
  scale_fill_manual(values = c("singlet" = "#BDBDBD", "multiplet" = "orange")) +
  scale_y_continuous(limits = c(0, 1),
  breaks = seq(0, 1, 0.2)) +
  labs(
    title = NULL,
    x = "Sample",
    y = "Proportion of cells",
    fill = "Singlet Multiplet proportion"
  ) +
  theme_bw(base_size = 30) +
  theme(
  text = element_text(family = "Arial"),
  axis.title = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text.x = element_text(angle = 0,
      hjust = 0.5,
      family = "Arial",
      size = 30),
    panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        legend.position = "none"
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_singlet_multiplet_proportion_no_legend.pdf")), plot = doublet_type_proportion , width = 12, height = 10, dpi = 600, device = cairo_pdf)

# histogram stacked bar plot (without the legend png)
doublet_type_proportion <- ggplot(singlet_multi_prop, aes(x = orig.ident, y = prop, fill = sgx_multiplet)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("singlet" = "#BDBDBD", "multiplet" = "orange")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = NULL,
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Singlet Multiplet proportion"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position="none"
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_singlet_multiplet_proportion_no_legend.png")), plot = doublet_type_proportion , device ="png", width = 12, height = 10, dpi = 300)
















##################################             ###################################
                                # Data formatting #
##################################              ##################################

# # Expected columns in adata.obs: # - "singlet_type" -> singlet identity (from "cell_type_simplified" in Seurat) # - "sgx_multiplet" -> labels like "singlet" or "multiplet" # - "doublet_type" -> type of doublet, e.g., "CD4T|Monocyte" (NA for singlets)


sc10x$multiplet_type <- sapply(sc10x$cluster_renames, function(x) {
  
  if (grepl("triplet", x, ignore.case = TRUE)) {
    # Take first 3 parts
    parts <- str_split(x, "_")[[1]]
    paste(parts[1:min(3, length(parts))], collapse = "_")
    
  } else if (grepl("doublet", x, ignore.case = TRUE)) {
    # Take first 2 parts
    parts <- str_split(x, "_")[[1]]
    paste(parts[1:min(2, length(parts))], collapse = "_")
    
  } else {
    # Singlet: take first part
    parts <- str_split(x, "_")[[1]]
    parts[1]
  }
})


# Creating a multiplet summary dataframe
multiplet_summary <- sc10x@meta.data %>%
  group_by(orig.ident, multiplet_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = count / sum(count))



#####################################" To plot each multiplet type proportion ####################
# to order individuals for the plotting  
  multiplet_summary$orig.ident <- factor(
  multiplet_summary$orig.ident,
  levels = sample_order
)

# histogram stacked bar plot
doublet_type_proportion <- ggplot(multiplet_summary, aes(x = orig.ident, y = prop, fill = multiplet_type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of Multiplet Types per Sample",
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Multiplet Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_multiplet_type_proportion.png")), plot = doublet_type_proportion , device ="png", width = 12, height = 10, dpi = 300)

# Same without legend
# histogram stacked bar plot
doublet_type_proportion <- ggplot(multiplet_summary, aes(x = orig.ident, y = prop, fill = multiplet_type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = NULL,
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Multiplet Type"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_multiplet_type_proportion_no_legend.png")), plot = doublet_type_proportion , device ="png", width = 12, height = 10, dpi = 300)

############## plotting simplified multiplet type by creating a new metadata column in sc10x

# --- Step 1: Extract simplified multiplet type
sc10x$multiplet_type_simplified <- sapply(sc10x$cluster_renames, function(x) {
  
  if (grepl("triplet", x, ignore.case = TRUE)) {
    # Take first 3 parts
    parts <- str_split(x, "_")[[1]]
    paste(parts[1:min(3, length(parts))], collapse = "_")
    
  } else if (grepl("doublet", x, ignore.case = TRUE)) {
    # Take first 2 parts
    parts <- str_split(x, "_")[[1]]
    paste(parts[1:min(2, length(parts))], collapse = "_")
    
  } else {
    # Singlet: take first part
    parts <- str_split(x, "_")[[1]]
    parts[1]
  }
})

# --- Step 2: Remove unwanted categories
remove_categories <- c(
  "Proliferating cell", 
  "T cell", 
  "Pro-inflammatory cell", 
  "No mod", 
  "Acutely inflammatory cell"
)

sc10x$multiplet_type_simplified[sc10x$multiplet_type_simplified %in% remove_categories] <- "Other"

# --- Step 3: Optional: convert to factor with desired order
# This ensures consistent order in plots
multiplet_levels <- unique(sc10x$multiplet_type_simplified)
sc10x$multiplet_type_simplified <- factor(sc10x$multiplet_type_simplified, levels = multiplet_levels)

# Creating a multipleet summary dataframe
multiplet_summary_simplified <- sc10x@meta.data %>%
  group_by(orig.ident, multiplet_type_simplified) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = count / sum(count))

# to order individuals for the plotting  
  multiplet_summary_simplified$orig.ident <- factor(
  multiplet_summary_simplified$orig.ident,
  levels = sample_order
)

# creating cluster numbers for each multiplet type to enhance visualization
# unique categories
types <- sort(unique(sc10x$multiplet_type_simplified))

# assign numbers 1, 2, 3...
type_ids <- seq_along(types)

# map category → number
type_lookup <- setNames(type_ids, types)

# create metadata column with numeric cluster ID
sc10x$mt_simplified_id <- as.numeric(factor(sc10x$multiplet_type_simplified, levels = names(type_lookup)))

palette_12 <- c(
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
  "#666666"  # dark gray
)

palette_use <- palette_12[seq_len(length(types))]
names(palette_use) <- types

palette_use[c(
  "Classical monocyte",
  "Non-classical monocyte",
  "Intermediate monocyte"
)] <- "#BDBDBD"

multiplet_summary_simplified$sample_id <- id_map[as.character(multiplet_summary_simplified$orig.ident)] 

# histogram stacked bar plot
doublet_type_proportion_simplified <- ggplot(multiplet_summary_simplified, aes(x = sample_id, y = prop, fill = multiplet_type_simplified)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = palette_use) +
  scale_y_continuous(limits = c(0, 1),
  breaks = seq(0, 1, 0.2)) +
  labs(
    title = NULL,
    x = "Sample",
    y = "Proportion of cells",
    fill = "Multiplet Type"
  ) +
  theme_bw(base_size = 30) +
  theme(
  text = element_text(family = "Arial"),
  axis.title = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text.x = element_text(angle = 0,
      hjust = 0.5,
      family = "Arial",
      size = 30),
      legend.title = element_blank(),
    legend.text = element_text(
      family = "Arial",
      size = 30
    ),
    legend.key.size = unit(1.2, "cm"),
    panel.grid.minor = element_blank(),
        plot.title = element_blank()
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_multiplet_type_proportion_simplified.pdf")), plot = doublet_type_proportion_simplified , width = 12, height = 10, dpi = 600, device = cairo_pdf)


# Same without the legend
# histogram stacked bar plot
doublet_type_proportion_simplified <- ggplot(multiplet_summary_simplified, aes(x = sample_id, y = prop, fill = multiplet_type_simplified)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.8) +
  scale_fill_manual(values = palette_use) +
  scale_y_continuous(limits = c(0, 1),
  breaks = seq(0, 1, 0.2)) +
  labs(
    title = NULL,
    x = "Sample",
    y = "Proportion of cells",
    fill = "Multiplet Type"
  ) +
  theme_bw(base_size = 30) +
  theme(
  text = element_text(family = "Arial"),
  axis.title = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text.x = element_text(angle = 0,
      hjust = 0.5,
      family = "Arial",
      size = 30),
    panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        legend.position ="none"
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_multiplet_type_proportion_simplified_no_legend.pdf")), plot = doublet_type_proportion_simplified , width = 12, height = 10, dpi = 600, device = cairo_pdf)


# Multiplet prediction UMAP split by condition
sgx_multiplet_umap_split <- DimPlot(
  sc10x,
  group.by = "multiplet_type_simplified",
  split.by = "condition",
  reduction = "umap",
  cols = palette_use,
  pt.size = 0.2
) +
  guides(
  color = guide_legend(
    override.aes = list(size = 5)
  )
  ) +
  theme(
    text = element_text(family = "Arial"),
    
    # Axis labels and numbers
    axis.title = element_text(
      family = "Arial",
      size = 22
    ),
    axis.text = element_text(
      family = "Arial",
      size = 22
    ),
    
    # Remove Chronic / Symptomatic labels above panels
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(
      family = "Arial",
      size = 25
    ),
    legend.key.size = unit(0.8, "cm"),
)
# Save plot
ggsave(
  plot = sgx_multiplet_umap_split,
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "_multiplet_simplified_umap.pdf"
    )
  ),
  width = 30,
  height = 13,
  units = "cm",
  dpi = 600, device = cairo_pdf
)

# Same without legend
# Multiplet prediction UMAP split by condition
sgx_multiplet_umap_split <- DimPlot(
  sc10x,
  group.by = "multiplet_type_simplified",
  split.by = "condition",
  reduction = "umap",
  cols = palette_use,
  pt.size = 0.2
) +
  guides(
  color = guide_legend(
    override.aes = list(size = 5)
  )
  ) +
  theme(
    text = element_text(family = "Arial"),
    
    # Axis labels and numbers
    axis.title = element_text(
      family = "Arial",
      size = 22
    ),
    axis.text = element_text(
      family = "Arial",
      size = 22
    ),
    
    # Remove Chronic / Symptomatic labels above panels
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
)
# Save plot
ggsave(
  plot = sgx_multiplet_umap_split,
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "_multiplet_simplified_umap_no_legend.pdf"
    )
  ),
  width = 30,
  height = 13,
  units = "cm",
  dpi = 600, device = cairo_pdf
)



# Same png

# Plot the map with clusters with ggplot
umap_multiplet_simplified <- DimPlot(
    sc10x,
    reduction = ifelse(exists("useReduction"), useReduction, "umap"),
    group.by = "multiplet_type_simplified",
    cols = palette_use,
    label = FALSE     # ensures NO labels appear on the UMAP
) + 
  ggtitle("Map of cells with clusters")

#save plot 
ggsave(plot = umap_multiplet_simplified, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_multiplet_simplified_umap.png"))
       , width = 23, height = 13, units ="cm")
       




# Same but split by condition

# Plot the map with clusters with ggplot
umap_multiplet_simplified_split <- DimPlot(
    sc10x,
    reduction = ifelse(exists("useReduction"), useReduction, "umap"),
    group.by = "multiplet_type_simplified",
    split.by = "condition",
    cols = palette_use,
    label = FALSE     # ensures NO labels appear on the UMAP
) + 
  ggtitle("Map of cells with clusters")


#save plot 
ggsave(plot = umap_multiplet_simplified_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_multiplet_simplified_umap_split.png"))
       , width = 30, height = 13, units ="cm")


# Same without the legend
umap_multiplet_simplified_split <- DimPlot(
    sc10x,
    reduction = ifelse(exists("useReduction"), useReduction, "umap"),
    group.by = "multiplet_type_simplified",
    split.by = "condition",
    cols = palette_use,
    label = FALSE     # ensures NO labels appear on the UMAP
) + 
  theme( plot.title = element_blank(), legend.position = "none")


#save plot 
ggsave(plot = umap_multiplet_simplified_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_multiplet_simplified_umap_split_no_legend.png"))
       , width = 30, height = 13, units ="cm")





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

# 2. Violin plot (stacked + flipped) for all marker genes by cluster
vln_plot <- VlnPlot(
  sc10x,
  features = marker_genes,
  group.by = "multiplet_type_simplified",
  pt.size = 0.1,
  stack = TRUE,
  flip = TRUE
) + 
  ggtitle("Immune Marker Expression by Cluster") +
  theme_minimal(base_size = 14) +
  theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "VlnPlot_ImmuneMarkers.png")),
  plot = vln_plot,
  width = 20,
  height = 12,
  bg = "white"
)




##################################             ###################################
                                # Data formatting to create a new metadata column to compare multiplet_number betwen conditions (singlet/doublet/triplet)#
##################################              ##################################

#Creating the new metaata column with multiplet number


sc10x$multiplet_number <- sapply(sc10x$cluster_renames, function(x) {
  
  if (grepl("triplet", x, ignore.case = TRUE)) {
    "triplet"
    
  } else if (grepl("doublet", x, ignore.case = TRUE)) {
    "doublet"
    
  } else {
    "singlet"
  }
})

# Creating a multiplet summary dataframe
multiplet_nbr_summary <- sc10x@meta.data %>%
  group_by(orig.ident, multiplet_number) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = count / sum(count))



##################################### To plot each multiplet_number proportion ####################


# to order individuals for the plotting  
  multiplet_nbr_summary$orig.ident <- factor(
  multiplet_nbr_summary$orig.ident,
  levels = sample_order
)

# To orders multiplet numbers for the plotting
multiplet_nbr_summary$multiplet_number <- factor(
  multiplet_nbr_summary$multiplet_number,
  levels = c("singlet", "doublet", "triplet")
)

multiplet_nbr_summary$sample_id <- id_map[as.character(multiplet_nbr_summary$orig.ident)] 

# histogram stacked bar plot
multiplet_nbr_proportion <- ggplot(multiplet_nbr_summary, aes(x = sample_id, y = prop, fill = multiplet_number)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("singlet" = "#BDBDBD", "doublet" = "green", "triplet" = "purple"))  +
  scale_y_continuous(limits = c(0, 1),
  breaks = seq(0, 1, 0.2)) +
  labs(
    title = NULL,
    x = "Sample",
    y = "Proportion of cells",
    fill = "Singlet Doublet Triplet proportion"
  ) +
  theme_bw(base_size = 30) +
  theme(
  text = element_text(family = "Arial"),
  axis.title = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text.x = element_text(angle = 0,
      hjust = 0.5,
      family = "Arial",
      size = 30),
      legend.title = element_blank(),
    legend.text = element_text(
      family = "Arial",
      size = 30
    ),
    legend.key.size = unit(1.2, "cm"),
    panel.grid.minor = element_blank(),
        plot.title = element_blank()
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_multiplet_number.pdf")), plot = multiplet_nbr_proportion , width = 12, height = 10, dpi = 600, device = cairo_pdf)


# histogram stacked bar plot without legend
multiplet_nbr_proportion <- ggplot(multiplet_nbr_summary, aes(x = sample_id, y = prop, fill = multiplet_number)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2, width = 0.8) +
  scale_fill_manual(values = c("singlet" = "#BDBDBD", "doublet" = "green", "triplet" = "purple"))  +
  scale_y_continuous(limits = c(0, 1),
  breaks = seq(0, 1, 0.2)) +
  labs(
    title = NULL,
    x = "Sample",
    y = "Proportion of cells",
    fill = "Singlet Doublet Triplet proportion"
  ) +
  theme_bw(base_size = 30) +
  theme(
  text = element_text(family = "Arial"),
  axis.title = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text = element_text(
      family = "Arial",
      size = 30
    ),
    axis.text.x = element_text(angle = 0,
      hjust = 0.5,
      family = "Arial",
      size = 30),
    panel.grid.minor = element_blank(),
        plot.title = element_blank(),
        legend.position = "none"
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_multiplet_number_no_legend.pdf")), plot = multiplet_nbr_proportion , width = 12, height = 10, dpi = 600, device = cairo_pdf)

# histogram stacked bar plot PNG
multiplet_nbr_proportion <- ggplot(multiplet_nbr_summary, aes(x = orig.ident, y = prop, fill = multiplet_number)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("singlet" = "#BDBDBD", "doublet" = "green", "triplet" = "purple")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of singlet/doublet/triplet per Sample",
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Singlet Doublet Triplet proportion"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_multiplet_number.png")), plot = multiplet_nbr_proportion , device ="png", width = 12, height = 10, dpi = 300)

# Same without the legend PNG
# histogram stacked bar plot
multiplet_nbr_proportion <- ggplot(multiplet_nbr_summary, aes(x = orig.ident, y = prop, fill = multiplet_number)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("singlet" = "#BDBDBD", "doublet" = "green", "triplet" = "purple")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = NULL,
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Singlet Doublet Triplet proportion"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    legend.position="none"
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_multiplet_number_no_legend.png")), plot = multiplet_nbr_proportion , device ="png", width = 12, height = 10, dpi = 300)





# Plot the corresponding UMAP
umap_multiplet_nbr <- DimPlot(sc10x, group.by = "multiplet_number", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "#BDBDBD", "doublet" = "green", "triplet" = "purple")) +
  ggtitle("Proportion of singlet/doublet/triplet per Sample")
  
#save plot 
ggsave(plot = umap_multiplet_nbr, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_multiplet_number_umap.png"))
       , width = 23, height = 13, units ="cm")



# UMAP PDF

# Multiplet prediction UMAP split by condition
umap_multiplet_nbr_split <- DimPlot(
  sc10x,
  group.by = "multiplet_number",
  split.by = "condition",
  reduction = "umap",
  pt.size = 0.2
) +
    scale_color_manual(values = c("singlet" = "#BDBDBD", "doublet" = "green", "triplet" = "purple"))+
  guides(
  color = guide_legend(
    override.aes = list(size = 5)
  )
  ) +
  theme(
    text = element_text(family = "Arial"),
    
    # Axis labels and numbers
    axis.title = element_text(
      family = "Arial",
      size = 22
    ),
    axis.text = element_text(
      family = "Arial",
      size = 22
    ),
    
    # Remove Chronic / Symptomatic labels above panels
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(
      family = "Arial",
      size = 25
    ),
    legend.key.size = unit(1.2, "cm"),
)
# Save plot
ggsave(
  plot = umap_multiplet_nbr_split,
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "_multiplet_number_umap_split.pdf"
    )
  ),
  width = 30,
  height = 13,
  units = "cm",
  dpi = 600, device = cairo_pdf
)



# Same without the legend

# Multiplet prediction UMAP split by condition
umap_multiplet_nbr_split <- DimPlot(
  sc10x,
  group.by = "multiplet_number",
  split.by = "condition",
  reduction = "umap",
  pt.size = 0.2
) +
    scale_color_manual(values = c("singlet" = "#BDBDBD", "doublet" = "green", "triplet" = "purple"))+
  guides(
  color = guide_legend(
    override.aes = list(size = 5)
  )
  ) +
  theme(
    text = element_text(family = "Arial"),
    
    # Axis labels and numbers
    axis.title = element_text(
      family = "Arial",
      size = 22
    ),
    axis.text = element_text(
      family = "Arial",
      size = 22
    ),
    
    # Remove Chronic / Symptomatic labels above panels
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
)
# Save plot
ggsave(
  plot = umap_multiplet_nbr_split,
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "_multiplet_number_umap_split_no_legend.pdf"
    )
  ),
  width = 30,
  height = 13,
  units = "cm",
  dpi = 600, device = cairo_pdf
)



## same but split by condition (chronic vs sympto) PNG

umap_multiplet_nbr_split <- DimPlot(sc10x, group.by = "multiplet_number", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "#BDBDBD", "doublet" = "green", "triplet" = "purple")) +
  ggtitle("Proportion of singlet/doublet/triplet per Sample between conditions")
  
#save plot 
ggsave(plot = umap_multiplet_nbr_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_multiplet_number_umap_split.png"))
       , width = 30, height = 13, units ="cm")

# Same without the legend PNG

umap_multiplet_nbr_split <- DimPlot(sc10x, group.by = "multiplet_number", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "#BDBDBD", "doublet" = "green", "triplet" = "purple")) +
  theme( plot.title = element_blank(), legend.position = "none")
  
#save plot 
ggsave(plot = umap_multiplet_nbr_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_multiplet_number_umap_split_no_legend.png"))
       , width = 30, height = 13, units ="cm")






##################################             ###################################
                                # Data formatting to create a new metadata column to compare interaction_capacity of classical monocytes and non-classical monocytes betwen conditions (CM/NCM/CM_multiplet/NCM_multiplet)#
##################################              ##################################


# Create mono_interaction_capacity column
sc10x$mono_interaction_capacity <- case_when(

  # Classical monocytes interacting (doublet/triplet or with platelet)
  grepl("Classical monocyte", sc10x$cluster_renames) &
    (grepl("doublet|triplet", sc10x$cluster_renames) | grepl("platelet", sc10x$cluster_renames)) ~
      "Classical monocyte_multiplet",

  # Non-classical monocytes interacting (doublet/triplet)
  grepl("Non-classical monocyte", sc10x$cluster_renames) &
    grepl("doublet|triplet", sc10x$cluster_renames) ~
      "Non-classical monocyte_multiplet",
      
       # Intermediate monocytes interacting (doublet/triplet)
  grepl("Intermediate monocyte", sc10x$cluster_renames) &
    grepl("doublet|triplet", sc10x$cluster_renames) ~
      "Intermediate monocyte_multiplet",

  # Classical singlets
  grepl("Classical monocyte", sc10x$cluster_renames) &
    !grepl("doublet|triplet|platelet", sc10x$cluster_renames) ~
      "Classical monocyte_singlet",

  # Non-classical singlets
  grepl("Non-classical monocyte", sc10x$cluster_renames) &
    !grepl("doublet|triplet", sc10x$cluster_renames) ~
      "Non-classical monocyte_singlet",

    # Intermediate singlets
    grepl("Intermediate monocyte", sc10x$cluster_renames) &
    !grepl("doublet|triplet", sc10x$cluster_renames) ~
      "Intermediate monocyte_singlet"
)


sc10x$mono_interaction_capacity <- factor(
  sc10x$mono_interaction_capacity,
  levels = c(
    "Classical monocyte_singlet",
    "Classical monocyte_multiplet",
    "Non-classical monocyte_singlet",
    "Non-classical monocyte_multiplet",
    "Intermediate monocyte_singlet",
    "Intermediate monocyte_multiplet"
  )
)

# keep only monocytes (classical / non-classical)
mono_only <- sc10x@meta.data %>%
  filter(mono_interaction_capacity %in% c(
    "Classical monocyte_singlet",
    "Classical monocyte_multiplet",
    "Non-classical monocyte_singlet",
    "Non-classical monocyte_multiplet",
    "Intermediate monocyte_singlet",
    "Intermediate monocyte_multiplet"
  ))

# summary with correct denominator
mono_inter_capacity_summary <- mono_only %>%
  group_by(orig.ident, mono_interaction_capacity) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%

  mutate(
    total_classical = sum(count[grepl("Classical", mono_interaction_capacity)]),
    total_nonclassical = sum(count[grepl("Non-classical", mono_interaction_capacity)]),
    total_intermediate = sum(count[grepl("Intermediate", mono_interaction_capacity)]),
    prop = case_when(
      grepl("Classical", mono_interaction_capacity) ~ count / total_classical,
      grepl("Non-classical", mono_interaction_capacity) ~ count / total_nonclassical,
      grepl("Intermediate", mono_interaction_capacity) ~ count / total_intermediate
    )
  )





##################################### To plot each inetracting monocyte proportion ####################
# to change sample order
sample_order <- c(
  "08_scRNASeq_A82_S_tp2",
  "12_scRNASeq_A86_S_tp2",
  "15_scRNASeq_A76_S_tp2",
  "04_scRNASeq_A81_AS_tp2",
  "16_scRNASeq_A83_AS_tp2",
  "13_scRNASeq_A62_AS_tp2",
  "14_scRNASeq_A75_AS_tp2"
)

# to order individuals for the plotting  
  mono_inter_capacity_summary$orig.ident <- factor(
  mono_inter_capacity_summary$orig.ident,
  levels = sample_order
)

# To orders mono_inter_capacity for the plotting
mono_inter_capacity_summary$mono_interaction_capacity <- factor(
  mono_inter_capacity_summary$mono_interaction_capacity,
  levels = c("Classical monocyte_singlet", "Classical monocyte_multiplet", "Non-classical monocyte_multiplet", "Non-classical monocyte_singlet")
)

# histogram stacked bar plot
mono_inter_proportion <- ggplot(mono_inter_capacity_summary, aes(x = orig.ident, y = prop, fill = mono_interaction_capacity)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("Classical monocyte_singlet" = "#88CCEE" , "Non-classical monocyte_singlet" = "#E69F00", "Classical monocyte_multiplet" = "#332288", "Non-classical monocyte_multiplet"="red")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of classical/non-classical monocyte interactions per Sample",
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Interacting monocyte proportion between conditions"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_mono_inter_capacity.png")), plot = mono_inter_proportion , device ="png", width = 12, height = 10, dpi = 300)


# Plot the corresponding UMAP
umap_inter_mono_prop <- DimPlot(sc10x, group.by = "mono_interaction_capacity", reduction = "umap") +
  scale_color_manual(values = c("Classical monocyte_singlet" = "#88CCEE" , "Non-classical monocyte_singlet" = "#E69F00", "Classical monocyte_multiplet" = "#332288", "Non-classical monocyte_multiplet"="red")) +
  ggtitle("Interacting monocyte proportion between conditions")
  
#save plot 
ggsave(plot = umap_inter_mono_prop, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_interacting_mono_umap.png"))
       , width = 23, height = 13, units ="cm")



# UMAP split by condition PDF
umap_inter_mono_prop_split <- DimPlot(
  sc10x,
  group.by = "mono_interaction_capacity",
  split.by = "condition",
  reduction = "umap",
  pt.size = 0.2
) +
    scale_color_manual(values = c("Classical monocyte_singlet" = "#88CCEE" , "Non-classical monocyte_singlet" = "#E69F00", "Classical monocyte_multiplet" = "#332288", "Non-classical monocyte_multiplet"="red", "Intermediate monocyte_singlet" = "#8DB600", "Intermediate monocyte_multiplet" = "#117733")
    )+
  guides(
  color = guide_legend(
    override.aes = list(size = 5)
  )
  ) +
  theme(
    text = element_text(family = "Arial"),
    
    # Axis labels and numbers
    axis.title = element_text(
      family = "Arial",
      size = 22
    ),
    axis.text = element_text(
      family = "Arial",
      size = 22
    ),
    
    # Remove Chronic / Symptomatic labels above panels
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.title = element_blank(),
    legend.text = element_text(
      family = "Arial",
      size = 25
    ),
    legend.key.size = unit(1.2, "cm"),
)
# Save plot
ggsave(
  plot = umap_inter_mono_prop_split,
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "_interacting_mono_umap_split.pdf"
    )
  ),
  width = 30,
  height = 13,
  units = "cm",
  dpi = 600, device = cairo_pdf
)



# Plot the corresponding UMAP
umap_inter_mono_prop <- DimPlot(sc10x, group.by = "mono_interaction_capacity", reduction = "umap") +
  scale_color_manual(values = c("Classical monocyte_singlet" = "#88CCEE" , "Non-classical monocyte_singlet" = "#E69F00", "Classical monocyte_multiplet" = "#332288", "Non-classical monocyte_multiplet"="red", "Intermediate monocyte_singlet" = "#8DB600", "Intermediate monocyte_multiplet" = "#117733")) +
  ggtitle("Interacting monocyte proportion between conditions")
  
#save plot 
ggsave(plot = umap_inter_mono_prop, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_interacting_mono_umap.png"))
       , width = 23, height = 13, units ="cm")



# UMAP split by condition PDF without the legend
umap_inter_mono_prop_split <- DimPlot(
  sc10x,
  group.by = "mono_interaction_capacity",
  split.by = "condition",
  reduction = "umap",
  pt.size = 0.2
) +
    scale_color_manual(values = c("Classical monocyte_singlet" = "#88CCEE" , "Non-classical monocyte_singlet" = "#E69F00", "Classical monocyte_multiplet" = "#332288", "Non-classical monocyte_multiplet"="red", "Intermediate monocyte_singlet" = "#8DB600", "Intermediate monocyte_multiplet" = "#117733")
    )+
  guides(
  color = guide_legend(
    override.aes = list(size = 5)
  )
  ) +
  theme(
    text = element_text(family = "Arial"),
    
    # Axis labels and numbers
    axis.title = element_text(
      family = "Arial",
      size = 22
    ),
    axis.text = element_text(
      family = "Arial",
      size = 22
    ),
    
    # Remove Chronic / Symptomatic labels above panels
    strip.text = element_blank(),
    plot.title = element_blank(),
    legend.position = "none"
)
# Save plot
ggsave(
  plot = umap_inter_mono_prop_split,
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(
      outputFilesPrefix,
      "_interacting_mono_umap_split_no_legend.pdf"
    )
  ),
  width = 30,
  height = 13,
  units = "cm",
  dpi = 600, device = cairo_pdf
)



## same but split by condition (chronic vs sympto) PNG

umap_inter_mono_prop_split <- DimPlot(sc10x, group.by = "mono_interaction_capacity", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("Classical monocyte_singlet" = "#88CCEE" , "Non-classical monocyte_singlet" = "#E69F00", "Classical monocyte_multiplet" = "#332288", "Non-classical monocyte_multiplet"="red", "Intermediate monocyte_singlet" = "#8DB600", "Intermediate monocyte_multiplet" = "#117733")) +
  ggtitle("Interacting monocyte proportion between conditions")
  
#save plot 
ggsave(plot = umap_inter_mono_prop_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_interacting_mono_umap_split.png"))
       , width = 30, height = 13, units ="cm")
       
         
# Same without legend PNG
umap_inter_mono_prop_split <- DimPlot(sc10x, group.by = "mono_interaction_capacity", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("Classical monocyte_singlet" = "#88CCEE" , "Non-classical monocyte_singlet" = "#E69F00", "Classical monocyte_multiplet" = "#332288", "Non-classical monocyte_multiplet"="red")) +
   theme(legend.position = "none", plot.title = element_blank())
  
#save plot 
ggsave(plot = umap_inter_mono_prop_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_interacting_mono_umap_split_no_legend.png"))
       , width = 30, height = 13, units ="cm")


# Save the seurat object
saveRDS( object = sc10x, file = PATH_SEURAT_RDS_FILE)


##############  Violin plot representing the % of mitochondrial genes in multiplet vs singlet

# Plot the cells percentage of mitochondrial genes in UMAP  
mitok_umap <- FeaturePlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
             features = "percent.mito") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) +
  ggtitle( "Map of cells with level of percentage of\n mitochondrial genes")
ggsave(
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    
    paste0(outputFilesPrefix, "umap_singlet_multi_mitok.png")
  ),
  plot = mitok_umap,
  width = 12,
  height = 8,
  dpi = 300
)

# if QC EXPLORATION MODE is active, plot if the cells are excluded by %of mitochondrial genes or not
umap_QC_mitok <- DimPlot( sc10x, reduction = ifelse(exists("useReduction"), useReduction, "umap"),
                  group.by = "outlier.percent.mito") + 
           ggtitle( "Map of cells by filter status\n on %mitochondrial gene")
ggsave(
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    
    paste0(outputFilesPrefix, "umap_QC_singlet_multi_mitok.png")
  ),
  plot = umap_QC_mitok,
  width = 12,
  height = 8,
  dpi = 300
)


# violin plot with boxplot:
mitok_violin_plot <- VlnPlot(
  sc10x,
  features = "percent.mito",
  group.by = "sgx_multiplet",
  pt.size = 0
) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle("% of mitochondrial genes in singlet vs multiplet") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    
    paste0(outputFilesPrefix, "vln_singlet_multi_mitok.png")
  ),
  plot = mitok_violin_plot,
  width = 12,
  height = 8,
  dpi = 300
)



# Table nb multiplet mono/singlet

table_singlet_multi = table(sc10x$orig.ident, sc10x$sgx_multiplet)


# Save to CSV file
write.csv(as.data.frame.matrix(table_singlet_multi), file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "nb_singlet_multiplet_per_indiv.csv")), row.names = TRUE)


