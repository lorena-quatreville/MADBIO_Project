####################  Data loading and Doublet Visualization ###############################

## Loading the sc10x merged object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

############ All multiplets taken together


# Create sgx_multiplet column
sc10x$sgx_multiplet <- ifelse(
  grepl("doublet|triplet", sc10x$sgx_clusters, ignore.case = TRUE),
  "multiplet",
  "singlet"
)

# Convert to factor for consistency
sc10x$sgx_multiplet <- factor(sc10x$sgx_multiplet, levels = c("singlet", "multiplet"))

# Check distribution
table(sc10x$sgx_multiplet)

sgx_multiplet_umap <- DimPlot(sc10x, group.by = "sgx_multiplet", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "orange", "multiplet" = "purple")) +
  ggtitle("Multiplet prediction based on SciGeneX signatures")
  
#save plot 
ggsave(plot = sgx_multiplet_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_SGX_multiplet_prediction_umap_all_cells.png"))
       , width = 23, height = 13, units ="cm")




## same but split by condition (chronic vs sympto)

sgx_multiplet_umap_split <- DimPlot(sc10x, group.by = "sgx_multiplet", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "orange", "multiplet" = "purple")) +
  ggtitle("Multiplet prediction based on SciGeneX signatures")
  
#save plot 
ggsave(plot = sgx_multiplet_umap_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_SGX_multiplet_prediction_umap_split_all_cells.png"))
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
  "08_scRNASeq_A82_S_tp2",
  "12_scRNASeq_A86_S_tp2",
  "15_scRNASeq_A76_S_tp2",
  "04_scRNASeq_A81_AS_tp2",
  "16_scRNASeq_A83_AS_tp2",
  "13_scRNASeq_A62_AS_tp2",
  "14_scRNASeq_A75_AS_tp2"
)

# to order individuals for the plotting  
  singlet_multi_prop$orig.ident <- factor(
  singlet_multi_prop$orig.ident,
  levels = sample_order
)

# histogram stacked bar plot
doublet_type_proportion <- ggplot(singlet_multi_prop, aes(x = orig.ident, y = prop, fill = sgx_multiplet)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("singlet" = "orange", "multiplet" = "purple")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of Multiplet Types per Sample",
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Singlet Multiplet proportion"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "stacked_barplot_singlet_multiplet_proportion_all_cells.png")), plot = doublet_type_proportion , device ="png", width = 12, height = 10, dpi = 300)





                                            ### B cells



# Subseting seurat object to get only B cells
B_cells <- rownames(sc10x@meta.data)[
  grepl("B cell", sc10x@meta.data$cluster_renames, ignore.case = TRUE)
]

sc10x_Bcell <- subset(sc10x, cells = B_cells)


## Plotting doublets based on SciGeneX clustering and identification

# Create sgx_multiplet column
sc10x_Bcell$sgx_multiplet <- ifelse(
  grepl("doublet", sc10x_Bcell$sgx_clusters, ignore.case = TRUE),
  "doublet",
  "singlet"
)

# Convert to factor for consistency
sc10x_Bcell$sgx_multiplet <- factor(sc10x_Bcell$sgx_multiplet, levels = c("singlet", "doublet"))

# Check distribution
table(sc10x_Bcell$sgx_multiplet)

sgx_multiplet_umap <- DimPlot(sc10x_Bcell, group.by = "sgx_multiplet", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "orange", "doublet" = "purple")) +
  ggtitle("Multiplet prediction based on SciGeneX signatures")
  
#save plot 
ggsave(plot = sgx_multiplet_umap, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_SGX_multiplet_prediction_B_cell_umap.png"))
       , width = 23, height = 13, units ="cm")


## same but split by condition (chronic vs sympto)

sgx_multiplet_umap_split <- DimPlot(sc10x_Bcell, group.by = "sgx_multiplet", split.by = "condition", reduction = "umap") +
  scale_color_manual(values = c("singlet" = "orange", "doublet" = "purple")) +
  ggtitle("Multiplet prediction based on SciGeneX signatures")
  
#save plot 
ggsave(plot = sgx_multiplet_umap_split, filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_SGX_multiplet_prediction_B_cell_umap_split.png"))
       , width = 30, height = 13, units ="cm")


#####################################" To plot singlet/multiplet proportion ####################
# Creating a multiplet summary dataframe
singlet_multi_prop <- sc10x_Bcell@meta.data %>%
  group_by(orig.ident, sgx_multiplet) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(prop = count / sum(count))

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
  singlet_multi_prop$orig.ident <- factor(
  singlet_multi_prop$orig.ident,
  levels = sample_order
)

# histogram stacked bar plot
doublet_type_proportion <- ggplot(singlet_multi_prop, aes(x = orig.ident, y = prop, fill = sgx_multiplet)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.2) +
  scale_fill_manual(values = c("singlet" = "orange", "doublet" = "purple")) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    title = "Proportion of Multiplet Types per Sample",
    x = "Sample",
    y = "Proportion of Cells",
    fill = "Singlet Multiplet proportion"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "stacked_barplot_singlet_multiplet_proportion_B_cell.png")), plot = doublet_type_proportion , device ="png", width = 12, height = 10, dpi = 300)



write.csv(singlet_multi_prop,
          file = file.path( PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "singlet_multi_prop.csv")),
          row.names = FALSE)
























