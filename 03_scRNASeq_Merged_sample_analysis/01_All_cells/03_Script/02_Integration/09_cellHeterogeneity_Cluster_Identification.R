

### To save the cell proportions in a CSV file

# Loading the seurat object
PATH_SEURAT_RDS_FILE = file.path( PATH_EXPERIMENT_OUTPUT,
                                  ANALYSIS_STEP_NAME,
                                 "Resolution_0.6",
                                  paste0( outputFilesPrefix, "seuratObject_Harmony.RDS"))
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

## remove line breaks

# Remove accidental line breaks in cluster names
sc10x$cluster_renames <- gsub("\n", " ", sc10x$cluster_renames)
sc10x$cluster_renames <- trimws(sc10x$cluster_renames)  # Remove leading/trailing spaces
sc10x$cluster_renames <- gsub("\\s+", " ", sc10x$cluster_renames)  # Replace multiple spaces with single space

# New metadata column
sc10x$monocyte_status <- ifelse(
  grepl("monocyte", sc10x$cluster_renames, ignore.case = TRUE),
  "monocyte",
  "other"
)

## split by condition (chronic vs sympto) and highligh monocytes

monocyte_umap_split <- DimPlot(sc10x, group.by = "monocyte_status", split.by = "condition", reduction = "umap") +
  ggtitle("Monocytes by condition") +
  scale_color_manual(values = c("monocyte" = "orange", "other" = "purple"))
  
#save plot 
ggsave(plot = monocyte_umap_split, filename = file.path(PATH_ANALYSIS_OUTPUT, "Resolution_0.6", paste0(outputFilesPrefix, "_monocyte_umap_split_condition.png"))
       , width = 30, height = 13, units ="cm")
