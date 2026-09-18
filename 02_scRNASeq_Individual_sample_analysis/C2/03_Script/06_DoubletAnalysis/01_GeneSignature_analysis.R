


     
########################## To rename clusters on the Seurat object created after the SciGeneX analysis ##########################

cellType = "mono"


# loading the Seurat object with the new clustering based on the SciGeneX analysis
sc10x <- readRDS(file = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, cellType, "_seuratObject_after_SciGeneX.RDS")))


# Retrieving gene names for the OxPhos pathway


# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Suppose your TSV file is 'reactome_oxphos.tsv'
oxphos_path <- read.delim(PATH_REFERENCE_FILE , header = TRUE, sep = "\t")

# 2nd column is UniProt ID
uniprot_ids <- unique(oxphos_path[,2])

# Map UniProt -> Gene Symbol
mapping <- getBM(
  attributes = c("uniprotswissprot", "hgnc_symbol"),
  filters = "uniprotswissprot",
  values = uniprot_ids,
  mart = mart
)

head(mapping)


genes_oxphos <- unique(mapping$hgnc_symbol)
genes_oxphos <- genes_oxphos[genes_oxphos != ""]   # remove empty entries

# This is now a vector you can use in AddModuleScore
head(genes_oxphos)
length(genes_oxphos)

# Ensure genes are present in Seurat object
available_genes <- rownames(sc10x@assays$RNA@data)  # Extract gene names from the RNA assay
genes_oxphos <- genes_oxphos[genes_oxphos %in% available_genes]

# This is now a vector you can use in AddModuleScore
head(genes_oxphos)
length(genes_oxphos)


# Example usage in Seurat
sc10x <- AddModuleScore(sc10x,
                         features = list(genes_oxphos),
                         name = "OxPhos_score")
                         
                         

# Result is stored in seurat@meta.data$OxPhos1


print("OxPhos Genes:")
print(list(genes_oxphos))

# Plot UMAP FeaturePlot for the scores
umap_oxphos_score <- FeaturePlot(sc10x, features = "OxPhos_score1", cols = c("lightgray", "red")) + ggtitle("OxPhos Score")

# Display plots
print(umap_oxphos_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, "umap_OxPhos_score.png")), plot = umap_oxphos_score, width = 12, height =8, dpi = 300)

# Generate Violin + Boxplot
violin_boxplot <- VlnPlot(
  sc10x,
  features = "OxPhos_score1",
  group.by = "cluster_renames"
) +
  ggtitle("OxPhos Score by cluster") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.shape = NA)  # Overlay boxplot

# Print and save the plot
print(violin_boxplot)
ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT, cellType, paste0(outputFilesPrefix, "vlnBoxPlot_OxPhos_score.png")),
       plot = violin_boxplot, width = 12, height = 8, dpi = 300)
















































