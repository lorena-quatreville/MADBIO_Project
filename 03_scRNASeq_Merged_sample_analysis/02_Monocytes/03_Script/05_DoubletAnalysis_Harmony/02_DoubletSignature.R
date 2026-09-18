



## -----------------------------------------
## Load all cluster-set (SciGeneX) objects
## -----------------------------------------

cluster_sets <- lapply(PATH_CLUSTERSET_RDS_FILE, readRDS)

names(cluster_sets) <- basename(unlist(PATH_CLUSTERSET_RDS_FILE))

cat("Loaded cluster-set objects:\n")
print(names(cluster_sets))


## -------------------------------------------------------------
## Define the SciGeneX module IDs per cell type
##    
## -------------------------------------------------------------

celltype_markers <- list(
  B_cell              = "CD79A",      
  T_cell              = "CD3D",
  Classical_monocyte  = "S100A8",
  NonClassical_mono   = "FCGR3A"
)

find_modules_by_marker <- function(cluster_set, marker_genes) {

  modules <- list()

  for (m in seq_along(cluster_set@gene_clusters)) {
    module_genes <- cluster_set@gene_clusters[[m]]

    if (length(intersect(module_genes, marker_genes)) > 0) {
      modules <- c(modules, m)
    }
  }

  unique(unlist(modules))
}

# Extract genes per cell type across all samples
extract_genes_by_celltype <- function(celltype) {

  all_genes <- c()

  for (obj in cluster_sets) {

    module_ids <- find_modules_by_marker(
      cluster_set  = obj,
      marker_genes = celltype_markers[[celltype]]
    )

    if (length(module_ids) == 0) {
      message("⚠ No module found for ", " in cluster-set: ", obj@name)
      next
    }

    # collect genes from modules
    for (m in module_ids) {
      all_genes <- c(all_genes, obj@gene_clusters[[m]])
    }
  }

  unique(all_genes)
}

# Generate gene lists for each cell type
B_cell_genes            <- extract_genes_by_celltype("B_cell")
T_cell_genes            <- extract_genes_by_celltype("T_cell")
ClassicalMono_genes     <- extract_genes_by_celltype("Classical_monocyte")
NonClassicalMono_genes  <- extract_genes_by_celltype("NonClassical_mono")

# Keeping only genes that are present in the SEurat object
available_genes <- rownames(sc10x)

filter_present <- function(x) intersect(x, available_genes)

B_cell_genes            <- filter_present(unique(B_cell_genes))
T_cell_genes            <- filter_present(unique(T_cell_genes))
ClassicalMono_genes     <- filter_present(unique(ClassicalMono_genes))
NonClassicalMono_genes  <- filter_present(unique(NonClassicalMono_genes))

# Removing genes that are in the intersection between lists
gene_lists <- list(
  B_cell = B_cell_genes,
  T_cell = T_cell_genes,
  ClassicalMono = ClassicalMono_genes,
  NonClassicalMono = NonClassicalMono_genes
)

# Get the names of the gene lists
names <- names(gene_lists)

# Loop through all unique pairs of lists (i and j)
for (i in 1:(length(names) - 1)) {
  for (j in (i + 1):length(names)) {

    list_name_i <- names[i]
    list_name_j <- names[j]

    # Find the genes common to both lists
    common_genes <- intersect(gene_lists[[list_name_i]], gene_lists[[list_name_j]])

    # Remove the common genes from list i
    gene_lists[[list_name_i]] <- setdiff(gene_lists[[list_name_i]], common_genes)

    # Remove the common genes from list j
    gene_lists[[list_name_j]] <- setdiff(gene_lists[[list_name_j]], common_genes)
  }
}

# Assign the filtered lists back to the original variable names
B_cell_genes <- gene_lists$B_cell
T_cell_genes <- gene_lists$T_cell
ClassicalMono_genes <- gene_lists$ClassicalMono
NonClassicalMono_genes <- gene_lists$NonClassicalMono

# Compute a score for each set of markers (each cell type)
sc10x <- AddModuleScore(sc10x, list(B_cell_genes),            name = "Bcell_Score")
sc10x <- AddModuleScore(sc10x, list(T_cell_genes),            name = "Tcell_Score")
sc10x <- AddModuleScore(sc10x, list(ClassicalMono_genes),     name = "ClassicalMono_Score")
sc10x <- AddModuleScore(sc10x, list(NonClassicalMono_genes),  name = "NonClassical_Score")

# Visualisation using violin plots

# Extracting the 4 scores: 
score_features <- grep("_Score[0-9]+$", colnames(sc10x@meta.data), value = TRUE)
print(score_features)



vln_plot <- VlnPlot(
  sc10x,
  features = score_features,
  group.by = "multiplet_type_simplified",
  pt.size = 0,          # avoid clutter
  stack = TRUE,         # stacked violins
  flip = TRUE           # horizontal orientation
) +
  facet_wrap(~ feature, scales = "free_y") +   # 🔥 different y-scale per signature
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +  # overlay boxplot
  ggtitle("Immune Marker Expression by Multiplet Type") +
  theme_minimal(base_size = 14) +
  theme(
  legend.position = "none", # remove legend
    panel.background  = element_rect(fill = "white", color = NA),
    plot.background   = element_rect(fill = "white", color = NA),
    axis.text.y       = element_text(size = 12),
    axis.text.x       = element_text(size = 12, angle = 45, hjust = 1),  # rotate x-axis labels
    strip.text        = element_text(face = "bold", size = 14),
    plot.title        = element_text(hjust = 0.5, face = "bold", size = 16),
    panel.grid.minor  = element_blank()
  )


ggsave(
  filename = file.path(PATH_ANALYSIS_OUTPUT, 
                       paste0(outputFilesPrefix, "VlnPlot_ModuleScores_doubletType.png")),
  plot = vln_plot,
  width = 20,
  height = 12,
  bg = "white"
)


## Adhesion and interaction signature
interaction_gene_list <- list("ITGAM", "ITGAL", "ITGB2", "ICAM1", "SELPLG", "VCAM1", "CD14", "FCGR1A", "FCGR3A", "CR1", "CD47", "SIRPA", "MMP9", "ADAM17")


available_genes <- rownames(sc10x@assays$RNA@data)  # Extract gene names from the RNA assay
interaction_gene_list <- interaction_gene_list[interaction_gene_list %in% available_genes]


# Compute module scores
sc10x <- AddModuleScore(sc10x, features = interaction_gene_list, name = "Interaction_Score")

# Plot UMAP FeaturePlot for the scores
umap_interaction_score <- FeaturePlot(sc10x, features = "Interaction_Score1", cols = c("lightgray", "red")) + ggtitle("Interaction Score")

# Display plots
print(umap_interaction_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "umap_Interaction_score.png")), plot = umap_interaction_score, width = 12, height =8, dpi = 300)

# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x, features = "Interaction_Score1", group.by = "sgx_multiplet", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle("Interaction score for multiplets and singlets")
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "vlnBoxPlot_interaction_score_doublet_singlet.png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)


# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x, features = "Interaction_Score1", group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle("Interaction score by condition")
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "vlnBoxPlot_interaction_score_condition.png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)



# Phagocytosis signature (same than the one compare classical monocytes from the different symptomatic and chronic individuals)

phagocytosis_gene_list <- list("IL1B","SRC", "CLEC7A", "SPHK1", "NCF2", "ITGAL", "ABCA1", "ADORA2A", "TNF", "P2RX7", "FCER1G", "PTX3", "RARA", "RAB5A", "PECAM1", "ARL8B", "TYROBP",  "TLR2" )


available_genes <- rownames(sc10x@assays$RNA@data)  # Extract gene names from the RNA assay
phagocytosis_gene_list <- phagocytosis_gene_list[phagocytosis_gene_list %in% available_genes]


# Compute module scores
sc10x <- AddModuleScore(sc10x, features = phagocytosis_gene_list, name = "Phagocytosis_Score")

# Plot UMAP FeaturePlot for the scores
umap_phagocytosis_score <- FeaturePlot(sc10x, features = "Phagocytosis_Score1", cols = c("lightgray", "red")) + ggtitle("Phagocytosis Score")

# Display plots
print(umap_phagocytosis_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "umap_Phagocytosis_score.png")), plot = umap_phagocytosis_score, width = 12, height =8, dpi = 300)

# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x, features = "Phagocytosis_Score1", group.by = "multiplet_number", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle("Phagocytosis score per singlet-doublet-triplet")
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "vlnBoxPlot_phagocytosis_score_multiplet_number.png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)


# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x, features = "Phagocytosis_Score1", group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle("Phagocytosis score by condition")
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "vlnBoxPlot_phagocytosis_score_condition.png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)






## Ag presentation signature (same than the one compare classical monocytes from the different symptomatic and chronic individuals)


antigen_gene_list <- list("FGL2",   "HAVCR2", "ICAM1",  "ITGAL",  "CTSS",   "CTSL",   "CD68",   "FCER1G", "LILRB2", "ARL8B" )


available_genes <- rownames(sc10x@assays$RNA@data)  # Extract gene names from the RNA assay
antigen_gene_list <- antigen_gene_list[antigen_gene_list %in% available_genes]


# Compute module scores
sc10x <- AddModuleScore(sc10x, features = antigen_gene_list, name = "Antigen_Score")

# Plot UMAP FeaturePlot for the scores
umap_antigen_score <- FeaturePlot(sc10x, features = "Antigen_Score1", cols = c("lightgray", "red")) + ggtitle("Antigen Score")

# Display plots
print(umap_antigen_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "umap_Antigen_score.png")), plot = umap_antigen_score, width = 12, height =8, dpi = 300)

# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x, features = "Antigen_Score1", group.by = "multiplet_number", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle("Antigen presentation score per singlet-doublet-triplet")
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "vlnBoxPlot_antigen_score_multiplet_number.png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)


# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x, features = "Antigen_Score1", group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle("Antigen presentation score by condition")
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "vlnBoxPlot_antigen_score_condition.png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)








## Cytokine signature to reflect CD4 T cell-monocyte interaction

cytokine_signature_T_M <- c(
  "IFNG", "TNF", "IL2", "CCL5", "IL21",   # T-cell induced
  "IL1B", "IL6", "CXCL8", "CCL3", "CCL4", "CCL2", # Monocyte cytokines induced by T cells
  "CXCL9", "CXCL10", "CXCL11"            # IFNγ-response
)


cytokine_signature_T_M <- c(
  "CXCL9", "CXCL10", "CXCL11",
  "IL1B", "CCL4", "IL16",
  "CSF2", "CSF3", "CSF1",
  "IFNG", "TNF",
  "CCL2", "CCL27",
  "NGF"
)

## If needed to load the merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x) <- sc10x$cluster_renames

# Subsetting doublets
sc10x_subset <- subset(
  sc10x,
  idents = grep("doublet", levels(Idents(sc10x)), value = TRUE)
)

Idents(sc10x_subset ) <- sc10x_subset$cluster_renames

# Subsetting Classical monocyte
sc10x_subset <- subset(
  sc10x_subset,
  idents = grep("Classical monocyte", levels(Idents(sc10x_subset)), value = TRUE)
)

# Subsetting T cell
sc10x_subset <- subset(
  sc10x_subset,
  idents = grep("T cell", levels(Idents(sc10x_subset)), value = TRUE)
)


cytokine_signature_T_M <- c("IER3", "IL1B", "MARCKS", "NFKBIA", "NLRP3", "C15orf48", "KDM6B", "ICAM1", "CCL3", "CXCL8", "CD83", "ZFP36", "DUSP1", "PPP1R15A", "CCL4", "CCL4L2", "CCL3L1", "MCL1", "FOSB", "TNFAIP2", "RASGEF1B", "NINJ1", "G0S2", "RNF19B", "APOBEC3A", "FFAR2", "FOS", "SERPINB9", "IFNGR2", "SPHK1", "PLEK", "PLAUR", "OLR1", "PTX3", "EGR1", "EGR3", "IER2", "CSRNP1", "PFKFB3", "NFKBIE", "IRAK2", "IL1A", "ZC3H12A", "CD44", "SRC", "ZEB2", "BCL3", "REL", "SMAD3", "KLF6", "CCRL2", "CXCL2", "ATF3", "DMXL2", "P2RX4", "NFKBIZ", "LCP2", "PRDM1", "LAMB3", "GPR35", "TNF", "C3AR1", "PIM3", "ZC3H12C", "ARL5B", "PLK3", "LITAF", "BTG2", "MFSD2A", "SLC35B2", "ZNF267", "HES1")

cytokine_signature_T_M <- c("ICAM1", "CD44", "PLAUR", "CCRL2")

cytokine_signature_T_M <- c("IL1B", "TNF", "CCL3/4", "CXCL8", "NFKBIA/NFKBIZ", "REL", "BCL3")

cytokine_signature_T_M <- c("VASP","ACTB","CFL1","TUBA1B","RHOA","MSN","HCK","RAB10","LCP1","BCL2A1","UBE2J1")
cytokine_signature_T_M <- c("IFI30","CD68","S100A11","FTL","GRINA","EMP3","FBP1","EIF4A1","LGALS3","SASH3","LAIR1")

available_genes <- rownames(sc10x@assays$RNA@data)
cytokine_signature_T_M <- cytokine_signature_T_M[cytokine_signature_T_M %in% available_genes]


# Compute module scores
sc10x <- AddModuleScore(
  sc10x,
  features = list(cytokine_signature_T_M),
  name = "T_M_interaction"
)

# Plot UMAP FeaturePlot for the scores
umap_cytokine_score <- FeaturePlot(sc10x, features = "T_M_interaction1", cols = c("lightgray", "red")) + ggtitle("Cytokine T-monocyte inetraction Score")

# Display plots
print(umap_cytokine_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, paste0( outputFilesPrefix, "umap_cytokine_Tmono_score.png")), plot = umap_cytokine_score, width = 12, height =8, dpi = 300)

# create a new metadata column just for the plotting 
sc10x$interaction_group <- NA

# Singlets
sc10x$interaction_group[sc10x$multiplet_number == "singlet" & sc10x$condition == "Chronic"] <- "Chronic_singlet"
sc10x$interaction_group[sc10x$multiplet_number == "singlet" & sc10x$condition == "Symptomatic"] <- "Symptomatic_singlet"

# Multiplets containing a T cell
sc10x$interaction_group[
  grepl("T cell", sc10x$multiplet_type_simplified, ignore.case = TRUE) &
  sc10x$condition == "Chronic"
] <- "Chronic_Tcell_multiplet"

sc10x$interaction_group[
  grepl("T cell", sc10x$multiplet_type_simplified, ignore.case = TRUE) &
  sc10x$condition == "Symptomatic"
] <- "Symptomatic_Tcell_multiplet"



# violin plot with boxplot:
violin_plot <- VlnPlot(
  sc10x,
  features = "T_M_interaction1",
  group.by = "interaction_group",
  pt.size = 0
) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle("CD4–Monocyte Interaction Signature by Condition and Multiplet Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    
    paste0(outputFilesPrefix, "vln_T_M_interaction_by_group.png")
  ),
  plot = violin_plot,
  width = 12,
  height = 8,
  dpi = 300
)

























