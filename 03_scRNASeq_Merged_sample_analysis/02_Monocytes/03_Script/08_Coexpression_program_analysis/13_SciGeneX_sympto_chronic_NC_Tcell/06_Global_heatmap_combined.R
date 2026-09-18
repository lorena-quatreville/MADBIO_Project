

# Building a global heatmap of all modules accross all cells


###############################################################################
# Extract NC-Tcell modules
###############################################################################
scigx_res_NC_Tcell <- readRDS(file.path("/home/quatreville/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2/05_Output/08_Coexpression_program_analysis/13_SciGeneX_sympto_chronic_NC_Tcell/01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_13_SciGeneX_sympto_chronic_NC_Tcell__clusterset_object_module_activity.RDS"))

modules_NC_Tcell <- scigx_res_NC_Tcell@top_genes
names(modules_NC_Tcell) <- paste0("NC_Tcell_M", seq_along(modules_NC_Tcell))

###############################################################################
# Extract NC modules
###############################################################################
scigx_res_NCmono <- readRDS(file.path("/home/quatreville/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2/05_Output/08_Coexpression_program_analysis/11_SciGeneX_sympto_chronic_NCM/01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_11_SciGeneX_sympto_chronic_NCM__clusterset_object_module_activity.RDS")) 

modules_NCmono <- scigx_res_NCmono@top_genes
names(modules_NCmono) <- paste0("NCmono_M", seq_along(modules_NCmono))

###############################################################################
# Merge all modules
###############################################################################

all_modules <- c(
    modules_NCmono,
    modules_NC_Tcell
)
# Reloading the seurat objects for each group
sc10x_NC_Tcell = readRDS(file.path("/home/quatreville/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2/05_Output/08_Coexpression_program_analysis/13_SciGeneX_sympto_chronic_NC_Tcell/01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_13_SciGeneX_sympto_chronic_NC_Tcell__seuratObject_module_activity.RDS"))
sc10x_NCmono = readRDS(file.path("/home/quatreville/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2/05_Output/08_Coexpression_program_analysis/11_SciGeneX_sympto_chronic_NCM/01_MADBIO_Cameroon_01_scRNASeq_42_scRNASeq_A76_A82_A86_A83_A81_A62_A75_mono_tp2_11_SciGeneX_sympto_chronic_NCM__seuratObject_module_activity.RDS"))

# Merging of the two seurat objects
sc10x = merge(sc10x_NC_Tcell, sc10x_NCmono)

# Keeping the information of which cell type each cluster comes from
sc10x$cell_origin <- ifelse(
    sc10x$multiplet_type == "Non-classical monocyte",
    "NCmono",
    "NCmono_Tcell"
)

sc10x$cluster_origin <- paste0(
    sc10x$cell_origin,
    "_",
    sc10x$heatmap_clusters
)

###############################################################################
# Computing the module score on all cells
###############################################################################

# Select the top genes per module
scigx_res <- top_genes(scigx)


# Compute the UCell score for each cell and each module
  sc10x <- UCell::AddModuleScore_UCell(
    sc10x,
    features = all_modules,
    name = ""
  )


# ###########################################
# Compute activity of gene modules per cell
# ###########################################

# Defining expression matrix from sc10x object and not SciGeneX
expr <- sc10x@assays$RNA@data

# Mean computing
gene_means <- Matrix::rowMeans(expr)

# Binary expression assessment
activity_bin <- expr > gene_means

# Module activity scorign
module_activity_score_names <- c()
module_activity_names <- c()

for (module_name in names(all_modules)) {

  genes <- all_modules[[module_name]]

  genes <- intersect(genes, rownames(activity_bin))

  score <- Matrix::colMeans(activity_bin[genes, , drop = FALSE])

  score_bin <- as.numeric(score >= 0.5)

  sc10x <- AddMetaData(
    sc10x,
    metadata = score,
    col.name = paste0("Activity.Score.", module_name)
  )

  sc10x <- AddMetaData(
    sc10x,
    metadata = score_bin,
    col.name = paste0("Activity.Binary.", module_name)
  )
}

# ###########################################
# Compute activity of gene modules per cluster as the percentage of cells activating this module
# ###########################################


# Loop over the cluster to study them one by one
Idents( sc10x) = "cluster_origin"

# Creating module list
module_names <- names(all_modules)

# Computing cluster × module activity (% cells active)
clusters_modules_activity <- data.frame()

for (cluster in levels(Idents(sc10x))) {

  cluster_cells <- WhichCells(sc10x, idents = cluster)

  cat("\nCluster:", cluster, "n =", length(cluster_cells))

  for (module_name in module_names) {

    bin_col <- paste0("Activity.Binary.", module_name)

    values <- sc10x@meta.data[cluster_cells, bin_col]
    # mean(binary) = fraction of active cells
    percentage_activity <- mean(values, na.rm = TRUE)

    clusters_modules_activity <- rbind(
      clusters_modules_activity,
      data.frame(
        cluster = cluster,
        module = module_name,
        activity.pct = percentage_activity
      )
    )
  }
}

# Factor ordering
clusters_modules_activity$module <- factor(
  clusters_modules_activity$module,
  levels = module_names
)
clusters_modules_activity$cluster <- factor(
  clusters_modules_activity$cluster,
  levels = levels(Idents(sc10x))
)

# Formatting for heatmap plotting
library(tidyr)

clusters_modules_activity_wide <- pivot_wider(
  clusters_modules_activity,
  names_from = cluster,
  values_from = activity.pct
)

clusters_modules_activity_wide_df <- as.data.frame(clusters_modules_activity_wide)
# so that module column become row names and so that all column are numeric and pheatmap doesn't crash
rownames(clusters_modules_activity_wide_df) <- clusters_modules_activity_wide_df$module
clusters_modules_activity_wide_df$module <- NULL





###############################################################################
# Plotting the gloabl heatmap
###############################################################################
# Heatmap expects a matrix with a module column, a cluster_rigin column and an activity score
library(tidyr)

# Heatmap
png(
    file.path(
        PATH_ANALYSIS_OUTPUT,
        paste0(outputFilesPrefix, "_NCmono_vs_NCmonoTcell_module_activity_heatmap.png")
    ),
    width = 3000,
    height = 2500,
    res = 300
)

pheatmap(
    as.matrix(clusters_modules_activity_wide_df),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(
        brewer.pal(7, "YlOrRd")
    )(100)
)


dev.off()



# Ordering clusters and modules
# Extracting cluster names
cluster_levels <- colnames(clusters_modules_activity_wide_df)
#Ordering
cluster_levels <- c(
  grep("^NCmono_", cluster_levels, value = TRUE),
  grep("^NCmono_Tcell_", cluster_levels, value = TRUE)
)

# Ordering module names
module_levels <- names(all_modules)

module_levels <- c(
  grep("^NCmono_M", module_levels, value = TRUE),
  grep("^NCmono_Tcell_M", module_levels, value = TRUE)
)

# numeric sort within each group
module_levels <- c(
  module_levels[order(as.numeric(sub(".*_M", "", module_levels[grep("^NCmono_M", module_levels)])))],
  module_levels[order(as.numeric(sub(".*_M", "", module_levels[grep("^NCmono_Tcell_M", module_levels)])))]
)

clusters_modules_activity_wide_df <- as.data.frame(clusters_modules_activity_wide)
# so that module column become row names and so that all column are numeric and pheatmap doesn't crash
rownames(clusters_modules_activity_wide_df) <- clusters_modules_activity_wide_df$module
clusters_modules_activity_wide_df$module <- NULL


png(
  file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "_NCmono_vs_NCmonoTcell_module_activity_heatmap_ordered.png")
  ),
  width = 3000,
  height = 2500,
  res = 300
)

pheatmap(
  as.matrix(clusters_modules_activity_wide_df),
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(
    RColorBrewer::brewer.pal(7, "YlOrRd")
  )(100)
)

dev.off()




# Module correlation matrix (as pseudo_balk
# 1. define module columns
module_scores <- names(all_modules)
modules_to_remove <- c("NCmono_M11", "NC_Tcell_M3", "NC_Tcell_M10", "NC_Tcell_M12")
module_scores <- setdiff(names(all_modules), modules_to_remove)



# 2. aggregate (pseudobulk per sample/condition)
patient_modules_activity_scoring <- sc10x@meta.data %>%
  group_by(orig.ident, condition) %>%
  summarise(across(all_of(module_scores), median), .groups = "drop")

# 3. matrix
mat <- as.matrix(patient_modules_activity_scoring[, module_scores])

# 4. remove bad columns
mat <- mat[, colSums(is.na(mat)) == 0, drop = FALSE]

# 5. correlation
cor_mat <- cor(mat, use = "pairwise.complete.obs", method = "spearman")

cor_mat[is.na(cor_mat)] <- 0

 # Save correlation matrix as CSV
write.csv(
  cor_mat,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "correlation_matrix_gloabl_module_activity_scoring.csv")
  ),
  row.names = TRUE
)


# 6. heatmap
global_corr_plot = Heatmap(
    cor_mat,
    name = "Spearman",
    col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
)

png(
    file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "correlation_heatmap_module_activity_scoring_global.png")),
    width = 3000,
    height = 2500,
    res = 300
)

draw(global_corr_plot)

dev.off()









