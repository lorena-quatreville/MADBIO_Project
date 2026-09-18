

############################################################
################### LOAD MODULES ###########################
############################################################

gene_modules_df <- read.csv(file= file.path(PATH_EXPERIMENT_OUTPUT,
                                  PREVIOUS_ANALYSIS_STEP_NAME,
                                 cellType, paste0(outputFilesPrefix, "_ModuleGenes_list.csv")),   stringsAsFactors = FALSE)



############################################################
############ IDENTIFY IL1B / ICAM1 MODULE ##################
############################################################

# Find modules containing BOTH IL1B and ICAM1

selected_module <- NULL
gene_modules <- gene_modules_df[,2]

for(i in seq_along(gene_modules)) {

  # Split genes in module
  genes <- unlist(strsplit(as.character(gene_modules[i]), split = ","))

  # Remove spaces
  genes <- trimws(genes)

  # Check if BOTH genes are present
  if(all(c("IL1B") %in% genes)) {

    selected_module <- genes

    cat("Module found at row:", i, "\n")
    cat("Number of genes:", length(genes), "\n")

    break
  }
}

# Inspect module
head(selected_module)
length(selected_module)


############################################################
################### ADD MODULE SCORE #######################
############################################################

sc10x <- AddModuleScore(
  object = sc10x,
  features = list(selected_module),
  name = "activation_module"
)


############################################################
################### VISUALIZE UMAP #########################
############################################################


FeaturePlot(
  sc10x,
  features = "activation_module1",
  reduction = "umap"
)


############################################################
################ TRANSFER TO MONOCLE #######################
############################################################

data_CDS_SciGeneX@colData$activation_module1 <-
  sc10x@meta.data[
    colnames(data_CDS_SciGeneX),
    "activation_module1"
  ]

############################################################
################ TRAJECTORY VISUALIZATION ##################
############################################################
# Trajectory inference step
data_CDS_SciGeneX <- learn_graph(data_CDS_SciGeneX, use_partition = FALSE) # use partitions calculated during cluster_cells and therefore to learn disjoint graph in each partition. When use_partition = FALSE, a single graph is learned across all partitions. Default is TRUE.

trajectory_activation <- plot_cells(
  data_CDS_SciGeneX,
  color_cells_by = "activation_module1",
  cell_size = 1,
  label_cell_groups = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE
)

print(trajectory_activation)

ggsave(
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    cellType,
    paste0(
      outputFilesPrefix,
      cellType,
      "_ActivationModule_Trajectory.png"
    )
  ),
  plot = trajectory_activation,
  device = "png",
  width = 23,
  height = 13,
  units = "cm"
)

############################################################
################ PSEUDOTIME DYNAMICS #######################
############################################################


# Learn graph
data_CDS_SciGeneX <- learn_graph(
    data_CDS_SciGeneX,
    use_partition = FALSE
)

# Order cells = compute pseudotime
data_CDS_SciGeneX <- order_cells(
    data_CDS_SciGeneX,
    reduction_method = "UMAP"
)



df <- data.frame(
  pseudotime = pseudotime(data_CDS_SciGeneX),
  module_score = colData(data_CDS_SciGeneX)$activation_module1
)

library(ggplot2)

p_pseudotime <- ggplot(
  df,
  aes(x = pseudotime, y = module_score)
) +
  geom_point(alpha = 0.2) +
  geom_smooth(
    method = "loess",
    span = 0.3,
    color = "red"
  ) +
  theme_classic() +
  labs(
    title = "Activation module dynamics along pseudotime",
    x = "Pseudotime",
    y = "Activation module score"
  )

print(p_pseudotime)

ggsave(
  filename = file.path(
    PATH_ANALYSIS_OUTPUT,
    cellType,
    paste0(
      outputFilesPrefix,
      cellType,
      "_ActivationModule_Pseudotime.png"
    )
  ),
  plot = p_pseudotime,
  device = "png",
  width = 12,
  height = 8
)

  
  
