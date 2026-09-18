################################                                                        #######################################################
                                         # Ordering cells along the trajectory #
################################                                                        #######################################################



# ...............sc10x from GlobalHeterogeneity step ........................ 

# Trajectory inference step
data_CDS <- learn_graph(data_CDS, use_partition = FALSE) # use partitions calculated during cluster_cells and therefore to learn disjoint graph in each partition. When use_partition = FALSE, a single graph is learned across all partitions. Default is TRUE.

# Ordering cells on the trajectory
data_CDS <- order_cells(data_CDS) # ordering cells along the trajectory

# Ploting the trajectory tree
trajectory_tree_CDS <- plot_cells(data_CDS, color_cells_by = "partition", cell_size = 1, group_label_size = 5)

dir.create(file.path(PATH_ANALYSIS_OUTPUT))
dir.create(file.path(PATH_ANALYSIS_OUTPUT, cellType))

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, cellType, "_Trajectory_Tree_before_SciGeneX.png")), plot = trajectory_tree_CDS, device ="png", width = 23, height = 13, units ="cm")


# ...............sc10x_SciGeneX from the SciGeneX analysis ..................


############# All SciGeneX genes present in modules

# Trajectory inference step
data_CDS_SciGeneX <- learn_graph(data_CDS_SciGeneX, use_partition = FALSE) # use partitions calculated during cluster_cells and therefore to learn disjoint graph in each partition. When use_partition = FALSE, a single graph is learned across all partitions. Default is TRUE.

# When ordering cells, specify your ordering genes:
data_CDS_SciGeneX <- order_cells(data_CDS_SciGeneX) # ordering cells along the trajectory

# Ploting the trajectory tree
trajectory_tree_CDS_SciGeneX <- plot_cells(data_CDS_SciGeneX, color_cells_by = "partition", cell_size = 1, group_label_size = 5)

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, cellType, paste0( outputFilesPrefix, cellType, "_Trajectory_Tree_after_SciGeneX.png")), plot = trajectory_tree_CDS_SciGeneX, device ="png", width = 23, height = 13, units ="cm")








