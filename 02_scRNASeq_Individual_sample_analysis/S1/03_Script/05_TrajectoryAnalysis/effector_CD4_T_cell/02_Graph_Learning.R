################################                                                        #######################################################
                                         # Ordering cells along the trajectory #
################################                                                        #######################################################



# ...............sc10x from GlobalHeterogeneity step ........................ 

# Trajectory inference step
data_CDS <- learn_graph(data_CDS, use_partition = FALSE) # use partitions calculated during cluster_cells and therefore to learn disjoint graph in each partition. When use_partition = FALSE, a single graph is learned across all partitions. Default is TRUE.

# When ordering cells, specify your ordering genes:
data_CDS <- order_cells(data_CDS) # ordering cells along the trajectory

# Ploting the trajectory tree
plot_cells(data_CDS, color_cells_by = "partition")




# ...............sc10x_SciGeneX from the SciGeneX analysis ..................


############# All SciGeneX genes present in modules

# Trajectory inference step
data_CDS_SciGeneX <- learn_graph(data_CDS_SciGeneX, use_partition = FALSE) # use partitions calculated during cluster_cells and therefore to learn disjoint graph in each partition. When use_partition = FALSE, a single graph is learned across all partitions. Default is TRUE.

# When ordering cells, specify your ordering genes:
data_CDS_SciGeneX <- order_cells(data_CDS_SciGeneX) # ordering cells along the trajectory

# Ploting the trajectory tree
plot_cells(data_CDS_SciGeneX, color_cells_by = "partition")










