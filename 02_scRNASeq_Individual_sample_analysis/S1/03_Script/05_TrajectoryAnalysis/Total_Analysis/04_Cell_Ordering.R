################################                                                        #######################################################
                                         # Ordering cells along the trajectory #
################################                                                        #######################################################



################################# Unsupervised (using Differential Gene Expression analysis between clusters)

# Orering the cells along the trajectory
data_CDS_DEG <- learn_graph(data_CDS_DEG, use_partition = FALSE) # trajectory inference step # use partitions calculated during cluster_cells and therefore to learn disjoint graph in each partition. When use_partition = FALSE, a single graph is learned across all partitions. Default is TRUE.
data_CDS_DEG <- order_cells(data_CDS_DEG) # ordering cells along the trajectory

# Ploting the trajectory tree
plot_cell_trajectory(data_CDS_DEG, color_by = "State")
plot_cells(data_CDS_DEG, color_cells_by = "partition")
plot_cells(data_CDS_DEG, color_cells_by = "partition") +
  facet_wrap(~partition, nrow = 1)

# Same but splitted by state
plot_cell_trajectory(data_CDS_DEG, color_by = "State") +
    facet_wrap(~State, nrow = 1)
    

# Ploting cells according to expression of genes of interest
gene_of_interest <- c()
meaningful_genes <- row.names(subset(fData(data_CDS_DEG ),gene_short_name %in% gene_of_interest))

# Expression vs state graph
plot_genes_jitter(data_CDS_DEG[meaningful_genes,], grouping = "State", min_expr = 0.1)




data_CDS <- learn_graph(data_CDS, use_partition = FALSE)
data_CDS <- order_cells(data_CDS)




# When ordering cells, specify your ordering genes:
data_CDS <- order_cells(data_CDS, ordering_genes = DEG_ordering_genes)




