################################                                                        #######################################################
                                               # Dimensionality Reduction #
################################                                                        #######################################################


################################# Unsupervised (using Differential Gene Expression analysis between clusters)
    
data_CDS_DEG <- preprocess_cds(data_CDS_DEG, num_dim = 50) # performs PCA 
data_CDS_DEG <- reduce_dimension(data_CDS_DEG, reduction_method = "UMAP", max_components = 50)  # Neded to run learn_graph() afterwards
    
    
################################## Semi-supervised (based on SciGeneX modules)

data_CDS_SciGeneX <- preprocess_cds(data_CDS_SciGeneX, num_dim = 50) # performs PCA
data_CDS_SciGeneX <- reduce_dimension(data_CDS_SciGeneX, reduction_method = "PCA") # Gets PCA coordinates


data_CDS_SciGeneX_sub <- preprocess_cds(data_CDS_SciGeneX_sub, num_dim = 10) # performs PCA
data_CDS_SciGeneX_sub <- reduce_dimension(data_CDS_SciGeneX_sub, reduction_method = "PCA") # Gets PCA coordinates
