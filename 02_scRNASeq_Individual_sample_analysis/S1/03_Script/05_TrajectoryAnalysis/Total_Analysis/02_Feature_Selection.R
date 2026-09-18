################################                                                        #######################################################
                                 # Informativegene selection for trajectory analysis #
################################                                                        #######################################################




################################# Unsupervised (using Differential Gene Expression analysis between clusters)

# Use DE genes between clusters as ordering genes 
DEG_markers <- FindAllMarkers(sc10x_subtype_SciGeneX, only.pos = FALSE) # both down and up-regulated along the process are informative

# Select top 50 most differentially expressed genes for each seurat cluster and removes duplicates
DEG_ordering_genes <- DEG_markers %>%
  filter(p_val_adj < 0.05) %>%  # Keep only significant DE genes
  group_by(cluster) %>%
  top_n(50, abs(avg_log2FC)) %>%  # Keep top up or down regulated
  pull(gene) %>%
  unique()

# Subseting the cds to only the ordering genes
data_CDS_DEG <- data_CDS[DEG_ordering_genes, ]

# To visualize results
# Select 10 top genes for each cluster
top_DEG_ordering_genes <- DEG_markers %>%
  filter(p_val_adj < 0.05) %>%  # Keep only significant DE genes
  group_by(cluster) %>%
  top_n(5, avg_log2FC) %>%  # Keep top up regulated
  pull(gene) %>%
  unique()
  
data_CDS_top_DEG <- data_CDS[top_DEG_ordering_genes, ]  

plot_genes_by_group(
    cds = data_CDS_top_DEG,
    markers = top_DEG_ordering_genes,
    group_cells_by = "cluster"
 ) # maybe using an already defined variable for the number of top genes( cf. FINDMARKERS_SHOWTOP_HEATMAP)
 
 # same with a heatmap (https://cole-trapnell-lab.github.io/monocle3/docs/differential/#gene-modules) (more to do a combination of modules)




################################## Semi-supervised (based on SciGeneX modules)

# All SciGeneX genes present in modules
SciGeneX_ordering_genes <- row_names(res_subtype_SciGeneX) 

# Subseting the cds to only the ordering genes
data_CDS_SciGeneX <- data_CDS[SciGeneX_ordering_genes, ]

# To visualize results? Creating a dataframe with a list of genes per module


# A subtype of genes based on specific modules
module_nbs = c(2, 3, 4, 5, 6, 14, 56)
SciGeneX_ordering_genes_sub <- row_names(res_subtype_SciGeneX[module_nbs,])

# Subseting the cds to only the ordering genes
data_CDS_SciGeneX_sub <- data_CDS[SciGeneX_ordering_genes_sub, ]


