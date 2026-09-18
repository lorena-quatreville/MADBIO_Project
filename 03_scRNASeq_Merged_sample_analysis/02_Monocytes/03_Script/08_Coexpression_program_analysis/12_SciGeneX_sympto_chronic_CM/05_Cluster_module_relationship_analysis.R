# ##############################################
# This script analyse the SciGeneX gene modules
# ##############################################


## .................................................................
## Analyze the content of genes in each SciGeneX modules
## and the activity of each SciGeneX Modules
## .................................................................

## @knitr heterogeneity_SGXGeneModules
options(warn = 1) # initially on 2
# Creating list of module name
sgx_module_name_set <- paste0(
    "Module_",
    seq_along(scigx_res@gene_clusters)
)

sgx_module_feature_list <- setNames(
    scigx_res@gene_clusters,
    sgx_module_name_set
)

# Creating a module activity name list
module_activity_names <- paste0(
    "Activity.SGXModule",
    seq_along(scigx_res@gene_clusters)
)

# Creating the list of genes for each module
sgx_module_name_set <- paste0(
    "Module_",
    seq_along(scigx_res@gene_clusters)
)

sgx_module_feature_list <- setNames(
    scigx_res@gene_clusters,
    sgx_module_name_set
)





cat(" \n \n")
cat("### Individual module analysis {.tabset .tabset-fade}")
cat(" \n \n")

for( module_index in 1: length( sgx_module_name_set)){
  

  # Get the name of the module
  module_name = sgx_module_name_set[ module_index]
  
    # Create a directory
  module_dir <- file.path(PATH_ANALYSIS_OUTPUT, module_name)
    dir.create(module_dir, showWarnings = FALSE, recursive = TRUE)
    
    
  cat(" \n \n")  
  cat("#### ", module_name, " {.tabset .tabset-fade}")
  cat(" \n \n")
  
  cat(" \n \n")
  cat("##### Module top genes list")
  cat(" \n \n")
  
  # Get the list of module genes
  gene_ids = sgx_module_feature_list[[ module_name]]
  
  # Display the genes as a list

  cat("<div style='width:48%; display:inline-block; vertical-align:top;'>")
  print( table( sc10x@meta.data[[ module_activity_names[ module_index]]]) %>% kable( caption = paste( "Cell count of", module_name)) %>% kable_styling( full_width = FALSE))
  
  # Saving the count table 
  cell_count_df <- as.data.frame(
  table(sc10x@meta.data[[module_activity_names[module_index]]]))
    write.csv(
  cell_count_df,
  file.path(module_dir, "cell_counts.csv"),
  row.names = FALSE
    )
  cat("</div>")
    
  cat("<div style='width:48%; display:inline-block; vertical-align:top;'>")
  # df <- matrix( sort( gene_ids), ncol = 4, byrow = TRUE) |> as.data.frame()
  genes_sorted <- sort(gene_ids)

remainder <- length(genes_sorted) %% 4

if (remainder != 0) {
  genes_sorted <- c(
    genes_sorted,
    rep(NA, 4 - remainder)
  )
}

df <- as.data.frame(
  matrix(
    genes_sorted,
    ncol = 4,
    byrow = TRUE
  )
)
    write.csv(
  data.frame(gene = sort(gene_ids)),
  file.path(module_dir, "top_genes_long.csv"),
  row.names = FALSE
)


  print( kable(df, col.names = NULL, caption = paste( "Top genes of", module_name)) %>%  kable_styling(full_width = FALSE, position = "center", font_size = 12))
  cat("</div>")

  cat(" \n \n")  
  cat("##### UMAP of Module Activity")
  cat(" \n \n")
  
  print( DimPlot( sc10x, group.by = module_activity_names[ module_index]) + ggtitle( paste( "Activity of", module_name)))  
  # Saving
  p1 <- DimPlot(
  sc10x,
  group.by = module_activity_names[module_index]
) +
ggtitle(paste("Activity of", module_name))

ggsave(
  filename = file.path(module_dir, "UMAP_activity.png"),
  plot = p1,
  width = 8,
  height = 6
)

print(p1)


  print( DimPlot( sc10x, group.by = module_activity_names[ module_index], split.by = module_activity_names[ module_index]) + ggtitle( paste( "Activity of", module_name))) 
 # Saving
  p2 <- DimPlot(
  sc10x,
  group.by = module_activity_names[module_index],
  split.by = module_activity_names[module_index]
) +
ggtitle(paste("Activity of", module_name))

ggsave(
  filename = file.path(module_dir, "UMAP_activity_split.png"),
  plot = p2,
  width = 12,
  height = 6
)

print(p2)
  
  cat(" \n \n")  
  cat("##### UMAP of Module Score")
  cat(" \n \n")
  
  print( FeaturePlot( sc10x, feature = module_name) + 
           ggtitle( paste( "Score of", module_name)))  
  # Saving
  p3 <- FeaturePlot(
  sc10x,
  feature = module_name
) +
ggtitle(paste("Score of", module_name))

ggsave(
  filename = file.path(module_dir, "UMAP_score.png"),
  plot = p3,
  width = 8,
  height = 6
)

print(p3)
  
  cat(" \n \n")
  cat("##### Functional enrichment analysis")
  cat(" \n \n")
  
  # Compute GO BP functionnal enrichment
  # ...................................
  ego_bp <- enrichGO(gene          = gene_ids,
                     OrgDb         = org.Hs.eg.db,
                     keyType       = "SYMBOL",
                     ont           = "BP",        # Biological Process
                     pAdjustMethod = "BH",
                     pvalueCutoff  = 0.05,
                     qvalueCutoff  = 0.05,
                     readable      = FALSE)

  # Convert result to dataframe
  ego_df <- as.data.frame( ego_bp)
  
  # Check if some result is found
  if( nrow( ego_df) > 1){
  
    # Display result dotplot
    print( dotplot(ego_bp, showCategory = 10, title = "GO Biological Process Enrichment"))
    
    # Saving
    p_go <- dotplot(
  ego_bp,
  showCategory = 10,
  title = "GO Biological Process Enrichment"
    )

    ggsave(
  filename = file.path(module_dir, "GO_BP_dotplot.png"),
  plot = p_go,
  width = 10,
  height = 7
    )

    print(p_go)

    
    # Use REVIGO to reduce the GOTerms
    # ...................................
    # 0. Prepare GO IDs and scores (-log10 p.adjust)
    go_ids <- ego_df$ID
    scores <- setNames( -log10(ego_df$p.adjust), go_ids)
    
    # 1. Compute the semantic similarity matrix
    simMatrix <- calculateSimMatrix( go_ids,
                                     orgdb    = "org.Hs.eg.db",
                                     ont      = "BP",
                                     method   = "Rel")  # or "Resnik", "Lin", etc.
    
    # 2. Reduce redundant GO terms
    reduced <- reduceSimMatrix(simMatrix,
                               scores    = scores,
                               threshold = 0.7,               # similarity threshold (0–1)
                               orgdb     = "org.Hs.eg.db")
    
    reduced = reduced[ which( !duplicated( reduced$go)), ]
    
    # 3. Display the results as Treemap
    treemapPlot(reduced)
    # Saving
    png(
  file.path(module_dir, "REVIGO_treemap.png"),
  width = 2000,
  height = 1600,
  res = 300
    )

    treemapPlot(reduced)

    dev.off()

    cat(" \n \n")  
  }else{
    cat("<BR>No functional enrichment for this module.")
  }
}










# plot_sgx_modules_combined_scores function
plot_module_overlap <- function(
    seurat_obj,
    module1,
    module2,
    cluster_col = "heatmap_clusters"
){

    m1 <- as.character(seurat_obj[[module1]][,1])
    m2 <- as.character(seurat_obj[[module2]][,1])

    state <- rep(NA, length(m1))

    state[m1=="0" & m2=="0"] <- "Neither"
    state[m1=="1" & m2=="0"] <- paste0(module1, "_only")
    state[m1=="0" & m2=="1"] <- paste0(module2, "_only")
    state[m1=="1" & m2=="1"] <- "Both"

    df <- data.frame(
        cluster = seurat_obj[[cluster_col]][,1],
        state = state
    )

    df_sum <- df %>%
        count(cluster, state) %>%
        group_by(cluster) %>%
        mutate(freq = n/sum(n))

    mat <- df_sum %>%
        select(cluster, state, freq) %>%
        pivot_wider(
            names_from = cluster,
            values_from = freq,
            values_fill = 0
        )

    mat2 <- as.matrix(mat[,-1])
    rownames(mat2) <- mat$state

    pheatmap(
        mat2,
        cluster_rows = FALSE,
        cluster_cols = FALSE,
        main = paste(module1, "vs", module2)
    )
}
png(
    filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix,"_module_comb_2_4.png")),
    width = 1800,
    height = 1200,
    res = 300
)

plot_module_overlap(
    sc10x,
    module1 = "Activity.SGXModule2",
    module2 = "Activity.SGXModule4"
)

dev.off()





cat("\n \n")


## ..........................................................................
## Analyze the relationship between clusters and sciegenex module activity
## ..........................................................................

## @knitr heterogeneity_SGXGeneModules_vs_clusters

cat(" \n \n")
cat("### Percentage of cells activating module in clusters")
cat(" \n \n")

# Loop over the cluster to study them one by one
Idents( sc10x) = "heatmap_clusters"
clusters_modules_activity = data.frame()
for( cluster in levels( Idents( sc10x))){

  # Get the list of cells in the cluster
  cluster_cells = names( Idents(sc10x))[ which( Idents(sc10x) == cluster)]
  cat("There are this number of cells in this cluster", length(cluster_cells))
  # For each SGX module, compute the number of cells in activity in the cluster
  for( module_index in 1:length( sgx_module_name_set)){
    
    module_activity_name = module_activity_names[ module_index]
    cluster_module_activity_table = table( sc10x@meta.data[ module_activity_name][ cluster_cells,])
    percentage_activity = cluster_module_activity_table[ "1"]/( cluster_module_activity_table[ "0"] + cluster_module_activity_table[ "1"])
    clusters_modules_activity = rbind( clusters_modules_activity, data.frame( cluster = cluster,
                                                                              module = sgx_module_name_set[ module_index],
                                                                              activity.pct = percentage_activity
                                                                              ))
  }
}

clusters_modules_activity$module = factor( clusters_modules_activity$module, levels = sgx_module_name_set)
clusters_modules_activity$cluster = factor( clusters_modules_activity$cluster, levels = levels( Idents( sc10x)))
clusters_modules_activity_wide = spread(clusters_modules_activity, cluster, activity.pct)
row.names( clusters_modules_activity_wide) = clusters_modules_activity_wide$module
clusters_modules_activity_wide$module = NULL

# heatmap of module activity by cluster
print( pheatmap( clusters_modules_activity_wide, 
          cluster_rows = FALSE, 
          cluster_cols = FALSE, 
          cellwidth = 20, 
          cellheight = 20,
          colorRampPalette(brewer.pal(n = 7, name = "YlOrRd"))(100))
       )

# Saving
# Open device
png(
    file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "module_activity_heatmap_by_cluster.png")),
    width = 2500,
    height = 1800,
    res = 300
)

# Draw heatmap
pheatmap(
    clusters_modules_activity_wide,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    cellwidth = 20,
    cellheight = 20,
    color = colorRampPalette(
        brewer.pal(n = 7, name = "YlOrRd")
    )(100)
)

# Close device
dev.off()

# -----------------------------
# Save PDF
# -----------------------------
pdf(
  file.path(PATH_ANALYSIS_OUTPUT,
            paste0(outputFilesPrefix, "module_activity_heatmap_by_cluster.pdf")),
  width = 10,
  height = 8
)

pheatmap(
  clusters_modules_activity_wide,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  cellwidth = 20,
  cellheight = 20,
  color = colorRampPalette(
    brewer.pal(n = 7, name = "YlOrRd")
  )(100)
)

dev.off()


# If I want to know for a specific cluster and 2 spcific modules if the same cells are expressing both modules
Idents(sc10x) <- "heatmap_clusters"

cluster1_cells <- WhichCells(sc10x, idents = "1")

m1 <- sc10x$Activity.SGXModule1[cluster1_cells]
m2 <- sc10x$Activity.SGXModule2[cluster1_cells]

table(m1, m2)



# Checking for cell composition between clusters (individuals and condition)

# Getting the percentage of cells belongign to each individual in each cluster
# Extract metadata
metadata <- sc10x@meta.data

# Count cells per cluster and individual
cluster_composition_indiv <- metadata %>%
  group_by(heatmap_clusters, orig.ident) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  ungroup()

# Calculate percentage of each individual within each cluster
cluster_composition_indiv <- cluster_composition_indiv %>%
  group_by(heatmap_clusters) %>%
  mutate(percent = 100 * cell_count / sum(cell_count)) %>%
  select(heatmap_clusters, orig.ident, percent) %>%
  pivot_wider(
    names_from = orig.ident,
    values_from = percent,
    values_fill = 0
  )

# Save
write.table(
  cluster_composition,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix, "_cell_cluster_composition.csv")
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)



