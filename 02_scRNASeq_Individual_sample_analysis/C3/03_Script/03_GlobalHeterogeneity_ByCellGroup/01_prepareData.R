# #########################################
# This script reads and filters sc10x  data
# #########################################




# READ DATA
# ---------

## @knitr loadData

# Create the seurat object (RNA only)
sc10x = readRDS( PATH_SEURAT_RDS_FILE)

cat("<BR><b>Loading data from: ", PATH_SEURAT_RDS_FILE, "</b>")

# Select the cell in the chosen cluster group 
cat("<BR><p style='color:red;'><b>WARNING: Subsetting cells to keep only those in cluster group", CLUSTER_GROUP, "(original clusters", paste( CLUSTER_GROUP_LIST[[ CLUSTER_GROUP]], collapse = ","), ")</b></p>")

Idents( sc10x) = "seurat_clusters"
cells_in_group <- WhichCells( sc10x, idents = CLUSTER_GROUP_LIST[[ CLUSTER_GROUP]])

# Highlight the cells in UMAP
print(
  DimPlot( sc10x, 
           label=TRUE, label.size = 10,
           cells.highlight= list(cells_in_group), 
           cols.highlight = "red", 
           cols= "grey") +
    ggtitle( paste( 'UMAP of Cells in', CLUSTER_GROUP)) +
    theme( axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           legend.position = "none",
           plot.margin = margin( 0, 0, 0, 0, "cm"),
           plot.title = element_text( face = "bold",
                                      size = rel( 16/14),
                                      hjust = 0.5,
                                      vjust = 1,
                                      margin = margin( b = 7)))
)

# Subset the Seurat object
sc10x = subset( sc10x, cells = cells_in_group )

# Memorize the original seurat clusters in new meta-data
sc10x = AddMetaData( sc10x, metadata = Idents( sc10x), col.name = "original_seurat_clusters")


# ### Remove eventual NULL (empty) list elements from list of genes to be monitored
# monitoredGroupEmpty = sapply( MONITORED_GENES, is.null);
# if(any( monitoredGroupEmpty)) warning( paste("Following group(s) of genes to be monitored will be ignored because empty:", paste( names(monitoredGroupEmpty)[monitoredGroupEmpty], collapse=" - ")));
# MONITORED_GENES = MONITORED_GENES[! monitoredGroupEmpty];
# 
# # Check whether genes in MONITORED_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
# #matchMonitoredGenes = match( toupper( unlist( MONITORED_GENES)), toupper( rownames( GetAssayData( sc10x))));
# matchMonitoredGenes = match( ( unlist( MONITORED_GENES)), ( rownames( GetAssayData( sc10x))));
# monitoredGenesNotFound = unique( unlist( MONITORED_GENES)[is.na( matchMonitoredGenes)]);
# if(any( is.na( matchMonitoredGenes))) warning( paste( "Following gene(s) to be monitored could not be found in experimental data and will be ignored:", paste( paste0("'", monitoredGenesNotFound, "'"), collapse=" - ")));
# # Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
# MONITORED_GENES = relist( rownames( GetAssayData( sc10x))[ matchMonitoredGenes ], skeleton = MONITORED_GENES); # Does not work with NULL list elements (removed earlier)
# # Finally remove names that did not match
# MONITORED_GENES = lapply( MONITORED_GENES, na.omit);
# 
# 
# ### Remove eventual NULL (empty) list elements from list of genes in modules
# modulesGroupEmpty = sapply( MODULES_GENES, is.null);
# if(any( modulesGroupEmpty)) warning( paste("Following module(s) of genes will be ignored because empty:", paste( names(modulesGroupEmpty)[modulesGroupEmpty], collapse=" - ")));
# MODULES_GENES = MODULES_GENES[! modulesGroupEmpty];
# 
# # Check whether genes in MODULES_GENES are actually found in assay object (case mismatch is not corrected anymore to avoid genes names confusion between species)
# #matchModulesGenes = match( toupper( unlist( MODULES_GENES)), toupper( rownames( GetAssayData( sc10x))));
# matchModulesGenes = match( ( unlist( MODULES_GENES)), ( rownames( GetAssayData( sc10x))));
# modulesGenesNotFound = unique( unlist( MODULES_GENES)[is.na( matchModulesGenes)]);
# if(any( is.na( matchModulesGenes))) warning( paste( "Following gene(s) from modules list could not be found in experimental data and will be ignored:", paste( paste0("'", modulesGenesNotFound, "'"), collapse=" - ")));
# # Replace names by names found in assay (correcting eventual case mismatch, NA for not found)
# MODULES_GENES = relist( rownames( GetAssayData( sc10x))[ matchModulesGenes ], skeleton = MODULES_GENES); # Does not work with NULL list elements (removed earlier)
# # Finally remove names that did not match
# MODULES_GENES = lapply( MODULES_GENES, na.omit);
# 
# 
# ### Transfer genes in very small modules (MODULES_GENES) to be analyzed individually (MONITORED_GENES)
# modulesToTransfer = sapply(MODULES_GENES, length) < MONITORED_GENES_SMALL_MODULE
# if(any(modulesToTransfer))
# {
#     warning( paste0( "Following Module(s) contained very few genes (<", MONITORED_GENES_SMALL_MODULE, "). These genes were transfered to 'Monitored genes' to be analyzed individually: ", paste( names(modulesToTransfer)[modulesToTransfer], collapse = " - ")));
#     # Alter name so they can be recognized in list of Monitored genes
#     names( MODULES_GENES)[modulesToTransfer] = paste0( "MOD_", names( MODULES_GENES)[modulesToTransfer]);
#     MONITORED_GENES = c( MONITORED_GENES, MODULES_GENES[modulesToTransfer]);
#     MODULES_GENES = MODULES_GENES[!modulesToTransfer];
# }
# 
# # ...................
# # Heat Shock data
# # ...................
# 
# # Load the list of 512 genes associated with heat shock and stress in article from Campbell et al.
# hs_stress_genes_df = read.table( file = PATH_HS_STRESS_MARKER_GENES_TABLE_FILE, header = TRUE, sep=",", quote = NULL)
# 
# # Convert list of genes from human symbol to mouse symbol
# converted_genes_df = convertHumanGeneList( hs_stress_genes_df$gene_symbol)
# hs_stress_genes_df = merge( hs_stress_genes_df, converted_genes_df, by.x = "gene_symbol", by.y="human_symbol")
# hs_stress_genes_df = hs_stress_genes_df[ order( hs_stress_genes_df$PValue, hs_stress_genes_df$logFC, decreasing = c(FALSE,TRUE)), ]
# 
# DT::datatable( hs_stress_genes_df, caption = "List of HS and Stress associated genes in human")
# 
# MODULES_GENES[[ "HeatShock"]] = unique( hs_stress_genes_df[ , "symbol"])
# MODULES_GENES[[ "HeatShock_Top40"]] = unique( hs_stress_genes_df[ 1:40, "symbol"])



