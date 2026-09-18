######### Data merging and normalization


# READ DATA
# ---------

## @knitr loadData

# Increase power
options(future.globals.maxSize = 16 * 1024^3)  # 16GB if needed
getOption("future.globals.maxSize")


sc10x_final <- NULL  # Initialize before the loop


# For each seurat object loaded (global per sample) from GlobalHeterogeneity
for (n in 1:nb_datasets) {
    PATH_SEURAT_RDS_FILE_GLOBAL <- get(paste0("PATH_SEURAT_RDS_FILE_", n))
    print(PATH_SEURAT_RDS_FILE_GLOBAL)
    sample_OutputFilesPrefix = get(paste0("OutputFilesPrefix_", n))
    print(sample_OutputFilesPrefix)
    sample_path_experiment_output = get(paste0("PATH_EXPERIMENT_OUTPUT_", n))
    sc10x <- readRDS(PATH_SEURAT_RDS_FILE_GLOBAL)
    
    # Remove line breaks in cluster names
    sc10x$cluster_renames <- gsub("\n", " ", sc10x$cluster_renames)
    sc10x$cluster_renames <- trimws(sc10x$cluster_renames)  # Remove leading/trailing spaces
    sc10x$cluster_renames <- gsub("\\s+", " ", sc10x$cluster_renames)  # Replace multiple spaces with single space
    sc10x$cluster_renames <- as.factor(sc10x$cluster_renames)  # Ensure it's a factor

    
    # Seurat Object loaded for each cell subtype
   # Initialize sgx_clusters as a copy of the global cluster_renames
    sc10x$sgx_clusters <- as.character(sc10x$cluster_renames)

    
    for (ct in cell_types) {
     # Load the cell-type specific Seurat object (SciGeneX)
     PATH_CT <- get(paste0("PATH_SGX_RDS_FILE_", ct, "_", n))
     print(PATH_CT)
     sc10x_ct <- readRDS(PATH_CT)
     
    # Extracting cell identity and cluster_renames label for each sub-dataset
    # Get the refined cluster labels from the SciGeneX object
    sgx_clusters <- as.character(sc10x_ct$cluster_renames) 
names(sgx_clusters) <- colnames(sc10x_ct)  # now each label has a cell barcode

    # Map these labels to the global object by name
sc10x$sgx_clusters[names(sgx_clusters)] <- sgx_clusters
  }
    
    # Remove line breaks in cluster names
    sc10x$sgx_clusters <- gsub("\n", " ", sc10x$sgx_clusters)
    sc10x$sgx_clusters <- trimws(sc10x$sgx_clusters)  # Remove leading/trailing spaces
    sc10x$sgx_clusters <- gsub("\\s+", " ", sc10x$sgx_clusters)  # Replace multiple spaces with single space
    sc10x$sgx_clusters <- as.factor(sc10x$sgx_clusters)  # Ensure it's a factor
   
    # Get all unique labels in the global object, including new ones
    all_labels <- unique(sc10x$sgx_clusters)
    sc10x$sgx_clusters <- factor(sc10x$sgx_clusters, levels = all_labels)

    # Set the Idents of the global object to the refined labels
  Idents(sc10x) <- sc10x$sgx_clusters
  
  # Optional: visualize
  umap_sample <- DimPlot(sc10x, reduction = "umap", group.by = "sgx_clusters", label = TRUE,       label.size = 4)
 print(umap_sample)
    ggsave(filename = file.path(sample_path_experiment_output, SEURAT_ANALYSIS_STEP_NAME,
    "Resolution_0.6", paste0(sample_OutputFilesPrefix, "_umap_refined.png")), plot = umap_sample, device ="png", width = 23, height = 13, units ="cm")
    
  # Saving the new version of the sc10x object
  print(sample_OutputFilesPrefix)
  saveRDS(sc10x, file = file.path(sample_path_experiment_output, SEURAT_ANALYSIS_STEP_NAME,
    "Resolution_0.6", paste0(sample_OutputFilesPrefix, "seuratObject_refined.RDS")))
 

  # Adding this dataset to the list of datasets
  if (is.null(sc10x_final)) { 
      sc10x_final <- sc10x  
    } else {
      sc10x_final <- merge(sc10x_final, y = sc10x)
    }
  }
sc10x <- sc10x_final
rm(sc10x_final)
gc()


# NORMALIZE DATA
# --------------

## @knitr normalizeData

sc10x = NormalizeData( object = sc10x ,
                       normalization.method = DATA_NORM_METHOD,
                       scale.factor = DATA_NORM_SCALEFACTOR,
                       verbose = .VERBOSE)

gc()

# Scaling of the sc10x object
sc10x <- ScaleData(object = sc10x, do.center = DATA_CENTER,
                   do.scale  = DATA_SCALE,
                   vars.to.regress = DATA_VARS_REGRESS,
                   verbose = .VERBOSE)
gc()





