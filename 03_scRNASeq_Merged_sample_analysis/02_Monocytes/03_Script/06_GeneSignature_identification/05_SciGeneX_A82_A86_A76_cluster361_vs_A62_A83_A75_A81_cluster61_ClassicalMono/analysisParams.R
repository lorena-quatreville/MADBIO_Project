

nb_datasets = 7
sample_id_A = c("08_scRNASeq_A82_S_tp2", "15_scRNASeq_A76_S_tp2", "12_scRNASeq_A86_S_tp2") # sample ID we want to compare (here symptomatic individuals)
sample_id_B = c("04_scRNASeq_A81_AS_tp2","13_scRNASeq_A62_AS_tp2","14_scRNASeq_A75_AS_tp2", "16_scRNASeq_A83_AS_tp2") # sample ID we want to compare to (here chronic individual)

# To create a folder per comparison
ANALYSIS_NAME = "05_SciGeneX_A82_A86_A76_cluster361_vs_A62_A83_A75_A81_cluster61_ClassicalMono"
# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_", SCIENTIFIC_SUBPROJECT_NAME, "_",
                            EXPERIMENT_NAME, "_", ANALYSIS_NAME, "_"
                            )
PATH_EXPERIMENT_OUTPUT        = file.path( PATH_EXPERIMENT, "05_Output")


ANALYSIS_STEP_NAME = "06_GeneSignature_identification"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME)
PATH_ANALYSIS_EXTRA_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME, "Extra")


# To load the merged Seurat Object that is the input file of this analysis step
                            
PREVIOUS_ANALYSIS_STEP_NAME = "02_Integration" # but we take the seuratObject without Harmony bacause we are not comparing clusters
previous_outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_", SCIENTIFIC_SUBPROJECT_NAME, "_",
                            EXPERIMENT_NAME, "_"
                            )                  
 # To load the seuratObject
previous_outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_", SCIENTIFIC_SUBPROJECT_NAME, "_",
                            EXPERIMENT_NAME, "_"
                            )                  
PATH_SEURAT_RDS_FILE = file.path( PATH_EXPERIMENT_OUTPUT,
                                  PREVIOUS_ANALYSIS_STEP_NAME,
                                 "Resolution_0.6",
                                  paste0( previous_outputFilesPrefix, "seuratObject_Harmony.RDS"))                        

# Condition information
condition_info <- c("04_scRNASeq_A81_AS_tp2" = "Chronic", "08_scRNASeq_A82_S_tp2"= "Symptomatic",  "12_scRNASeq_A86_S_tp2" = "Symptomatic",  "13_scRNASeq_A62_AS_tp2" = "Chronic", "14_scRNASeq_A75_AS_tp2" = "Chronic", "15_scRNASeq_A76_S_tp2" = "Symptomatic", "16_scRNASeq_A83_AS_tp2" = "Chronic")

# Sequencing batch information
batch_info <- c("04_scRNASeq_A81_AS_tp2" = "batch_1", "08_scRNASeq_A82_S_tp2"= "batch_2", "12_scRNASeq_A86_S_tp2" = "batch_2",  "13_scRNASeq_A62_AS_tp2" = "batch_2", "14_scRNASeq_A75_AS_tp2" = "batch_2", "15_scRNASeq_A76_S_tp2" = "batch_2", "16_scRNASeq_A83_AS_tp2" = "batch_2")                                 
















#### General

# Seed for pseudo-random numbers
SEED = 42;

# Number of cores to use when possible (for Seurat3 using 'future')
NBCORES = 4;

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 20000;


#### Analysis parameters

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report

# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION     = 0.6;
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
nbPC_findclusters = FINDCLUSTERS_USE_PCA_NBDIMS
FINDCLUSTERS_ALGORITHM      = 4;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden

# RENAME_CLUSTERS = list()
# RENAME_CLUSTERS[[ SAMPLE_ALL]] = list()
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "0"]] = 0
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "1"]] = 5
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "2"]] = 1
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "3"]] = 3
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "4"]] = 6
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "5"]] = 2
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "6"]] = 4
# RENAME_CLUSTERS[[ SAMPLE_ALL]][[ "7"]] = 7
  

# PCA parameters
PCA_NPC              = 50;  # Default number of dimensions to use for PCA (see Seurat::RunPCA())
PCA_PLOTS_NBDIMS     = 3;   # Number of dimensions to show in PCA-related plots
PCA_PLOTS_NBFEATURES = 15;  # Number of'top' features to show when plotting PCA loadings

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
useReduction = "umap" # the reduction method to use in DimPlot()

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP_TABLE     = 100;       # Number of marker genes to show in report tables (NULL for all)
FINDMARKERS_SHOWTOP_HEATMAP   = 5;       # Number of marker genes to show in repot heatmaps (NULL for all)

# Parameter for enrichment analysis in GO Terms
ENRICHMENT_GO_PVALUECUTOFF = 0.05
ENRICHMENT_GO_QVALUECUTOFF = 0.05


