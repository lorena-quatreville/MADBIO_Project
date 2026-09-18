###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis
#

ANALYSIS_STEP_NAME = "03_Diff_compo_analysis"
PATH_ANALYSIS_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME)
PATH_ANALYSIS_EXTRA_OUTPUT = file.path( PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, "Extra")


# To load the merged Seurat Object that is the input file of this analysis step
                            
PREVIOUS_ANALYSIS_STEP_NAME = "02_Integration" # but we take the seurat object without integration  
                            
PATH_SEURAT_RDS_FILE = file.path( PATH_EXPERIMENT_OUTPUT,
                                  PREVIOUS_ANALYSIS_STEP_NAME,
                                 "Resolution_0.6",
                                  paste0( outputFilesPrefix, "seuratObject_Harmony.RDS"))
                                  
# Condition information
condition_info <- c("04_scRNASeq_A81_AS_tp2" = "Chronic", "08_scRNASeq_A82_S_tp2"= "Symptomatic",  "12_scRNASeq_A86_S_tp2" = "Symptomatic",  "13_scRNASeq_A62_AS_tp2" = "Chronic", "14_scRNASeq_A75_AS_tp2" = "Chronic", "15_scRNASeq_A76_S_tp2" = "Symptomatic", "16_scRNASeq_A83_AS_tp2" = "Chronic")

# Sequencing batch information
batch_info <- c("04_scRNASeq_A81_AS_tp2" = "batch_1", "08_scRNASeq_A82_S_tp2"= "batch_2", "12_scRNASeq_A86_S_tp2" = "batch_2",  "13_scRNASeq_A62_AS_tp2" = "batch_2", "14_scRNASeq_A75_AS_tp2" = "batch_2", "15_scRNASeq_A76_S_tp2" = "batch_2", "16_scRNASeq_A83_AS_tp2" = "batch_2")


sample_order <- c(
  "16_scRNASeq_A83_AS_tp2",
  "04_scRNASeq_A81_AS_tp2",
  "13_scRNASeq_A62_AS_tp2",
  "14_scRNASeq_A75_AS_tp2",
  "15_scRNASeq_A76_S_tp2",
  "12_scRNASeq_A86_S_tp2",
    "08_scRNASeq_A82_S_tp2" 
)

#### General

# Seed for pseudo-random numbers
SEED = 42;
set.seed(42)

# Number of cores to use when possible (for Seurat3 using 'future')
NBCORES = 4;

# Number of cells above which use ggplot instead of interactive plotly
PLOT_RASTER_NBCELLS_THRESHOLD = 20000;



#### Filtering / Normalization

# Switch from QC exploration mode and QC filtering
# In exploration mode, the cells that would be filtered by QC threshold are not filtered but marked
# They are shown in the later analysis with a different color
QC_EXPLORATION_MODE = TRUE

# Filters for loading seurat object
LOAD_MIN_CELLS     = 3;    # Retain cells with at least this many features (annotations)
LOAD_MIN_FEATURES  = 200;  # Retain annotations appearing in at least this many cells

# Cells with number of UMIs outside the range will be excluded
FILTER_UMI_MIN     = 0;
FILTER_UMI_MAX     = 40000;

# Cells with number of genes outside the range will be excluded
FILTER_FEATURE_MIN = 0;
FILTER_FEATURE_MAX = 6500;

# Cells with percentage of mitocohondrial genes above threshold will be excluded
FILTER_MITOPCT_MAX = 10;

# Cells with percentage of ribosomal genes below threshold will be excluded
FILTER_RIBOPCT_MIN = 0;

# Normalization parameters (see Seurat::NormalizeData())
DATA_NORM_METHOD      = "LogNormalize";
DATA_NORM_SCALEFACTOR = 10000;

# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE;
DATA_SCALE        = FALSE;
DATA_VARS_REGRESS = NULL;  # c("nCount_RNA") for UMIs (NULL to ignore)




#### Analysis parameters

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report

# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION     = 0.6;
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
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

# Parametre for DoubletFinder 
DELETE_MULTIPET              = FALSE
DOUBLET_FORMATION_RATE       = 0.055; # Assuming 8% doublet formation rate 
DOUBLETFINDER_USE_PCA_NBDIMS = 30;    # Number of dimensions to use from PCA results
pN=0.25;                              #  Default number of pN from DoubletFinder

# Dimensionality reduction parameters (TSNE/UMAP)
DIMREDUC_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE ;     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)
FINDMARKERS_MINPCT    = 0.1;      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25;     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001;    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP_TABLE     = 100;       # Number of marker genes to show in report tables (NULL for all)
FINDMARKERS_SHOWTOP_HEATMAP   = 5;       # Number of marker genes to show in repot heatmaps (NULL for all)

# Parameter for enrichment analysis in GO Terms
ENRICHMENT_GO_PVALUECUTOFF = 0.05
ENRICHMENT_GO_QVALUECUTOFF = 0.05



#### List of group of clusters of interest

CLUSTER_GROUP_LIST = list( T_NK_cell = c( 0, 3, 5, 7, 11, 14, 17, 19),
                			     Neutrophil_cell = c( 12),
                			     B_cell = c( 4,6,13,18),
                			     Myeloid_cell = c( 1,2, 8,9, 10, 16, 20, 15))


#### Lists of genes of interest

# Number of genes under which genes of a module (MODULES_GENES) are transfered to be analyzed individually (MONITORED_GENES)
MONITORED_GENES_SMALL_MODULE = 5;
MODULES_CONTROL_SIZE = 100;


## Contamination-related genes (txt file, one gene name by line)
#CONTAMINATION_GENES = readLines( file.path( PATH_PROJECT_EXTERNALDATA, "ContaminationGenes.txt"))
CONTAMINATION_GENES=NULL

## Genes monitored individually (tsv file, one column for each group of genes)
# MONITORED_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
#                                                   "03_MonitoredGenes",
#                                                   "Monitored.csv"),
#                                        sep = ",",
#                                        header = TRUE,
#                                        stringsAsFactors = FALSE,
#                                        row.names = NULL, fill = TRUE));
# MONITORED_GENES = Map('[', MONITORED_GENES, lapply(MONITORED_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
# 
# ## Genes monitored individually (tsv file, one column for each group of genes)
# CHOSEN_MARKER_GENES = read.table( file.path( PATH_EXPERIMENT_REFERENCE,
#                                                   "03_MonitoredGenes",
#                                                   "ChosenMarkerGenes.csv"),
#                                        sep = ",",
#                                        header = TRUE,
#                                        stringsAsFactors = FALSE,
#                                        row.names = NULL, fill = TRUE);
# CHOSEN_MARKER_GENES = CHOSEN_MARKER_GENES$Gene
# 
# ## Genes groups to analyze (csv file, one column for each group of genes)
# GROUP_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
#                                               "04_GeneGroups",
#                                               "GeneGroups.csv"),
#                                    sep = ",", quote = '"',
#                                    header = TRUE,
#                                    stringsAsFactors = FALSE,
#                                    row.names = NULL, fill = TRUE));
# GROUP_GENES = Map('[', GROUP_GENES, lapply(GROUP_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
# 
# 
# ## Genes monitored as modules (tsv file, one column for each group of genes)
# MODULES_GENES = as.list( read.table( file.path( PATH_EXPERIMENT_REFERENCE,
#                                                 "05_Modules",
#                                                 "Modules.csv"),
#                                         sep = ",",
#                                         header = TRUE,
#                                         stringsAsFactors = FALSE,
#                                         row.names = NULL, fill = TRUE));
# MODULES_GENES = Map('[', MODULES_GENES, lapply(MODULES_GENES, function(x){ which( nchar( x)>0)})); # Remove empty strings
# MODULES_GENES[[ "S_PHASE"]] = CELL_CYCLE_SPHASE_GENELIST
# MODULES_GENES[[ "G2M_PHASE"]] = CELL_CYCLE_G2MPHASE_GENELIST
# 
# 
