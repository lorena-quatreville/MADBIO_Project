###############################################################################
# This file defines PROJECT parameters as global variables that will be loaded
# before analysis starts. It should define common parameters shared by several
# samples. 
#


#### General

GLOBAL_DESCRIPTION = "MADBIO_Project"

# SCIENTIFIC_GROUP = "BMMLAB"
SCIENTIFIC_PROJECT_NAME = "01_MADBIO_Cameroon"
SCIENTIFIC_SUBPROJECT_NAME = "01_scRNASeq"
EXPERIMENT_NAME = "scRNASeq_S3"

#### Additional custom variables (for a specific analysis step, tools, ...)

#SAMPLE_ID = ""



#### Input / Output (PROJECT Folder)

# Output folder name in data folder (for R session object, lists of cells/genes) 
PATH_PROJECT = file.path( "/home/quatreville", 
                                        # SCIENTIFIC_GROUP,
                                        # "BIOINFO", 
                                        "03_BIOINFO_PROJECT",
                                        SCIENTIFIC_PROJECT_NAME,
                                        SCIENTIFIC_SUBPROJECT_NAME)

PATH_EXPERIMENT = file.path( PATH_PROJECT, EXPERIMENT_NAME)

PATH_EXPERIMENT_RAWDATA       = file.path( PATH_EXPERIMENT, "00_InputData")
PATH_EXPERIMENT_REFERENCE     = file.path( PATH_EXPERIMENT, "01_Reference")
PATH_EXPERIMENT_OUTPUT        = file.path( PATH_EXPERIMENT, "05_Output")



# Create a 'safe' unique prefix for output files
outputFilesPrefix = paste0( SCIENTIFIC_PROJECT_NAME, "_", SCIENTIFIC_SUBPROJECT_NAME, "_",
                            EXPERIMENT_NAME, "_"
                            )



#### Debug

.SHOWFLEXBORDERS = FALSE;
.VERBOSE = FALSE;



