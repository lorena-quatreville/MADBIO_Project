

Script for analyzing samples C1 to C4 and S1 to S3 :

01_Demultiplexing, 02_scRNASeq_Individual_sample_analysis and 03_scRNASeq_Merged_sample_analysis are for analyzing single-cell RNA-Seq data
04_Cytometry is for analizing Cytometry data
00_Container contains all Dockerfile and readme used to create container with pinned package versions

00_Container contains Docker files and readme files for three different containers:
    - Container_CellRanger: Dockerfile and readme for creating the Docker container with pinned package versions used for sample demultiplexing
    - rstudio_seurat_v4.2.2 is the container used for all analysis except Trajectory analysis and Compositional analysis
    - rstudio_statistics_v4.4.0 is the container used for Trajectory analysis (Monocle3) and Compositional analysis (DCATS)
    - rstudio_cytometry_v4.6.0 si the continer used for compositional analysis from cytometry data (sccomp)
    
   
01_Demultiplexing: Contains
- 01_Config_files: the configuration files to give to CellRanger for demultiplexing with probe barcodes used for each sample
- 02_Script: the slurm script to demultiplexe samples on a distant server

02_scRNASeq_Individual_sample_analysis contains all files for individual analysis of each sample (one file per sample).

03_scRNASeq_Merged_sample_analysis contains all script for analyzing merged samples alltogether:
01_All_cells contains scripts that have been used for global analyses of all cells (and especially differential composition analysis). 
02_Monocytes contains scripts that have been used for  analyses of monocytes from all individual (and especially differential composition analysis). 

04_Cytometry contains script for differential composition analysis



