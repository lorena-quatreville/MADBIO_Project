
# ##############################################################
            # Doublet Detection Script: DoubletFinder    #
# ##############################################################


# Load the seurat object


sc10x <- readRDS(PATH_SEURAT_RDS_FILE)



# pK Identification (no ground-truth) 
sweep.list <- paramSweep_v2(sc10x, PCs = 1:DOUBLETFINDER_USE_PCA_NBDIMS, sct = FALSE);
sweep.stats <- summarizeSweep(sweep.list, GT = FALSE);
bcmvn<- find.pK(sweep.stats);

bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),];
optimal.pk <- bcmvn.max$pK;
optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk];

# Homotypic Doublet Proportion Estimate 
annotations <- sc10x@meta.data$seurat_clusters;
homotypic.prop <- modelHomotypic(annotations); 
nExp_poi <- round(DOUBLET_FORMATION_RATE*nrow(sc10x@meta.data));
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop));

# Run DoubletFinder with varying classification stringencies
sc10x <- doubletFinder_v2(seu = sc10x, 
                          PCs = 1:DOUBLETFINDER_USE_PCA_NBDIMS, 
                          pN = pN, 
                          pK = optimal.pk, 
                          nExp = nExp_poi.adj, 
                          reuse.pANN = FALSE, 
                          sct = FALSE);
                          
                          
# Rename the name of metadata create by Doublet Finder in DF_prediction                          
colnames(sc10x@meta.data)[which(colnames(sc10x@meta.data) == 
                         paste("DF.classifications", pN, optimal.pk, nExp_poi.adj, sep = "_"))] <- "DF_prediction"
                         
                          
cat("<br>Doublet formation rate:", DOUBLET_FORMATION_RATE)
cat("<br>Expected doublets:", nExp_poi)
cat("<br>Adjusted doublets:", nExp_poi.adj)
cat("<br>Optimal pK:", optimal.pk, "\n")
              

