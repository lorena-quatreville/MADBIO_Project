# #########################################
# This script compute and filters multiplet
# #########################################


# ..........................................................................................................
## @knitr heterogeneity_identificationdoublet_DoubletFinder
# ..........................................................................................................

options(future.globals.maxSize = 8000 * 1024^2)
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

# Rename the name of metadata create by Doublet Finder in multiplet_prediction
colnames(sc10x@meta.data)[which(colnames(sc10x@meta.data) == 
                         paste("DF.classifications", pN, optimal.pk, nExp_poi.adj, sep = "_"))] <- "multiplet_prediction"


cat("<br>Doublet formation rate:", DOUBLET_FORMATION_RATE );
cat("<br>Limit of cells considered as doublet:", nExp_poi.adj );
cat("<br>Optimal pK:", optimal.pk );
cat("\n<br>\n");


# ..........................................................................................................
## @knitr filterDoublet_filterObject
# ..........................................................................................................

sc10x_Filtered_withMultiplet = sc10x;
# Filter the excluded cells in the Seurat object
Idents(object = sc10x) = "multiplet_prediction"
if( DELETE_MULTIPET == TRUE){
  sc10x <- subset(sc10x, idents = "Singlet");
}


nbPC=PCA_NPC
if(PCA_NPC>length(Cells(sc10x))){
  warning( paste0( "Number of cells in object (", length(Cells(sc10x)), ") smaller than requested number of PCs (", PCA_NPC,"), setting lower PC number..." ))
  nbPC = length(Cells(sc10x))
}           

sc10x <- RunPCA( object   = sc10x,
                 features = VariableFeatures( sc10x),
                 npcs     = nbPC,
                 verbose  = .VERBOSE)



