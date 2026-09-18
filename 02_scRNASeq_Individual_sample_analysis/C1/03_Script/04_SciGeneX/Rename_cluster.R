
################################################################################
#                               Rename cluster                                 #
################################################################################

# Create the seurat object (RNA only)
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x)<-"seurat_clusters"

sc10x[["old.seurat_clusters"]] <- Idents(sc10x)

sc10x <- RenameIdents(object = sc10x,
                      "1" = "Classical monocytes", "2" = "Naive and memory CD4 T cells", "3" = "Naive and memory CD8 T cells", "4" = "?", "5" = "Effector CD4 T cells", "6" = "Class-switched memory B cells", "7" = "Neutrophils", "8" = "Non-classical monocytes", "9" = "Naive B cells", "10" = "Plasmacitoid dendritic cells", "11" = "Plasma cells", "12" = "?", "13" = "Conventional dendritic cells", "14" = "Pletelets", "15" = "Proliferating cells", "16"= "Proliferating cells", "17" = "Plasmacitoid dendritic cells", "18" = "?")
                      
# atB cell : FCRL5, ITGAX, TNFRSF1B, LILRB1, CD19, and MS4A1
# proliferation : MKI67 TOP2A  
# pre-DC : AXL SIGLEC6
# cDC1 CLEC9A CADM1 XCR1 



sc10x[["seurat_clusters"]] <- Idents(sc10x)



###############################################################################
