# ##################################################
# Global declarations and libraries for the analysis
# ##################################################


######## R Libraries (BiocManager version 3.16, R version 4.2.2)

library( digest)        # digest (hashing)
library( DT)            # datatable
library( forcats)       # fct_inorder (factors)
library( fs)            # path_sanitize
library( future)        # plan (multicore)
library( ggplot2)
library( ggpubr)
library( ggrepel)
library( gridExtra)
library( htmltools)     # browsable
library( htmlwidgets)   # JS (for datatable)
library( iheatmapr)     # iheatmap
library( kableExtra)
library( knitr)
library( pander)        # pander
library( patchwork)     # +/ (ggplot layout)
library( pheatmap)      # pheatmap
library( plotly)
library( dplyr)
library( rmarkdown)
library( scales)        # hue_pal
library( RColorBrewer)
library(KernSmooth)   #needed for summarizeSweep()
library(ROCR)       #needed for summarizeSweep()
library(leiden)     #to load the leidenalg() package

# Single-cell technology
library( Seurat)
library( Matrix) # Force reload Matrix 1.5.3 instead of 1.5.1 that cause a bug in FindNeighbors function

# Functional Enrichment analysis
library( biomaRt)
#library( babelgene)
library( clusterProfiler)
library( org.Hs.eg.db)
#library( rrvgo) # no package of this name and doesn't seem necessary for this analysis
library( rstatix)
library(DoubletFinder)


######## DEFINITION OF FUNCTIONS #########################

# ########################################################
# Convert a list of Human genes to their homolog in Mouse
# ########################################################

convertHumanGeneList <- function(x){
  
  genesV2 = babelgene::orthologs( genes = x, human = TRUE, species = "mouse")
  
  # human = useEnsembl( biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  # mouse = useEnsembl( biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  # genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  
  return( genesV2)
}

######## Violin/Jitter plot using plotly

# Function that creates a violin plot with jitter points
# To be used for showing #UMIs, #Genes, and %Mito
plotViolinJitter = function(plotData, xAxisFormula, yAxisFormula, colorFormula = NULL, pointSize = 3, pointsColorFormula = "#444444", pointsOpacity = 0.5, xAxisTitle = "", yAxisTitle = "", hoverText = NULL, traceName = "", xTicklabels = c(""), panelWidth = 300, panelHeight = 600, thresholdHigh = NULL, thresholdLow = NULL)
{
  # Generate a jitter that will not be recomputed at each redraw
  constantJitter = jitter(rep(0, nrow(plotData)), amount = 0.15);
  # Add it to data so it can be accessed with updated formula (markers)
  plotData = cbind(plotData, jitterCoords = constantJitter + 0.17);
  
  p = plot_ly( plotData,
               x = xAxisFormula,
               y = yAxisFormula,
               name = traceName,
               fillcolor = colorFormula, # 'color' argument does not work with custom values (?)
               height = panelHeight,
               width = panelWidth) %>%
    add_trace( type='violin',
               showlegend=FALSE,
               points = FALSE,
               line = list( width = 1.5, 
                            color = "#44444444"), # impossible to specify individual colors (?)
               meanline = list(visible = TRUE),
               scalemode = 'width',
               width = 0.6,
               side = 'negative',
               hoveron = 'kde',
               hoverinfo = 'skip') %>%
    add_trace( type = "scatter",
               mode = "markers",
               x = update(xAxisFormula, ~ . + jitterCoords), # Update formula adding jitter values
               marker = list(size = pointSize,
                             color = pointsColorFormula,
                             opacity = pointsOpacity),
               text = hoverText,
               hoverinfo = 'text+name') %>%
    layout( yaxis = list( title = yAxisTitle,
                          rangemode = "nonnegative"),
            xaxis = list( title = xAxisTitle,
                          ticktext = xTicklabels,
                          tickvals = ((1:length( xTicklabels)) + 0.01)),
            margin = list( l = 0, r = 0, t = 0, b = 0, pad = 5));
  
  
  # Function creating a horizontal line as a shape to be added to layout (use 'paper' coords ref for x axis)
  hline = function(y, color = "blue") 
  {
    if (is.null( y)) return( NULL);
    list( list( type = "line", 
                x0 = 0, 
                x1 = 1, 
                xref = "paper",
                y0 = y, 
                y1 = y, 
                line = list( color = color)))
  }
  
  # Shapes to be added to layout (for eventual threshold lines)
  shapes = c( hline( thresholdHigh, "red"), hline( thresholdLow, "blue"));
  
  if(length( shapes)) return(p %>% layout( shapes = shapes));
  return(p);
}




######## Highlight cells from a specific cluster on a dimreduc plot
highlightClusterPlot = function(clusterName, seuratObject, reduction = "umap")
{
  clusterCells = which( Idents( seuratObject) == clusterName);
  # Create a dimreduc plot with current cluster highlighted (useful for large number of clusters)
  print( DimPlot( seuratObject, 
                  reduction = reduction, 
                  cols="#44444422", 
                  cells.highlight = clusterCells, 
                  cols.highlight = "#FF000088", 
                  sizes.highlight = Seurat:::AutoPointSize(seuratObject)*1.5, 
                  order = clusterCells,  # Plot highlighted cells last
                  group.by=NULL) + 
           ggtitle( paste( "Cluster", clusterName)) +
           theme( axis.title.x = element_blank(),
                  axis.title.y = element_blank(),
                  plot.title = element_text( face = "bold",
                                             size = rel( 16/14), 
                                             hjust = 0.5,
                                             vjust = 1, 
                                             margin = margin( b = 7)),
                  legend.position = "none"));

}




######## Create a violinplot showing feature (gene or metadata) values for each cluster (+ number of 'zero' and 'not zero' cells)
violinFeatureByCluster = function(currentFeature, seuratObject, slot = "counts", clustersColor = NULL, yLabel = "Counts", addStats = TRUE, addTitle = TRUE, trimTitle = FALSE) 
{
  if(length( currentFeature)>1) warning( "Several features in 'currentFeature' argument, corresponding values will be averaged...");

  # Extract feature values from Seurat object (average if several features given)
  featureValues = apply( FetchData( object = seuratObject, vars = currentFeature, slot = slot), 1, mean);

  # Plot expression values of current annotation for each cluster as violin + jitter
  ggFigure = ggplot( data.frame( Counts = featureValues, Cluster = Idents( seuratObject)), 
                                 aes( x = Cluster, y = Counts)) +
               geom_jitter( width = 0.2, height = 0.3, size = 0.4, colour = "#44444444") +
               geom_violin( aes( col = Cluster, fill = Cluster), scale = "width", alpha = 0.4, draw_quantiles = 0.5) +
               theme_minimal() +
               ylab( label = yLabel) +
               theme( legend.position = "None");

  # Change ggplot default clusters colors ('clustersColor' should be a named vector of colors)
  if(!is.null(clustersColor))
  {
    ggFigure = ggFigure + 
      scale_color_manual( values = clustersColor) +
      scale_fill_manual( values = clustersColor);
  }

  if(addStats)
  {
    # Get the number of cells in which score is zero (or not) for this annotation
    statsZero = tapply( featureValues, 
                        Idents( seuratObject), 
                        function(x) { paste( table( factor( as.logical( x), levels = c( TRUE, FALSE))), collapse="\n") });

    # Create the text layer with stats to be added
    ggFigure = ggFigure + geom_text( data = data.frame( Stat = statsZero, Cluster = names( statsZero)), 
                                     aes( x = Cluster, y = max( featureValues)*1.1, label = Stat), 
                                     size = 2.5) +
                          geom_text( data = data.frame( text = "Pos.\nNull"), 
                                     aes( label = text), 
                                     x = -0.4, y = max( featureValues)*1.1, hjust = 0, size = 2.5) +
                          coord_cartesian( clip = 'off', xlim = c( 0, length(statsZero))); # Allow drawing outside plot region
  } 

  if(addTitle)
  {
    # Remove 'trimTitle' number of ending character(s) (module names have a number suffix added automatically)
    title = if(trimTitle) substr( currentFeature, 1, nchar( currentFeature)-trimTitle) else currentFeature;

    # Create the title layer to be added (copy 'FeaturePlot' title style from cowplot)
    ggFigure = ggFigure + ggtitle( label = title) + 
                          theme( plot.title =  element_text( face = "bold",
                                                             size = rel( 16/14), 
                                                             hjust = 0.5,
                                                             vjust = 1, 
                                                             margin = margin( b = 7)));
  }

  # Render figure
  #print( ggFigure)
  suppressWarnings(print( ggFigure));
}




######## Datatable default buttons
# Configure a set of buttons for exporting data from datatables
# (NOTE: an orthogonal named 'export' must be defined in the concerned columns for proper formatting)
exportButtonsListDT = list( 
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save full table
    extend = "collection",
    text = "<center>Export<br>Full</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index', # give result in the same order as original table
                                                                page     = 'all',   # do not restrict to currently displayed page
                                                                search   = 'none'), # do not restrict to current filtered rows
                                                # do not define 'selected' here as 'undefined' gives expected behavior (selection ignored, not sure how to set to null or undefined using R)
                                                orthogonal = 'export')), # use 'export' orthogonal source (no formatting, defined in 'columnsDef')
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none'),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none'),
                                                orthogonal = 'export')))),
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save filtered table
    extend = "collection",
    text = "<center>Export<br>Filtered</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index',      # give result in the same order as original table
                                                                page     = 'all',        # do not restrict to currently displayed page
                                                                search   = 'applied'),   # restrict to current filtered rows
                                                # do not define 'selected' here as 'undefined' gives expected behavior (selection ignored, not sure how to set to null or undefined using R)
                                                orthogonal = 'export')),
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'applied'),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'applied'),
                                                orthogonal = 'export')))),
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save rows from current page only
    extend = "collection",
    text = "<center>Export<br>Page</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index',      # give result in the same order as original table
                                                                page     = 'current',    # restrict to currently displayed page
                                                                search   = 'none'),      # do not restrict to current filtered rows
                                                # do not define 'selected' here as 'undefined' gives expected behavior (selection ignored, not sure how to set to null or undefined using R)
                                                orthogonal = 'export')),
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'current',
                                                                search   = 'none'),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'current',
                                                                search   = 'none'),
                                                orthogonal = 'export')))),
  list( # Collection containing export buttons with altered 'exportOptions' to use non-formatted/non-ordered data and save selected rows only (independently of the filtering)
    extend = "collection",
    text = "<center>Export<br>Selected</center>",
    buttons = list( list( extend  = "excel",
                          exportOptions = list( modifier = list(order    = 'index',      # give result in the same order as original table
                                                                page     = 'all',        # do not restrict to currently displayed page
                                                                search   = 'none',       # do not restrict to current filtered rows
                                                                selected = TRUE),        # extension Select
                                                orthogonal = 'export')),
                    list( extend  = "csv",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none',
                                                                selected = TRUE),
                                                orthogonal = 'export')),
                    list( extend  = "copy",
                          exportOptions = list( modifier = list(order    = 'index',
                                                                page     = 'all',
                                                                search   = 'none',
                                                                selected = TRUE),
                                                orthogonal = 'export')))),
  # # -- Print buttons
  list( # Collection containing print buttons with altered 'exportOptions' to render visible columns only in displayed rows order, and select desired rows (see comments in export buttons, here it uses default orthogonal that gives formatted data)
    extend = "collection",
    text = "<center>Print<br>(as seen)</center>",
    buttons = list( list( extend  = "print",
                          text = "Full",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'all',
                                                                search   = 'none'),
                                                columns = ':visible')),
                    list( extend  = "print",
                          text = "Filtered",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'all',
                                                                search   = 'applied'),
                                                columns = ':visible')),
                    list( extend  = "print",
                          text = "Page",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'current',
                                                                search   = 'none'),
                                                columns = ':visible')),
                    list( extend  = "print",
                          text = "Selected",
                          exportOptions = list( modifier = list(order    = 'current',
                                                                page     = 'all',
                                                                search   = 'none',
                                                                selected = TRUE),
                                                columns = ':visible')))),
  list( extend = "selectNone",
        text = "<center>Deselect<br>all</center>"));




######## Generate a datatable summarizing values 
# For environments (parameters), all values/variables are shown

showSimpleDT = function( dataToShow, rownames = TRUE, colnames = "Value")
{
  valuesDF = NULL;
  
  if(is.environment( dataToShow))
  {
    # Extract values from environment as character strings
    envValues = sapply(lapply(dataToShow, function(x) {ifelse(is.null(x), "NULL", x)}), paste, collapse = ", ");
    # Sort them by name and convert to data.frame
    valuesDF = data.frame("Value" = envValues[order(names(envValues))]);
  } else
  {
    valuesDF = dataToShow;
  }
  
  # Create a datatable with collected information
  # Create datatable
  datatable( as.data.frame(valuesDF), 
             class = "compact",
             rownames = rownames,
             colnames = colnames,
             options = list(dom = "<'row'rt>", # Set elements for CSS formatting ('<Blf><rt><ip>')
                            autoWidth = FALSE,
                            columnDefs = list( # Center all columns
                              list( targets = "_all",
                                    className = 'dt-center')),
                            orderClasses = FALSE, # Disable flag for CSS to highlight columns used for ordering (for performance)
                            ordering = FALSE,
                            paging = FALSE, # Disable pagination (show all)
                            processing = TRUE, 
                            scrollCollapse = TRUE,
                            scroller = TRUE,  # Only load visible data
                            scrollX = TRUE,
                            scrollY = "525px",
                            stateSave = TRUE));
}


######## 

#Fuction of DoubletFinder for Seurat4

parallel_paramSweep_v2 <- function(n, n.real.cells, real.cells, pK, pN, data, orig.commands, PCs, sct)  {
  
  sweep.res.list = list()
  list.ind = 0
  
  ## Make merged real-artifical data
  # print(paste("Creating artificial doublets for pN = ", pN[n]*100,"%",sep=""))
  n_doublets <- round(n.real.cells/(1 - pN[n]) - n.real.cells)
  real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
  real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
  doublets <- (data[, real.cells1] + data[, real.cells2])/2
  colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
  data_wdoublets <- cbind(data, doublets)
  
  ## Pre-process Seurat object
  if (sct == FALSE) {
    # print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    
    # print("Normalizing Seurat object...")
    seu_wdoublets <- NormalizeData(seu_wdoublets,
                                   normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                   scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                   margin = orig.commands$NormalizeData.RNA@params$margin)
    
    # print("Finding variable genes...")
    seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                          selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                          loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                          clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                          mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                          dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                          num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                          binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                          nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                          mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                          dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
    
    # print("Scaling data...")
    seu_wdoublets <- ScaleData(seu_wdoublets,
                               features = orig.commands$ScaleData.RNA$features,
                               model.use = orig.commands$ScaleData.RNA$model.use,
                               do.scale = orig.commands$ScaleData.RNA$do.scale,
                               do.center = orig.commands$ScaleData.RNA$do.center,
                               scale.max = orig.commands$ScaleData.RNA$scale.max,
                               block.size = orig.commands$ScaleData.RNA$block.size,
                               min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
    
    # print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets,
                            features = orig.commands$ScaleData.RNA$features,
                            npcs = length(PCs),
                            rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                            weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                            verbose=FALSE)
  }
  
  if (sct == TRUE) {
    require(sctransform)
    # print("Creating Seurat object...")
    seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
    
    # print("Running SCTransform...")
    seu_wdoublets <- SCTransform(seu_wdoublets)
    
    # print("Running PCA...")
    seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
  }
  
  ## Compute PC distance matrix
  # print("Calculating PC distance matrix...")
  nCells <- nrow(seu_wdoublets@meta.data)
  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
  rm(seu_wdoublets)
  gc()
  dist.mat <- fields::rdist(pca.coord)[,1:n.real.cells]
  
  ## Pre-order PC distance matrix prior to iterating across pK for pANN computations
  # print("Defining neighborhoods...")
  for (i in 1:n.real.cells) {
    dist.mat[,i] <- order(dist.mat[,i])
  }
  
  ## Trim PC distance matrix for faster manipulations
  ind <- round(nCells * max(pK))+5
  dist.mat <- dist.mat[1:ind, ]
  
  ## Compute pANN across pK sweep
  # print("Computing pANN across all pK...")
  for (k in 1:length(pK)) {
    # print(paste("pK = ", pK[k], "...", sep = ""))
    pk.temp <- round(nCells * pK[k])
    pANN <- as.data.frame(matrix(0L, nrow = n.real.cells, ncol = 1))
    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells
    list.ind <- list.ind + 1
    
    for (i in 1:n.real.cells) {
      neighbors <- dist.mat[2:(pk.temp + 1),i]
      pANN$pANN[i] <- length(which(neighbors > n.real.cells))/pk.temp
    }
    
    sweep.res.list[[list.ind]] <- pANN
    
  }
  
  return(sweep.res.list)
}

paramSweep_v2 <- function(seu, PCs=1:10, sct = FALSE, num.cores=1) {
  require(Seurat); require(fields); require(parallel)
  ## Set pN-pK param sweep ranges
  pK <- c(0.0005, 0.001, 0.005, seq(0.01,0.3,by=0.01))
  pN <- seq(0.05,0.3,by=0.05)
  
  ## Remove pK values with too few cells
  min.cells <- round(nrow(seu@meta.data)/(1-0.05) - nrow(seu@meta.data))
  pK.test <- round(pK*min.cells)
  pK <- pK[which(pK.test >= 1)]
  
  ## Extract pre-processing parameters from original data analysis workflow
  orig.commands <- seu@commands
  
  ## Down-sample cells to 10000 (when applicable) for computational effiency
  if (nrow(seu@meta.data) > 10000) {
    real.cells <- rownames(seu@meta.data)[sample(1:nrow(seu@meta.data), 10000, replace=FALSE)]
    data <-seu@assays$RNA@counts[ , real.cells]
    n.real.cells <- ncol(data)
  }
  
  if (nrow(seu@meta.data) <= 10000){
    real.cells <- rownames(seu@meta.data)
    data <-seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }
  
  ## Iterate through pN, computing pANN vectors at varying pK
  #no_cores <- detectCores()-1
  if(num.cores>1){
    require(parallel)
    cl <- makeCluster(num.cores)
    output2 <- mclapply(as.list(1:length(pN)),
                        FUN = parallel_paramSweep_v2,
                        n.real.cells,
                        real.cells,
                        pK,
                        pN,
                        data,
                        orig.commands,
                        PCs,
                        sct,mc.cores=num.cores)
    stopCluster(cl)
  }else{
    output2 <- lapply(as.list(1:length(pN)),
                      FUN = parallel_paramSweep_v2,
                      n.real.cells,
                      real.cells,
                      pK,
                      pN,
                      data,
                      orig.commands,
                      PCs,
                      sct)
  }
  
  ## Write parallelized output into list
  sweep.res.list <- list()
  list.ind <- 0
  for(i in 1:length(output2)){
    for(j in 1:length(output2[[i]])){
      list.ind <- list.ind + 1
      sweep.res.list[[list.ind]] <- output2[[i]][[j]]
    }
  }
  
  ## Assign names to list of results
  name.vec <- NULL
  for (j in 1:length(pN)) {
    name.vec <- c(name.vec, paste("pN", pN[j], "pK", pK, sep = "_" ))
  }
  names(sweep.res.list) <- name.vec
  return(sweep.res.list)
  
}

doubletFinder_v2 <- function(seu, PCs, pN = 0.25, pK, nExp, reuse.pANN = FALSE, sct = FALSE, annotations = NULL) {
  require(Seurat); require(fields); require(KernSmooth)
  
  ## Generate new list of doublet classificatons from existing pANN vector to save time
  if (reuse.pANN != FALSE ) {
    pANN.old <- seu@meta.data[ , reuse.pANN]
    classifications <- rep("Singlet", length(pANN.old))
    classifications[order(pANN.old, decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    return(seu)
  }
  
  if (reuse.pANN == FALSE) {
    ## Make merged real-artifical data
    real.cells <- rownames(seu@meta.data)
    data <- seu@assays$RNA@counts[, real.cells]
    n_real.cells <- length(real.cells)
    n_doublets <- round(n_real.cells/(1 - pN) - n_real.cells)
    # print(paste("Creating",n_doublets,"artificial doublets...",sep=" "))
    real.cells1 <- sample(real.cells, n_doublets, replace = TRUE)
    real.cells2 <- sample(real.cells, n_doublets, replace = TRUE)
    doublets <- (data[, real.cells1] + data[, real.cells2])/2
    colnames(doublets) <- paste("X", 1:n_doublets, sep = "")
    data_wdoublets <- cbind(data, doublets)
    # Keep track of the types of the simulated doublets
    if(!is.null(annotations)){
      stopifnot(typeof(annotations)=="character")
      stopifnot(length(annotations)==length(Cells(seu)))
      stopifnot(!any(is.na(annotations)))
      annotations <- factor(annotations)
      names(annotations) <- Cells(seu)
      doublet_types1 <- annotations[real.cells1]
      doublet_types2 <- annotations[real.cells2]
    }
    ## Store important pre-processing information
    orig.commands <- seu@commands
    
    ## Pre-process Seurat object
    if (sct == FALSE) {
      # print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      # print("Normalizing Seurat object...")
      seu_wdoublets <- NormalizeData(seu_wdoublets,
                                     normalization.method = orig.commands$NormalizeData.RNA@params$normalization.method,
                                     scale.factor = orig.commands$NormalizeData.RNA@params$scale.factor,
                                     margin = orig.commands$NormalizeData.RNA@params$margin)
      
      # print("Finding variable genes...")
      seu_wdoublets <- FindVariableFeatures(seu_wdoublets,
                                            selection.method = orig.commands$FindVariableFeatures.RNA$selection.method,
                                            loess.span = orig.commands$FindVariableFeatures.RNA$loess.span,
                                            clip.max = orig.commands$FindVariableFeatures.RNA$clip.max,
                                            mean.function = orig.commands$FindVariableFeatures.RNA$mean.function,
                                            dispersion.function = orig.commands$FindVariableFeatures.RNA$dispersion.function,
                                            num.bin = orig.commands$FindVariableFeatures.RNA$num.bin,
                                            binning.method = orig.commands$FindVariableFeatures.RNA$binning.method,
                                            nfeatures = orig.commands$FindVariableFeatures.RNA$nfeatures,
                                            mean.cutoff = orig.commands$FindVariableFeatures.RNA$mean.cutoff,
                                            dispersion.cutoff = orig.commands$FindVariableFeatures.RNA$dispersion.cutoff)
      
      # print("Scaling data...")
      seu_wdoublets <- ScaleData(seu_wdoublets,
                                 features = orig.commands$ScaleData.RNA$features,
                                 model.use = orig.commands$ScaleData.RNA$model.use,
                                 do.scale = orig.commands$ScaleData.RNA$do.scale,
                                 do.center = orig.commands$ScaleData.RNA$do.center,
                                 scale.max = orig.commands$ScaleData.RNA$scale.max,
                                 block.size = orig.commands$ScaleData.RNA$block.size,
                                 min.cells.to.block = orig.commands$ScaleData.RNA$min.cells.to.block)
      
      # print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets,
                              features = orig.commands$ScaleData.RNA$features,
                              npcs = length(PCs),
                              rev.pca =  orig.commands$RunPCA.RNA$rev.pca,
                              weight.by.var = orig.commands$RunPCA.RNA$weight.by.var,
                              verbose=FALSE)
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc() # Free up memory
    }
    
    if (sct == TRUE) {
      require(sctransform)
      # print("Creating Seurat object...")
      seu_wdoublets <- CreateSeuratObject(counts = data_wdoublets)
      
      # print("Running SCTransform...")
      seu_wdoublets <- SCTransform(seu_wdoublets)
      
      # print("Running PCA...")
      seu_wdoublets <- RunPCA(seu_wdoublets, npcs = length(PCs))
      pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[ , PCs]
      cell.names <- rownames(seu_wdoublets@meta.data)
      nCells <- length(cell.names)
      rm(seu_wdoublets); gc()
    }
    
    ## Compute PC distance matrix
    # print("Calculating PC distance matrix...")
    dist.mat <- fields::rdist(pca.coord)
    
    ## Compute pANN
    # print("Computing pANN...")
    pANN <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = 1))
    if(!is.null(annotations)){
      neighbor_types <- as.data.frame(matrix(0L, nrow = n_real.cells, ncol = length(levels(doublet_types1))))
    }
    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"
    k <- round(nCells * pK)
    for (i in 1:n_real.cells) {
      neighbors <- order(dist.mat[, i])
      neighbors <- neighbors[2:(k + 1)]
      pANN$pANN[i] <- length(which(neighbors > n_real.cells))/k
      if(!is.null(annotations)){
        for(ct in unique(annotations)){
          neighbors_that_are_doublets = neighbors[neighbors>n_real.cells]
          if(length(neighbors_that_are_doublets) > 0){
            neighbor_types[i,] <-
              table( doublet_types1[neighbors_that_are_doublets - n_real.cells] ) +
              table( doublet_types2[neighbors_that_are_doublets - n_real.cells] )
            neighbor_types[i,] <- neighbor_types[i,] / sum( neighbor_types[i,] )
          } else {
            neighbor_types[i,] <- NA
          }
        }
      }
    }
    # print("Classifying doublets..")
    classifications <- rep("Singlet",n_real.cells)
    classifications[order(pANN$pANN[1:n_real.cells], decreasing=TRUE)[1:nExp]] <- "Doublet"
    seu@meta.data[, paste("pANN",pN,pK,nExp,sep="_")] <- pANN[rownames(seu@meta.data), 1]
    seu@meta.data[, paste("DF.classifications",pN,pK,nExp,sep="_")] <- classifications
    if(!is.null(annotations)){
      colnames(neighbor_types) = levels(doublet_types1)
      for(ct in levels(doublet_types1)){
        seu@meta.data[, paste("DF.doublet.contributors",pN,pK,nExp,ct,sep="_")] <- neighbor_types[,ct]
      }
    }
    return(seu)
  }
}



###################################### SciGeneX ###########################################################

library(Seurat)
library(scigenex)
library(dplyr)
library(dbplyr)
library(patchwork)
library(ggplot2)
library(viridis)
library(scCustomize)
library(DT)
library( scales) 

#install.packages("dbplyr")
#devtools::install_github("RobinHankin/Brobdingnag")
set_verbosity(1)


#################################### Monocle3 (Trajectory analysis) ###########################################

library(monocle3)
library(R.utils)
library(SeuratWrappers)






