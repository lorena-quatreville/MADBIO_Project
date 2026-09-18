# ##################################################
# Global declarations and libraries for the analysis
# ##################################################

######## R Libraries

library(digest)
library(DT)
library(forcats)
library(fs)
library(future)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(htmltools)
library(htmlwidgets)
library(kableExtra)
library(knitr)
library(pander)
library(patchwork)
library(plotly)
library(dplyr)
library(tidyr)
library(rmarkdown)
library(scales)
library(RColorBrewer)
library(KernSmooth)
#library(ROCR)
#library(ggbeeswarm)
library(stringr)
library(ggsignif)

######## sccomp / CmdStan

library(cmdstanr)
library(sccomp)


# ##################################################
# DEFINITION OF FUNCTIONS
# ##################################################


# ########################################################
# Convert a list of Human genes to their homolog in Mouse
# ########################################################

convertHumanGeneList <- function(x) {

  genesV2 <- babelgene::orthologs(
    genes = x,
    human = TRUE,
    species = "mouse"
  )

  return(genesV2)
}


# ########################################################
# Violin/Jitter plot using plotly
# ########################################################

plotViolinJitter <- function(
    plotData,
    xAxisFormula,
    yAxisFormula,
    colorFormula = NULL,
    pointSize = 3,
    pointsColorFormula = "#444444",
    pointsOpacity = 0.5,
    xAxisTitle = "",
    yAxisTitle = "",
    hoverText = NULL,
    traceName = "",
    xTicklabels = c(""),
    panelWidth = 300,
    panelHeight = 600,
    thresholdHigh = NULL,
    thresholdLow = NULL
) {

  constantJitter <- jitter(
    rep(0, nrow(plotData)),
    amount = 0.15
  )

  plotData <- cbind(
    plotData,
    jitterCoords = constantJitter + 0.17
  )

  p <- plot_ly(
    plotData,
    x = xAxisFormula,
    y = yAxisFormula,
    name = traceName,
    fillcolor = colorFormula,
    height = panelHeight,
    width = panelWidth
  ) %>%
    add_trace(
      type = "violin",
      showlegend = FALSE,
      points = FALSE,
      line = list(
        width = 1.5,
        color = "#44444444"
      ),
      meanline = list(visible = TRUE),
      scalemode = "width",
      width = 0.6,
      side = "negative",
      hoveron = "kde",
      hoverinfo = "skip"
    ) %>%
    add_trace(
      type = "scatter",
      mode = "markers",
      x = update(
        xAxisFormula,
        ~ . + jitterCoords
      ),
      marker = list(
        size = pointSize,
        color = pointsColorFormula,
        opacity = pointsOpacity
      ),
      text = hoverText,
      hoverinfo = "text+name"
    ) %>%
    layout(
      yaxis = list(
        title = yAxisTitle,
        rangemode = "nonnegative"
      ),
      xaxis = list(
        title = xAxisTitle,
        ticktext = xTicklabels,
        tickvals = ((1:length(xTicklabels)) + 0.01)
      ),
      margin = list(
        l = 0,
        r = 0,
        t = 0,
        b = 0,
        pad = 5
      )
    )

  hline <- function(y, color = "blue") {

    if (is.null(y))
      return(NULL)

    list(
      list(
        type = "line",
        x0 = 0,
        x1 = 1,
        xref = "paper",
        y0 = y,
        y1 = y,
        line = list(color = color)
      )
    )
  }

  shapes <- c(
    hline(thresholdHigh, "red"),
    hline(thresholdLow, "blue")
  )

  if (length(shapes))
    return(p %>% layout(shapes = shapes))

  return(p)
}


# ########################################################
# Highlight cells from a specific cluster
# ########################################################

highlightClusterPlot <- function(
    clusterName,
    seuratObject,
    reduction = "umap"
) {

  clusterCells <- which(
    Idents(seuratObject) == clusterName
  )

  print(
    DimPlot(
      seuratObject,
      reduction = reduction,
      cols = "#44444422",
      cells.highlight = clusterCells,
      cols.highlight = "#FF000088",
      sizes.highlight = Seurat:::AutoPointSize(seuratObject) * 1.5,
      order = clusterCells,
      group.by = NULL
    ) +
      ggtitle(paste("Cluster", clusterName)) +
      theme(
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        plot.title = element_text(
          face = "bold",
          size = rel(16 / 14),
          hjust = 0.5,
          vjust = 1,
          margin = margin(b = 7)
        ),
        legend.position = "none"
      )
  )
}


# ########################################################
# Violin plot by cluster
# ########################################################

violinFeatureByCluster <- function(
    currentFeature,
    seuratObject,
    slot = "counts",
    clustersColor = NULL,
    yLabel = "Counts",
    addStats = TRUE,
    addTitle = TRUE,
    trimTitle = FALSE
) {

  if (length(currentFeature) > 1)
    warning(
      "Several features in 'currentFeature' argument, corresponding values will be averaged..."
    )

  featureValues <- apply(
    FetchData(
      object = seuratObject,
      vars = currentFeature,
      slot = slot
    ),
    1,
    mean
  )

  ggFigure <- ggplot(
    data.frame(
      Counts = featureValues,
      Cluster = Idents(seuratObject)
    ),
    aes(
      x = Cluster,
      y = Counts
    )
  ) +
    geom_jitter(
      width = 0.2,
      height = 0.3,
      size = 0.4,
      colour = "#44444444"
    ) +
    geom_violin(
      aes(
        col = Cluster,
        fill = Cluster
      ),
      scale = "width",
      alpha = 0.4,
      draw_quantiles = 0.5
    ) +
    theme_minimal() +
    ylab(label = yLabel) +
    theme(legend.position = "none")

  if (!is.null(clustersColor)) {

    ggFigure <- ggFigure +
      scale_color_manual(values = clustersColor) +
      scale_fill_manual(values = clustersColor)
  }

  if (addStats) {

    statsZero <- tapply(
      featureValues,
      Idents(seuratObject),
      function(x) {
        paste(
          table(
            factor(
              as.logical(x),
              levels = c(TRUE, FALSE)
            )
          ),
          collapse = "\n"
        )
      }
    )

    ggFigure <- ggFigure +
      geom_text(
        data = data.frame(
          Stat = statsZero,
          Cluster = names(statsZero)
        ),
        aes(
          x = Cluster,
          y = max(featureValues) * 1.1,
          label = Stat
        ),
        size = 2.5
      ) +
      geom_text(
        data = data.frame(
          text = "Pos.\nNull"
        ),
        aes(label = text),
        x = -0.4,
        y = max(featureValues) * 1.1,
        hjust = 0,
        size = 2.5
      ) +
      coord_cartesian(
        clip = "off",
        xlim = c(0, length(statsZero))
      )
  }

  if (addTitle) {

    title <- if (trimTitle)
      substr(
        currentFeature,
        1,
        nchar(currentFeature) - trimTitle
      )
    else
      currentFeature

    ggFigure <- ggFigure +
      ggtitle(label = title) +
      theme(
        plot.title = element_text(
          face = "bold",
          size = rel(16 / 14),
          hjust = 0.5,
          vjust = 1,
          margin = margin(b = 7)
        )
      )
  }

  suppressWarnings(
    print(ggFigure)
  )
}


# ########################################################
# Datatable default buttons
# ########################################################

exportButtonsListDT <- list(

  list(
    extend = "collection",
    text = "<center>Export<br>Full</center>",
    buttons = list(
      list(
        extend = "excel",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "all",
            search = "none"
          ),
          orthogonal = "export"
        )
      ),
      list(
        extend = "csv",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "all",
            search = "none"
          ),
          orthogonal = "export"
        )
      ),
      list(
        extend = "copy",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "all",
            search = "none"
          ),
          orthogonal = "export"
        )
      )
    )
  ),

  list(
    extend = "collection",
    text = "<center>Export<br>Filtered</center>",
    buttons = list(
      list(
        extend = "excel",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "all",
            search = "applied"
          ),
          orthogonal = "export"
        )
      ),
      list(
        extend = "csv",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "all",
            search = "applied"
          ),
          orthogonal = "export"
        )
      ),
      list(
        extend = "copy",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "all",
            search = "applied"
          ),
          orthogonal = "export"
        )
      )
    )
  ),

  list(
    extend = "collection",
    text = "<center>Export<br>Page</center>",
    buttons = list(
      list(
        extend = "excel",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "current",
            search = "none"
          ),
          orthogonal = "export"
        )
      ),
      list(
        extend = "csv",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "current",
            search = "none"
          ),
          orthogonal = "export"
        )
      ),
      list(
        extend = "copy",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "current",
            search = "none"
          ),
          orthogonal = "export"
        )
      )
    )
  ),

  list(
    extend = "collection",
    text = "<center>Export<br>Selected</center>",
    buttons = list(
      list(
        extend = "excel",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "all",
            search = "none",
            selected = TRUE
          ),
          orthogonal = "export"
        )
      ),
      list(
        extend = "csv",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "all",
            search = "none",
            selected = TRUE
          ),
          orthogonal = "export"
        )
      ),
      list(
        extend = "copy",
        exportOptions = list(
          modifier = list(
            order = "index",
            page = "all",
            search = "none",
            selected = TRUE
          ),
          orthogonal = "export"
        )
      )
    )
  ),

  list(
    extend = "collection",
    text = "<center>Print<br>(as seen)</center>",
    buttons = list(
      list(
        extend = "print",
        text = "Full",
        exportOptions = list(
          modifier = list(
            order = "current",
            page = "all",
            search = "none"
          ),
          columns = ":visible"
        )
      ),
      list(
        extend = "print",
        text = "Filtered",
        exportOptions = list(
          modifier = list(
            order = "current",
            page = "all",
            search = "applied"
          ),
          columns = ":visible"
        )
      ),
      list(
        extend = "print",
        text = "Page",
        exportOptions = list(
          modifier = list(
            order = "current",
            page = "current",
            search = "none"
          ),
          columns = ":visible"
        )
      ),
      list(
        extend = "print",
        text = "Selected",
        exportOptions = list(
          modifier = list(
            order = "current",
            page = "all",
            search = "none",
            selected = TRUE
          ),
          columns = ":visible"
        )
      )
    )
  ),

  list(
    extend = "selectNone",
    text = "<center>Deselect<br>all</center>"
  )
)


# ########################################################
# Generate a datatable summarizing values
# ########################################################

showSimpleDT <- function(
    dataToShow,
    rownames = TRUE,
    colnames = "Value"
) {

  valuesDF <- NULL

  if (is.environment(dataToShow)) {

    envValues <- sapply(
      lapply(
        dataToShow,
        function(x) {
          ifelse(is.null(x), "NULL", x)
        }
      ),
      paste,
      collapse = ", "
    )

    valuesDF <- data.frame(
      "Value" = envValues[order(names(envValues))]
    )

  } else {

    valuesDF <- dataToShow
  }

  datatable(
    as.data.frame(valuesDF),
    class = "compact",
    rownames = rownames,
    colnames = colnames,
    options = list(
      dom = "<'row'rt>",
      autoWidth = FALSE,
      columnDefs = list(
        list(
          targets = "_all",
          className = "dt-center"
        )
      ),
      orderClasses = FALSE,
      ordering = FALSE,
      paging = FALSE,
      processing = TRUE,
      scrollCollapse = TRUE,
      scroller = TRUE,
      scrollX = TRUE,
      scrollY = "525px",
      stateSave = TRUE
    )
  )
}


# ########################################################
# DoubletFinder custom functions
# ########################################################

parallel_paramSweep_v2 <- function(
    n,
    n.real.cells,
    real.cells,
    pK,
    pN,
    data,
    orig.commands,
    PCs,
    sct
) {

  sweep.res.list <- list()
  list.ind <- 0

  n_doublets <- round(
    n.real.cells / (1 - pN[n]) - n.real.cells
  )

  real.cells1 <- sample(
    real.cells,
    n_doublets,
    replace = TRUE
  )

  real.cells2 <- sample(
    real.cells,
    n_doublets,
    replace = TRUE
  )

  doublets <- (
    data[, real.cells1] +
      data[, real.cells2]
  ) / 2

  colnames(doublets) <- paste0(
    "X",
    1:n_doublets
  )

  data_wdoublets <- cbind(
    data,
    doublets
  )

  if (!sct) {

    seu_wdoublets <- CreateSeuratObject(
      counts = data_wdoublets
    )

    seu_wdoublets <- NormalizeData(
      seu_wdoublets,
      normalization.method =
        orig.commands$NormalizeData.RNA@params$normalization.method,
      scale.factor =
        orig.commands$NormalizeData.RNA@params$scale.factor,
      margin =
        orig.commands$NormalizeData.RNA@params$margin
    )

    seu_wdoublets <- FindVariableFeatures(
      seu_wdoublets,
      selection.method =
        orig.commands$FindVariableFeatures.RNA$selection.method,
      loess.span =
        orig.commands$FindVariableFeatures.RNA$loess.span,
      clip.max =
        orig.commands$FindVariableFeatures.RNA$clip.max,
      mean.function =
        orig.commands$FindVariableFeatures.RNA$mean.function,
      dispersion.function =
        orig.commands$FindVariableFeatures.RNA$dispersion.function,
      num.bin =
        orig.commands$FindVariableFeatures.RNA$num.bin,
      binning.method =
        orig.commands$FindVariableFeatures.RNA$binning.method,
      nfeatures =
        orig.commands$FindVariableFeatures.RNA$nfeatures,
      mean.cutoff =
        orig.commands$FindVariableFeatures.RNA$mean.cutoff,
      dispersion.cutoff =
        orig.commands$FindVariableFeatures.RNA$dispersion.cutoff
    )

    seu_wdoublets <- ScaleData(
      seu_wdoublets,
      features =
        orig.commands$ScaleData.RNA$features,
      model.use =
        orig.commands$ScaleData.RNA$model.use,
      do.scale =
        orig.commands$ScaleData.RNA$do.scale,
      do.center =
        orig.commands$ScaleData.RNA$do.center,
      scale.max =
        orig.commands$ScaleData.RNA$scale.max,
      block.size =
        orig.commands$ScaleData.RNA$block.size,
      min.cells.to.block =
        orig.commands$ScaleData.RNA$min.cells.to.block
    )

    seu_wdoublets <- RunPCA(
      seu_wdoublets,
      features =
        orig.commands$ScaleData.RNA$features,
      npcs = length(PCs),
      rev.pca =
        orig.commands$RunPCA.RNA$rev.pca,
      weight.by.var =
        orig.commands$RunPCA.RNA$weight.by.var,
      verbose = FALSE
    )
  }

  if (sct) {

    requireNamespace(
      "sctransform"
    )

    seu_wdoublets <- CreateSeuratObject(
      counts = data_wdoublets
    )

    seu_wdoublets <- SCTransform(
      seu_wdoublets
    )

    seu_wdoublets <- RunPCA(
      seu_wdoublets,
      npcs = length(PCs)
    )
  }

  nCells <- nrow(
    seu_wdoublets@meta.data
  )

  pca.coord <- seu_wdoublets@reductions$pca@cell.embeddings[, PCs]

  rm(seu_wdoublets)
  gc()

  dist.mat <- fields::rdist(
    pca.coord
  )[, 1:n.real.cells]

  for (i in 1:n.real.cells) {
    dist.mat[, i] <- order(
      dist.mat[, i]
    )
  }

  ind <- round(
    nCells * max(pK)
  ) + 5

  dist.mat <- dist.mat[1:ind, ]

  for (k in 1:length(pK)) {

    pk.temp <- round(
      nCells * pK[k]
    )

    pANN <- as.data.frame(
      matrix(
        0L,
        nrow = n.real.cells,
        ncol = 1
      )
    )

    colnames(pANN) <- "pANN"
    rownames(pANN) <- real.cells

    list.ind <- list.ind + 1

    for (i in 1:n.real.cells) {

      neighbors <- dist.mat[
        2:(pk.temp + 1),
        i
      ]

      pANN$pANN[i] <-
        length(
          which(neighbors > n.real.cells)
        ) / pk.temp
    }

    sweep.res.list[[list.ind]] <- pANN
  }

  return(sweep.res.list)
}


paramSweep_v2 <- function(
    seu,
    PCs = 1:10,
    sct = FALSE,
    num.cores = 1
) {

  requireNamespace("Seurat")
  requireNamespace("fields")
  requireNamespace("parallel")

  pK <- c(
    0.0005,
    0.001,
    0.005,
    seq(0.01, 0.3, by = 0.01)
  )

  pN <- seq(
    0.05,
    0.3,
    by = 0.05
  )

  min.cells <- round(
    nrow(seu@meta.data) /
      (1 - 0.05) -
      nrow(seu@meta.data)
  )

  pK.test <- round(
    pK * min.cells
  )

  pK <- pK[
    which(pK.test >= 1)
  ]

  orig.commands <- seu@commands

  if (nrow(seu@meta.data) > 10000) {

    real.cells <- rownames(
      seu@meta.data
    )[
      sample(
        1:nrow(seu@meta.data),
        10000,
        replace = FALSE
      )
    ]

    data <- seu@assays$RNA@counts[, real.cells]
    n.real.cells <- ncol(data)

  } else {

    real.cells <- rownames(
      seu@meta.data
    )

    data <- seu@assays$RNA@counts
    n.real.cells <- ncol(data)
  }

  if (num.cores > 1) {

    output2 <- mclapply(
      as.list(1:length(pN)),
      FUN = parallel_paramSweep_v2,
      n.real.cells,
      real.cells,
      pK,
      pN,
      data,
      orig.commands,
      PCs,
      sct,
      mc.cores = num.cores
    )

  } else {

    output2 <- lapply(
      as.list(1:length(pN)),
      FUN = parallel_paramSweep_v2,
      n.real.cells,
      real.cells,
      pK,
      pN,
      data,
      orig.commands,
      PCs,
      sct
    )
  }

  sweep.res.list <- list()
  list.ind <- 0

  for (i in 1:length(output2)) {

    for (j in 1:length(output2[[i]])) {

      list.ind <- list.ind + 1

      sweep.res.list[[list.ind]] <-
        output2[[i]][[j]]
    }
  }

  name.vec <- NULL

  for (j in 1:length(pN)) {

    name.vec <- c(
      name.vec,
      paste(
        "pN",
        pN[j],
        "pK",
        pK,
        sep = "_"
      )
    )
  }

  names(sweep.res.list) <- name.vec

  return(sweep.res.list)
}


doubletFinder_v2 <- function(
    seu,
    PCs,
    pN = 0.25,
    pK,
    nExp,
    reuse.pANN = FALSE,
    sct = FALSE,
    annotations = NULL
) {

  requireNamespace("Seurat")
  requireNamespace("fields")
  requireNamespace("KernSmooth")

  if (reuse.pANN != FALSE) {

    pANN.old <- seu@meta.data[, reuse.pANN]

    classifications <- rep(
      "Singlet",
      length(pANN.old)
    )

    classifications[
      order(
        pANN.old,
        decreasing = TRUE
      )[1:nExp]
    ] <- "Doublet"

    seu@meta.data[
      ,
      paste(
        "DF.classifications",
        pN,
        pK,
        nExp,
        sep = "_"
      )
    ] <- classifications

    return(seu)
  }

  if (reuse.pANN == FALSE) {

    real.cells <- rownames(
      seu@meta.data
    )

    data <- seu@assays$RNA@counts[, real.cells]

    n_real.cells <- length(
      real.cells
    )

    n_doublets <- round(
      n_real.cells /
        (1 - pN) -
        n_real.cells
    )

    real.cells1 <- sample(
      real.cells,
      n_doublets,
      replace = TRUE
    )

    real.cells2 <- sample(
      real.cells,
      n_doublets,
      replace = TRUE
    )

    doublets <- (
      data[, real.cells1] +
        data[, real.cells2]
    ) / 2

    colnames(doublets) <- paste0(
      "X",
      1:n_doublets
    )

    data_wdoublets <- cbind(
      data,
      doublets
    )

    if (!is.null(annotations)) {

      stopifnot(
        typeof(annotations) == "character"
      )

      stopifnot(
        length(annotations) ==
          length(Cells(seu))
      )

      stopifnot(
        !any(is.na(annotations))
      )

      annotations <- factor(
        annotations
      )

      names(annotations) <- Cells(seu)

      doublet_types1 <-
        annotations[real.cells1]

      doublet_types2 <-
        annotations[real.cells2]
    }

    orig.commands <- seu@commands

    if (!sct) {

      seu_wdoublets <- CreateSeuratObject(
        counts = data_wdoublets
      )

      seu_wdoublets <- NormalizeData(
        seu_wdoublets,
        normalization.method =
          orig.commands$NormalizeData.RNA@params$normalization.method,
        scale.factor =
          orig.commands$NormalizeData.RNA@params$scale.factor,
        margin =
          orig.commands$NormalizeData.RNA@params$margin
      )

      seu_wdoublets <- FindVariableFeatures(
        seu_wdoublets,
        selection.method =
          orig.commands$FindVariableFeatures.RNA$selection.method,
        loess.span =
          orig.commands$FindVariableFeatures.RNA$loess.span,
        clip.max =
          orig.commands$FindVariableFeatures.RNA$clip.max,
        mean.function =
          orig.commands$FindVariableFeatures.RNA$mean.function,
        dispersion.function =
          orig.commands$FindVariableFeatures.RNA$dispersion.function,
        num.bin =
          orig.commands$FindVariableFeatures.RNA$num.bin,
        binning.method =
          orig.commands$FindVariableFeatures.RNA$binning.method,
        nfeatures =
          orig.commands$FindVariableFeatures.RNA$nfeatures,
        mean.cutoff =
          orig.commands$FindVariableFeatures.RNA$mean.cutoff,
        dispersion.cutoff =
          orig.commands$FindVariableFeatures.RNA$dispersion.cutoff
      )

      seu_wdoublets <- ScaleData(
        seu_wdoublets,
        features =
          orig.commands$ScaleData.RNA$features,
        model.use =
          orig.commands$ScaleData.RNA$model.use,
        do.scale =
          orig.commands$ScaleData.RNA$do.scale,
        do.center =
          orig.commands$ScaleData.RNA$do.center,
        scale.max =
          orig.commands$ScaleData.RNA$scale.max,
        block.size =
          orig.commands$ScaleData.RNA$block.size,
        min.cells.to.block =
          orig.commands$ScaleData.RNA$min.cells.to.block
      )

      seu_wdoublets <- RunPCA(
        seu_wdoublets,
        features =
          orig.commands$ScaleData.RNA$features,
        npcs = length(PCs),
        rev.pca =
          orig.commands$RunPCA.RNA$rev.pca,
        weight.by.var =
          orig.commands$RunPCA.RNA$weight.by.var,
        verbose = FALSE
      )

      pca.coord <-
        seu_wdoublets@reductions$pca@cell.embeddings[, PCs]

      nCells <-
        nrow(seu_wdoublets@meta.data)

      rm(seu_wdoublets)
      gc()
    }

    if (sct) {

      requireNamespace("sctransform")

      seu_wdoublets <- CreateSeuratObject(
        counts = data_wdoublets
      )

      seu_wdoublets <- SCTransform(
        seu_wdoublets
      )

      seu_wdoublets <- RunPCA(
        seu_wdoublets,
        npcs = length(PCs)
      )

      pca.coord <-
        seu_wdoublets@reductions$pca@cell.embeddings[, PCs]

      nCells <-
        nrow(seu_wdoublets@meta.data)

      rm(seu_wdoublets)
      gc()
    }

    dist.mat <- fields::rdist(
      pca.coord
    )

    pANN <- as.data.frame(
      matrix(
        0L,
        nrow = n_real.cells,
        ncol = 1
      )
    )

    if (!is.null(annotations)) {

      neighbor_types <- as.data.frame(
        matrix(
          0L,
          nrow = n_real.cells,
          ncol = length(
            levels(doublet_types1)
          )
        )
      )
    }

    rownames(pANN) <- real.cells
    colnames(pANN) <- "pANN"

    k <- round(
      nCells * pK
    )

    for (i in 1:n_real.cells) {

      neighbors <- order(
        dist.mat[, i]
      )

      neighbors <- neighbors[
        2:(k + 1)
      ]

      pANN$pANN[i] <-
        length(
          which(neighbors > n_real.cells)
        ) / k

      if (!is.null(annotations)) {

        neighbors_that_are_doublets <-
          neighbors[
            neighbors > n_real.cells
          ]

        if (
          length(neighbors_that_are_doublets) > 0
        ) {

          neighbor_types[i, ] <-
            table(
              doublet_types1[
                neighbors_that_are_doublets -
                  n_real.cells
              ]
            ) +
            table(
              doublet_types2[
                neighbors_that_are_doublets -
                  n_real.cells
              ]
            )

          neighbor_types[i, ] <-
            neighbor_types[i, ] /
            sum(neighbor_types[i, ])

        } else {

          neighbor_types[i, ] <- NA
        }
      }
    }

    classifications <- rep(
      "Singlet",
      n_real.cells
    )

    classifications[
      order(
        pANN$pANN[1:n_real.cells],
        decreasing = TRUE
      )[1:nExp]
    ] <- "Doublet"

    seu@meta.data[
      ,
      paste(
        "pANN",
        pN,
        pK,
        nExp,
        sep = "_"
      )
    ] <- pANN[
      rownames(seu@meta.data),
      1
    ]

    seu@meta.data[
      ,
      paste(
        "DF.classifications",
        pN,
        pK,
        nExp,
        sep = "_"
      )
    ] <- classifications

    if (!is.null(annotations)) {

      colnames(neighbor_types) <-
        levels(doublet_types1)

      for (ct in levels(doublet_types1)) {

        seu@meta.data[
          ,
          paste(
            "DF.doublet.contributors",
            pN,
            pK,
            nExp,
            ct,
            sep = "_"
          )
        ] <- neighbor_types[, ct]
      }
    }

    return(seu)
  }
}
