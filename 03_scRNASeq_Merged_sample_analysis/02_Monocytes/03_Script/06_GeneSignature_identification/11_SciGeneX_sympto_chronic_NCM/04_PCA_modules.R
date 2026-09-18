

# PCA on modules
module_cols <- grep("^Module_", colnames(sc10x_subset@meta.data), value = TRUE)

library(dplyr)

patient_modules <- sc10x_subset@meta.data %>%
  group_by(orig.ident, condition) %>%
  summarise(across(all_of(module_cols), median))
  
mat <- as.matrix(patient_modules[, module_cols])

mat_filtered <- mat[, colSums(is.na(mat)) == 0]


pca <- prcomp(mat_filtered, scale. = TRUE)
pca$rotation
sort(pca$rotation[,3], decreasing = TRUE)
sort(pca$rotation[,3], decreasing = FALSE)


# Create dataframe for plotting
pca_df <- data.frame(
  PC3 = pca$x[,3],
  PC4 = pca$x[,4],
  individual = patient_modules$orig.ident,
  condition = patient_modules$condition
)
pca_df$condition <- factor(pca_df$condition)

my_colors <- c(
  "Chronic" = "orange",
  "Symptomatic" = "red"
)

# Plot
library(ggplot2)

pca_plot <- ggplot(pca_df, aes(x = PC3, y = PC4,
                   color = condition,
                   label = individual)) +
  geom_point(size = 5) +
  geom_text(vjust = -1, size = 4) +
  scale_color_manual(values = my_colors) +
  theme_classic() +
  labs(
    title = "PCA of module scores",
    x = paste0("PC3 (", round(summary(pca)$importance[2,3] * 100, 1), "%)"),
    y = paste0("PC4 (", round(summary(pca)$importance[2,4] * 100, 1), "%)")
  )
ggsave(pca_plot , file = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_plot_PCA_PC3_PC4.png")), width = 23, height = 15, units ="cm" )


# Plotting PCA loadings

library(ggplot2)
library(ggrepel)

loadingsMatrix <- pca$rotation

for (dims in combn(1:4, 2, simplify = FALSE)) {

  pc_names <- paste0("PC", dims)

  loadings_sub <- loadingsMatrix[, pc_names]

  top_features <- head(
    rownames(loadings_sub)[
      order(
        rowMeans(abs(loadings_sub)),
        decreasing = TRUE
      )
    ],
    20
  )

  loadingsDF <- data.frame(
    loadings_sub[top_features, ],
    feature = top_features
  )

  axesLimit <- max(abs(loadingsDF[, pc_names]))

  p <- ggplot(
    loadingsDF,
    aes_string(x = pc_names[1], y = pc_names[2])
  ) +
    coord_cartesian(
      xlim = c(-axesLimit, axesLimit),
      ylim = c(-axesLimit, axesLimit)
    ) +
    geom_segment(
      aes_string(
        xend = pc_names[1],
        yend = pc_names[2]
      ),
      x = 0,
      y = 0,
      arrow = arrow(length = unit(2, "mm")),
      colour = "grey50"
    ) +
    geom_text_repel(aes(label = feature)) +
    theme_classic() +
    ggtitle(
      paste(pc_names[1], "vs", pc_names[2])
    )

  print(p)
}

# Extracting top contributing modules for each PC
top_pc_features <- lapply(
  seq_len(ncol(loadingsMatrix)),
  function(pc){

    ord <- order(
      abs(loadingsMatrix[,pc]),
      decreasing = TRUE
    )

    data.frame(
      PC = paste0("PC", pc),
      Feature = rownames(loadingsMatrix)[ord],
      Loading = loadingsMatrix[ord, pc]
    )
  }
)

top_pc_features <- do.call(rbind, top_pc_features)

write.table(
  top_pc_features,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix,
           "_PC_loadings.tsv")
  ),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

# Getting the angle between each module and a PC:
loadings <- pca$rotation

eig <- pca$sdev^2

cos_theta <- sweep(
  loadings,
  2,
  sqrt(eig),
  "*"
)

angle_deg <- acos(
  pmax(
    pmin(cos_theta,1),
    -1
  )
) * 180/pi

table_cos_theta <- lapply(
  seq_len(ncol(loadings)),
  function(pc){

    ord <- order(
      abs(cos_theta[,pc]),
      decreasing = TRUE
    )

    data.frame(
      PC = paste0("PC", pc),
      Module = rownames(loadings)[ord],
      Loading = loadings[ord, pc],
      CosTheta = cos_theta[ord, pc],
      Angle_deg = angle_deg[ord, pc]
    )
  }
)

table_cos_theta <- do.call(rbind, table_cos_theta)

write.table(
  table_cos_theta,
  file = file.path(
    PATH_ANALYSIS_OUTPUT,
    paste0(outputFilesPrefix,
           "_PC_loadings_cos_theta.tsv")
  ),
  sep = "\t",
  row.names = TRUE,
  quote = FALSE
)









# To plot module score by condition
pc4 <- pca$rotation[,4]

sort(pc4)
sort(pc4, decreasing = TRUE)
pc4 <- pca$rotation[,4]

sort(pc4)
sort(pc4, decreasing = TRUE)

library(ggplot2)

for(m in module_cols){

  print(
    ggplot(patient_modules,
           aes(condition, .data[[m]])) +
      geom_boxplot() +
      geom_jitter(width=0.1)
  )
}






















# Correlation between modules
 cor_mat <- cor(mat_filtered, use="pairwise.complete.obs", method="spearman") # option to handle NAs
  
library(pheatmap)
cor_mat[is.na(cor_mat)] <- 0
correlation_plot <- pheatmap(
  cor_mat,
  silent = TRUE
)
png(file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "correlation_heatmap.png")), width = 1200, height = 1000, res = 150)
grid::grid.draw(correlation_plot$gtable)
dev.off()

# Test with only the two most contributing modules
ggplot(patient_modules, aes(Module_11, Module_5, color = condition)) +
  geom_point(size = 5)
  
  


# To chanage labels for the paper
pca_df$individual <- gsub(".*_", "scRNAseq_", pca_df$individual)
pca_df$individual <- paste0("scRNAseq_", sprintf("%02d", seq_len(nrow(pca_df))))



