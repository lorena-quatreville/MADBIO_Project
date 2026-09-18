


# Loading sc10x object

## If needed to load the merged sc10x object
sc10x = readRDS(PATH_SEURAT_RDS_FILE)

Idents(sc10x) <- "cluster_renames"

# Subset to only include cells from Individual A and Individual B
sc10x_subset <- subset(
  sc10x,
  subset = multiplet_type == "Non-classical monocyte" &
           orig.ident %in% c(sample_id_A, sample_id_B)
)

sc10x_subset_up <- sc10x_subset


# Loading scigeneX gene modules
# Retrieving the list of gene from SciGeneX
scigenex_module_list <- read.csv(file= file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "_ModuleGenes_list.csv")))
# Extract all genes from column 2
scigenex_gene_list <- unlist(
    strsplit(scigenex_module_list[[2]], split = ",")
)
# Check result
length(scigenex_gene_list)

# Remove extra spaces
scigenex_gene_list <- trimws(scigenex_gene_list)

# Remove duplicates
scigenex_gene_list <- unique(scigenex_gene_list)


# Loading the gene ontology file
ontology_file_path <- file.path(PATH_EXPERIMENT_OUTPUT, "04_Diff_gene_exp_analysis", "06_sympto_vs_chronic_NCmono", "GO_Biological_Process_2026_table_up.csv")

ontology_file <- read.csv(ontology_file_path, header = TRUE)
# Filtering ontology file
ontology_file <- ontology_file[
  ontology_file$Adjusted.P.value < 0.05,]


extract_signature_genes <- function(
    GO_df,
    scigenex_gene_list,
    patterns
){

    selected_GO <- GO_df[
        Reduce(`|`, lapply(
            patterns,
            function(x)
                grepl(x, GO_df$Term, ignore.case = TRUE)
        )),
    ]
    unique(selected_GO)
    print(selected_GO)
    genes <- unique(
    unlist(strsplit(selected_GO$Genes, ";"))
)
    print(length(genes))
    print(genes)
    #genes <- intersect(scigenex_gene_list, genes)
    genes <- genes[!grepl("^[0-9.]+$", genes)]

    return(genes)
}

acute_inflammatory_activation_patterns <- c(
  "inflammatory response",
  "tumor necrosis factor",
  "toll-like receptor",
  "pattern recognition receptor",
  "response to lipopolysaccharide",
  "macrophage activation",
  "innate immune",
   # specific IL1 family only
  "interleukin-1([^0-9]|$)",
  "interleukin-6",
  "nf-kappab",
  "cellular defense response",
  "cytokine production",
  "defense response to bacterium"
)


inflammatory_activation_patterns <- c(
  "inflammatory response",
  "tumor necrosis factor",
  "toll-like receptor",
  "pattern recognition receptor",
  "response to lipopolysaccharide",
  "macrophage activation",
  "innate immune",
  "interleukin-1([^0-9]|$)",
  "interleukin-6",
  "nf-kappab",
  "cytokine production",
  "cellular defense response"
)


inflammatory_activation_genes <- extract_signature_genes(
    GO_df = ontology_file,
    scigenex_gene_list = scigenex_gene_list,
    patterns = inflammatory_activation_patterns
)
length(inflammatory_activation_genes)


chemotaxis_patterns <- c(
  "chemotaxis",
  "migration",
  "extravasation",
  "leukocyte adhesion",
  "cell motility",
  "chemokine-mediated signaling",
  "monocyte chemotaxis",
  "neutrophil migration",
  "granulocyte chemotaxis",
  "cell-cell adhesion"
) # trafficking

chemotaxis_genes <- extract_signature_genes(
    GO_df = ontology_file,
    scigenex_gene_list = scigenex_gene_list,
    patterns = chemotaxis_patterns
)
length(chemotaxis_genes )

inflammasome_patterns <- c(
  "inflammasome",
  "interleukin-1([^0-9]|$)",
  "reactive oxygen",
  "superoxide",
  "respiratory burst",
  "toll-like receptor",
  "pattern recognition receptor",
  "lipopolysaccharide",
  "bacterial lipopeptide",
  "nf-kappab"
)

oxidative_antimicrobial_patterns <- c(
  "reactive oxygen",
  "superoxide",
  "respiratory burst",
  "defense response to bacterium",
  "antimicrobial",
  "cell killing",
  "phagocytosis",
  "myeloid leukocyte mediated immunity"
)

oxidative_antimicrobial_genes <- extract_signature_genes(
    GO_df = ontology_file,
    scigenex_gene_list = scigenex_gene_list,
    patterns = oxidative_antimicrobial_patterns
)
length(oxidative_antimicrobial_genes )


intersect(chemotaxis_genes,inflammatory_activation_genes)  
length(intersect(
  inflammatory_activation_genes,
  oxidative_antimicrobial_genes
))      
        
length(inflammatory_activation_genes)        
length(chemotaxis_genes )       
length(oxidative_antimicrobial_genes )
     
BP_name = "chemotaxis"
#BP_name = "inflammatory_activation"
#BP_name= "oxidative_antimicrobial"# NOT USED
#BP_name = "inflammasome_ros"# NOT USED
final_filtered_genes = chemotaxis_genes

score_name = paste0(BP_name, "_Score")
features_name = paste0(score_name, "1")


# Create the output directory
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME, BP_name))

# Saving the list in a csv file
write.csv(data.frame(gene = final_filtered_genes),
          file = file.path(
        PATH_ANALYSIS_OUTPUT, BP_name,
        paste0(outputFilesPrefix, BP_name, "gene_signature_up.csv")), row.names = FALSE)
        

# Compute module scores
sc10x_subset_up <- AddModuleScore(sc10x_subset_up, features = list(final_filtered_genes), name = score_name)

Seurat::DoHeatmap(
  sc10x_subset_up,
  features = final_filtered_genes,
  group.by = "cluster_renames"
)

# to change sample order 
sample_order <- c("12_scRNASeq_A86_S_tp2", "15_scRNASeq_A76_S_tp2", "04_scRNASeq_A81_AS_tp2", "16_scRNASeq_A83_AS_tp2", "13_scRNASeq_A62_AS_tp2", "14_scRNASeq_A75_AS_tp2" ) 
# to order individuals for the plotting 
sc10x_subset_up$orig.ident <- factor(
  sc10x_subset_up$orig.ident,
  levels = sample_order
)


# Global heatmap of the genes contained in the signature

heatmap_plot <- DoHeatmap(sc10x_subset_up, features = final_filtered_genes, group.by = "orig.ident") +
  ggtitle(paste0("Heatmap of Genes contained in ", BP_name, " signature")) +
  theme(axis.text.y = element_text(size = 8))  # smaller text for readability

print(heatmap_plot)

# Save to file
ggsave(
  plot = heatmap_plot,
  filename = file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, "_", BP_name, "_signature_Heatmap_up_reg_sympto.png")),
  width = 10, height = 8, dpi = 300
)



## Displaying the UMAP with the scores

# Plot UMAP FeaturePlot for the scores
umap_score <- FeaturePlot(
  sc10x_subset_up,
  features = features_name
) + 
  scale_color_gradientn(
    colors = c("green", "blue", "purple", "red", "yellow")
  ) +
  ggtitle(score_name)


# Display plots
print(umap_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,   paste0( outputFilesPrefix, "umap_", score_name, "_up.png")), plot = umap_score, width = 12, height =8, dpi = 300)



# Visualisation using violin plots and violin boxplot


# Visualisation of the score accross clusters:

# violin plot:
violin_plot <- VlnPlot(sc10x_subset_up, features = features_name, group.by = "condition") +
  ggtitle(paste0(score_name, " by Condition")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,  paste0( outputFilesPrefix, "vlnPlot_", score_name, "_up.png")), plot = violin_plot, width = 12, height =8, dpi = 300)
print(violin_plot)


# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x_subset_up, features = features_name, group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle(paste0(score_name, " by Condition"))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,   paste0( outputFilesPrefix, "vlnBoxPlot_", score_name, "_up.png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)
print(violin_box_plot)



################## statistical test to compare means (Welsh Anova which assumes normal distribution but unequal variances) ############################

# Test
features_scores <- sc10x_subset_up@meta.data %>%
  filter(condition == "Symptomatic") %>%
  pull(.data[[features_name]])
  # Shapiro–Wilk is very sensitive; With large numbers of cells, it will almost always reject normality; A significant p-value does not automatically invalidate Welch ANOVA
#shapiro.test(features_scores) Not doable because observations shuold not go up to 5000

# Histogram
hist(features_scores, breaks = 50,
     main = "Histogram of the score accross cells",
     xlab = score_name)

# QQplot
qqnorm(features_scores)
qqline(features_scores, col = "red")
# Mild deviations in the tails → usually acceptable; Strong S-shape or heavy skew → non-normal

###### Test for equal variance
#leveneTest(
#  .data[[features_name]] ~ condition,
#  data = sc10x_subset@meta.data,
# center = median   # more robust than mean
#)
# p < 0.05 → variances are significantly different; p ≥ 0.05 → no evidence of unequal variance
# If not significant, we have to use: Wesh anova (correction for unequal variance) and Games-Howell post-hoc test (for unequal variance and sample size


# Visualisation of the score accross  specific clusters:

# Wilcoxon rank sum test instead

wilcox_test <- wilcox.test(
  formula = reformulate("condition", response = features_name),
  data = sc10x_subset_up@meta.data,
  exact = FALSE
)

print(wilcox_test)

p_val <- wilcox_test$p.value
p_label <- paste0("p = ", signif(p_val, 3))

wilcox_df <- data.frame(
  group1 = "Chronic",
  group2 = "Symptomatic",
  p.adj = p_val,
  p_label = p_label
)

wilcox_df <- wilcox_df %>%
  mutate(label = paste0("p = ", scales::pvalue(p.adj)))
  
max_y <- max(sc10x_subset_up@meta.data[[features_name]], na.rm = TRUE)
wilcox_df$y.position <- max_y + 0.05

# Define custom colors
condition_colors <- c(
  "Chronic" = "#E69F00",      # orange
  "Symptomatic" = "#D55E00"   # red
)

# Plot (PAPER)
violin_boxplot_specific_clusters <- ggplot(sc10x_subset_up@meta.data, 
                                           aes(x = condition, y = .data[[features_name]])) +
    geom_violin(aes(fill = condition, color = condition), trim = FALSE, alpha = 0.5, width = 0.8) +
    geom_boxplot(aes(fill = condition, color = condition), 
                 width = 0.2, position = position_dodge(1.5), outlier.shape = NA) +  
    # Mean bar
    stat_summary(
        fun = mean,
        geom = "crossbar",
        width = 0.3,
        color = "black",
        linewidth = 0.6
    ) +
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    stat_pvalue_manual(wilcox_df, label = "p_label", y.position = "y.position", size = 12, step.increase = 0.1) +  
    theme_classic(base_size = 22) +   # removes grey background & grids
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 25),
    axis.text = element_text(size = 25), legend.position = "none") +
    labs(y = score_name, x = "Condition")

print(violin_boxplot_specific_clusters)

# Print and save the plot
print(violin_boxplot_specific_clusters)
ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT, BP_name, paste0(outputFilesPrefix, "vlnBoxPlot_", score_name, "_chronic_sympto_wilcox_up_reg_sympto.png")),
       plot = violin_boxplot_specific_clusters, width = 12, height = 8, dpi = 300)



# to change sample order 
sample_order <- c("12_scRNASeq_A86_S_tp2", "15_scRNASeq_A76_S_tp2", "04_scRNASeq_A81_AS_tp2", "16_scRNASeq_A83_AS_tp2", "13_scRNASeq_A62_AS_tp2", "14_scRNASeq_A75_AS_tp2" ) 
# to order individuals for the plotting 
sc10x_subset_down$orig.ident <- factor(
  sc10x_subset_down$orig.ident,
  levels = sample_order
)


# Global heatmap of the genes contained in the signature

heatmap_plot <- DoHeatmap(sc10x_subset_up, features = final_filtered_genes, group.by = "orig.ident") +
  ggtitle(paste0("Heatmap of Genes contained in ", BP_name, " signature")) +
  theme(axis.text.y = element_text(size = 8))  # smaller text for readability

print(heatmap_plot)

# Save to file
ggsave(
  plot = heatmap_plot,
  filename = file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, "_", BP_name, "_signature_Heatmap.png")),
  width = 10, height = 8, dpi = 300
)












########################## Same for down-regulated genes ################


# Loading the gene ontology file
ontology_file_path <- file.path(PATH_EXPERIMENT_OUTPUT, "04_Diff_gene_exp_analysis", "06_sympto_vs_chronic_NCmono", "GO_Biological_Process_2026_table_down.csv")


ontology_file <- read.csv(ontology_file_path, header = TRUE)

# Filtering ontology file
ontology_file <- ontology_file[
  ontology_file$Adjusted.P.value < 0.05,]

sc10x_subset_down <- sc10x_subset

stress_adaptation_patterns <- c(
  "unfolded protein",
  "endoplasmic reticulum stress",
  "heat",
  "protein folding",
  "integrated stress response",
  "ubiquitin",
  "protein catabolic process",
  "mrna catabolic",
  "response to stress",
  "cellular response to unfolded protein"
) # NOT USED ANYMORE

stress_adaptation_genes <- extract_signature_genes(
    GO_df = ontology_file,
    scigenex_gene_list = scigenex_gene_list,
    patterns = stress_adaptation_patterns
)
length(stress_adaptation_genes)

immune_regulation_patterns <- c(
  "negative regulation of cytokine",
  "negative regulation of inflammatory",
#  "negative regulation of immune",
  "negative regulation of toll-like receptor",
  "negative regulation of nf-kappab",
  "immune response-inhibiting",
#  "regulation of inflammatory response",
#  "regulation of cytokine production",
#  "regulation of tumor necrosis factor",
#  "regulation of interleukin",
#  "regulation of leukocyte degranulation",
#  "regulation of immune system process",
#  "regulation of t cell activation",
#  "regulation of lymphocyte proliferation",
  "tolerance",
  "interleukin-10"
)


immune_regulation_genes <- extract_signature_genes(
    GO_df = ontology_file,
    scigenex_gene_list = scigenex_gene_list,
    patterns = immune_regulation_patterns
)
length(immune_regulation_genes)


endothelial_interaction_patterns <- c(
  "angiogenesis",
  "vasculature",
  "endothelial",
  "blood vessel",
  "vasculogenesis",
  "vasodilation",
  "leukocyte adhesion to vascular endothelial",
  "sprouting angiogenesis"
)

endothelial_interaction_genes <- extract_signature_genes(
    GO_df = ontology_file,
    scigenex_gene_list = scigenex_gene_list,
    patterns = endothelial_interaction_patterns
)
length(endothelial_interaction_genes)





intersect(immune_regulation_genes, endothelial_interaction_genes)
        
        
        
        
length(immune_regulation_genes)        
length(endothelial_interaction_genes)         
        


#BP_name = "endothelial_interaction"
BP_name = "immune_regulation"
score_name = paste0(BP_name, "_Score")
features_name = paste0(score_name, "1")
final_filtered_genes = immune_regulation_genes


# Create the output directory
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME, BP_name))

# Saving the list in a csv file
write.csv(data.frame(gene = final_filtered_genes),
          file = file.path(
        PATH_ANALYSIS_OUTPUT, BP_name,
        paste0(outputFilesPrefix, BP_name, "gene_signature_down.csv")), row.names = FALSE)
        
       
# Compute module scores
sc10x_subset_down <- AddModuleScore(sc10x_subset_down, features = list(final_filtered_genes), name = score_name)


## Displaying the UMAP with the scores

# Plot UMAP FeaturePlot for the scores
umap_score <- FeaturePlot(
  sc10x_subset_down,
  features = features_name
) + 
  scale_color_gradientn(
    colors = c("green", "blue", "purple", "red", "yellow")
  ) +
  ggtitle(score_name)


# Display plots
print(umap_score)

# Save the plot with appropriate size
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,  paste0( outputFilesPrefix, "umap_", score_name, "_down.png")), plot = umap_score, width = 12, height =8, dpi = 300)



# Visualisation using violin plots and violin boxplot


# Visualisation of the score accross clusters:

# violin plot:
violin_plot <- VlnPlot(sc10x_subset_down, features = features_name, group.by = "condition") +
  ggtitle(paste0(score_name, " by Condition")) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name, paste0( outputFilesPrefix, "vlnPlot_", score_name, "_down.png")), plot = violin_plot, width = 12, height =8, dpi = 300)
print(violin_plot)


# violin plot with boxplot:
violin_box_plot <- VlnPlot(sc10x_subset_down, features = features_name, group.by = "condition", pt.size = 0) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  ggtitle(paste0(score_name, " by Condition"))
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,  paste0( outputFilesPrefix, "vlnBoxPlot_", score_name, "_down.png")), plot = violin_box_plot, width = 12, height =8, dpi = 300)
print(violin_box_plot)



################## statistical test to compare means (Welsh Anova which assumes normal distribution but unequal variances) ############################

# Test
features_scores <- sc10x_subset_down@meta.data %>%
  filter(condition == "Symptomatic") %>%
  pull(.data[[features_name]])
  # Shapiro–Wilk is very sensitive; With large numbers of cells, it will almost always reject normality; A significant p-value does not automatically invalidate Welch ANOVA
#shapiro.test(features_scores) Not doable because observations shuold not go up to 5000

# Histogram
hist(features_scores, breaks = 50,
     main = "Histogram of the score accross cells",
     xlab = score_name)

# QQplot
qqnorm(features_scores)
qqline(features_scores, col = "red")
# Mild deviations in the tails → usually acceptable; Strong S-shape or heavy skew → non-normal

###### Test for equal variance
#leveneTest(
#  .data[[features_name]] ~ condition,
#  data = sc10x_subset@meta.data,
# center = median   # more robust than mean
#)
# p < 0.05 → variances are significantly different; p ≥ 0.05 → no evidence of unequal variance
# If not significant, we have to use: Wesh anova (correction for unequal variance) and Games-Howell post-hoc test (for unequal variance and sample size


# Visualisation of the score accross  specific clusters:


# Welch's ANOVA
welch_anova <- oneway.test(
  formula = reformulate("condition", response = features_name),
  data = sc10x_subset_down@meta.data,
  var.equal = FALSE
)
print(welch_anova)


# Ensure custom_group exists
sc10x_subset_down@meta.data$condition <- as.factor(sc10x_subset_down@meta.data$condition)
    
    
# Perform pairwise Games-Howell test # Handles unequal sample size and variances)
games_howell_results <- gamesHowellTest(
  formula = reformulate("condition", response = features_name),
  data = sc10x_subset_down@meta.data
)


# Convert p-value matrix to tidy data frame
pvals <- games_howell_results$p.value

games_howell_df <- as.data.frame(as.table(pvals))
colnames(games_howell_df) <- c("group1", "group2", "p.adj")


# Convert to character before filtering NA
games_howell_df <- games_howell_df %>%
  mutate(group1 = as.character(group1),
         group2 = as.character(group2))%>%
  filter(!is.na(p.adj))


# Add significance label
games_howell_df <- games_howell_df %>%
  mutate(p.adj.signif = p_format(p.adj))

# Remove comparisons that failed or produce NA
games_howell_df <- games_howell_df %>% 
    filter(!is.na(p.adj), !is.na(group1), !is.na(group2))

# Add the y.position column
max_y <- max(sc10x_subset_down@meta.data[[features_name]], na.rm = TRUE)


games_howell_df <- games_howell_df %>%
  mutate(
    p_label = paste0(
      "p = ",
      format(pmax(p.adj, .Machine$double.xmin),
             scientific = TRUE,
             digits = 3)
    )
  )
  
games_howell_df$y.position <- max_y + seq(0.05, by = 0.05, length.out = nrow(games_howell_df))

# Define custom colors
condition_colors <- c(
  "Chronic" = "#E69F00",      # orange
  "Symptomatic" = "#D55E00"   # red
)

# Plot
violin_boxplot_specific_clusters <- ggplot(sc10x_subset_down@meta.data, 
                                           aes(x = condition, y = .data[[features_name]])) +
    geom_violin(aes(fill = condition, color = condition), trim = FALSE, alpha = 0.5, width = 0.8) +
    geom_boxplot(aes(fill = condition, color = condition), 
                 width = 0.2, position = position_dodge(1.5), outlier.shape = NA) +  
    # Mean bar
    stat_summary(
        fun = mean,
        geom = "crossbar",
        width = 0.3,
        color = "black",
        linewidth = 0.6
    ) +
    scale_fill_manual(values = condition_colors) +
    scale_color_manual(values = condition_colors) +
    stat_pvalue_manual(
    games_howell_df,
    label = "p_label",
    y.position = "y.position",
    size = 12,
    step.increase = 0.1
    ) +
    theme(plot.title = element_text(hjust = 0.5), axis.title = element_text(size = 25),
    axis.text = element_text(size = 25), legend.position = "none") +
    theme_classic(base_size = 22) +   # removes grey background & grids
    labs(y = score_name, x = "Condition")

print(violin_boxplot_specific_clusters)

# Print and save the plot
print(violin_boxplot_specific_clusters)
ggsave(filename = file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, "vlnBoxPlot_", score_name, "_ANOVA_chronic_sympto_down.png")),
       plot = violin_boxplot_specific_clusters, width = 12, height = 8, dpi = 300)
       


# to change sample order 
sample_order <- c("12_scRNASeq_A86_S_tp2", "15_scRNASeq_A76_S_tp2", "04_scRNASeq_A81_AS_tp2", "16_scRNASeq_A83_AS_tp2", "13_scRNASeq_A62_AS_tp2", "14_scRNASeq_A75_AS_tp2" ) 
# to order individuals for the plotting 
sc10x_subset_down$orig.ident <- factor(
  sc10x_subset_down$orig.ident,
  levels = sample_order
)


# Global heatmap of the genes contained in the signature

heatmap_plot <- DoHeatmap(sc10x_subset_down, features = final_filtered_genes, group.by = "orig.ident") +
  ggtitle(paste0("Heatmap of Genes contained in ", BP_name, " signature")) +
  theme(axis.text.y = element_text(size = 8))  # smaller text for readability

print(heatmap_plot)

# Save to file
ggsave(
  plot = heatmap_plot,
  filename = file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, "_", BP_name, "_signature_Heatmap.png")),
  width = 10, height = 8, dpi = 300
)




     
     
     
     
     
     
     

   
    ############ En plus
    
    


phago_ros_patterns <- c(
  "phagocytosis",
  "endocytosis",
  "pinocytosis",
  "receptor-mediated endocytosis",
  "lysosomal",
  "superoxide",
  "respiratory burst",
  "reactive oxygen",
  "neutrophil degranulation",
  "antimicrobial",
  "vacuolar acidification",
  "pattern recognition receptor"
)

inflammasome_ros_patterns <- c(
  "inflammasome",
  "interleukin-1",
  "reactive oxygen",
  "superoxide",
  "respiratory burst",
  "toll-like receptor",
  "pattern recognition receptor",
  "lipopolysaccharide",
  "bacterial lipopeptide",
  "myeloid leukocyte mediated immunity",
  "nf-kappab"
)

phago_genes <- extract_signature_genes(
    GO_df = ontology_file,
    scigenex_gene_list = scigenex_gene_list,
    patterns = phago_patterns
)


length(cytokine_genes)
length(chemotaxis_genes)
length(phago_genes)


