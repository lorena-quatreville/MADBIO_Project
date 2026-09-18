BP_name = "adhesion"
score_name = paste0(BP_name, "_Score")
features_name = paste0(score_name, "1")

# Create the output directory
dir.create(file.path(PATH_EXPERIMENT_OUTPUT, ANALYSIS_STEP_NAME, ANALYSIS_NAME, BP_name))


# load the signatures up and down
gene_signature_down <- read.csv(file.path(
        PATH_ANALYSIS_OUTPUT, BP_name,
        paste0(outputFilesPrefix, BP_name, "gene_signature_down.csv")), header = TRUE)

gene_signature_up <- read.csv(file.path(
        PATH_ANALYSIS_OUTPUT, BP_name,
        paste0(outputFilesPrefix, BP_name, "gene_signature_up.csv")), header = TRUE)
        
        
# Up and down separated
# Initialize STRINGdb
string_db <- STRINGdb$new(version="11.5", 
            species=9606, # Homo sapiens
            score_threshold=700, # high confidence interaction
            input_directory="")

# Map your genes to STRING IDs
mapped_up <- string_db$map(data.frame(gene=gene_signature_up), "gene", removeUnmappedRows = TRUE)
mapped_down <- string_db$map(data.frame(gene=gene_signature_down), "gene", removeUnmappedRows = TRUE)

# get the best 200 hits
hits_up = mapped_up$STRING_id[1:200]  
hits_down = mapped_down$STRING_id[1:200] 

# plot the STRING network png 
string_db$plot_network( hits_up)
png(file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "network_up.png")), width = 1000, height = 800, res = 150)
string_db$plot_network(hits_up)
dev.off()
pdf(file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "network_up.pdf")), width = 1000, height = 800)
string_db$plot_network(hits_up)
dev.off()
string_db$plot_network( hits_down)
png(file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "network_down.png")), width = 1000, height = 800, res = 150)
string_db$plot_network(hits_down)
dev.off()
pdf(file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "network_down.pdf")), width = 1000, height = 800)
string_db$plot_network(hits_down)
dev.off()

# plot a protein-protein enrichment graph of the best 1000 hits in order to see how the ppi signal is distributed along the sorted list
string_db$plot_ppi_enrichment(mapped_up$STRING_id[1:1000])
string_db$plot_ppi_enrichment(mapped_down$STRING_id[1:1000])

######## compute enrichment in GO annotations  ########
enrichmentGO_up = string_db$get_enrichment( hits_up, category = "Process", methodMT = "fdr", iea = TRUE )
enrichmentKEGG_up = string_db$get_enrichment( hits_up, category = "KEGG", methodMT = "fdr", iea = TRUE )
head(enrichmentGO_up, n=7)
head(enrichmentKEGG_up, n=7)

# Saving
write.csv(
  enrichmentGO_up,
  file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "enrich_GO_up.csv")),
  row.names = FALSE
)

write.csv(
  enrichmentKEGG_up,
  file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "enrich_KEGG_up.csv"))
)

enrichmentGO_down = string_db$get_enrichment( hits_down, category = "Process", methodMT = "fdr", iea = TRUE )
enrichmentKEGG_down = string_db$get_enrichment( hits_down, category = "KEGG", methodMT = "fdr", iea = TRUE )
head(enrichmentGO_down, n=7)
head(enrichmentKEGG_down, n=7)

# Saving
write.csv(
  enrichmentGO_down,
  file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "enrich_GO_down.csv")),
  row.names = FALSE
)

write.csv(
  enrichmentKEGG_down,
  file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "enrich_KEGG_down.csv"))
)






                                    ### By combining both up and down regulated genes ###

# Combining up and down gene lists                                            
gene_signature_all <- rbind(gene_signature_up, gene_signature_down)

# Map your genes to STRING IDs
mapped_all <- string_db$map(
  data.frame(gene=gene_signature_all),
  "gene",
  removeUnmappedRows=TRUE
)

# First visualization
hits_all = mapped_all$STRING_id[1:200] 
string_db$plot_network( hits_all)
png(file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "network_all.png")), width = 1000, height = 800, res = 150)
string_db$plot_network(hits_all)
dev.off()
pdf(file.path(PATH_ANALYSIS_OUTPUT, BP_name,  paste0(outputFilesPrefix, BP_name, "network_all.pdf")), width = 1000, height = 800)
string_db$plot_network(hits_all)
dev.off()

# Retrieve the network
all_genes_edge <- string_db$get_interactions(mapped_all$STRING_id)


# Create igraph object
all_genes_graph <- graph_from_data_frame(all_genes_edge[,c("from","to")], directed=FALSE)


# Color nodes based on gene list
V(all_genes_graph)$color <- ifelse(
  V(all_genes_graph)$name %in% mapped_down$STRING_id, "Upregulated in chronic",
  ifelse(
    V(all_genes_graph)$name %in% mapped_up$STRING_id, "Upregulated in sympto",
    "gray"
  )
)

# Cet node labels to gene symbols
V(all_genes_graph)$label <- mapped_all$gene[match(V(all_genes_graph)$name, mapped_all$STRING_id)]

# Plot with ggraph
ntwk_all <- ggraph(all_genes_graph, layout = "fr") +
  geom_edge_link(alpha = 0.5) +
  geom_node_point(aes(color=color), size=6) +
  geom_node_text(aes(label=label), repel=TRUE, size=4) +
  theme_void()
  
 
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,  paste0( outputFilesPrefix, BP_name, "network_all_colored.png")), plot = ntwk_all, width = 12, height =8, dpi = 300)
  
######## Find clusters of interactions ########

# get clusters
clustersList = string_db$get_clusters(mapped_all$STRING_id[1:600])

# plot first 4 clusters
png(
  file.path(PATH_ANALYSIS_OUTPUT, BP_name, paste0(outputFilesPrefix, BP_name, "top4_clusters_all.png")),
  width = 2000,
  height = 1600,
  res = 300
)

par(mfrow = c(2,2), mar = c(1,1,1,1))

for(i in seq(1:4)){
  string_db$plot_network(clustersList[[i]])
}

dev.off()



                                    ### Same but by filtering ###

# Loading the dataset fom STRING
string_up = read_tsv(file.path(
        PATH_ANALYSIS_OUTPUT, BP_name,
        "string_interactions_migration_up.tsv"), col_names = TRUE)
string_down = read_tsv(file.path(
        PATH_ANALYSIS_OUTPUT, BP_name,
        "string_interactions_migration_down.tsv"), col_names = TRUE)
string_all = read_tsv(file.path(
        PATH_ANALYSIS_OUTPUT, BP_name,
        "string_interactions_migration_all.tsv"), col_names = TRUE)


# Filtering genes:
filtered_edges <- string_all %>%
  filter(
    (experimentally_determined_interaction >= 0 | database_annotated >= 0)
  )

# Generating the graph
all_genes_graph <- graph_from_data_frame(
  filtered_edges %>% select(from = `#node1`, to = node2,
      experimentally_determined_interaction,
      database_annotated,
      coexpression,
      automated_textmining,
      combined_score
))

# Get gene ID
up_gene_id <- unique(mapped_up$gene)
down_gene_id <- unique(mapped_down$gene)

# Remove genes that are not part of down or up regulated gene list
all_genes_graph <- induced_subgraph(
  all_genes_graph,
  vids = V(all_genes_graph)[V(all_genes_graph)$name %in% c(up_gene_id, down_gene_id)]
)

# Color nodes based on gene list
V(all_genes_graph)$group <- ifelse(
  V(all_genes_graph)$name %in% down_gene_id, "Chronic",
  ifelse(
    V(all_genes_graph)$name %in% up_gene_id, "Symptomatic",
    "Others"
  )
)

# Set node labels to gene symbols
V(all_genes_graph)$label <- V(all_genes_graph)$name


#Plotting
ntwk_all_database <- ggraph(all_genes_graph, layout = "fr") +
  geom_edge_link(aes(alpha = database_annotated), color = "gray50") +
  geom_node_point(aes(color = group), size = 5) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3) +
   scale_color_manual(values = c(
    "Symptomatic" = "red",
    "Chronic" = "orange",
    "Others" = "gray"
  )) +
  theme_void()
 
ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,  paste0( outputFilesPrefix, BP_name, "network_all_database_colored.png")), plot = ntwk_all_database, width = 12, height =8, dpi = 300)
   
 #Plotting
ntwk_all_experiment <- ggraph(all_genes_graph, layout = "fr") +
  geom_edge_link(aes(color = experimentally_determined_interaction), width = 0.5) +
  geom_node_point(aes(color = group), size = 5) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3) +
  scale_edge_color_gradient(
    low = "lightgray",
    high = "red"
  ) +
  scale_color_manual(values = c(
    "Symptomatic" = "red",
    "Chronic" = "orange",
    "Others" = "gray"
  )) +
  theme_void() 

ggsave(filename = file.path( PATH_ANALYSIS_OUTPUT, BP_name,  paste0( outputFilesPrefix, BP_name, "network_all_experiment_colored.png")), plot = ntwk_all_experiment, width = 12, height =8, dpi = 300)         
            
            
            
            
            
            
            
            
                          
  


                                            
