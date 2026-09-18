   #######  Visualise results with a plot
                
                

                                    # all genes

# BP
if (!is.null(go_results_BP) && nrow(as.data.frame(go_results_BP)) > 0) {
barplot <- barplot(go_results_BP, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in BP")
print(barplot)
dotplot <- dotplot(go_results_BP, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in BP") 
print(dotplot)
networkplot <- cnetplot(go_results_BP, showCategory = 5)
print(networkplot)

ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "barplot_GO_analysis_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "dotplot_GO_analysis_BP.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "cnetplot_GO_analysis_BP.png")), width = 12, height = 8, dpi=300)
}

# MF
if (!is.null(go_results_MF) && nrow(as.data.frame(go_results_MF)) > 0) {
barplot <- barplot(go_results_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF")
print(barplot)
dotplot <- dotplot(go_results_MF, showCategory = 10, title = "GO Enrichment (Top 10 Terms) in MF")
print(dotplot)
networkplot <- cnetplot(go_results_MF, showCategory = 5)
print(networkplot)


ggsave(plot = barplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "barplot_GO_analysis_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = dotplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "dotplot_GO_analysis_MF.png")), width = 12, height = 8, dpi=300)
ggsave(plot = networkplot, file= file.path(PATH_ANALYSIS_OUTPUT, paste0("Module_", module_nb), paste0(outputFilesPrefix, "cnetplot_GO_analysis_MF.png")), width = 12, height = 8, dpi=300)
}

