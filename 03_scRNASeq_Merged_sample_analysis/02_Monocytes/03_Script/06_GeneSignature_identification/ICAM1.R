
plot = VlnPlot(sc10x_subset, "ICAM1", group.by="condition")
sc10x@meta.data$ICAM1 <- FetchData(sc10x, vars = "ICAM1")[,1]
summary(sc10x@meta.data$ICAM1)
t_test_result <- t.test(
     ICAM1 ~ condition,
     data = sc10x@meta.data)
 
t_test_result

plot = VlnPlot(sc10x_subset, "ICAM1", group.by="condition")
ggsave(
filename = file.path(PATH_ANALYSIS_OUTPUT, paste0(outputFilesPrefix, "VlnPlot_ICAM1_t_test.png")),
  plot = plot,
  width = 20,
  height = 12
)

