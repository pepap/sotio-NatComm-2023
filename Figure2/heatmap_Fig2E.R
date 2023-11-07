library(ComplexHeatmap)
library(circlize)
#gene expression upload.Genes in rows, samples in columns
my_data=read.table(file="heatmap_Fig2E.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
my_data=log2(1+my_data)
mat = as.matrix(my_data)
mat_scaled = t(apply(mat, 1, scale))
colnames(mat_scaled)=colnames(mat)
set.seed(123456)
hmap <- Heatmap(mat_scaled,name = "expression",col = colorRamp2(c(-3, 0, 3),c("green", "black", "red")),show_row_names = T,show_column_names =F,cluster_rows =F,cluster_columns = TRUE,show_column_dend = TRUE,show_row_dend = TRUE,row_dend_reorder = TRUE,column_dend_reorder = TRUE,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",row_names_gp = grid::gpar(fontsize = 12),na_col = "black")
draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")


