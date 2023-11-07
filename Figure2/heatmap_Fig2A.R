library(ComplexHeatmap)
library(circlize)
#gene expression upload.Genes in rows, samples in columns
my_data=read.table(file="heatmap_Fig2A.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
mat = as.matrix(my_data)
mat_scaled = t(apply(mat, 1, scale))
colnames(mat_scaled)=colnames(mat)
#information about annotation.Each annotation in columnns
pheno=read.table(file="col_annot_Fig2A.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
pheno_row=read.table(file="row_annot_Fig2A.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
#Set annotation
ann <- data.frame(pheno$Categories)
colnames(ann) <- c("Categories")
colours <- list("Categories"=c("1"="black","2"="green","3"="blue","4"="red"))
colAnn <- HeatmapAnnotation(df=ann, which="col", col=colours, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"),simple_anno_size = unit(0.5, "cm"))

ann_row <- data.frame(pheno_row$Immunotype)
colnames(ann_row) <- c("Immunotype")
colours_IM <- list(Immunotype=c("B cells"="#CCFF33","cytotoxicity"="#333300","DC"="#FFFF00","Immune population"="#990000","Immunosupression"="#FF0000","NK cells"="#660066","T cell activation"="#66FFFF","TEM"="royal blue","Th1"="#3300FF","Th2"="#6633FF"))
rowAnn <- HeatmapAnnotation(df=ann_row, which="row",col=colours_IM)

set.seed(123456)
hmap <- Heatmap(mat_scaled,name = "expression",col = colorRamp2(c(-3, 0, 3),c("green", "black", "red")),show_row_names = T,show_column_names =F,cluster_rows =F,cluster_columns =F,show_column_dend = TRUE,show_row_dend = TRUE,row_dend_reorder = TRUE,column_dend_reorder = TRUE,clustering_method_rows = "ward.D2",clustering_method_columns = "ward.D2",top_annotation=colAnn,row_names_gp = grid::gpar(fontsize = 5),na_col = "grey",left_annotation = rowAnn)
draw(hmap, heatmap_legend_side="left", annotation_legend_side="right")
