library(ComplexHeatmap)
library(circlize)

#gene expression upload.Genes in rows, samples in columns
my_data=read.table(file="Immune_population.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
mat = as.matrix(my_data)
mat_scaled = t(apply(mat, 1, scale))

#do clustering & name columns by the 
col_annot=read.table(file="col_annot.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
COLsplit <- factor( paste0(col_annot$Categories),levels=c("cat.1","cat.3","cat.2") )

colnames(mat_scaled)=colnames(mat)
pheno=read.table(file="Immune_population_pheno.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
ROWsplit <- factor( paste0(pheno$group),levels=c("Cell population","Cell function","Inhibitory receptor","Plasma cell") )

col_fun=colorRamp2(c(-3, 0, 3),c("green", "black", "red"))
lgd = Legend(col_fun = col_fun,direction="horizontal",title_position="topleft")

set.seed(123456)

hmap <-
 Heatmap(
  mat_scaled,
  row_split=ROWsplit,column_split=COLsplit,
  show_row_names=T,show_column_names=F,
  cluster_rows=F,cluster_columns=T,
  show_column_dend=T,show_row_dend=T,
  row_dend_reorder=T,column_dend_reorder=T,
  clustering_method_rows="ward.D2",clustering_method_columns="ward.D2",
  show_heatmap_legend=F,
  show_parent_dend_line=F,
  row_names_side="left",
  col=col_fun,
  row_names_gp=grid::gpar( fontsize=10 ),
column_names_gp=grid::gpar( fontsize=5 ),na_col="grey"
 )

draw( hmap,heatmap_legend_side="top",annotation_legend_side="left" )
draw(lgd, x = unit(0.6, "cm"), y = unit(17, "cm"), just = c("left", "top"))


