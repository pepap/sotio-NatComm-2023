library(ComplexHeatmap)
library(circlize)

#gene expression upload.Genes in rows, samples in columns
my_data=read.table(file="heatmap.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
mat = as.matrix(my_data)
mat_scaled = t(apply(mat, 1, scale))

#do clustering & name columns by the 
COLsplit <- cluster::pam( hclust( d=dist( x=t(mat_scaled),method="euclidean" ),method="ward.D2" ),k=3 )
COLsplit <- as.character(as.vector( COLsplit$clustering ))
COLsplit[ COLsplit=="3" ] <- "3"; COLsplit[ COLsplit=="2" ] <- "2"; COLsplit[ COLsplit=="1" ] <- "1"
COLsplit <- factor( x=COLsplit,levels=c("1","3","2") )

colnames(mat_scaled)=colnames(mat)
pheno=read.table(file="row_annot.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
ROWsplit <- factor( paste0(pheno$Group),levels=c("Cell population","Cell function","Inhibitory receptor") )

col_fun=colorRamp2(c(-3, 0, 3),c("green", "black", "red"))
lgd = Legend(col_fun = col_fun,direction="horizontal",title_position="topleft")
pheno_col=read.table(file="col_annot.txt",sep="\t",header=TRUE,row.names=1,stringsAsFactors=FALSE,check.names=FALSE)
ann <- data.frame(pheno_col$TMB_med,pheno_col$TLS_med)
colnames(ann) <- c("TMB_med","TLS_med")
colours <- list("TMB_med"=c("High"="red2","Low"="royalblue"),"TLS_med"=c("High"="red2","Low"="royalblue"))
colAnn <- HeatmapAnnotation(df=ann, which="col", col=colours, annotation_width=unit(c(1, 4), "cm"), gap=unit(1, "mm"))

#Set annotation

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
  top_annotation=colAnn,
  show_heatmap_legend=F,
  show_parent_dend_line=F,
  row_names_side="left",
  col=col_fun,
  row_names_gp=grid::gpar( fontsize=10 ),na_col="grey"
 )

draw( hmap,heatmap_legend_side="top",annotation_legend_side="right" )
draw(lgd, x = unit(0.6, "cm"), y = unit(17, "cm"), just = c("left", "top"))

