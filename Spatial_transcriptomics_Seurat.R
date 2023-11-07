library(ggplot2)
library(Matrix)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(Seurat)
library(dplyr)
library(hdf5r)
library(data.table)
library(patchwork)
library(limma)
sample2 <- Seurat::Load10X_Spatial(
  # The directory contains the read count matrix H5 file and the image data in a subdirectory called `spatial`. 
  data.dir ="M:/PROJECTS/OvCA/snRNA-sequencing/analysis/sample2/my_analysis" , 
  filename ="filtered_feature_bc_matrix.h5",
  assay = "Spatial", # specify name of the initial assay
  slice = "slice1", # specify name of the stored image
  filter.matrix = TRUE, 
  to.upper = FALSE
)
sample2
dim(x = sample2)#genes_vs_spots
#Extract the feature and sample names using rownames or colnames
head(x = rownames(sample2), n = 5)
tail(x = colnames(sample2), n = 5)
class(sample2[[]])
class(sample2[[]])
colnames(sample2[[]])
colnames(sample2[[]])
# nFeature_Spatial: the number of unique genes in each sample
sum(sample2$nFeature_Spatial ==  colSums(sample2@assays$Spatial@counts > 0))
# nCount_Spatial: the total number of detected molecules in each sample
sum(sample2$nCount_Spatial ==  colSums(sample2@assays$Spatial@counts))
# A vector of names of associated objects can be had with the names function
# These can be passed to the double [[ extract operator to pull them from the Seurat object
names(x = sample2)
sample2[['Spatial']] # equivalent to: sample2@assays$Spatial
sample2[['slice1']] # equivalent to: sample2@images$slice1
sample2[['Spatial']]@meta.features
head(sample2[['Spatial']][[]])

#Quality control
#The number of unique genes detected in each sample (nFeature_Spatial).
#The total number of molecules detected within a sample (nCount_Spatial).
#The percentage of reads that map to the mitochondrial genome.
sample2[["percent.mt"]] <- PercentageFeatureSet(sample2, 
                                                   pattern = "^mt-")
VlnPlot(
  sample2, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt"), 
  pt.size = 0.1, ncol = 3) & 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())


#Jointly (rather than separately) consider the QC metrics when filtering
plot1 <- FeatureScatter(
  sample2, feature1 = "nCount_Spatial", feature2 = "percent.mt") + NoLegend()
plot2 <- FeatureScatter(
  sample2, feature1 = "nCount_Spatial", feature2 = "nFeature_Spatial") +
  NoLegend()
plot1 + plot2

sample2[["percent.mt"]] <- PercentageFeatureSet( sample2,pattern="^mt-" )

sample2_subset <- subset(
  sample2, 
  subset = nFeature_Spatial > 300 & 
    nCount_Spatial < 50000 & percent.mt < 30)

print(paste("Filter out", ncol(sample2) - ncol(sample2_subset), 
            "samples because of the outlier QC metrics, with", ncol(sample2_subset),
            "samples left."))
#Normalization
SpatialFeaturePlot(
  sample2_subset, features = c("nFeature_Spatial", "nCount_Spatial", "percent.mt")) &
  theme(legend.position = "bottom")

sample2_norm <- SCTransform(sample2_subset, assay = "Spatial", verbose = FALSE)
names(sample2_norm)
#dimension reduction, data visualization, cluster annotation, differential expression
sample2_obj <- RunPCA(sample2_norm, assay = "SCT", verbose = FALSE)
# compute K nearest neighbors (KNN)
sample2_obj <- FindNeighbors(sample2_obj, reduction = "pca", dims = 1:30)
# Leiden algorithm for community detection
sample2_obj <- FindClusters(sample2_obj, verbose = FALSE)
#PCA result is the default UMAP input, use dimensions 1:30 as input features
sample2_obj <- RunUMAP(sample2_obj, reduction = "pca", dims = 1:30)
plot3 <- DimPlot(sample2_obj, reduction = "umap", label = TRUE) + NoLegend()
plot4 <- SpatialDimPlot(sample2_obj, label = TRUE, label.size = 3) + NoLegend()
plot3 + plot4
SpatialDimPlot(sample2_obj, cells.highlight = CellsByIdentities(object = sample2_obj, idents = c(0,1,2,3,4,5,6,7,8,9,10,11)), facet.highlight = TRUE, ncol = 4)
#identity class of each sample
table(sample2_obj@active.ident)
#find all markers of cluster 9
clusters_markers <- FindMarkers(sample2_obj, ident.1 = 11, min.pct = 0.25)
write.table(clusters_markers,file="cluster11_markers.txt",sep="\t",col.names=NA)

VlnPlot(sample2_obj, features = c("PTPRC", "CCL19", "CXCL13","MS4A1"))
VlnPlot(sample2_obj, features = c("MS4A1","CXCL13","CD8A","LAMP3","NCAM1","CD68"))

SpatialFeaturePlot(object = sample2_obj, 
                   features = rownames(cluster0_markers)[1:6], 
                   alpha = c(0.1, 1), ncol = 3)


SpatialFeaturePlot(object = sample2_obj, 
                   features = c("MS4A1","CXCL13","CD8A","LAMP3","NCAM1","CD68"), 
                   alpha = c(0.2, 1),n=3)
SpatialFeaturePlot(object = sample2_obj, 
                   features = c("FOXP3","PDCD1","CD274","HAVCR2","LAG3","CTLA4","TIGIT","LILRB1"), 
                   alpha = c(0.2, 1),n=4)

#find markers for every cluster compared to all remaining cells
sample2_obj_markers <- FindAllMarkers(sample2_obj, only.pos = TRUE, min.pct = 0.25, 
                                    logfc.threshold =0.25)
write.table(sample2_obj_markers,file="sample2_obj_markers_log.fc.treshold0.25.txt",sep="\t",col.names=NA)

top10 <- sample2_obj_markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)
DoHeatmap(sample2_obj, features = top10$gene) +theme(axis.text.y = element_text(size = 5))

sample2_moransi <- FindSpatiallyVariableFeatures(
  sample2_obj, assay = "SCT", 
  features = VariableFeatures(sample2_obj)[1:10],
    selection.method = "moransi") 
top_features_moransi <- head(
  SpatiallyVariableFeatures(sample2_moransi, 
                            selection.method = "moransi"), 3)
SpatialFeaturePlot(sample2_moransi, 
                   features = top_features_moransi, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
  title = "Top 3 genes with the largest Moran's I",
  subtitle = "among 10 top variable genes for illustration purposes")


sample2_variogram <- FindSpatiallyVariableFeatures(
  sample2_obj, assay = "SCT", 
  features = VariableFeatures(sample2_obj)[1:10],
    selection.method = "markvariogram")  
variogram_output_df <- sample2_variogram@assays$SCT@meta.features %>%
  na.exclude # there are NA rows b/c we only calculated the variogram for 10 genes
head(variogram_output_df[order(variogram_output_df$r.metric.5), ])

top_features_variogram <- head(
  SpatiallyVariableFeatures(sample2_variogram, 
                            selection.method = "markvariogram"), 3)
SpatialFeaturePlot(sample2_variogram, 
                   features = top_features_variogram, ncol = 3, alpha = c(0.1, 1)) + 
  plot_annotation(
  title = "3 genes with the top spatially variable rank (by mark-variogram)",
  subtitle = "among 10 top variable genes for illustration purposes")
