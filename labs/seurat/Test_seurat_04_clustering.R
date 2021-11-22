library(Seurat)
library(cowplot)
library(ggplot2)
library(pheatmap)
library(rafalib)
library(clustree)
setwd("C:/Users/LYdia/Documents/GitHub/workshop-scRNAseq/labs/seurat")



## Input
alldata <- readRDS("../data/results/covid_qc_dr_int.rds")


## Graph clustering

# Building kNN / SNN graph
# check that CCA is still the active assay
alldata@active.assay


alldata <- FindNeighbors(alldata,
                         dims = 1:30,
                         k.param = 60,
                         prune.SNN = 1/15)

alldata_graph <- FindNeighbors(alldata,
                         dims = 1:30,
                         k.param = 60,
                         prune.SNN = 1/15,
                         return.neighbor = TRUE, compute.SNN = TRUE)

# check the names for graphs in the object.
names(alldata@graphs)


hist(alldata@graphs$CCA_snn@x, xlim=c(0,1))


pheatmap(alldata@graphs$CCA_nn[1:200,1:200],
         col=c("white","black"),border_color = "grey90",
         legend = F,cluster_rows = F,cluster_cols = F,fontsize = 2)


# Clustering on a graph

# Clustering with louvain (algorithm 1)
for (res in c( 0.1 , 0.25 , .5 , 1 , 1.5 , 2 )){
  alldata <- FindClusters(alldata, graph.name = "CCA_snn", resolution = res , algorithm = 1)
}

# each time you run clustering, the data is stored in meta data columns:
# seurat_clusters - lastest results only
# CCA_snn_res.XX - for each different resolution you test.

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.0.5")+ggtitle("louvain_0.5"),
  DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.1")+ggtitle("louvain_1"),
  DimPlot(alldata, reduction = "umap", group.by = "CCA_snn_res.2")+ggtitle("louvain_2")
)


clustree(alldata@meta.data, prefix = "CCA_snn_res.")



## K-means clustering
for (k in c( 5 , 7 , 10 , 12 , 15 , 17 , 20)){
  alldata@meta.data[,paste0("kmeans_",k)] <- kmeans(x = alldata@reductions[["pca"]]@cell.embeddings,
                                          centers = k,nstart = 100)$cluster
}

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "umap", group.by = "kmeans_5")+ggtitle("kmeans_5"),
  DimPlot(alldata, reduction = "umap", group.by = "kmeans_10")+ggtitle("kmeans_10"),
  DimPlot(alldata, reduction = "umap", group.by = "kmeans_15")+ggtitle("kmeans_15")
)



clustree(alldata@meta.data, prefix = "kmeans_")



## Hierarchical clustering
# Defining distance between cells
d <- dist( alldata@reductions[["pca"]]@cell.embeddings,
           method="euclidean")

#Compute sample correlations
sample_cor <- cor( Matrix::t(alldata@reductions[["pca"]]@cell.embeddings) )

#Transform the scale from correlations
sample_cor <- (1 - sample_cor) / 2

#Convert it to a distance object
d2 <- as.dist(sample_cor)


d3<-as.dist(1 - cor( Matrix::t(alldata@reductions[["pca"]]@cell.embeddings) ))
range(2*d2-d3)


# Clustering
#euclidean
h_euclidean <- hclust(d, method="ward.D2")

#correlation
h_correlation <- hclust(d2, method="ward.D2")

#euclidean distance
alldata$hc_euclidean_5 <- cutree(h_euclidean,k = 5)
alldata$hc_euclidean_10 <- cutree(h_euclidean,k = 10)
alldata$hc_euclidean_15 <- cutree(h_euclidean,k = 15)

#correlation distance
alldata$hc_corelation_5 <- cutree(h_correlation,k = 5)
alldata$hc_corelation_10 <- cutree(h_correlation,k = 10)
alldata$hc_corelation_15 <- cutree(h_correlation,k = 15)


plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "umap", group.by = "hc_euclidean_5")+ggtitle("hc_euc_5"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_euclidean_10")+ggtitle("hc_euc_10"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_euclidean_15")+ggtitle("hc_euc_15"),
  
  DimPlot(alldata, reduction = "umap", group.by = "hc_corelation_5")+ggtitle("hc_cor_5"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_corelation_10")+ggtitle("hc_cor_10"),
  DimPlot(alldata, reduction = "umap", group.by = "hc_corelation_15")+ggtitle("hc_cor_15")
)



saveRDS(alldata,"../data/results/covid_qc_dr_int_cl.rds")


sessionInfo()

save.image("~/GitHub/workshop-scRNAseq/labs/seurat/Test_seurat_04_clustering.RData")
