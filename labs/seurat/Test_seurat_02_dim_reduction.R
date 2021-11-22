library(Seurat)
library(cowplot)
library(ggplot2)
library(scran)
setwd("C:/Users/LYdia/Documents/GitHub/workshop-scRNAseq/labs/seurat")


alldata <- readRDS("../data/results/seurat_covid_qc.rds")
dim(alldata)
# 18121  5532


## Data preparation
# Feature selection
alldata <- FindVariableFeatures(alldata, selection.method = "vst", nfeatures = 2000 ,verbose = FALSE,assay = "RNA")
range(apply(alldata@assays$RNA@counts, 1, mean)- alldata@assays$RNA@meta.features$vst.mean)
# -8.881784e-16  4.440892e-16
range(apply(alldata@assays$RNA@counts, 1, var)- alldata@assays$RNA@meta.features$vst.variance)
# -4.220055e-10  1.164153e-09
plot(log(alldata@assays$RNA@meta.features$vst.mean), log(alldata@assays$RNA@meta.features$vst.variance))


tapply(alldata@assays$RNA@meta.features$vst.variance.standardized, alldata@assays$RNA@meta.features$vst.variable, range)
# $`FALSE`
# [1] 0.000000 1.216165
# $`TRUE`
# [1]  1.216551 18.493540
sum(alldata@assays$RNA@meta.features$vst.variable)
# 2000

top20 <- head(VariableFeatures(alldata), 20)

LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, repel = TRUE)


# Z-score transformation
alldata <- ScaleData(alldata, vars.to.regress = c("percent_mito", "nFeature_RNA"), assay = "RNA")
# Regressing out percent_mito, nFeature_RNA
# Centering and scaling data matrix
identical(sort(rownames(alldata@assays$RNA@scale.data)), sort(alldata@assays$RNA@var.features))

alldata_scale2<-ScaleData(alldata, vars.to.regress = c("percent_mito", "nFeature_RNA"), features = rownames(alldata), assay = "RNA")
identical(alldata_scale2@assays$RNA@scale.data[rownames(alldata@assays$RNA@scale.data),], alldata@assays$RNA@scale.data)



## PCA
alldata <- RunPCA(alldata, npcs = 50, verbose = F)

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "pca", group.by = "orig.ident",dims = 1:2),
  DimPlot(alldata, reduction = "pca", group.by = "orig.ident",dims = 3:4),
  DimPlot(alldata, reduction = "pca", group.by = "orig.ident",dims = 5:6) )

VizDimLoadings(alldata, dims = 1:5, reduction = "pca",ncol = 5,balanced = T)

ElbowPlot(alldata, reduction = "pca",ndims = 50)



## tSNE
alldata <- RunTSNE(alldata, reduction = "pca", dims = 1:30, 
                   perplexity=30,
                   max_iter=1000,
                   theta=0.5,
                   eta=200,
                   num_threads=0 )

alldata_tsne2<-alldata
alldata_tsne2@reductions<-list()
alldata_tsne2 <- RunTSNE(alldata_tsne2, reduction = "ica", dims = 1:30, 
                   perplexity=30,
                   max_iter=1000,
                   theta=0.5,
                   eta=200,
                   num_threads=0 )

plot_grid(ncol = 3,DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"))



## UMAP
alldata <- RunUMAP(alldata, reduction = "pca", dims = 1:30,
                   n.components=2,
                   n.neighbors=30,
                   n.epochs=200,
                   min.dist=0.3,
                   learning.rate=1,
                   spread=1 )

alldata <- RunUMAP(alldata, reduction.name = "UMAP10_on_PCA",
                   reduction = "pca", 
                   dims = 1:30,
                   n.components=10,
                   n.neighbors=30,
                   n.epochs=200,
                   min.dist=0.3,
                   learning.rate=1,
                   spread=1 )


plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_PCA"),
  DimPlot(alldata, reduction = "UMAP10_on_PCA", group.by = "orig.ident",dims = 1:2)+ ggplot2::ggtitle(label ="UMAP10_on_PCA"),
  DimPlot(alldata, reduction = "UMAP10_on_PCA", group.by = "orig.ident",dims = 3:4)+ ggplot2::ggtitle(label ="UMAP10_on_PCA")
)

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "pca", group.by = "orig.ident"),
  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident"),
  DimPlot(alldata, reduction = "umap", group.by = "orig.ident")
)



## Using ScaledData and graphs for DR
# Using ScaledData for UMAP
alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_ScaleData",
                   features = alldata@assays$RNA@var.features,
                   assay = "RNA",
                   n.components=2,
                   n.neighbors=30,
                   n.epochs=200,
                   min.dist=0.3,
                   learning.rate=1,
                   spread=1 )

plot_grid(ncol = 3,
          DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_PCA"),
          DimPlot(alldata, reduction = "UMAP10_on_PCA", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP10_on_PCA"),
          DimPlot(alldata, reduction = "UMAP_on_ScaleData", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_ScaleData")
)

#Using a Graph for UMAP
alldata <- FindNeighbors(alldata,
                         reduction = "pca",
                         graph.name = "SNN",
                         assay = "RNA",
                         k.param = 20,
                         features = alldata@assays$RNA@var.features)
# doesn't work since R cannot find miniconda installed

alldata <- RunUMAP(alldata, reduction.name = "UMAP_on_Graph",
                   graph = "SNN",
                   assay = "RNA" )


p1 <- DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_PCA")
p2 <- DimPlot(alldata, reduction = "UMAP_on_ScaleData", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_ScaleData")
# p3 <- DimPlot(alldata, reduction = "UMAP_on_Graph", group.by = "orig.ident")+ ggplot2::ggtitle(label ="UMAP_on_Graph")
p3 <- DimPlot(alldata, reduction = "pca", group.by = "orig.ident")+ ggplot2::ggtitle(label ="PCA")
leg <- get_legend(p1)

gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    p1 + NoLegend() + NoAxes(),
    p2 + NoLegend() + NoAxes(),
    p3 + NoLegend() + NoAxes(), 
    leg,nrow=2),
  ncol=1,widths=c(1)
)



## Ploting genes of interest
myfeatures <- c("CD3E","CD4","CD8A","NKG7","GNLY","MS4A1","CD14","LYZ","MS4A7","FCGR3A","CST3","FCER1A")
FeaturePlot(alldata, reduction = "umap",dims = 1:2,
            features = myfeatures,ncol = 3,order = T)




saveRDS(alldata,"../data/results/covid_qc_dr.rds")

sessionInfo()


save.image("~/GitHub/workshop-scRNAseq/labs/seurat/Test_seurat_02_dim_reduction.RData")
