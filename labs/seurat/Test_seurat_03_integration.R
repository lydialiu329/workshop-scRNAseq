library(Seurat)
library(cowplot)
library(ggplot2)
setwd("C:/Users/LYdia/Documents/GitHub/workshop-scRNAseq/labs/seurat")


alldata <- readRDS("../data/results/covid_qc_dr.rds")
print(names(alldata@reductions))


alldata.list <- SplitObject(alldata, split.by = "orig.ident")
str(alldata.list, max.level = 3)

for (i in 1:length(alldata.list)) {
    alldata.list[[i]] <- NormalizeData(alldata.list[[i]], verbose = FALSE)
    alldata.list[[i]] <- FindVariableFeatures(alldata.list[[i]], selection.method = "vst", nfeatures = 2000,verbose = FALSE)
}

hvgs_per_dataset <- lapply(alldata.list, function(x) { x@assays$RNA@var.features })
venn::venn(hvgs_per_dataset,opacity = .4,zcolor = scales::hue_pal()(3),cexsn = 1,cexil = 1,lwd=1,col="white",frame=F,borders = NA)
venn::venn(hvgs_per_dataset,opacity = .4,zcolor = scales::hue_pal()(3),sncs = 1,ilcs = 1,lwd=1,col="white",box=F,borders = NA)
venn::venn(hvgs_per_dataset)



alldata.anchors <- FindIntegrationAnchors(object.list = alldata.list, dims = 1:30,reduction = "cca")

alldata.anchorsFeatures <- SelectIntegrationFeatures(object.list = alldata.list)
identical(alldata.anchorsFeatures, alldata.anchors@anchor.features)


alldata.int <- IntegrateData(anchorset = alldata.anchors, dims = 1:30, new.assay.name = "CCA")
identical(alldata.int@assays$CCA@var.features, alldata.anchors@anchor.features)

length(paste(alldata.anchors@anchors$dataset1, alldata.anchors@anchors$cell1, sep="-"))
# 60680
length(unique(paste(alldata.anchors@anchors$dataset1, alldata.anchors@anchors$cell1, sep="-")))
# 5413



names(alldata.int@assays)

# by default, Seurat now sets the integrated assay as the default assay, so any operation you now perform will be on the ingegrated data.
alldata.int@active.assay



#Run Dimensionality reduction on integrated space
alldata.int <- ScaleData(alldata.int, verbose = FALSE)
alldata.int <- RunPCA(alldata.int, npcs = 30, verbose = FALSE)
alldata.int <- RunUMAP(alldata.int, dims = 1:30)
alldata.int <- RunTSNE(alldata.int, dims = 1:30)

plot_grid(ncol = 3,
  DimPlot(alldata, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA raw_data"),
  DimPlot(alldata, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE raw_data"),
  DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP raw_data"),
  
  DimPlot(alldata.int, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA integrated"),
  DimPlot(alldata.int, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE integrated"),
  DimPlot(alldata.int, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP integrated")
)


FeaturePlot(alldata.int, 
            reduction = "umap", 
            features = c("CD3E","CD4","CD8A","NKG7",
                         "GNLY","MS4A1","CD14","LYZ",
                         "MS4A7","FCGR3A","CST3","FCER1A"),
            order = T,slot = "data",combine = T)




## another integration method
library(harmony)

alldata.harmony <- RunHarmony(
  alldata,
  group.by.vars = "orig.ident",
  reduction = "pca",
  dims.use = 1:50,
  assay.use = "RNA")



alldata.harmony <- RunUMAP(alldata.harmony, dims = 1:50, reduction = "harmony", reduction.name = "umap_harmony")




#Here we use all PCs computed from Harmony for UMAP calculation
alldata.int[["harmony"]] <- alldata.harmony[["harmony"]]
alldata.int <- RunUMAP(alldata.int, dims = 1:50, reduction = "harmony", reduction.name = "umap_harmony")

identical(alldata.int@reductions$umap_harmony, alldata.harmony@reductions$umap_harmony)


hvgs <- unique(unlist(hvgs_per_dataset))

assaylist <- list()
genelist <- list()
for(i in 1:length(alldata.list)) {
  assaylist[[i]] <- t(as.matrix(GetAssayData(alldata.list[[i]], "data")[hvgs,]))
  genelist[[i]] <- hvgs
}

lapply(assaylist,dim)



## didn't run for the following part
library(reticulate)
scanorama <- import("scanorama")

integrated.data <- scanorama$integrate(datasets_full = assaylist,
                                       genes_list = genelist )

intdimred <- do.call(rbind, integrated.data[[1]])
colnames(intdimred) <- paste0("PC_", 1:100)
rownames(intdimred) <- colnames(alldata.int)

# Add standard deviations in order to draw Elbow Plots in Seurat
stdevs <- apply(intdimred, MARGIN = 2, FUN = sd)

alldata.int[["scanorama"]] <- CreateDimReducObject(
  embeddings = intdimred,
  stdev      = stdevs,
  key        = "PC_",
  assay      = "RNA")

#Here we use all PCs computed from Scanorama for UMAP calculation
alldata.int <- RunUMAP(alldata.int, dims = 1:100, reduction = "scanorama",reduction.name = "umap_scanorama")


p1 <- DimPlot(alldata, reduction = "umap", group.by = "orig.ident")+ggtitle("UMAP raw_data")
p2 <- DimPlot(alldata.int, reduction = "umap", group.by = "orig.ident")+ggtitle("UMAP CCA")
p3 <- DimPlot(alldata.int, reduction = "umap_harmony", group.by = "orig.ident")+ggtitle("UMAP Harmony")
p4 <- DimPlot(alldata.int, reduction = "umap_scanorama", group.by = "orig.ident")+ggtitle("UMAP Scanorama")
leg <- get_legend(p1)

gridExtra::grid.arrange(
  gridExtra::arrangeGrob(
    p1 + NoLegend() + NoAxes(),
    p2 + NoLegend() + NoAxes(),
    p3 + NoLegend() + NoAxes(),
    p4 + NoLegend() + NoAxes(), nrow=2),
  leg, ncol=2,widths=c(8,2)
)


saveRDS(alldata.int,"../data/results/covid_qc_dr_int.rds")



sessionInfo()


save.image("~/GitHub/workshop-scRNAseq/labs/seurat/Test_seurat_03_integration.RData")
