library(Seurat)
library(venn)
library(dplyr)
library(cowplot)
library(ggplot2)
library(pheatmap)
library(rafalib)
library(scPred)
setwd("C:/Users/LYdia/Documents/GitHub/workshop-scRNAseq/labs/seurat")



## Load and process data
# Covid-19 data
#load the data and select 'ctrl_13` sample
alldata <- readRDS("../data/results/covid_qc_dr_int_cl.rds")
ctrl = alldata[, alldata$orig.ident == 'ctrl_13']

# set active assay to RNA and remove the CCA assay
ctrl@active.assay = 'RNA' 
ctrl[['CCA']] = NULL


# Reference data
reference <- scPred::pbmc_1
reference
unique(reference$orig.ident)
table(reference$cell_type)
# B cell  CD4 T cell  CD8 T cell         cDC       cMono      ncMono     NK cell         pDC Plasma cell 
#    280        1620         945          26         212          79         312          20           6


# Rerun analysis pipeline
reference <- reference %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30)


DimPlot(reference, group.by = "cell_type", label = TRUE, repel = TRUE) + NoAxes()



#Set the identity as louvain with resolution 0.3
ctrl <- SetIdent(ctrl, value = "CCA_snn_res.0.5")
  
ctrl <- ctrl %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData() %>%
  RunPCA(verbose = F) %>%
  RunUMAP(dims = 1:30)

DimPlot(ctrl,  label = TRUE, repel = TRUE) + NoAxes()



## Seurat label transfer
transfer.anchors <- FindTransferAnchors(reference = reference, query = ctrl, 
    dims = 1:30)
predictions <- TransferData(anchorset = transfer.anchors, refdata = reference$cell_type, 
    dims = 1:30)
ctrl <- AddMetaData(object = ctrl, metadata = predictions)

DimPlot(ctrl, group.by = "predicted.id", label = T, repel = T) + NoAxes()


ggplot(ctrl@meta.data, aes(x=CCA_snn_res.0.5, fill = predicted.id)) + geom_bar() + theme_classic()



## scPred
reference <- getFeatureSpace(reference, "cell_type")
identical(reference@misc$scPred@cell_embeddings, reference@reductions$pca@cell.embeddings)
identical(reference@misc$scPred@feature_loadings, reference@reductions$pca@feature.loadings)

reference <- trainModel(reference)
get_probabilities(reference) %>% head()

get_scpred(reference)
plot_probabilities(reference)


ctrl <- scPredict(ctrl, reference)
table(ctrl@meta.data$scpred_prediction)
# B cell  CD4 T cell  CD8 T cell         cDC       cMono      ncMono     NK cell         pDC Plasma cell  unassigned 
#    105         200         233          24         212         109         171           2           3          79 
table(ctrl@meta.data$scpred_no_rejection)
# B cell  CD4 T cell  CD8 T cell         cDC       cMono      ncMono     NK cell         pDC Plasma cell 
#    105         208         240          26         269         111         174           2           3

DimPlot(ctrl, group.by = "scpred_prediction", label = T, repel = T) + NoAxes()

ggplot(ctrl@meta.data, aes(x=CCA_snn_res.0.5, fill = scpred_prediction)) + geom_bar() + theme_classic()




## Compare results
crossTab(ctrl, "predicted.id", "scpred_prediction")




## GSEA with celltype markers
# DEG overlap
# run differential expression in our dataset, using clustering at resolution 0.3
alldata <- SetIdent(alldata,value = "CCA_snn_res.0.5")
DGE_table <- FindAllMarkers(alldata,
                               logfc.threshold = 0,
                               test.use = "wilcox",
                               min.pct = 0.1,
                               min.diff.pct = 0,
                               only.pos = TRUE,
                               max.cells.per.ident = 20,
                               return.thresh = 1,
                               assay = "RNA")

# split into a list
DGE_list <- split(DGE_table, DGE_table$cluster)

unlist(lapply(DGE_list, nrow))


# Compute differential gene expression in reference dataset (that has cell annotation)
reference <- SetIdent( reference, value = "cell_type")
reference_markers <- FindAllMarkers( reference , min.pct = .1 , 
                                     min.diff.pct = .2, only.pos = T, 
                                     max.cells.per.ident = 20 ,
                                     return.thresh = 1 )

# Identify the top cell marker genes in reference dataset
# select top 50 with hihgest foldchange among top 100 signifcant genes.
reference_markers <- reference_markers [ order(reference_markers$avg_logFC,decreasing = T), ]
reference_markers %>% 
  group_by(cluster) %>% 
  top_n(-100, p_val) %>% 
  top_n(50, avg_log2FC) -> top50_cell_selection

# Transform the markers into a list
ref_list = split(top50_cell_selection$gene, top50_cell_selection$cluster)

unlist(lapply(ref_list, length))



suppressPackageStartupMessages(library(fgsea))

# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x){
  fgseaRes <- fgsea( pathways=ref_list, stats=setNames(x$avg_log2FC, x$gene),nperm=10000)
})
names(res) <- names(DGE_list)
unlist(lapply(res, nrow))
# 0 1 2 3 4 5 6 7 8 9 
# 9 7 8 8 9 8 7 6 9 8

# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {x[ x$pval < 0.1 , ]} )
res <- lapply(res, function(x) {x[ x$size > 2 , ]} )
res <- lapply(res, function(x) {x[order(x$NES,decreasing = T), ]} )
unlist(lapply(res, nrow))
# 0 1 2 3 4 5 6 7 8 9 
# 6 3 5 5 3 6 2 2 4 6
res



new.cluster.ids <- unlist(lapply(res,function(x){as.data.frame(x)[1,1]}))
alldata$ref_gsea <- new.cluster.ids[as.character(alldata@active.ident)]

cowplot::plot_grid( ncol = 2,
DimPlot(alldata,label = T,group.by = "CCA_snn_res.0.5") + NoAxes(),
DimPlot(alldata,label = T, group.by = "ref_gsea") + NoAxes())


ctrl$ref_gsea = alldata$ref_gsea[alldata$orig.ident == "ctrl_13"]

cowplot::plot_grid( ncol = 3,
DimPlot(ctrl,label = T,group.by = "ref_gsea") + NoAxes() + ggtitle("GSEA"),
DimPlot(ctrl,label = T, group.by = "predicted.id") + NoAxes() + ggtitle("LabelTransfer"),
DimPlot(ctrl,label = T, group.by = "scpred_prediction") + NoAxes() + ggtitle("scPred")
)



## With annotated gene sets
# Download gene marker list
if(!dir.exists("../data/CellMarker_list/")) {
  dir.create("../data/CellMarker_list")
  download.file(url="http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt",
                destfile = "../data/CellMarker_list/Human_cell_markers.txt")
  download.file(url="http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Mouse_cell_markers.txt",
                destfile = "../data/CellMarker_list/Mouse_cell_markers.txt")
}


# Load the human marker table
markers <- read.delim("../data/CellMarker_list/Human_cell_markers.txt")
markers <- markers [ markers$speciesType == "Human", ]
markers <- markers [ markers$cancerType == "Normal", ]

#Filter by tissue (to reduce computational time and have tissue-specific classification)
# sort(unique(markers$tissueType))
# grep("blood",unique(markers$tissueType),value = T)
# markers <- markers [ markers$tissueType %in% c("Blood","Venous blood",
#                                                "Serum","Plasma",
#                                                "Spleen","Bone marrow","Lymph node"), ]

# remove strange characters etc.
celltype_list <- lapply( unique(markers$cellName) , function(x){
  x <- paste(markers$geneSymbol[markers$cellName == x],sep=",")
  x <- gsub("[[]|[]]| |-",",",x)
  x <- unlist(strsplit( x , split = ","))
  x <- unique(x [ ! x %in% c("","NA","family") ])
  x <- casefold(x,upper = T)
})
names(celltype_list) <- unique(markers$cellName)
# celltype_list <- lapply(celltype_list , function(x) {x[1:min(length(x),50)]} )
celltype_list <- celltype_list[ unlist(lapply(celltype_list,length)) < 100 ]
celltype_list <- celltype_list[ unlist(lapply(celltype_list,length)) > 5 ]


# run fgsea for each of the clusters in the list
res <- lapply(DGE_list, function(x){
  gene_rank <- setNames(x$avg_log2FC, x$gene)
  fgseaRes <- fgsea( pathways=celltype_list, stats=gene_rank,nperm=10000)
  return(fgseaRes)
})
names(res) <- names(DGE_list)

# You can filter and resort the table based on ES, NES or pvalue
res <- lapply(res, function(x) {x[ x$pval < 0.01 , ]} )
res <- lapply(res, function(x) {x[ x$size > 5 , ]} )
res <- lapply(res, function(x) {x[order(x$NES,decreasing = T), ]} )

# show top 3 for each cluster.
lapply(res,head,3)




new.cluster.ids <- unlist(lapply(res,function(x){as.data.frame(x)[1,1]}))
alldata$cellmarker_gsea <- new.cluster.ids[as.character(alldata@active.ident)]

cowplot::plot_grid( ncol = 2,
DimPlot(alldata,label = T,group.by = "ref_gsea") + NoAxes(),
DimPlot(alldata,label = T, group.by = "cellmarker_gsea") + NoAxes()
)




## Save data
saveRDS(ctrl,"../data/results/ctrl13_qc_dr_int_cl_celltype.rds")


sessionInfo()
save.image("~/GitHub/workshop-scRNAseq/labs/seurat/Test_seurat_06_celltype.RData")
