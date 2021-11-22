setwd("C:/Users/LYdia/Documents/GitHub/workshop-scRNAseq/labs/trajectory")


### Downloading dataset

```{bash, eval=F}
#Create data folder
mkdir data

#Download file from NCBI server into data folder
curl ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE72nnn/GSE72857/suppl/GSE72857%5Fumitab%2Etxt%2Egz -o data/GSE72857_umitab.txt.gz

#Decompress it
gunzip GSE72857_series_matrix.txt.gz
```


### Loading matrix into R and converting it to compressed `.rds` format.
data <- read.delim("data/umitab.txt", header = T, row.names = 1)
comp_matrix <- Matrix::Matrix(as.matrix(data),sparse = T)
saveRDS(comp_matrix,"data/GSE72857_umitab.rds")



### Loading data
umi_counts <- readRDS("data/GSE72857_umitab.rds")
dim(umi_counts)
# 27297 10368
umi_counts <- umi_counts[ , c(T,F,F,F,F)]
dim(umi_counts)
# 27297  2074


#Define a color pallete to use
pal <- c( RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Set2"))



### Basic processing with Seurat pipeline
suppressPackageStartupMessages({
library(Seurat)
library(cowplot)
})

#Data analysis with Seurat pipeline
data <- CreateSeuratObject(counts = umi_counts)
data <- NormalizeData(data)
data <- FindVariableFeatures(data, nfeatures = 2000)
data <- ScaleData(data)
data <- RunPCA(data)
data <- FindNeighbors(data)
data <- FindClusters(data,resolution = 1)
data <- RunUMAP(data, n.neighbors = 10, dims = 1:50,spread = 2,min.dist = 0.3 )

#Plot the clusters
DimPlot(data, group.by = "RNA_snn_res.1")
FeaturePlot(data, reduction = "umap", c("Adgr6","Ctsg","Cebpe","Irf8","Hba-a2"),order = T)

#Save the objects as separate matrices for input in slingshot
dimred <- data@reductions$umap@cell.embeddings
clustering <- data$RNA_snn_res.1
counts <- as.matrix( data@assays$RNA@counts[ data@assays$RNA@var.features , ] )


#### Trajectory inference with Slingshot

#### Defining cell lineages with Slingshot
suppressPackageStartupMessages({
library(slingshot)})

#Run default Slingshot lineage identification
set.seed(1)
lineages <- getLineages(data = dimred,
                     clusterLabels = clustering)
identical(lineages@elementMetadata@listData$reducedDim, dimred)

lineages
lineages@metadata$lineages
# $Lineage1
# [1] "11" "2"  "5"  "10" "12" "0"  "1"  "8"  "4" 
# 
# $Lineage2
# [1] "11" "2"  "5"  "10" "12" "0"  "1"  "8"  "7" 
# 
# $Lineage3
# [1] "11" "2"  "5"  "10" "12" "0"  "1"  "9" 
# 
# $Lineage4
# [1] "11" "2"  "5"  "10" "12" "6" 
# 
# $Lineage5
# [1] "11" "2"  "5"  "3" 


#Plot the lineages
par(mfrow=c(1,2))
plot(dimred[,1:2], col = pal[clustering],  cex=.5,pch = 16)
for(i in levels(clustering)){ 
text( mean(dimred[clustering==i,1]),
     mean(dimred[clustering==i,2]), labels = i,font = 2) }
plot(dimred[,1:2], col = pal[clustering],cex=.5, pch = 16)
lines(as.SlingshotDataSet(lineages), lwd = 3, col = 'black')


#Run default Slingshot
set.seed(1)
lineages <- getLineages(data = dimred,
                     clusterLabels = clustering,
                     #end.clus = c("11","7","10","9","5"), #define how many branches/lineages to consider
                     start.clus = "0") #define where to start the trajectories

lineages
lineages@metadata$lineages
# $Lineage1
# [1] "0"  "12" "10" "5"  "2"  "11"
# 
# $Lineage2
# [1] "0"  "12" "10" "5"  "3" 
# 
# $Lineage3
# [1] "0" "1" "8" "4"
# 
# $Lineage4
# [1] "0" "1" "8" "7"
# 
# $Lineage5
# [1] "0" "1" "9"
# 
# $Lineage6
# [1] "0"  "12" "6" 


#Plot the lineages
par(mfrow=c(1,2))
plot(dimred[,1:2], col = pal[clustering],  cex=.5,pch = 16)
for(i in levels(clustering)){ 
text( mean(dimred[clustering==i,1]),
     mean(dimred[clustering==i,2]), labels = i,font = 2) }
plot(dimred, col = pal[clustering],cex=.5,  pch = 16)
lines(as.SlingshotDataSet(lineages), lwd = 3, col = 'black')



#### Defining Principal Curves
curves <- getCurves(lineages, approx_points = 300, thresh = 0.01, stretch = .8, allow.breaks = FALSE, shrink=.99)
curves

plot(dimred, col = pal[clustering], asp = 1,cex=.5, pch = 16)
lines(as.SlingshotDataSet(curves), lwd = 3, col = 'black')



### Finding differentially expressed genes
BiocParallel::register(BiocParallel::SerialParam())
library(tradeSeq)

#Removing some genes to speed up the computations for this tutorial
filt_counts <- counts [ rowSums(counts > 5) > ncol(counts)/100, ] 
dim(filt_counts)
# 215 2074

sce <- fitGAM(  counts = as.matrix(filt_counts),
             sds = curves )

plotGeneCount(curves, filt_counts, clusters = clustering, models = sce)


#Define function to plot
library(dplyr)
# plot_differential_expression <- function(feature_id) {
# feature_id <- pseudotime_association %>% filter(pvalue < 0.05) %>% top_n(1, -waldStat) %>% pull(feature_id)
# cowplot::plot_grid(
# plotGeneCount(curves, filt_counts, gene=feature_id[1], clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
# plotSmoothers(sce, as.matrix(counts), gene = feature_id[1])
# )}




#### Genes that change with pseudotime
pseudotime_association <- associationTest(sce)
pseudotime_association$fdr <- p.adjust(pseudotime_association$pvalue, method = "fdr")
pseudotime_association <- pseudotime_association[ order(pseudotime_association$pvalue), ]
pseudotime_association$feature_id <- rownames(pseudotime_association)


feature_id <- pseudotime_association %>%
filter(pvalue < 0.05) %>%
top_n(1, -waldStat) %>%
pull(feature_id)
# "Gfi1b"
# plot_differential_expression(feature_id)

cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene=feature_id, clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = feature_id)
)


#### Genes that change between two pseudotime points
pseudotime_start_end_association <- startVsEndTest(sce, pseudotimeValues = c(0, 1))
pseudotime_start_end_association$feature_id <- rownames(pseudotime_start_end_association)

feature_id <- pseudotime_start_end_association %>% 
filter(pvalue < 0.05) %>% 
top_n(1, waldStat) %>% 
pull(feature_id)
# "Ermap"

# plot_differential_expression(feature_id)
cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene=feature_id, clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = feature_id)
)


pseudotime_start_end_association2 <- startVsEndTest(sce)
pseudotime_start_end_association2$feature_id <- rownames(pseudotime_start_end_association2)

feature_id2 <- pseudotime_start_end_association2 %>% 
  filter(pvalue < 0.05) %>% 
  top_n(1, waldStat) %>% 
  pull(feature_id)
# "Mpo"

cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene=feature_id2, clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = feature_id2)
)


pseudotime_start_end_association3 <- startVsEndTest(sce, lineages=TRUE, pseudotimeValues = c(0, 1))
pseudotime_start_end_association3$feature_id <- rownames(pseudotime_start_end_association3)

feature_id <- pseudotime_start_end_association3 %>% 
  filter(pvalue_lineage6 < 0.05) %>% 
  top_n(-1, pvalue_lineage6) %>% 
  pull(feature_id)
# "Cd52"

cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene=feature_id, clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = feature_id)
)



#### Genes that are different between lineages
# Different at the end points, using diffEndTest
different_end_association <- diffEndTest(sce)
different_end_association$feature_id <- rownames(different_end_association)
feature_id <- different_end_association %>% 
filter(pvalue < 0.05) %>% 
arrange(desc(waldStat)) %>% 
dplyr::slice(1) %>% 
pull(feature_id)
# "Prtn3"
cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene=feature_id, clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = feature_id)
)


# Different at the branching point, using earlyDETest
branch_point_association <- earlyDETest(sce)
branch_point_association$feature_id <- rownames(branch_point_association)

feature_id <- branch_point_association %>% 
filter(pvalue < 0.05) %>% 
arrange(desc(waldStat)) %>% 
dplyr::slice(1) %>% 
pull(feature_id)
# "Mpo"

cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene=feature_id, clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = feature_id)
)



## Discovering genes with different expression patterns
patternRes <- patternTest(sce)
oPat <- order(patternRes$waldStat, decreasing = TRUE)
head(rownames(patternRes)[oPat])

feature_id<-rownames(patternRes)[oPat][4]
# "Car2"

cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene=feature_id, clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = feature_id)
)



### Example on combining `patternTest` with `diffEndTest` results
library(ggplot2)
different_pattern<-patternRes
patternRes$Gene <- rownames(patternRes)
patternRes$pattern <- patternRes$waldStat
patternRes <- patternRes[, c("Gene", "pattern")]

endRes<-different_end_association
endRes$Gene <- rownames(different_end_association)
endRes$end <- endRes$waldStat
endRes <- endRes[, c("Gene", "end")]

compare <- merge(patternRes, endRes, by = "Gene", all = FALSE)
compare$transientScore <- 
  rank(-compare$end, ties.method = "min")^2 + rank(compare$pattern, ties.method = "random")^2
ggplot(compare, aes(x = log(pattern), y = log(end))) +
  geom_point(aes(col = transientScore)) +
  labs(x = "patternTest Wald Statistic (log scale)",
       y = "diffEndTest Wald Statistic (log scale)") +
  scale_color_continuous(low = "yellow", high = "red") +
  theme_classic()


topTransient <- compare[which.max(compare$transientScore), "Gene"]
# "Atpase6"
cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene=topTransient, clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = topTransient)
)


head(
  compare[order(compare$transientScore, decreasing = TRUE), "Gene"],
  n = 5
)
cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene="Sphk1", clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = "Sphk1")
)



## Early drivers of differentiation
plotGeneCount(curve = curves, counts = counts,
              clusters = apply(slingClusterLabels(curves), 1, which.max),
              models = sce)
earlyDERes <- earlyDETest(sce, knots = c(1, 2))
oEarly <- order(earlyDERes$waldStat, decreasing = TRUE)
head(rownames(earlyDERes)[oEarly])

cowplot::plot_grid(
  plotGeneCount(curves, filt_counts, gene=rownames(earlyDERes)[oEarly][2], clusters = clustering, models = sce)+ ggplot2::theme(legend.position = "none"),
  plotSmoothers(sce, as.matrix(counts), gene = rownames(earlyDERes)[oEarly][2])
)



# Differential expression in large datasets

# testing against fold change threshold of 2
start2 <- startVsEndTest(sce, l2fc = log2(2))
# testing against fold change threshold of 1.5
pat2 <- patternTest(sce, l2fc = log2(1.5))


# Clustering of genes according to their expression pattern

## Extracting fitted values to use with any clustering method
yhat <- predictCells(models = sce, gene = "Irf8")
ysmooth <- predictSmooth(models = sce, gene = "Irf8", nPoints = 40)
table(ysmooth$lineage)
# 1  2  3  4  5  6 
# 40 40 40 40 40 40 
plot(ysmooth[ysmooth$lineage==1,2])
ysmooth2 <- predictSmooth(models = sce, gene = "Irf8", nPoints = 40, tidy=FALSE)



## Clustering using RSEC, clusterExperiment
library(clusterExperiment)
clusterExperiment::listBuiltInFunctions()

nPointsClus <- 20
clusPat <- clusterExpressionPatterns(sce, nPoints = nPointsClus,
                                     genes = rownames(filt_counts)[1:100])

clusterLabels <- primaryCluster(clusPat$rsec)


cUniq <- unique(clusterLabels)
cUniq <- cUniq[!cUniq == -1] # remove unclustered genes
for (xx in cUniq[1:4]) {
  cId <- which(clusterLabels == xx)
  p <- ggplot(data = data.frame(x = 1:nPointsClus,
                                y = rep(range(clusPat$yhatScaled[cId, ]),
                                        nPointsClus / 2)),
              aes(x = x, y = y)) +
    geom_point(alpha = 0) +
    labs(title = paste0("Cluster ", xx),  x = "Pseudotime", y = "Normalized expression") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5))
  for (ii in 1:length(cId)) {
    geneId <- rownames(clusPat$yhatScaled)[cId[ii]]
    p <- p +
      geom_line(data = data.frame(x = rep(1:nPointsClus, 6),
                                  y = clusPat$yhatScaled[geneId, ],
                                  lineage = rep(0:5, each = nPointsClus)),
                aes(col = as.character(lineage), group = lineage), lwd = 1.5)
  }
  # p <- p + guides(color = FALSE) +
  #   scale_color_manual(values = c("orange", "darkseagreen3"),
  #                      breaks = c("0", "1"))  
  print(p)
}




save.image("~/GitHub/workshop-scRNAseq/labs/trajectory/Test_slingshot2021.RData")
