library(Seurat)
library(Matrix)
library(DoubletFinder)
setwd("C:/Users/LYdia/Documents/GitHub/workshop-scRNAseq/labs/seurat")


## Get data
cov.15 <- Seurat::Read10X_h5(
  filename = "../data/covid_data_GSE149689/sub/nCoV_PBMC_15.h5",
  use.names = T)
cov.1 <- Seurat::Read10X_h5(
  filename = "../data/covid_data_GSE149689/sub/nCoV_PBMC_1.h5",
  use.names = T)
cov.17 <- Seurat::Read10X_h5(
  filename = "../data/covid_data_GSE149689/sub/nCoV_PBMC_17.h5",
  use.names = T)

ctrl.5 <- Seurat::Read10X_h5(
  filename = "../data/covid_data_GSE149689/sub/Normal_PBMC_5.h5",
  use.names = T)
ctrl.13 <- Seurat::Read10X_h5(
  filename = "../data/covid_data_GSE149689/sub/Normal_PBMC_13.h5",
  use.names = T)
ctrl.14 <- Seurat::Read10X_h5(
  filename = "../data/covid_data_GSE149689/sub/Normal_PBMC_14.h5",
  use.names = T)


# cov15_matrix<-as.matrix(cov.15)



## Create one merged object
sdata.cov15 <- CreateSeuratObject(cov.15,  project = "covid_15")
sdata.cov1 <- CreateSeuratObject(cov.1,  project = "covid_1")
sdata.cov17 <- CreateSeuratObject(cov.17,  project = "covid_17")
sdata.ctrl5 <- CreateSeuratObject(ctrl.5,  project = "ctrl_5")
sdata.ctrl13 <- CreateSeuratObject(ctrl.13,  project = "ctrl_13")
sdata.ctrl14 <- CreateSeuratObject(ctrl.14,  project = "ctrl_14")

# add metadata
sdata.cov1$type = "Covid"
sdata.cov15$type = "Covid"
sdata.cov17$type = "Covid"
sdata.ctrl5$type = "Ctrl"
sdata.ctrl13$type = "Ctrl"
sdata.ctrl14$type = "Ctrl"



# Merge datasets into one single seurat object
alldata <- merge(sdata.cov15, c(sdata.cov1, sdata.cov17, sdata.ctrl5, sdata.ctrl13, sdata.ctrl14), add.cell.ids=c("covid_15","covid_1","covid_17","ctrl_5","ctrl_13", "ctrl_14"))


# remove all objects that will not be used.
rm(cov.15, cov.1, cov.17, ctrl.5, ctrl.13, ctrl.14, sdata.cov15, sdata.cov1, sdata.cov17, sdata.ctrl5, sdata.ctrl13, sdata.ctrl14)

# run garbage collect to free up memory
gc()




## Calculate QC
# mitochondrial reads 
# Way1: Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^MT-", col.name = "percent_mito")

# Way2: Doing it manually
total_counts_per_cell <- colSums(  alldata@assays$RNA@counts  )
mito_genes <- rownames(alldata)[grep("^MT-",rownames(alldata))]
alldata$percent_mito <- colSums(  alldata@assays$RNA@counts[mito_genes,]  ) / total_counts_per_cell



# ribosomal proteins
# Way1: Doing it using Seurat function
alldata <- PercentageFeatureSet(alldata, "^RP[SL]", col.name = "percent_ribo")

# Way2: Doing it manually
ribo_genes <- rownames(alldata)[grep("^RP[SL]",rownames(alldata))]
alldata$percent_ribo <- colSums(  alldata@assays$RNA@counts[ribo_genes,]  ) / total_counts_per_cell


# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
alldata <- PercentageFeatureSet(alldata, "^HB[^(P)]", col.name = "percent_hb")
rownames(alldata)[grep("^HB[^(P)]",rownames(alldata))]

alldata <- PercentageFeatureSet(alldata, "PECAM1|PF4", col.name = "percent_plat")
rownames(alldata)[grep("PECAM1|PF4",rownames(alldata))]



## Plot QC
feats <- c("nFeature_RNA","nCount_RNA","percent_mito","percent_ribo", "percent_hb")
VlnPlot(alldata, group.by= "orig.ident", features = feats, pt.size = 0.1,ncol = 3) + NoLegend()

FeatureScatter(alldata, "nCount_RNA"  , "nFeature_RNA", group.by = "orig.ident", pt.size = .5)
FeatureScatter(alldata, "percent_mito"  , "nFeature_RNA", group.by = "orig.ident", pt.size = .5)
FeatureScatter(alldata, "percent_mito"  , "nCount_RNA", group.by = "orig.ident", pt.size = .5)


## Filtering

# Detection-based filtering
selected_c <- WhichCells(alldata, expression = nFeature_RNA > 200)
selected_f <- rownames(alldata)[ Matrix::rowSums(alldata) > 3] # in the demo but wrong
# selected_f_1 <- rownames(alldata)[apply(as.matrix(alldata@assays$RNA@counts), 1, function(x) sum(x>0))>3]


data.filt <- subset(alldata, features=selected_f, cells=selected_c)
dim(data.filt)
# 18147  7973

# skip for now and run DoubletFinder first!

#high.det.v3 <- WhichCells(data.filt, expression = nFeature_RNA > 4100)
#high.det.v2 <- WhichCells(data.filt, expression = nFeature_RNA > 2000 & orig.ident == "v2.1k")

# remove these cells
#data.filt <- subset(data.filt, cells=setdiff(WhichCells(data.filt),c(high.det.v2,high.det.v3)))

# check number of cells
# ncol(data.filt)


#Compute the relative expression of each gene per cell
#Use sparse matrix operations, if your dataset is large, doing matrix devisions the regular way will take a very long time.
par(mar=c(4,8,2,1))
C <- data.filt@assays$RNA@counts
C <-  Matrix::t( Matrix::t(C) / Matrix::colSums(C) ) * 100
most_expressed <- order(apply(C,1,median),decreasing = T)[20:1]
boxplot( as.matrix(t(C[most_expressed,])),cex=.1, las=1, xlab="% total count per cell",col=scales::hue_pal()(20)[20:1],horizontal=TRUE)


# Mito/Ribo filtering
selected_mito <- WhichCells(data.filt, expression = percent_mito < 20)
selected_ribo <- WhichCells(data.filt, expression = percent_ribo > 5)

# and subset the object to only keep those cells
data.filt <- subset(data.filt, cells = selected_mito)
data.filt <- subset(data.filt, cells = selected_ribo)

dim(data.filt)
# 18147  5762

table(data.filt$orig.ident)
# covid_1 covid_15 covid_17  ctrl_13  ctrl_14   ctrl_5 
# 878      585     1042     1154     1063     1040 


# Plot filtered QC

feats <- c("nFeature_RNA","nCount_RNA","percent_mito","percent_ribo", "percent_hb")
VlnPlot(data.filt, group.by= "orig.ident", features = feats, pt.size = 0.1,ncol = 3) + NoLegend()


# Filter genes
dim(data.filt)

# Filter MALAT1
data.filt <- data.filt[ ! grepl("MALAT1", rownames(data.filt)), ]

# Filter Mitocondrial
data.filt <- data.filt[ ! grepl("^MT-", rownames(data.filt)), ]

# Filter Ribossomal gene (optional if that is a problem on your data)
# data.filt <- data.filt[ ! grepl("^RP[SL]", rownames(data.filt)), ]

# Filter Hemoglobin gene (optional if that is a problem on your data)
data.filt <- data.filt[ ! grepl("^HB[^(P)]", rownames(data.filt)), ]

dim(data.filt)
# 18121  5762



## Sample sex
genes.file = "../data/results/genes.table.csv"

if (!file.exists(genes.file)){
  suppressMessages(require(biomaRt))

  # initialize connection to mart, may take some time if the sites are unresponsive.
  mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

  # fetch chromosome info plus some other annotations
  genes.table <- try(biomaRt::getBM(attributes = c("ensembl_gene_id",
              "external_gene_name", "description","gene_biotype", "chromosome_name","start_position"),
              mart = mart, useCache = F))
  
  if(!dir.exists("../data/results")){dir.create("../data/results")}
  if(is.data.frame(genes.table)){write.csv(genes.table, file = genes.file)}
  
  if (!file.exists(genes.file)){
  download.file("https://raw.githubusercontent.com/NBISweden/workshop-scRNAseq/master/labs/misc/genes.table.csv",destfile = "data/results/genes.table.csv")
    genes.table = read.csv(genes.file)
    }

} else {
  genes.table = read.csv(genes.file)
}

genes.table <- genes.table[genes.table$external_gene_name %in% rownames(data.filt),]



chrY.gene = genes.table$external_gene_name[genes.table$chromosome_name == "Y"]

data.filt$pct_chrY = colSums(data.filt@assays$RNA@counts[chrY.gene,]) / colSums(data.filt@assays$RNA@counts)

FeatureScatter(data.filt, feature1 = "XIST", feature2 = "pct_chrY")

VlnPlot(data.filt, features = c("XIST", "pct_chrY"))

# remove sex-related genes



## Calculate cell-cycle scores
# Before running CellCycleScoring the data need to be normalized and logtransformed.
data.filt = NormalizeData(data.filt)
data.filt_2 = NormalizeData(data.filt, scale.factor=10000)
table(apply(exp(data.filt_2@assays$RNA@data), 2, sum))
# 28121 
#  5762 
table(apply(expm1(data.filt_2@assays$RNA@data), 2, sum))
# 9999.99999999999            10000 
#                1             5761
data.filt_3 = NormalizeData(data.filt, scale.factor=28121)
table(apply(exp(data.filt_3@assays$RNA@data), 2, sum))
# 46242 
#  5762
# the sum of exp() of each cell is 18121 (# of cells) plus the scale.factor, instead of scale.factor itself
# that's because log1p() instead of log() function is used, to avoid of the disaster of log(0)
data.filt_5 = NormalizeData(data.filt, scale.factor=-8121)
# does not work
data.filt_6 = NormalizeData(data.filt, scale.factor=0)
table(apply(exp(data.filt_6@assays$RNA@data), 2, sum))
# 18121 
#  5762
sum(data.filt_6@assays$RNA@data)
# 0

total_counts_per_cell <- colSums(  data.filt@assays$RNA@counts  )
data.filt_2_normedData <-t( t(data.filt@assays$RNA@counts) / total_counts_per_cell) *10000


data.filt_4 = NormalizeData(data.filt, normalization.method = "RC", scale.factor=28121)
table(apply(data.filt_4@assays$RNA@data, 2, sum))
# 28121 
#  5762
# the sum of each cell is the scale.factor
data.filt_7 = NormalizeData(data.filt, normalization.method = "RC")
table(apply(data.filt_7@assays$RNA@data, 2, sum))
# 10000 
#  5762
identical(data.filt_2_normedData, data.filt_7@assays$RNA@data)
# TRUE
identical(log1p(data.filt_2_normedData), data.filt_2@assays$RNA@data)
# TRUE


data.filt <- CellCycleScoring(object = data.filt,
                              g2m.features = cc.genes$g2m.genes,
                              s.features = cc.genes$s.genes)
# Warning: The following features are not present in the object: MLF1IP, not searching for symbol synonyms
# Warning: The following features are not present in the object: FAM64A, HN1, not searching for symbol synonyms
table(data.filt$Phase)
#   G1  G2M    S 
# 2536 1265 1961

data.filt_cellCycle2 <- CellCycleScoring(object = data.filt,
                                         g2m.features = cc.genes$g2m.genes,
                                         s.features = cc.genes$s.genes)
object.cc <- AddModuleScore(object = data.filt, features = list(S.Score = cc.genes$s.genes, G2M.Score = cc.genes$g2m.genes), 
                            name = "Cell.Cycle", ctrl = 42)
cc.columns <- grep(pattern = "Cell.Cycle", x = colnames(x = object.cc[[]]), 
                   value = TRUE)
cc.scores <- object.cc[[cc.columns]]
assignments <- apply(X = cc.scores, MARGIN = 1, FUN = function(scores, first = "S", second = "G2M", null = "G1") {
  if (all(scores < 0)) {
    return(null)
  }
  else {
    if (length(which(x = scores == max(scores))) > 1) {
      return("Undecided")
    }
    else {
      return(c(first, second)[which(x = scores == max(scores))])
    }
  }
})
cc.scores <- merge(x = cc.scores, y = data.frame(assignments), 
                   by = 0)
colnames(x = cc.scores) <- c("rownames", "S.Score", 
                             "G2M.Score", "Phase")
rownames(x = cc.scores) <- cc.scores$rownames
cc.scores <- cc.scores[, c("S.Score", "G2M.Score", 
                           "Phase")]

data.filt_cellCycle3 <- CellCycleScoring(object = data.filt_cellCycle2,
                                         g2m.features = cc.genes$g2m.genes,
                                         s.features = cc.genes$s.genes)


VlnPlot(data.filt, features = c("S.Score","G2M.Score"), group.by= "orig.ident",ncol = 4, pt.size = .1)



## Predict doublets
suppressMessages(require(DoubletFinder))

data.filt = FindVariableFeatures(data.filt, verbose = F)
data.filt = ScaleData(data.filt, vars.to.regress = c("nFeature_RNA", "percent_mito"), verbose = F)
data.filt = RunPCA(data.filt, verbose = F, npcs = 20)
data.filt = RunUMAP(data.filt, dims = 1:10, verbose = F)
# Warning: The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
# To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
# This message will be shown once per session




# Can run parameter optimization with paramSweep, but skip for now.

sweep.res <- paramSweep_v3(data.filt)
sweep.stats <- summarizeSweep(sweep.res, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
barplot(bcmvn$BCmetric, names.arg = bcmvn$pK, las=2)

# define the expected number of doublet cellscells.
nExp <- round(ncol(data.filt)* 0.04) # expect 4% doublets
data.filt <- doubletFinder_v3(data.filt, pN=0.25, pK = 0.09, nExp = nExp, PCs = 1:10)


# name of the DF prediction can change, so extract the correct column name.
# DF.name = colnames(data.filt@meta.data)[grepl("DF.classification", colnames(data.filt@meta.data))]
DF.name = grep("DF.classification", colnames(data.filt@meta.data), value=TRUE)
table(data.filt$DF.classifications_0.25_0.09_230)
# Doublet Singlet 
#     230    5532 

cowplot::plot_grid( ncol = 2,
DimPlot(data.filt, group.by = "orig.ident") + NoAxes(),
DimPlot(data.filt, group.by = DF.name) + NoAxes()
)


VlnPlot(data.filt, features = "nFeature_RNA", group.by = DF.name, pt.size = .1)
VlnPlot(data.filt, features = "pANN_0.25_0.09_230", group.by = DF.name, pt.size = .1)


dim(data.filt)
# 18121  5762
data.filt = data.filt[,data.filt@meta.data[,DF.name] == "Singlet"]
dim(data.filt)
# 18121  5532


## Save data
dir.create('../data/results', showWarnings = F)

saveRDS(data.filt,"../data/results/seurat_covid_qc.rds")


#SESSION_INFO:
sessionInfo()



save.image("~/GitHub/workshop-scRNAseq/labs/seurat/Test_seurat_01_qc.R.RData")
