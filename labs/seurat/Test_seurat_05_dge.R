library(Seurat)
library(venn)
library(dplyr)
library(cowplot)
library(ggplot2)
library(pheatmap)
library(enrichR)
library(rafalib)
setwd("C:/Users/LYdia/Documents/GitHub/workshop-scRNAseq/labs/seurat")


alldata <- readRDS("../data/results/covid_qc_dr_int_cl.rds")



## Differential gene expression
#Set the identity as louvain with resolution 0.5
sel.clust = "CCA_snn_res.0.5"

alldata <- SetIdent(alldata, value = sel.clust)
table(alldata@active.ident)

# plot this clustering
plot_grid(ncol = 3,
  DimPlot(alldata, label = T) + NoAxes(),
  DimPlot(alldata, group.by = "orig.ident") + NoAxes(),
  DimPlot(alldata, group.by = "type") + NoAxes() )

plot_grid(ncol = 3,
          DimPlot(alldata, group.by = "CCA_snn_res.0.5", label = T) + NoAxes(),
          DimPlot(alldata, group.by = "orig.ident", label = T) + NoAxes(),
          DimPlot(alldata, group.by = "type", label = T) + NoAxes() )



# Cell marker genes
#Compute differential expression
markers_genes <- FindAllMarkers(alldata,
                               logfc.threshold = 0.2,
                               test.use = "wilcox",
                               min.pct = 0.1,
                               min.diff.pct = 0.2,
                               only.pos = TRUE,
                               max.cells.per.ident = 50,
                               assay = "RNA")
range(markers_genes$avg_log2FC)
range(pmax(markers_genes$pct.1, markers_genes$pct.2))
range(abs(apply(markers_genes[,3:4],1,diff)))

tapply(markers_genes$gene, markers_genes$cluster, length)
#    0    1    2    3    4    5    6    7    8    9 
# 1055  114  167  117  133  106  118   63  164  717
rle(as.character(markers_genes$cluster))
# Run Length Encoding
# lengths: int [1:10] 1055 114 167 117 133 106 118 63 164 717
# values : chr [1:10] "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"


markers_genes %>% group_by(cluster)  %>% top_n(-25, p_val_adj) -> top25
top25
rle(as.character(top25$cluster))
# Run Length Encoding
# lengths: int [1:10] 25 25 25 25 133 25 25 25 25 25
# values : chr [1:10] "0" "1" "2" "3" "4" "5" "6" "7" "8" "9"
markers_genes[markers_genes$cluster=="4",]
# the gene number in cluster 4 is more than 25 due to ties



mypar(2,5,mar=c(4,6,3,1))
for(i in unique(top25$cluster)){
  barplot( sort( setNames(top25$avg_log2FC, top25$gene) [top25$cluster == i])[1:25],
           horiz = T,las=1 ,main=paste0(i," vs. rest"),border = "white", yaxs="i" )
  abline(v=c(0,0.25),lty=c(1,2))
}



markers_genes %>% group_by(cluster)  %>% top_n(-5, p_val_adj) -> top5

# create a scale.data slot for the selected genes
alldata <- ScaleData(alldata, features = as.character(unique(top5$gene)), assay = "RNA")
DoHeatmap(alldata, features = as.character(unique(top5$gene)),group.by = sel.clust, assay = "RNA")

DotPlot(alldata, features = rev(as.character(unique(top5$gene))),group.by = sel.clust,assay = "RNA")+coord_flip()


# take top 3 genes per cluster/
top5 %>% group_by(cluster)  %>% top_n(-3, p_val) -> top3


# set pt.size to zero if you do not want all the points to hide the violin shapes, or to a small value like 0.1
VlnPlot(alldata, features = as.character(unique(top3$gene)), ncol = 5, group.by = sel.clust, assay = "RNA", pt.size = 0)



## Differential expression across conditions
# select all cells in cluster 1
cell_selection <- subset(alldata, cells = colnames(alldata)[ alldata@meta.data[,sel.clust] == 2])
cell_selection <- SetIdent(cell_selection, value = "type")
#Compute differential expression
DGE_cell_selection <- FindAllMarkers(cell_selection,
                               logfc.threshold = 0.2,
                               test.use = "wilcox",
                               min.pct = 0.1,
                               min.diff.pct = 0.2,
                               only.pos = TRUE,
                               max.cells.per.ident = 50,
                               assay = "RNA")

DGE_cell_selection %>% group_by(cluster)  %>% top_n(-5, p_val) -> top5_cell_selection

VlnPlot(cell_selection, features = as.character(unique(top5_cell_selection$gene)),
        ncol = 5,group.by = "type",assay = "RNA", pt.size = .1)


VlnPlot(alldata, features = as.character(unique(top5_cell_selection$gene)),
        ncol = 5, split.by = "type",assay = "RNA", pt.size = 0)



## Gene Set Analysis
# Hypergeometric enrichment test
# Load additional packages
library(enrichR)

# Check available databases to perform enrichment (then choose one)
enrichR::listEnrichrDbs()

# Perform enrichment
enrich_results <- enrichr(
 genes     = DGE_cell_selection$gene[DGE_cell_selection$cluster == "Covid"],
 databases = "GO_Biological_Process_2017b" )[[1]]

# Some databases of interest:
# GO_Biological_Process_2017b
# KEGG_2019_Human
# KEGG_2019_Mouse
# WikiPathways_2019_Human
# WikiPathways_2019_Mouse



par(mfrow=c(1,1),mar = c(3, 25, 2, 1))
barplot( height    = -log10(enrich_results$P.value)[10:1],
        names.arg = enrich_results$Term[10:1],
        horiz     = TRUE,
        las       = 1,
        border    = FALSE,
        cex.names = .6 )
abline(v = c(-log10(0.05)), lty = 2)
abline(v = 0, lty = 1)


# Gene Set Enrichment Analysis (GSEA)
DGE_cell_selection <- FindMarkers(cell_selection,
                                  ident.1 = "Covid",
                               logfc.threshold = -Inf,
                               test.use = "wilcox",
                               min.pct = 0.1,
                               min.diff.pct = 0,
                               only.pos = FALSE,
                               max.cells.per.ident = 50,
                               assay = "RNA")

# Create a gene rank based on the gene expression fold change
gene_rank <- setNames( DGE_cell_selection$avg_log2FC, casefold(rownames(DGE_cell_selection),upper=T) )


#Download gene sets
msigdbgmt <- msigdbr::msigdbr("Homo sapiens")
msigdbgmt <- as.data.frame(msigdbgmt)

#List available gene sets
unique(msigdbgmt$gs_subcat)
# "MIR:MIR_Legacy"  "TFT:TFT_Legacy"  "CGP"             "TFT:GTRD"        ""               
# "VAX"             "CP:BIOCARTA"     "CGN"             "GO:BP"           "GO:CC"          
# "IMMUNESIGDB"     "GO:MF"           "HPO"             "CP:KEGG"         "MIR:MIRDB"      
# "CM"              "CP"              "CP:PID"          "CP:REACTOME"     "CP:WIKIPATHWAYS"

#Subset which gene set you want to use.
msigdbgmt_subset <- msigdbgmt[msigdbgmt$gs_subcat == "CP:WIKIPATHWAYS",]
gmt <- lapply( unique(msigdbgmt_subset$gs_name),function(x){msigdbgmt_subset [msigdbgmt_subset$gs_name == x ,"gene_symbol"]} )
names(gmt) <- unique(paste0(msigdbgmt_subset$gs_name,"_",msigdbgmt_subset$gs_exact_source))



library(fgsea)

# Perform enrichemnt analysis
fgseaRes <- fgsea( pathways=gmt, stats=gene_rank, minSize=15, maxSize=500,nperm = 10000)
fgseaRes <- fgseaRes[ order(fgseaRes$ES, decreasing = T) ,]

# Filter the results table to show only the top 10 UP or DOWN regulated processes (optional)
top10_UP <- fgseaRes$pathway [1:10]

# Nice summary table (shown as a plot)
dev.off()
plotGseaTable(gmt[top10_UP], gene_rank, fgseaRes, gseaParam = 0.5)



saveRDS(alldata,"../data/3pbmc_qc_dr_int_cl_dge.rds")
write.csv(markers_genes,"../data/3pbmc_qc_dr_int_cl_dge_markerGenes.csv")




sessionInfo()

save.image("~/GitHub/workshop-scRNAseq/labs/seurat/Test_seurat_05_dge.RData")
