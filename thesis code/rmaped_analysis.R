library(Seurat)
library(ggplot2)
library(dplyr)
setwd("H:/R/R 09Apr2020/Data/")
########################################
## Load in the RNA UMI matrix
JL1.rna <- as.sparse(read.csv(file = "E:/scRNA data Remap/Data/180426_Jane_1_S1_ExpressionMatrix.csv", sep = ",", 
                               header = TRUE, row.names = 1))
JL2.rna <- as.sparse(read.csv(file = "E:/scRNA data Remap/Data/180426_Jane_2_S1_ExpressionMatrix.csv", sep = ",", 
                              header = TRUE, row.names = 1))
JL3.rna <- as.sparse(read.csv(file = "E:/scRNA data Remap/Data/180426_Jane_3_S2_ExpressionMatrix.csv", sep = ",", 
                              header = TRUE, row.names = 1))
JL4.rna <- as.sparse(read.csv(file = "E:/scRNA data Remap/Data/180426_Jane_4_S2_ExpressionMatrix.csv", sep = ",", 
                              header = TRUE, row.names = 1))
STATIC_DMSO <- CreateSeuratObject(counts = JL1.rna)
STATIC_SBCHIR <- CreateSeuratObject(counts = JL2.rna)
SHEAR_DMSO <- CreateSeuratObject(counts = JL3.rna)
SHEAR_SBCHIR <- CreateSeuratObject(counts = JL4.rna)
########################################
# write experimental group information into cell id
#head(x = colnames(x = STATIC_DMSO))
#STATIC_DMSO<-RenameCells(STATIC_DMSO,add.cell.id = "C1")
#head(x = colnames(x = STATIC_DMSO))

#head(x = colnames(x = STATIC_SBCHIR))
#STATIC_SBCHIR<-RenameCells(STATIC_SBCHIR,add.cell.id = "C2")
#head(x = colnames(x = STATIC_SBCHIR))

#head(x = colnames(x = SHEAR_DMSO))
#SHEAR_DMSO<-RenameCells(SHEAR_DMSO,add.cell.id = "C3")
#head(x = colnames(x = SHEAR_DMSO))

#head(x = colnames(x = SHEAR_SBCHIR))
#SHEAR_SBCHIR<-RenameCells(SHEAR_SBCHIR,add.cell.id = "C4")
#head(x = colnames(x = SHEAR_SBCHIR))
########################################
## add metadata with the experimental group for each cell
GroupID<-rep("STATIC_DMSO",length(colnames(STATIC_DMSO)))
Sample<-data.frame(GroupID,row.names=colnames(STATIC_DMSO))
STATIC_DMSO<-AddMetaData(STATIC_DMSO,Sample,col.name="GroupID")
head(STATIC_DMSO@meta.data)

GroupID<-rep("STATIC_SBCHIR",length(colnames(STATIC_SBCHIR)))
Sample<-data.frame(GroupID,row.names=colnames(STATIC_SBCHIR))
STATIC_SBCHIR<-AddMetaData(STATIC_SBCHIR,Sample,col.name="GroupID")
head(STATIC_SBCHIR@meta.data)

GroupID<-rep("SHEAR_DMSO",length(colnames(SHEAR_DMSO)))
Sample<-data.frame(GroupID,row.names=colnames(SHEAR_DMSO))
SHEAR_DMSO<-AddMetaData(SHEAR_DMSO,Sample,col.name="GroupID")
head(SHEAR_DMSO@meta.data)

GroupID<-rep("SHEAR_SBCHIR",length(colnames(SHEAR_SBCHIR)))
Sample<-data.frame(GroupID,row.names=colnames(SHEAR_SBCHIR))
SHEAR_SBCHIR<-AddMetaData(SHEAR_SBCHIR,Sample,col.name="GroupID")
head(SHEAR_SBCHIR@meta.data)

aggr <- merge(STATIC_DMSO,y = c(STATIC_SBCHIR,SHEAR_DMSO,SHEAR_SBCHIR),add.cell.ids = c("C1","C2","C3","C4"),project = "all")
saveRDS(aggr,file="remaped_data/RemapedAggregate.Robj")
########################################
## remove
########################################
## General data process pipeline
load("remaped_data/ReMapAggregate.Robj")
## General data visualisation by vlnPlot & GenePlot
# The number of genes and UMIs (nGene and nUMI) are automatically calculated for every object by Seurat.  For non-UMI
# data, nUMI represents the sum of the non-normalized values within a cell We calculate the percentage of mitochondrial
# genes here and store it in percent.mito using AddMetaData.  We use aggr@counts (the raw data) since this represents
# non-transformed and non-log-normalized counts The % of UMI mapping to MT-genes is a common scRNA-seq QC metric.
mito.genes <- grep(pattern = "^MT-", x = rownames(aggr), value = TRUE)
percent.mito <- Matrix::colSums(aggr[mito.genes, ])/Matrix::colSums(aggr)
# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats.  This also allows us to plot the
# metadata values using the Seurat's VlnPlot().
head(aggr@meta.data)  # Before adding
aggr <- AddMetaData(object = aggr, metadata = percent.mito, col.name = "percent.mito")
head(aggr@meta.data) # After adding
save(aggr,file="remaped_data/ReMapAggregate.Robj")
VlnPlot(object = aggr, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"),group.by = "GroupID")
# nCount_RNA is equavalent to nUMI; nFeature_RNA is equivalent to nGene
FeatureScatter (object = aggr, "nCount_RNA",  "percent.mito",group.by = "GroupID")
FeatureScatter (object = aggr, "nCount_RNA", "nFeature_RNA",group.by = "GroupID")
########################################
## Filter out dying cells with high percentage of mitochondrial gene and doublets with outside range of No.Feature genes(No.Gene)
# We filter out cells that have unique gene counts over 5000 or less than
# 200 Note that low.thresholds and high.thresholds are used to define a
# 'gate'.  -Inf and Inf should be used if you don't want a lower or upper
# threshold.
aggr <- subset(x = aggr, subset = nFeature_RNA > 50 & nFeature_RNA < 5000 & percent.mito >  -Inf & percent.mito < 0.05 )
aggr <- NormalizeData(object = aggr, normalization.method = "LogNormalize", 
                         scale.factor = 10000)
########################################
## Detection of variable genes across the single cells
# In S3 FindVariableGenes has been changed to FindVariableFeature, with a optimised n.feature=2000 in data scaled to 10E4
# and the cutoff are not used.
# New default method for FindVariableFeatures is VST. Seurat 2 use MVP
# vst: First, fits a line to the relationship of log(variance) and log(mean) using local polynomial regression (loess). 
#      Then standardizes the feature values using the observed mean and expected variance (given by the fitted line). 
#      Feature variance is then calculated on the standardized values after clipping to a maximum (see clip.max parameter).
# mean.var.plot (mvp): First, uses a function to calculate average expression (mean.function) and dispersion (dispersion.function) for each feature.
#      Next, divides features into num.bin (deafult 20) bins based on their average expression, and calculates z-scores for dispersion within each bin.
#      The purpose of this is to identify variable features while controlling for the strong relationship between variability and average expression.

#aggr <- FindVariableFeatures(object = aggrNew,Selection.method="mean.var.plot", mean.function = ExpMean, dispersion.function = LogVMR, 
#                             x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)

aggr <- FindVariableFeatures(object = aggr,selection.method = "vst",nfeatures = 5000) # nfeatures found from S2 work was 4612 genes

# calculate how many varibal genes identified by cutoff filtering
head(x = HVFInfo(object = aggr))

# scaling data by repression
aggr <- ScaleData(object = aggr, vars.to.regress = c("nCount_RNA", "percent.mito"))

VlnPlot(object = aggr, features = c("mEGFP"),group.by = "GroupID")
VlnPlot(object = aggr, features = c("SOX17"),group.by = "GroupID")
# perfom PCA
aggr<- RunPCA(object = aggr, features = VariableFeatures(object = aggr), ndims.print = 1:10,
                  nfeatures.print = 30)

DimPlot(object = aggr, reduction = "pca", group.by = "GroupID")
VariableFeaturePlot(object = aggr)
aggr <- ProjectDim(object = aggr)
DimHeatmap(object = aggr, reduction = "pca", cells = 500, balanced = TRUE)
PCHeatmap(object = aggr, reduction = "pca", dims = 1:12, cells = 500, balanced = TRUE)
ElbowPlot(aggr,ndims = 50)
aggr <- JackStraw(object = aggr, num.replicate = 100,reduction = "pca")
aggr <- ScoreJackStraw(aggr,dims = 1:20)
JackStrawPlot(object = aggr, dim = 1:20, reduction = "pca")
########################################
## Cluster the cells
# save.SNN = T saves the SNN so that the clustering algorithm can be rerun
# using the same graph but with a different resolution value 

aggr <- FindNeighbors(aggr,reduction = "pca", dims = 1:20)
aggr <- FindClusters(aggr, reduction = "pca", dims = 1:20, resolution = 0.6)

## Old color set
#Col.set <- c("#fb8500", "#ff0000", "#008e17", "#fde800", "#4ffc00", "#bc9000", "#0099cc", "#00bcac",
#             "#00eefd", "#5f777f", "#cf6bd6", "#D35400", "#99cc00","#aa00ff", "#ff00ff", "#0053c8",
#             "#f2a287", "#800000")
########################################
## Run and plot TSNE
aggr <- RunTSNE(object = aggr, dims.use = 1:20, do.fast = TRUE, do.label=TRUE)
DimPlot(object = aggr, reduction = "tsne", label = T,label.size = 4,pt.size = 2, split.by="GroupID")
########################################
## Run UMAP
aggr <- RunUMAP(aggr, reduction = "pca", dims = 1:20)
DimPlot(aggr, reduction = "umap",label = T, label.size = 4,pt.size = 2, split.by = "GroupID")
DimPlot(aggr, reduction = "umap",label = T, label.size = 10, pt.size = 2)

##
# Quick inventery of interest genes expression
FeaturePlot(aggr,c("SOX17","mEGFP","RUNX1","CD34"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("SOX17"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("mEGFP"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("CD34"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("PTPRC"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("SPN"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("CXCR4"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("HOXA10"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("CD34"),label = T)

FeaturePlot(aggr,c("NCKAP1L"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("FN1"),label = T, split.by = "GroupID")
FeaturePlot(aggr,c("COL1A2"),label = T, split.by = "GroupID")

VlnPlot(object = aggr, features = c("HOXA11","HOXA10","HOXA9"),group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("DNM3OS"),group.by = "ClusterID_res0.6",,split.by = "GroupID")

FeaturePlot(aggr,c("CPE"),label = T,split.by = "GroupID")
FeaturePlot(aggr,c("DNM3OS","GPC3"),label = T,split.by = "GroupID")
########################################
## Save current cluster identites in object@meta.data under 'clusterID'
# Run only if Seurat::FindClusters() was executed
aggr[["ClusterID_res0.6"]] <- Idents(object = aggr)
# create new ID combine clusters and conditions
MixID_0.6<- as.character(paste(aggr@meta.data$ClusterID_res0.6, aggr@meta.data$GroupID, sep = "_"))
# Save new ID to seurat object @Meta.data and @ident
Sample<-data.frame(MixID_0.6,row.names=colnames(aggr))
aggr<-AddMetaData(aggr,Sample,col.name="MixID_0.6")
head(aggr@meta.data)
########################################
## Annotation
#cluster 0,5, 6,8,1 Mesenchymal stromal cell?Fibroblast? Vascular smooth muscle cells? pericytes? 
FeaturePlot(aggr,c("PDGFRA","HOXA10"),label = T)
########################################
#cluster 0 mesodermal pericytes HOXA+PDGFRa+ (Old X0)
FeaturePlot(aggr,c("DNM3OS","IGFBP5","LUM","MEST"),label = T)
FeaturePlot(aggr,c("GPC3","MME", "POSTN","IGFBP3"),label = T)
FeaturePlot(aggr,c("MDK","PTN", "WNT5A"),label = T)
FeaturePlot(aggr,c("PDGFRA","HOXA10","THY1"),label = T)
FeaturePlot(aggr,c("GFAP","VIM","ACTA2","PALLD"),label = T)
# "VIM" vimintin a genral mesenchymal marker
# "ACTA2" Smooth muscle actin A
# "PALLD" cytoskeletal actin
########################################
#cluster 1 Mesenchymal stromal cell HOXA+PDGFRa+ highly proliferative (Old X6)
FeaturePlot(aggr,c("TOP2A","MKI67","CENPF","HMGB2"),label = T)
FeaturePlot(aggr,c("NUSAP1","HIST1H4C","TPX2","PRC1"),label = T)
FeaturePlot(aggr,c("DLGAP5","ASPM","SMC4","PTTG1"),label = T)
FeaturePlot(aggr,c("MCAM","PDGFRB","CSPG4","DLK1"),label = T)
FeaturePlot(aggr,c("APLNR","PDGFRA","KDR","PECAM1"),label = T)
# DLK1+ arteriolar PCs https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6428685/ 
########################################
#cluster 3  Myofibroblast C7+ HOXA-PDGFRa+ (Old X1)
FeaturePlot(aggr,c("PDGFRA","CALD1"),label = T)
FeaturePlot(aggr,c("FBN1","LAMB1","MYH10","FN1"),label = T)
FeaturePlot(aggr,c("COL1A1","VCAN","BGN","CGNL1"),label = T)
FeaturePlot(aggr,c("MGP","C7","HAPLN1","CDH11"),label = T)
FeaturePlot(aggr,c("BEX1","GUCY1A3","LRRN3","ITGA1"),label = T)

########################################
#cluster 16 Cardiomyocyte (Old X2) 
FeaturePlot(aggr,c("TTN","MT-RNR2","MAP1B","VCAN"),label = T)
FeaturePlot(aggr,c("MYL7","TNNT2","MYH7","TTN"),label = T)
########################################
#cluster 6 Somitic/Vascular Smooth Muscle HOXA+  (Old X3) 
FeaturePlot(aggr,c("LUM","IGF2","TAGLN","MYH10"),label = T)
FeaturePlot(aggr,c("TGFBI","MME","PTN","PDGFRA"),label = T)
FeaturePlot(aggr,c("HOXA10","VIM","ACTA2","DNM3OS"),label = T)
FeaturePlot(aggr,c("CAV1", "KDR"),label = T)
# subaortic mesenchymal population that expresses smooth muscle cell markers 
# and that contributes to blood formation via an endothelial intermediate
#TBX18 gene : this gene may allow some progenitor cell develop into 
#smooth muscel progenitor. https://doi.org/10.1113/JP270033
# DNM30S  https://doi.org/10.1002/dvdy.21787            
########################################
#cluster 5 Somitic/Vascular Smooth Muscle HOXA- (Old X1)
FeaturePlot(aggr,c("IGFBP7","CTHRC1","SDC2","ACTA2"),label = T)
FeaturePlot(aggr,c("HOXA10","VIM","ACTA2","DNM3OS"),label = T)
FeaturePlot(aggr,c("MYH10","MYH11","IGF2","PTN"),label = T)
########################################
#Cluster 4  Megakaryocyte (old X5)
FeaturePlot(aggr,c("PPBP","PF4","SPN","PLEK"),label = T)
FeaturePlot(aggr,c("CLEC1B","TUBB1","THBS1","RGX18"),label = T)
FeaturePlot(aggr,c("RAB27B","GPX1","GNG11","ITGA2B"),label = T)
FeaturePlot(aggr,c("TLN1","ACTB","LIMS1","VCL"),label = T)
# ITGA2B CD41 Platelet marker
########################################

#Cluster 2  MEP (old X7)
FeaturePlot(aggr,c("HBG1","HBG2","CNST","CMIP"),label = T)
FeaturePlot(aggr,c("PPBP","PF4","TXNIP","MTURN"),label = T)
FeaturePlot(aggr,c("BTK","NEAT1","ARHGAP6","MPP1"),label = T)
########################################
#Cluster 10 Macrophages CD14 (old X8)
FeaturePlot(aggr,c("SPP1","CD74","C1QB","CTSB"),label = T)
FeaturePlot(aggr,c("HLA-DRA","AIF1","RAG1","KIT"),label = T)
FeaturePlot(aggr,c("FCGR3A","CD14","C1QA","CYBB"),label = T)
FeaturePlot(aggr,c("CSF1R"),label = T)
# "FCGR3A" CD16 https://www.genecards.org/cgi-bin/carddisp.pl?gene=FCGR3A
########################################
#cluster 13 Aortic endothelium / Haematopoietic progenitor cells(Old X9)
FeaturePlot(aggr,c("CD93","ERG","CD34","KDR"),label = T)
FeaturePlot(aggr,c("SOX17","GJA5","IGFBP5","CDH5"),label = T)
FeaturePlot(aggr,c("SOX17","NOTCH1","DLL4","CDH5"),label = T)
FeaturePlot(aggr,c("CAV1","MLLT3","IGFBP5","IGFBP4"),label = T)
FeaturePlot(aggr,c("NEURL3"),label = T)
FeaturePlot(aggr,c("IGFBP5","IGFBP4","KDR"),label = T)
VlnPlot(object = aggr, features = c("CD34","SOX17","NOTCH1","DLL4","CDH5"),group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("PTGER2"),group.by = "ClusterID_res0.6",split.by = "GroupID")


VlnPlot(object = aggr, features = c("IGFBP5"),group.by = "ClusterID_res0.6",split.by = "GroupID")
VlnPlot(object = aggr, features = c("IGFBP4"),group.by = "ClusterID_res0.6",split.by = "GroupID")
VlnPlot(object = aggr, features = c("KDR"),group.by = "ClusterID_res0.6",split.by = "GroupID")
VlnPlot(object = aggr, features = c("IGFBP5","IGFBP4","KDR"),group.by = "ClusterID_res0.6")
#CD144 CDH5
###########################################
#Cluster 14 Lymphatic Endothlium 
# HAPLN1 LEC https://pubmed.ncbi.nlm.nih.gov/18084611/
FeaturePlot(aggr,c("CD93","ERG","CD34","KDR"),label = T)
FeaturePlot(aggr,c("SOX17","GJA5","COL4A1","CDH5"),label = T)
FeaturePlot(aggr,c("HAPLN1","CDH11","CD34","KDR"),label = T)
VlnPlot(object = aggr, features = c("ITGB2"),group.by = "ClusterID_res0.6",split.by = "GroupID")


#Cluster 6 has some CD73 expression as well as Cluster 13, Cluster14 
# KDR is high in Cluster 2, 13, 14
########################################
#cluster 16  Cardiomyocytes (Old X10)
FeaturePlot(aggr,c("TTN","NEBL","ACTC1","TNNT2"),label = T)
FeaturePlot(aggr,c("MYL7","EMC10","ALPK2","MYH7"),label = T)
FeaturePlot(aggr,c("CDH2","PCDH7","SERINC5","SORBS2"),label = T)
########################################
#cluster 17 myeloid Progenitors  
FeaturePlot(aggr,c("SRGN","CD24","RNASE2","BPI"),label = T)
FeaturePlot(aggr,c("PRG2","PRG3","EPX","LYZ"),label = T)
########################################
#cluster 18  myeloid Progenitors-(MPO) (Old X12)
FeaturePlot(aggr,c("ITGB2","MPO"),label = T)
FeaturePlot(aggr,c("SRGN","CD24","RNASE2","BPI"),label = T)
########################################
#cluster 9   Common myeloid progenitors (CPA3, HDC) (Old X13)
#https://doi.org/10.1182/blood.V126.23.2768.2768 
FeaturePlot(aggr,c("CPA3","HDC","GATA2","CLC"),label = T)
FeaturePlot(aggr,c("SAMSN1","S100A4","LCP1","SRGN"),label = T)
FeaturePlot(aggr,c("PRG2","PRG3","EPX","LYZ"),label = T)
########################################
#cluster 7  MPPs (CD34 & CD43, TMPO,CD3E, ) (Old X14)
FeaturePlot(aggr,c("AIF1","HMGA1","DUT","S100A4"),label = T)
FeaturePlot(aggr,c("HMGB2","VAMP8","H2AFZ","MKI67"),label = T)
FeaturePlot(aggr,c("HSPE1","ABRACL","RAN","C1QBP"),label = T)
FeaturePlot(aggr,c("ANP32B","CENPF","CACYBP","SNRPG"),label = T)
FeaturePlot(aggr,c("CD34","PTPRC","mEGFP","SPN"),label = T)
FeaturePlot(aggr,c("TMPO"),label = T)
FeaturePlot(aggr,c("CD3E"),label = T)
# from Ali's paper, TMPO is main marker for one of the haematopoietic progenitor marker
# it is expressed in Cluster 7, 17, 18, 9 and a lot in 1.
########################################
#cluster 12  Erythroid  (Old X15)
FeaturePlot(aggr,c("GYPA","HBZ","HBA2","HBE1"),label = T)
FeaturePlot(aggr,c("GYPB","HBM","ALAS2","HBG1"),label = T)
FeaturePlot(aggr,c("HBA1","HBZ"),label = T,split.by = "GroupID")
#######################################
#cluster 15   lymphocytes or NK progenitor  (Old X16)
FeaturePlot(aggr,c("KLF10","CREBRF","CLK1","TXNIP"),label = T)
FeaturePlot(aggr,c("NFKBIA","TANK","HES1","DDIT3"),label = T)
FeaturePlot(aggr,c("CD69","TRA2A","KIFAP3","EPB41"),label = T)
FeaturePlot(aggr,c("CD69","NFKBIA","KLRB1","FCGR3A"),label = T)
#CD69:Involved in lymphocyte proliferation and 
# functions as a signal transmitting receptor in lymphocytes, 
# natural killer (NK) cells, and platelets 
# https://www.genecards.org/cgi-bin/carddisp.pl?gene=CD69
#########################################
#cluster 11  sympathoadrenal progenitors 
# doi: 10.1097/PAS.0b013e3182a0218f   https://cordis.europa.eu/project/id/302739/reporting
FeaturePlot(aggr,c("ERBB4","MSX2","TBR1","PRTG","GATA3"),label = T)
FeaturePlot(aggr,c("GATA3","PDGFRA","HAND1","TH"),label = T)
#GATA3 is very high in Cluster 11, Ali's paper use it also for NK and T cell progenitor.
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3510442/ Simon's paper show GATA3 not express in haemaotpoietic cells
# in E11.5 AGMs but paraortic cells, but in mesonephric ducts in the developing kidney and
# sympathoadrenal progenitors (Coexpressed HAND1 and TH) which form the developing adrenal gland and sympathetic ganglia of the SNS.
# GATA3 deficiency decrease HSCs production independent from blood flow but via TH and hence catecholamine production in cells of the SNS.
########################################



FeaturePlot(aggr,c("RIPOR3","TSC22D1","MTURN","RAP2B"),label = T)
FeaturePlot(aggr,c("MALAT1"),label = T)
########################################
# for (gene in c10$gene) {
#   FeaturePlot(abcd,gene,label = T)
#   ggsave(filename = paste("c10/",gene,".png",sep = ""))
# }
######################################## 

########################################
#cluster 14 renal stromal progenitors
#doi: 10.1097/PAS.0b013e3182a0218f   https://cordis.europa.eu/project/id/302739/reporting
FeaturePlot(aggr,c("ERBB4","MSX2","TBR1","PRTG","GATA3"),label = T)
#maybe kidney
#abcd_anotate <- RenameIdents(abcd_anotate,"14" = "14.renal stromal progenitors")
#DimPlot(abcd_anotate,label = T)



########################################
#intra paraxial mesoidem
#cluster 0 Somitic https://doi.org/10.1002/dvdy.21787           AGM stroma cell
#TBX18 gene : this gene may allow some progenitor cell develop into smooth muscel progenitor. https://doi.org/10.1113/JP270033
FeaturePlot(aggr,c("TBX18","ACTA2"),label = T)
FeaturePlot(aggr,c("DNM3OS","LUM","MEST","HOXA11"),label = T)
abcd <- RenameIdents(abcd,"0" = "Skeletal system?")
abcd <- RenameIdents(abcd,"Skeletal system?" = "Somite")
DimPlot(abcd,label = T)

########################################
#cluster 4 
FeaturePlot(aggr,c("TGFBI","COL12A1","ACTA2","COL11A1"),label = T)
#TGFBI: tranforming growth factor, May resultin endochondrial ossification.
#ACTA2 gene : aortic smooth muscle gene, mutation cause variety of cascular diseases.
FeaturePlot(aggr,c("LUM","IGFBP5","PARM1","MEST"),label = T)
########################################
## cluster annotation Plotting
Idents(object = aggr)<-aggr[["ClusterID_res0.6"]]
aggr <- RenameIdents(aggr,"0"="Stromal HOXA+PDGFRA+","1"="Stromal HOXA-PDGFRA+","2"="Mesenchymal PDGFRA+","3"="MEMP-1","4"="Megakaryocyte","5"="Stromal",
                          "6"="Endothelium","7"="Myeloid-2","8"="Macrophages","9"="MPP-1","10"="Cardiac Smooth muscle cells","11"="Erythroid",
                           "12"="Cardiomyocytes","13"="Myeloid-3","14"="Macrophage Progenitor","15"="MPP-2","16"="Myeloid-1")
Col.set <- c("#fb8500", "#ff0000", "#008e17", "#fde800", "#4ffc00", "#bc9000", "#0099cc", "#00bcac",
             "#00eefd", "#5f777f", "#cf6bd6", "#D35400", "#99cc00","#aa00ff", "#ff00ff", "#0053c8",
              "#f2a287", "#800000","#001100")
  
#  c("#fb8500", "#ff0000", "#0099cc", "#00bcac", "#bc9000", "#fde800", "#5f777f", "#aa00ff",
#             "#00eefd", "#ff00ff", "#008e17", "#008e17", "#0053c8","#cf6bd6", "#f2a287", "#ff00ff",
#             "#99cc00")
# saveRDS(aggr,file="remaped_data/RemapedAggregate_clustered.Robj")
DimPlot(aggr, reduction = "umap",cols=Col.set,label = T, label.size = 7, pt.size = 1.5)
DimPlot(aggr, reduction = "umap",cols=Col.set,label = T, label.size = 5, pt.size = 1.5,split.by = "GroupID")
##############################################################################################
setwd("D:/Work/R/R 09Apr2020/Data/")
aggr <- readRDS("remaped_data/RemapedAggregate_clustered.Robj")
##############################
####### DPA analysis##########
##############################
setwd("D:/Work/R/R 09Apr2020/Data/DPA/")
source("diffprop_functions.R");
source("CountCellNum.R")
# use VlnPlot function to extract feature expressing cell counts
# library(dplyr)
# p <- VlnPlot(object = aggr, features =c("mEGFP"))
# GFP_counts_Group <- p$data %>% group_by(aggr@meta.data$GroupID) %>% summarize(counts = sum(mEGFP, na.rm = TRUE))
# GFP_counts_Clusters <- p$data %>% group_by(ident) %>% summarize(counts = sum(mEGFP, na.rm = TRUE))
# GFP_counts_MixID <- p$data %>% group_by(aggr@meta.data$MixID) %>% summarize(counts = sum(mEGFP, na.rm = TRUE))
setwd("H:/R/R 09Apr2020/Data/")
aggr <- readRDS("remaped_data/RemapedAggregate_clustered.Robj")

# load("updated_data/updatedAggregate.Robj")
# use function CountCellNum to extract feature expressing cell counts
genes <- c("mEGFP")
GFP_counts_Clusters <-PrctCellExpringGene(aggr,genes,group.by = "ident")
GFP_counts_Group <- PrctCellExpringGene(aggr,genes,group.by = "GroupID")
GFP_Counts_MixID <- PrctCellExpringGene(aggr,genes,group.by = "MixID")
write.csv(GFP_counts_Clusters,("GFP_counts_cluster_res0.6.csv"), quote= FALSE)
write.csv(GFP_counts_Group,("GFP_counts_Group_res0.6.csv"), quote= FALSE)
write.csv(GFP_Counts_MixID,("GFP_counts_MixID.csv"), quote= FALSE)

#use housekeeping genes for all cells. but may depends on cell types
genes <- c("MALAT1")
Cell_counts_Clusters <-PrctCellExpringGene(aggr,genes,group.by = "ident")
write.csv(Cell_counts_Clusters,("Cell_counts_cluster_res0.6.csv"), quote= FALSE)

## count cell number in each cluster
library(data.table)
library(magrittr)
# extract meta data
md <- aggr@meta.data %>% as.data.table
# the resulting md object has one "row" per cell
## count the number of cells per unique combinations of "Sample" and "seurat_clusters"
clustercount<-md[, .N, by = c("GroupID", "ClusterID_res0.6")] %>% dcast(., GroupID ~ ClusterID_res0.6, value.var = "N")

###############################
## Hiararchycal clustering####
###############################
# Find differentially expessed genes as biomarkers
cluster17.markers <- FindMarkers(aggr, ident.1 = 17, min.pct = 0.25)
head(cluster17.markers, n = 5)
cluster0.markers <- FindMarkers(aggr, ident.1 = 0, min.pct = 0.25)
head(cluster0.markers, n = 5)
cluster1.markers <- FindMarkers(aggr, ident.1 = 1, min.pct = 0.25)
cluster2.markers <- FindMarkers(aggr, ident.1 = 2, min.pct = 0.25)
cluster3.markers <- FindMarkers(aggr, ident.1 = 3, min.pct = 0.25)
Endo.markers <- FindMarkers(aggr, ident.1 = 13,ident.2 =14,  min.pct = 0.25)
cluster1.markers <- FindMarkers(aggr, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
fibroblast.markers <- FindMarkers(aggr, ident.1 = 0,ident.2 =8,  min.pct = 0.25)
# find all markers
aggr.markers <- FindAllMarkers(aggr, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# DoHeatmap(aggr, features = aggr.markers$gene,group.colors = Col.set,size = 3) + NoLegend()
# find markers for each cluster and test set as ROC
top10 <- aggr.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(aggr, features = top10$gene,group.colors = Col.set,size = 3) + NoLegend()
Top10DEG<-as.data.table(top10)
top50 <- aggr.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(aggr, features = top50$gene,group.colors = Col.set,size = 3) + NoLegend()
# write differentially expressed genes 

write.csv(aggr.markers, ("Differnetially expressed genes.csv"), quote= FALSE)
write.csv(top10, ("TOP10 Differnetially expressed genes.csv"), quote= FALSE)
write.csv(Endo.markers, ("Endothelial cluster Differnetially expressed genes.csv"), quote= FALSE)
write.csv(fibroblast.markers, ("fibroblast cluster Differnetially expressed genes.csv"), quote= FALSE)
################################################
## Find all DEGs between all clusters 23Jun2020
#################################################
library(plyr)
library(tidygraph)
library(igraph)
library(ggraph)
library(scales)
library(Seurat)
library(progress)
library(gplots)
setwd("H:/R/R 09Apr2020/Data/Plotting and Expression functions")
source("plotting_functions_JL.R")
# source("plotting_functions.R")
source("expression_functions.R")

writeAllDifferentialExpressionTests(aggr,"H:/R/R 09Apr2020/Data/Plotting and Expression functions/Differentially expressed genes_Remapped/")


#############################################
# Create a rule based model to identify clusters using gene markers using C5.0
#############################################
library(C50)
setwd("H:/R//R 09Apr2020/Data/remaped_data/")
source("plotting_functions_JL.R")
# combine all differentially expressed genes
Gene.Names<-rownames(aggr.markers)
expression.data<-t(aggr@assays$RNA@data)
ndx<-colnames(expression.data) %in% Gene.Names 
Predictor<-as.matrix(expression.data[,ndx])
cell.clusters<-aggr@active.ident
DT<-C5.0(Predictor,cell.clusters)
genes.importance<-C5imp(DT) # get importance of genes used for classification

# select 300 most important genes
Gene.Names.Selection<-rownames(genes.importance)[1:300]
ndx<-colnames(expression.data) %in% Gene.Names.Selection 
Predictor<-as.matrix(expression.data[,ndx])
DT<-C5.0(Predictor,cell.clusters)
summary(DT)
doHC2DotPlot(aggr,Gene.Names.Selection)

FeaturePlot(aggr,c("NCKAP1L","FN1","COL1A2","COL3A1"),label = F)

write.csv(cluster.heatmap, ("HC Heatmap of 300 genes.csv"), quote= FALSE)

####################################################
## Calculating total expression of HOXA genes
setwd("H:/R/R 09Apr2020/Data/DPA/")
# source("diffprop_functions.R");
source("CountCellNum.R")
Idents(object = aggr) <- 'GroupID'
Group.averages<-AverageExpression(aggr,features = "HOXA11")
Group.averages$HOXA13<-AverageExpression(aggr,features = "HOXA13")
Group.averages$HOXA11<-AverageExpression(aggr,features = "HOXA11")
Group.averages$HOXA10<-AverageExpression(aggr,features = "HOXA10")
Group.averages$HOXA9<-AverageExpression(aggr,features = "HOXA9")
Group.averages$HOXA7<-AverageExpression(aggr,features = "HOXA7")
Group.averages$HOXA5<-AverageExpression(aggr,features = "HOXA5")
Group.averages$HOXA3<-AverageExpression(aggr,features = "HOXA4")
Group.averages$HOXA3<-AverageExpression(aggr,features = "HOXA3")
Group.averages$HOXA2<-AverageExpression(aggr,features = "HOXA2")
Group.averages$HOXA1<-AverageExpression(aggr,features = "HOXA1")
HOXA13_counts_Group <- PrctCellExpringGene(aggr,"HOXA13",group.by = "GroupID")

library(data.table)
HOXAgenes<-c("HOXA1","HOXA2","HOXA3","HOXA4","HOXA5","HOXA6","HOXA7","HOXA9","HOXA10","HOXA11","HOXA13")
genes<-c("HOXA11")
#i=1;
#for(genes in HOXAgenes){
 # a<- VlnPlot(object = aggr, features = genes,group.by = "GroupID")
#  aa<-as.data.table(a$data)
#  Group.Total<-xtabs(genes~ident, aa)
#  i<-i+1
#}
remove(a)
remove(aa)
#
a<-VlnPlot(object = aggr, features = c("HOXA1"),group.by = "GroupID")
aa<-as.data.table(a$data)
xtabs(HOXA1~ident, aa)


############################################################
## GO Pathway over-representation test with GSOAP package ##
############################################################
library()
devtools::install_github(c("GuangchuangYu/DOSE", "GuangchuangYu/clusterProfiler"))


library(gsoap)

####################################################################################
## Find differentially expessed genes between experimental condition##
####################################################################################
Data<-aggr
Idents(object = Data) <- 'GroupID'
# find all markers for different experimental condition
Conditions.markers <- FindAllMarkers(Data, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Conditions.markers, ("DEGs between condtions.csv"), quote= FALSE)
top10.conditions <- Conditions.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Data, features = top10.conditions$gene,size = 5) + NoLegend()
write.csv(top10.conditions, ("Top10 DEGs between condtions.csv"), quote= FALSE)




###################################################################################
## HOXA genes plot######################################
######################################################
VlnPlot(object = aggr, features = c("HOXA9","HOXA10","HOXA11"),group.by = "GroupID")
VlnPlot(object = aggr, features = HOXAgenes,cols=Col.set,group.by = "ClusterID_res0.6")
col = rgb(red = 1, green = 0, blue = 0, alpha = 0.5)
FeaturePlot(aggr,features = c("HOXA9","HOXA10","HOXA11"),pt.size=1.2,label = T)

###Need to custom ggplot for color transparency
# p <- ggplot(a@data$HOXA9, axes(UMAP1, UMAP2))
# p + geom_point()
VlnPlot(object = aggr, features = c("PDGFRA"),cols=Col.set,group.by = "GroupID")
VlnPlot(object = aggr, features = c("PDGFRA"),cols=Col.set,group.by = "ClusterID_res0.6")
FeaturePlot(aggr,features = c("PDGFRA"),pt.size=0.8,label = T, split.by=("GroupID"))
FeaturePlot(aggr,features = c("SPN"),pt.size=0.8,label = T, split.by=("GroupID"))
FeaturePlot(aggr,features = c("CD44"),pt.size=0.8,label = T, split.by=("GroupID"))
VlnPlot(object = aggr, features = c("CD34","SPN","CD44","PTPRC","MYB","mEGFP"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("NFKBIA"),cols=Col.set,group.by = "ClusterID_res0.6",split.by=("GroupID"))
FeaturePlot(aggr,features = c("HIF1A"),pt.size=0.8,label = T, split.by=("GroupID"))
VlnPlot(object = aggr, features = c("HIF1A"),cols=Col.set,group.by = "ClusterID_res0.6",split.by=("GroupID"))


#
VlnPlot(object = aggr, features = c("ERBB4","ALPK2"),cols=Col.set,group.by = "ClusterID_res0.6")
FeaturePlot(aggr,features = c("GATA3","HAND1","ERBB4","ALPK2"),pt.size=0.8,label = T)
FeaturePlot(aggr,features = c("RAMP2"),pt.size=0.8,label = T)
VlnPlot(object = aggr, features = c("Podxl","Aldh1a2", "Ppargc1a"),cols=Col.set, group.by = "ClusterID_res0.6",split.by=("GroupID"))
VlnPlot(object = aggr, features = c("PODXL"),cols=Col.set, group.by = "ClusterID_res0.6",split.by=("GroupID"))
VlnPlot(object = aggr, features = c("ALDH1A2"),cols=Col.set, group.by = "ClusterID_res0.6",split.by=("GroupID"))
VlnPlot(object = aggr, features = c("PPARGC1A"),cols=Col.set, group.by = "ClusterID_res0.6",split.by=("GroupID"))
VlnPlot(object = aggr, features = c("ADM"),cols=Col.set, group.by = "ClusterID_res0.6",split.by=("GroupID"))
VlnPlot(object = aggr, features = c("SVEP1"),cols=Col.set, group.by = "ClusterID_res0.6",split.by=("GroupID"))

#
FeaturePlot(aggr,c("CPA3","HDC","GATA2","CLC"),label = T)
VlnPlot(object = aggr, features = c("CPA3","HDC"),cols=Col.set,group.by = "ClusterID_res0.6",split.by=("GroupID"))

# Cluster 12 Erythroid
VlnPlot(object = aggr, features = c("HBA1","HBA2","HBZ","HBG1","HBG2","HBE1","HBB","HBM","GYPA","GYPB","GYPC"),cols=Col.set,group.by = "ClusterID_res0.6")
FeaturePlot(aggr,features = c("HBA1","HBZ","HBG1","HBE1"),pt.size=0.8,label = T,split.by=("GroupID"))

#
FeaturePlot(aggr,c("DEFA4","MPO","RETN","PRTN3"),label = T)
VlnPlot(object = aggr, features = c("RETN","PRTN3"),cols=Col.set,,group.by = "ClusterID_res0.6")

#
VlnPlot(object = aggr, features = c("CD3E","CD3D"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("CD7","CD79A"),cols=Col.set,group.by = "ClusterID_res0.6")

#
FeaturePlot(aggr,c("PPBP","PF4","CLEC1B","PLEK","ITGA2B"),label = T)
FeaturePlot(aggr,c("PPBP","PF4","GYPA","HBA1"),label = T)
FeaturePlot(aggr,c("CLEC1B","TUBB1","THBS1","RGX18"),label = T)
VlnPlot(object = aggr, features = c("PPBP","PF4"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("GYPA","HBA1"),cols=Col.set,group.by = "ClusterID_res0.6")

FeaturePlot(aggr,c("ITGB2"),label = T,split.by=("GroupID"))
VlnPlot(object = aggr, features = c("PPBP","PF4","ITGA2B"),cols=Col.set,group.by = "ClusterID_res0.6")

#
#Cluster 10 Macrophages/Monocytes CD14
FeaturePlot(aggr,c("SPP1","CD74","C1QB","CTSB"),label = T)
FeaturePlot(aggr,c("HLA-DRA","AIF1","RAG1","KIT"),label = T)
FeaturePlot(aggr,c("FCGR3A","CD14","C1QA","CYBB"),label = T)
FeaturePlot(aggr,c("CD68","FCGR3A","CD14","ITGAM"),label = T)
FeaturePlot(aggr,c("SPP1","CD74","C1QB","C1QA","HLA-DRA","HLA-DRB5"),label = T)
VlnPlot(object = aggr, features = c("CD68","ITGAM"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("FCGR3A","CD14"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("SPP1","C1QA"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("CD74","HLA-DRA","HLA-DRB5","SPP1","C1QA","C1QB"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("TOP2A","CENPF","MKI67"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("mEGFP"),cols=Col.set,group.by = "ClusterID_res0.6",split.plot =TRUE,split.by = "GroupID")
VlnPlot(object = aggr, features = c("JAG1"),cols=Col.set,group.by = "ClusterID_res0.6",split.plot =TRUE,split.by = "GroupID")

#
FeaturePlot(aggr,c("POSTN"),label = T)
FeaturePlot(aggr,c("CCL5"),label = T)
FeaturePlot(aggr,c("KDR"),label = T)

# Endothlial
FeaturePlot(aggr,c("KDR","CDH5","CD34","GJA4","VWF"),label = T)
FeaturePlot(aggr,c("CALCRL","CD93","ENG","IGFBP4"),label = T)
VlnPlot(object = aggr, features = c("KDR","CD34"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("CDH5","GJA4"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("CD93","IGFBP4","HAPLN1"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("KDR","CD34", "NOTCH1","DLL4","GJA4","CXCR4","CDH5","SOX17","IGFBP5"),cols=Col.set,group.by = "ClusterID_res0.6")
FeaturePlot(aggr,c("JAG1","CLASP2","SVEP1","RAMP2","ADM","CALCRL"),label = T)
VlnPlot(object = aggr, features = c("CALCRL"),cols=Col.set,group.by = "ClusterID_res0.6",split.by = "GroupID")
VlnPlot(object = aggr, features = c("RAMP2"),cols=Col.set,group.by = "ClusterID_res0.6",split.by = "GroupID")
VlnPlot(object = aggr, features = c("SVEP1"),cols=Col.set,group.by = "ClusterID_res0.6",split.by = "GroupID")
VlnPlot(object = aggr, features = c("CLASP2"),cols=Col.set,group.by = "ClusterID_res0.6",split.by = "GroupID")
VlnPlot(object = aggr, features = c("NOTCH1"),cols=Col.set,group.by = "ClusterID_res0.6",split.by = "GroupID")
VlnPlot(object = aggr, features = c("BMP4"),cols=Col.set,group.by = "GroupID")
VlnPlot(object = aggr, features = c("GATA2"),cols=Col.set,group.by = "GroupID")
# BoxPlot(object = aggr, features = c("GATA2"),cols=Col.set,group.by = "GroupID")

#
VlnPlot(object = aggr, features = c("TNNT2","NKX2-5"),group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("MYH7","MYH7"),group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("NKX2-5"),cols=Col.set,group.by = "ClusterID_res0.6",split.by = "GroupID")
##Notch signalling for Cardiacmyocyte proliferation


#

VlnPlot(object = aggr, features = c("RAMP2"),cols=Col.set,group.by = "ClusterID_res0.6",split.by = "GroupID")
# cluster 15 & 10 & 17
VlnPlot(object = aggr, features = c("CD74","CD14","ITGAM","FCGR3A","CD58","CD68","HES1","ICAM1","HLA-DRB5"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("NFKBIA","CREBRF","CD58","HES1","CXCL5","CXCL3","ICAM1"),cols=Col.set,group.by = "ClusterID_res0.6")
VlnPlot(object = aggr, features = c("IL3","IL6","TPO"),cols=Col.set,group.by = "ClusterID_res0.6")

###########################################################
## Averge expression plot for clusters of split condition##
###########################################################
setwd("H:/R/R 09Apr2020/Data/remaped_data/")
source("plotting_functions_JL.R")
source("expression_functions.R")
cellnames=names(aggr@active.ident)
aggr <- AddMetaData(object = aggr, metadata = cellnames, col.name = "cell.names")
MarkerGenes=c("HBA1","HBA2","HBZ","HBG1","HBG2","HBE1","HBB","HBD","GYPA","GYPB","GYPC")
cluster.set=c("1", "12")
labels<- unique(as.character(aggr@meta.data$GroupID))
col.set=c("#ff0000","#ff0000","#ff0000","#ff0000")
col.set=NULL
doDotPlotSplit(aggr,MarkerGenes,labels,cluster.set,figure.title = NULL, do.plot = TRUE, 
               return.plot = FALSE, flip.plot = FALSE, dot.scale = 6, pal.num = 20,
               col.low = "#d9d9d9", col.set, col.max = 3, col.min = -3)

MarkerGenes=c("NOTCH1","JAG1","DLL4","NR2F2","FOXC2","ANG2","EFNB2")
MarkerGenes=c("PPBP","PF4","CLEC1B","PLEK","ITGA2B")

VlnPlot(object = aggr, features = c("HBA1","HBA2"),group.by = "ClusterID_res0.6",split.plot =TRUE,split.by = "GroupID")
VlnPlot(object = aggr, features = c("HBG1","HBG2"),group.by = "ClusterID_res0.6",split.plot =TRUE,split.by = "GroupID")
VlnPlot(object = aggr, features = c("HBZ","HBE1"),group.by = "ClusterID_res0.6",split.plot =TRUE,split.by = "GroupID")
VlnPlot(object = aggr, features = c("HBE1"),group.by = "ClusterID_res0.6",split.plot =TRUE,split.by = "GroupID")
VlnPlot(object = aggr, features = c("TIMP3","ESAM","RHOJ","DLL4","GATA2","ERG","mEGFP","RAG1"),group.by = "ClusterID_res0.6",split.plot =TRUE,split.by = "GroupID")
VlnPlot(object = aggr, features = c("CD34","SPN","CD44","LMO4","CD33","KIT","GATA2","PTPRC","RUNX1"),cols=Col.set,group.by = "ClusterID_res0.6")

###########################################################
## Monocle for haematoppoietic clusters ###################
###########################################################
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

# BiocManager::install("monocle")

# library(monocle)
# import Seurat object
importCDS("remaped_data/RemapedAggregate_clustered.Robj", import_all = TRUE)


