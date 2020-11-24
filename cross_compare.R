library("Seurat")
library("plyr")
library('ggplot2')
library("devtools")
library("CSO")
install("../Rpackage//CSO")

setwd("../../THESIS/")


A <- readRDS("./original_data/GFP/STATIC_DMSO.Robj")
B <- readRDS("./original_data/GFP/SHEAR_DMSO.Robj")
C <- readRDS("./original_data/GFP/STATIC_SBCHIR.Robj")
D <- readRDS("./original_data/GFP/SHEAR_SBCHIR.Robj")

#a c1 static #b c3 shear #c c2 static+sbchir #d c4 Shear+ sbchir


GFP_Jgroup <- merge(SD,y=c(SS,STD,STS) ,add.cell.ids = c("C3","C4","C1","C2"), project = "GFP_with_jane_group")


GFP_Jgroup <- NormalizeData(GFP_Jgroup,normalization.method = "LogNormalize",scale.factor = 10000)
GFP_Jgroup <- FindVariableFeatures(GFP_Jgroup, selection.method = "vst",nfeatures = 2000)
all.genes <- rownames(GFP_Jgroup)
GFP_Jgroup <- ScaleData(GFP_Jgroup,features = all.genes)
GFP_Jgroup <- RunPCA(GFP_Jgroup, features = VariableFeatures(object = GFP_Jgroup))
ElbowPlot(GFP_Jgroup,ndims = 50) # suggest maybe 17
# GFP_Jgroup <- JackStraw(GFP_Jgroup, num.replicate = 100)
# GFP_Jgroup <- ScoreJackStraw(GFP_Jgroup,dims = 1:13)
# JackStrawPlot(GFP_Jgroup,dims = 1:13)
GFP_Jgroup <- FindNeighbors(GFP_Jgroup,dims = 1:17)
GFP_Jgroup <- FindClusters(GFP_Jgroup,resolution = 0.7)
# GFP_Jgroup <- RunTSNE(GFP_Jgroup,dims = 1:10)
GFP_Jgroup <- RunUMAP(GFP_Jgroup,dims = 1:17)
DimPlot(GFP_Jgroup,label = TRUE)
saveRDS(GFP_Jgroup,"./updated_data/GFP_Jgroup.Robj")
  
  
DimPlot(GFP,label = T)


#a c1 static #b c3 shear #c c2 static+sbchir #d c4 Shear+ sbchir

A <- RenameCells(A,add.cell.id = "C1")
A<-NormalizeData(A,normalization.method = "LogNormalize",scale.factor = 10000)
A <- FindVariableFeatures(A, selection.method = "vst",nfeatures = 2000)
all.genes <- rownames(A)
A <- ScaleData(A,features = all.genes)
A <- RunPCA(A, features = VariableFeatures(object = A))
ElbowPlot(A,ndims = 50) # suggest maybe 11
# A <- JackStraw(A, num.replicate = 100)
# A <- ScoreJackStraw(A,dims = 1:13)
# JackStrawPlot(A,dims = 1:13)
A <- FindNeighbors(A,dims = 1:11)
A <- FindClusters(A,resolution = 0.7)
# A <- RunTSNE(A,dims = 1:10)
A <- RunUMAP(A,dims = 1:11)
DimPlot(A,label = TRUE)
saveRDS(A,"./updated_data/A_GFP_j.Robj")
A <- readRDS("./updated_data/A_GFP_j.Robj")

B <- RenameCells(B,add.cell.id = "C3")
B<-NormalizeData(B,normalization.method = "LogNormalize",scale.factor = 10000)
B <- FindVariableFeatures(B, selection.method = "vst",nfeatures = 2000)
all.genes <- rownames(B)
B <- ScaleData(B,features = all.genes)
B <- RunPCA(B, features = VariableFeatures(object = B))
ElbowPlot(B,ndims = 50) # suggest maybe 11
# B <- JackStraw(B, num.replicate = 100)
# B <- ScoreJackStraw(B,dims = 1:13)
# JackStrawPlot(B,dims = 1:13)
B <- FindNeighbors(B,dims = 1:11)
B <- FindClusters(B,resolution = 0.7)
# B <- RunTSNE(B,dims = 1:10)
B <- RunUMAP(B,dims = 1:11)
DimPlot(B,label = TRUE)
saveRDS(B,"./updated_data/B_GFP.Robj")

C <- RenameCells(C,add.cell.id = "C2")
C<-NormalizeData(C,normalization.method = "LogNormalize",scale.factor = 10000)
C <- FindVariableFeatures(C, selection.method = "vst",nfeatures = 2000)
all.genes <- rownames(C)
C <- ScaleData(C,features = all.genes)
C <- RunPCA(C, features = VariableFeatures(object = C))
ElbowPlot(C,ndims = 50) # suggest maybe 12
# C <- JackStraw(C, num.replicate = 100)
# C <- ScoreJackStraw(C,dims = 1:13)
# JackStrawPlot(C,dims = 1:13)
C <- FindNeighbors(C,dims = 1:12)
C <- FindClusters(C,resolution = 0.7)
# C <- RunTSNE(C,dims = 1:10)
C <- RunUMAP(C,dims = 1:12)
DimPlot(C,label = TRUE)
saveRDS(C,"./updated_data/C_GFP.Robj")

D <- RenameCells(D,add.cell.id = "C4")
D<-NormalizeData(D,normalization.method = "LogNormalize",scale.factor = 10000)
D <- FindVariableFeatures(D, selection.method = "vst",nfeatures = 2000)
all.genes <- rownames(D)
D <- ScaleData(D,features = all.genes)
D <- RunPCA(D, features = VariableFeatures(object = D))
ElbowPlot(D,ndims = 50) # suggest maybe 13
# D <- JackStraw(D, num.replicate = 100)
# D <- ScoreJackStraw(D,dims = 1:13)
# JackStrawPlot(D,dims = 1:13)
D <- FindNeighbors(D,dims = 1:12)
D <- FindClusters(D,resolution = 0.7)
# D <- RunTSNE(D,dims = 1:10)
D <- RunUMAP(D,dims = 1:12)
DimPlot(D,label = TRUE)
saveRDS(D,"./updated_data/D_GFP.Robj")

A_anotated <- A
DimPlot(A_anotated,label = TRUE)
#cluster 0 Stromal
FeaturePlot(A_anotated,c("COL8A1","IGFBP7","ID4","ACTA2"),label = T,pt.size = 1)
A_anotated <- RenameIdents(A_anotated,"0" = "0.Stromal")

#cluster 1 MEP
FeaturePlot(A_anotated,c("KEL","TSC22D1","PPP1R15A"),label = T)
A_anotated <- RenameIdents(A_anotated,"1" = "1.MEP")

#cluster 2 Megakaryocyte
FeaturePlot(A_anotated,c("PPBP","PF4","SPN","GFI1B"),label = T)
A_anotated <- RenameIdents(A_anotated,"2" = "2.Megakaryocyte")

#cluster 3 Common myeloid progenitor
FeaturePlot(A_anotated,c("SAMSN1","S100A4","CLC","SRGN"),label = T)
A_anotated <- RenameIdents(A_anotated,"3" = "3.CMP")

#cluster 4 c7+ Myofibroblast
A_anotated <- RenameIdents(A_anotated,"4" = "4.Myofibroblast")

#cluster 5 CM
FeaturePlot(A,c("MYL7"),label = T)
A_anotated <- RenameIdents(A_anotated,"5" = "5.Cardiac Muscel")

#cluster 6 Mesenchymal stromal cell HOXA+PDGFRa+ highly proliferative
FeaturePlot(A_anotated,c("TOP2A","MKI67","HOXA10","PDGFRA"),label = T)
FeaturePlot(A_anotated,c("NUSAP1","HIST1H4C","TPX2","PRC1"),label = T)
A_anotated <- RenameIdents(A_anotated,"6" = "6.HP Stromal")

#cluster 7 MEP
FeaturePlot(A_anotated,c("RP3-395M20.9","KLF10","RP11-62J1.4","CTA-29F11.1","CREBRF","MYL6"),label = T)

#cluster 8 Endothelial
FeaturePlot(A_anotated,c("CD93","ERG","CD34","KDR"),label = T)
A_anotated <- RenameIdents(A_anotated,"8" = "8.Endothelial")

#cluster 9 
FeaturePlot(A_anotated,c("PTN","TGFBI"),label = T)

#cluster 10 sympathoadrenal progenitors
FeaturePlot(A_anotated,c("ERBB4","MSX2","LINC00458","PRTG","GATA3"),label = T)
A_anotated <- RenameIdents(A_anotated,"10" = "10.sympathoadrenal")

#cluster 11 MEP

#cluster 12 Erythroid
FeaturePlot(A_anotated,c("HBZ"))
A_anotated <- RenameIdents(A_anotated,"12" = "12.Erythroid")

#cluster 13 Macrophage
FeaturePlot(A_anotated,c("SPP1","CD74","C1QB","CTSB","CSF1R","CD14"),label = T)
A_anotated <- RenameIdents(A_anotated,"13" = "13.Macrophage")

DimPlot(A_anotated,label = T)
saveRDS(A_anotated,"./updated_data/A_GFP_anotated_j.Robj")


A_markers <- FindAllMarkers(A)

DimPlot(A,label = TRUE,label.size = 8)
DimPlot(B,label = TRUE,label.size = 8)
DimPlot(C,label = TRUE,label.size = 8)
DimPlot(D,label = TRUE,label.size = 8)








