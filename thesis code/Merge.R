library("Seurat")

# load the data 
A <- readRDS("updated_data/A_Shear_culture+drug_treatment.Robj")
B <- readRDS("updated_data/B_Shear_culture.Robj")
C <- readRDS("updated_data/C_Static_culture+drug_treatment.Robj")
D <- readRDS("updated_data/D_Static_culture.Robj")



# set the meta original data into the seurat object, make things easier later on.
A$original <- A$orig.ident
A$orig.ident <- NULL
A <- AddMetaData(object = A, metadata = "D.Shear+SB/CHIR", col.name = "orig.ident")
head(A[[]])

B$original <- B$orig.ident
B$orig.ident <- NULL
B <- AddMetaData(object = B, metadata = "B.Shear", col.name = "orig.ident")
head(B[[]])

C$original <- C$orig.ident
C$orig.ident <- NULL
C <- AddMetaData(object = C, metadata = "C.static+SB/CHIR", col.name = "orig.ident")
head(C[[]])

D$original <- D$orig.ident
D$orig.ident <- NULL
D <- AddMetaData(object = D, metadata = "A.Static", col.name = "orig.ident")
head(D[[]])


abcd <- merge(A,y = c(B,C,D),add.cell.ids = c("C4","C3","C2","C1"),project = "all")

#save the raw unclusterd data
saveRDS(abcd,"updated_data/ABCD.Robj")


#pipeline for umap
abcd <- FindVariableFeatures(abcd, selection.method = "vst",nfeatures = 2000)
all.genes <- rownames(abcd)
abcd <- ScaleData(abcd,features = all.genes)
abcd <- RunPCA(abcd, features = VariableFeatures(object = abcd))
ElbowPlot(abcd,ndims = 50)
abcd <- JackStraw(abcd, num.replicate = 100)
abcd <- ScoreJackStraw(abcd,dims = 1:20)
JackStrawPlot(abcd,dims = 1:20)
abcd <- FindNeighbors(abcd,dims = 1:10)
abcd <- FindClusters(abcd,resolution = 0.7)
# abcd <- RunTSNE(abcd,dims = 1:10)
abcd <- RunUMAP(abcd,dims = 1:10)
DimPlot(abcd,label = T,split.by = "orig.ident")
saveRDS(abcd,"updated_data/ABCD.Robj")

DimPlot(abcd,label = T,coord.fixed = T)

# ct <- Idents(abcd)
# abcd$replicate <- sample()
# Idents(abcd) <- "replicate"
# DimPlot(abcd,split.by = "orig.ident")
