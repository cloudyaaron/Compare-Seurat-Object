library("Seurat")
c5 <- subset(abcd,idents = 5)
c5 <- FindVariableFeatures(c5, selection.method = "vst",nfeatures = 2000)
all.genes <- rownames(c5)
c5 <- ScaleData(c5,features = all.genes)
c5 <- RunPCA(c5, features = VariableFeatures(object = c5))
ElbowPlot(c5,ndims = 50)
c5 <- JackStraw(c5, num.replicate = 100)
c5 <- ScoreJackStraw(c5,dims = 1:20)
JackStrawPlot(c5,dims = 1:20)
c5 <- FindNeighbors(c5,dims = 1:10)
c5 <- FindClusters(c5,resolution = 0.5)
c5 <- RunUMAP(c5,dims = 1:10)
DimPlot(c5,label = T,split.by = "orig.ident")
# PC_ 1 
# Positive:  FN1, HAPLN1, CGNL1, LAMB1, ==="FBN1"=== embryonic semilunar valves, S100A10, EFNB2, SPARC, ==="RGS5"===, DSTN 
# MAP1B, IFI16, ITGAV, MYH10, COL4A1, EDN1, BGN, GUCY1A3, PCDH7, OPTN 
# LTBP1, SERPINE2, TPM4, S100A6, IGFBP7, COL1A2, COL1A1, VIM, COL4A2, GBP4 
# Negative:  ==="IGFBP5"===https://pubmed.ncbi.nlm.nih.gov/15010534/, CRABP1, PARM1, NEFM, PDE1A, NR2F1, FRZB, SFRP2, MEST, GRIA4 
# SIX1, MME, THSD7A, MOXD1, COL9A2, COL2A1, ZFHX4, COL24A1, TSHZ2, WNT5A 
# NKX2-3, CD24, IGDCC3, NEFL, KIF26B, CXXC4, MAB21L2, CLC, FIBIN, SLC26A7 
# PC_ 2 
# Positive:  LGALS1, COL4A1, F2R, COL1A2, COL1A1, SPRY1, PTN, ITGA1, IGFBP7, TGFBI 
# COL11A1, BGN, CTGF, TIMP1, ITGA4, S100A6, CAVIN1, COL4A2, CYTOR, FLT1 
# LAMA4, SPARC, SEMA5A, PCDH18, FBN1, TM4SF1, FILIP1L, STMN2, CTHRC1, LRRC17 
# Negative:  ERBB4, NEBL, PI15, ALPK2, MSX2, DYNC1I1, TBR1, PRTG, DAB1, LINC00458 
# CYSLTR2, DOK4, KRT19, NRP2, MYL7, GATA3, ARL6IP1, PCDH10, AP003716.1, CENPF 
# RHOBTB3, TMEM88, TRDN, MYEF2, LEF1, LINC00937, SERINC5, TOP2A, NRP1, PKP2 


c9 <- subset(abcd,idents = 9)
c9 <- FindVariableFeatures(c9, selection.method = "vst",nfeatures = 2000)
all.genes <- rownames(c9)
c9 <- ScaleData(c9,features = all.genes)
c9 <- RunPCA(c9, features = VariableFeatures(object = c9))
ElbowPlot(c9,ndims = 20)
c9 <- JackStraw(c9, num.replicate = 100)
c9 <- ScoreJackStraw(c9,dims = 1:20)
JackStrawPlot(c9,dims = 1:20)
c9 <- FindNeighbors(c9,dims = 1:5)
c9 <- FindClusters(c9,resolution = 0.5)
c9 <- RunUMAP(c9,dims = 1:5)
DimPlot(c9,label = T)
# PC_ 1 
# Positive:  COL1A2, COL3A1, HAPLN1, GBP4, TPM1, EDN1, CGNL1, SLIT3, NPM1, RPL14 
# RPL12, LTBP1, VCAN, FBN1, DDR2, CXCL12, ADGRG6, PCDH7, NAP1L1, MT-ND5 
# RPS25, RGS5, C7, RPS6, LAMB1, RPL6, HMGA2, ECE1, ADGRL2, BEX1 
# Negative:  IGFBP5, LYVE1, CD34, CLDN5, TMSB4X, LPL, A2M, CD74, ADGRF5, STC1 
# CYP1B1, COL15A1, HLA-DRA, AC245060.5, CMKLR1, SCN9A, SLC9A3R2, NEAT1, CLEC14A, THBD 
# C8orf4, HLA-DRB1, NQO1, PMP22, NNMT, TM4SF18, MN1, MALL, LHX6, CYP1A1 
# PC_ 2 
# Positive:  CD34, ADGRF5, TIMP3, TMSB4X, PECAM1, GNG11, CALM1, S100A6, CLEC14A, COL15A1 
# LAMA4, CAVIN1, PMP22, TM4SF18, HOXB7, CD59, ZFHX3, THSD7A, NOTCH4, C8orf4 
# CAV1, CRIM1, MT-ATP8, IGFBP5, LYVE1, MCAM, CLDN5, AFAP1L1, MT-ND2, PLK2 
# Negative:  COL3A1, COL1A2, VCAN, CXCL12, TPM1, SLIT3, GBP4, COL1A1, CGNL1, HAND2 
# TGFB2, EDN1, LRRC17, LTBP1, CDH11, FBN1, PTN, NCOA7, GATA4, HAPLN1 
# CDH2, EMILIN1, COL21A1, SYT4, CCDC80, C7, PCOLCE, RND3, SEMA5A, ITGAV 
