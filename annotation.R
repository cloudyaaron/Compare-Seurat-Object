library("Seurat")
library("enrichR")
library("ggplot2")
abcd <- readRDS("updated_data/ABCD.Robj")

#ABCD
#paraxial mesoderm##
##
OLD_anotate <- abcd
########################################
#cluster 0 Mesochymal stromal cell, Stromal-1 HOXA+PDGFRa+
FeaturePlot(OLD_anotate,c("WNT5A"),label = T)
FeaturePlot(OLD_anotate,c("DNM3OS","LUM","MEST","HOXA11"))
OLD_anotate <- RenameIdents(OLD_anotate,"0" = "0.Stromal HOXA+PDGFRa+")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 1 Stromal or early cardiac muscle (multipotent secondary-heart-field (SHF) progenitors)
FeaturePlot(OLD_anotate,c("TTN","ISL1")) #https://www.nature.com/articles/news.2009.522
OLD_anotate <- RenameIdents(OLD_anotate,"1" = "1.Early Cardiac muscle")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 2 Stromal-2 HOXA-PDGFRa+
FeaturePlot(OLD_anotate,c("C7"))
OLD_anotate <- RenameIdents(OLD_anotate,"2" = "2.Stromal HOXA-PDGFRa+")
DimPlot(OLD_anotate,label = T)


########################################
#cluster 3 Mesoderm stromal Stromal HOXA+PDGFRa+
FeaturePlot(OLD_anotate,c("HOXA10","PDGFRA"),label = T,split.by = "orig.ident",pt.size = 1)
OLD_anotate <- RenameIdents(OLD_anotate,"3" = "3.Stromal HOXA+PDGFRa+")
DimPlot(OLD_anotate,label = T)


########################################
#cluster 4 Stromal PDGRA+
FeaturePlot(OLD_anotate,c("TGFBI","COL12A1","ACTA2","COL11A1"))
#TGFBI: tranforming growth factor, May resultin endochondrial ossification.
#ACTA2 gene : aortic smooth muscle gene, mutation cause variety of cascular diseases.
OLD_anotate <- RenameIdents(OLD_anotate,"4" = "4.Stromal HOXA+PDGFRa+")
DimPlot(OLD_anotate,label = T)


########################################
#cluster 5 assumming some sort of multipotent cell close related to c10
# 
#SIX2 gene is necessary for nephrogenesis during murine development
FeaturePlot(abcd,c("MKI67","PBK"),label = T)
#MKI67: for cell proliferation
OLD_anotate <- RenameIdents(OLD_anotate,"5" = "5.PDGFRa+ Endothelial to Mesenchymal transition")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 6 MEP
FeaturePlot(abcd,c("KEL","TSC22D1","PPP1R15A"),label = T)
OLD_anotate <- RenameIdents(OLD_anotate,"6" = "6.Megakaryocyte and erthroid progenitor")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 7 Megakarypcyte    #very important GFI1B:hematopoietic lineage, control development and maturation of erythrocytes and megakaryocytes GENECARD
FeaturePlot(abcd,c("PPBP","PF4","SPN","GFI1B"),label = T)
OLD_anotate <- RenameIdents(OLD_anotate,"7" = "7.Megakarypcyte")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 8 Macrophage
FeaturePlot(OLD,c("MRC1"))
OLD_anotate <- RenameIdents(OLD_anotate,"8" = "8.Macrophage")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 9 endothelium
readRDS("updated_data/ABCD.Robj")
FeaturePlot(abcd,c("CD93","ERG"))
OLD_anotate <- RenameIdents(OLD_anotate,"9" = "9.Endothelium")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 10 Myeloid-12
FeaturePlot(abcd,c("LYZ","CPA3","CD24","CLC","S100A8","PRG3"),label = T)
OLD_anotate <- OLD_anotate(OLD_anotate,"10" = "10.Myeloid")
DimPlot(OLD_anotate,label = T)


########################################
#cluster 11 Erythroid
FeaturePlot(abcd,c("HBZ"))
OLD_anotate <- RenameIdents(OLD_anotate,"11" = "11.Erythroid")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 12 AGM MP cell
#lot's of gene in c10 is also express in c5
FeaturePlot(abcd,c("SPINK2","AIF1","BUB1","C1QBP","CCNA2"))
#CACYB  #MYB from jane
#AIF1: induced by cytokines and interferon and may promote macrophage activation and growth of vascular smooth muscle cells and T-lymphocytes.https://doi.org/10.1634/stemcells.2005-0185
#SPINK2: find a lot of express in Cordblood, high expressed will results in acute myeloid leukemia patient death https://doi.org/10.3892/ol.2019.10665
#BUB1: Budding Uninhibited By Benzimidazoles, (Mitotic Checkpoint Serine)
#C1QBP: the first component of the serum complement system
#CCNA2:  binds and activates cyclin-dependent kinase 2 and thus promotes transition through G1/S and G2/M. 
OLD_anotate <- RenameIdents(OLD_anotate,"12" = "12.AGM_MPP")
DimPlot(OLD_anotate,label = T)

# for (gene in c10$gene) {
#   FeaturePlot(abcd,gene,label = T)
#   ggsave(filename = paste("c10/",gene,".png",sep = ""))
# }

######################################## 
#cluster 13 Common myeloid progenitors
#https://doi.org/10.1182/blood.V126.23.2768.2768
FeaturePlot(abcd,c("GATA2","CPA3","HDC"))
OLD_anotate <- RenameIdents(OLD_anotate,"13"="13.CMP")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 14 renal stromal progenitors
#doi: 10.1097/PAS.0b013e3182a0218f   https://cordis.europa.eu/project/id/302739/reporting
FeaturePlot(OLD_anotate,c("ERBB4","MSX2","TBR1","PRTG","GATA3"))
#maybe kidney
OLD_anotate <- RenameIdents(OLD_anotate,"14" = "14.renal stromal progenitors")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 15 Myeloid-3
# FeaturePlot(abcd,c("MYL7"))
OLD_anotate <- RenameIdents(OLD_anotate,"15"="15.Myeloid")
DimPlot(OLD_anotate,label = T)

########################################
#cluster 16 CM
FeaturePlot(OLD_anotate,c("MYL7"))
OLD_anotate <- RenameIdents(OLD_anotate,"16" = "16.Cardiac Muscel")
DimPlot(OLD_anotate,label = T)


########################################
#intra paraxial mesoidem
#cluster 0 Somitic https://doi.org/10.1002/dvdy.21787           AGM stroma cell
#TBX18 gene : this gene may allow some progenitor cell develop into smooth muscel progenitor. https://doi.org/10.1113/JP270033
FeaturePlot(abcd,c("TBX18","ACTA2"))
FeaturePlot(abcd,c("DNM3OS","LUM","MEST","HOXA11"))
abcd <- RenameIdents(abcd,"0" = "Skeletal system?")
abcd <- RenameIdents(abcd,"Skeletal system?" = "Somite")
DimPlot(abcd,label = T)



########################################
DimPlot(OLD_anotate,label = T)
DimPlot(abcd,label = T,split.by = "orig.ident")
saveRDS(OLD_anotate,"updated_data/OLD_anotate.Robj")
OLD_anotate <- readRDS("updated_data/OLD_anotate.Robj")
########################################
abcd_marker <- FindAllMarkers(abcd)
write.csv(abcd_marker,"Results/abcd_marker.csv")


v512 <- FindMarkers(abcd,ident.1 = "5",ident.2 = "12.AGM_MPP",logfc.threshold = 0.5)
print(head(v512))
VlnPlot(abcd,row.names(v512))

# change of original data
# D - > A
# B -> B
# C -> C   sbchir
# A -> D

# find difference between 2 & 4
d2_4 <- FindMarkers(abcd,ident.1 = "2",ident.2 = "4")
d12_5 <- FindMarkers(abcd,ident.1 = "12", ident.2 = "5")
# c12 has col type gene
# FN1 less express in D

# for (gene in row.names(d12_5)) {
#   f <- paste(gene,".png",sep = "")
#   t <- FeaturePlot(abcd,gene)
#   ggsave(f,t)
# }

#GYPA & BLVRB gene express in the red blood cell type