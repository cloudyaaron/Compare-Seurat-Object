DimPlot(GFP_anotated,label = T)

FeaturePlot(GFP_anotated,c("PPBP","PF4","SPN","GFI1B"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"7" = "7.Megakarypcyte")
DimPlot(GFP_anotated,label = T)


FeaturePlot(GFP_anotated,c("WNT5A","HOXA10","PDGFRA"),label = T,pt.size = 1)



########################################
#cluster 0 Megakarypcyte    #very important GFI1B:hematopoietic lineage, control development and maturation of erythrocytes and megakaryocytes GENECARD
FeaturePlot(GFP_anotated,c("PPBP","PF4","SPN","GFI1B"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"0" = "0.Megakaryocyte")
DimPlot(GFP_anotated,label = T)
########################################
#cluster 1 mesodermal pericytes,  HOXA+PDGFRa+WNT5+
FeaturePlot(GFP_anotated,c("WNT5A"),label = T)
FeaturePlot(GFP_anotated,c("DNM3OS","LUM","MEST","HOXA10"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"1" = "1.mesodermal pericytes ")
DimPlot(GFP_anotated,label = T)
########################################
#cluster 2 Stromal-1 HOXA-PDGFRa+C7+
FeaturePlot(GFP_anotated,c("C7","PDGFRA","HOXA10","ACTA2"),label = T,pt.size = 1)
GFP_anotated <- RenameIdents(GFP_anotated,"2" = "2.Myofibroblast")
DimPlot(GFP_anotated,label = T)
########################################
#cluster 3 Stromal-1 HOXA+PDGFRa++
FeaturePlot(GFP_anotated,c("TGFBI","PDGFRB","HOXA10","PDGFRA","CRABP1"),label = T,pt.size = 1)
GFP_anotated <- RenameIdents(GFP_anotated,"3" = "3.Stromal-(CRABP1+)")

########################################
#cluster 4 Stromal-2 HOXA-PDGFRb++
#has smooth muscel signal
FeaturePlot(GFP_anotated,c("TGFBI","PDGFRB","HOXA10","PDGFRA","TBX18"),label = T,pt.size = 1)
GFP_anotated <- RenameIdents(GFP_anotated,"4" = "4.Stromal-(smooth muscle)")
########################################
# cluster 5 MEP Megakaryocyte and erthroid progenitor
FeaturePlot(GFP_anotated,c("KEL","TSC22D1","PPP1R15A"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"5" = "5.MEP")
DimPlot(GFP_anotated,label = T)
########################################
#cluster 6 Stromal SCG2+
FeaturePlot(GFP_anotated,c("SCG2","HOXA10","ACTA2","PDGFRA"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"6" = "6.Stromal-(SCG2+)")
DimPlot(GFP_anotated,label = T)

########################################
# cluster 7 AGM MP
FeaturePlot(GFP_anotated,c("CD34","PTPRC","MKI67","SPN"),label = T)
FeaturePlot(GFP_anotated,c("SPINK2","AIF1","BUB1","C1QBP","CCNA2"))
GFP_anotated <- RenameIdents(GFP_anotated,"7" = "7.AGM MPP")
DimPlot(GFP_anotated,label = T)
########################################
# cluster 8 Macrophages 
FeaturePlot(GFP_anotated,c("SPP1","CD74","C1QB","CTSB","CSF1R","CD14"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"8" = "8.Macrophage")
DimPlot(GFP_anotated,label = T)
########################################
# cluster 9 Mesenchymal stromal cell HOXA+PDGFRa+ highly proliferative
FeaturePlot(GFP_anotated,c("TOP2A","MKI67","HOXA10","PDGFRA"),label = T)
FeaturePlot(GFP_anotated,c("NUSAP1","HIST1H4C","TPX2","PRC1"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"9" = "9.Mesenchymal stromal-1")
DimPlot(GFP_anotated,label = T)
########################################
# cluster 10 Aortic endothelium / Haematopoietic progenitor cells
FeaturePlot(GFP_anotated,c("CD93","ERG","CD34","KDR"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"10" = "10.Aortic endothelium")
DimPlot(GFP_anotated,label = T)
########################################
# cluster 11 Mesenchymal stromal cell HOXA-PDGFRa+ highly proliferative
FeaturePlot(GFP_anotated,c("TOP2A","MKI67","HOXA10","PDGFRA"),label = T)
FeaturePlot(GFP_anotated,c("NUSAP1","HIST1H4C","TPX2","PRC1"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"11" = "11.Mesenchymal stromal-2")
DimPlot(GFP_anotated,label = T)
########################################
# cluster 12 Common myeloid progenitors
FeaturePlot(GFP_anotated,c("SAMSN1","S100A4","CLC","SRGN"),label = T)
FeaturePlot(GFP_anotated,c("GATA2","CPA3","HDC"))
GFP_anotated <- RenameIdents(GFP_anotated,"12" = "12.Common myeloid progenitors")
DimPlot(GFP_anotated,label = T)
########################################
#cluster 13 CM
FeaturePlot(GFP_anotated,c("MYL7"))
GFP_anotated <- RenameIdents(GFP_anotated,"13" = "13.Cardiac Muscle")
GFP_anotated <- RenameIdents(GFP_anotated,"13.Cardiac Muscel" = "13.Cardiac Muscle")

DimPlot(GFP_anotated,label = T)
########################################
#cluster 14 sympathoadrenal progenitors 
FeaturePlot(GFP_anotated,c("ERBB4","MSX2","TBR1","PRTG","GATA3"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"14" = "14.sympathoadrenal progenitors")
DimPlot(GFP_anotated,label = T)
########################################
#cluster 15 esoderm stromal Stromal HOXA+PDGFRa+
FeaturePlot(GFP_anotated,c("HOXA10","PDGFRA","IGJ","C20orf112"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"15" = "15.Stromal")
DimPlot(GFP_anotated,label = T)
########################################
#cluster 16 Erythroid
FeaturePlot(GFP_anotated,c("HBZ"))
GFP_anotated <- RenameIdents(GFP_anotated,"16" = "16.Erythroid")
DimPlot(GFP_anotated,label = T)
########################################
#cluster 17 Macrophage
FeaturePlot(GFP_anotated,c("SPP1","CD74","C1QB","CTSB","CSF1R","CD14"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"17" = "17.Macrophage")
DimPlot(GFP_anotated,label = T)
########################################
#cluster 18 Myeloid
FeaturePlot(GFP_anotated,c("LYZ","CPA3","CD24","CLC","S100A8","PRG3"),label = T)
GFP_anotated <- RenameIdents(GFP_anotated,"18" = "18.Myeloid")
DimPlot(GFP_anotated,label = T)


DimPlot(GFP_anotated,label = T,pt.size = 1,order = c("0.Megakaryocyte","1.mesodermal pericytes","2.Myofibroblast","3.Stromal-(CRABP1+)","4.Stromal-(smooth muscle)","5.MEP","6.Stromal-(SCG2+)","7.AGM MPP",
                                                     "8.Macrophage","9.Mesenchymal stromal-1","10.Aortic endothelium","11.Mesenchymal stromal-2","12.Common myeloid progenitors","13.Cardiac Muscel",
                                                     "14.sympathoadrenal progenitors","15.Stromal","16.Erythroid","17.Macrophage","18.Myeloid"))




# PC_ 1 
# Positive:  RGS10, GPX1, PLEK, RGS18, PPBP, BTK, FERMT3, PF4, ITGA2B, CLEC1B 
# TUBB1, CD84, GP9, BIN2, RAB27B, MPP1, RHAG, SELP, ITGB3, HBG2 
# CNST, C2orf88, PSTPIP2, FCER1G, FYB, GP1BA, HBG1, ARHGAP6, GRAP2, MMRN1 
# Negative:  COL1A2, COL3A1, TPM1, MYH10, VCAN, COL1A1, FN1, COL4A1, LAMB1, COL4A2 
# SPARC, MAP1B, COL5A1, FBN1, CDH2, IGFBP7, DSTN, TPM2, CSRP2, DSP 
# PALLD, PCDH7, PLS3, CCDC80, HAND2, CGNL1, CYR61, NES, PRSS23, IGF2 
# PC_ 2 
# Positive:  CYBB, CTSS, C1QB, CD74, MS4A6A, C1QC, LCP1, C1QA, HLA-DRA, IGSF6 
# SPP1, CD36, S100B, ITGB2, SAMHD1, PTPRC, PLA2G7, HLA-DRB1, RNASE6, MS4A7 
# IFI30, HCK, CSF1R, HLA-DRB5, LY86, CCR1, CAPG, RASSF4, FCGR3A, FTL 
# Negative:  CLEC1B, PPBP, TUBB1, NEXN, PF4, DGKI, ITGB3, GP9, THBS1, SELP 
# PRKAR2B, CNST, GP1BA, ITGA2B, GRAP2, GNG11, NRGN, RAB27B, ARHGAP6, C2orf88 
# MTURN, CD226, CCL5, MYLK, MYOM1, TMEM40, SDPR, RGS18, MMRN1, RHAG 
# PC_ 3 
# Positive:  COL1A1, COL3A1, COL1A2, SPARC, COL6A3, FBN1, COL11A1, FN1, CTHRC1, LUM 
# TGFBI, PTN, SULF1, COL5A1, CXCL12, BGN, DLK1, CYR61, LAMB1, COL4A1 
# VCAN, MYH10, ITGA1, IGFBP3, POSTN, LRRC17, SERPINE2, MME, IGFBP7, SDC2 
# Negative:  SMPX, TNNC1, ACTN2, CTD-2545M3.8, TNNT2, MYL7, MYH7, MYBPC3, TRIM55, MYL3 
# CSRP3, ACTC1, NEBL, SYNPO2L, UNC45B, TNNI1, MYH6, MYL4, SORBS2, TRDN 
# SMYD1, CKM, LDB3, MASP1, ENO3, ALPK2, TTN, CMYA5, FABP3, PLN 
# PC_ 4 
# Positive:  C1QB, C1QA, C1QC, MAF, PLA2G7, SPP1, FCGR3A, MS4A7, RASSF4, CD36 
# MS4A6A, SLCO2B1, CD14, LGMN, IGSF6, CTSB, FOLR2, VSIG4, NPL, SIGLEC1 
# CTSS, CD163, CMKLR1, CTSD, TREM2, CTSL, MS4A4A, LY86, HMOX1, CCR1 
# Negative:  TOP2A, MKI67, TPX2, NUSAP1, CENPF, PRC1, HMGB2, CDK1, KIF11, ASPM 
# PTTG1, DLGAP5, NDC80, BUB1, NUF2, UBE2C, BUB1B, CCNA2, CENPE, HMMR 
# CCNB1, SMC4, KIF4A, SPC25, KIF2C, KIAA0101, ANLN, CDCA8, ARHGAP11A, KIF14 
# PC_ 5 
# Positive:  DAB2, FN1, GNG11, SERPINE2, LGMN, SPARC, C1QB, TM4SF1, HAPLN1, LAMB1 
# BGN, C1QC, C1QA, CTSB, CGNL1, FRMD4B, CD59, FBN1, ITGAV, CTSZ 
# TLN1, F13A1, LIMS1, EFNB2, PDLIM3, TIMP1, ANXA2, SPP1, ECE1, ITGA1 
# Negative:  HBZ, HBA2, ALAS2, GYPA, HDC, HBM, HBA1, CLC, AHSP, GYPB 
# MYB, CPA3, GAPT, CA1, MS4A3, HBE1, MYBPC3, GATA2, TNNC1, MYL3 
# SPTA1, ACTN2, PRSS57, SMPX, MYL4, C1orf186, CTD-2545M3.8, TRIM55, MYH7, S100P 