library("Seurat")
library("enrichR")
library("ggplot2")
abcd <- readRDS("updated_data/ABCD.Robj")
load("original_data/Aggregate_MIX.RData")

tsne <- UpdateSeuratObjedct(aggr)
rm(aggr)

DimPlot(tsne,label = T)


tsnecells <- as.data.frame(tsne@active.ident)
umapcells <- as.data.frame(abcd@active.ident)

rowname <- row.names(tsnecells)
tsnecells$cell_id <- rowname
write.csv(tsnecells,"tsne_cells_ID.csv")

rowname <- row.names(umapcells)
umapcells$cell_id <- rowname
write.csv(umapcells,"umap_cells_ID.csv")



##PYTHON EDIT EXTRA TEXT##



tsnecells <- read.csv("tsne_cells_ID.csv") #8851 cell
umapcells <- read.csv("umap_cells_ID.csv") #8893 cell


c <- 0
c1 <- 0
tsnemap <- data.frame(matrix(nrow = 0,ncol = 3))
umaponly <- data.frame(matrix(nrow = 0,ncol = 2))
colnames(tsnemap) <- c("Cell_id","umap_ident","tsne_ident")
colnames(umaponly) <- c("Cell_id","umap_ident")
for (i in 1:nrow(umapcells)) {
  
  if (umapcells$cell_id[i] %in% tsnecells$cell_id){
    index <- which(grepl(umapcells$cell_id[i],tsnecells$cell_id))
    newrow <- data.frame(Cell_id = umapcells$cell_id[i], umap_ident = umapcells$abcd.active.ident[i],tsne_ident = tsnecells$tsne.active.ident[index])
    tsnemap <- rbind(tsnemap,newrow)
    c <- c + 1
  }else{
    newrow <- data.frame(Cell_id = umapcells$cell_id[i], umap_ident = umapcells$abcd.active.ident[i])
    umaponly <- rbind(umaponly,newrow)
    c1 <- c1 + 1
  }
}
umaponly_summary <- table(umaponly$umap_ident)
#umaponly result c1 5 cells, c3 34 cells, c6 1 cell, c10 2 cells

write.csv(tsnemap,"tsne_umap_map.csv")

