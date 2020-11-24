library("Seurat")
library("plyr")
library('ggplot2')
library(reshape2)
files <- list.files("./Results/GeneDiscovery/", pattern="*.csv", full.names=TRUE)


# for(file in files){
#   text <- readLines(file)
#   write.csv(text,paste(file,".csv",sep = ""))
# } 


all_gene <- list()
for(file in files){
  text <- read.csv(file,header = F)
  all_gene <- c(all_gene,list(text$V2))
} 
print(all_gene)

nlength <- sapply(all_gene,length)
seqmax <- seq_len(max(nlength))
dataf <- t(sapply(all_gene,"[",i = seqmax))
write.csv(dataf,"./Results/summary.csv")

all_gene <- unique(all_gene)
print(all_gene)







files <- list.files("./Results/GeneDiscovery/", pattern="*bd.txt$", full.names=TRUE)

all_gene2 <- c()
for(file in files){
  text <- readLines(file)
  if(length(text)!= 1 ){
    print(file)
    print(text)
    all_gene2 <- c(all_gene2,text)
    
  }
} 
print(all_gene2)
unique2 <- unique(all_gene2)
print(all_gene2)
print(unique2)
t<-table(all_gene2)
t <- as.data.frame(t)
t <- t[order(t$Freq,decreasing = T),]

write.csv(t,"Gene_frequency(bd).csv")

sheargene <- intersect(all_gene,all_gene2)
VlnPlot(GFP,sheargene,group.by = "orig.ident")

genelist <- c()

files <- list.files("./Results/GeneDiscovery/", pattern="^c0_.*csv$", full.names=TRUE)
for (file in files){
  t <- read.csv(file,header = T)
  genelist <- c(genelist, as.character(t$x))
  
}
featuresgene <- unique(genelist)
DoHeatmap(subset(GFP,idents = "0"),group.by = "orig.ident",features = featuresgene)
VlnPlot(subset(GFP,idents = "0"),group.by = "orig.ident", features = featuresgene)
DotPlot(subset(GFP,idents = "0"),group.by = "orig.ident",features = featuresgene) + coord_flip()

file <- read.csv('Results/cluster 0.csv')
ma <- subset(file,file$X %in% featuresgene)
ma <- ma[,1:5]
mma <- melt(ma)
heatmap(as.matrix(ma))
ggplot(data = mma, aes(x=variable,y=X,fill=value)) + geom_tile() + scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                                                                                        midpoint = 0, limit = c(-2,2), space = "Lab", 
                                                                                        name="Gene\nExpression")
