library("Seurat")
library("plyr")
library('ggplot2')
library(reshape2)
files <- list.files("./analysis_jane/raw/", pattern="*.csv", full.names=TRUE)
files <- sort(files)

# for(file in files){
#   text <- readLines(file)
#   write.csv(text,paste(file,".csv",sep = ""))
# } 


all_genelist <- list()
all_gene <- c()
for(f in files){
  print(f)
  if(file.size(f )!= 4 ){
    
    text <- read.csv(f,header = T)
    all_genelist <- c(all_genelist,list(text$Gene))
    print(head(text$Gene))
    all_gene <- c(all_gene,as.character( text$Gene))
  }else{
    all_genelist <- c(all_genelist,list(as.factor(c("NA"))))
  }

} 
print(all_gene)
print(all_genelist)

nlength <- sapply(all_genelist,length)
seqmax <- seq_len(max(nlength))
dataf <- t(sapply(all_genelist,"[",i = seqmax))
write.csv(dataf,"./analysis_jane/summary.csv")

all_gene <- unique(all_gene)
print(all_gene)
print(length(all_gene))



c13gene <- read.csv("./analysis_jane/c13_shearVSshearSBchir.csv")
f <- c13gene[order(c13gene$Pvalue),]
f <- f[1:10,]
DotPlot(subset(jane,idents = "13"),features = f$Gene,group.by = "GroupID",) +coord_flip()
VlnPlot(subset(jane,idents = "13"),group.by = "GroupID", features = f$Gene)
table(Idents(jane),jane$GroupID)
files <- list.files("./Results/GeneDiscovery/", pattern="*bd.txt$", full.names=TRUE)



t<-table(all_gene)
t <- as.data.frame(t)
t <- t[order(t$Freq,decreasing = T),]
