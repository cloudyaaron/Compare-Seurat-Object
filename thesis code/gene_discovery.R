library("Seurat")
library(datasets)
library("futile.matrix")

# gfp_csv <- as.data.frame(as.matrix(GFP@scale.data))
# write.table(as.matrix(GetAssayData(GFP,slot = "scale.data")),'./gfp.csv',sep = ",",row.names = T,col.names = T,quote = F)


# temp <- FeaturePlot(GFP,"IGFBP5",split.by = "orig.ident",combine = T)

gdbg <- AverageExpression(GFP,add.ident = "orig.ident",slot = "scale.data")
g <-  gdbg$RNA
write.csv(g,file = "raw_average_expression(scaled).csv")
genes <- rownames(g)

maxn <- 19
i <- 0
colindex <- colnames(g)
while(i < 19){
  i2 <- 1
  df <- data.frame(matrix(nrow = 23597,ncol = 0))
  row.names(df) <- genes
  while (i2 < length(colindex)) {
    if(grepl(paste("^",i,"_",sep = ""),colindex[i2])){
      print(paste(i,colindex[i2]))
      dft <- as.data.frame( g[,colindex[i2]] )
      print(head(dft ))
      colnames(dft) <- colindex[i2]
      df <- cbind(df,dft)
    }
    i2 <- i2 + 1
  }
  write.csv(df,file = paste("./Results/cluster ",i,".csv",sep = ""))
  i <- i + 1
  print("")
}


f <- read.csv("Results/cluster 0.csv",header = TRUE,row.names = 1,sep = ",")


rown <- 1
maxn <- length(rownames(f))
while (rown <= maxn) {
  coln <- 1
  maxcol <- length(colnames(f[rown,]))
  total <- 0
  avg <- 0
  while(coln <= maxcol){
    total <- total + f[rown,coln]
    coln <- coln + 1
  }
  avg <- total/maxcol
  if (total > 1){
    coln <- 1
    while(coln <= maxcol){
      
      # if( abs(f[rown,coln]-avg)/avg > 1.5){
      #   cat(rownames(f)[rown],colnames(f)[coln],"Avg","\n")
      # }
      # if(abs(f[rown,coln]-f[rown,1])/avg > 1.5){
      #   cat(rownames(f)[rown],colnames(f)[coln],"B","\n")      
      # }
      if(abs(f[rown,coln]-f[rown,2])/avg > 1.5){
        cat(rownames(f)[rown],colnames(f)[coln],"D",f[rown,1],f[rown,2],f[rown,3],f[rown,4],"\n")
      }
      # if(abs(f[rown,coln]-f[rown,3])/avg > 1.5){
      #   cat(rownames(f)[rown],colnames(f)[coln],"A","\n")      
      # }
      # if(abs(f[rown,coln]-f[rown,4])/avg > 1.5){
      #   cat(rownames(f)[rown],colnames(f)[coln],"C","\n")      
      # }
      coln <- coln + 1
    }
  }
  rown <- rown + 1
}

c0 <- subset(GFP,idents = "18")
table(c0$orig.ident)
df <- c0@assays$RNA@scale.data
dfa <- df[,grepl(paste("^",'A',"_",sep = ""),colnames(df))]
dfb <- df[,grepl(paste("^",'B',"_",sep = ""),colnames(df))]
dfc <- df[,grepl(paste("^",'C',"_",sep = ""),colnames(df))]
dfd <- df[,grepl(paste("^",'D',"_",sep = ""),colnames(df))]
# multiple ocrrections by bonfaronni any p value under 3.531522e-07 will be consider as significant.
# 0.05/(23857*6) = 3.531522e-07
genes <- rownames(c0)
i <- 1 
ab <- c()
ac <- c()
ad <- c()
bc <- c()
bd <- c()
cd <- c()
while ( i<= length(genes)) {
  
  # ccheck group A,B
  if (mean(dfa[i,]) != mean(dfb[i,])){
    r <- t.test(dfa[i,],dfb[i,])
    if(r$p.value <= 3.531522e-07){
      ab <- c(ab,genes[i])
      #print(r)
    }
  }
  
  # check group A,C
  if (mean(dfa[i,]) != mean(dfc[i,])){
    r <- t.test(dfa[i,],dfc[i,])
    if(r$p.value <= 3.531522e-07){
      ac <- c(ac,genes[i])
      #print(r)
    }
  }
  # check group A,D
  if (mean(dfa[i,]) != mean(dfd[i,])){
    r <- t.test(dfa[i,],dfd[i,])
    if(r$p.value <= 3.531522e-07){
      ad <- c(ad,genes[i])
      #print(r)
    }
  }
  # check group B,C
  if (mean(dfb[i,]) != mean(dfc[i,])){
    r <- t.test(dfb[i,],dfc[i,])
    if(r$p.value <= 3.531522e-07){
      bc <- c(bc,genes[i])
      #print(r)
    }
  }
  # check group B,D
  if (mean(dfb[i,]) != mean(dfd[i,])){
    r <- t.test(dfb[i,],dfd[i,])
    if(r$p.value <= 3.531522e-07){
      bd <- c(bd,genes[i])
      #print(r)
    }
  }
  # check group C,D  
  if (mean(dfc[i,]) != mean(dfd[i,])){
  r <- t.test(dfc[i,],dfd[i,])
  if(r$p.value <= 3.531522e-07){
    cd <- c(cd,genes[i])
    #print(r)
  }
}
  
  i<- i + 1
}

VlnPlot(c0,"ABRACL",group.by = "orig.ident")
intersect(ab,cd)
VlnPlot(c0,c("AMD1","KIF2A","MDM1","NEXN","PRKAR2B"),group.by = "orig.ident")




general <- c(0,1,3,5,6,7,8,9,10,11,12,15,16,17,18)
specific <- c(2,4,13,14)

#############################################################################
for (filenumber in general){
  
  
  c0 <- subset(GFP,idents = filenumber)
  table(c0$orig.ident)
  df <- c0@assays$RNA@scale.data
  dfa <- df[,grepl(paste("^",'A',"_",sep = ""),colnames(df))]
  dfb <- df[,grepl(paste("^",'B',"_",sep = ""),colnames(df))]
  dfc <- df[,grepl(paste("^",'C',"_",sep = ""),colnames(df))]
  dfd <- df[,grepl(paste("^",'D',"_",sep = ""),colnames(df))]
  # multiple ocrrections by bonfaronni any p value under 3.531522e-07 will be consider as significant.
  # 0.05/(23857*6) = 3.531522e-07
  genes <- rownames(c0)
  i <- 1 
  ab <- c()
  ac <- c()
  ad <- c()
  bc <- c()
  bd <- c()
  cd <- c()
  while ( i<= length(genes)) {
    
    # ccheck group A,B
    if (mean(dfa[i,]) != mean(dfb[i,])){
      r <- t.test(dfa[i,],dfb[i,])
      if(r$p.value <= 3.531522e-07){
        ab <- c(ab,genes[i])
        #print(r)
      }
    }
    
    # check group A,C
    if (mean(dfa[i,]) != mean(dfc[i,])){
      r <- t.test(dfa[i,],dfc[i,])
      if(r$p.value <= 3.531522e-07){
        ac <- c(ac,genes[i])
        #print(r)
      }
    }
    #check group A,D
    if (mean(dfa[i,]) != mean(dfd[i,])){
      r <- t.test(dfa[i,],dfd[i,])
      if(r$p.value <= 3.531522e-07){
        ad <- c(ad,genes[i])
        #print(r)
      }
    }
    # check group B,C
    if (mean(dfb[i,]) != mean(dfc[i,])){
      r <- t.test(dfb[i,],dfc[i,])
      if(r$p.value <= 3.531522e-07){
        bc <- c(bc,genes[i])
        #print(r)
      }
    }
    # check group B,D
    if (mean(dfb[i,]) != mean(dfd[i,])){
      r <- t.test(dfb[i,],dfd[i,])
      if(r$p.value <= 3.531522e-07){
        bd <- c(bd,genes[i])
        #print(r)
      }
    }
    # check group C,D
    if (mean(dfc[i,]) != mean(dfd[i,])){
      r <- t.test(dfc[i,],dfd[i,])
      if(r$p.value <= 3.531522e-07){
        cd <- c(cd,genes[i])
        #print(r)
      }
    }
    
    i<- i + 1
  }
  
  write(ab,paste("./Results/GeneDiscovery/c",filenumber,"_ab.txt",sep = ""))
  write(ac,paste("./Results/GeneDiscovery/c",filenumber,"_ac.txt",sep = ""))
  write(ad,paste("./Results/GeneDiscovery/c",filenumber,"_ad.txt",sep = ""))
  write(bc,paste("./Results/GeneDiscovery/c",filenumber,"_bc.txt",sep = ""))
  write(bd,paste("./Results/GeneDiscovery/c",filenumber,"_bd.txt",sep = ""))
  write(cd,paste("./Results/GeneDiscovery/c",filenumber,"_cd.txt",sep = ""))
  


}


VlnPlot(c0,ab[1:9],group.by = "orig.ident")
intersect(ab,cd)
