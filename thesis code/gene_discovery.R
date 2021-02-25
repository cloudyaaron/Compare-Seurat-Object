library("Seurat")
library(datasets)
library("futile.matrix")

jane <- readRDS("./analysis_jane/RemapedAggregate_clustered.Robj")
gdbg <- AverageExpression(jane,add.ident = "GroupID")
g <-  gdbg$RNA
write.csv(g,file = "./analysis_jane/raw_average_expression(not_scaled).csv")
table(Idents(jane),jane$GroupID)


genes <- rownames(g)

maxn <- 19
i <- 0
colindex <- colnames(g)
while(i < 19){
  i2 <- 1
  df <- data.frame(matrix(nrow = length(genes),ncol = 0))
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
  write.csv(df,file = paste("./analysis_jane/(no_scale)cluster ",i,".csv",sep = ""))
  i <- i + 1
  print("")
}



general <- c(0,1,2,4,5,6,7,8,9,10,12,13,15,18)
specific <- c(3,11,14,16)
sp2 <- c(17)

#############################################################################
for (filenumber in specific){
  
  #a c1 static #b c3 shear #c c2 static+sbchir #d c4 Shear+ sbchir
  c0 <- subset(jane,idents = filenumber)
  table(c0$GroupID)
  df <- GetAssayData(c0,"data")
  dfa <- df[,grepl(paste("^",'C1',"_",sep = ""),colnames(df))]
  dfb <- df[,grepl(paste("^",'C3',"_",sep = ""),colnames(df))]
  dfc <- df[,grepl(paste("^",'C2',"_",sep = ""),colnames(df))]
  dfd <- df[,grepl(paste("^",'C4',"_",sep = ""),colnames(df))]
  # multiple ocrrections by bonfaronni any p value under 3.531522e-07 will be consider as significant.
  # 0.05/(23597*6) = 3.531522e-07
  th = 0.05/(23597*6)
  genes <- rownames(df)
  i <- 1 
  ab <- c()
  ac <- c()
  ad <- c()
  bc <- c()
  bd <- c()
  cd <- c()
  abp <- c()
  acp <- c()
  adp <- c()
  bcp <- c()
  bdp <- c()
  cdp <- c()
  while ( i<= length(rownames(df))) {
    
    # ccheck group A,B
    if (mean(dfa[i,]) != mean(dfb[i,])){
      r <- t.test(dfa[i,],dfb[i,])
      if(r$p.value <= th){
        ab <- c(ab,genes[i])
        abp <- c(abp,r$p.value)
        #print(r)
      }
    }

    # check group A,C
    if (mean(dfa[i,]) != mean(dfc[i,])){
      r <- t.test(dfa[i,],dfc[i,])
      if(r$p.value <=th){
        ac <- c(ac,genes[i])
        acp <- c(acp,r$p.value)

        #print(r)
      }
    }
    #check group A,D
    if (mean(dfa[i,]) != mean(dfd[i,])){
      r <- t.test(dfa[i,],dfd[i,])
      if(r$p.value <= th){
        ad <- c(ad,genes[i])
        adp <- c(adp,r$p.value)

        #print(r)
      }
    }
    # check group B,C
    if (mean(dfb[i,]) != mean(dfc[i,])){
      r <- t.test(dfb[i,],dfc[i,])
      if(r$p.value <= th){
        bc <- c(bc,genes[i])
        bcp <- c(bcp,r$p.value)

        #print(r)
      }
    }
    # check group B,D
    if (mean(dfb[i,]) != mean(dfd[i,])){
      r <- t.test(dfb[i,],dfd[i,])
      if(r$p.value <= th){
        bd <- c(bd,genes[i])
        bdp <- c(bdp,r$p.value)

        #print(r)
      }
    }
    # check group C,D
    if (mean(dfc[i,]) != mean(dfd[i,])){
      r <- t.test(dfc[i,],dfd[i,])
      if(r$p.value <= th){
        cd <- c(cd,genes[i])
        cdp <- c(cdp,r$p.value)

        #print(r)
      }
    }
    
    i<- i + 1
  }
  #a c1 static #b c3 shear #c c2 static+sbchir #d c4 Shear+ sbchir
  
  temp <- do.call(rbind,Map(data.frame,Gene=ab,Pvalue=abp))
  write.csv(temp,paste("./analysis_jane/c",filenumber,"_staticVSshear.csv",sep = ""))
  temp <- do.call(rbind,Map(data.frame,Gene=ac,Pvalue=acp))
  write.csv(temp,paste("./analysis_jane/c",filenumber,"_staticVSstaticSBchir.csv",sep = ""))
  temp <- do.call(rbind,Map(data.frame,Gene=ad,Pvalue=adp))
  write.csv(temp,paste("./analysis_jane/c",filenumber,"_staticVSshearSBchir.csv",sep = ""))
  temp <- do.call(rbind,Map(data.frame,Gene=bc,Pvalue=bcp))
  write.csv(temp,paste("./analysis_jane/c",filenumber,"_shearVSstaticSBchir.csv",sep = ""))
  temp <- do.call(rbind,Map(data.frame,Gene=bd,Pvalue=bdp))
  write.csv(temp,paste("./analysis_jane/c",filenumber,"_shearVSshearSBchir.csv",sep = ""))
  temp <- do.call(rbind,Map(data.frame,Gene=cd,Pvalue=cdp))
  write.csv(temp,paste("./analysis_jane/c",filenumber,"_staticSBchirVSshearSBchir.csv",sep = ""))
  
}
