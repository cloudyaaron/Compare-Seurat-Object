#' compare two cluster between two different seurat object
#' @description This function take 2 seurat object, and their corresponding clusternumber as input, return a set of genes that make the difference
#' @return a list of genes
#' ClusterDifference()

ClusterDifference <- function(so1,c1,so2,c2){
  library(Seurat)
  d1 <- subset(so1,idents = c1)
  d2 <- subset(so2,idents = c2)
  df1 <- GetAssayData(d1,"data")
  df2 <- GetAssayData(d2,"data")
  df1 <- df1[order(row.names(df1)),]
  df2 <- df2[order(row.names(df2)),]
  genes1 <- rownames(df1)
  print("sample 1 has")
  cat(length(genes1))
  print("genes")
  genes2 <- rownames(df2)
  print("sample 2 has")
  cat(length(genes2))
  print("genes")
  ana_gene <- intersect(genes1,genes2)
  print("only analyze:")
  print(length(ana_gene))
  print(" numbers of genes")
  returngenes <- c()
  plist<-c()
  i<- 1
  threshold <-0.05/(length(ana_gene)*2)
  for (g in ana_gene) {
    if (range(df1[g,])[1]!=range(df1[g,])[2] || range(df2[g,])[1]!=range(df2[g,])[2] )    {
      r <- t.test(df1[g,],df2[g,])
      if(r$p.value <= threshold){
        returngenes <- c(returngenes,g)
        plist <- c(plist,r$p.value)
        print(r)
      }
    }

  }
  dfr <- do.call(rbind, Map(data.frame, Gene=returngenes, p.value=plist))
  return(dfr)
}
