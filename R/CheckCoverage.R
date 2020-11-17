#' CheckCellCoverage
#'@description This function take 2 Seurat3 Object as input, and return the intersection of two object by Cell ID.
#'
#'
#' CheckCellCoverage()
CheckCellCoverage <- function(so1,so2){

  so1cells <- as.data.frame(so1@active.ident)
  so2cells <- as.data.frame(so2@active.ident)

  cells1 <- row.names(so1cells)
  cells2 <- row.names(so2cells)

  temp <- Reduce(intersect,list(v1 = cells1,
                        v2 = cells2))

  print(paste("Seurat Obeject 1 has: ",length(cells1)," cells"))
  print(paste((length(temp)/length(so1@active.ident)*100),"%"))

  print(paste("Seurat Obeject 2 has: ",length(cells2)," cells"))
  print(paste( (length(temp)/length(so2@active.ident)*100),"%"))

  print(paste("Intersection: ",length(temp)))
  return(as.data.frame( temp ))
}
