#' Make a confusion matrix
#' @description This function take 2 seurat object as input, return a confusion matrix
#' @return a dataframe contain the confusion matrix
#' MakeConfusionMatrix()
MakeConfusionMatrix <- function(so1,so2){
  library(plyr)
  library(ggplot2)
  if(length(so1@reductions) == 0 || length(so2@reductions) == 0){
    warning("One of the Seurat Object has not been cluster yet")
    return()
  }
  intersect_cell <- as.array( length(CheckCellCoverage(so1,so2)$temp))

  if((intersect_cell/length(so1@active.ident))<=0.5 ||(intersect_cell/length(so2@active.ident))<=0.5){
    s <- readline(prompt = "Cell coverage is less than 50%, Forced start?(Y|N)")

  }else{

    s <- readline("Start making ConfusionMatrix? (Y|N)")

  }
  if(toupper(s) == "Y"){
    print("Start making confusion matrix....")

    # remember to change back to so1 and so2
    df1 <-as.data.frame(so1@active.ident)
    df2 <-as.data.frame(so2@active.ident)

    df1$cell <- rownames(df1)
    df2$cell <- rownames(df2)


    initial_length_df1 <- sort(as.array( unique( df1[,1])))
    initial_length_df2 <- sort(as.array( unique( df2[,1])))
    confusion_matrix <- data.frame(matrix(0,nrow = length(initial_length_df1) ,ncol = length(initial_length_df2) ))
    colnames(confusion_matrix)<-initial_length_df2
    rownames(confusion_matrix)<- initial_length_df1

    df <- data.frame(matrix(nrow = 0,ncol = 3))
    colnames(df) <- c("Cell_id","so1_ident","so2_ident")
    for (i in 1:nrow(df1)) {
      if(df1$cell[i]%in% df2$cell){
        index <- which(grepl(df1$cell[i],df2$cell))
        newrow <- data.frame(cell_id=df1$cell[i],so1_ident=df1[i,1],so2_ident=df2[index,1])
        df <- rbind(df,newrow)
      }
    }

    # use DF to create a 2d table
    for (i in 1:nrow(df)) {
      confusion_matrix[as.numeric( df$so1_ident[i]),as.numeric(df$so2_ident[i])] <- confusion_matrix[as.numeric(df$so1_ident[i]),as.numeric(df$so2_ident[i])]+1
    }
    heatmap(data.matrix(confusion_matrix),ylab = paste("1st data",so1@project.name) ,xlab = paste("2nd data",so2@project.name))
    return(confusion_matrix)
  }else{

    return()
  }

  # if(){
  #
  # }

}
