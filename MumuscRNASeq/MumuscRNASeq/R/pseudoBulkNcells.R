

#' This function create the pseudobulk using N cells
#' @param mat matrix of single cells
#' @param N number of cells
#' @export
PseudoBulkNcells<-function(mat, 
                           N.cells=10,
                           N.pseudoBulk=10,
                           suffix.pseudoBulk=NULL,
                           replace=TRUE){
  
  require(dplyr)
  
  mat<-mat %>% as.matrix()

  out<-vector()
  
  for(i in 1:N.pseudoBulk){
  
  id.sampled.i<-sample(ncol(mat), size=N.cells, replace = replace, prob = NULL)
  
  pseudo.bulk.i<-mat[, id.sampled.i]
  
  pseudo.bulk.i<-rowSums(pseudo.bulk.i)
  
  out<-cbind(out, pseudo.bulk.i)
  
  }
  
  if(!is.null(suffix.pseudoBulk)){colnames(out)<-paste0( "suffix.pseudoBulk", "pseudo.bulk.",  1:N.pseudoBulk)}else{
  
  colnames(out)<-paste0("pseudo.bulk.",  1:N.pseudoBulk)
  
    }
  
  return(out)
  
}
