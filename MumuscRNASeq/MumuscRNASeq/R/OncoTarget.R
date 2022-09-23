


#' This function returns the OncoTarget genes
#'
#'
#' @export
OncoTargetGenes<-function(){

  require(drugbank)
  require(n1database)
  data(drugbase)
  pos <- sapply(names(dbase),function(x) dbase[[x]][['oncology']] & !dbase[[x]][['topical']])
  dbase2 <- dbase[pos]
  targets <- unique(unlist(sapply(names(dbase2),function(x) names(dbase2[[x]][['targets']]),USE.NAMES = FALSE)))

  names<-NameConvertor()

  tmp<-match(targets, names$entrez)

  OncoTarget<-names$hgnc_symbol[tmp]

  OncoTarget<-OncoTarget[is.na(OncoTarget)==FALSE]

  return(OncoTarget)

}







#' This funtion filter the matrix and get the OncoTarget mat
#' @export
OncoTarget<-function(vp.mat){

require(drugbank)
require(n1database)
data(drugbase)
pos <- sapply(names(dbase),function(x) dbase[[x]][['oncology']] & !dbase[[x]][['topical']])
dbase2 <- dbase[pos]
targets <- unique(unlist(sapply(names(dbase2),function(x) names(dbase2[[x]][['targets']]),USE.NAMES = FALSE)))

names<-NameConvertor()

tmp<-match(targets, names$entrez)

OncoTarget<-names$hgnc_symbol[tmp]

OncoTarget<-OncoTarget[is.na(OncoTarget)==FALSE]

targets<-OncoTarget

tmp<-intersect(targets, rownames(vp.mat))
tmp<-match(tmp,  rownames(vp.mat))

vp.mat<-vp.mat[tmp, ]

return(vp.mat)

}







