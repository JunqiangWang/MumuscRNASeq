


# reference
# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/



convertHumanGeneList2Mouse <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}



convertMouseGeneList2Human <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}





# This function convert the mouse genes to human
#' @export
GenesMouse2Human_Mat<- function(mat){

  x<-rownames(mat)
  require("biomaRt")

#
#
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

   # mouse = useEnsembl("ensembl",dataset="mmusculus_gene_ensembl", mirror = "uswest")
   # human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror="uswest")

  tx2gene = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

  tx2gene <- tx2gene[!is.na(tx2gene$HGNC.symbol),]

  tx2gene <- tx2gene[!is.na(tx2gene$MGI.symbol),]



  tx2gene<-tx2gene %>% distinct(MGI.symbol, .keep_all= TRUE)
  tx2gene<-tx2gene %>% distinct(HGNC.symbol, .keep_all= TRUE)

  #----------------------------------

  rownames(mat)

  mat<-mat[which(match(rownames(mat), tx2gene$MGI.symbol)!= "NA"),]

  rownames(mat)<-tx2gene$HGNC.symbol[match(rownames(mat), tx2gene$MGI.symbol)]

  return(mat)

}



# This function convert the human genes to mouse
#' @export
GenesHuman2Mouse_Mat<- function(mat){

  x<-rownames(mat)
  require("biomaRt")

  #human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  #mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")



  tx2gene = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", values = x , mart = human, attributesL = c("mgi_symbol"), martL = mouse, uniqueRows=T)

  tx2gene <- tx2gene[!is.na(tx2gene$HGNC.symbol),]

  tx2gene <- tx2gene[!is.na(tx2gene$MGI.symbol),]



  tx2gene<-tx2gene %>% distinct(MGI.symbol, .keep_all= TRUE)
  tx2gene<-tx2gene %>% distinct(HGNC.symbol, .keep_all= TRUE)

  #----------------------------------

  rownames(mat)

  mat<-mat[which(match(rownames(mat), tx2gene$HGNC.symbol)!= "NA"),]

  rownames(mat)<-tx2gene$MGI.symbol[match(rownames(mat), tx2gene$HGNC.symbol)]

  return(mat)

}









