

#' This function generate a one-to-one name list including entrez and hgnc_symbol
#' @author Junqiang Wang
#' @export
NameConvertor<-function(){

tx2gene <- NA


if(is.na(tx2gene)){

  require(dplyr)
  require(biomaRt)

  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
  t2g <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id","ensembl_gene_id","entrezgene_id","external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, entrez = entrezgene_id)
  tx2gene <- t2g[,c(1,4)]
  tx2gene <- tx2gene[!is.na(tx2gene$entrez),]

  tx2gene<-tx2gene %>% distinct(entrez, .keep_all= TRUE)
  tx2gene<-tx2gene %>% distinct(hgnc_symbol, .keep_all= TRUE)

  #save(tx2gene,file = "tx2gene.rda")
}

return(tx2gene)

}



#' This function convert the human Entrez names to symbl
#' @param mat matrix, for example viper matrix
#' @return A matrix associated with gene symbol
#' @author  Junqiang Wang
#' @export
Entrez2symbolMat<-function(mat){

  require(dplyr)
  require(biomaRt)

  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
  t2g <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id","ensembl_gene_id","entrezgene_id","external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, entrez = entrezgene_id)
  tx2gene <- t2g[,c(1,4)]
  tx2gene <- tx2gene[!is.na(tx2gene$entrez),]


  tx2gene<-tx2gene %>% distinct(entrez, .keep_all= TRUE)
  tx2gene<-tx2gene %>% distinct(hgnc_symbol, .keep_all= TRUE)

  #----------------------------------

  rownames(mat)

  mat<-mat[which(match(rownames(mat), tx2gene$entrez)!= "NA"),]

  rownames(mat)<-tx2gene$hgnc_symbol[match(rownames(mat), tx2gene$entrez)]

  return(mat)

}






#' This function convert the human symbl to Entrez names
#' @param mat matrix, for example viper matrix
#' @return A matrix associated with gene symbol
#' @author  Junqiang Wang
#' @export
Symbol2EntrezMat<-function(mat){

  require(dplyr)
  require(biomaRt)

  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
  t2g <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id","ensembl_gene_id","entrezgene_id","external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, entrez = entrezgene_id)
  tx2gene <- t2g[,c(1,4)]
  tx2gene <- tx2gene[!is.na(tx2gene$entrez),]


  tx2gene<-tx2gene %>% distinct(entrez, .keep_all= TRUE)
  tx2gene<-tx2gene %>% distinct(hgnc_symbol, .keep_all= TRUE)

  #----------------------------------

  rownames(mat)

  mat<-mat[which(match(rownames(mat), tx2gene$hgnc_symbol)!= "NA"),]

  rownames(mat)<-tx2gene$entrez[match(rownames(mat), tx2gene$hgnc_symbol)]

  return(mat)

}











#' This function convert the human Entrez names to symbl
#' @param mat matrix, for example viper matrix
#' @return A matrix associated with gene symbol
#' @author  Junqiang Wang
#' @export
Ensembl2SymbolMat_Human<-function(mat){

  require(dplyr)
  require(biomaRt)

  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")

  ref<- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id","ensembl_gene_id","entrezgene_id","external_gene_name"), mart = mart)


  ref <- ref[,c(1,3)]

 ref<-ref %>%
     filter(hgnc_symbol!="") %>%
   distinct(ensembl_gene_id, .keep_all= TRUE) %>%
   distinct(hgnc_symbol, .keep_all= TRUE)

  rownames(mat)

  mat<-mat[which(match(rownames(mat), ref$ensembl_gene_id)!= "NA"),]

  rownames(mat)<-ref$hgnc_symbol[match(rownames(mat), ref$ensembl_gene_id)]

  return(mat)

}










#' This function convert the human Entrez names to symbl
#' @param mat matrix, for example viper matrix
#' @return A matrix associated with gene symbol
#' @author  Junqiang Wang
#' @export
Entrez2symbolIndex<-function(){

  require(dplyr)
  require(biomaRt)

  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host = "ensembl.org")
  t2g <- biomaRt::getBM(attributes = c("hgnc_symbol", "ensembl_transcript_id","ensembl_gene_id","entrezgene_id","external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, entrez = entrezgene_id)
  tx2gene <- t2g[,c(1,4)]
  tx2gene <- tx2gene[!is.na(tx2gene$entrez),]


  tx2gene<-tx2gene %>% distinct(entrez, .keep_all= TRUE)
  tx2gene<-tx2gene %>% distinct(hgnc_symbol, .keep_all= TRUE)

  #----------------------------------

 # rownames(mat)

  # mat<-mat[which(match(rownames(mat), tx2gene$entrez)!= "NA"),]

  # rownames(mat)<-tx2gene$hgnc_symbol[match(rownames(mat), tx2gene$entrez)]


  #return(mat)

  Entrez2symbolIndex<-tx2gene


  return(Entrez2symbolIndex)

}









































