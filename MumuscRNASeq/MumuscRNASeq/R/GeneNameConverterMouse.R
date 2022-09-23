

#' this function convert the mouse gene names
#' @param mat The matrix with the rownames of
#' @author Junqing Wang
#' @export
Ens2SymblMat_Mouse<- function(mat){

  x<-rownames(mat)

  require(dplyr)

  require(biomaRt)

  #mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  #tx2gene = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "mmusculus_gene_ensembl", host = "ensembl.org")

  t2g <- biomaRt::getBM(attributes = c("mgi_symbol", "ensembl_transcript_id","ensembl_gene_id","entrezgene_id","external_gene_name"), mart = mart)

 # t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id, ens_gene = ensembl_gene_id, ext_gene = external_gene_name, entrez = entrezgene_id)

  tx2gene <- t2g[,c(3,5)]

  tx2gene <- tx2gene[!is.na(tx2gene$external_gene_name),]

  tx2gene<-tx2gene %>% distinct(ensembl_gene_id, .keep_all= TRUE)

  tx2gene<-tx2gene %>% distinct(external_gene_name, .keep_all= TRUE)

  rownames(mat)

  mat<-mat[which(match(rownames(mat), tx2gene$ensembl_gene_id)!= "NA"),]

  rownames(mat)<-tx2gene$external_gene_name[match(rownames(mat), tx2gene$ensembl_gene_id)]

  return(mat)

}
















