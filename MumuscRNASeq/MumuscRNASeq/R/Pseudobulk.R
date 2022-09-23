

#' PseudobulkCountGen
#' @description  generate the Pseudobulk counts from scRNA-Seq counts
#' @param  seu The seurat object
#' @param cluster.by the colnames of the seu@meta.data
#' @return A count matrix
#' @export
seuPseudoBulk<-function(seu, cluster.by){


  tmp<-match(cluster.by, colnames(seu@meta.data))
  clusters<-as.vector(sort(unique(seu@meta.data[,tmp])))


  out<-matrix(0, nrow(seu@assays$RNA@counts), length(clusters))

  rownames(out)<-rownames(seu@assays$RNA@counts)
  colnames(out)<-clusters

  for (j in 1:length(clusters)){

    reg_pattern <- FetchData(object = seu, vars = cluster.by)

    seu.i<-seu[, which(x=reg_pattern==clusters[j])]

    exp<-as.matrix(seu.i@assays$RNA@counts)

    exp<-rowSums(exp)

    # convert to RPKM
    #exp<-log2(1e6*exp/sum(exp)+1)

    out[,j]<-exp

  }

  out<-round(out, digits=2)

  return(out)

}






















