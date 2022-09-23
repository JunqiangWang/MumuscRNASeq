



#' This function generate the differential signatures
#' @param seu seurat object
#' @param cluster.by the ID in seu@meata.data
#' @param reverse.clusters whether the order of the cluster id shoud be reversed
#' @param cluster.order the order of the cluster
#' @param groups.use the groups that are used run the presto
#' @param top.n return the top.n dif features
#' @param assay The assay that is used
#' @param slot the slot of the assay that is used
#' @param logFC.cutoff the cutoff of the presto
#' @return return the markers and store the markers in a new assay
#' @author  Junqiang Wang
#'@export
Presto1Signature<-function(seu, cluster.by="seurat_clusters",  signature.cluster= NULL, groups.use=NULL, reverse.clusters=c(TRUE, FALSE),  assay="RNA", slot=c("counts", "scale.data", "data")){

  require(dplyr)
  require(Seurat)
  require(presto)

seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]


if(slot=="scale.data"){

  tmp.mat<-seu@assays$RNA@scale.data

}else if(slot=="data"){

  tmp.mat<-seu@assays$RNA@data

}else if(slot=="counts"){

  tmp.mat<-seu@assays$RNA@counts
}


tmp.mat<-GetAssayData(object = seu[[assay]], slot = slot)


tmp.cluster.vector<-unname(unlist(subset(seu@meta.data, select=cluster.by)))

clusters<-sort(unique(tmp.cluster.vector))

if(is.null(groups.use)==TRUE){

markers<- presto::wilcoxauc(as.matrix(tmp.mat), tmp.cluster.vector)

} else if(is.null(groups.use)==FALSE){
  markers<- presto::wilcoxauc(as.matrix(tmp.mat), tmp.cluster.vector,  groups.use = groups.use)
  clusters<-groups.use
}




if(reverse.clusters==FALSE){

  clusters<-clusters

} else if (reverse.clusters==TRUE){

  clusters<-rev(clusters)

}


  markers %>%
    group_by(group) %>%
    arrange(desc(logFC)) -> all.markers


  signature<-all.markers%>%
    dplyr::filter(group==signature.cluster) %>%
    arrange(desc(logFC))

 return(signature)

}



#







#' This function generate the differential signatures
#' @param seu seurat object
#' @param cluster.by the ID in seu@meata.data
#' @param reverse.clusters whether the order of the cluster id shoud be reversed
#' @param cluster.order the order of the cluster
#' @param groups.use the groups that are used run the presto
#' @param top.n return the top.n dif features
#' @param assay The assay that is used
#' @param slot the slot of the assay that is used
#' @param logFC.cutoff the cutoff of the presto
#' @return return the markers and store the markers in a new assay
#' @author  Junqiang Wang
#'@export
PrestoMultiSignatures<-function(seu, cluster.by="seurat_clusters", groups.use=NULL, reverse.clusters=c(TRUE, FALSE),  assay="RNA", slot=c("counts", "scale.data", "data")){

  require(dplyr)
  require(Seurat)
  require(presto)

  seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]


  if(slot=="scale.data"){

    tmp.mat<-seu@assays$RNA@scale.data

  }else if(slot=="data"){

    tmp.mat<-seu@assays$RNA@data

  }else if(slot=="counts"){

    tmp.mat<-seu@assays$RNA@counts
  }


  tmp.mat<-GetAssayData(object = seu[[assay]], slot = slot)


  tmp.cluster.vector<-unname(unlist(subset(seu@meta.data, select=cluster.by)))

  clusters<-sort(unique(tmp.cluster.vector))

  if(is.null(groups.use)==TRUE){

    markers<- presto::wilcoxauc(as.matrix(tmp.mat), tmp.cluster.vector)

  } else if(is.null(groups.use)==FALSE){
    markers<- presto::wilcoxauc(as.matrix(tmp.mat), tmp.cluster.vector,  groups.use = groups.use)
    clusters<-groups.use
  }




  if(reverse.clusters==FALSE){

    clusters<-clusters

  } else if (reverse.clusters==TRUE){

    clusters<-rev(clusters)

  }


  markers %>%
    group_by(group) %>%
    arrange(desc(logFC)) -> all.markers


  signature<-all.markers%>%
    dplyr::filter(group==signature.cluster)

  return(signature)

}












