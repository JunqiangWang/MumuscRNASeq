

#' This function do Stouffer integration of the PA in each cluster by rows
#' @param  seu The seurat object
#'@cluster.by The colname of the meta used as clustering info
#'@return Integreated NES
#' @export
StoufferIntegrateSeuClusterRows<-function(seu,
                                         # assay=c("Activity", "RNA"),
                                         # slot=c("data", "scale.data", "counts"),

                                          assay="Activity",
                                          slot="data",
                                          top_n=NULL,
                                          cluster.by=NULL,
                                           order.cluster=NULL,
                                          remove.duplicates=TRUE,
                                         z.cutoff=NULL,
                                         p.cutoff=NULL,
                                         output.format=c("NES", "log10.p.adj"),
                                         log10.p.adj.cutoff=NULL
                                         ){


  require(Seurat)
  require(dplyr)


  Idents(seu) <- cluster.by


  if(cluster.by=="seurat_clusters"){

    clusters<-paste0("C", seu@meta.data$seurat_clusters)

    seu@meta.data$tmp<-clusters

    Idents(seu) <- "tmp"

    cluster.by<-"tmp"

  }


  if(cluster.by=="PA.seurat_clusters"){

    clusters<-paste0("C", seu@meta.data$PA.seurat_clusters)

    seu@meta.data$tmp<-clusters

    Idents(seu) <- "tmp"

    cluster.by<-"tmp"

  }



  tmp<-match(cluster.by, colnames(seu@meta.data))

  clusters<-unique(seu@meta.data[,tmp])

  clusters<-sort(as.vector(clusters))


 if(is.null(order.cluster)==FALSE){clusters<-order.cluster}



   nes.stouffer<-vector()


for (i in clusters){

  seu.i <- subset(seu, idents =i)

   mat<-GetAssayData(object = seu.i[[assay]], slot = slot) %>% as.matrix()

  #mat<-as.matrix(seu.i@assays$Activity@data)

  nes.stouffer.i<-rowSums(mat)/sqrt(ncol(mat))

  nes.stouffer<-cbind(nes.stouffer, nes.stouffer.i)

}

   nes.stouffer<-as.matrix(nes.stouffer)

rownames(nes.stouffer)<-rownames(seu)
colnames(nes.stouffer)<-clusters







#--------------------------
# big loop
#---------------------------

if(output.format=="NES"){


  if(is.null(p.cutoff)==FALSE){

 # convert to p_values and adjust the value

    nes<-as.matrix(nes.stouffer)

    pval<-pnorm(nes, lower.tail = FALSE)

    #adjust the p-value by FDR

    unlist(pval)

    pval.fdr<-matrix(as.numeric(p.adjust(unlist(pval), method="fdr")), nrow=nrow(pval))

    rownames(pval.fdr)<-rownames(pval)

    colnames(pval.fdr)<-colnames(pval)


    tmp<-apply(nes, 1, min)

    nes.stouffer<-nes.stouffer[which(tmp<p.cutoff), ]

  }


if(is.null(top_n)==TRUE){

return(nes.stouffer)

}else{


markers.stouffer<-nes.stouffer

top.n.markers<-vector()

cluster.names<-colnames(markers.stouffer)


if(is.null(order.cluster)==FALSE){

cluster.names<-order.cluster

}


for (i in cluster.names){

  print(i)

  nes.i<-markers.stouffer[,i]

  nes.i<-sort(nes.i, decreasing=TRUE)

  mrs<-sort(nes.i, decreasing=TRUE)[1:top_n]  %>% as.data.frame()


  mrs<-cbind(mrs, i)

  colnames(mrs)<-c("Nes.stoufferIntegrated", "clusters")

  mrs$gene<-rownames(mrs)

  top.n.markers<-rbind(top.n.markers, mrs)

}


# remove the duplicate genes

if(remove.duplicates==TRUE && is.null(order.cluster)==TRUE){

top.n.markers %>%
  arrange(-Nes.stoufferIntegrated) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(clusters, -Nes.stoufferIntegrated)->top.n.markers

} else if(remove.duplicates==TRUE && is.null(order.cluster)==FALSE){

  top.n.markers %>%
    arrange(-Nes.stoufferIntegrated) %>%
    distinct(gene, .keep_all = TRUE) ->top.n.markers


 tmp<-vector()

for (i in order.cluster){

tmp.i<-top.n.markers %>%
  dplyr::filter(clusters==i)

tmp<-rbind(tmp, tmp.i)

}

 top.n.markers<-tmp


}

# z scores


tmp<-match(top.n.markers$gene, rownames(markers.stouffer))

markers.stouffer<-markers.stouffer[tmp, ]

return(markers.stouffer)

}





}else if (output.format=="log10.p.adj"){

#----------------------------------
# adj output format
#----------------------------------


  # convert to p values and adjust the p values

  nes<-as.matrix(nes.stouffer)

  pval<-pnorm(nes, lower.tail = FALSE)

  unlist(pval)

  pval.fdr<-matrix(as.numeric(p.adjust(unlist(pval), method="fdr")), nrow=nrow(pval))

  rownames(pval.fdr)<-rownames(pval)

  colnames(pval.fdr)<-colnames(pval)

  #

  pval.fdr<--log10(pval.fdr)

  pval.fdr<-as.data.frame(pval.fdr)

  tmp<-apply(pval.fdr,2, as.numeric)
  rownames(tmp)<-rownames(pval.fdr)

  pval.fdr<-tmp

  #pval.fdr<-pval.fdr[,-1]

  mat.p<-as.data.frame(pval.fdr)


#------------------ p cutoff


    if(is.null(log10.p.adj.cutoff)==FALSE){

      tmp<-apply(mat.p, 1, max)

      mat.p<-mat.p[which(tmp>log10.p.adj.cutoff), ]

    }




#----convert nes to log adj p values, end


  if(is.null(top_n)==TRUE){

    return(mat.p)

  }else{


    markers.stouffer<-mat.p

    top.n.markers<-vector()

    cluster.names<-colnames(markers.stouffer)


    if(is.null(order.cluster)==FALSE){

      cluster.names<-order.cluster

    }


    for (i in cluster.names){

      print(i)

      nes.i<-markers.stouffer[,i]

      nes.i<-sort(nes.i, decreasing=TRUE)

      mrs<-sort(nes.i, decreasing=TRUE)[1:top_n]  %>% as.data.frame()


      mrs<-cbind(mrs, i)

      colnames(mrs)<-c("Nes.stoufferIntegrated", "clusters")

      mrs$gene<-rownames(mrs)

      top.n.markers<-rbind(top.n.markers, mrs)

    }


    # remove the duplicate genes

    if(remove.duplicates==TRUE && is.null(order.cluster)==TRUE){

      top.n.markers %>%
        arrange(-Nes.stoufferIntegrated) %>%
        distinct(gene, .keep_all = TRUE) %>%
        arrange(clusters, -Nes.stoufferIntegrated)->top.n.markers

    } else if(remove.duplicates==TRUE && is.null(order.cluster)==FALSE){

      top.n.markers %>%
        arrange(-Nes.stoufferIntegrated) %>%
        distinct(gene, .keep_all = TRUE) ->top.n.markers


      tmp<-vector()

      for (i in order.cluster){

        tmp.i<-top.n.markers %>%
          dplyr::filter(clusters==i)

        tmp<-rbind(tmp, tmp.i)

      }

      top.n.markers<-tmp


    }

    # z scores


    tmp<-match(top.n.markers$gene, rownames(markers.stouffer))

    markers.stouffer<-markers.stouffer[tmp, ]

    return(markers.stouffer)

}

}

}








#' This function do Stouffer integration of the PA in 1 cluster by rows
#' @param  seu The seurat object
#'@cluster.by The colname of the meta used as clustering info
#'@return Integreated NES
#' @export
StoufferIntegrateSeu1ClusterRows<-function(seu, nn=25,
                                           # assay=c("Activity", "RNA"),
                                           # slot=c("data", "scale.data", "counts"),

                                           assay="Activity",
                                           slot="data",

                                           cluster.by="PA.seurat_clusters",

                                           sub_clusters="Null"


                                           ){



  require(Seurat)
  require(dplyr)


  Idents(seu) <- cluster.by


  if(cluster.by=="seurat_clusters"){

    clusters<-paste0("C", seu@meta.data$seurat_clusters)

    seu@meta.data$tmp<-clusters

    Idents(seu) <- "tmp"

    cluster.by<-"tmp"

  }


  if(cluster.by=="PA.seurat_clusters"){

    clusters<-paste0("C", seu@meta.data$PA.seurat_clusters)

    seu@meta.data$tmp<-clusters

    Idents(seu) <- "tmp"

    cluster.by<-"tmp"

  }



  tmp<-match(cluster.by, colnames(seu@meta.data))

  clusters<-unique(seu@meta.data[,tmp])

  clusters<-sort(as.vector(clusters))



  nes.stouffer<-vector()


    seu <- subset(seu, idents =sub_clusters)

    mat<-GetAssayData(object = seu[[assay]], slot = slot) %>% as.matrix()

    #mat<-as.matrix(seu.i@assays$Activity@data)


    NES.Stouffer<-rowSums(mat)/sqrt(ncol(mat))

    names(NES.Stouffer)<-rownames(mat)

    NES.Stouffer<-sort(NES.Stouffer, decreasing=TRUE)

    MRs.postive<-sort(NES.Stouffer, decreasing=TRUE)[1:nn]


    MRs.negative<-sort(NES.Stouffer, decreasing=FALSE)[1:nn]

    MRs<-names(c(MRs.postive, MRs.negative))

    tmp<-match(MRs, rownames(mat))

    mat.MR<-mat[tmp, ]

    mat.MR.postive<-mat[match(names(MRs.postive), rownames(mat)), ]


    NES.Stouffer.2<-colSums(mat.MR.postive)/sqrt(nrow(mat.MR.postive))

    NES.Stouffer.2<-sort(NES.Stouffer.2, decreasing=TRUE)

    tmp<-match(names(NES.Stouffer.2), colnames(mat.MR))

    mat.MR<-mat.MR[,tmp]



    seu.markers<-CreateSeuratObject(counts = mat.MR, min.cells = 0, min.features = 0)

    seu.markers@meta.data<-seu@meta.data[match(colnames(mat.MR), rownames(seu@meta.data)), ]

    seu.markers@assays$RNA@misc$diffMarkers<-rownames(mat.MR)

    seu.markers<-AddMarkerAssay(seu=seu.markers, marker.mat=mat.MR)

    seu.markers@assays$Marker@misc$diffMarkers<-rownames(mat.MR)


    return(seu.markers)



}









#' This function do Stouffer integration and identify the markers
#' @param  seu The seurat object
#'@cluster.by The colname of the meta used as clustering info
#'@return Integreated NES
#' @aliases StoufferIntegrateSeuClusterMarkers
#' @export
GetActivityMarkersStouffer<-function(seu, top.n=25,
                                           # assay=c("Activity", "RNA"),
                                           # slot=c("data", "scale.data", "counts"),
                                           assay="Activity",
                                           slot="data",
                                           output.format=c("seurat", "matrix"),
                                           reverse.clusters=FALSE,
                                         order.cluster=NULL,
                                           cluster.by=NULL){


  require(dplyr)

  require(Seurat)

  DefaultAssay(object = seu) <- assay

  Idents(seu)<-cluster.by

  seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

  seu@meta.data$NES.Stouffer<-NA

  tmp.mat<-GetAssayData(object = seu[[assay]], slot = slot)

  tmp.cluster.vector<-unname(unlist(subset(seu@meta.data, select=cluster.by)))

  clusters<-as.vector(sort(unique(tmp.cluster.vector)))



message("clusters are found:")

print(clusters)




  if (reverse.clusters==TRUE){

    clusters<-rev(clusters)

  }


if(is.null(order.cluster)==FALSE){

  clusters<-order.cluster

}





# Get the markers of each cluster seperately

top.n.markers<-vector()


for (i in clusters){

  print(i)

  seu.i <- subset(seu, idents =i)

  mat<-GetAssayData(object = seu.i[[assay]], slot = slot) %>% as.matrix()


  NES.Stouffer<-rowSums(mat)/sqrt(ncol(mat))

  names(NES.Stouffer)<-rownames(mat)

  NES.Stouffer<-sort(NES.Stouffer, decreasing=TRUE)

  mrs<-sort(NES.Stouffer, decreasing=TRUE)[1:top.n] %>% as.data.frame()


 mrs<-cbind(mrs, i)

colnames(mrs)<-c("Nes.stoufferIntegrated", "clusters")

mrs$gene<-rownames(mrs)

top.n.markers<-rbind(top.n.markers, mrs)

}



# remove duplciate markers across multiple clusters


top.n.markers %>%
  arrange(-Nes.stoufferIntegrated) %>%
  distinct(gene, .keep_all = TRUE) %>%
  arrange(clusters, -Nes.stoufferIntegrated)->top.n.markers




  message("Here are the top markers:")

  print(top.n.markers)


  if (output.format=="matrix"){

    return (top.n.markers)

  }else{


  all.cluster.markers<-vector()


    for (cluster.i in clusters){

      cluster.markers<-top.n.markers%>%
        dplyr::filter(clusters==cluster.i)

      markers.id<-match(cluster.markers$gene, rownames(tmp.mat))

      cells<-as.data.frame(seu@meta.data) %>%
        dplyr::filter(Clusters==cluster.i)

      cells.id<-match(rownames(cells), colnames(tmp.mat))

      mat<-tmp.mat[markers.id, cells.id ] %>% as.matrix()

      NES.Stouffer<-colSums(mat)/sqrt(ncol(mat))

      seu@meta.data$NES.Stouffer[match(names(NES.Stouffer), rownames(seu@meta.data))]<-NES.Stouffer

      all.cluster.markers<-c(all.cluster.markers, cluster.markers$gene)

    }





  message("Have got all cluster markers...")

  print(all.cluster.markers)

  all.markers.id<-match(all.cluster.markers, rownames(tmp.mat))

  mat<-tmp.mat[all.markers.id, ] %>% as.matrix()


if(is.null(seu@assays$RNA@data)==FALSE){

tmp.exp<-seu@assays$RNA@data %>% as.matrix()

mat.exp<-tmp.exp[all.markers.id, ] %>% as.matrix()

}


  if(is.null(seu@assays$RNA@scale.data)==FALSE){

    tmp.exp.scaled<-seu@assays$RNA@scale.data %>% as.matrix()

    all.markers.id<-match(all.cluster.markers, rownames(tmp.exp.scaled))

    mat.exp.scaled<-tmp.exp.scaled[all.markers.id, ] %>% as.matrix()

  }


  if(is.null(seu@assays$RNA@counts)==FALSE){

    tmp.counts<-seu@assays$RNA@counts %>% as.matrix()

    all.markers.id<-match(all.cluster.markers, rownames(tmp.counts))

    mat.counts<-tmp.counts[all.markers.id, ] %>% as.matrix()

  }



  # order the columns by group and nes,stouffer

  tmp<-seu@meta.data%>%
    as.data.frame() %>%
    arrange(Clusters, -NES.Stouffer) %>%
    rownames()


  mat<-mat[, match(tmp, colnames(mat))]

  mat.exp<-mat.exp[, match(tmp, colnames(mat))]

  mat.exp.scaled<-mat.exp.scaled[, match(tmp, colnames(mat.exp.scaled))]

 mat.counts<-mat.counts[, match(tmp, colnames(mat.counts))]


  seu.markers<-CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

  seu.markers@meta.data<-seu@meta.data[match(colnames(mat), rownames(seu@meta.data)), ]

  seu.markers@assays$RNA@misc$diffMarkers<-all.cluster.markers

  seu.markers<-AddMarkerAssay(seu=seu.markers, marker.mat=mat)

  seu.markers@assays$Marker@misc$diffMarkers<-all.cluster.markers



  seu.markers@assays$RNA@counts<-mat.counts


  seu.markers@assays$RNA@data<-mat.exp

  seu.markers@assays$RNA@scale.data<-mat.exp.scaled
#


  return(seu.markers)
}


}























