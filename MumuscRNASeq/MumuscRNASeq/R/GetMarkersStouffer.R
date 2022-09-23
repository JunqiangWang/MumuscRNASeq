
#' This function generates the markers of the Top Differentially Activate Proteins by using the Stouffer integration method
#' @param seu seurat object
#' @param cluster.by the ID in seu@meata.data
#' @param reverse.clusters whether the order of the cluster id shoud be reversed
#' @param cluster.order the order of the cluster
#' @param top.n return the top.n dif features
#' @param assay The assay that is used
#' @param slot the slot of the assay that is used
#' @return return the markers and store the markers in a new assay
#' @author  Junqiang Wang
#' @export
GetMarkersStouffer<-function(seu, cluster.by="seurat_clusters",
                          reverse.clusters=c(FALSE, TRUE),
                          groups.use=NULL,
                          top.n=25,
                          marker.type=c("both", "top", "bottom"),
                          assay="Activity",
                          slot=c("data", "scale.data", "counts")){

    require(dplyr)

    require(Seurat)

    Idents(seu)<-cluster.by


    seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

    seu@meta.data$NES.Stouffer<-NA

    tmp.mat<-GetAssayData(object = seu[[assay]], slot = slot)

    tmp.cluster.vector<-unname(unlist(subset(seu@meta.data, select=cluster.by)))

    clusters<-sort(unique(tmp.cluster.vector))



    if(reverse.clusters==FALSE){

      clusters<-clusters

    } else if (reverse.clusters==TRUE){

      clusters<-rev(clusters)

    }



    all.cluster.markers<-vector()



    for (cluster.i in clusters){


      seu.i<-subset(seu, idents=cluster.i)

      tmp.mat.i<-GetAssayData(object = seu.i[[assay]], slot = slot)

      nes.stouffer.rows<-rowSums(tmp.mat.i)/sqrt(ncol(tmp.mat.i))

      names(nes.stouffer.rows)<-rownames(tmp.mat.i)

      nes.stouffer.rows<-sort(nes.stouffer.rows, decreasing = TRUE)

      cluster.markers.top<-nes.stouffer.rows%>%
        as.data.frame() %>%
        top_n(n = top.n) %>%
        rownames()


      nes.stouffer.rows<-sort(nes.stouffer.rows, decreasing = FALSE)

      cluster.markers.bottom<-nes.stouffer.rows%>%
        as.data.frame() %>%
        top_n(n = -top.n) %>%
        rownames()


      cluster.markers<-c(cluster.markers.top, cluster.markers.bottom)

      markers.id.top<-match(cluster.markers.top, rownames(tmp.mat))

      markers.id.bottom<-match(cluster.markers.bottom, rownames(tmp.mat))



      cells<-as.data.frame(seu@meta.data) %>%
        dplyr::filter(Clusters==cluster.i)

      cells.id<-match(rownames(cells), colnames(tmp.mat))

      mat<-tmp.mat[markers.id.top, cells.id ] %>% as.matrix()

      NES.Stouffer<-colSums(mat)/sqrt(ncol(mat))

      seu@meta.data$NES.Stouffer[match(names(NES.Stouffer), rownames(seu@meta.data))]<-NES.Stouffer

      all.cluster.markers<-c(all.cluster.markers, cluster.markers)




    }


    all.markers.id<-match(all.cluster.markers, rownames(tmp.mat))

    mat<-tmp.mat[all.markers.id, ] %>% as.matrix()


    # order the columns by group and nes,stouffer

    tmp<-seu@meta.data%>%
      as.data.frame() %>%
      arrange(Clusters, desc(NES.Stouffer)) %>%
      rownames()


    mat<-mat[, match(tmp, colnames(mat))]


    seu.markers<-CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

    seu.markers@meta.data<-seu@meta.data[match(colnames(mat), rownames(seu@meta.data)), ]

    seu.markers@assays$RNA@misc$diffMarkers<-all.cluster.markers

    seu.markers<-AddMarkerAssay(seu=seu.markers, marker.mat=mat)

    seu.markers@assays$Marker@misc$diffMarkers<-all.cluster.markers



    return(seu.markers)

}
























