


#' This function generate an accustomed Heatmap Plot of the Top Differentially Activate Proteins
#' @param seu seurat object
#' @param cluster.by the ID in seu@meata.data
#' @param reverse.clusters whether the order of the cluster id shoud be reversed
#' @param cluster.order the order of the cluster
#' @param groups.use the groups that are used run the presto
#' @param top.n return the top.n dif features
#' @param assay The assay that is used
#' @param slot the slot of the assay that is used
#' @return return the markers and store the markers in a new assay
#' @author  Junqiang Wang
#' @export
GetActivityMarkers<-function(seu,
                               cluster.by="seurat_clusters",
                        reverse.clusters=c(FALSE, TRUE),
                        order.cluster=NULL,
                        ident.query=NULL,
                        ident.ref=NULL,
                        top.n=25,
                        assay="Activity",
                        slot=c("data", "scale.data","counts"),
                     only.pos = TRUE,
                     test.use=c("roc","wilcox"),
                     min.pct.threshold = 0.25,
                     myAUC.cutoff=0.5,
                     avg_logfc.threshold = 0,
                     avg_diff.threshold=0,
                     p_val_adj.threshold=0.01){

  require(dplyr)

  require(Seurat)

  DefaultAssay(object = seu) <- assay

   Idents(seu)<-cluster.by

seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

seu@meta.data$NES.Stouffer<-NA

tmp.mat<-GetAssayData(object = seu[[assay]], slot = slot)

tmp.cluster.vector<-unname(unlist(subset(seu@meta.data, select=cluster.by)))

clusters<-sort(unique(tmp.cluster.vector))


if(is.null(ident.query)==TRUE){


  markers <-FindAllMarkers(seu,  test.use = test.use, only.pos = TRUE, slot=slot, assay =assay, min.pct = min.pct.threshold)


} else if(is.null(ident.query)==FALSE){

 markers <- FindMarkers(seu, ident.1=ident.query, ident.2=ident.ref, slot=slot, assay=assay, test.use = test.use, only.pos = TRUE, min.pct = min.pct.threshold)

  clusters<-c(ident.query, ident.ref)

}


print(markers)


 if (reverse.clusters==TRUE){

  clusters<-rev(clusters)

 }


if(is.null(order.cluster)==FALSE){

clusters<-order.cluster

}





if(is.null(ident.query)==TRUE){

if(test.use =="wilcox"){

  markers %>%
    group_by(cluster) %>%
    filter(p_val_adj<p_val_adj.threshold) %>%
    filter(avg_log2FC> avg_logfc.threshold) %>%
    top_n(n = top.n, wt=-p_val_adj) -> top.n.markers

  top.n.markers<-top.n.markers %>%
    as.data.frame() %>%
    arrange(cluster, p_val_adj)

# remove duplicate

top.n.markers %>%
arrange(p_val_adj) %>%
  distinct(gene, .keep_all = TRUE) %>%
arrange(cluster)->top.n.markers



}else if (test.use =="roc") {

  markers %>%
    group_by(cluster) %>%
    filter(myAUC>myAUC.cutoff) %>%
    # filter(p_val_adj<p_val_adj.threshold) %>%
    filter(avg_diff>avg_diff.threshold) %>%
    top_n(n = top.n, wt=avg_diff) -> top.n.markers

  top.n.markers<-top.n.markers %>%
    as.data.frame() %>%
 arrange(cluster, -avg_diff)

  # remove duplicate

  top.n.markers %>%
    arrange(p_val_adj) %>%
    distinct(gene, .keep_all = TRUE) %>%
    arrange(cluster)->top.n.markers



}
}






if(is.null(ident.query)==FALSE){

  if(test.use =="wilcox"){

    markers %>%
     # group_by(cluster) %>%
      filter(p_val_adj<p_val_adj.threshold) %>%
      filter(avg_log2FC>avg_logfc.threshold) %>%
      top_n(n = top.n, wt=-p_val_adj) -> top.n.markers
     # top_n(n = top.n, wt=-avg_log2FC) -> top.n.markers

  }else if (test.use =="roc") {

    markers %>% as.data.frame() %>%
      arrange(desc(avg_diff)) %>%
      group_by(cluster) %>%

      filter(myAUC>myAUC.cutoff) %>%
      # filter(p_val_adj<p_val_adj.threshold) %>%
      filter(avg_diff>avg_diff.threshold) %>%
      top_n(n = top.n, wt=avg_diff) -> top.n.markers

    top.n.markers<-as.data.frame(top.n.markers)

  }
}



message("Here are the top markers:")

print(top.n.markers)





all.cluster.markers<-vector()





if(is.null(ident.query)==TRUE){

for (cluster.i in clusters){

  cluster.markers<-top.n.markers%>%
    dplyr::filter(cluster==cluster.i)

  markers.id<-match(cluster.markers$gene, rownames(tmp.mat))

  cells<-as.data.frame(seu@meta.data) %>%
    dplyr::filter(Clusters==cluster.i)

  cells.id<-match(rownames(cells), colnames(tmp.mat))

  mat<-tmp.mat[markers.id, cells.id ] %>% as.matrix()

  NES.Stouffer<-colSums(mat)/sqrt(ncol(mat))

  seu@meta.data$NES.Stouffer[match(names(NES.Stouffer), rownames(seu@meta.data))]<-NES.Stouffer

  all.cluster.markers<-c(all.cluster.markers, cluster.markers$gene)

}
}




if(is.null(ident.query)==FALSE){

    cluster.markers<-top.n.markers

    markers.id<-match(rownames(cluster.markers), rownames(tmp.mat))

    cells<-as.data.frame(seu@meta.data)

    cells.id<-match(rownames(cells), colnames(tmp.mat))

    mat<-tmp.mat[markers.id, cells.id ] %>% as.matrix()

    NES.Stouffer<-colSums(mat)/sqrt(ncol(mat))

    seu@meta.data$NES.Stouffer[match(names(NES.Stouffer), rownames(seu@meta.data))]<-NES.Stouffer

    all.cluster.markers<- rownames(cluster.markers)

  }





message("Have got all cluster markers...")

print(all.cluster.markers)

all.markers.id<-match(all.cluster.markers, rownames(tmp.mat))

 mat<-tmp.mat[all.markers.id, ] %>% as.matrix()




# order the columns by group and nes,stouffer

tmp<-seu@meta.data%>%
  as.data.frame() %>%
  arrange(Clusters, -NES.Stouffer) %>%
  rownames()


mat<-mat[, match(tmp, colnames(mat))]


seu.markers<-CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)

seu.markers@meta.data<-seu@meta.data[match(colnames(mat), rownames(seu@meta.data)), ]

seu.markers@assays$RNA@misc$diffMarkers<-all.cluster.markers

seu.markers<-AddMarkerAssay(seu=seu.markers, marker.mat=mat)

seu.markers@assays$Marker@misc$diffMarkers<-all.cluster.markers


return(seu.markers)

}















