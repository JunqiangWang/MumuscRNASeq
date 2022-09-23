


#' This function generate an accustomed Heatmap Plot of the Top Differentially Activate Proteins
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
#' @export
PrestoMarkers<-function(seu, cluster.by="seurat_clusters",
                        reverse.clusters=c(FALSE, TRUE),
                        groups.use=NULL,
                        top.n=25,
                        assay="Activity",
                        slot=c("data", "scale.data", "counts"),
                        pval.cutoff=1,
                        padj.cutoff=0.01,
                         logFC.cutoff=0.25){

  require(dplyr)

  require(Seurat)

  require(presto)


seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

seu@meta.data$NES.Stouffer<-NA

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
  dplyr::filter(logFC>logFC.cutoff)%>%
  dplyr::filter(padj<padj.cutoff)%>%
  dplyr::filter(pval<pval.cutoff)%>%

   arrange(desc(logFC)) %>%

  top_n(n = top.n, wt = logFC) -> all.markers

all.cluster.markers<-vector()



for (cluster.i in clusters){

 # cluster.i<-clusters[1]

  cluster.markers<-all.markers%>%
    dplyr::filter(group==cluster.i)

  markers.id<-match(cluster.markers$feature, rownames(tmp.mat))

  cells<-as.data.frame(seu@meta.data) %>%
    dplyr::filter(Clusters==cluster.i)

  cells.id<-match(rownames(cells), colnames(tmp.mat))

  mat<-tmp.mat[markers.id, cells.id ] %>% as.matrix()

  NES.Stouffer<-colSums(mat)/sqrt(ncol(mat))

  seu@meta.data$NES.Stouffer[match(names(NES.Stouffer), rownames(seu@meta.data))]<-NES.Stouffer

  all.cluster.markers<-c(all.cluster.markers, cluster.markers$feature)

}





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

seu.markers<-AddActivityAssay(seu=seu.markers, vp_mat=mat)

seu.markers@assays$Activity@misc$diffMarkers<-all.cluster.markers

#
#  seu.markers<-CreateSeuratObject(counts = mat, min.cells = 0, min.features = 0)
#
#  seu.markers@meta.data<-seu@meta.data
#
#
# NES.Stouffer<-seu.markers@meta.data$NES.Stouffer
#
# names(NES.Stouffer)<-rownames(seu.markers@meta.data)
#
# NES.Stouffer<-sort(-NES.Stouffer)
#
# id<-match(names(NES.Stouffer), rownames(seu.markers@meta.data))
#
# seu.markers@meta.data<-seu.markers@meta.data[id,]
#
#
# seu.markers@assays$RNA@counts<-seu.markers@assays$RNA@counts[,id]
#
# seu.markers@assays$RNA@misc$diffMarkers<-all.cluster.markers
#
# # Add assay to seurat markers
#
# dat<-as.matrix(seu.markers@assays$RNA@counts)
#
# seu.markers<-AddActivityAssay(seu=seu.markers, vp_mat=dat)
#
#
# seu.markers@assays$Activity@misc$diffMarkers<-all.cluster.markers

return(seu.markers)

}



















