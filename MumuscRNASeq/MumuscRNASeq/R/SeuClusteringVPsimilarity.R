


#' This function do clustering analysis using viper similarity as distance
#' @param  seu The seu object of protein activity
#' @return  The seurat object associated with cluster info
#' @author Junqiang Wang
#' @export
SeuClusteringVPsimilarity<-function(seu, nn=50,
                                           method=c("two.sided", "greater", "less"),
                                           my.dims=1:50,
                                           features = rownames(seu),
                                           k.param=20,
                                           prune.SNN=1/15,
                                           resolution=0.3,
                                           graph.name="snn"){

  require(Seurat)
  require(viper)

# vp has negative values, vst cannot be used here!

#seu<- FindVariableFeatures(seu, selection.method = "dispersion", nfeatures = 2000)

#get variable features

seu@assays$RNA@scale.data<-as.matrix(seu@assays$RNA@counts)

seu <- RunPCA(seu, features = features)

seu <- RunUMAP(seu, reduction = "pca", dims = my.dims)
seu <- RunTSNE(seu, reduction = "pca", dims = my.dims)

vpmat<-as.matrix(seu@assays$RNA@counts)

vp_dist <- as.dist(viperSimilarity(vpmat, nn=nn, method=method))

#seu@graphs <-FindNeighbors(as.matrix(vp_dist), k.param = k.param, prune.SNN = prune.SNN)


seu@graphs <-FindNeighbors(as.matrix(vp_dist), distance.matrix = TRUE, k.param = k.param, prune.SNN = prune.SNN)


seu <- FindClusters(seu, resolution =resolution, graph.name=graph.name)

seu@meta.data$seurat_clusters<-as.factor(as.numeric(seu@meta.data$seurat_clusters))

seu@meta.data$Clusters<-as.factor(as.numeric(seu@meta.data$seurat_clusters))

return(seu)

}







