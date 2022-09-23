

#' This function do Clustering Analysis. Different distance methods are included, Euclidean, Viper Similairty and Correlation distance. Please choose the assay abd the slot to do the analysis.
#' @param  seu The seu expression object
#' res=0.6 for PA
#' @return  The seu results
#' @author Junqiang Wang
#' @export
SeuPAClustering<-function(seu, ndims=30, dims=10,
                          resolution=0.8,
                          assay=c("Activity","RNA"),
                          slot=c("data", "counts", "scaled.data"),
                          feature.sd.cutoff=0,
                         # top.n.features=1000,
                          distance=c("default", "viper_similarity", "cor"),
                         k.param = 20,
                         prune.SNN = 1/15,
                          nn=50 ){

  require(Seurat)
  require(ggplot2)




  if (is.null(seu@assays$Activity)==TRUE){

  message("Activity assay is NULL, add assay from counts")

  seu<-AddActivityAssay(seu, vp_mat=as.matrix(seu@assays$RNA@counts))

  }




  seu@active.assay <- 'Activity'


  message ("clustering using seurat default method...")

  seu@assays$Activity@scale.data<- as.matrix(seu@assays$Activity@data)


#tmp.sd<-apply(as.matrix(seu@assays$Activity@data), 1, sd)

#tmp.sd<-sort(tmp.sd, decreasing=TRUE)

# features<-names(tmp.sd[tmp.sd>feature.sd.cutoff])

 #features<-names(tmp.sd[1:top.n.features])


  features<-rownames(seu)

  seu <- RunPCA(seu, features = features)



  # set the k.param and prune.SNN parameter
  seu <- FindNeighbors(seu, k.param = k.param, prune.SNN = prune.SNN)

  seu <- FindClusters(seu, resolution = resolution)

  seu <- RunUMAP(seu, dims=1:dims)
  seu <- RunTSNE(seu, dims=1:dims)



  seu@meta.data$Clusters_PA_EuclideanDist<-as.factor(as.numeric(seu@meta.data$seurat_clusters))


  if(distance=="viper_similarity"){


  message ("Clustering using viper similarity as distance...")


  vpmat<-as.matrix(seu@assays$Activity@data)

  vp_dist <- as.dist(viperSimilarity(vpmat, nn=nn)) %>% as.matrix()

  seu@graphs <-FindNeighbors(vp_dist, distance.matrix=TRUE, k.param = k.param, prune.SNN = prune.SNN)

 # seu@graphs <-FindNeighbors(vp_dist, distance.matrix=TRUE, k.param = 20, prune.SNN = 1/15)

 # seu@graphs <-FindNeighbors(as.matrix(vp_dist), distance.matrix=TRUE, k.param = 50, prune.SNN = 1/15)

  #seu@graphs <-FindNeighbors(as.matrix(vp_dist), k.param = 50, prune.SNN = 1/15)

  seu <- FindClusters(seu, resolution =resolution, graph.name="snn")

  seu@meta.data$Clusters_PA_ViperDist<-as.factor(as.numeric(seu@meta.data$seurat_clusters))

  seu@meta.data$Clusters_PA_ViperDist<-paste0("C", seu@meta.data$Clusters_PA_ViperDist)

  }



  if(distance=="cor"){

    message ("clustering using viper similarity as distance...")

    vpmat<-as.matrix(seu@assays$Activity@data)

    vp_dist <- 1-cor(vpmat)

    seu@graphs <-FindNeighbors(vp_dist, distance.matrix=TRUE, k.param = k.param, prune.SNN = prune.SNN)

    #seu@graphs <-FindNeighbors(as.matrix(vp_dist), distance.matrix=TRUE, k.param = 20, prune.SNN = 1/15)

    #seu@graphs <-FindNeighbors(as.matrix(vp_dist), k.param = 50, prune.SNN = 1/15)

    seu <- FindClusters(seu, resolution =resolution, graph.name="snn")

    seu@meta.data$Clusters_PA_CorDist<-as.factor(as.numeric(seu@meta.data$seurat_clusters))

    seu@meta.data$Clusters_PA_CorDist<-paste0("C", seu@meta.data$Clusters_PA_CorDist)

  }


  seu@meta.data$seurat_clusters<- NULL

  return(seu)

}











SeuPAClusteringDefault<-function(seu, dims=50, resolution=0.6, assay=c("Activity","RNA"), slot=c("data", "counts", "scaled.data") ){

  require(Seurat)
  require(ggplot2)


  message ("clustering using seurat default method...")

  seu@active.assay <- 'Activity'

  seu@assays$Activity@scale.data<- as.matrix(seu@assays$Activity@data)

  seu<- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)


  #seu <- RunPCA(seu, features = rownames(seu))

  seu <- RunPCA(seu, features = VariableFeatures(seu))

  seu <- FindNeighbors(seu, dims=1:10)

  seu <- FindClusters(seu, resolution = resolution)


  seu <- RunUMAP(seu, dims=1:dims)
  seu <- RunTSNE(seu, dims=1:dims)


  seu@meta.data$PA.Clusters<-as.factor(as.numeric(seu@meta.data$seurat_clusters))

  return(seu)

}








