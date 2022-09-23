

#' This function run the seurat do SCTransform
#' @param  seu The seu expression object
#' @species The species
#' @return  The seu results
#' @author Junqiang Wang
#' @export
SeuExpClustering<-function(seu, species=c("Human","Mouse"),  ndims=50, dims=50, resolution=0.8){

  require(Seurat)
  require(sctransform)
  require(ggplot2)

  seu@active.assay<-'RNA'



  message("Procressing expresion data sample using the standarded Seurat Pipeline!")

  if(species=="Human"){

    seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")

  }else if(species=="Mouse"){

    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^mt-")

    #seu <- PercentageFeatureSet(seu, pattern = "^Mt-", col.name = "percent.mt")

  }else{

    stop("Please select the species: species:Human or Mouse!")

  }

  seu <- CellCycleScoring(seu,g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
  seu@meta.data$Phase <- factor(seu@meta.data$Phase,levels=c("G1","S","G2M"))
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score

  # VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0, ncol = 3)

  # store mitochondrial percentage in object meta data

  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

  seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)

  seu <- ScaleData(seu, features = rownames(seu))

  seu <- RunPCA(seu, features = VariableFeatures(seu))


  ElbowPlot(seu, ndims = ndims)

  seu <- FindNeighbors(seu, dims=1:10)

  seu <- FindClusters(seu, resolution = resolution)


  seu <- RunUMAP(seu, dims=1:dims)
  seu <- RunTSNE(seu, dims=1:dims)


  seu@meta.data$Clusters_expression<-paste0("C", as.factor(as.numeric(seu@meta.data$seurat_clusters)))

  return(seu)

}





















