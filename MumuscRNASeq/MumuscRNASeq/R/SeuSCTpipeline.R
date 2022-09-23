

#' This function run the seurat do SCTransform
#' @param  seu The seu expression object
#' @param method The method used for SCTransform
#' @param vars.to.regress The variables used in SCTransform, for example, "percent.mt", "nCount_RNA"
#' @return  The seu results
#' @author Junqiang Wang
#' @export
SeuSCTpipeline<-function(seu, species=c("Human","Mouse"), method= c("default", "glmGamPoi"), vars.to.regress = NULL, ndims=50, dims=50, resolution=0.8){

  require(Seurat)
  require(sctransform)
  require(ggplot2)

  message("Procressing expresion data sample using Seurat Pipeline with SCT transformation!")

if(species=="Human"){

  seu <- PercentageFeatureSet(seu, pattern = "^MT-", col.name = "percent.mt")

}else if(species=="Mouse"){

  seu <- PercentageFeatureSet(seu, pattern = "^Mt-", col.name = "percent.mt")

}else{

  stop("Please select the species: species:Human or Mouse!")

}

  # VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),pt.size = 0, ncol = 3)

  # store mitochondrial percentage in object meta data

if (method=="default"){

  seu <- SCTransform(seu, vars.to.regress = vars.to.regress, return.only.var.genes=FALSE)

} else if (method=="glmGamPoi") {

  seu <- SCTransform(seu, method = "glmGamPoi", vars.to.regress = vars.to.regress, return.only.var.genes=FALSE)

} else{

  stop("Please select one method used in SCtransform!")

}


  seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)

  seu <- ScaleData(pbmc, features = rownames(seu))

  seu <- RunPCA(seu, features = VariableFeatures(seu))
  seu <- CellCycleScoring(seu,g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes, set.ident = TRUE)
  seu@meta.data$Phase <- factor(seu@meta.data$Phase,levels=c("G1","S","G2M"))
  seu$CC.Difference <- seu$S.Score - seu$G2M.Score


  ElbowPlot(seu, ndims = ndims)
  seu <- FindNeighbors(seu, dims = 1:dims)
  seu <- FindClusters(seu, resolution = resolution)
  seu <- RunUMAP(seu, dims=1:dims)
  seu <- RunTSNE(seu, dims=1:dims)


  seu@meta.data$seurat_clusters<-as.factor(as.numeric(seu@meta.data$seurat_clusters))

  return(seu)

}


































