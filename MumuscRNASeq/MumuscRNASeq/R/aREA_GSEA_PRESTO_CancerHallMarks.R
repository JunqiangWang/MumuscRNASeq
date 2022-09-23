




#' this function do gsea enrichment analysis
#' @author Junqing Wang
#' @export
aREA_GSEA <- function(current.geneset,current.sig,uniform.wts = TRUE,tcm.size = NULL){


  require(ggplot2)
  require(viper)


  current.regul <- list(tfmode = sign(current.geneset), likelihood = as.numeric(abs(current.geneset)))



  if(!is.null(tcm.size)){
    subset.targets <- current.regul$tfmode*current.regul$likelihood
    subset.targets <- c(sort(subset.targets,decreasing = TRUE)[seq(from = 1, to = round(tcm.size/2), by = 1)],rev(sort(subset.targets,decreasing = FALSE)[seq(from = 1, to = round(tcm.size/2), by = 1)]))
    # modify the geneset if desired
    if(uniform.wts){
      subset.targets <- sign(subset.targets)
    }
    new.regul <- current.regul
    new.regul$tfmode <- sign(subset.targets)
    new.regul$likelihood <- as.numeric(abs(subset.targets))
    current.int <- list(new.regul)
    class(current.int) <- "regulon"
    current.geneset <- new.regul$tfmode
  } else {
    subset.targets <- current.regul$tfmode*current.regul$likelihood
    # modify the geneset if desired
    if(uniform.wts){
      subset.targets <- sign(subset.targets)
    }
    new.regul <- current.regul
    new.regul$tfmode <- sign(subset.targets)
    new.regul$likelihood <- as.numeric(abs(subset.targets))
    current.int <- list(new.regul)
    class(current.int) <- "regulon"
    current.geneset <- new.regul$tfmode
  }

  current.enrich.area <- aREA(eset = current.sig, regulon = current.int)

  area.es.value <- signif(current.enrich.area$es,3)
  area.nes.value <- signif(current.enrich.area$nes,3)

  #area.log10.pvalue <- (log(2) + pnorm(q = abs(current.enrich.area$nes), lower.tail = FALSE, log.p = TRUE))/log(10)

  area.pvalue <-pnorm(q = abs(current.enrich.area$nes), lower.tail = FALSE, log.p = FALSE)

  area.log10.pvalue<-log10(area.pvalue)*(-1)

  area.pvalue<-signif(area.pvalue, 3)


  #area.pvalue.power <- floor(area.log10.pvalue)
 # area.pvalue.num <- signif((10^(area.log10.pvalue - floor(area.log10.pvalue))),3)


  out<-list("aREA.ES"=area.es.value, "aREA.NES"=area.nes.value,"aREA.P"=area.pvalue, "aREA.log10.P"=area.log10.pvalue)

  return(out)

}






#' this function do gsea enrichment analysis
#' @author Junqing Wang
#' @export
aREA_PRESTO_CancerHallMarkers<-function(seu,
                                        cluster.by="seurat_clusters",
                                        species=c("human", "mouse"),
                                        signature.cluster= NULL,
                                        groups.use=NULL,
                                        reverse.clusters=c(TRUE, FALSE),
                                        assay="RNA",
                                        slot=c("scale.data","counts", "data")){

   require(presto)
   require(tidyverse)


  seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

  nc<-length(unique(seu@meta.data$Clusters))

  clusters<-as.vector(unique(seu@meta.data$Clusters))

#nc<-length(unique(seu@meta.data$seurat_clusters))
#clusters<-as.vector(unique(seu@meta.data$seurat_clusters))



presto.out <- wilcoxauc(X=seu, group_by=cluster.by, assay=slot, seurat_assay = assay)

head(presto.out)
# we have all the genes for each cluster
# BiocManager::install(c("msigdbr", "fgsea"))
library(msigdbr)
library(fgsea)
library(dplyr)
library(ggplot2)

msigdbr_show_species()
m_df<- msigdbr(species = species, category = "H")

#m_df<- msigdbr(species = "human", category = "H")

head(m_df)
gene_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)


p<-list()
nes<-matrix(0,50,nc)
aREA.p<-matrix(0,50,nc)

colnames(nes)<-clusters
colnames(aREA.p)<-clusters

for (j in 1:nc){

  #j<-1

  cj<-clusters[j]
  print (cj)
  # select only the feature and auc columns for fgsea, which statistics to use is an open question
  cluster.genes<- presto.out %>%
    dplyr::filter(group == cj) %>%
    arrange(desc(logFC)) %>%
    dplyr::select(feature, logFC)
  ranks<- deframe(cluster.genes)

  head(ranks)



  for (i in 1:length(gene_sets)){

  #  i<-1

  tmp.id<-match(gene_sets[[i]], names(ranks))

  tmp.id<-tmp.id[!is.na(tmp.id)]

  tmp.gene.set<-ranks[tmp.id]

  tmp.gene.set[tmp.gene.set<1e-100]<-1

  aREA.gsea.ij<- aREA_GSEA(current.geneset=tmp.gene.set,current.sig=ranks,uniform.wts = TRUE, tcm.size = NULL)

  print( aREA.gsea.ij)

  nes[i,j]<-aREA.gsea.ij$aREA.NES

  aREA.p[i,j]<-aREA.gsea.ij$aREA.P

  }
}

rownames(nes)<-names(gene_sets)
rownames(aREA.p)<-names(gene_sets)


aREA.p.fdr<-matrix(p.adjust(as.vector(aREA.p), method="fdr"), ncol=ncol(aREA.p))

colnames(aREA.p.fdr)<-colnames(aREA.p)

rownames(aREA.p.fdr)<-rownames(aREA.p)



out<-list("nes"=nes, "aREA.p"=aREA.p, "aREA.p.fdr"=aREA.p.fdr)

return(out)

}





#' This function do heatmap plot of the gsea analysis result
#'@author Junqiang
#'@export
PlotHeatmapGSEA<-function(mat, col_fun_range=c(-10, 0, 10)){

  require(ComplexHeatmap)

# log10.FDR.signed<-aREA.out$aREA.p.fdr*sign(aREA.out$nes)

# mat<-log10.FDR.signed

# mat<-as.data.frame(mat) %>% arrange(desc(SF8628_Radiation_1w)) %>% as.matrix()

#q<-quantile(mat, c(0.01, 0.99))

#brewer.pal(11,"RdYlBu")

colors = rev(c("#D53E4F", "white","#3288BD"))

#col_fun = circlize::colorRamp2(c(q[1], 0, q[2]), colors)


col_fun = circlize::colorRamp2(col_fun_range, colors)


p<-Heatmap(mat, name = "log10FDR_Signed",

           cell_fun = function(j, i, x, y, width, height, fill) {
             grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 10))
           },

           cluster_columns = FALSE,
           show_column_dend = FALSE,
           cluster_column_slices = TRUE,
           column_title_gp = gpar(fontsize = 7),
           column_gap = unit(0.5, "mm"),
           cluster_rows = FALSE,
           show_row_dend = TRUE,
           col = col_fun,

           column_names_rot = 60,

           row_names_gp = gpar(fontsize = 7),
           column_title_rot = 0,
           #  top_annotation = HeatmapAnnotation(foo = anno_block(gp = gpar(fill = scales::hue_pal()(9)))),
           show_column_names = TRUE,
           use_raster = FALSE,
           # raster_device = 'png',
           raster_quality = 4)

p

}






