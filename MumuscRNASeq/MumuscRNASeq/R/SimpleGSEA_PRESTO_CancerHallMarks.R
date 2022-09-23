

#' this function do gsea enrichment analysis using fixed reference
#' @param ref.cluster The reference cluster
#' @author Junqing Wang
#' @export
SimpleGSEA_PRESTO_CancerHallMarkers_Ref<-function(seu,
                                                  cluster.by="seurat_clusters",
                                                  species=c("human", "mouse"),
                                                  signature.cluster= NULL,
                                                  groups.use=NULL,
                                                  ref.cluster=NULL,
                                                  reverse.clusters=c(TRUE, FALSE),
                                                  assay="RNA",
                                                  slot=c("scale.data","counts", "data")){

  require(presto)
  require(tidyverse)

  # we have all the genes for each cluster
  # BiocManager::install(c("msigdbr", "fgsea"))

  library(msigdbr)
  library(fgsea)
  library(dplyr)
  library(ggplot2)

  #--------------
  # get the gene sets

  msigdbr_show_species()

  m_df<- msigdbr(species = species, category = "H")

  #m_df<- msigdbr(species = "human", category = "H")

  head(m_df)

  gene_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)

  #----------------------------

  seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

  nc<-length(unique(seu@meta.data$Clusters))

  clusters<-as.vector(unique(seu@meta.data$Clusters))

  clusters<-sort(setdiff(clusters, ref.cluster))

  p<-list()

  nes<-matrix(0,50,length(clusters))

  p.value<-matrix(0,50, length(clusters))

  colnames(nes)<-paste0("C", clusters)

  colnames(p.value)<-paste0("C", clusters)

  for (j in 1:length(clusters)){

    #j<-1

    cj<-clusters[j]

    print (cj)

    presto.out <- wilcoxauc(X=seu, group_by=cluster.by, groups_use=c(cj, ref.cluster),  assay=slot, seurat_assay = assay)

    head(presto.out)

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

      gsea.out<- gsea(signature=ranks, geneset=tmp.gene.set, score=1, twoTails=FALSE, pout=FALSE, per=1000)


      nes[i,j]<-gsea.out$nes

      p.value[i,j]<-gsea.out$p.value

    }
  }

  rownames(nes)<-names(gene_sets)

  rownames(p.value)<-names(gene_sets)

  p.value.fdr<-matrix(p.adjust(as.vector(p.value), method="fdr"), ncol=ncol(p.value))

  colnames(p.value.fdr)<-colnames(p.value)

  rownames(p.value.fdr)<-rownames(p.value)

  out<-list("nes"=nes, "p.value"=p.value, "p.value.fdr"=p.value.fdr)

  return(out)

}


#' this function do gsea enrichment analysis
#' @author Junqing Wang
#' @export
SimpleGSEA_PRESTO_CancerHallMarkers<-function(seu,
                                        cluster.by="seurat_clusters",
                                        species=c("human", "mouse"),
                                        signature.cluster= NULL,
                                        groups.use=NULL,
                                        reverse.clusters=c(TRUE, FALSE),
                                        assay="RNA",
                                        slot=c("scale.data","counts", "data")){

   require(presto)
   require(tidyverse)

  # we have all the genes for each cluster
  # BiocManager::install(c("msigdbr", "fgsea"))
  library(msigdbr)
  library(fgsea)
  library(dplyr)
  library(ggplot2)


  seu@meta.data$Clusters<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

  nc<-length(unique(seu@meta.data$Clusters))

  clusters<-as.vector(unique(seu@meta.data$Clusters))

presto.out <- wilcoxauc(X=seu, group_by=cluster.by, groups_use=groups.use,  assay=slot, seurat_assay = assay)

head(presto.out)


msigdbr_show_species()
m_df<- msigdbr(species = species, category = "H")

#m_df<- msigdbr(species = "human", category = "H")

head(m_df)
gene_sets<- m_df %>% split(x = .$gene_symbol, f = .$gs_name)


p<-list()
nes<-matrix(0,50,nc)
p.value<-matrix(0,50,nc)

colnames(nes)<-clusters
colnames(p.value)<-clusters

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

  gsea.out<- gsea(signature=ranks, geneset=tmp.gene.set, score=1, twoTails=FALSE, pout=FALSE, per=1000)


  nes[i,j]<-gsea.out$nes

  p.value[i,j]<-gsea.out$p.value

  }
}

rownames(nes)<-names(gene_sets)
rownames(p.value)<-names(gene_sets)


p.value.fdr<-matrix(p.adjust(as.vector(p.value), method="fdr"), ncol=ncol(p.value))

colnames(p.value.fdr)<-colnames(p.value)

rownames(p.value.fdr)<-rownames(p.value)



out<-list("nes"=nes, "p.value"=p.value, "p.value.fdr"=p.value.fdr)

return(out)

}






#' This function do heatmap plot of the gsea analysis result
#'@author Junqiang
#'@export
PlotHeatmapGSEA<-function(mat, col_fun_range=c(-10, 0, 10)){

# log10.FDR.signed<-aREA.out$p.value.fdr*sign(aREA.out$nes)

# mat<-log10.FDR.signed

# mat<-as.data.frame(mat) %>% arrange(desc(SF8628_Radiation_1w)) %>% as.matrix()

#q<-quantile(mat, c(0.01, 0.99))

#brewer.pal(11,"RdYlBu")

colors = rev(c("#D53E4F", "white","#3288BD"))

#col_fun = circlize::colorRamp2(c(q[1], 0, q[2]), colors)


col_fun = circlize::colorRamp2(col_fun_range, colors)


p<-Heatmap(mat, name = "+/- -log10FDR",

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






