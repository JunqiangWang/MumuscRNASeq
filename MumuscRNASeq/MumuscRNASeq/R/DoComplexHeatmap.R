

#' Accustomed Heatmap Plot of the Top Differentially Activate Proteins
#' @param seu the seurat markers
#' @param cluster.by the column names indicating the column that is used in the analysis
#' @auther Junqiang Wang
#' @export
DoComplexHeatmap<-function(seu,
                           cluster.by="seurat_clusters",
                           column_split_option=c("YES", "NO"),
                           assay="Activity",
                           slot=c("data", "scale.data", "data"),
                           reverse.clusters=c(FALSE, TRUE),
                           ann.top.df.columns=c("seurat_clusters", "Phase"),
                           data.type=c("expression", "activity"),
                           color_quantile=0.01,
                           fontsize.column.title.gp=7,
                           fontsize.row.names.gp=7
                           ){

  require(ComplexHeatmap)

  require(Seurat)

  set.seed(77)

 tmp.cluster.vector<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

 ann<-seu@meta.data

# get the mat

mat<-GetAssayData(object = seu[[assay]], slot = slot)

mat<-as.matrix(mat)

#mat<-as.matrix(seu@assays$RNA@counts)

rownames(mat)<-seu@assays$Activity@misc$diffMarkers

column_split<-tmp.cluster.vector


if(reverse.clusters==FALSE){
  column_split_factor<-factor(column_split)
} else if (reverse.clusters==TRUE){
  column_split_factor<-factor(column_split, levels=rev(levels(factor(column_split))))
}


ann.top.df<-subset(ann,  select= ann.top.df.columns)

ann.top<-HeatmapAnnotation(df=ann.top.df,
                           name=ann.top.df.columns,
                           annotation_name_rot = 0,
                           col=list(Phase=c("G1"="purple", "G2M"="yellow", "S"="green"))
                           )

q<-quantile(mat, c(color_quantile, 1-color_quantile))

#brewer.pal(11,"RdYlBu")
#colors = rev(c("#D53E4F", "white","#3288BD"))

colors = rev(c("red", "white","blue"))

col_fun = circlize::colorRamp2(c(q[1], 0, q[2]), colors)



if(column_split_option=="YES"){

if(data.type=="expression"){

  set.seed(77)

p<-Heatmap(mat, name = "Expression",

        top_annotation = ann.top,
        column_split = column_split_factor,
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,

        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun,

        column_title_gp = gpar(fontsize = fontsize.column.title.gp),

        row_names_gp = gpar(fontsize = fontsize.row.names.gp),
        column_title_rot = 0,
        show_column_names = FALSE,
        use_raster = FALSE,
       # raster_device = 'png',
        raster_quality = 4)

}else{

          set.seed(77)
          p<-Heatmap(mat, name = "Protein Activity",

                     top_annotation = ann.top,
                     column_split = column_split_factor,
                     cluster_columns = FALSE,
                     show_column_dend = FALSE,
                     cluster_column_slices = TRUE,

                     column_gap = unit(0.5, "mm"),
                     cluster_rows = FALSE,
                     show_row_dend = FALSE,
                       col = col_fun,

                     column_title_gp = gpar(fontsize = fontsize.column.title.gp),

                     row_names_gp = gpar(fontsize = fontsize.row.names.gp),
                     column_title_rot = 0,
                     show_column_names = FALSE,
                     use_raster = FALSE,
                     # raster_device = 'png',
                     raster_quality = 4)

        }

return(p)


}else  {

  if(data.type=="expression"){

    set.seed(77)

    p<-Heatmap(mat, name = "Expression",

               top_annotation = ann.top,
               #column_split = column_split_factor,
               cluster_columns = FALSE,
               show_column_dend = FALSE,
               cluster_column_slices = TRUE,

               column_gap = unit(0.5, "mm"),
               cluster_rows = FALSE,
               show_row_dend = FALSE,
               col = col_fun,
               column_title_gp = gpar(fontsize = fontsize.column.title.gp),

               row_names_gp = gpar(fontsize = fontsize.row.names.gp),
               column_title_rot = 0,
               show_column_names = FALSE,
               use_raster = FALSE,
               # raster_device = 'png',
               raster_quality = 4)
  } else{

    set.seed(77)
    p<-Heatmap(mat, name = "Protein Activity",

               top_annotation = ann.top,
               #column_split = column_split_factor,
               cluster_columns = FALSE,
               show_column_dend = FALSE,
               cluster_column_slices = TRUE,

               column_gap = unit(0.5, "mm"),
               cluster_rows = FALSE,
               show_row_dend = FALSE,
               col = col_fun,
               column_title_gp = gpar(fontsize = fontsize.column.title.gp),

               row_names_gp = gpar(fontsize = fontsize.row.names.gp),

               column_title_rot = 0,
               show_column_names = FALSE,
               use_raster = FALSE,
               # raster_device = 'png',
               raster_quality = 4)

  }

  return(p)
}
}


















