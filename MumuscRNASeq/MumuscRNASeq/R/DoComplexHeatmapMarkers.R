


#' Accustomed Heatmap Plot of the Top Differentially Activate Proteins
#' @param seu the seurat markers
#' @param cluster.by the column names indicating the column that is used in the analysis
#' @auther Junqiang Wang
#' @alises DoComplexHeatmapMarkers
#' @export
PlotSeuMarkers<-function(seu,
                           cluster.by="seurat_clusters",
                           column_split_option=c("YES", "NO"),
                           assay="Marker",
                           slot=c("data", "scale.data", "count"),
                           reverse.clusters=c(FALSE, TRUE),
                           ann.top.df.columns=c("seurat_clusters", "Phase"),
                           legend.name=c("Expression", "Activity"),
                           color_quantile=0.01,
                           fontsize.column.title.gp=7,
                           fontsize.row.names.gp=7,
                         column_names_rot = 0,
                         column_title_rot = 0,

                         order.cluster=NULL,
                          rowAnnotation=NULL
                           ){

  require(ComplexHeatmap)

  require(Seurat)

  set.seed(77)

 tmp.cluster.vector<-seu@meta.data[,match(cluster.by, colnames(seu@meta.data))]

 ann<-seu@meta.data

# get the mat

#mat<-GetAssayData(object = seu[["Marker"]], slot = "data")

 mat<-GetAssayData(object = seu[[assay]], slot = slot)

mat<-as.matrix(mat)

rownames(mat)<-seu@assays$Marker@misc$diffMarkers

column_split<-tmp.cluster.vector




if(is.null(order.cluster)==TRUE && reverse.clusters==FALSE){
  column_split_factor<-factor(column_split)
} else if (is.null(order.cluster)==TRUE && reverse.clusters==TRUE){
  column_split_factor<-factor(column_split, levels=rev(levels(factor(column_split))))
} else if (is.null(order.cluster)==FALSE){

  column_split_factor<-factor(column_split, levels=order.cluster)
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

  set.seed(77)

p<-Heatmap(mat, name = legend.name,

        top_annotation = ann.top,
        column_split = column_split_factor,
        cluster_columns = FALSE,
        show_column_dend = FALSE,
        cluster_column_slices = TRUE,

        column_gap = unit(0.5, "mm"),
        cluster_rows = FALSE,
        show_row_dend = FALSE,
        col = col_fun,
        column_names_rot=column_names_rot,
        column_title_gp = gpar(fontsize = fontsize.column.title.gp),

        row_names_gp = gpar(fontsize = fontsize.row.names.gp),
        column_title_rot = column_title_rot,
        show_column_names = FALSE,
        use_raster = FALSE,
       # raster_device = 'png',
        raster_quality = 4)


}else  {

    set.seed(77)

    p<-Heatmap(mat, name = legend.name,

               top_annotation = ann.top,
               #column_split = column_split_factor,
               cluster_columns = FALSE,
               show_column_dend = FALSE,
               cluster_column_slices = TRUE,

               column_gap = unit(0.5, "mm"),
               cluster_rows = FALSE,
               show_row_dend = FALSE,
               col = col_fun,
               column_names_rot=column_names_rot,

               column_title_gp = gpar(fontsize = fontsize.column.title.gp),

               row_names_gp = gpar(fontsize = fontsize.row.names.gp),
               column_title_rot = column_title_rot,
               show_column_names = FALSE,
               use_raster = FALSE,
               # raster_device = 'png',
               raster_quality = 4,
               )


}



if(is.null(rowAnnotation)==FALSE){


  p<-p+rowAnnotation
  draw(p, auto_adjust = FALSE)

}else{

 return(p)

}

}





















