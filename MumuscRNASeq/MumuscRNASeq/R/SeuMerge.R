

#' This function merge the seurat object
#' @param  seu.list list of seurat object
#' @export
SeuMerge<-function(seu.list, add.cell.ids = NULL, project = "Seu"){

require(Seurat)



if(is.null(add.cell.ids)==TRUE){

if(is.null(names(seu.list))==TRUE){

  add.cell.ids<-paste0("Set", 1:length(seu.list))}else {

    add.cell.ids<-names(seu.list)
  }

}



obj<-vector()


for (i in 1:length(seu.list)){

  print(i)

    obj <- c(obj, seu.list[[i]])

}


# for (i in 1:length(seu.list)){
#
#   print(i)
#
#   if(length(obj)==0){obj<-seu.list[[1]]}else{
#
#     obj <- merge(obj, y=seu.list[[i]])
#
#   }
# }

seu.merged <- merge(seu.list[[1]], y = unlist(seu.list[-1]), add.cell.ids = add.cell.ids, project = project, merge.data = TRUE)




require(tidyr)


if(is.null(seu.merged@assays$RNA)==FALSE){


  message("remvoing NAs from RNA assay")


counts<-as.data.frame(seu.merged@assays$RNA@counts) %>% drop_na()

seu.merged@assays$RNA@counts<-as.matrix(counts)
}




if(is.null(seu.merged@assays$Activity)==FALSE){

  message("remvoing NAs from Activity assay")

data<-as.data.frame(seu.merged@assays$Activity@data) %>% drop_na()

seu.merged@assays$Activity@data<-as.matrix(data)

}

return(seu.merged)

}









