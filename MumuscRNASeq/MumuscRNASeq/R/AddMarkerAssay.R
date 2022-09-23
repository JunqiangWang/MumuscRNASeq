

#' This function adds a named assay 'activity' to a Seurat object.
#' @export
AddMarkerAssay <- function(seu, marker.mat) {
  seu.marker <- CreateAssayObject(data = as.matrix(marker.mat))
  seu[['Marker']] <- seu.marker
  seu@active.assay <- 'Marker'
  return(seu)
}


