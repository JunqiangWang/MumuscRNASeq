

#' This function adds a named assay 'activity' to a Seurat object.
#' @export
AddActivityAssay <- function(seu, vp_mat) {
  seu.pa <- CreateAssayObject(data = as.matrix(vp_mat))
  seu[['Activity']] <- seu.pa
  seu@active.assay <- 'Activity'
  return(seu)
}


