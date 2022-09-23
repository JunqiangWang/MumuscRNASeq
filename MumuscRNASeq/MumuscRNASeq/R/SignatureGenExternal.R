
#'This function generates the signatue of an expression mat using the external expression mat

#' @param exp The matrix
#' @param exp.ref The external matrix as centroid reference
#' @param useLogTransform Whether process log transformation on the input datasets
#' @author Junqiang Wang
#' @export
SignaureGenExternal<-function(exp, exp.ref, useLogTransform=TRUE){


 if(useLogTransform==TRUE){
   message("Process log transformation!")

   exp<-Log2CPM(exp)
   exp.ref<-Log2CPM(exp.ref)

 } else{}


  # Get the overlapped genes
  tmp <- intersect(rownames(exp), rownames(exp.ref))
  tmp1<-match(tmp, rownames(exp)) ; tmp2<-match(tmp, rownames(exp.ref))
  exp <- exp[tmp1,];  exp.ref<- exp.ref[tmp2,]



  # GTex caudate centroid reference
  tmp.mean <- apply(exp.ref, 1, mean)
  tmp.sd <- apply(exp.ref, 1, sd)

  exp.ges <- apply(exp, 2, function(x){(x - tmp.mean) / tmp.sd})

  exp.ges <- exp.ges[is.finite(rowSums(exp.ges)),]
  exp.ges<-round(exp.ges, digits=2)

return(exp.ges)

}


