

#' This function generate the arance mat files from matrix
#' @param  log2CPM if on, use log2CPM transformation
#' @param  write.output if TRUE, output will write to files
#' @export
Mat2Aracne<-function( mat, outfile.name.suffix,
                      log2CPM=c("On", "Off"),
                      write.output=c(TRUE, FALSE),
                      sample.use="On",
                      r.cutoff=0.1,
                      sample.size=500){


  if(log2CPM=="On"){

    require(MumuscRNASeq)

mat<-Log2CPM(as.matrix(mat))

  }


  if (sample.use=="On"){

id<-sample(1:ncol(mat), size=min(sample.size, ncol(mat)))

mat<-mat[,id]
}

exprs<-mat

exprs<-round(exprs, 3)

# remove genes do not express in 95% samples

nonzeors<-apply(exprs, 1, FUN=function(x){
  length(which(x != 0))
  }
  )

r<-nonzeors/ncol(exprs)

exprs<-exprs[which(r>r.cutoff), ]


# genes<-rownames(exprs)[which(rowSums(exprs)>0)]
#
# r.id<-which(rownames(exprs)%in%genes==TRUE)
#
# exprs <- exprs[r.id, ]

exprs<-as.data.frame(as.matrix(exprs))

exprs<-round(exprs, digits=2)

colnames.exprs<-c("gene", colnames(exprs))
exprs<-cbind(rownames(exprs), exprs)
exprs<-rbind(colnames.exprs, exprs)
exprs<-as.data.frame(exprs)

rownames(exprs)<-NULL
colnames(exprs)<-NULL


out.file=paste("log2CPM_aracne_", outfile.name.suffix, ".dat", sep="")

if(write.output==TRUE){

write.table(exprs, file=out.file, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)

}

return(exprs)

}



