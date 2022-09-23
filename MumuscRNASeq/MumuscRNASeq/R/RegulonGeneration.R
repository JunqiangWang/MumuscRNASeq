

#' This function generate the regulon from the output of the ARACNe
#' @param lgo2mat_dir the directory of the expression data
#' @param network_dir the directory of the ARACNe network
#' @param cutoff the param in function pruneRegulon
#' @return pruned regulon
#' @author  Junqiang Wang
#' @export
RegulGeneration<-function(log2mat_dir, network_dir, cutoff=50){

  require(viper)


  networks<-list.files(network_dir, full.names=TRUE, include.dirs = TRUE)


  net<-NULL

  for (i in networks){

    net.i<-read.table(file=i, head=TRUE)

    if(is.null(net)==TRUE){

      net<-net.i

    }else{

      net<-rbind(net, net.i)

    }

  }

  colnames(net)<-NULL

  write.table(net, file="tmp_aracne_net.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)

  exp.dat<-list.files(log2mat_dir, full.names=TRUE, include.dirs = TRUE)

  log2mat<- read.table(file=exp.dat, sep="\t", head=TRUE)

  log2mat.rownames<-log2mat[,1]

  log2mat<-log2mat[,-1]

  rownames(log2mat)<-log2mat.rownames


  sdat<-as.matrix(log2mat)

  regul <- aracne2regulon("tmp_aracne_net.txt", sdat, gene = FALSE, format ="3col")


  regul.pruned<-pruneRegulon(regul, cutoff = cutoff, adaptive = FALSE, eliminate = TRUE, wm = NULL)

  message<-paste0("processing ", i, " done!", sep="")
  message(message)

  name.reg.i<-paste0("regul.pruned.", i)

  return(regul.pruned)

}







