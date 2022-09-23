

#' This function generate the regulon
#' @author Junqiang Wang
#' @export
Aracne2Regulon<-function(log2mat_path, network_dir_path, cutoff = 50, adaptive = FALSE, eliminate = TRUE, wm = NULL){

library(viper)


  networks<-list.files(network_dir_path, full.names=TRUE, include.dirs = TRUE)


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

  log2mat.file<-list.files(log2mat_path, full.names=TRUE, include.dirs = TRUE)

  log2mat<- read.table(file=log2mat.file, sep="\t", head=TRUE)

  log2mat.rownames<-log2mat[,1]

  log2mat<-log2mat[,-1]

  rownames(log2mat)<-log2mat.rownames

  sdat<-as.matrix(log2mat)

  regul <- aracne2regulon(afile="tmp_aracne_net.txt", eset=sdat, gene = FALSE, format ="3col")

  regul.pruned<-pruneRegulon(regul, cutoff = cutoff, adaptive = adaptive, eliminate = eliminate, wm = wm)

  unlink("tmp_aracne_net.txt")

  return(regul.pruned)

}





#' This function generate the regulon
#' @author Junqiang Wang
#' @export
Aracne2RegulonMono<- function(log2mat.file, network.file, cutoff = 50, adaptive = FALSE, eliminate = TRUE, wm = NULL){


  library(viper)

  log2mat<- read.table(file=log2mat.file, sep="\t", head=TRUE)

  log2mat.rownames<-log2mat[,1]

  log2mat<-log2mat[,-1]

  rownames(log2mat)<-log2mat.rownames

  sdat<-as.matrix(log2mat)

  regul <- aracne2regulon(afile=network.file, eset=sdat, gene = FALSE, format ="3col")

  regul.pruned<-pruneRegulon(regul, cutoff = cutoff, adaptive = adaptive, eliminate = eliminate, wm = wm)

  return(regul.pruned)

}






#' This function generate the regulon
#' @author Junqiang Wang
#' @export
Aracne2Regulon2<-function(network_dir_path, cutoff = 50, adaptive = FALSE, eliminate = TRUE, wm = NULL){

  library(viper)


  files<-list.files(network_dir_path, full.names=TRUE, include.dirs = TRUE)

  files.dat<-files[grep("dat", files)]

  files.networks<-setdiff(files, files.dat)


  net<-NULL

  for (i in files.networks){

    net.i<-read.table(file=i, head=TRUE)

    if(is.null(net)==TRUE){

      net<-net.i

    }else{

      net<-rbind(net, net.i)

    }

  }


  colnames(net)<-NULL

  write.table(net, file="tmp_aracne_net.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE, append=FALSE)


  log2mat<- read.table(file=files.dat, sep="\t", head=TRUE)

  log2mat.rownames<-log2mat[,1]

  log2mat<-log2mat[,-1]

  rownames(log2mat)<-log2mat.rownames

  sdat<-as.matrix(log2mat)

  regul <- aracne2regulon(afile="tmp_aracne_net.txt", eset=sdat, gene = FALSE, format ="3col")

  regul.pruned<-pruneRegulon(regul, cutoff = cutoff, adaptive = adaptive, eliminate = eliminate, wm = wm)

  unlink("tmp_aracne_net.txt")

  return(regul.pruned)

}







#' This function generate the regulon
#' @author Junqiang Wang
#' @export
Aracne2RegulonList<-function(network_dir_path, cutoff = 50, adaptive = FALSE, eliminate = TRUE, wm = NULL){

  require(viper)

  network_dir_paths<-list.files(network_dir_path, full.names=TRUE, include.dirs = TRUE)

  network_dir_names<-list.files(network_dir_path, full.names=FALSE, include.dirs = TRUE)

  regul.list<-list()


  for (i in 1:length(network_dir_paths)){

    print(i)

    regul.i<-Aracne2Regulon2(network_dir_path=network_dir_paths[i], cutoff = cutoff, eliminate = eliminate)

    message<-paste0("processing ", network_dir_names[i], " done!", sep="")

    name.regul.pruned.i<-network_dir_names[i]

    assign(name.regul.pruned.i, regul.i)

    # regul.list<-append(regul.list, regul.i)

    regul.list<-append(regul.list,  list( assign(name.regul.pruned.i, regul.i)))

  }

  names(regul.list)<-network_dir_names

  return(regul.list)

}


















