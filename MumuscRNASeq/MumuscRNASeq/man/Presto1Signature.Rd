% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrestoSignature.R
\name{Presto1Signature}
\alias{Presto1Signature}
\title{This function generate the differential signatures}
\usage{
Presto1Signature(
  seu,
  cluster.by = "seurat_clusters",
  signature.cluster = NULL,
  groups.use = NULL,
  reverse.clusters = c(TRUE, FALSE),
  assay = "RNA",
  slot = c("counts", "scale.data", "data")
)
}
\arguments{
\item{seu}{seurat object}

\item{cluster.by}{the ID in seu@meata.data}

\item{groups.use}{the groups that are used run the presto}

\item{reverse.clusters}{whether the order of the cluster id shoud be reversed}

\item{assay}{The assay that is used}

\item{slot}{the slot of the assay that is used}

\item{cluster.order}{the order of the cluster}

\item{top.n}{return the top.n dif features}

\item{logFC.cutoff}{the cutoff of the presto}
}
\value{
return the markers and store the markers in a new assay
}
\description{
This function generate the differential signatures
}
\author{
Junqiang Wang
}
