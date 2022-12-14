% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/StoufferIntegrateSeuClusterRow.R
\name{StoufferIntegrateSeuClusterRows}
\alias{StoufferIntegrateSeuClusterRows}
\title{This function do Stouffer integration of the PA in each cluster by rows}
\usage{
StoufferIntegrateSeuClusterRows(
  seu,
  assay = "Activity",
  slot = "data",
  top_n = NULL,
  cluster.by = NULL,
  order.cluster = NULL,
  remove.duplicates = TRUE,
  z.cutoff = NULL,
  output.format = c("NES", "log10.p.adj"),
  log10.p.adj.cutoff = NULL
)
}
\arguments{
\item{seu}{The seurat object}
}
\value{
Integreated NES
}
\description{
This function do Stouffer integration of the PA in each cluster by rows
}
