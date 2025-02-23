% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\name{FindAllClusterMarkers}
\alias{FindAllClusterMarkers}
\alias{FindAllClusterMarkers.SingleCellExperiment}
\alias{FindAllClusterMarkers,SingleCellExperiment-method}
\title{Identification of feature markers for all clusters}
\usage{
FindAllClusterMarkers.SingleCellExperiment(
  object,
  clustering.label,
  test,
  log2fc.threshold,
  min.pct,
  min.diff.pct,
  min.cells.group,
  max.cells.per.cluster,
  return.thresh,
  only.pos
)

\S4method{FindAllClusterMarkers}{SingleCellExperiment}(
  object,
  clustering.label,
  test = "wilcox",
  log2fc.threshold = 0.25,
  min.pct = 0.1,
  min.diff.pct = NULL,
  min.cells.group = 3,
  max.cells.per.cluster = NULL,
  return.thresh = 0.01,
  only.pos = FALSE
)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{clustering.label}{A variable name (of class \code{character}) available 
in the cell metadata \code{colData(object)} with the clustering labels 
(\code{character} or \code{factor}) to use.}

\item{test}{Which test to use. Only "wilcox" (the Wilcoxon rank-sum test,
AKA Mann-Whitney U test) is supported at the moment.}

\item{log2fc.threshold}{Filters out features that have log2 fold-change of the
averaged feature expression values below this threshold. Default is \code{0.25}.}

\item{min.pct}{Filters out features that have dropout rate (fraction of cells
expressing a feature) below this threshold in both comparison groups. Default is 
\code{0.1}.}

\item{min.diff.pct}{Filters out features that do not have this minimum
difference in the dropout rates (fraction of cells expressing a feature)
between the two comparison groups. Default is \code{NULL}.}

\item{min.cells.group}{The minimum number of cells in the two comparison
groups to perform the DE analysis. If the number of cells is below the
threshold, then the DE analysis of this cluster is skipped. Default is \code{3}.}

\item{max.cells.per.cluster}{The maximum number of cells per cluster if
downsampling is performed to speed up the DE analysis. Default is \code{NULL}, 
i.e., no downsampling.}

\item{return.thresh}{If \code{only.pos=TRUE}, then return only features that have the
adjusted p-value (adjusted by the Bonferroni method) below or equal to this
threshold. Default is \code{0.01}.}

\item{only.pos}{Whether to return only features that have an adjusted p-value
(adjusted by the Bonferroni method) below or equal to the threshold. Default 
is \code{FALSE}.}
}
\value{
A data frame of the results if positive results were found, else \code{NULL}.
}
\description{
\code{FindAllClusterMarkers} enables identifying feature markers for all 
clusters at once. This is done by differential expresission analysis where 
cells from one cluster are compared against the cells from the rest of the 
clusters. Feature and cell filters can be applied to accelerate the analysis, 
but this might lead to missing weak signals.
}
\examples{
# Import package
suppressPackageStartupMessages(library("SingleCellExperiment"))

# Create toy SCE data
batches <- c("b1", "b2")
set.seed(239)
batch <- sample(x = batches, size = nrow(iris), replace = TRUE)
sce <- SingleCellExperiment(assays = list(logcounts = t(iris[,1:4])),  
                            colData = DataFrame("Species" = iris$Species, 
                                               "Batch" = batch))
colnames(sce) <- paste0("samp", 1:ncol(sce))

# Markers
dge <- FindAllClusterMarkers(sce, clustering.label = "Species")
dge

}
\keyword{DE}
\keyword{analysis}
\keyword{differential}
\keyword{expression}
\keyword{feature}
\keyword{markers}
