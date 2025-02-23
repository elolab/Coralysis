% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\name{RunPCA}
\alias{RunPCA}
\alias{RunPCA.SingleCellExperiment}
\alias{RunPCA,SingleCellExperiment-method}
\title{Principal Component Analysis}
\usage{
RunPCA.SingleCellExperiment(
  object,
  assay.name,
  p,
  scale,
  center,
  threshold,
  pca.method,
  return.model,
  select.icp.tables,
  features,
  dimred.name
)

\S4method{RunPCA}{SingleCellExperiment}(
  object,
  assay.name = "joint.probability",
  p = 50,
  scale = TRUE,
  center = TRUE,
  threshold = 0,
  pca.method = "irlba",
  return.model = FALSE,
  select.icp.tables = NULL,
  features = NULL,
  dimred.name = "PCA"
)
}
\arguments{
\item{object}{A \code{SingleCellExperiment} object.}

\item{assay.name}{Name of the assay to compute PCA. One of \code{assayNames(object)}
or \code{joint.probability}. By default \code{joint.probability} is used. Use 
\code{joint.probability} to obtain an integrated embedding after running 
\code{RunParallelDivisiveICP}. One of the assays in \code{assayNames(object)}
can be provided before performing integration to assess if data requires 
integration.}

\item{p}{A positive integer denoting the number of principal components to 
calculate and select. Default is \code{50}.}

\item{scale}{A logical specifying whether the probabilities should be
standardized to unit-variance before running PCA. Default is \code{TRUE}.}

\item{center}{A logical specifying whether the probabilities should be
centered before running PCA. Default is \code{TRUE}.}

\item{threshold}{A threshold for filtering out ICP runs before PCA with the 
lower terminal projection accuracy below the threshold. Default is \code{0}.}

\item{pca.method}{A character specifying the PCA method. One of \code{"irlba"}
(default), \code{"RSpectra"} or \code{"stats"}. Set seed before, if the method 
is \code{"irlba"} to ensure reproducibility.}

\item{return.model}{A logical specifying if the PCA model should or not be 
retrieved. By default \code{FALSE}. Only implemented for \code{pca.method = "stats"}. 
If \code{TRUE}, the \code{pca.method} is coerced to \code{"stats"}.}

\item{select.icp.tables}{Select the ICP cluster probability tables to perform 
PCA. By default \code{NULL}, i.e., all are used, except if the ICP tables were
obtained with the function \code{RunParallelDivisiveICP}, in which the ICP 
tables correspond to the last round of divisive clustering for every epoch.  
A vector of \code{integers} should be given otherwise.}

\item{features}{A character of feature names matching \code{row.names(object)} 
to select from before computing PCA. Only used if \code{assay.name} is one of
the assays in \code{assayNames(object)}, otherwise it is ignored.}

\item{dimred.name}{Dimensional reduction name given to the returned PCA. By 
default \code{"PCA"}.}
}
\value{
object of \code{SingleCellExperiment} class
}
\description{
Perform principal component analysis using assays or the joint 
probability matrix as input.
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

# Prepare SCE object for analysis
sce <- PrepareData(sce)

# Multi-level integration (just for highlighting purposes; use default parameters)
set.seed(123)
sce <- RunParallelDivisiveICP(object = sce, batch.label = "Batch", 
                              k = 2, L = 25, C = 1, train.k.nn = 10, 
                              train.k.nn.prop = NULL, use.cluster.seed = FALSE,
                              build.train.set = FALSE, ari.cutoff = 0.1, 
                             threads = 2)

# Integrated PCA
set.seed(125) # to ensure reproducibility for the default 'irlba' method
sce <- RunPCA(object = sce, assay.name = "joint.probability", p = 10)

# Plot result 
cowplot::plot_grid(PlotDimRed(object = sce, color.by = "Batch", 
                              legend.nrow = 1),
                   PlotDimRed(object = sce, color.by = "Species", 
                             legend.nrow = 1), ncol = 2)

}
\keyword{PCA}
\keyword{eigendecomposition}
