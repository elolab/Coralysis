% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoreMethods.R
\name{RunUMAP}
\alias{RunUMAP}
\alias{RunUMAP.SingleCellExperiment}
\alias{RunUMAP,SingleCellExperiment-method}
\title{Uniform Manifold Approximation and Projection (UMAP)}
\usage{
RunUMAP.SingleCellExperiment(
  object,
  dims,
  dimred.type,
  return.model,
  umap.method,
  dimred.name,
  ...
)

\S4method{RunUMAP}{SingleCellExperiment}(
  object,
  dims = NULL,
  dimred.type = "PCA",
  return.model = FALSE,
  umap.method = "umap",
  dimred.name = "UMAP",
  ...
)
}
\arguments{
\item{object}{An object of \code{SingleCellExperiment} class.}

\item{dims}{Dimensions to select from \code{dimred.type}. By default \code{NULL}, 
i.e., all the dimensions are selected. Provide a numeric vector to select a 
specific range, e.g., \code{dims = 1:10} to select the first 10 dimensions.}

\item{dimred.type}{Dimensional reduction type to use. By default \code{"PCA"}.}

\item{return.model}{Return UMAP model. By default \code{FALSE}.}

\item{umap.method}{UMAP method to use: \code{"umap"} or \code{"uwot"}. 
By default \code{"umap"}.}

\item{dimred.name}{Dimensional reduction name given to the returned UMAP. 
By default \code{"UMAP"}.}

\item{...}{Parameters to be passed to the \code{umap} function. The parameters 
given should match the parameters accepted by the \code{umap} function depending
on the \code{umap.method} given. Check possible parameters with \code{?umap::umap} 
or \code{?uwot::umap} depending if \code{umap.method} is \code{"umap"} or 
\code{"uwot"}.}
}
\value{
A \code{SingleCellExperiment} object.
}
\description{
Run nonlinear dimensionality reduction using UMAP with a dimensional 
reduction as input.
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

# Run UMAP
set.seed(123)
sce <- RunUMAP(sce, dimred.type = "PCA")

# Plot results
# Plot result 
cowplot::plot_grid(PlotDimRed(object = sce, color.by = "Batch", 
                              legend.nrow = 1),
                   PlotDimRed(object = sce, color.by = "Species", 
                             legend.nrow = 1), ncol = 2)
                             
}
\keyword{Approximation}
\keyword{Manifold}
\keyword{Projection}
\keyword{UMAP}
\keyword{Uniform}
\keyword{and}
