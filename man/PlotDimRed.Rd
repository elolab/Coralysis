% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VizMethods.R
\name{PlotDimRed}
\alias{PlotDimRed}
\alias{PlotDimRed.SingleCellExperiment}
\alias{PlotDimRed,SingleCellExperiment-method}
\title{Plot dimensional reduction categorical variables}
\usage{
PlotDimRed.SingleCellExperiment(
  object,
  color.by,
  dimred,
  dims,
  use.color,
  point.size,
  point.stroke,
  legend.nrow,
  seed.color,
  label,
  plot.theme,
  rasterise,
  rasterise.dpi,
  legend.justification,
  legend.size,
  legend.title
)

\S4method{PlotDimRed}{SingleCellExperiment}(
  object,
  color.by,
  dimred = tail(reducedDimNames(object), n = 1),
  dims = 1:2,
  use.color = NULL,
  point.size = 1,
  point.stroke = 1,
  legend.nrow = 2,
  seed.color = 123,
  label = FALSE,
  plot.theme = theme_classic(),
  rasterise = (ncol(object) <= 30000),
  rasterise.dpi = 300,
  legend.justification = "center",
  legend.size = 10,
  legend.title = color.by
)
}
\arguments{
\item{object}{An object of \code{SingleCellExperiment} class.}

\item{color.by}{Categorical variable available in \code{colData(object)} to 
plot.}

\item{dimred}{Dimensional reduction available in \code{ReducedDimNames(object)}
to plot. By default the last dimensional reduction in the object is used.}

\item{dims}{Dimensions from the dimensional reduction embedding to plot.}

\item{use.color}{Character specifying the colors. By default \code{NULL}, i.e., 
colors are randomly chosen based on the seed given at \code{seed.color}.}

\item{point.size}{Size of points. By default \code{1}.}

\item{point.stroke}{Size of stroke. By default \code{1}.}

\item{legend.nrow}{Display legend items by this number of rows. By default \code{2}.}

\item{seed.color}{Seed to randomly select colors. By default \code{123}.}

\item{label}{Logical to add or not categorical labels to the centroid categories. 
By default \code{FALSE}, i.e., labels are not added.}

\item{plot.theme}{Plot theme available in \code{ggplot2}. By default \code{theme_classic()}.}

\item{rasterise}{Logical specifying if points should be rasterised or not. By 
default \code{TRUE}, if more than 3e4 cells, otherwise \code{FALSE}.}

\item{rasterise.dpi}{In case \code{rasterise = TRUE}, DPI to use. By default 
\code{300}.}

\item{legend.justification}{Legend justification. By default \code{"center"}.}

\item{legend.size}{Legend size. By default \code{10}}

\item{legend.title}{Legend title. By default the same as given at \code{color.by}.}
}
\value{
A plot of class \code{ggplot}.
}
\description{
Plot categorical variables in dimensional reduction.
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

# Compute dimensional reduction 
sce <- RunPCA(object = sce, assay.name = "logcounts", p = 4, 
            pca.method = "stats")

# Plot batch 
PlotDimRed(object = sce, color.by = "Batch", dimred = "PCA", legend.nrow = 1)

# Plot cell type annotations
PlotDimRed(object = sce, color.by = "Species", legend.nrow = 1, 
          dimred = "PCA", label = TRUE)

}
\keyword{Dimensional}
\keyword{reduction}
\keyword{visualization}
