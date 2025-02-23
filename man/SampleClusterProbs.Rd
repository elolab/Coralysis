% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{SampleClusterProbs}
\alias{SampleClusterProbs}
\title{Sample cells based on cluster probabilities distribution}
\usage{
SampleClusterProbs(cluster, probs, q.split = 0.5)
}
\arguments{
\item{cluster}{Clustering cell labels predicted by ICP (factor).}

\item{probs}{Clustering probabilities predicted by ICP (matrix).}

\item{q.split}{Split (cell) batch principal component distribution by this 
quantile (numeric). By default {0.5}, i.e., median.}
}
\value{
A factor with cell cluster identities.
}
\description{
Samples cells based on cluster probabilities distribution
}
\keyword{internal}
