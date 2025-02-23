% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Utils.R
\name{Scale}
\alias{Scale}
\title{Scale a sparse matrix by row or column}
\usage{
Scale(x, center = TRUE, scale = TRUE, scale.by = "col")
}
\arguments{
\item{x}{A matrix of class `dgCMatrix`.}

\item{center}{A logical. By default \code{TRUE}. Subtract the values by the 
row or column mean depending on the `scale.by` parameter.}

\item{scale}{A logical. By default \code{TRUE}. Divide the values by the row 
or column standard deviation depending on the `scale.by` parameter}

\item{scale.by}{Scale by `row` or `col` (=column), i.e., use the row or 
column mean and/or standard deviations to center and /or scale the data. 
Default is \code{col}.}
}
\value{
A matrix of class `dgCMatrix`.
}
\description{
Faster implementation of \code{scale} function. It diverges from 
the \code{scale} function by not performing the root-mean-square when 
\code{scale=TRUE} and \code{center=FALSE}. In this case it divides the values
by the standard deviation of the column or row used (depending on \code{scale.by}).
}
\keyword{internal}
