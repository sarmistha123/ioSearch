% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f.R
\name{f}
\alias{f}
\title{Runs iosearch algorithm to find top omics variables for significant set(s)/pathway(s)}
\usage{
f(x, mat1, mat2, phen, sets, K, ntopOmics2, ntopOmics1)
}
\arguments{
\item{x}{A character vector with one element.}

\item{mat1}{Omic1 matrix with pathway information.}

\item{mat2}{Omic2 matrix with pathway information.}

\item{phen}{Phenotype matrix.}

\item{sets}{Array of pathways.}

\item{K}{Default is 1. Tuning parameter specific to elasticnet package.}

\item{ntopOmics2}{Number of top omics1 variables to be selected from each pathway (user-defined).}

\item{ntopOmics1}{Number of top omics1 variables to be selected from each pathway (user-defined).}
}
\value{
A list.
}
\description{
Runs iosearch algorithm to find top omics variables for significant set(s)/pathway(s)
}
