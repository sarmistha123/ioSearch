% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pathwayOmic2_expr.R
\name{pathwayOmic2_expr}
\alias{pathwayOmic2_expr}
\title{Collates pathway information with omic2 data values for all subjects}
\usage{
pathwayOmic2_expr(pathway.genes, omic2.expr_Group1, omic2.expr_Group2)
}
\arguments{
\item{pathway.genes}{Read the information obtained from find_pathwayGenes().}

\item{omic2.expr_Group1}{Read omic1 data matrix for all Group 1 subjects obtained from omic2().}

\item{omic2.expr_Group2}{Read omic1 data matrix for all Group 2 subjects obtained from omic2().}
}
\value{
A matrix with pathway information and omic2 expression values for all subjects.
}
\description{
Collates pathway information with omic2 data values for all subjects
}
