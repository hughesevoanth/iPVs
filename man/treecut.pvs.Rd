% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/treecut.pvs.R
\name{treecut.pvs}
\alias{treecut.pvs}
\title{A function to identify the principal variables in your variable data set, given a defined tree cut height.}
\usage{
treecut.pvs(tree, variabledata, cutheight, cormat)
}
\arguments{
\item{tree}{a hclust object of you variables, built with the cormat and the hclust method "complete"}

\item{variabledata}{your data.frame of variables used to generate your tree}

\item{cutheight}{the height at which you wish to cut your tree}

\item{cormat}{your correlation matrix}
}
\description{
This function to identify the principal variables in your variable data set, given a defined tree cut height.
}
\examples{
treecut.pvs()
}
\keyword{cut,}
\keyword{principal}
\keyword{tree}
\keyword{variables}
