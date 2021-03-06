% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indicators.R
\name{build_Gd}
\alias{build_Gd}
\title{Build double-edge adjacency matrix}
\usage{
build_Gd(adj_mat)
}
\arguments{
\item{adj_mat}{The original adjacency matrix.}
}
\value{
A double-edge adjacency matrix \code{Gd} in sparse form.
}
\description{
Build the sparse double-edge adjacency matrix \code{Gd} from a
graph adjacency matrix.
}
\keyword{internal}
