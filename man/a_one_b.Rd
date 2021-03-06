% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{a_one_b}
\alias{a_one_b}
\title{Compute a left-multiplication with the ones matrix}
\usage{
a_one_b(a, b)
}
\arguments{
\item{a, b}{Square matrices.}
}
\value{
The square matrix \code{a * (one_mat \%*\% b)}.
}
\description{
Compute \code{a * (one_mat \%*\% b)} where \code{a}, \code{b},
\code{ones_mat} are square matrices of the same size,
and \code{ones_mat} contains all entries equal to one.
The product \code{*} is an entry-wise (Hadamard) product,
while \code{\%*\%} represents matrix multiplication.
This method is more efficient than the naive approach
when \code{a} or \code{b} are sparse.
}
\keyword{internal}
