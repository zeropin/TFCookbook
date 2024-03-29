% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matrix_functions.R
\name{as.PEM}
\alias{as.PEM}
\alias{as.PEM,matrix-method}
\alias{as.PEM,lm-method}
\alias{as.PEM,elnet-method}
\alias{getEnergyMatrix}
\title{Convert defined subject to position energy matrix}
\usage{
as.PEM(subject, ...)

\S4method{as.PEM}{matrix}(subject)

\S4method{as.PEM}{lm}(subject)

\S4method{as.PEM}{elnet}(subject)

as.PEM(subject, ...)
}
\description{
Convert defined subject to position energy matrix

Derive energy matrix from position-based energy model including parameters like 1CG, 1CT, 2CA, 3CG, etc
}
\section{Methods (by class)}{
\itemize{
\item \code{as.PEM(matrix)}: PFM

\item \code{as.PEM(lm)}: energy model

\item \code{as.PEM(elnet)}: energy model

}}
