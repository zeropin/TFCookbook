% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fit_models.R
\name{predictEnergyMW}
\alias{predictEnergyMW}
\title{Predict the binding energy of input sequences based on energy model including methylation parameters}
\usage{
predictEnergyMW(sequences, model, encoding = "4L+1")
}
\arguments{
\item{sequences}{Input sequences, including M and W as methylated cytosines at upper and lower strands respectively}

\item{model}{Linear model containing position-based energy coefficients, e.g., 1CA, 1CG, 1CT, 1MW, 2CA, etc.}

\item{encoding}{Encoding scheme. By default 4L+1}
}
\value{
Predicted binding energy values of \code{sequences} based on \code{model}
}
\description{
Predict the binding energy of input sequences based on energy model including methylation parameters
}
\examples{
predictEnergyMW(sequences, model, encoding = "4L+1")
}
