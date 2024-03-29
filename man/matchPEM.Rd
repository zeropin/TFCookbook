% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matchPEM.R
\name{matchPEM,matrix,GenomicRanges-method}
\alias{matchPEM,matrix,GenomicRanges-method}
\alias{matchPEM,matrix,character-method}
\alias{matchPEM,matrix,DNAStringSet-method}
\title{find binding sites within some GenomicRanges or DNA sequences that match some PEM}
\usage{
\S4method{matchPEM}{matrix,GenomicRanges}(PEM, subject, genome = GenomeInfoDb::genome(subject), E.cutoff = -3)

\S4method{matchPEM}{matrix,character}(PEM, subject, E.cutoff = -3)

\S4method{matchPEM}{matrix,DNAStringSet}(PEM, subject, E.cutoff = -3)
}
\arguments{
\item{PEM}{a PEM matrix}

\item{subject}{a list of DNA sequences}

\item{genome}{the genome name, e.., hg38, mm10, etc}

\item{E.cutoff}{the threshold to select hit}
}
\description{
find binding sites within some GenomicRanges or DNA sequences that match some PEM
}
\section{Methods (by class)}{
\itemize{
\item \code{matchPEM(PEM = matrix, subject = GenomicRanges)}: matrix/GenomicRanges

\item \code{matchPEM(PEM = matrix, subject = character)}: matrix/character

\item \code{matchPEM(PEM = matrix, subject = DNAStringSet)}: matrix/DNAStringSet

}}
