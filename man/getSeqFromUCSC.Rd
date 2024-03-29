% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utility_functions.R
\name{getSeqFromUCSC}
\alias{getSeqFromUCSC}
\title{Get DNA sequences from UCSC genome database directly through REST API}
\usage{
getSeqFromUCSC(genome, names, start = NA, end = NA, strand = "+")
}
\arguments{
\item{genome}{genome assembly name, hg38, mm10, panPan3, etc}

\item{names}{either a GRanges object or a character of the chromosome name}

\item{start, end}{Vector of integers (eventually with NAs) specifying the locations of the subsequences to extract}

\item{strand}{Vector of integers (eventually with NAs) specifying the locations of the subsequences to extract}
}
\description{
Get DNA sequences from UCSC genome database directly through REST API
}
\examples{
getSeqFromUCSC("hg38", "chr1", start=55, end=65, strand="-")

ranges = GenomicRanges::GRanges(seqnames = "chrM", ranges = IRanges::IRanges(15466, 15555), strand="+")
getSeqFromUCSC("panPan3", names=ranges)
}
