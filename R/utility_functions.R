#' Convert DNA sequences to binary code based on various encoding schemes
#' 
#' @param sequences List of sequences.
#' @param encoding The encoding scheme used to convert sequences to binary vector, including 3L+1, 4L+1, 5L+1. 3L+1 as default.
#' @return Binary codes converted from sequences
#' @examples
#' SeqToByte("GTAC")
#' SeqToByte("GTAAN", encoding = "3L+1“)
SeqToByte <- function(sequences, encoding = "3L+1"){
  seqs = strsplit(sequences, "")

  if(encoding == "3L+1"){
    Bases=cbind('A', 'C', "G", "T", "N", "-")
    Codes=cbind("0 1 0 ", "0 0 0 ", "1 0 0 ", "0 0 1 ", "0.25 0.25 0.25 ", "0 0 0 ")}
  else if(encoding == "5L+1"){
    Bases=cbind('A', 'C', "G", "T", "M", "W", "N", "-")
    Codes=cbind("0 1 0 0 0 ", "0 0 0 0 0 ", "1 0 0 0 0 ", "0 0 1 0 0 ", "0 0 0 1 0 ", "1 0 0 0 1 ", "0.25 0.25 0.25 0 0 ", "0 0 0 0 0 ")
  }
  else if(encoding == "4L+1"){
    Bases=cbind('A', 'C', "G", "T", "M", "W", "N", "-")
    Codes=cbind("0 1 0 0 ", "0 0 0 0 ", "1 0 0 0 ", "0 0 1 0 ", "0 0 0 1 ", "1 0 0 0 ", "0.25 0.25 0.25 0 ", "0 0 0 0 ")
  }
  purrr::map_chr(seqs,
                 function(seq){
                   Codes[match(seq, Bases)] %>%
                     stringr::str_c(collapse = "")
                 }
  )
}

#' Convert DNA sequences to binary code only including methylation-related bytes
#' 
#' @param sequences List of sequences.
#' @param encoding Either (3+1)L+1 or (3+1)L+1
#' @return Binary codes converted from sequences, only M or W positions have non-zero bytes.
#' @examples
#' SeqToMethylByte("GTAMW")
#' SeqToByte("GTAMWCGN", encoding = "(3+2)L+1“)
SeqToMethylByte <- function(sequences, encoding = "(3+2)L+1"){
  seqs = strsplit(sequences, "")
  
  if(encoding == "(3+2)L+1"){
    Bases=cbind("M", "W", "A", "C", "G", "T", "N", "-")
    Codes=cbind("1 0 ", "0 1 ", "0 0 ", "0 0 ", "0 0 ", "0 0 ", "0 0 ", "0 0 ")
  }
  else if(encoding == "(3+1)L+1"){
    Bases=cbind("M", "W", "A", "C", "G", "T", "N", "-")
    Codes=cbind("1 ", "0 ", "0 ", "0 ", "0 ", "0 ", "0 ", "0 ")
  }
  purrr::map_chr(seqs,
                 function(seq){
                   Codes[match(seq, Bases)] %>%
                     stringr::str_c(collapse = "")
                 }
  )
}


#' Convert energy matrix to energy model coefficients
#' Encoding scheme: 3L+1
#' @param PEM
#' @return energy model coefficients
#' @examples
#' toCoeffs("GTAC")
toCoeffs <- function(PEM) {
  coeff <- array(0L, dim = 3 * ncol(PEM))
  names(coeff) <- paste0(rep(1:ncol(PEM), each = 3), c("CG", "CA", "CT"))

  for (i in 1:ncol(PEM)) {
    coeff[paste0(i, "CG")] <- PEM["G", i] - PEM["C", i]
    coeff[paste0(i, "CA")] <- PEM["A", i] - PEM["C", i]
    coeff[paste0(i, "CT")] <- PEM["T", i] - PEM["C", i]
  }

  return(coeff)
}

#' Enumerate all k-mer sequences
#' 
#' @param k
#' @return all \code{k}-mer sequences
#' @examples
#' kmer(6)
kmer <- function(k){
  if(k==1)
    return(c("A", "C", "G", "T"))
  else(
    return(paste0(kmer(k-1),
                  cbind(rep("A", 4^(k-1)), rep("C", 4^(k-1)), rep("G", 4^(k-1)), rep("T", 4^(k-1))
                  )))
  )
}

#' Select sequence variants based on the max number of mismatches to some reference site
#' 
#' @param data The data frame includes Sequence variable
#' @param reference Reference site
#' @param maxMismatches The max number of mismatches of each sequence to the reference site
#' @return The data frame containing all sequences with no more than \code{maxMismatches} to the reference
#' @examples
#' selectVariants(data, "GTACGAC", maxMismatches=1)
selectVariants <- function(data, reference, maxMismatches = 1) {
  return(dplyr::filter(data, countMismatch(Sequence, reference) <= maxMismatches))
}

#' Count the number of mismatches of sequence(s) to some reference site
#' 
#' @param sequences List of sequences.
#' @param reference Reference sequence.
#' @return The number of mismatches between \code{sequences} and \code{reference}.
#' @examples
#' countMismatch("GTAC", "GAAC")
#' countMismatch(c("GTAC", "GAAC"), "GAAC")
countMismatch <- function(sequences, reference){
  seqs = strsplit(sequences, "")
  ref  = strsplit(reference, "")[[1]]
  
  purrr::map_int(seqs,
                 function(seq) {sum(seq!=ref, na.rm=TRUE)})
}


#' Highlight the mismatches of sequence(s) to some reference site with '-'
#' 
#' @param sequences List of sequences.
#' @param reference Reference sequence.
#' @return Sequences with those mismatched portions to reference changed to '-'
#' @examples
#' highlightMismatch("GTAC", "GAAC")
#' highlightMismatch(c("GTAC", "GAAC"), "GAAC")
highlightMismatch<- function(sequences, reference){
  seqs = strsplit(sequences, "")
  ref  = strsplit(reference, "")[[1]]
  
  purrr::map_chr(seqs,
                 function(seq){
                   seq[seq==ref] <- '-'
                   stringr::str_c(seq,collapse = "")
                 })
}


## validate genome input

setGeneric(".validate_genome_input",
           function(genome) standardGeneric(".validate_genome_input"))

setMethod(".validate_genome_input", signature(genome = "BSgenome"),
          function(genome) {
            return(genome)
          })

setMethod(".validate_genome_input", signature(genome = "DNAStringSet"),
          function(genome) {
            return(genome)
          })

setMethod(".validate_genome_input", signature(genome = "character"),
          function(genome) {
            if (any(is.na(genome)))
              stop("No genome provided")
            if (length(genome) > 1){
              stopifnot(all(genome == genome[[1]]))
              genome <- genome[[1]]
            }
            return(BSgenome::getBSgenome(genome))
          })

setMethod(".validate_genome_input", signature(genome = "ANY"),
          function(genome) {
            stop("genome input must be a BSgenome, DNAStringSet, or FaFile",
                 "object or a string recognized by getBSgenome")
          })
