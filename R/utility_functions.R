countMismatch <- function(sequences, reference){
  seqs = strsplit(sequences, "")
  ref  = strsplit(reference, "")[[1]]

  purrr::map_int(seqs,
                 function(seq) {sum(seq!=ref, na.rm=TRUE)})
}

highlightMismatch<- function(sequences, reference){
  seqs = strsplit(sequences, "")
  ref  = strsplit(reference, "")[[1]]

  purrr::map_chr(seqs,
                 function(seq){
                   seq[seq==ref] <- '-'
                   stringr::str_c(seq,collapse = "")
                 })
}

SeqToByte <- function(sequences,
                      Bases=cbind('A', 'C', "G", "T", "N"),
                      Codes=cbind("0 1 0 ", "0 0 0 ", "1 0 0 ", "0 0 1 ", "0.25 0.25 0.25 ") ){
  seqs = strsplit(sequences, "")

  purrr::map_chr(seqs,
                 function(seq){
                   Codes[match(seq, Bases)] %>%
                     stringr::str_c(collapse = "")
                 }
  )
}


toCoeffs <- function(energyMatrix) {
  coeff <- array(0L, dim = 3 * ncol(energyMatrix))
  names(coeff) <- paste0(rep(1:ncol(energyMatrix), each = 3), c("CG", "CA", "CT"))

  for (i in 1:ncol(energyMatrix)) {
    coeff[paste0(i, "CG")] <- energyMatrix["G", i] - energyMatrix["C", i]
    coeff[paste0(i, "CA")] <- energyMatrix["A", i] - energyMatrix["C", i]
    coeff[paste0(i, "CT")] <- energyMatrix["T", i] - energyMatrix["C", i]
  }

  return(coeff)
}


kmer <- function(k){
  if(k==1)
    return(c("A", "C", "G", "T"))
  else(
    return(paste0(kmer(k-1),
                  cbind(rep("A", 4^(k-1)), rep("C", 4^(k-1)), rep("G", 4^(k-1)), rep("T", 4^(k-1))
                  )))
  )
}

selectVariants <- function(data, reference, maxMismatches = 1) {
  return(dplyr::filter(data, countMismatch(Sequence, reference) <= maxMismatches))
}
