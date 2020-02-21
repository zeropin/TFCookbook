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

getEnergyMatrix <- function(model) {
  coeff <- model$coefficients
  num <- (length(model$coefficients) - 1) / 3

  energyMatrix <- matrix(nrow = 4, ncol = num, dimnames = list(c("A", "C", "G", "T")))

  for (i in 1:num) {
    energyMatrix["C", i] <- 0
    energyMatrix["G", i] <- coeff[paste0("`", i, "CG`")]
    energyMatrix["A", i] <- coeff[paste0("`", i, "CA`")]
    energyMatrix["T", i] <- coeff[paste0("`", i, "CT`")]

    energyMatrix[, i] <- energyMatrix[, i] - 0.25 * sum(energyMatrix[, i])
  }

  return(energyMatrix)
}

toEnergyMatrix <- function(PFM){
  num <- dim(PFM)[2]
  
  ## If some position equals to zero, add pseudocount for calculation.
  PFM[PFM==0] <- 1
  
  energyMatrix <- matrix(0L, nrow=4, ncol=num, dimnames = list(c("A", "C", "G", "T")))

  for (i in 1:num){
    energyMatrix[,i] <- log(sum(PFM[,i])/PFM[,i] -1)
    energyMatrix[,i] <- energyMatrix[,i]-0.25*sum(energyMatrix[,i])
  }
  return(energyMatrix)
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

reverseComplement <- function(energyMatrix){
  reverseMatrix <- matrix(0L,
              nrow = 4, ncol = dim(energyMatrix)[2], byrow = TRUE,
              dimnames = list(c("A", "C", "G", "T")))
              
  reverseMatrix["A", ] <- rev(energyMatrix["T",])
  reverseMatrix["C", ] <- rev(energyMatrix["G",])
  reverseMatrix["G", ] <- rev(energyMatrix["C",])
  reverseMatrix["T", ] <- rev(energyMatrix["A",])
  
  return(reverseMatrix)
}


anchorMatrix <- function(anchor, upstream = 0, downstream = 0) {
  m <- matrix(0L,
    nrow = 4, ncol = upstream + nchar(anchor) + downstream, byrow = TRUE,
    dimnames = list(c("A", "C", "G", "T"))
  )

  m["A", upstream + which(strsplit(anchor, "")[[1]] == "A")] <- 1
  m["C", upstream + which(strsplit(anchor, "")[[1]] == "C")] <- 1
  m["G", upstream + which(strsplit(anchor, "")[[1]] == "G")] <- 1
  m["T", upstream + which(strsplit(anchor, "")[[1]] == "T")] <- 1

  return(m)
}

addAnchorSite <- function(energyMatrix, anchor, position=1, height=0.2) {
  energyMatrix[c("A","C","G","T"), position:(position+nchar(anchor)-1)] <- 0
  
  energyMatrix["A", position - 1 + which(strsplit(anchor, "")[[1]] == "A")] <- -height
  energyMatrix["C", position - 1 + which(strsplit(anchor, "")[[1]] == "C")] <- -height
  energyMatrix["G", position - 1 + which(strsplit(anchor, "")[[1]] == "G")] <- -height
  energyMatrix["T", position - 1 + which(strsplit(anchor, "")[[1]] == "T")] <- -height
  
  return(energyMatrix)
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
