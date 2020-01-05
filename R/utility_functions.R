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
                      Bases=cbind('A', 'C', "G", "T"),
                      Codes=cbind("0 1 0 ", "0 0 0 ", "1 0 0 ", "0 0 1 ") ){
  seqs = strsplit(sequences, "")

  purrr::map_chr(seqs,
                 function(seq){
                   Codes[match(seq, Bases)] %>%
                     stringr::str_c(collapse = "")
                 }
  )
}

getEnergyMatrix <- function(model){
  coeff <- model$coefficients
  num <- (length(model$coefficients)-1)/3

  energyMatrix = matrix(nrow=4, ncol=num, dimnames=list(c('A', 'C', 'G', 'T')))

  for(i in 1:num){
    energyMatrix["C", i]<-0
    energyMatrix["G", i]<-coeff[paste0("`", i, "CG`")]
    energyMatrix["A", i]<-coeff[paste0("`", i, "CA`")]
    energyMatrix["T", i]<-coeff[paste0("`", i, "CT`")]

    energyMatrix[,i]<-energyMatrix[,i]-0.25*sum(energyMatrix[,i])
  }

  return(energyMatrix)
}

getCoeffs <- function(matrix) {
  coeff <- array(0L, dim = 3 * ncol(matrix))
  names(coeff) <- paste0(rep(1:ncol(matrix), each = 3), c("CG", "CA", "CT"))

  for (i in 1:ncol(matrix)) {
    coeff[paste0(i, "CG")] <- matrix["G", i] - matrix["C", i]
    coeff[paste0(i, "CA")] <- matrix["A", i] - matrix["C", i]
    coeff[paste0(i, "CT")] <- matrix["T", i] - matrix["C", i]
  }

  return(coeff)
}

plotEnergyLogo <- function(energyMatrix){
  ggplot2::ggplot() +
    ggseqlogo::geom_logo(-energyMatrix, method='custom', seq_type='dna') +
    ggseqlogo::theme_logo() +
    ggplot2::xlab("Position") +
    ggplot2::ylab("-Energy (kT)") %>%
    return()
}

anchorMatrix <- function(anchor, upstream=0,downstream=0){
  m <- matrix(0L, nrow=4, ncol=upstream+nchar(anchor)+downstream, byrow=TRUE,
              dimnames=list(c("A" ,"C" ,"G" , "T")))

  m["A",upstream + which(strsplit(anchor,"")[[1]]=="A")] = 1
  m["C",upstream + which(strsplit(anchor,"")[[1]]=="C")] = 1
  m["G",upstream + which(strsplit(anchor,"")[[1]]=="G")] = 1
  m["T",upstream + which(strsplit(anchor,"")[[1]]=="T")] = 1

  return(m)
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
