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
                      Bases=cbind('C', "G", "A", "T"),
                      Codes=cbind("0 0 0 ", "1 0 0 ", "0 1 0 ", "0 0 1 ") ){
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
  num <- (model$rank-1)/3

  energyMatrix = matrix(nrow=4, ncol=num, dimnames=list(c('C', 'G', 'A', 'T')))

  for(i in 1:num){
    energyMatrix["C", i]<-0
    energyMatrix["G", i]<-coeff[paste0("`", i, "CG`")]
    energyMatrix["A", i]<-coeff[paste0("`", i, "CA`")]
    energyMatrix["T", i]<-coeff[paste0("`", i, "CT`")]

    energyMatrix[,i]<-energyMatrix[,i]-0.25*sum(energyMatrix[,i])
  }

  return(energyMatrix)
}


plotEnergyLogo <- function(energyMatrix){
  ggplot2::ggplot() +
    ggseqlogo::geom_logo(-energyMatrix, method='custom', seq_type='dna') + ggseqlogo::theme_logo() +
    ggplot2::xlab("Position") +
    ggplot2::ylab("-Energy (kT)") %>%
    return()
}
