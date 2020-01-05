SeqToByte_MW <- function(sequences){
  seqs = strsplit(sequences, "")
  Bases=cbind('C', "G", "A", "T", "M", "W")
  Codes=cbind("0 0 0 0 0 ", "1 0 0 0 0 ", "0 1 0 0 0 ", "0 0 1 0 0 ",
              "0 0 0 1 0 ", "1 0 0 0 1 ")

  purrr::map_chr(seqs,
                 function(seq){
                   Codes[match(seq, Bases)] %>%
                     stringr::str_c(collapse = "")
                 }
  )
}


getEnergyMatrix_MW <- function(model){
  coeff <- model$coefficients
  num <- (model$rank-1)/3

  energyMatrix = matrix(nrow=6, ncol=num, dimnames=list(c('C', 'G', 'A', 'T', 'M', 'W')))

  for(i in 1:num){
    energyMatrix["C", i]<-0
    energyMatrix["G", i]<-coeff[paste0("`", i, "CG`")]
    energyMatrix["A", i]<-coeff[paste0("`", i, "CA`")]
    energyMatrix["T", i]<-coeff[paste0("`", i, "CT`")]
    energyMatrix["M", i]<-coeff[paste0("`", i, "CM`")]
    energyMatrix["W", i]<-coeff[paste0("`", i, "GW`")]

    energyMatrix[,i]<-energyMatrix[1:4,i]-0.25*sum(energyMatrix[2:4,i])
  }

  return(energyMatrix)
}

plotEnergyLogo_MW <- function(energyMatrix){

  ggplot2::ggplot() +
    ggseqlogo::geom_logo(-energyMatrix,
                         method='custom',
                         namespace <- c('C', 'G', 'A', 'T', 'M', 'W')) +
    ggseqlogo::theme_logo() +
    ggplot2::xlab("Position") +
    ggplot2::ylab("-Energy (kT)") %>%
    return()
}
