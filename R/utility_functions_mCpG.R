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
  num <- (length(model$coefficients) - 1) / 5

  energyMatrix = matrix(nrow=6, ncol=num, dimnames=list(c('A', 'C', 'G', 'T', 'M', 'W')))

  for(i in 1:num){
    energyMatrix["C", i]<-0
    energyMatrix["G", i]<-coeff[paste0("`", i, "CG`")]
    energyMatrix["A", i]<-coeff[paste0("`", i, "CA`")]
    energyMatrix["T", i]<-coeff[paste0("`", i, "CT`")]
    energyMatrix[1:4,i]<-energyMatrix[1:4,i]-0.25*sum(energyMatrix[1:4,i])
    
    energyMatrix["M", i]<-coeff[paste0("`", i, "CM`")]
    energyMatrix["W", i]<-coeff[paste0("`", i, "GW`")]

    
  }

  return(energyMatrix)
}

plotEnergyLogo_MW <- function(energyMatrix){

  ggplot2::ggplot() +
    ggseqlogo::geom_logo(-energyMatrix,
                         method='custom',
                         namespace = 'ACGTMW') +
    ggseqlogo::theme_logo() +
    ggplot2::xlab("Position") +
    ggplot2::ylab("-Energy (kT)") %>%
    return()
}
