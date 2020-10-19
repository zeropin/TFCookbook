#' Derive energy matrix from position-based energy model including parameters like 1CG, 1CT, 2CA, 3CG, etc
#' 
#' @param model position-based energy model, returned by buildEnergymodel function.
#' @return Energy matrix
#' @examples
#' getEnergyMatrix(model)
getEnergyMatrix <- function(model) {
  coeff <- model$coefficients
  num <- (length(coeff) - 1) / 3
  num <- model$seqLength
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


addMethylMatrix <- function(energyMatrix, MethylModel, encoding = "(3+2)L+1"){
  coeff <- MethylModel$coefficients
  MethylMatrix <- matrix(nrow = 2, ncol = ncol(energyMatrix), dimnames = list(c("M", "W")))
  
  if(encoding == "(3+2)L+1")
    for (i in 1:ncol(energyMatrix)){
      MethylMatrix["M", i] <- coeff[paste0("`", i, "CM`")]
      MethylMatrix["W", i] <- coeff[paste0("`", i, "GW`")] 
    }
  else if(encoding == "(3+1)L+1")
    for (i in 1:(ncol(energyMatrix)-1)){
      MethylMatrix["M", i] <- 0.5*coeff[paste0("`", i, "MW`")]
      MethylMatrix["W", i+1] <- 0.5*coeff[paste0("`", i, "MW`")] 
    }

  return(rbind(energyMatrix, MethylMatrix))
}

#' Add some prefixed sequence to existing energy matrix
#' 
#' @param energyMatrix existing Energy Matrix
#' @param anchor The prefixed sequence
#' @param position The starting position of the prefixed sequence
#' @param height The height of prefixed sequence, which determines how it looks at energy logo.
#' @return Energy matrix with some prefixed sequence at designated position
#' @examples
#' addAnchorMatrix(energyMatrix, anchor = "TAGC", position = 2, height = 0.25)
addAnchorMatrix <- function(energyMatrix, anchor, position=1, height=0.2) {
  energyMatrix[c("A","C","G","T"), position:(position+nchar(anchor)-1)] <- 0
  
  energyMatrix["A", position - 1 + which(strsplit(anchor, "")[[1]] == "A")] <- -height
  energyMatrix["C", position - 1 + which(strsplit(anchor, "")[[1]] == "C")] <- -height
  energyMatrix["G", position - 1 + which(strsplit(anchor, "")[[1]] == "G")] <- -height
  energyMatrix["T", position - 1 + which(strsplit(anchor, "")[[1]] == "T")] <- -height

  energyMatrix["A", position - 1 + which(strsplit(anchor, "")[[1]] == "M")] <- -height*0.5
  energyMatrix["C", position - 1 + which(strsplit(anchor, "")[[1]] == "M")] <- -height*0.5
  
  energyMatrix["T", position - 1 + which(strsplit(anchor, "")[[1]] == "K")] <- -height*0.5
  energyMatrix["G", position - 1 + which(strsplit(anchor, "")[[1]] == "K")] <- -height*0.5
 
  energyMatrix["A", position - 1 + which(strsplit(anchor, "")[[1]] == "R")] <- -height*0.5
  energyMatrix["G", position - 1 + which(strsplit(anchor, "")[[1]] == "R")] <- -height*0.5
  
  energyMatrix["T", position - 1 + which(strsplit(anchor, "")[[1]] == "Y")] <- -height*0.5
  energyMatrix["C", position - 1 + which(strsplit(anchor, "")[[1]] == "Y")] <- -height*0.5
  
  energyMatrix["A", position - 1 + which(strsplit(anchor, "")[[1]] == "N")] <- -height*0.25
  energyMatrix["G", position - 1 + which(strsplit(anchor, "")[[1]] == "N")] <- -height*0.25
  energyMatrix["T", position - 1 + which(strsplit(anchor, "")[[1]] == "N")] <- -height*0.25
  energyMatrix["C", position - 1 + which(strsplit(anchor, "")[[1]] == "N")] <- -height*0.25
  return(energyMatrix)
}

#' Reverse some existing energy matrix
#' 
#' @param energyMatrix
#' @return reversed energy matrix
#' @examples
#' reverseComplement(energyMatrix)
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

toPFM <- function(energyMatrix, miu=0){
  num <- dim(energyMatrix)[2]
  
  #PFM <- matrix(0L, nrow=4, ncol=num, dimnames = list(c("A", "C", "G", "T")))
  
  PFM <- 1000/(exp(energyMatrix-miu)+1)
  
  return(PFM)
}


#' Build composite motif based on two existing motifs with specified configurations
#' 
#' @param motifA The first motif (energyMatrix by default)
#' @param motifB The second motif (energyMatrix by default)
#' @param configuration "F0F" by default, i.e., Forward-0bp gap-Forward orientation
#' @return compositeMotif constructed from motifA and motifB
#' @examples
#' buildCompositeMotif(SOX2.matrix, Oct4.matrix, configuration = "F0F")
composeMotif <- function(motifA, motifB, configuration = "F0F", type = "energy") {
  conf <- strsplit(configuration, split="")[[1]]
  names(conf) <- c("A", "gap", "B")
  
  gap <- as.integer(conf["gap"])
  
  compositeMotif <- matrix(0L,
                           nrow = 4,
                           ncol = ncol(motifA) + ncol(motifB) + gap, 
                           dimnames = list(c("A", "C", "G", "T")))
  
  if(conf["A"] == "F")
    compositeMotif[,1:ncol(motifA)] <- motifA
  else if(conf["A"] == "R")
    compositeMotif[,1:ncol(motifA)] <- TFCookbook::reverseComplement(motifA)
  
  if(conf["B"] == "F")
    compositeMotif[,(ncol(motifA)+gap+1):(ncol(motifA)+gap+ncol(motifB))] <- motifB
  else if(conf["B"] == "R")
    compositeMotif[,(ncol(motifA)+gap+1):(ncol(motifA)+gap+ncol(motifB))] <- TFCookbook::reverseComplement(motifB)
  
  
  return(compositeMotif)
}


#' Decompose some composite motif into two separate submotifs
#' 
#' @param motifA The first motif (energyMatrix by default)
#' @param motifB The second motif (energyMatrix by default)
#' @param configuration "F0F" by default, i.e., Forward-0bp gap-Forward orientation
#' @return compositeMotif constructed from motifA and motifB
#' @examples
#' buildCompositeMotif(SOX2.matrix, Oct4.matrix, configuration = "F0F")
decomposeMotif <- function(compositeMotif, configuration = "F0F", type = "energy") {
  
}
