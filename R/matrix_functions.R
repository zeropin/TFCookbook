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