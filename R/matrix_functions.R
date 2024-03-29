

addMethylMatrix <- function(PEM, MethylModel, encoding = "(3+2)L+1"){
  coeff <- MethylModel$coefficients
  MethylMatrix <- matrix(nrow = 2, ncol = ncol(PEM), dimnames = list(c("M", "W")))
  
  if(encoding == "(3+2)L+1")
    for (i in 1:ncol(PEM)){
      MethylMatrix["M", i] <- coeff[paste0("`", i, "CM`")]
      MethylMatrix["W", i] <- coeff[paste0("`", i, "GW`")] 
    }
  else if(encoding == "(3+1)L+1")
    for (i in 1:(ncol(PEM)-1)){
      MethylMatrix["M", i] <- 0.5*coeff[paste0("`", i, "MW`")]
      MethylMatrix["W", i+1] <- 0.5*coeff[paste0("`", i, "MW`")] 
    }

  return(rbind(PEM, MethylMatrix))
}

#' Add some prefixed sequence to existing energy matrix
#' 
#' @param PEM existing Energy Matrix
#' @param anchor The prefixed sequence
#' @param position The starting position of the prefixed sequence
#' @param height The height of prefixed sequence, which determines how it looks at energy logo.
#' @return Energy matrix with some prefixed sequence at designated position
#' @examples
#' addAnchorMatrix(PEM, anchor = "TAGC", position = 2, height = 0.25)
addAnchorMatrix <- function(PEM = NULL, anchor, position=1, height=0.2) {
  if(is.null(PEM))
    PEM <- matrix(nrow = 4, ncol = nchar(anchor), dimnames = list(c("A", "C", "G", "T")))
  
  
  PEM[c("A","C","G","T"), position:(position+nchar(anchor)-1)] <- 0
  
  PEM["A", position - 1 + which(strsplit(anchor, "")[[1]] == "A")] <- -height
  PEM["C", position - 1 + which(strsplit(anchor, "")[[1]] == "C")] <- -height
  PEM["G", position - 1 + which(strsplit(anchor, "")[[1]] == "G")] <- -height
  PEM["T", position - 1 + which(strsplit(anchor, "")[[1]] == "T")] <- -height

  PEM["A", position - 1 + which(strsplit(anchor, "")[[1]] == "M")] <- -height
  PEM["C", position - 1 + which(strsplit(anchor, "")[[1]] == "M")] <- -height
  
  PEM["T", position - 1 + which(strsplit(anchor, "")[[1]] == "K")] <- -height
  PEM["G", position - 1 + which(strsplit(anchor, "")[[1]] == "K")] <- -height
 
  PEM["A", position - 1 + which(strsplit(anchor, "")[[1]] == "R")] <- -height
  PEM["G", position - 1 + which(strsplit(anchor, "")[[1]] == "R")] <- -height
  
  PEM["T", position - 1 + which(strsplit(anchor, "")[[1]] == "Y")] <- -height
  PEM["C", position - 1 + which(strsplit(anchor, "")[[1]] == "Y")] <- -height
  
  PEM["A", position - 1 + which(strsplit(anchor, "")[[1]] == "N")] <- -height
  PEM["G", position - 1 + which(strsplit(anchor, "")[[1]] == "N")] <- -height
  PEM["T", position - 1 + which(strsplit(anchor, "")[[1]] == "N")] <- -height
  PEM["C", position - 1 + which(strsplit(anchor, "")[[1]] == "N")] <- -height
  return(PEM)
}

#' Reverse some existing position energy matrix
#' 
#' @param PEM position energy matrix
#' @return reversed energy matrix
#' @examples
#' reverseComplement(PEM)
reverseComplement <- function(PEM){
  PEM = .validate_PEM_input(PEM)
  reverseMatrix <- matrix(0L,
                          nrow = 4, ncol = dim(PEM)[2], byrow = TRUE,
                          dimnames = list(c("A", "C", "G", "T")))
  
  reverseMatrix["A", ] <- rev(PEM["T",])
  reverseMatrix["C", ] <- rev(PEM["G",])
  reverseMatrix["G", ] <- rev(PEM["C",])
  reverseMatrix["T", ] <- rev(PEM["A",])
  
  return(reverseMatrix)
}

#' Convert defined subject to position energy matrix
setGeneric("as.PEM",
           function(subject, ...) standardGeneric("as.PEM"))

#' @describeIn as.PEM PFM
#' @export
#' @example as.PEM(PFMatrix)
setMethod("as.PEM",
          signature(subject = "matrix"),
          function(subject) {
            num <- dim(subject)[2]
            
            ## If some position equals to zero, add pseudocount for calculation.
            subject[subject == 0] <- 1
            
            PEM <-
              matrix(0L,
                     nrow = 4,
                     ncol = num,
                     dimnames = list(c("A", "C", "G", "T")))
            
            for (i in 1:num) {
              PEM[, i] <- log(sum(subject[, i]) / subject[, i] - 1)
              PEM[, i] <-
                PEM[, i] - 0.25 * sum(PEM[, i])
            }
            return(PEM)
          })


#' @describeIn as.PEM energy model
#' @export
#' @example as.PEM(model)
setMethod("as.PEM",
          signature(subject = "lm"),
          function(subject) {
            coeff <- subject$coefficients

            if(is.null(subject$seqLength))
                num <- (length(coeff) - 1) / 3
            else
                num <- subject$seqLength
  
            PEM <- matrix(nrow = 4, ncol = num, dimnames = list(c("A", "C", "G", "T")))
    
            for (i in 1:num) {
                PEM["C", i] <- 0
                PEM["G", i] <- coeff[paste0("`", i, "CG`")]
                PEM["A", i] <- coeff[paste0("`", i, "CA`")]
                PEM["T", i] <- coeff[paste0("`", i, "CT`")]
    
                PEM[, i] <- PEM[, i] - 0.25 * sum(PEM[, i])
            }
  
            return(PEM)
          })

#' @describeIn as.PEM energy model
#' @export
#' @example as.PEM(model)
setMethod("as.PEM",
          signature(subject = "elnet"),
          function(subject) {
            PEM <- matrix(nrow = 4, ncol = subject$seqLength, dimnames = list(c("A", "C", "G", "T")))
            
            for (i in 1:subject$seqLength) {
              PEM["C", i] <- 0
              PEM["G", i] <- subject$beta[paste0(i, "CG"),1]
              PEM["A", i] <- subject$beta[paste0(i, "CA"),1]
              PEM["T", i] <- subject$beta[paste0(i, "CT"),1]
              
              PEM[, i] <- PEM[, i] - 0.25 * sum(PEM[, i])
            }
            
            return(PEM)
          })

#' Derive energy matrix from position-based energy model including parameters like 1CG, 1CT, 2CA, 3CG, etc
getEnergyMatrix <- as.PEM


#' Build composite motif based on two existing motifs with specified configurations
#' 
#' @param motifA The first motif (PEM by default)
#' @param motifB The second motif (PEM by default)
#' @param configuration "F0F" by default, i.e., Forward-0bp gap-Forward orientation
#' @return compositeMotif constructed from motifA and motifB
#' @examples
#' buildCompositeMotif(SOX2.matrix, Oct4.matrix, configuration = "F0F")
composeMotif <- function(motifA, motifB, configuration = "F0F", type = "PEM") {
  motifA <- .validate_PEM_input(motifA)
  motifB <- .validate_PEM_input(motifB)
  
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
#' @param motifA The first motif (PEM by default)
#' @param motifB The second motif (PEM by default)
#' @param configuration "F0F" by default, i.e., Forward-0bp gap-Forward orientation
#' @return compositeMotif constructed from motifA and motifB
#' @examples
#' buildCompositeMotif(SOX2.matrix, Oct4.matrix, configuration = "F0F")
decomposeMotif <- function(compositeMotif, configuration = "F0F", type = "PEM") {
  
}
