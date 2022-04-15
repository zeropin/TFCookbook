#' Build energy model based on Sequence-Energy relationship
#' 
#' @param data Data frame containing at least two variables -- Sequence and Energy
#' @param encoding The encoding scheme used to build model, including 3L+1, 4L+1, 5L+1.
#' @return Linear model containing position-based energy coefficients, e.g., 1CA, 1CG, 1CT, 2CA, etc.
#' @examples
#' buildEnergyModel(data)
buildEnergyModel <- function(data, encoding = "3L+1") {
  num <- nchar(data$Sequence[1])

  if(encoding == "3L+1") scheme = paste0(rep(1:num, each = 3), c("CG", "CA", "CT"))
  else if(encoding == "5L+1") scheme = paste0(rep(1:num, each = 5), c("CG", "CA", "CT", "CM", "GW"))
  else if(encoding == "4L+1") scheme = paste0(rep(1:num, each = 4), c("CG", "CA", "CT", "MW"))
  
  model <- data %>%
    dplyr::mutate(ByteCode = SeqToByte(Sequence, encoding)) %>%
    dplyr::select(Energy, ByteCode) %>%
    tidyr::separate(ByteCode,
      into = scheme,
      convert = TRUE
    ) %>%
    lm(Energy ~ . - Energy, data = .)
  
  model$seqLength <- num
  return(model)
}


#' Build model of methylation effect based on pair-wise compared data between methylated and unmethylated sites
#' 
#' @param data Data frame containing at least two variables -- Sequence and Energy
#' @param encoding The encoding scheme used for modeling, either (3+2)L+1 or (3+1)L+1.
#' @param withIntercept Working with pairwise compaired data doesn't require intercept parameter, so it is FALSE by default.
#' @return Linear model containing position-based energy coefficients, e.g., 1CM, 1GW, 2CM, 2GW, etc.
#' @examples
#' buildMehylationModel(data, encoding = "(3+1)L+1")
buildMethylationModel <- function(data, encoding = "(3+2)L+1", withIntercept = FALSE) {
  num <- nchar(data$Sequence[1])
  
  if(encoding == "(3+2)L+1") scheme = paste0(rep(1:num, each = 2), c("CM", "GW"))
  else if(encoding == "(3+1)L+1") scheme = paste0(rep(1:num, each = 1), c("MW"))
  
  data %>%
    dplyr::mutate(MethylCode = SeqToMethylByte(Sequence, encoding)) %>%
    dplyr::select(Energy, MethylCode) %>%
    tidyr::separate(MethylCode,
                    into = scheme,
                    convert = TRUE
    ) -> input
  
  if(withIntercept == TRUE)
    return(lm(Energy ~ . - Energy, data = input))
  else if(withIntercept == FALSE)
    return(lm(Energy ~ . - Energy + 0, data = input))
}


#' Predict the binding energy of input sequences based on energy model
#' @export
#' @param sequences Input sequences
#' @param model Linear model containing position-based energy coefficients, e.g., 1CA, 1CG, 1CT, 2CA, etc.
#' @return Predicted binding energy values of \code{sequences} based on \code{model} or Position energy matrix (PEM)
#' @examples
#' predictEnergy(sequences, model)
setGeneric("predictEnergy",
           function(sequences, model, ...) standardGeneric("predictEnergy"))

#' @describeIn predictEnergy Sequences/EnergyModel
#' @export
setMethod("predictEnergy", signature(sequences = "character",
                                     model     = "lm"),
          function(sequences, model) {
            #num <- (model$rank - 1) / 3  
            num <- nchar(sequences[[1]]) 
            
            #num shoulbe be the number of nucleotide positions used for prediction，model$rank is the number of coefficients in the models.
            # In some instances, there are NA values for some parameters, which are not included in total rank, so nchar(sequence) is used instead.
            
            tibble(ByteCode = SeqToByte(sequences)) %>%
              tidyr::separate(ByteCode,
                              into = paste0(rep(1:num, each = 3), c("CG", "CA", "CT")),
                              convert = TRUE
              ) %>%
              predict.lm(model, newdata = ., type = "response")
          })


#' @describeIn predictEnergy Sequences/PEM
#' @export
setMethod("predictEnergy", signature(sequences = "character",
                                     model     = "matrix"),
          function(sequences, model) {
            model <- .validate_PEM_input(model)
            num <- ncol(model) 
            
            purrr::map_dbl(sequences, function(seq){
              bases = strsplit(seq, "")[[1]][1:num]
              
              return(sum(model['A', which(bases=="A")])+
                     sum(model['C', which(bases=="C")])+
                     sum(model['T', which(bases=="T")])+
                     sum(model['G', which(bases=="G")]))
              })
      
          })

 

#' Predict the binding energy of input sequences based on energy model including methylation parameters
#' 
#' @param sequences Input sequences, including M and W as methylated cytosines at upper and lower strands respectively
#' @param model Linear model containing position-based energy coefficients, e.g., 1CA, 1CG, 1CT, 1MW, 2CA, etc.
#' @param encoding Encoding scheme. By default 4L+1
#' @return Predicted binding energy values of \code{sequences} based on \code{model}
#' @examples
#' predictEnergyMW(sequences, model, encoding = "4L+1")
predictEnergyMW <- function(sequences, model, encoding = "4L+1") {
  #num <- (model$rank - 1) / 3  
  num <- nchar(sequences[[1]]) 
  
  #num shoulbe be the number of nucleotide positions used for prediction，model$rank is the number of coefficients in the models.
  # In some instances, there are NA values for some parameters, which are not included in total rank, so nchar(sequence) is used instead.
  
  tibble(ByteCode = SeqToByte(sequences, encoding = "4L+1")) %>%
    tidyr::separate(ByteCode,
                    into = paste0(rep(1:num, each = 4), c("CG", "CA", "CT", "MW")),
                    convert = TRUE
    ) %>%
    predict.lm(model, newdata = ., type = "response")
}