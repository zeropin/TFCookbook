#' Build energy model based on Sequence-Energy relationship
#' 
#' @param data Data frame containing at least two variables -- Sequence and Energy
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


#' Build energy model based on Sequence-Energy relationship
#' 
#' @param data Data frame containing at least two variables -- Sequence and Energy
#' @return Linear model containing position-based energy coefficients, e.g., 1CA, 1CG, 1CT, 2CA, etc.
#' @examples
#' buildEnergyModel(data)
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
#' 
#' @param sequences Input sequences
#' @param model Linear model containing position-based energy coefficients, e.g., 1CA, 1CG, 1CT, 2CA, etc.
#' @return Predicted binding energy values of \code{sequences} based on \code{model}
#' @examples
#' predictEnergy(sequences, model)
predictEnergy <- function(sequences, model) {
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
}


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