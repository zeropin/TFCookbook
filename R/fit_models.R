buildEnergyModel <- function(data) {
  num <- nchar(data$Sequence[1])

  data %>%
    dplyr::mutate(ByteCode = SeqToByte(Sequence)) %>%
    dplyr::select(Energy, ByteCode) %>%
    tidyr::separate(ByteCode,
      into = paste0(rep(1:num, each = 3), c("CG", "CA", "CT")),
      convert = TRUE
    ) %>%
    lm(Energy ~ . - Energy, data = .) %>%
    return()
}

buildEnergyModel_MW <- function(data) {
  num <- nchar(data$Sequence[1])

  data %>%
    dplyr::mutate(ByteCode = SeqToByte_MW(Sequence)) %>%
    dplyr::select(Energy, ByteCode) %>%
    tidyr::separate(ByteCode,
      into = paste0(rep(1:num, each = 5), c("CG", "CA", "CT", "CM", "GW")),
      convert = TRUE
    ) %>%
    lm(Energy ~ . - Energy, data = .) %>%
    return()
}

predictEnergy <- function(sequences, model) {
  num <- (model$rank - 1) / 3

  tibble(ByteCode = SeqToByte(sequences)) %>%
    tidyr::separate(ByteCode,
      into = paste0(rep(1:num, each = 3), c("CG", "CA", "CT")),
      convert = TRUE
    ) %>%
    predict.lm(model, newdata = ., type = "response")
}
