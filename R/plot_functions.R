plotEnergyLogo <- function(energyMatrix) {
  ggplot2::ggplot() +
    ggseqlogo::geom_logo(-energyMatrix, method = "custom", seq_type = "dna") +
    ggseqlogo::theme_logo() +
    ggplot2::xlab("Position") +
    ggplot2::ylab("-Energy (kT)") %>%
    return()
}
