#' Visualize some energy matrix by plotting logo
#' 
#' @param energyMatrix
#' @param col_scheme The color scheme used to illustrate different nucleotides, which can be created by function ggseqlogo::make_col_scheme
#' @return A ggplot object illustrating the energy logo.
#' @examples
#' plotEnergyLogo(energyMatrix)
plotEnergyLogo <- function(energyMatrix, col_scheme = NA) {
  if(is.na(col_scheme))
    col_scheme = ggseqlogo::make_col_scheme(chars=c('A', 'C', 'G', 'T', 'M', 'W', 'N'),
                             cols=c("#0E927B", "#59A9D8", "#DC9514", "#1A1A1A", "#0001FC", "#FB0605", "#5C5C5D"))
  
  ggplot2::ggplot() +
    ggseqlogo::geom_logo(-energyMatrix, method = "custom", seq_type = "dna", col_scheme = col_scheme) +
    ggseqlogo::theme_logo() +
    ggplot2::xlab("Position") +
    ggplot2::ylab("-Energy (kT)") %>%
    return()
}