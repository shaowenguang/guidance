#' DIA-guidance: A package to infer protein abundance from DIA/SWATH-MS 
#' peptide measurements
#'
#' The DIA-guidance package contains functions to normalize peptide/fragment-level intensities, 
#' estimate protein abundances and conduct statistical tests on protein measurements. Although 
#' the package has primarily been developed for mass spectrometry proteomics (DIA/SWATH),
#' it should also be applicable to most proteomics data with minor adaptations.
#'
#' To learn more about proBatch, start with the vignettes:
#' \code{browseVignettes(package = "guidance")}
#'
#' @importFrom data.table as.data.table
#' @importFrom data.table fread
#' @importFrom data.table dcast
#' @importFrom data.table copy
#' @importFrom MASS lda



