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
#' @import dplyr
#' @import ggfortify
#' @import ggplot2
#' @import pheatmap
#' @import reshape2
#' @import tibble
#' @import lazyeval
#' @import WGCNA
#' @import viridis
#' @import RColorBrewer
#' @import readr
#' @importFrom corrplot corrplot.mixed
#' @importFrom magrittr %>%
#' @importFrom purrr map
#' @importFrom rlang UQ sym syms
#' @importFrom tidyr complete nest unnest
#' @importFrom data.table setDT
#' @importFrom data.table transpose
#' @importFrom data.table :=
#' @importFrom data.table IDateTime
#' @importFrom purrr map2
#' 
#' @docType package
#' @name proBatch
if(getRversion() >= "2.15.1")  utils::globalVariables(c(".", 
                                                        "batch_size", "batch_the_same", "batch_total", "category", "dateTime", 
                                                        "fit", "label", "mean_fit", "median_batch", "median_global", 
                                                        "median_run", "optimise_bw", "optimise_df", "peptide_col_name", 
                                                        "sample_annotatation_col", "Step", "tipping.poings", "Var1", "Var2"))
NULL


#' @importFrom data.table, as.data.table
#' @importFrom data.table, fread
#' @importFrom data.table, dcast
#' @importFrom data.table, copy
#' @importFrom MASS lda


