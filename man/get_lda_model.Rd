% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selecter.R
\name{get_lda_model}
\alias{get_lda_model}
\title{Generate linear discriminant analysis (LDA) model}
\usage{
get_lda_model(input_dt, input_features)
}
\arguments{
\item{input_dt}{data table or data frame in wide representation. The data typically 
contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
measurements of each peptide or precursor ions in rows.}

\item{input_features}{a vector of features to be used for training. The examples 
include \code{scaled_mean_intensity_all}, \code{scaled_cv_intensity_all}, 
\code{scaled_numNA_intensity_all}, \code{scaled_averaged_score_all}, 
\code{scaled_median_PCC}, \code{scaled_sd_width_all} and \code{label}.}
}
\description{
A function to generate linear discriminant analysis (LDA) model 
by using a subset of feature (peptide) statistics. This function utilizes 
\code{lda()} function in \code{MASS}.
}
\examples{
global_level = "PeptideIon"
d_feature <- calc_features(peptideIon_st)

ecoli_std <- c(2, 3, 4, 6, 8)
index_mean_int <- which(grepl("^mean_intensity", names(d_feature)))
d_feature <- d_feature[, cor_std := 0]
d_feature$cor_std <- apply(d_feature[, index_mean_int, with=F], 1,
 function(x) cor(x, ecoli_std, use="p"))
d_feature$cor_std[apply(d_feature[, index_mean_int, with=F], 1, 
 function(x) count_pairwise_number(x, ecoli_std)) < 4] <- NA
d_feature[, label := "bad"]
d_feature[ cor_std > 0.95, ]$label <- "good"

index_features <- c("scaled_mean_intensity_all", "scaled_cv_intensity_all", 
                    "scaled_numNA_intensity_all", "scaled_averaged_score_all", 
                     "scaled_median_PCC", "scaled_sd_width_all", "label")
model_lda <- get_lda_model(d_feature, index_features)

}
