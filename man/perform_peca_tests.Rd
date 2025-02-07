% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statisticalTest.R
\name{perform_peca_tests}
\alias{perform_peca_tests}
\title{Perform PECA test}
\usage{
perform_peca_tests(input_file, sample_annotation, input_test = "modt",
  input_bool_paired = FALSE)
}
\arguments{
\item{input_file}{data table or data frame in wide representation. The data typically 
contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
protein measurements in rows. The input data is usually an output from 
\code{pept2prot()} and \code{merge_replicates()}}

\item{sample_annotation}{data matrix with \code{SampleName}, biological covariates 
(biological replicates) and technical covariates (technical replicates, batches, etc)}

\item{input_test}{the type of t-test to be conducted either ordinary 
\code{"t"} or modified \code{"modt"} t-test}

\item{input_bool_paired}{a logical indicating whether a paired test is performed}
}
\description{
A function to perform probe-level expression change averaging (PECA)
to identify expression change between two groups of samples for each protein.
}
