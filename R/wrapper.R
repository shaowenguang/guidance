#' Wrapper function from importing peptide data to infering protein abundance
#' 
#' @param data A data frame containing the SWATH-MS data. This data typically
#' contains peptide precursors in each row with corresponding \code{ProteinName}, 
#' \code{Intensity}, \code{RT}, \code{Score} and etc in columns. The data can be loaded
#' from local directory.
#' @param data_fromEuler A data frame containing the SWATH-MS data. This data typically
#' contains peptide ions in each row with corresponding \code{ProteinName}, 
#' \code{Intensity}, \code{RT}, \code{Score} and etc in columns. The data can be loaded
#' from local directory.
#' @param sample_annotation data matrix with \code{SampleName}, biological covariates 
#' (biological replicates) and technical covariates (technical replicates, batches, etc)
#' @param level the protein-level of imported data. The options include 
#' \code{"PeptideIon"} (by default), \code{"Transition"}, \code{"Peptide"} or 
#' \code{"PeptideWithMod"}
#' @param replaceNA whether to treat missing values. The options include to \code{"remove"}, 
#' \code{"keep"}, replace them with \code{"zero"}, or minimum intensity (\code{"min_intensity"})
#' @param bool_NA_means_requant boolean value (\code{TRUE} or \code{FALSE}) determining if 
#' the missing values correspond to requants. \code{bool_NA_means_requant = TRUE} will 
#' use the number of requant values (m_score = 2) as the number of NAs.
#' @param averageFun method to compute mean peptide or precursor ion intensity 
#' of biological replicates. Options include \code{"mean"} and \code{"median"}. 
#' @param normalization different methods of normalization. The options include 
#' median-centering (\code{"mediancenter"}), quantile normalization (\code{"quantile"},
#' and normalized based on total ion current (\code{"TIC"}) and indexed retention 
#' time (iRT) standards \code{"iRT"}. The default is \code{"mediancenter"} and denote 
#' \code{"none"} if normalization is not necessary.
#' @param filter_prob a numeric value in range of 0 to 1 denoting posterior 
#' probability threshold to filter peptides by.
#' @param input_rank_index 
#' @param topN number of peptides utilized to infer protein abundance for each protein 
#' @param aggfun method to aggregate peptide measurements to estimate protein abundance.
#' Options include \code{"mean"} and \code{"sum"}
#' @param bool_weighted_by_prob boolean value (\code{TRUE} or \code{FALSE}) determining 
#' whether to weight the intensity by the posterior probability of being a representative 
#' peptide
#' @param bool.removeDecoy boolean value (\code{TRUE} or \code{FALSE}) determining 
#' if the decoy peptides should be removed from the imported data 
#' @param remove_prefixInFileName boolean value (\code{TRUE} or \code{FALSE}) whether to 
#' remove unnecessary prefix from euler portal file name. For example, a typical 
#' filename will look like \code{"/scratch/71239421.tmpdir/xuep_J180621_SW_3.mzXML.gz"} 
#' and \code{remove_prefixInFileName = TRUE} will result in \code{"xuep_J180621_SW_3.mzXML.gz"}
#'
#' @example
#' donprot_table <- dia_guidance(data= "S:/SWATH-guidance/feature_alignment.csv", 
#'                 sample_annotation="S:/SWATH-guidance/sample_annotation", 
#'                 level="PeptideIon") 
#'                 
#' @export
dia_guidance <- function(data = NULL, data_fromEuler = NULL, sample_annotation = NULL, level = "PeptideIon", 
                         replaceNA="keep", bool_NA_means_requant = FALSE, 
                         averageFun = "mean",
                         normalization="mediancenter", filter_prob = 0.25, 
                         input_rank_index = "prob", topN = 3, aggfun = "sum", 
                         bool_weighted_by_prob = TRUE, bool.removeDecoy = T, 
                         remove_prefixInFileName = FALSE){
  
  if(!is.null(data)){
    peptide <- import_openswath(search_results= data, 
                                    sample_annotation = sample_annotation, 
                                    level=level) 
    all_peptide <- long2wide(peptide)
  }
  
  if(!is.null(data_fromEuler)){
    all_peptide <- import_openswath_matrix_fromEulerPortal(
      search_results=data_fromEuler, sample_annotation=sample_annotation) 
  }

  if(!is.null(data) & !is.null(data_fromEuler)){
    stop("Input one dataset only")
  }
  
  # normalize data (if necessary)
  all_peptide_normalized <- normalize_data(all_peptide, replaceNA= replaceNA, 
                                               normalization=normalization)
  
  # merge replicates - CHECK!! THE STEP, SAMPLENAME NOT VALID 
  cons_peptide <- merge_replicates(wide = all_peptide_normalized, sample_annotation = anno, 
                                       bool_NA_means_requant = bool_NA_means_requant, averageFun = averageFun)

  # filter proteotypic peptides 
  cons_peptide <- cons_peptide[which(grepl("^1/", cons_peptide$ProteinName)), ]
  
  # calculate features 
  cons_peptide_features <- calc_features(cons_peptide)
  
  # calculate posterior probability of being RARE using LDA model 
  index_mean_int <- which(grepl("^mean_intensity", names(cons_peptide_features)))
  test <- perform_selection(cons_peptide_features)
  
  # filter by probability, impute missing value 
  test_yesFiltered <- test[prob > filter_prob, ]

  columns <- which(grepl("^Intensity_", colnames(test_yesFiltered)))
  test_yesFiltered_yesImputated <- impute_missing_values(test_yesFiltered, columns)
  
  # peptide to protein inference 
  prot_inf_table <- merge_replicates(
    pept2prot(test_yesFiltered_yesImputated, input_rank_index = input_rank_index, 
              topN = topN, aggfun=aggfun, bool_weighted_by_prob=bool_weighted_by_prob), anno)
  
  return(prot_inf_table)

}



#' @export
prepare_matrix <- function(data, data_type = "openswath"
                               , sample_annotation = NULL
                               , level = "PeptideIon"
                               , bool.keepProteotypic = T
                               , normalization="mediancenter"
                               , replaceNA="keep"
                               , remove_prefixInFileName = FALSE
                               , bool.removeDecoy = T) {

if(data_type=="openswath") {
  peptideIons <- import_openswath(search_results=data, sample_annotation=sample_annotation, level=level
                                                     , bool.removeDecoy=bool.removeDecoy, remove_prefixInFileName=remove_prefixInFileName) 
} else if(data_type=="openswath_fromEulerPortal") {
  peptideIons <- import_openswath_matrix_fromEulerPortal(search_results=data, sample_annotation=sample_annotation)
} else if(data_type=="spectronaut") {
  peptideIons <- import_spectronaut_matrix(search_results=data, sample_annotation=sample_annotation)
} else {
  stop("Please select a valid type for imported data to be kept in Prom. Options:  \"openswath(default)\", \"spectronaut\", \"openswath_fromEulerPortal\"")
}

  all_peptideIons <- long2wide(peptideIons)

  all_peptideIons_normalized <- normalize_data(all_peptideIons, replaceNA=replaceNA, normalization=normalization)

  cons_peptideIons <- merge_replicates(all_peptideIons_normalized, anno)
   
  if(bool.keepProteotypic == T) {
    cons_peptideIons <- keep_proteotypic_only(cons_peptideIons)
  }

  cons_peptideIons_features <- calc_features(cons_peptideIons)

  return(cons_peptideIons_features)

}

