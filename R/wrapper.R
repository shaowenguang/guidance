#' wrapper function
#' 
#' @example
#'  prot_table <- dia_guidance(data= "S:/SWATH-guidance/feature_alignment.csv", 
#'                 sample_annotation="S:/SWATH-guidance/sample_annotation", 
#'                 level="PeptideIon") 
#'                 
#' @export
dia_guidance <- function(data = NULL, data_fromEuler = NULL, sample_annotation = NULL, level = "PeptideIon", 
                         replaceNA="keep", bool_NA_means_requant = FALSE, 
                         averageFun = "mean",
                         normalization="mediancenter", filter_prob = 0.2, 
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





