#' wrapper function
#' example
#'  dia_guidance(data= "S:/SWATH-guidance/feature_alignment.csv", 
#'                 sample_annotation="S:/SWATH-guidance/sample_annotation", 
#'                 level="PeptideIon") 



#' @export
dia_guidance <- function(data, sample_annotation, level){
  
  # import data
  peptideIons <- import_openswath(search_results= data, 
                                  sample_annotation = sample_annotation, 
                                  level=level) 
  
  # normalize data (if necessary)
  all_peptideIons <- long2wide(peptideIons)
  all_peptideIons_normalized <- normalize_data(all_peptideIons, replaceNA="keep", 
                                               normalization="none")
  
  # merge replicates
  cons_peptideIons <- merge_replicates(all_peptideIons_normalized, anno)
  
  # filter proteotypic peptides 
  cons_peptideIons <- cons_peptideIons[which(grepl("^1/", cons_peptideIons$ProteinName)), ]
  
  # calculate features 
  cons_peptideIons_features <- calc_features(cons_peptideIons)
  
  # calculate posterior probability of being RARE using LDA model 
  index_mean_int <- which(grepl("^mean_intensity", names(cons_peptideIons_features)))
  test <- perform_selection(cons_peptideIons_features)
  
  # filter by probability, impute missing value 
  test_yesFiltered <- test[prob > 0.2, ]
  test_yesFiltered_yesImputated <- imputate_missing_values(test_yesFiltered, c(3:17))
  
  # peptide to protein inference 
  cons_prot_test_yesFiltered_top3_sum_yesImputated_yesWeighted <- merge_replicates(
    pept2prot(test_yesFiltered_yesImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=T), anno)
  
  return(cons_prot_test_yesFiltered_top3_sum_yesImputated_yesWeighted)
}





