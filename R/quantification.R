#' Infer protein abundance from peptide measurements  
#' 
#' @description A function infer protein abundance from peptide/fragment-level 
#' measurements. The users can denote parameters such as the number of peptides 
#' to be utilized for protein inference (\code{topN}), the method to aggregate 
#' peptide/fragment intensities (\code{aggfun}) and whether to weight intensity 
#' of each peptide by their posterior probability in estimating protein abundance 
#' (\code{bool_weighted_by_prob}). This function then outputs estimated protein 
#' intensities for each sample. 
#' 
#' @param input_dt data table or data frame in wide representation. The data typically 
#' contains \code{PeptideIon}, \code{ProteinName} and sample names in columns and 
#' measurements of each peptide or precursor ions in rows
#' @param input_rank_index name of a column to rank the protein by. The default 
#' column is \code{"prob} which depicts the probability of being peptide 
#' representative of protein 
#' @param topN number of peptides utilized to infer protein abundance for each protein 
#' @param aggfun method to aggregate peptide measurements to estimate protein abundance.
#' Options include \code{"mean"} and \code{"sum"}
#' @param bool_weighted_by_prob boolean value (\code{TRUE} or \code{FALSE}) determining 
#' whether to weight the intensity by the posterior probability of being a representative 
#' peptide
#'
#' @export
#' 
#' @examples 
#' peptideIons_features <- calc_features(all_peptideIons)
#' peptideIons_features_select <- perform_selection(peptideIons_features)
#' 
#' peptide_to_protein <- pept2prot(peptideIons_features_select, 
#' "prob", 3, aggfun="sum", bool_weighted_by_prob=T)
#' 
pept2prot <- function(input_dt, input_rank_index = "prob", 
                      topN = 3, aggfun = "mean", bool_weighted_by_prob = TRUE) {
  
  select <- copy(input_dt)
  
  for(i in 1:dim(anno)[1]) {
#    select[which(is.na(select[, anno$Injection[i], with=F])), anno$Injection[i] := get(paste0("mean_intensity_", anno$SampleName[i])) ]
  }
  
  select <- select[, rank := rank( - as.numeric(get(input_rank_index)), na.last = T), by=ProteinName]
  
  select <- select[rank <= topN, ]
  
  select[, numPerProt := length( get(names(select)[1]) ), by=ProteinName]
  select[, protQuantProb := as.numeric(mean(prob)), by=ProteinName]
  
  #select[, ProbPerProt = as.numeric(mean(Prob)), by=ProteinName]
  
  long <- melt(select, id.vars = c(names(input_dt)[1], "ProteinName", "numPerProt", "protQuantProb", "prob")
                     #, measure.vars = paste0("Intensity_", anno$InjectionName)
                     , measure.vars = names(select)[which(grepl("^Intensity_", names(select)))]
                     , variable.name = "run_id"
                     , value.name = "Intensity" )
  
  if(aggfun=="mean") {
  #wenguang: here, as.numeric is necessary, otherwise errors will occur "Column 1 of result for group 2 is type 'logical' but expecting type 'double'. Column types must be consistent for each group." This is because mean_na will return a number or NA, and they belong to two different classes...
    if(bool_weighted_by_prob == TRUE) {
      long_combined <- long[, .(Quant = as.numeric(sum_na(Intensity * prob) / sum_na( sign(Intensity) * prob) )), by=.(ProteinName, numPerProt, protQuantProb, run_id)]
    } else {
      long_combined <- long[, .(Quant = as.numeric(mean_na(Intensity))), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    }
    
    
  } else if (aggfun=="sum") {
    #wenguang: please note that using sum function, na will be automatically replaced with zero.
    #long_combined <- long[, .(Quant = as.numeric(sum_na(Intensity))), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    #long_combined <- long[, .(Quant = as.numeric(2^mean_na(log2(Intensity)))), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    #long_combined <- long[, .(Quant = as.numeric(sum_na(Intensity * prob) / sum_na(prob))), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    
    if(bool_weighted_by_prob == TRUE) {
      long_combined <- long[, .(Quant = as.numeric( sum_na(Intensity * prob) )), by=.(ProteinName, numPerProt, protQuantProb, run_id)]
    } else {
      long_combined <- long[, .(Quant = as.numeric( sum_na(Intensity) )), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    }

  }

  combined <- dcast(long_combined, ProteinName + numPerProt + protQuantProb ~ run_id, value.var="Quant")
  
  return(combined)
  
}




#' Log2 tranform peptide intensity measurements and infer protein abunance 
#' 
#' @description A function infer protein abundance from peptide/fragment-level 
#' measurements after log2 tranformation. The users can denote parameters 
#' such as the number of peptides to be utilized for protein inference 
#' (\code{topN}), the method to aggregate peptide/fragment intensities 
#' (\code{aggfun}) and whether to weight intensity of each peptide by 
#' their posterior probability in estimating protein abundance
#' (\code{bool_weighted_by_prob}). This function then outputs estimated protein 
#' intensities for each sample. 
#' 
#' @param input_dt data table or data frame in wide representation. The data typically 
#' contains \code{PeptideIon}, \code{ProteinName} and sample names in columns and 
#' measurements of each peptide or precursor ions in rows
#' @param input_rank_index name of a column to rank the protein by. The default 
#' column is \code{"prob} which depicts the probability of being peptide 
#' representative of protein 
#' @param topN number of peptides utilized to infer protein abundance for each protein 
#' @param aggfun method to aggregate peptide measurements to estimate protein abundance.
#' Options include \code{"mean"} and \code{"sum"}
#' @param bool_weighted_by_prob boolean value (\code{TRUE} or \code{FALSE}) determining 
#' whether to weight the intensity by the posterior probability of being a representative 
#' peptide
#'
#' @export
#' 
#' @examples 
#' peptideIons_features <- calc_features(all_peptideIons)
#' peptideIons_features_select <- perform_selection(peptideIons_features)
#' 
#' peptide_to_protein <- pept2prot(peptideIons_features_select, 
#' "prob", 3, aggfun="sum", bool_weighted_by_prob=T)
#' 
pept2prot_log2 <- function(input_dt, input_rank_index = "prob", 
                           topN = 3, aggfun = "mean", bool_weighted_by_prob = TRUE) {
  
  select <- copy(input_dt)
  
  for(i in 1:dim(anno)[1]) {
    #    select[which(is.na(select[, anno$Injection[i], with=F])), anno$Injection[i] := get(paste0("mean_intensity_", anno$SampleName[i])) ]
  }
  
  select <- select[, rank := rank( - as.numeric(get(input_rank_index)), na.last = T), by=ProteinName]
  
  select <- select[rank <= topN, ]
  
  select[, numPerProt := length( get(names(select)[1]) ), by=ProteinName]
  select[, protQuantProb := as.numeric(mean(prob)), by=ProteinName]
  
  #select[, ProbPerProt = as.numeric(mean(Prob)), by=ProteinName]
  
  long <- melt(select, id.vars = c(names(input_dt)[1], "ProteinName", "numPerProt", "protQuantProb", "prob")
               #, measure.vars = paste0("Intensity_", anno$InjectionName)
               , measure.vars = names(select)[which(grepl("^Intensity_", names(select)))]
               , variable.name = "run_id"
               , value.name = "Intensity" )
  
  long$Intensity <- as.numeric(mean_na(log2(long$Intensity + 1)))
  
  if(aggfun=="mean") {
    #wenguang: here, as.numeric is necessary, otherwise errors will occur "Column 1 of result for group 2 is type 'logical' but expecting type 'double'. Column types must be consistent for each group." This is because mean_na will return a number or NA, and they belong to two different classes...
    if(bool_weighted_by_prob == TRUE) {
      long_combined <- long[, .(Quant = as.numeric(sum_na(Intensity * prob) / sum_na( sign(Intensity) * prob) )), by=.(ProteinName, numPerProt, protQuantProb, run_id)]
    } else {
      long_combined <- long[, .(Quant = as.numeric(mean_na(Intensity))), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    }
    
    
  } else if (aggfun=="sum") {
    #wenguang: please note that using sum function, na will be automatically replaced with zero.
    #long_combined <- long[, .(Quant = as.numeric(sum_na(Intensity))), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    #long_combined <- long[, .(Quant = as.numeric(2^mean_na(log2(Intensity)))), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    #long_combined <- long[, .(Quant = as.numeric(sum_na(Intensity * prob) / sum_na(prob))), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    
    if(bool_weighted_by_prob == TRUE) {
      long_combined <- long[, .(Quant = as.numeric( sum_na(Intensity * prob) )), by=.(ProteinName, numPerProt, protQuantProb, run_id)]
    } else {
      long_combined <- long[, .(Quant = as.numeric( sum_na(Intensity) )), by=.(ProteinName, numPerProt, protQuantProb, run_id)] 
    }
    
  }
  
  combined <- dcast(long_combined, ProteinName + numPerProt + protQuantProb ~ run_id, value.var="Quant")

  return(combined)
  
}



#' Impute missing values 
#' 
#' @description A function to impute for missing value by fitting data to a 
#' uniform distribution between minimum $0.1 \times CV$ and minimum 
#' $1.1 \times CV$ and infer the intensity of missing values. 
#' 
#' @param input_dt data table or data frame in wide representation. The data typically 
#' contains \code{PeptideIon}, \code{ProteinName} and sample names in columns and 
#' measurements of each peptide or precursor ions in rows. 
#' @param input_index a integer vector denoting numeric columns to apply imputation
#'
#' @return  data.table data.frame  
#' 
#' @export
#' 
#' @examples 
#' peptideIons_features <- calc_features(all_peptideIons)
#' peptideIons_features_select <- calc_features(peptideIons_features)
#' Imputated <- impute_missing_values(test_yesFiltered, c(3:17))
#' 
impute_missing_values <- function(input_dt, input_index) {
  
  message("It starts to impute missing values...")
  
  output_dt <- copy(input_dt)
  
  for(i in 1:dim(output_dt)[1]) {
    
    #wenguang: this is a very slow way to imputate, needs to be improved using data.table functions.
    
    if( (i %% 1000) == 0 ) {
      message("  imputing missing values for ",  (i %/% 1000), "000 out of ", dim(output_dt)[1], " peptide ions...")
    }
    
    index_na <- which(is.na(output_dt[i, input_index, with=F]))
    if(length(index_na) > 0) {
      #wenguang: this will imputate a value following the uniform distribution between min_per_peptide * (0.1*cv_per_peptide, 1.1*cv_per_peptide).
      imputated_values <- (1 - runif(length(index_na), -0.1, 0.9) * sd(output_dt[i, input_index[-index_na], with=F])/mean(unlist(output_dt[i, input_index[-index_na], with=F])) ) * min(output_dt[i, input_index[-index_na], with=F])
      #imputated_values <- (0.6 - runif(length(index_na), 0.2, 0.5) * sd(output_dt[i, input_index[-index_na], with=F])/mean(unlist(output_dt[i, input_index[-index_na], with=F])) ) * min(output_dt[i, input_index[-index_na], with=F])
      
      
      for(j in 1:length(index_na)) {
        
        set(output_dt, i, input_index[index_na[j]], imputated_values[j])
        
      }
      
    }
    
  }
  
  message("Done with imputating missing values...")
  
  return(output_dt)
  
}
