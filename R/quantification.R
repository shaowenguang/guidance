#' @export
generate_protein_table <- function(input_dt, input_rank_index = "prob", topN = 3, aggfun = "sum", bool_weighted_by_prob = TRUE, bool_imputation = TRUE,
                                   prob_threshold = 0.2, bool_keep_low_confident_prot = FALSE) {
  
  
  test_filtered <- input_dt[prob > prob_threshold, ]
  
  if(bool_imputation == TRUE) {
    
    index_intensity <- which(grepl("^Intensity", names(test)))
    
    test_filtered_imputated <- imputate_missing_values(test_filtered, index_intensity)
    
    cons_prot_table <- merge_replicates(pept2prot(test_filtered_imputated, input_rank_index, topN, aggfun, bool_weighted_by_prob), anno)
    
  } else{
    
    cons_prot_table <- merge_replicates(pept2prot(test_filtered, input_rank_index, topN, aggfun, bool_weighted_by_prob), anno)
    
  }
  
  
  if(bool_keep_low_confident_prot == TRUE) {

    low_confident_prot_list <- setdiff(unique(input_dt$ProteinName), unique(test_filtered$ProteinName))
    lowSignal <- input_dt[which(input_dt$ProteinName %in% low_confident_prot_list), ]
    
    if(bool_imputation == TRUE) {
      
      lowSignal_imputated <- imputate_missing_values(lowSignal, index_intensity)
      cons_prot_table_lowSignal <- merge_replicates(pept2prot(lowSignal_imputated, input_rank_index, topN, aggfun, bool_weighted_by_prob), anno)
      
    } else{
      
      cons_prot_table_lowSignal <- merge_replicates(pept2prot(lowSignal, input_rank_index, topN, aggfun, bool_weighted_by_prob), anno)
      
    }
    
    output_prot_table <- rbind(cons_prot_table, cons_prot_table_lowSignal)
    
  } else {
    
    output_prot_table <- copy(cons_prot_table)
    
  }
  
  return(output_prot_table)
  
  
}






#' @export
pept2prot <- function(input_dt, input_rank_index = "prob", topN = 3, aggfun = "sum", bool_weighted_by_prob = TRUE) {
  
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
      long_combined <- long[, .(Quant = as.numeric(sum_na(Intensity * prob) / sum_na(prob))), by=.(ProteinName, numPerProt, protQuantProb, run_id)]
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




pept2prot_log2 <- function(input_dt, input_index, topN) {
  
  select <- copy(input_dt)
  
  for(i in 1:dim(anno)[1]) {
    #    select[which(is.na(select[, anno$Injection[i], with=F])), anno$Injection[i] := get(paste0("mean_intensity_", anno$SampleName[i])) ]
  }
  
  select <- select[, rank := rank( - as.numeric(get(input_index)), na.last = T), by=ProteinName]
  
  select <- select[rank <= topN, ]
  
  select[, numPerProt := length( get(names(select)[1]) ), by=(ProteinName)]
  select[, protQuantProb := as.numeric(mean(prob)), by=ProteinName]
  
  #select[, ProbPerProt = as.numeric(mean(Prob)), by=ProteinName]
  
  long <- melt(select, id.vars=c(names(input_dt)[1], "ProteinName", "numPerProt", "protQuantProb"), measure.vars=anno$InjectionName, variable.name="run_id", value.name="Intensity")
  
  long_combined <- long[, .(Quant = as.numeric(mean_na(log2(Intensity)))), by=.(ProteinName, numPerProt, protQuantProb, run_id)] #wenguang: here, as.numeric is necessary, otherwise errors will occur "Column 1 of result for group 2 is type 'logical' but expecting type 'double'. Column types must be consistent for each group." This is because mean_na will return a number or NA, and they belong to two different classes...
  
  combined <- dcast(long_combined, ProteinName + numPerProt + protQuantProb ~ run_id, value.var="Quant")
  
  return(combined)
  
}



#' @export
imputate_missing_values <- function(input_dt, input_index) {
  
  message("start to imputate missing values...")
  
  output_dt <- copy(input_dt)
  
  for(i in 1:dim(output_dt)[1]) {
    
    #wenguang: this is a very slow way to imputate, needs to be improved using data.table functions.
    
    if( (i %% 1000) == 0 ) {
      message("imputating missing values for ",  (i %/% 1000), "000 out of ", dim(output_dt)[1], " peptide ions...")
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
  
  message("done with imputating missing values...")
  
  return(output_dt)
  
}
