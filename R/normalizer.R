#' @importFrom data.table dcast copy

#' @export
long2wide <- function(input_dt) {

  if (!global_level %in% c("PeptideIon","Transition", "Peptide", "PeptideWithMod")) {
    stop("Please select a valid type for imported data to be kept in Prom. Options:  \"Transition\", \"PeptideIon(default)\", \"PeptideWithMod\", \"Peptide\"")
  }
  
  if(global_level=="PeptideIon") {
    
    if("InjectionName" %in% names(input_dt)) {
      
      #annotated data table
      
      #wide <- dcast(input_dt, PeptideIon + ProteinName ~ InjectionName, value.var="Intensity", fun.aggregate=sum) #wenguang: for duplicate peptideIon, but this will change na to 0.
      wide <- dcast(input_dt, PeptideIon + ProteinName ~ InjectionName, value.var=c("Intensity", "Score", "RT", "Width"), fun.aggregate=function(x) max(x, na.rm=T)) #wenguang:  for duplicate peptideIon, but this will change na to -inf, then replace_inf was applied to solve this...
    } else {
      wide <- dcast(input_dt, PeptideIon + ProteinName ~ filename, value.var=c("Intensity", "Score", "RT", "Width"), fun.aggregate=function(x) max(x, na.rm=T))
    }
  
  } else if(global_level=="Transition") {
      
    if("InjectionName" %in% names(input_dt)) {
      #annotated data table
      wide <- dcast(input_dt, aggr_Fragment_Annotation + PeptideIon + ProteinName ~ InjectionName, value.var=c("aggr_Peak_Area", "Score", "RT")) 
    } else {
      wide <- dcast(input_dt, aggr_Fragment_Annotation + PeptideIon + ProteinName ~ filename, value.var=c("aggr_Peak_Area", "Score", "RT"))
    }
    
  }
  
  wide <- replace_inf(wide)
  
  return(wide)
  
}

#' @export
summarize_data <- function(input_dt) {
  
  message("This table contains ", length(unique(input_dt$PeptideIon)), " peptide ions; " 
                                , length(unique(input_dt$ProteinName)), " proteins." )
  message("Among them, ", length(which(complete.cases(input_dt))), " are completed rows (no NAs)." )
  
}

#' @export
normalize_data <- function(input_dt, replaceNA="keep", normalization="mediancenter") {
  
  if (!replaceNA %in% c("remove","keep","zero","min_intensity")) {
    stop("Please select a valid method for NA replacement. Options:  \"remove\", \"keep\", \"zero\", \"min_intensity\"")
  }
  
  if (!normalization %in% c("mediancenter","quantile","TIC", "none", "iRT")) {
    stop("Please select a valid method for normalization. Options:  \"mediancenter\", \"quantile\", \"TIC\", \"none\", \"iRT\"")
  }
  
  if (replaceNA=="keep" && normalization=="quantile") {
    #warning("quantile normalization is not compatile with NA-contained data!! Now NAs will be replaced with min_intensity")
    #replaceNA <- "min_intensity"
  }
  
  normalized_dt <- copy(input_dt)

#  if( length(which( grepl("Intensity", names(normalized_dt))) ) > 0 ) { 
#    #wenguang: this means that input table is the matrix from Euler Portal
#    index_intensity <- which( grepl("Intensity", names(normalized_dt)) )
#  } else {
#    index_intensity <- which(!names(normalized_dt) %in% c("aggr_Fragment_Annotation", "PeptideIon", "ProteinName"))
#  }
 
  if( length(which( grepl("Intensity", names(normalized_dt))) ) > 0 ) { 
    #wenguang: peptide level data
    index_intensity <- which( grepl("Intensity", names(normalized_dt)) )
  } else {
    #wenguang: transition level data
    index_intensity <- which( grepl("aggr_Peak_Area", names(normalized_dt)) )
  }
   
#  if(log2transform) {
#    normalized_dt[, index_intensity := log2(normalized_dt[, index_intensity, with=F]), with=F]
#  }
  
  if(replaceNA=="remove") {
    normalized_dt <- normalized_dt[which(complete.cases(normalized_dt)), ]
  }  else if(replaceNA=="zero") {
    normalized_dt[is.na(normalized_dt)] <- 0
  }  else if(replaceNA=="min_intensity") {
    normalized_dt[is.na(normalized_dt)] <- min(normalized_dt[, index_intensity, with=F], na.rm = T)
  }
    
  
  if(normalization=="mediancenter") {
    
    #normalized_dt[, index_intensity := log2(normalized_dt[, index_intensity, with=F]), with=F]
    
    for(i in 1:length(index_intensity)) {
      normalized_dt[, names(normalized_dt)[index_intensity[i]] ] <- log2(normalized_dt[, index_intensity[i], with=F])
    }
    
    run_median <- sapply(normalized_dt[, index_intensity, with=F], median, na.rm=T)
    for(i in 1:length(index_intensity)) {
      #normalized_dt[, index_intensity[i] := normalized_dt[, index_intensity[i], with=F] - run_median[i] + median(run_median), with=F]
      normalized_dt[, names(normalized_dt)[index_intensity[i]]] <- normalized_dt[, index_intensity[i], with=F] - run_median[i] + median(run_median)
      normalized_dt[, names(normalized_dt)[index_intensity[i]]] <- 2^normalized_dt[, index_intensity[i], with=F]
    }
  } else if(normalization=="iRT") {  
    normalized_dt[, index_intensity := log2(normalized_dt[, index_intensity, with=F]), with=F]
    run_median <- sapply(normalized_dt[ProteinName=="1/iRT", index_intensity, with=F], median, na.rm=T)
    for(i in 1:length(index_intensity)) {
      #normalized_dt[, index_intensity[i] := normalized_dt[, index_intensity[i], with=F] - run_median[i] + median(run_median), with=F]
      normalized_dt[, names(normalized_dt)[index_intensity[i]]] <- normalized_dt[, index_intensity[i], with=F] - run_median[i] + median(run_median)
      normalized_dt[, names(normalized_dt)[index_intensity[i]]] <- 2^normalized_dt[, index_intensity[i], with=F]
    }
  } else if(normalization=="TIC") {
    run_TIC <- sapply(normalized_dt[, index_intensity, with=F], sum, na.rm=T)
    for(i in 1:length(index_intensity)) {
      normalized_dt[, names(normalized_dt)[index_intensity[i]]] <- normalized_dt[, index_intensity[i], with=F] / run_TIC[i] * median(run_TIC)
    }
  } else if(normalization=="quantile") {

    temp_rank <- apply(normalized_dt[, index_intensity, with=F], 2, function(x) rank(x, na.last = F))
    temp_sort <- apply(normalized_dt[, index_intensity, with=F], 2, function(x) sort(x, na.last = F))
    temp_mean <- apply(temp_sort, 1, mean_na)
    
    for(i in 1:length(index_intensity)) {
      normalized_dt[, names(normalized_dt)[index_intensity[i]]] <- temp_mean[temp_rank[,i]]
    }
    
  }
  
  normalized_dt[, numPerProt := length( get(names(normalized_dt)[1]) ), by=(ProteinName)]
  
  return(normalized_dt)  
  
}


#' @export
merge_replicates <- function(wide, sample_annotation = NULL, bool_NA_means_requant = FALSE, averageFun = "mean") {
  
  #annotated_long <- annotate_sample(search_result, sample_annotation)
  #annotated_wide <- long2wide(annotated_long)
  #annotated_wide_normalized <- normalize_data(annotated_wide, replaceNA="keep", normalization="mediancenter", log2transform=FALSE)
  
  cons_dt <- copy(wide)
  
  list_samples <- unique(sort(sample_annotation$SampleName)) 
  
  message("start to merge replicates...")
  
  for(i in 1:length(list_samples)) {
  
    message("processing: sample_", i, ": ", list_samples[i])
    
    
    if(global_level=="PeptideIon") {
      
      if(averageFun=="mean") {
        
        #wenguang: some short comments for using "mean" as the default function of averaging. 
        #Using HEM benchmarking datasets on the protein level, it seems that mean performs consistently better median, possible because median is good for outliers, which have already been removed at this stage.
        
        cons_dt[, paste0("mean_intensity_", list_samples[i]) := apply(.SD, 1, mean_na), 
              .SDcols = paste0("Intensity_", sample_annotation[sample_annotation$SampleName %in% list_samples[i], ]$InjectionName) ]
      } else if (averageFun=="median") {
        cons_dt[, paste0("mean_intensity_", list_samples[i]) := apply(.SD, 1, median_na), 
                .SDcols = paste0("Intensity_", sample_annotation[sample_annotation$SampleName %in% list_samples[i], ]$InjectionName) ]
      }
      cons_dt[,   paste0("cv_intensity_", list_samples[i]) := apply(.SD, 1, cv_na), 
              .SDcols = paste0("Intensity_", sample_annotation[sample_annotation$SampleName %in% list_samples[i], ]$InjectionName) ] 
      
    } else if(global_level=="Transition") {
      
      if(averageFun=="mean") {
        cons_dt[, paste0("mean_intensity_", list_samples[i]) := apply(.SD, 1, mean_na), 
              .SDcols = paste0("aggr_Peak_Area_", sample_annotation[sample_annotation$SampleName %in% list_samples[i], ]$InjectionName) ]
      } else if (averageFun=="median") {
        cons_dt[, paste0("mean_intensity_", list_samples[i]) := apply(.SD, 1, median_na), 
                .SDcols = paste0("aggr_Peak_Area_", sample_annotation[sample_annotation$SampleName %in% list_samples[i], ]$InjectionName) ]
      }
      cons_dt[,   paste0("cv_intensity_", list_samples[i]) := apply(.SD, 1, cv_na), 
              .SDcols = paste0("aggr_Peak_Area_", sample_annotation[sample_annotation$SampleName %in% list_samples[i], ]$InjectionName) ] 
      
    }
  
    
    if(bool_NA_means_requant==FALSE) {
      
      if(global_level=="PeptideIon") {
        cons_dt[, paste0("numNA_intensity_", list_samples[i]) := apply(.SD, 1, function(x) length(which(is.na(x)))), 
                .SDcols = paste0("Intensity_", sample_annotation[sample_annotation$SampleName %in% list_samples[i], ]$InjectionName) ]
      } else if(global_level=="Transition") {
        cons_dt[, paste0("numNA_intensity_", list_samples[i]) := apply(.SD, 1, function(x) length(which(is.na(x)))), 
                .SDcols = paste0("aggr_Peak_Area_", sample_annotation[sample_annotation$SampleName %in% list_samples[i], ]$InjectionName) ]
      }
      
    } else {
      
      #wenguang: this will use the number of requant values (m_score==2) as the number of NAs. 
      cons_dt[, paste0("numNA_intensity_", list_samples[i]) := apply(.SD, 1, function(x) length(which(x == 2.0))), 
              .SDcols = paste0("Score_", sample_annotation[sample_annotation$SampleName %in% list_samples[i], ]$InjectionName) ]
      
    }
         
  }
  
  message("done with merging replicates...")
  
  
  
  for(i in 1:length(list_samples)) {
  
#    message("Globally, removing ", length(which(is.na(cons_dt[, paste0("mean_intensity_", list_samples[i]), with=F ]))), " empty features existing in sample ", list_samples[i])
    if(length(which(is.na(cons_dt[, paste0("mean_intensity_", list_samples[i]), with=F ])))) {
#      cons_dt <- cons_dt[-which(is.na(cons_dt[, paste0("mean_intensity_", list_samples[i]), with=F ])), ]
    }
  }
  
  #cons_dt[, numPerProt := length( get(names(cons_dt)[1]) ), by=(ProteinName)]
  
  return(cons_dt)
  
}





#' @export
keep_proteotypic_only <- function(input_dt) {

  output_dt <- copy(input_dt)

  output_dt <- output_dt[which(grepl("^1/", output_dt$ProteinName)), ]
  
  return(output_dt)

}
