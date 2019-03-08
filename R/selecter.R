#' @importFrom MASS lda

#Wenguang: the loop in this function is slow, the results of which, however, were well validated. An updated/faster version to calculate features was implemented below...
#' Title
#'
#' @param input_dt 
#' @param level 
#'
#' @return
#' @export
#'
#' @examples
calc_features_old_just_backup <- function(input_dt, level="PeptideIon") {
  
  message("start to calculate features per proteins...")
  
  output_dt <- copy(input_dt)
  
  output_dt[, feature_mean_intensity_all := apply(.SD, 1, mean_na), .SDcols=names(output_dt)[which(grepl("^mean_intensity", names(output_dt)))]]
  output_dt[, feature_cv_intensity_all := apply(.SD, 1, mean_na), .SDcols=names(output_dt)[which(grepl("^cv_intensity", names(output_dt)))]]
  output_dt[, feature_numNA_intensity_all := apply(.SD, 1, mean_na), .SDcols=names(output_dt)[which(grepl("^numNA_intensity", names(output_dt)))]]
  output_dt[, feature_median_PCC := -2]
  output_dt[, feature_median_SCC := -2]
  output_dt[, feature_MAD_dist := -2]
  
  #if cv is na (that means no replicates), then the median cv overall is assigned. 
  output_dt[is.na(output_dt$feature_cv_intensity_all), ]$feature_cv_intensity_all <- median_na(output_dt$feature_cv_intensity_all)
  
  output_dt[, median_PCC_PerProt := -2]
  
  
  prot_list <- unique(output_dt$ProteinName)
  
  for (i in 1:length(prot_list)) {
    
    if(i%%500 == 0) {
      message("processing: ", i, "th proteins: ", prot_list[i])
    }
    
    case <- output_dt[output_dt$ProteinName==prot_list[i], ]
    
    if(dim(case)[1] > 1) {
      
      index_injections <- which(names(output_dt) %in% anno$InjectionName)
      
      mat_pcc <- cor(t(case[, index_injections, with=F]), use="p", method="p")
      mat_scc <- cor(t(case[, index_injections, with=F]), use="p", method="s")
      
      # To exclude the self pairwise comparison. Otherwise median_PCC of proteins with two peptides will be overestimated (because of at least there is a ONE between two compairons).
      diag(mat_pcc) <- NA
      diag(mat_scc) <- NA
      
      # at least three completed paires to calculate the cor
      mat_pairwise_complete <-  count_pairwise_number_matrix(t(case[, index_injections, with=F]))
      
      mat_pcc[mat_pairwise_complete < 3] <- NA
      mat_scc[mat_pairwise_complete < 3] <- NA
      
      output_dt[which(output_dt$PeptideIon %in% case$PeptideIon), ]$feature_median_PCC <- apply(mat_pcc, 1, function(x) median_na(x))
      output_dt[which(output_dt$PeptideIon %in% case$PeptideIon), ]$feature_median_SCC <- apply(mat_scc, 1, function(x) median_na(x))
      
      output_dt[which(output_dt$PeptideIon %in% case$PeptideIon), ]$median_PCC_PerProt <- median_na(apply(mat_pcc, 1, function(x) median_na(x)))
      
      #ratio_case <- log2(case$mean_intensity_a / case$mean_intensity_b)
      #output_dt[which(output_dt$PeptideIon %in% case$PeptideIon), ]$feature_MAD_dist <- abs(ratio_case - median(ratio_case)) / mad(ratio_case, constant=1)
      
    }
    
  }
  
  
  
  output_dt[, `:=`(scaled_mean_intensity_all = scale(log2(feature_mean_intensity_all))
                  , scaled_cv_intensity_all = scale(feature_cv_intensity_all)
                  , scaled_numNA_intensity_all = scale(feature_numNA_intensity_all)
                  , scaled_median_PCC = scale(feature_median_PCC)
                  , scaled_MAD_dist = scale(feature_MAD_dist) )]
  
  
  #if median_PCC/median_SCC is na (that means many NAs), then the median median_PCC/median_SCC overall is assigned. 
  output_dt[is.na(output_dt$scaled_median_PCC), ]$scaled_median_PCC <- median_na(output_dt$scaled_median_PCC)

  
  message("done with calculating features...")
  
  return(output_dt)
  
}





#' @param input_dt 
#'
#' @export
calc_features <- function(input_dt) {
  
  message("start to calculate features per protein...")

  output_dt <- copy(input_dt)
  
  output_dt[, feature_mean_intensity_all := apply(.SD, 1, mean_na), .SDcols=names(output_dt)[which(grepl("^mean_intensity", names(output_dt)))]]
  output_dt[, feature_cv_intensity_all := apply(.SD, 1, mean_na), .SDcols=names(output_dt)[which(grepl("^cv_intensity", names(output_dt)))]]
  output_dt[, feature_numNA_intensity_all := apply(.SD, 1, mean_na), .SDcols=names(output_dt)[which(grepl("^numNA_intensity", names(output_dt)))]]

  #output_dt[, feature_averaged_score_all := apply(.SD, 1, function(x) mean_na( - log2(x)) ), .SDcols=names(output_dt)[which(grepl("^Score_", names(output_dt)))]]
  #wenguang: some extreme low-score (high confident) cases could be score==zero, resulting in INF after log2. Thus, a tiny number was constantly added...
  output_dt[, feature_averaged_score_all := apply(.SD, 1, function(x) mean_na( - log2(x + 4.744947e-16)) ), .SDcols=names(output_dt)[which(grepl("^Score_", names(output_dt)))]]
  
  output_dt[, feature_sd_width_all := apply(.SD, 1, function(x) sd_na(x) ), .SDcols=names(output_dt)[which(grepl("^Width_", names(output_dt)))]]

  #if cv_int is na (that means no replicates), then the median cv_int overall is assigned. 
  #output_dt[is.na(output_dt$feature_cv_intensity_all), ]$feature_cv_intensity_all <- median_na(output_dt$feature_cv_intensity_all)
  #wenguang: now it was changed to 75% of overall, to give some penalty...
  output_dt[is.na(output_dt$feature_cv_intensity_all), ]$feature_cv_intensity_all <- quantile(output_dt$feature_cv_intensity_all, breaks=c(0,0.25,0.5,0.75,1), na.rm=T)[4]
  
  #if sd_width is na (that means no replicates), then sd_width overall is assigned. 
  output_dt[is.na(output_dt$feature_sd_width_all), ]$feature_sd_width_all <- quantile(output_dt$feature_sd_width_all, breaks=c(0,0.25,0.5,0.75,1), na.rm=T)[4]
    
  if(global_level=="PeptideIon") {
    injections = paste0("Intensity_", anno$InjectionName)
  } else if(global_level=="Transition") {
    injections = paste0("aggr_Peak_Area_", anno$InjectionName)
  }
  
  result_cor <- input_dt[, {
    
    # for protein with multiple peptides
    if(dim(t(.SD[, injections, with=F]))[2] > 1) {
      
      mat_pcc <- cor(t(.SD[, injections, with=F]), use = "p", method="p")
      mat_scc <- cor(t(.SD[, injections, with=F]), use = "p", method="s")
      
      # To exclude the self pairwise comparison. 
      # Otherwise median_PCC of proteins with two peptides will be overestimated (because of at least there is a ONE between two compairons).
      diag(mat_pcc) <- NA
      diag(mat_scc) <- NA
      
      # at least three completed paires to calculate the cor
      mat_pairwise_complete <-  count_pairwise_number_matrix(t(.SD[, injections, with=F]))
      
      mat_pcc[mat_pairwise_complete < 3] <- NA
      mat_scc[mat_pairwise_complete < 3] <- NA
      
      
    } else {
      # this is for protein with only one peptide
      mat_pcc <- data.frame(x = -2)
      mat_scc <- data.frame(x = -2)
    }
    
    if(global_level=="PeptideIon") {
      .( PeptideIon = PeptideIon, 
         feature_median_PCC = apply(mat_pcc, 1, median_na),
         feature_median_SCC = apply(mat_scc, 1, median_na), 
         feature_MAD_dist = -2,
         median_PCC_PerProt = median_na(apply(mat_pcc, 1, function(x) median_na(x))) 
        )
    } else if(global_level=="Transition") {
      .( aggr_Fragment_Annotation = aggr_Fragment_Annotation,
         PeptideIon = PeptideIon, 
         feature_median_PCC = apply(mat_pcc, 1, median_na),
         feature_median_SCC = apply(mat_scc, 1, median_na), 
         feature_MAD_dist = -2,
         median_PCC_PerProt = median_na(apply(mat_pcc, 1, function(x) median_na(x))) 
        )
    }
    
  }, by=.(ProteinName)]
  
  if(global_level=="PeptideIon") {
    output_dt <- merge(output_dt, result_cor, by.x=c("PeptideIon", "ProteinName"), 
                                              by.y=c("PeptideIon", "ProteinName")  )
  } else if(global_level=="Transition") {
    output_dt <- merge(output_dt, result_cor, by.x=c("aggr_Fragment_Annotation", "PeptideIon", "ProteinName"), 
                                              by.y=c("aggr_Fragment_Annotation", "PeptideIon", "ProteinName")  )
  }
  
  output_dt[, `:=`(scaled_mean_intensity_all = scale(log2(feature_mean_intensity_all))
                   , scaled_cv_intensity_all = scale(feature_cv_intensity_all)
                   , scaled_numNA_intensity_all = scale(feature_numNA_intensity_all)
                   , scaled_averaged_score_all = scale(feature_averaged_score_all)
                   , scaled_sd_width_all = scale(feature_sd_width_all)
                   , scaled_median_PCC = scale(feature_median_PCC)
                   , scaled_MAD_dist = scale(feature_MAD_dist) )]
  
  
  #if median_PCC/median_SCC is na (that means many NAs), then the median median_PCC/median_SCC overall is assigned. 
  output_dt[is.na(output_dt$scaled_median_PCC), ]$scaled_median_PCC <- median_na(output_dt$scaled_median_PCC)
  
  # if all CVs were NAs (as all the samples have no replicates), scaled_cv_intensity would be set to zero.
  if( length(which(is.na(output_dt$scaled_cv_intensity_all))) == dim(output_dt)[1] ) {
    message("******WARNING******: It seems that all the CVs are NAs, probably due to no replicates available. ")
    message("******WARNING******: The feature of CV would not be used (i.e. all set to be zeros) in this analysis. \n")
    output_dt$scaled_cv_intensity_all <- 0
  }
  
  
  #the same as above. scaled_sd_width_all could be missing (i.e. they are all NAs), as the input matrix was from euler portal. 
  #In this case, this feature will not be used any more (i.e. all set to zero).
  if( length(which(is.na(output_dt$scaled_sd_width_all))) == dim(output_dt)[1] ) {
    message("******WARNING******: It seems that all the sd_width are NAs, probably because the input matrix was from euler portal. ")
    message("******WARNING******: The feature of sd_width would not be used (i.e. all set to be zeros) in this analysis. \n")
    output_dt$scaled_sd_width_all <- 0
  }
  
  
  message("done with calculating features...")
  
  return(output_dt)
  
}





#' @param input_dt 
#'
#' @param input_features 
#'
#' @export
get_lda_model <- function(input_dt, input_features) {
  
  model_lda <- lda(label ~., input_dt[, c(input_features), with=F])
  
  return(model_lda)
}


#' @param input_dt 
#'
#' @export
perform_selection <- function(input_dt) {
  
  #Wenugang: these two rules still need to be checked... as I think there is room to improve...
  output_dt_withCor <- input_dt[median_PCC_PerProt > 0.4 & numPerProt > 4, ]
  output_dt_withoutCor <- input_dt[-which(median_PCC_PerProt > 0.4 & numPerProt > 4), ] #wenguang: these include numPerProt<=4; median_PCC_PerProt==NA and median_PCC_PerProt <= 0.3 !
  
  
  # process data withCor first

  pred_lda <- predict(model_lda_withCor, output_dt_withCor[, c("scaled_mean_intensity_all", 
                                                               "scaled_cv_intensity_all", 
                                                               "scaled_numNA_intensity_all", 
                                                               "scaled_averaged_score_all",
                                                               "scaled_sd_width_all",
                                                               "scaled_median_PCC"), with=F])
  
  output_dt_withCor[, prob := 0]
  output_dt_withCor$prob <- pred_lda$posterior[,2]
  
  
  # then process data withoutCor
  pred_lda_1 <- predict(model_lda_withoutCor, output_dt_withoutCor[, c("scaled_mean_intensity_all", 
                                                                       "scaled_cv_intensity_all", 
                                                                       "scaled_numNA_intensity_all", 
                                                                       "scaled_averaged_score_all",
                                                                       "scaled_sd_width_all"), with=F])
  
  output_dt_withoutCor[, prob := 0]
  output_dt_withoutCor$prob <- pred_lda_1$posterior[,2]
  
  result <- rbind(output_dt_withCor, output_dt_withoutCor)
  
  return(result)

}