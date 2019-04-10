context("selecter")


test_that("calc_features", {
  d_feature <- calc_features(peptideIon_st)
  
  expect_equivalent(colnames(d_feature)[1:2], c("PeptideIon", "ProteinName" ))
  expect_equivalent(colnames(d_feature)[79:94], 
                    c( "feature_mean_intensity_all",  "feature_cv_intensity_all", 
                       "feature_numNA_intensity_all", "feature_averaged_score_all",
                       "feature_sd_width_all",        "feature_median_PCC" , 
                       "feature_median_SCC",          "feature_MAD_dist" , 
                       "median_PCC_PerProt", "scaled_mean_intensity_all",
                       "scaled_cv_intensity_all", "scaled_numNA_intensity_all", 
                       "scaled_averaged_score_all", "scaled_sd_width_all",  
                       "scaled_median_PCC","scaled_MAD_dist"))
  
  expect_equivalent(summary(as.matrix(d_feature$feature_mean_intensity_all))[1], "Min.   : 2383  ")
  expect_equivalent(summary(as.matrix(d_feature$feature_mean_intensity_all))[2],  "1st Qu.: 6347  ")
  expect_equivalent(summary(as.matrix(d_feature$feature_mean_intensity_all))[4], "Mean   :16927  ")
  
})


test_that("perform_prediction_and_filtering", {
  d_feature <- calc_features(peptideIon_st)
  d_select <- perform_prediction_and_filtering(d_feature)

  expect_equivalent(colnames(d_select)[1:2], c("PeptideIon", "ProteinName" ))
  expect_equivalent(colnames(d_select)[79:95], 
                    c( "feature_mean_intensity_all",  "feature_cv_intensity_all", 
                       "feature_numNA_intensity_all", "feature_averaged_score_all",
                       "feature_sd_width_all",        "feature_median_PCC" , 
                       "feature_median_SCC",          "feature_MAD_dist" , 
                       "median_PCC_PerProt", "scaled_mean_intensity_all",
                       "scaled_cv_intensity_all", "scaled_numNA_intensity_all", 
                       "scaled_averaged_score_all", "scaled_sd_width_all",  
                       "scaled_median_PCC","scaled_MAD_dist", "prob"))
  
  expect_equivalent(summary(as.matrix(d_select$prob))[1], "Min.   :0.2026  ")
  expect_equivalent(summary(as.matrix(d_select$prob))[2],  "1st Qu.:0.4008  ")
  expect_equivalent(summary(as.matrix(d_select$prob))[4],  "Mean   :0.6021  ")
  
})


test_that("get_lda_model", {
  d_feature <- calc_features(peptideIon_st)
  
  ecoli_std <- c(2, 3, 4, 6, 8)
  index_mean_int <- which(grepl("^mean_intensity", names(d_feature)))
  d_feature[, cor_std := 0]
  d_feature$cor_std <- apply(d_feature[, index_mean_int, with=F], 1, function(x) cor(x, ecoli_std, use="p"))
  d_feature$cor_std[apply(d_feature[, index_mean_int, with=F], 1, function(x) count_pairwise_number(x, ecoli_std)) < 4] <- NA
  d_feature[, label := "bad"]
  d_feature[ cor_std > 0.95, ]$label <- "good"
  
  index_features <- c("scaled_mean_intensity_all", "scaled_cv_intensity_all", 
                      "scaled_numNA_intensity_all", "scaled_averaged_score_all", 
                      "scaled_median_PCC", "scaled_sd_width_all", "label")
  model_lda <- get_lda_model(d_feature, index_features)
  
  expect_equivalent(names(model_lda$prior), c( "bad","good"))
  expect_equivalent(round(model_lda$prior[1], 2), 0.56)
  expect_equivalent(round(model_lda$prior[2],2), 0.44)
  
  expect_equivalent(dimnames(model_lda$means)[[2]], 
                    c("scaled_mean_intensity_all", "scaled_cv_intensity_all",
                      "scaled_numNA_intensity_all", "scaled_averaged_score_all", 
                      "scaled_median_PCC", "scaled_sd_width_all")) 
  
  expect_equivalent(round(model_lda$means[1],2), -0.5)
  expect_equivalent(round(model_lda$scaling[2], 2), -0.08)
  expect_equivalent(round(model_lda$scaling[3], 2), -0.03)
 
})
