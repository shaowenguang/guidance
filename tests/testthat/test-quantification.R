context("quantification")

### need to change the data!! calc_features require global_level == "PeptideIon" which
### does not get saved if only the processed data is saved. 
d_feature <- calc_features(peptideIon_st)
d_feature_select <- perform_selection(d_feature)

test_that("pept2prot", {
  protein_m <- pept2prot(d_feature_select, "prob", 3, aggfun="sum", bool_weighted_by_prob=T)
  
  expect_equivalent(dim(protein_m)[1], 9)
  expect_equivalent(dim(protein_m)[2], 18)
  
  colnames <- colnames(protein_m)
  expect_equivalent(colnames[1:2], c("ProteinName", "numPerProt"))
  expect_equivalent(colnames[4],  "Intensity_lgillet_J170408_001")
  
  expect_equivalent(round(mean(as.matrix(protein_m[,3:18]), na.rm=TRUE),2), 33409.88)

})



test_that("pept2prot", {
  protein_log <- pept2prot_log2(d_feature_select, 
                                input_rank_index = "prob", 3, aggfun="sum", bool_weighted_by_prob=T)  
  expect_equivalent(dim(protein_log)[1], 9)
  expect_equivalent(dim(protein_log)[2], 18)
  
  colnames <- colnames(protein_log)
  expect_equivalent(colnames[1:2], c("ProteinName", "numPerProt"))
  expect_equivalent(colnames[4],  "Intensity_lgillet_J170408_001")
  
  expect_equivalent(round(mean(as.matrix(protein_log[,3:18]), na.rm=TRUE),2), 18.35)
  
})



test_that("impute_missing_values", {
  imputated <- impute_missing_values(d_feature_select, c(3:17))
  
  expect_equivalent(round(mean(as.matrix(d_feature_select[,3:18]), na.rm=TRUE),2), 18879.38)
  expect_equivalent(round(mean(is.na(d_feature_select)), 5), 0.16559)
  
  expect_equivalent(round(mean(as.matrix(imputated[,3:18]), na.rm=TRUE),2), 15915.01)
  expect_equivalent(round(mean(is.na(imputated)), 5), 0.13441)
  
})


