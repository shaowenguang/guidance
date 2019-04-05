context("quantification")

peptideIons <- import_openswath(search_results= "data/QGS_SWATH_data", 
                                sample_annotation="data/QGS_sample_annotation", 
                                level="PeptideIon")
all_peptideIons <- long2wide(peptideIons)
peptideIon_n <- normalize_data(all_peptideIons, replaceNA="keep", normalization="none")
d <- merge_replicates(peptideIon_n, anno)
d_feature <- calc_features(d)
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


test_that("impute_missing_values", {
  imputated <- impute_missing_values(d_feature_select, c(3:17))
  
  expect_equivalent(round(mean(as.matrix(d_feature_select[,3:18]), na.rm=TRUE),2), 18879.38)
  expect_equivalent(round(mean(is.na(d_feature_select)), 5), 0.16559)
  
  expect_equivalent(round(mean(as.matrix(imputated[,3:18]), na.rm=TRUE),2), 15915.01)
  expect_equivalent(round(mean(is.na(imputated)), 5), 0.13441)
  
})


