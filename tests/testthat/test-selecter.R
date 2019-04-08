context("selecter")

peptideIons <- import_openswath(search_results= "data/QGS_SWATH_data", 
                                sample_annotation="data/QGS_sample_annotation", 
                                level="PeptideIon")
all_peptideIons <- long2wide(peptideIons)
peptideIon_n <- normalize_data(all_peptideIons, replaceNA="keep", normalization="none")
d <- merge_replicates(peptideIon_n, anno)
d_feature <- calc_features(d)
d_feature_select <- perform_selection(d_feature)


test_that("calc_features", {
  peptideIons_features <- calc_features(all_peptideIons)
  
  colnames <- colnames(peptideIons_features)
  expect_equivalent(colnames[1:2], c("PeptideIon", "ProteinName" ))
  expect_equivalent(colnames[63:65],  c("feature_mean_intensity_all",
                    "feature_cv_intensity_all", "feature_numNA_intensity_all"))
  
  expect_equivalent(round(mean(as.matrix(peptideIons_features$), na.rm=TRUE),2), 33409.88)
  
})

peptideIons_features <- calc_features(all_peptideIons)
peptideIons_features_select <- perform_selection(peptideIons_features)

