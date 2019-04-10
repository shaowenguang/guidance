context("wrapper")

roundup_to <- function(x, to = 10, up = FALSE){
  if(up) round(.Machine$double.eps^0.5 + x/to)*to else round(x/to)*to
}

test_that("dia_guidance", {
  guidance <- dia_guidance(data= "data/QGS_SWATH_data", data_type = "openswath",
                           sample_annotation="data/QGS_sample_annotation", level="PeptideIon",
                           replaceNA="keep", bool_NA_means_requant = FALSE, 
                           averageFun = "mean", normalization="none", filter_prob = 0.25, 
                           input_rank_index = "prob", topN = 3, aggfun = "sum", 
                           bool_weighted_by_prob = TRUE, bool.removeDecoy = T, remove_prefixInFileName = FALSE)
  
  expect_equivalent(colnames(guidance)[1:2], c( "ProteinName", "numPerProt" ))
  expect_equivalent(colnames(guidance)[4], "Intensity_lgillet_J170408_001")
  expect_equivalent(  colnames(guidance)[31:33], c("mean_intensity_E","cv_intensity_E", "numNA_intensity_E" ))

  expect_equivalent(roundup_to(mean(as.matrix(guidance$mean_intensity_D)), to = 100), 44500)
  expect_equivalent(roundup_to(mean(as.matrix(guidance$cv_intensity_D)), to = 0.1), 0.1)
  expect_equivalent(round(mean(as.matrix(guidance$numNA_intensity_C))), 0)

})

