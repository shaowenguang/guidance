context("Importer")

test_that("import_openswath", {
  peptideIons <- import_openswath(search_results= "data/QGS_SWATH_data", 
                                  sample_annotation="data/QGS_sample_annotation", 
                                  level="PeptideIon")
  
  expect_equivalent(dim(peptideIons)[1], 612)
  expect_equivalent(dim(peptideIons)[2], 8)
  expect_equivalent(colnames(peptideIons)[1], "PeptideIon")
  
  expect_equivalent(dim(anno)[1], 15)
  expect_equivalent(dim(anno)[2], 8)
  
})


test_that("import_openswath_matrix_fromEulerPortal", {
  all_peptideIons <- import_openswath_matrix_fromEulerPortal(
    search_results="data/YS_SWATH_data", 
    sample_annotation="data/YS_sample_annotation")
  
  expect_equivalent(dim(all_peptideIons)[1], 234)
  expect_equivalent(dim(all_peptideIons)[2], 65)
  expect_equivalent(colnames(all_peptideIons)[1],  "PeptideIons")
  expect_equivalent(colnames(all_peptideIons)[3:5], 
                    c("Intensity_uT01","RT_uT01","Score_uT01"))

  expect_equivalent(dim(anno)[1], 20)
  expect_equivalent(dim(anno)[2], 8)
  
})




test_that("import_spectronaut_matrix", {
  peptideIons <- import_spectronaut_matrix(search_results= "data/QGS_SWATH_data", 
                                           sample_annotation="data/QGS_sample_annotation")
  
  expect_equivalent(dim(peptideIons)[1], 612)
  expect_equivalent(dim(peptideIons)[2], 58)
  expect_equivalent(colnames(peptideIons)[1], "PeptideIon")
  expect_equivalent(colnames(peptideIons)[10:15], 
                    c("m/z", "Intensity_Intensity",
                      "ProteinName", "decoy", "assay_rt", "delta_rt"))
  
  expect_equivalent(dim(anno)[1], 15)
  expect_equivalent(dim(anno)[2], 8)
  
})

