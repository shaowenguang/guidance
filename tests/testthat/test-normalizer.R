context("normalizer")

test_that("long2wide", {
  peptideIons <- import_openswath(search_results= "data/QGS_SWATH_data", 
                                  sample_annotation="data/QGS_sample_annotation", 
                                  level="PeptideIon")
  all_peptideIons <- long2wide(peptideIons)
                                  
  expect_equivalent(dim(all_peptideIons)[1], 52)
  expect_equivalent(dim(all_peptideIons)[2], 62)
  
  colnames <- colnames(all_peptideIons)
  
  expect_equivalent(colnames[1:2], c("PeptideIon", "ProteinName"))
  expect_equivalent(colnames[3], "Intensity_lgillet_J170408_001")
  expect_equivalent(colnames[18], "Score_lgillet_J170408_001")
  expect_equivalent(colnames[33], "RT_lgillet_J170408_001")
  expect_equivalent(colnames[48], "Width_lgillet_J170408_001" )
  
})

test_that("normalize_data", {
  data("all_peptideIons")
  all_peptideIons_normalized <- normalize_data(all_peptideIons, replaceNA="keep", 
                                               normalization="none")
  
  expect_equivalent(dim(all_peptideIons_normalized)[1], 52)
  expect_equivalent(dim(all_peptideIons_normalized)[2], 63)
  
  colnames <- colnames(all_peptideIons_normalized)
  
  expect_equivalent(colnames[1:2], c("PeptideIon", "ProteinName"))
  expect_equivalent(colnames[3], "Intensity_lgillet_J170408_001")
  expect_equivalent(colnames[18], "Score_lgillet_J170408_001")
  expect_equivalent(colnames[33], "RT_lgillet_J170408_001")
  expect_equivalent(colnames[48], "Width_lgillet_J170408_001" )
  expect_equivalent(colnames[63], "numPerProt" )
  
})


test_that("merge_replicates", {
  data("all_peptideIons")
  all_peptideIons_normalized <- normalize_data(all_peptideIons, replaceNA="keep", 
                                               normalization="none")
  cons_peptideIons <- merge_replicates(all_peptideIons_normalized, anno)
  
  expect_equivalent(dim(cons_peptideIons)[1], 52)
  expect_equivalent(dim(cons_peptideIons)[2], 78)
  
  colnames <- colnames(cons_peptideIons)
  
  expect_equivalent(colnames[1:2], c("PeptideIon", "ProteinName"))
  expect_equivalent(colnames[3], "Intensity_lgillet_J170408_001")
  expect_equivalent(colnames[18], "Score_lgillet_J170408_001")
  expect_equivalent(colnames[33], "RT_lgillet_J170408_001")
  expect_equivalent(colnames[48], "Width_lgillet_J170408_001" )
  expect_equivalent(colnames[63], "numPerProt" )
  expect_equivalent(colnames[64:66], c("mean_intensity_A", "cv_intensity_A", "numNA_intensity_A"))

})


test_that("merge_replicates", {
  data("all_peptideIons")
  proteotypic_peptide_matrix <- keep_proteotypic_only(all_peptideIons)
  
  
  expect_equivalent(dim(proteotypic_peptide_matrix)[1], 52)
  expect_equivalent(dim(proteotypic_peptide_matrix)[2], 62)
  
  expect_equivalent(length(which(grepl("1/sp", proteotypic_peptide_matrix$ProteinName))), 
                    dim(proteotypic_peptide_matrix)[1]) 
  
})




