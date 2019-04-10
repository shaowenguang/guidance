context("visualizer")

library(gridExtra)
library(grid)

peptideIons <- import_openswath_matrix_fromEulerPortal(search_results= "data/YS_SWATH_data", 
                                sample_annotation="data/YS_sample_annotation")
peptideIon_n <- normalize_data(peptideIons, replaceNA="keep", normalization="none")
d <- merge_replicates(peptideIon_n, anno)
d_feature <- calc_features(d)
d_feature_select <- perform_selection(d_feature)
prot_name <- "1/sp|P45578|LUXS_ECOLI"
test_prot <- test[test$ProteinName==prot_name, ]



test_that("plot_protein_profile", {
  p <- plot_protein_profile(test_prot, cutoff_prob = 0.3)
  
  expect_equivalent(p$layout$t, c(1,2)) 
  expect_equivalent(p$layout$l, c(1, 1)) 
  
  expect_equivalent(class(p$grobs[[2]]), c("gtable", "gTree",  "grob",   "gDesc" ))
  expect_equivalent(p$grobs[[2]]$layout$z, c(1, 2, 3))
})


test_that("plot_peptide_intensity", {
  p <- plot_peptide_intensity(test_prot)
  
  expect_equivalent(summary(p$plot_env$a_reorder)[1],  "Min.   :0.2800  ")
  expect_equivalent(summary(p$plot_env$a_reorder)[2],   "1st Qu.:0.7650  ")
  expect_equivalent(summary(p$plot_env$a_reorder)[4],  "Mean   :0.7255  ")
  
  expect_equivalent( p$plot_env$bool_isPeptideIon, TRUE)
  expect_equivalent( p$plot_env$cutoff_prob, 0.3)
  expect_equivalent(p$plot_env$tmp_injections[1], "lgillet_J170408_001")
  
  expect_equivalent(p$labels$x, "Injections")
  expect_equivalent( p$labels$y, "log2(Intensity)")
  expect_equivalent(p$labels$colour, "PeptideIon")
  expect_equivalent( p$labels$group,  "PeptideIon")
  
})


test_that("plot_peptide_intensity", {
  p <- plot_cor_heatmap(test_prot)

  expect_equivalent(summary(p$plot_env$a_reorder)[1],  "Min.   :0.2800  ")
  expect_equivalent(summary(p$plot_env$a_reorder)[2],   "1st Qu.:0.7650  ")
  expect_equivalent(summary(p$plot_env$a_reorder)[4],  "Mean   :0.7255  ")
  
  expect_equivalent(p$layers[[2]]$aes_params$colour, "white")
  expect_equivalent(p$layers[[2]]$aes_params$size, 2.5)
  
  expect_equivalent(p$labels$x, "peptideIon")
  expect_equivalent(p$labels$fill, "Correlation")
  expect_equivalent(p$labels$label, "Correlation")
  
})


test_that("plot_peptide_intensity", {
  p <- plot_bar_intensity_n_probability(test_prot)
  
  expect_equivalent(p$layout$t, c(1,2))
  expect_equivalent(p$layout$l, c(1,1))
  expect_equivalent(p$name, "arrange")
  
})


test_that("plot_density", {
  p <- plot_density(d_feature_select, feature = "feature_mean_intensity_all")
  
  expect_equivalent(p$labels$x, "feature_mean_intensity_all")
  expect_equivalent(p$labels$y, "density")
  expect_equivalent(p$plot_env$feature,"feature_mean_intensity_all")
  
  expect_equivalent(p$layers[[1]]$aes_params$alpha, 0.8)
  
})




