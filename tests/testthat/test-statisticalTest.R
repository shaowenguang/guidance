context("statisticalTest")

library(PECA)
library(genefilter)
library(limma)

library(data.table)
library(MASS)
library(ggplot2)
library(grid)
library(gplots)
library(plyr)
library(GGally)
library(ggfortify)
library(gridExtra)
library(Prom)
require(plotrix)


test_that("perform_peca_tests", {
  
  # peca_result <- perform_peca_tests(protein_Filtered_top3_sum_ImputedWeighted, anno)
  # this function is currently buggy fl
  
  
})



test_that("perform_t_tests", {
  t_test_result <- perform_t_tests(protein_Filtered_top3_sum_ImputedWeighted, anno)
  
  colnames(t_test_result)[1] <- "ProteinName" 
  colnames(t_test_result)[7:8] <- c("p.value", "pval_adj" )

  expect_equivalent(round(mean(t_test_result$statistic, na.rm = TRUE), 2), -4.39)
  expect_equivalent(round(mean(t_test_result$p.value, na.rm = TRUE), 2), 0.23)
  expect_equivalent(round(mean(t_test_result$pval_adj, na.rm = TRUE), 2), 0.53)
  
})



test_that("perform_modt_tests", {
  modt_test_result <- perform_modt_tests(protein_Filtered_top3_sum_ImputedWeighted, anno)
  
  colnames(modt_test_result)[1] <- "ProteinName" 
  colnames(modt_test_result)[7:8] <- c("p.value", "pval_adj" )
  
  expect_equivalent(round(mean(modt_test_result$statistic, na.rm = TRUE), 2), 3.14)
  expect_equivalent(round(mean(modt_test_result$p.value, na.rm = TRUE), 2), 0.2)
  expect_equivalent(round(mean(modt_test_result$pval_adj, na.rm = TRUE), 2), 0.25)
  
})


test_that("perform_anova", {
  anova_results <- perform_anova(protein_Filtered_top3_sum_ImputedWeighted, anno)
  
  colnames(anova_results)[1] <- "ProteinName" 
  colnames(anova_results)[34:37] <- c("parametric_Fvalue", "parametric_pvalue",
                                      "KruskalWallis_chiSquared", "KruskalWallis_pvalue" )
  
  expect_equivalent(round(mean(anova_results$parametric_Fvalue, na.rm = TRUE), 2), 34.94)
  expect_equivalent(round(mean(anova_results$parametric_pvalue, na.rm = TRUE), 2), 0.18)
  expect_equivalent(round(mean(anova_results$KruskalWallis_chiSquared, na.rm = TRUE), 2), 9.05)
  expect_equivalent(round(mean(anova_results$KruskalWallis_pvalue, na.rm = TRUE), 2), 0.19)
  
})

