
# load library 
library(Prom)
library(data.table)
library(MASS)

library(ggplot2)
library(grid)
library(gplots)
library(plyr)
library(GGally)
library(ggfortify)
library(gridExtra)
library(qvalue) # not available for R-devel nor R 3.5.1 
                # Error: Bioconductor version '3.8' requires R version '3.5'; see https://bioconductor.org/install

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("qvalue", version = "3.8")


#wenguang: it will creat a global variable anno for sample_annotation table...
#peptideIons <- import_openswath(search_results="feature_alignment.csv", sample_annotation="sample_annotation", level="PeptideIon") 
peptideIons <- import_openswath(search_results="E1508100902_feature_alignment.tsv", sample_annotation="sample_annotation", level="PeptideIon")
peptideIons <- import_openswath(search_results= "D:/SWATH-guidance/E1508100902_feature_alignment.tsv", 
                                sample_annotation="D:/SWATH-guidance/sample_annotation", level="PeptideIon")


all_peptideIons <- long2wide(peptideIons)

#all_peptideIons_normalized <- normalize_data(all_peptideIons, replaceNA="keep", normalization="mediancenter")
all_peptideIons_normalized <- normalize_data(all_peptideIons, replaceNA="keep", normalization="none")
  # no normlization to adjust values, but added numPerProt column 

cons_peptideIons <- merge_replicates(all_peptideIons_normalized, anno)

cons_peptideIons <- cons_peptideIons[which(grepl("^1/", cons_peptideIons$ProteinName)), ]
cons_peptideIons_features <- calc_features(cons_peptideIons)

index_mean_int <- which(grepl("^mean_intensity", names(cons_peptideIons_features)))


test <- perform_selection(cons_peptideIons_features)

test_noFiltered <- copy(test)
test_yesFiltered <- test[prob > 0.2, ]

test_yesFiltered_yesImputated <- imputate_missing_values(test_yesFiltered, c(3:17))
 test_noFiltered_yesImputated <- imputate_missing_values( test_noFiltered, c(3:17))

test_noFiltered_noImputated <- copy(test_noFiltered)
test_yesFiltered_noImputated <- copy(test_yesFiltered)



cons_prot_test_yesFiltered_top1_sum_yesImputated_noWeighted <- merge_replicates(pept2prot(test_yesFiltered_yesImputated, "prob", 1, aggfun="sum", bool_weighted_by_prob=F), anno)
cons_prot_test_yesFiltered_top3_sum_yesImputated_noWeighted <- merge_replicates(pept2prot(test_yesFiltered_yesImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=F), anno)
cons_prot_test_yesFiltered_top99999_sum_yesImputated_noWeighted <- merge_replicates(pept2prot(test_yesFiltered_yesImputated, "prob", 99999, aggfun="sum", bool_weighted_by_prob=F), anno)


cons_prot_test_noFiltered_top1_sum_yesImputated_noWeighted <- merge_replicates(pept2prot(test_noFiltered_yesImputated, "prob", 1, aggfun="sum", bool_weighted_by_prob=F), anno)
cons_prot_test_noFiltered_top3_sum_yesImputated_noWeighted <- merge_replicates(pept2prot(test_noFiltered_yesImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=F), anno)
cons_prot_test_noFiltered_top99999_sum_yesImputated_noWeighted <- merge_replicates(pept2prot(test_noFiltered_yesImputated, "prob", 99999, aggfun="sum", bool_weighted_by_prob=F), anno)


cons_prot_test_yesFiltered_top1_sum_yesImputated_yesWeighted <- merge_replicates(pept2prot(test_yesFiltered_yesImputated, "prob", 1, aggfun="sum", bool_weighted_by_prob=T), anno)
cons_prot_test_yesFiltered_top3_sum_yesImputated_yesWeighted <- merge_replicates(pept2prot(test_yesFiltered_yesImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=T), anno)
cons_prot_test_yesFiltered_top99999_sum_yesImputated_yesWeighted <- merge_replicates(pept2prot(test_yesFiltered_yesImputated, "prob", 99999, aggfun="sum", bool_weighted_by_prob=T), anno)


cons_prot_test_noFiltered_top1_sum_yesImputated_yesWeighted <- merge_replicates(pept2prot(test_noFiltered_yesImputated, "prob", 1, aggfun="sum", bool_weighted_by_prob=T), anno)
cons_prot_test_noFiltered_top3_sum_yesImputated_yesWeighted <- merge_replicates(pept2prot(test_noFiltered_yesImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=T), anno)
cons_prot_test_noFiltered_top99999_sum_yesImputated_yesWeighted <- merge_replicates(pept2prot(test_noFiltered_yesImputated, "prob", 99999, aggfun="sum", bool_weighted_by_prob=T), anno)






cons_prot_test_yesFiltered_top1_sum_noImputated_noWeighted <- merge_replicates(pept2prot(test_yesFiltered_noImputated, "prob", 1, aggfun="sum", bool_weighted_by_prob=F), anno)
cons_prot_test_yesFiltered_top3_sum_noImputated_noWeighted <- merge_replicates(pept2prot(test_yesFiltered_noImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=F), anno)
cons_prot_test_yesFiltered_top99999_sum_noImputated_noWeighted <- merge_replicates(pept2prot(test_yesFiltered_noImputated, "prob", 99999, aggfun="sum", bool_weighted_by_prob=F), anno)


cons_prot_test_noFiltered_top1_sum_noImputated_noWeighted <- merge_replicates(pept2prot(test_noFiltered_noImputated, "prob", 1, aggfun="sum", bool_weighted_by_prob=F), anno)
cons_prot_test_noFiltered_top3_sum_noImputated_noWeighted <- merge_replicates(pept2prot(test_noFiltered_noImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=F), anno)
cons_prot_test_noFiltered_top99999_sum_noImputated_noWeighted <- merge_replicates(pept2prot(test_noFiltered_noImputated, "prob", 99999, aggfun="sum", bool_weighted_by_prob=F), anno)


cons_prot_test_yesFiltered_top1_sum_noImputated_yesWeighted <- merge_replicates(pept2prot(test_yesFiltered_noImputated, "prob", 1, aggfun="sum", bool_weighted_by_prob=T), anno)
cons_prot_test_yesFiltered_top3_sum_noImputated_yesWeighted <- merge_replicates(pept2prot(test_yesFiltered_noImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=T), anno)
cons_prot_test_yesFiltered_top99999_sum_noImputated_yesWeighted <- merge_replicates(pept2prot(test_yesFiltered_noImputated, "prob", 99999, aggfun="sum", bool_weighted_by_prob=T), anno)


cons_prot_test_noFiltered_top1_sum_noImputated_yesWeighted <- merge_replicates(pept2prot(test_noFiltered_noImputated, "prob", 1, aggfun="sum", bool_weighted_by_prob=T), anno)
cons_prot_test_noFiltered_top3_sum_noImputated_yesWeighted <- merge_replicates(pept2prot(test_noFiltered_noImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=T), anno)
cons_prot_test_noFiltered_top99999_sum_noImputated_yesWeighted <- merge_replicates(pept2prot(test_noFiltered_noImputated, "prob", 99999, aggfun="sum", bool_weighted_by_prob=T), anno)










prot <- fread("../mapDIA_v2.0/top99999/protein_level.txt")

prot[prot==0.0] <- NA #wenguang: Attention!!! This will make a huge difference, as protein table contains lots of zeros.

names(prot)[1] <- "ProteinName"

prot[, `:=`( mean_intensity_A = apply(.SD, 1, mean_na), cv_intensity_A = apply(.SD, 1, cv_na)), .SDcols=names(prot)[2:4]]
prot[, `:=`( mean_intensity_B = apply(.SD, 1, mean_na), cv_intensity_B = apply(.SD, 1, cv_na)), .SDcols=names(prot)[5:7]]
prot[, `:=`( mean_intensity_C = apply(.SD, 1, mean_na), cv_intensity_C = apply(.SD, 1, cv_na)), .SDcols=names(prot)[8:10]]
prot[, `:=`( mean_intensity_D = apply(.SD, 1, mean_na), cv_intensity_D = apply(.SD, 1, cv_na)), .SDcols=names(prot)[11:13]]
prot[, `:=`( mean_intensity_E = apply(.SD, 1, mean_na), cv_intensity_E = apply(.SD, 1, cv_na)), .SDcols=names(prot)[14:16]]


merged <- merge(cons_prot_test_filter_top99999_sum_imputated, prot, by.x="ProteinName", by.y="ProteinName")


mouse <- merged[which(grepl("MOUSE", merged$ProteinName)), ]
temp <- mouse[numPerProt > 1, ]

temp[, error.x := 0]
temp[, error.y := 0]

temp$error.x <- abs(log2(temp$mean_intensity_E.x) - log2(temp$mean_intensity_C.x)  + log2(3))
temp$error.y <- abs(log2(temp$mean_intensity_E.y) - log2(temp$mean_intensity_C.y)  + log2(3))



############################### train model ######################################
model_lda_ecoli <- get_lda_model(ecoli[numPerProt > 4 , ], index_feature_selected)

# data preparation 
ecoli_std <- c(2,3,4,6,8)

ecoli <- cons_peptideIons_features[grepl("ECOL", cons_peptideIons_features$ProteinName), ]

ecoli[, cor_std := 0]
ecoli$cor_std <- apply(ecoli[, index_mean_int, with=F], 1, function(x) cor(x, ecoli_std, use="p"))
ecoli$cor_std[apply(ecoli[, index_mean_int, with=F], 1, function(x) count_pairwise_number(x, ecoli_std)) < 4] <- NA


ecoli[, error := 0]
ecoli$error <- apply( cbind(abs(log2(ecoli$mean_intensity_A) - log2(ecoli$mean_intensity_B) - log2(2/3)), 
                            abs(log2(ecoli$mean_intensity_B) - log2(ecoli$mean_intensity_C) - log2(3/4)),
                            abs(log2(ecoli$mean_intensity_C) - log2(ecoli$mean_intensity_D) - log2(4/6)),
                            abs(log2(ecoli$mean_intensity_D) - log2(ecoli$mean_intensity_E) - log2(6/8)),
                            abs(log2(ecoli$mean_intensity_E) - log2(ecoli$mean_intensity_A) - log2(8/2)) ), 1, mean_na)

ecoli[, label := "bad"]
#ecoli[ cor_std > 0.95, ]$label <- "good"

ecoli[ error < 0.28, ]$label <- "good"

ecoli$label <- as.factor(ecoli$label)

index_feature_selected <- c("scaled_mean_intensity_all", "scaled_cv_intensity_all", "scaled_numNA_intensity_all", "scaled_averaged_score_all", "scaled_median_PCC", "scaled_sd_width_all", "label")

model_lda_ecoli <- get_lda_model(ecoli[numPerProt > 4 , ], index_feature_selected)


#library(devtools)
#model_lda_withCor <- copy(model_lda_ecoli)
#use_data(model_lda_withCor, pkg="~/project/rep/R_package/v2/Prom/R/", overwrite=T)