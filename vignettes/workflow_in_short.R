###importing file from guidance/data
peptideIons <- import_spectronaut_matrix(search_results= "data/QGS_SWATH_data", 
                                          sample_annotation="data/QGS_sample_annotation"
                                         )

peptideIons <- import_openswath(search_results= "data/QGS_SWATH_data", 
                                sample_annotation="data/QGS_sample_annotation", 
                                level="PeptideIon")

all_peptideIons <- long2wide(peptideIons)


# import data
peptideIons <- import_openswath(search_results= "S:/SWATH-guidance/feature_alignment.csv", 
                                sample_annotation="S:/SWATH-guidance/sample_annotation", 
                                level="PeptideIon") 

# normalize data (if necessary)
all_peptideIons <- long2wide(peptideIons)
all_peptideIons_normalized <- normalize_data(all_peptideIons, replaceNA="keep", 
                                             normalization="none")

# merge replicates
cons_peptideIons <- merge_replicates(all_peptideIons_normalized, anno)

# filter proteotypic peptides 
cons_peptideIons <- cons_peptideIons[which(grepl("^1/", cons_peptideIons$ProteinName)), ]

# calculate features 
cons_peptideIons_features <- calc_features(cons_peptideIons)

# calculate posterior probability of being RARE using LDA model 
index_mean_int <- which(grepl("^mean_intensity", names(cons_peptideIons_features)))
test <- perform_selection(cons_peptideIons_features)

# filter by probability, impute missing value 
test_yesFiltered <- test[prob > 0.2, ]
test_yesFiltered_yesImputated <- impute_missing_values(test_yesFiltered, c(3:17))

# peptide to protein inference 
cons_prot_test_yesFiltered_top3_sum_yesImputated_yesWeighted <- merge_replicates(
  pept2prot(test_yesFiltered_yesImputated, "prob", 3, aggfun="sum", bool_weighted_by_prob=T), anno)

# statistical tests 
t_test_result <- perform_t_tests(cons_prot_test_yesFiltered_top3_sum_yesImputated_yesWeighted)
modt_result <- perform_modt_tests(cons_prot_test_yesFiltered_top3_sum_yesImputated_yesWeighted)
peca_result <- perform_peca_tests(cons_prot_test_yesFiltered_top3_sum_yesImputated_yesWeighted)

# visualization 
prot_name <- "1/sp|P37108|SRP14_HUMAN" 
test_prot <- test[test$ProteinName==prot_name, ]

pdf("example_protein_profiles.pdf", width=7.5*3, height=4.1*2)
plot_a_heatmap_include_prob_update(test_prot)
dev.off()

plot_bar_intensity_prob(test_prot)
plot_heatmap(test_prot)
plot_peptide_trend(test_prot)