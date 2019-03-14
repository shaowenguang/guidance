#' Perform PECA test
#' 
#' @param input_file data table or data frame in wide representation. The data typically 
#' contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
#' protein measurements in rows. The input data is usually an output from 
#' \code{pept2prot()} and \code{merge_replicates()}
#' @param sample_annotation data matrix with \code{SampleName}, biological covariates 
#' (biological replicates) and technical covariates (technical replicates, batches, etc)
#' @param input_test the type of t-test to be conducted either ordinary 
#' \code{"t"} or modified \code{"modt"} t-test
#' @param input_bool_paired a logical indicating whether a paired test is performed
#'
#' @export
perform_peca_tests <- function(input_file, sample_annotation, input_test="modt",
                               input_bool_paired=FALSE) {
  
  numSamples <- length(unique(anno$SampleName)) 
  
  output_df <- data.frame()
  
  for(i in 1:numSamples) {
    
    for(j in 1:numSamples) {
      
      if(i > j) {
      
        message("Now performing ", input_test, " tests using PECA: ", unique(anno$SampleName)[i], " vs ", unique(anno$SampleName)[j])
        
        g1_peca <- paste0("Intensity_", anno[anno$SampleName==unique(anno$SampleName)[i], ]$Injection)
        g2_peca <- paste0("Intensity_", anno[anno$SampleName==unique(anno$SampleName)[j], ]$Injection)
      
        one_comp <- data.frame(PECA_tsv(file=input_file, samplenames1=g1_peca, samplenames2=g2_peca, test=input_test, paired=input_bool_paired))
      
        one_comp["ProteinName"] <- rownames(one_comp)
        one_comp["Label"] <- paste0(unique(anno$SampleName)[i], "/", unique(anno$SampleName)[j])
        
        #one_comp <- one_comp[c("ProteinName", "numPerProt", "Label", "log2fc", "t", "score", "pval", "padj")]
        one_comp <- one_comp[c("ProteinName", "n", "Label", "slr", "t", "score", "p", "p.fdr")]
        
        output_df <- rbind(output_df, one_comp)
      
       }
    
    }
    
  }
  
  return(output_df)
  
}



#' Perform t.test 
#' 
#' @param input_dt data table or data frame in wide representation. The data typically 
#' contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
#' protein measurements in rows. The input data is usually an output from 
#' \code{pept2prot()} and \code{merge_replicates()}
#'
#' @param sample_annotation data matrix with \code{SampleName}, biological covariates 
#' (biological replicates) and technical covariates (technical replicates, batches, etc)
#' @param input_bool_paired a logical indicating whether a paired test is performed
#' @param input_mtc_method multiple testing correction method. Options include
#' \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, 
#' \code{"BH"}, \code{"BY"}, \code{"fdr"} and \code{"none"}. Refer to 
#' \code{p.adjust()} for details. 
#'
#' @export
perform_t_tests <- function(input_dt, sample_annotation, input_bool_paired=FALSE, 
                            input_mtc_method="bonferroni") {
  
  numSamples <- length(unique(anno$SampleName)) 
  
  output_df <- data.frame()
  
  for(i in 1:numSamples) {
    
    for(j in 1:numSamples) {
      
      if(i > j) {
        
        message("Now performing normal t tests using genefilter package: ", unique(anno$SampleName)[i], " vs ", unique(anno$SampleName)[j])
        
        g1 <- paste0("Intensity_", anno[anno$SampleName==unique(anno$SampleName)[i], ]$Injection)
        g2 <- paste0("Intensity_", anno[anno$SampleName==unique(anno$SampleName)[j], ]$Injection)
        
        one_comp <- rowttests( as.matrix(cbind(log2(input_dt[, g1, with=F]+1), log2(input_dt[, g2, with=F]+1))), 
                   fac=factor(c(rep(unique(anno$SampleName)[i], length(g1)),
                                rep(unique(anno$SampleName)[j], length(g2)))) )
        
        one_comp["ProteinName"] <- input_dt$ProteinName
        one_comp["Label"] <- paste0(unique(anno$SampleName)[i], "/", unique(anno$SampleName)[j])
        one_comp["numPerProt"] <- input_dt$numPerProt
        one_comp["protQuantProb"] <- input_dt$protQuantProb
        
        one_comp["pval_adj"] <- p.adjust(one_comp$p.value, method=input_mtc_method)
        
        one_comp <- one_comp[c("ProteinName", "protQuantProb", "numPerProt", "Label", "dm", "statistic", "p.value", "pval_adj")]
        
        output_df <- rbind(output_df, one_comp)
      }
      
    }
  
  }
  
  return(output_df)
  
}


#' Performed modified t.test
#' 
#' @param input_dt data table or data frame in wide representation. The data typically 
#' contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
#' protein measurements in rows. The input data is usually an output from 
#' \code{pept2prot()} and \code{merge_replicates()}
#'
#' @param sample_annotation data matrix with \code{SampleName}, biological covariates 
#' (biological replicates) and technical covariates (technical replicates, batches, etc)
#' @param input_bool_paired a logical indicating whether a paired test is performed
#' @param input_mtc_method multiple testing correction method. Options include
#' \code{"holm"}, \code{"hochberg"}, \code{"hommel"}, \code{"bonferroni"}, 
#' \code{"BH"}, \code{"BY"}, \code{"fdr"} and \code{"none"}. Refer to 
#' \code{p.adjust()} for details. 
#'
#' @export
perform_modt_tests <- function(input_dt, sample_annotation, input_bool_paired=FALSE, input_mtc_method="BH") {
  
  numSamples <- length(unique(anno$SampleName)) 
  
  output_df <- data.frame()
  
  for(i in 1:numSamples) {
    
    for(j in 1:numSamples) {
      
      if(i > j) {
        
        message("Now performing modt tests using limma package: ", 
                unique(anno$SampleName)[i], " vs ", unique(anno$SampleName)[j])
        
        g1 <- paste0("Intensity_", anno[anno$SampleName==unique(anno$SampleName)[i], ]$Injection)
        g2 <- paste0("Intensity_", anno[anno$SampleName==unique(anno$SampleName)[j], ]$Injection)
        
        design <- cbind( G1=1, 
                         G1vsG2=c( rep(1,length(g1)), rep(0,length(g2))) )
        
        #one_comp <- rowttests( as.matrix(cbind(log2(input_dt[, g1, with=F]+1), log2(input_dt[, g2, with=F]+1))), 
        #                       fac=factor(c(rep(unique(anno$SampleName)[i], length(g1)), rep(unique(anno$SampleName)[j], length(g2)))) )
        
        one_comp <- data.frame(ProteinName = input_dt$ProteinName,
                               Label = paste0(unique(anno$SampleName)[i], 
                                              "/", unique(anno$SampleName)[j]),
                               numPerProt = input_dt$numPerProt,
                               protQuantProb = input_dt$protQuantProb )
        
        limma_input <- as.matrix(cbind(log2(input_dt[, g1, with=F]+1), 
                                       log2(input_dt[, g2, with=F]+1)))
        fit <- lmFit(limma_input, design)
        efit <- eBayes(fit)
        
        one_comp["dm"] <- efit$coefficients[,2]
        one_comp["statistic"] <- efit$t[,2]
        one_comp["p.value"] <- efit$p.value[, 2]
        one_comp["pval_adj"] <- p.adjust(one_comp$p.value, method=input_mtc_method)
        
        
        #one_comp["ProteinName"] <- input_dt$ProteinName
        #one_comp["Label"] <- paste0(unique(anno$SampleName)[i], "/", unique(anno$SampleName)[j])
        #one_comp["numPerProt"] <- input_dt$numPerProt
        #one_comp["protQuantProb"] <- input_dt$protQuantProb
        

        
        one_comp <- one_comp[c("ProteinName", "protQuantProb", "numPerProt", "Label", "dm", "statistic", "p.value", "pval_adj")]
        
        output_df <- rbind(output_df, one_comp)
        
      }
      
    }
    
  }
  
  return(output_df)
  
}


#' Perform ANOVA
#'
#' @param input_dt @param input_dt data table or data frame in wide representation. The data typically 
#' contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
#' protein measurements in rows. The input data is usually an output from 
#' \code{pept2prot()} and \code{merge_replicates()}
#'
#' @param sample_annotation data matrix with \code{SampleName}, biological covariates 
#' (biological replicates) and technical covariates (technical replicates, batches, etc)
#'
#' @return data.table data.frame containing statistics of parametric and 
#' non-parametric anova results
#' 
#' @export 
#'
perform_anova <- function(input_dt, sample_annotation) {
  
  proteinName <- input_dt$ProteinName
  sampleName <- sample_annotation$SampleName
  
  cols <- which(grepl("^Intensity_", names(input_dt)))
  data_t <- cbind(data.frame(sampleName), data.frame(t(data.frame(input_dt)[,cols])))
  
  baseformula <- " ~ sampleName"
  output <- matrix(nrow = ncol(data_t)-1, ncol = 4, 
                   dimnames = list(proteinName, c("parametric_Fvalue", "parametric_pvalue", 
                                                  "KruskalWallis_chiSquared", "KruskalWallis_pvalue")))
  
  message("start to compute anova...")
  for (i in 2:ncol(data_t)){
    
    formula <- paste(colnames(data_t)[i], baseformula, sep="")
    
    protein_row <- i - 1
    if( (protein_row %% 1000) == 0 ) {
      message("computing anova for ",  (protein_row %/% 1000),
              "000 out of ", dim(input_dt)[1], " peptide ions...")
    }
    
    aov_summary <- summary(aov(as.formula(formula), data=data_t))[[1]]
    output[protein_row,"parametric_Fvalue"] <- signif(aov_summary[["F value"]][1], 4)
    output[protein_row,"parametric_pvalue"] <- signif(aov_summary[["Pr(>F)"]][1], 4)
    
    kruskal_summary <- kruskal.test(as.formula(formula), data=data_t)
    output[protein_row,"KruskalWallis_chiSquared"] <- signif(kruskal_summary$statistic, 4)
    output[protein_row,"KruskalWallis_pvalue"] <- signif(kruskal_summary$p.value, 4)
  }
  
  output_df <- cbind(input_dt, data.frame(output))
  return(output_df)
  
}


