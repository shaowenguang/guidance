#' @param input_file 
#'
#' @param anno 
#' @param input_test 
#' @param input_bool_paired 
#'
#' @export
perform_peca_tests <- function(input_file, anno, input_test="modt", input_bool_paired=FALSE) {
  
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



#' @param input_dt 
#'
#' @param anno 
#' @param input_bool_paired 
#' @param input_mtc_method 
#'
#' @export
perform_t_tests <- function(input_dt, anno, input_bool_paired=FALSE, input_mtc_method="bonferroni") {
  
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


#' @param input_dt 
#'
#' @param anno 
#' @param input_bool_paired 
#' @param input_mtc_method 
#'
#' @export
perform_modt_tests <- function(input_dt, anno, input_bool_paired=FALSE, input_mtc_method="BH") {
  
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


perform_anova <- function(
  
  
)



