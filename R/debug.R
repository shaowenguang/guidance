#' @export
## not finised...


#' @importFrom data.table melt



calcRunTime <- function(input_function, input_data) {

  FUN <- match.fun(input_function) 
  
  a <- Sys.time()  
  FUN(input_data)
  b <- Sys.time()
  
  cat("it takes", b-a, "\n")
  
}




#' @export
reorder_cormat <- function(cormat){ # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

#' @export
plot_a_heatmap <- function(case, index_intensity, index_cv=NULL, bool_reorder) {

  mean_int <- apply(case[, index_intensity, with=F], 1, function(x) mean(x, na.rm = TRUE))
  if(is.null(index_cv)) {
    mean_cv <- 0
  } else {
    mean_cv <- apply(case[, index_cv, with=F], 1, function(x) mean(x, na.rm = TRUE))
  }
  
  case[, mean_int := mean_int]
  case[, mean_cv := mean_cv]
    
  a <- round(cor(t(case[, index_intensity, with=F]), use="p", method="p"), 2)
  rownames(a) <- case$PeptideIon
  colnames(a) <- case$PeptideIon

  a[is.na(a)] <- 0
  
  if(bool_reorder==1) {
    a_reorder <- reorder_cormat(a)
    p_heatmap <- ggplot(melt(a_reorder), aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient2(limits=c(-1, 1), low="red", mid="black", high="green", space = "Lab", midpoint=0) + geom_text(aes(Var1, Var2, label=value), col="white", size=2.5)
    
    case$PeptideIon <- factor(case$PeptideIon, levels=case[colnames(a_reorder), on="PeptideIon"]$PeptideIon)
    #p_bar_int <- ggplot(case, aes(x=PeptideIon, y=mean_int, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip()
    
    limits <- aes(ymax = case$mean_int + case$mean_cv * case$mean_int, ymin = case$mean_int - case$mean_cv*case$mean_int) 
    p_bar_error <- ggplot(case, aes(x=PeptideIon, y=mean_int, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5)
   
    
    quant_mean_all <- apply(case[, index_intensity, with=F], 2, function(x) mean(x, na.rm = TRUE))
    
    #wenguang: order by mean_all, or without any orders
    case_long <- melt(case, id.vars="PeptideIon", measure.vars=names(case)[index_intensity], variable.name="index", value.name="int")
    #case_long <- melt(case, id.vars="PeptideIon", measure.vars=names(case)[index_intensity][order(quant_mean_all, decreasing = T)], variable.name="index", value.name="int")
    
    p_matplot <- ggplot(case_long, aes(x=index, y=log2(int))) + geom_line(aes(colour=PeptideIon, group=PeptideIon)) + theme(axis.text.x = element_text(angle = 90))
    
    #p_plot <- arrangeGrob(p_heatmap, p_bar_int, ncol=2)
    p_plot <- grid.arrange( p_matplot, arrangeGrob(p_heatmap, p_bar_error, ncol=2), nrow=2 )

    
  } else {
    p_plot <- ggplot(melt(a), aes(x=Var1, y=Var2, fill=value)) + geom_tile() 
                                                               + scale_fill_gradient2(limits=c(-1, 1), low="red", mid="black", high="green", space = "Lab", midpoint=0) 
                                                               + geom_text(aes(Var1, Var2, label=value), col="white", size=2.5)
    
  }
  
  
  return(p_plot)
  
}



plot_a_heatmap_include_prob_old <- function(case, index_intensity, index_cv=NULL, bool_reorder=TRUE) {
  
  if("aggr_Fragment_Annotation" %in% names(case)) {
    # the input is transition level data
    bool_isPeptideIon <- 0 
  } else {
    # the input is peptide level data
    bool_isPeptideIon <- 1
  }
  
  mean_int <- apply(case[, index_intensity, with=F], 1, function(x) mean(x, na.rm = TRUE))
  
  if(is.null(index_cv)) {
    mean_cv <- 0
  } else {
    mean_cv <- apply(case[, index_cv, with=F], 1, function(x) mean(x, na.rm = TRUE))
  }
  
  case[, mean_int := mean_int]
  case[, mean_cv := mean_cv]
  
  a <- round(cor(t(case[, index_intensity, with=F]), use="p", method="p"), 2)
  
  if(bool_isPeptideIon == TRUE) {
    rownames(a) <- case$PeptideIon
    colnames(a) <- case$PeptideIon
  } else {
    rownames(a) <- case$aggr_Fragment_Annotation
    colnames(a) <- case$aggr_Fragment_Annotation
  }
  
  a[is.na(a)] <- 0
  
  if(bool_reorder==1) {
    
    a_reorder <- reorder_cormat(a)
    p_heatmap <- ggplot(melt(a_reorder), aes(x=Var1, y=Var2, fill=value)) + geom_tile() + scale_fill_gradient2(limits=c(-1, 1), low="red", mid="black", high="green", space = "Lab", midpoint=0) + geom_text(aes(Var1, Var2, label=value), col="white", size=2.5)
    
    if(bool_isPeptideIon == TRUE) {
      
      case$PeptideIon <- factor(case$PeptideIon, levels=case[colnames(a_reorder), on="PeptideIon"]$PeptideIon)
      #p_bar_int <- ggplot(case, aes(x=PeptideIon, y=mean_int, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip()
      
      limits <- aes(ymax = case$mean_int + case$mean_cv * case$mean_int, ymin = case$mean_int - case$mean_cv*case$mean_int) 
      p_bar_error <- ggplot(case, aes(x=PeptideIon, y=mean_int, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5)
      p_bar_prob <- ggplot(case, aes(x=PeptideIon, y=prob, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip()
      
      
      quant_mean_all <- apply(case[, index_intensity, with=F], 2, function(x) mean(x, na.rm = TRUE))
      
      #wenguang: order by mean_all, or without any orders
      case_long <- melt(case, id.vars="PeptideIon", measure.vars=names(case)[index_intensity], variable.name="index", value.name="int")
      #case_long <- melt(case, id.vars="PeptideIon", measure.vars=names(case)[index_intensity][order(quant_mean_all, decreasing = T)], variable.name="index", value.name="int")
      
      p_matplot <- ggplot(case_long, aes(x=index, y=log2(int))) + geom_line(aes(colour=PeptideIon, group=PeptideIon)) + theme(axis.text.x = element_text(angle = 90))
      
    } else {
      
      case$aggr_Fragment_Annotation <- factor(case$aggr_Fragment_Annotation, levels=case[colnames(a_reorder), on="aggr_Fragment_Annotation"]$aggr_Fragment_Annotation)
      #p_bar_int <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=mean_int, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip()
      
      limits <- aes(ymax = case$mean_int + case$mean_cv * case$mean_int, ymin = case$mean_int - case$mean_cv*case$mean_int) 
      p_bar_error <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=mean_int, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5)
      p_bar_prob <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=prob, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip()
      
      
      quant_mean_all <- apply(case[, index_intensity, with=F], 2, function(x) mean(x, na.rm = TRUE))
      
      #wenguang: order by mean_all, or without any orders
      case_long <- melt(case, id.vars="aggr_Fragment_Annotation", measure.vars=names(case)[index_intensity], variable.name="index", value.name="int")
      #case_long <- melt(case, id.vars="aggr_Fragment_Annotation", measure.vars=names(case)[index_intensity][order(quant_mean_all, decreasing = T)], variable.name="index", value.name="int")
      
      p_matplot <- ggplot(case_long, aes(x=index, y=log2(int))) + geom_line(aes(colour=aggr_Fragment_Annotation, group=aggr_Fragment_Annotation)) + theme(axis.text.x = element_text(angle = 90))
      
      
    }
    
    #p_plot <- arrangeGrob(p_heatmap, p_bar_int, ncol=2)
    p_plot <- grid.arrange( p_matplot, arrangeGrob(p_heatmap, p_bar_error, p_bar_prob, ncol=3), nrow=2 )
    
    
  } else {
    
    p_plot <- ggplot(melt(a), aes(x=Var1, y=Var2, fill=value)) + geom_tile() 
    + scale_fill_gradient2(limits=c(-1, 1), low="red", mid="black", high="green", space = "Lab", midpoint=0) 
    + geom_text(aes(Var1, Var2, label=value), col="white", size=2.5)
    
  }
  
  
  return(p_plot)
  
}



#' @export
plot_a_heatmap_include_prob_update <- function(case, cutoff_prob=0.3) {
  
  if("aggr_Fragment_Annotation" %in% names(case)) {
    # the input is transition level data
    bool_isPeptideIon <- 0 
    index_intensity <- which(grepl("^aggr_Peak_Area_", names(case)))
  } else {
    # the input is peptide level data
    bool_isPeptideIon <- 1
    index_intensity <- which(grepl("^Intensity_", names(case)))
  }
  
  case[, Selection := "Removed"]
  case[case$prob > cutoff_prob, ]$Selection <- "Kept"
  
  a <- round(cor(t(case[, index_intensity, with=F]), use="p", method="p"), 2)
  
  if(bool_isPeptideIon == TRUE) {
    rownames(a) <- case$PeptideIon
    colnames(a) <- case$PeptideIon
  } else {
    rownames(a) <- case$aggr_Fragment_Annotation
    colnames(a) <- case$aggr_Fragment_Annotation
  }
  
  a[is.na(a)] <- 0
  
  
  a_reorder <- reorder_cormat(a)
  p_heatmap <- ggplot(melt(a_reorder, value.name = "Correlation", varnames=c("peptideIon", "PeptideIon")), aes(x=peptideIon , y=PeptideIon, fill=Correlation)) + geom_tile() + scale_fill_gradient2(limits=c(-1, 1), low="red", mid="black", high="green", space = "Lab", midpoint=0) + geom_text(aes(peptideIon , PeptideIon, label=Correlation), col="white", size=2.5) + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
  
  if(bool_isPeptideIon == TRUE) {
    
    case$PeptideIon <- factor(case$PeptideIon, levels=case[colnames(a_reorder), on="PeptideIon"]$PeptideIon)
    
    limits <- aes(ymax = case$feature_mean_intensity_all + case$feature_cv_intensity_all * case$feature_mean_intensity_all, ymin = case$feature_mean_intensity_all - case$feature_cv_intensity_all*case$feature_mean_intensity_all) 
    
    #p_bar_int <- ggplot(case, aes(x=PeptideIon, y=feature_mean_intensity_all, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) + guides(fill=FALSE)
    p_bar_int <- ggplot(case, aes(x=PeptideIon, y=feature_mean_intensity_all, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) 
    
    #p_bar_prob <- ggplot(case, aes(x=PeptideIon, y=prob, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() 
    #wenguang: please note that coord_flip() and coord_cartesian() are exclusive.
    p_bar_prob <- ggplot(case, aes(x=PeptideIon, y=prob, fill=Selection)) + geom_bar(stat="identity") + coord_flip(ylim=c(0,1)) 
    
    tmp_injections <- gsub("Intensity_", "", names(case)[grepl("^Intensity", names(case))])
    
    case_long <- melt(case, id.vars=c("PeptideIon", "Selection"), measure.vars = patterns("^Intensity_", "^Score_"), value.name=c("Intensity", "Score"), variable.name="Injections")
    
    setattr(case_long$Injections, "levels", tmp_injections)
    
    case_long[, Extraction := "Original"]
    case_long[case_long$Score == 2.0, ]$Extraction <- "Requant"
    
    
    p_matplot <- ggplot(case_long, aes(x=Injections, y=log2(Intensity))) + geom_line(aes(colour=PeptideIon, group=PeptideIon, linetype=Selection)) + geom_point(aes(shape=Extraction, solid=F), size=3) + geom_text(data=case_long[case_long$Injections==tail(case_long$Injections, n=1), ], aes(colour=PeptideIon, label=PeptideIon, vjust=-2, hjust=1), size=3) + theme(axis.text.x = element_text(angle = 45))
    
    
  } else {
    
    case$aggr_Fragment_Annotation <- factor(case$aggr_Fragment_Annotation, levels=case[colnames(a_reorder), on="aggr_Fragment_Annotation"]$aggr_Fragment_Annotation)
    
    limits <- aes(ymax = case$feature_mean_intensity_all + case$feature_cv_intensity_all * case$feature_mean_intensity_all, ymin = case$feature_mean_intensity_all - case$feature_cv_intensity_all*case$feature_mean_intensity_all) 
    
    #p_bar_int <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=feature_mean_intensity_all, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) + guides(fill=FALSE) 
    p_bar_int <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=feature_mean_intensity_all, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) 
    
    #p_bar_prob <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=prob, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip()
    p_bar_prob <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=prob, fill=Selection)) + geom_bar(stat="identity") + coord_flip(ylim=c(0,1)) 
    
    tmp_injections <- gsub("aggr_Peak_Area_", "", names(case)[grepl("^aggr_Peak_Area", names(case))])
    
    case_long <- melt(case, id.vars=c("aggr_Fragment_Annotation", "PeptideIon", "Selection"), measure.vars = patterns("^aggr_Peak_Area_", "^Score_"), value.name=c("aggr_Peak_Area", "Score"), variable.name="Injections")
    
    setattr(case_long$Injections, "levels", tmp_injections)
    
    case_long[, Extraction := "Original"]
    case_long[case_long$Score == 2.0, ]$Extraction <- "Requant"
    
    
    p_matplot <- ggplot(case_long, aes(x=Injections, y=log2(aggr_Peak_Area))) + geom_line(aes(colour=PeptideIon, group=aggr_Fragment_Annotation, linetype=Selection)) + geom_point(aes(shape=Extraction, solid=F), size=3) + geom_text(data=case_long[case_long$Injections==tail(case_long$Injections, n=1), ], aes(colour=PeptideIon, label=aggr_Fragment_Annotation, vjust=-2, hjust=1), size=3) + theme(axis.text.x = element_text(angle = 45)) 
    
  }
  
  
  p_plot <- grid.arrange( p_matplot, arrangeGrob(p_heatmap, p_bar_int, p_bar_prob, ncol=3), nrow=2 )
  
  
  return(p_plot)
  
}



#' @export
plot_a_few_proteins <- function(input_prot_list, input_dt, input_cutoff_prob=0.3) {
  
  pdf(paste0(substitute(input_prot_list), ".pdf"), width=7.5*3, height=4.1*2)
  
  for(i in 1:length(input_prot_list)) {
    
    case <- input_dt[ProteinName==input_prot_list[i], ]
    
    plot_a_heatmap_include_prob_update(case, input_cutoff_prob)
    
  }
  
  dev.off()
  
}


check_one_protein <- function(case, index_intensity) {
  
  case_long <- melt(case, id.vars="PeptideIon", measure.vars=names(case)[index_intensity], variable.name="index", value.name="intensity")
  p_1 <- ggplot(case_long, aes(x=index, y=log2(intensity))) + geom_line(aes(colour=PeptideIon, group=PeptideIon))
  
  print(p_1)
  
}


check_prot_list <- function(input_prot_list) {

  one_prot_list_all <- vector()
  
  for(i in 1:length(input_prot_list)) {
    temp_num_prot <- as.numeric(sapply(strsplit(input_prot_list[i], "/"), "[[", 1))
    while(temp_num_prot > 0) {
#      cat(i, input_prot_list[i], "\n")
      one_prot_list_all <- c(one_prot_list_all, sapply(strsplit(input_prot_list[i], "/"), "[[", temp_num_prot+1) )
      temp_num_prot <- temp_num_prot - 1
    }
  }
  
  return(unique(one_prot_list_all))
  
}




#' @export

validate_mixed_data <- function(input_dt) {
  
  index_human <- which(grepl("HUMA", input_dt$ProteinName))
  index_yeast <- which(grepl("YEAS", input_dt$ProteinName))
  index_ecoli <- which(grepl("ECOL", input_dt$ProteinName))
  
  cat("human:", length(index_human), length(unique(input_dt[index_human, ]$ProteinName))
      , mean_na( abs(log2(input_dt[index_human, ]$mean_intensity_A / input_dt[index_human, ]$mean_intensity_B)    ))
      , mean_na(input_dt[index_human, ]$cv_intensity_A), mean_na(input_dt[index_human, ]$cv_intensity_B), "\n")
  cat("yeast:", length(index_yeast), length(unique(input_dt[index_yeast, ]$ProteinName))
      , mean_na( abs(log2(input_dt[index_yeast, ]$mean_intensity_A / input_dt[index_yeast, ]$mean_intensity_B) -1 ))
      , mean_na(input_dt[index_yeast, ]$cv_intensity_A), mean_na(input_dt[index_yeast, ]$cv_intensity_B), "\n")
  cat("ecoli:", length(index_ecoli), length(unique(input_dt[index_ecoli, ]$ProteinName))
      , mean_na( abs(log2(input_dt[index_ecoli, ]$mean_intensity_A / input_dt[index_ecoli, ]$mean_intensity_B) +2 ))
      , mean_na(input_dt[index_ecoli, ]$cv_intensity_A), mean_na(input_dt[index_ecoli, ]$cv_intensity_B), "\n")
  
  
  input_dt[, error := 0]
  input_dt[index_human, ]$error <- log2(input_dt[index_human, ]$mean_intensity_A / input_dt[index_human, ]$mean_intensity_B)
  input_dt[index_yeast, ]$error <- log2(input_dt[index_yeast, ]$mean_intensity_A / input_dt[index_yeast, ]$mean_intensity_B) - 1
  input_dt[index_ecoli, ]$error <- log2(input_dt[index_ecoli, ]$mean_intensity_A / input_dt[index_ecoli, ]$mean_intensity_B) + 2
  
  return(input_dt)
  
}

#' @export

validate_HEM_data <- function(input_dt) {
  
  index_human <- which(grepl("HUMA", input_dt$ProteinName))
  index_mouse <- which(grepl("MOUS", input_dt$ProteinName))
  index_ecoli <- which(grepl("ECOL", input_dt$ProteinName))
  
  cat("human:", length(index_human), length(unique(input_dt[index_human, ]$ProteinName))
      , median_na( abs(log2(input_dt[index_human, ]$mean_intensity_C / input_dt[index_human, ]$mean_intensity_E)    ))
      , median_na(input_dt[index_human, ]$cv_intensity_C), median_na(input_dt[index_human, ]$cv_intensity_E), "\n")
  cat("mouse:", length(index_mouse), length(unique(input_dt[index_mouse, ]$ProteinName))
      , median_na( abs(log2(input_dt[index_mouse, ]$mean_intensity_C / input_dt[index_mouse, ]$mean_intensity_E) - log2(3) ))
      , median_na(input_dt[index_mouse, ]$cv_intensity_C), median_na(input_dt[index_mouse, ]$cv_intensity_E), "\n")
  cat("ecoli:", length(index_ecoli), length(unique(input_dt[index_ecoli, ]$ProteinName))
      , median_na( abs(log2(input_dt[index_ecoli, ]$mean_intensity_C / input_dt[index_ecoli, ]$mean_intensity_E) + log2(2) ))
      , median_na(input_dt[index_ecoli, ]$cv_intensity_C), median_na(input_dt[index_ecoli, ]$cv_intensity_E), "\n")
  
  
  input_dt[, error1 := 0]
  input_dt[index_human, ]$error1 <- log2(input_dt[index_human, ]$mean_intensity_C / input_dt[index_human, ]$mean_intensity_E)
  input_dt[index_mouse, ]$error1 <- log2(input_dt[index_mouse, ]$mean_intensity_C / input_dt[index_mouse, ]$mean_intensity_E) - log2(3)
  input_dt[index_ecoli, ]$error1 <- log2(input_dt[index_ecoli, ]$mean_intensity_C / input_dt[index_ecoli, ]$mean_intensity_E) + log2(2)
  
  cat("#########################################\n")
  
  cat("human:", length(index_human), length(unique(input_dt[index_human, ]$ProteinName))
      , median_na( abs(log2(input_dt[index_human, ]$mean_intensity_A / input_dt[index_human, ]$mean_intensity_D)    ))
      , median_na(input_dt[index_human, ]$cv_intensity_A), median_na(input_dt[index_human, ]$cv_intensity_D), "\n")
  cat("mouse:", length(index_mouse), length(unique(input_dt[index_mouse, ]$ProteinName))
      , median_na( abs(log2(input_dt[index_mouse, ]$mean_intensity_A / input_dt[index_mouse, ]$mean_intensity_D) - log2(2) ))
      , median_na(input_dt[index_mouse, ]$cv_intensity_A), median_na(input_dt[index_mouse, ]$cv_intensity_D), "\n")
  cat("ecoli:", length(index_ecoli), length(unique(input_dt[index_ecoli, ]$ProteinName))
      , median_na( abs(log2(input_dt[index_ecoli, ]$mean_intensity_A / input_dt[index_ecoli, ]$mean_intensity_D) + log2(3) ))
      , median_na(input_dt[index_ecoli, ]$cv_intensity_A), median_na(input_dt[index_ecoli, ]$cv_intensity_D), "\n")
  
  
  input_dt[, error2 := 0]
  input_dt[index_human, ]$error2 <- log2(input_dt[index_human, ]$mean_intensity_A / input_dt[index_human, ]$mean_intensity_D)
  input_dt[index_mouse, ]$error2 <- log2(input_dt[index_mouse, ]$mean_intensity_A / input_dt[index_mouse, ]$mean_intensity_D) - log2(2)
  input_dt[index_ecoli, ]$error2 <- log2(input_dt[index_ecoli, ]$mean_intensity_A / input_dt[index_ecoli, ]$mean_intensity_D) + log2(3)
  
  cat("#########################################\n")
  
  cat("human:", length(index_human), length(unique(input_dt[index_human, ]$ProteinName))
      , median_na( abs(log2(input_dt[index_human, ]$mean_intensity_A / input_dt[index_human, ]$mean_intensity_E)    ))
      , median_na(input_dt[index_human, ]$cv_intensity_A), median_na(input_dt[index_human, ]$cv_intensity_E), "\n")
  cat("mouse:", length(index_mouse), length(unique(input_dt[index_mouse, ]$ProteinName))
      , median_na( abs(log2(input_dt[index_mouse, ]$mean_intensity_A / input_dt[index_mouse, ]$mean_intensity_E) - log2(4) ))
      , median_na(input_dt[index_mouse, ]$cv_intensity_A), median_na(input_dt[index_mouse, ]$cv_intensity_E), "\n")
  cat("ecoli:", length(index_ecoli), length(unique(input_dt[index_ecoli, ]$ProteinName))
      , median_na( abs(log2(input_dt[index_ecoli, ]$mean_intensity_A / input_dt[index_ecoli, ]$mean_intensity_E) + log2(4) ))
      , median_na(input_dt[index_ecoli, ]$cv_intensity_A), median_na(input_dt[index_ecoli, ]$cv_intensity_E), "\n")
  
  
  input_dt[, error3 := 0]
  input_dt[index_human, ]$error3 <- log2(input_dt[index_human, ]$mean_intensity_A / input_dt[index_human, ]$mean_intensity_E)
  input_dt[index_mouse, ]$error3 <- log2(input_dt[index_mouse, ]$mean_intensity_A / input_dt[index_mouse, ]$mean_intensity_E) - log2(4)
  input_dt[index_ecoli, ]$error3 <- log2(input_dt[index_ecoli, ]$mean_intensity_A / input_dt[index_ecoli, ]$mean_intensity_E) + log2(4)
  
  return(input_dt)
  
}




#' @export

validate_HEM_protein_table <- function(input_dt) {
  
  result <- data.frame(organism=character(), known_log2FC=numeric(), sample1=numeric(), sample2=numeric(), size=numeric(), median_error=numeric(), stringsAsFactors = FALSE)
  
  known_ecoli <- c(2,3,4,6,8)
  known_mouse <- c(8,7,6,4,2)
  known_ecoli_log2FC <- matrix(0, 5, 5)
  known_mouse_log2FC <- matrix(0, 5, 5)

  for (i in 1:5) {
    for (j in 1:5) {
      known_ecoli_log2FC[i, j] <- log2(known_ecoli[i]) - log2(known_ecoli[j])
      known_mouse_log2FC[i, j] <- log2(known_mouse[i]) - log2(known_mouse[j])
    }
  } 
  
  input_dt[, cor_std := NA]
  
  human <- input_dt[which(grepl("HUMA", input_dt$ProteinName)), ]
  mouse <- input_dt[which(grepl("MOUS", input_dt$ProteinName)), ]
  ecoli <- input_dt[which(grepl("ECOL", input_dt$ProteinName)), ]
  
  sample_list <- c("mean_intensity_A", "mean_intensity_B", "mean_intensity_C", "mean_intensity_D", "mean_intensity_E")
  
  #cat(cor( mouse[, sample_list, with=F],  known_mouse, use="p"))
  cat(median(apply(mouse[, sample_list, with=F], 1, function(x) cor(x, known_mouse, use="p")), na.rm = T), "\n")
  #cat(cor( ecoli[, sample_list, with=F],  known_ecoli, use="p"))
  cat(median(apply(ecoli[, sample_list, with=F], 1, function(x) cor(x, known_ecoli, use="p")), na.rm = T), "\n")
  
  cat( mean(1 - abs(apply(mouse[, sample_list, with=F], 1, function(x) cor(x, known_mouse, use="p"))), na.rm = T), "\n")
  cat( mean(1 - abs(apply(ecoli[, sample_list, with=F], 1, function(x) cor(x, known_ecoli, use="p"))), na.rm = T), "\n")
  
  input_dt$cor_std[which(input_dt$ProteinName %in% mouse$ProteinName)] <- apply(mouse[, sample_list, with=F], 1, function(x) cor(x, known_mouse, use="p"))
  input_dt$cor_std[which(input_dt$ProteinName %in% ecoli$ProteinName)] <- apply(ecoli[, sample_list, with=F], 1, function(x) cor(x, known_ecoli, use="p"))
  
  #input_dt[input_dt$ProteinName %in% ecoli$ProteinName, cor_std := apply(ecoli[, sample_list, with=F], 1, function(x) cor(x, known_ecoli, use="p"))]
  #input_dt[, cor_std := 9999]
  
  for (i in 1:5) {
    for (j in 1:5) {
      
      if(i > j ) {

#        cat("human", i, j, dim(human)[1], median_na( unlist(abs(log2(human[, sample_list[i], with=F]) -
#                                                        log2(human[, sample_list[j], with=F])))), "\n")
#        cat("mouse", known_mouse_log2FC[i, j], i, j, dim(mouse)[1], median_na(unlist(abs(log2(mouse[, sample_list[i], with=F]) -
#                                                                                         log2(mouse[, sample_list[j], with=F]) - known_mouse_log2FC[i, j]))), "\n")
#        cat("ecoli", known_ecoli_log2FC[i, j], i, j, dim(ecoli)[1], median_na(unlist(abs(log2(ecoli[, sample_list[i], with=F]) -
#                                                                                         log2(ecoli[, sample_list[j], with=F]) - known_ecoli_log2FC[i, j]))), "\n")
        
#        result <- rbind(result, c("human", 0, i, j, dim(human)[1], median_na( unlist(abs(log2(human[, sample_list[i], with=F]) -
#                                                                                           log2(human[, sample_list[j], with=F])))) ))
        
        result <- rbind(result, data.frame(organism="human", known_log2FC=0, sample1=i, sample2=j, size=dim(human)[1]
                                           , median_error=median_na( unlist(abs(log2(human[, sample_list[i], with=F]) - log2(human[, sample_list[j], with=F])))) ))
        
        result <- rbind(result, data.frame(organism="mouse", known_log2FC=known_mouse_log2FC[i, j], sample1=i, sample2=j, size=dim(mouse)[1]
                                           , median_error=median_na( unlist(abs(log2(mouse[, sample_list[i], with=F]) - log2(mouse[, sample_list[j], with=F]) - known_mouse_log2FC[i, j]))) ))
        
        result <- rbind(result, data.frame(organism="ecoli", known_log2FC=known_ecoli_log2FC[i, j], sample1=i, sample2=j, size=dim(ecoli)[1]
                                           , median_error=median_na( unlist(abs(log2(ecoli[, sample_list[i], with=F]) - log2(ecoli[, sample_list[j], with=F]) - known_ecoli_log2FC[i, j]))) ))
        
        
      }
    }
  }
  
  result <- result[order(result$organism), ]
  
  #return(result)
  return(input_dt)
  
}







#' @export

validate_HEM_protein_t_test <- function(input_cons_prot_table) {
  
  result <- data.frame(organism=character(), known_log2FC=numeric(), sample1=numeric(), sample2=numeric(), size=numeric(), median_error=numeric(), stringsAsFactors = FALSE)
  
  known_ecoli <- c(2,3,4,6,8)
  known_mouse <- c(8,7,6,4,2)
  known_ecoli_log2FC <- matrix(0, 5, 5)
  known_mouse_log2FC <- matrix(0, 5, 5)
  
  for (i in 1:5) {
    for (j in 1:5) {
      known_ecoli_log2FC[i, j] <- log2(known_ecoli[i]) - log2(known_ecoli[j])
      known_mouse_log2FC[i, j] <- log2(known_mouse[i]) - log2(known_mouse[j])
    }
  } 
  
  index_human <- which(grepl("HUMA", input_cons_prot_table$ProteinName))
  index_mouse <- which(grepl("MOUS", input_cons_prot_table$ProteinName))
  index_ecoli <- which(grepl("ECOL", input_cons_prot_table$ProteinName))
  
  
  sample_list <- c("mean_intensity_A", "mean_intensity_B", "mean_intensity_C", "mean_intensity_D", "mean_intensity_E")
  numNA_list <- c("numNA_intensity_A", "numNA_intensity_B", "numNA_intensity_C", "numNA_intensity_D", "numNA_intensity_E")
  
  for (i in 1:5) {
    for (j in 1:5) {
      
      if(i > j ) {
        
        v_pval <- rep(NA, dim(input_cons_prot_table)[1])
        
        for (k in 1:dim(input_cons_prot_table)[1]) {
          if(input_cons_prot_table[k,numNA_list[i],with=F] < 2 & input_cons_prot_table[k,numNA_list[j],with=F] < 2) {
            v_pval[k] <- t.test(input_cons_prot_table[k,anno[anno$SampleName==unique(anno$SampleName)[i], ]$Injection,with=F], 
                            input_cons_prot_table[k,anno[anno$SampleName==unique(anno$SampleName)[j], ]$Injection,with=F], paired=F)$p.value
          }
        }
        
        #v_pval_adj <- p.adjust(v_pval, method="BH")
        #index_sig <- which(v_pval < 0.01)
        
        v_qval <- qvalue(v_pval[which(v_pval > -9999)])
        index_sig <- which(v_pval > -9999)[which(v_qval$qvalues < 0.01)]
        
        
        result <- rbind(result, data.frame(organism="human", known_log2FC=0, sample1=i, sample2=j, size=length(intersect(index_human, index_sig))
                                           , median_error=median_na( unlist(abs(log2(input_cons_prot_table[intersect(index_human, index_sig), sample_list[i], with=F]) 
                                                                              - log2(input_cons_prot_table[intersect(index_human, index_sig), sample_list[j], with=F])   )) )))
        
        result <- rbind(result, data.frame(organism="mouse", known_log2FC=known_mouse_log2FC[i, j], sample1=i, sample2=j, size=length(intersect(index_mouse, index_sig))
                                           , median_error=median_na( unlist(abs(log2(input_cons_prot_table[intersect(index_mouse, index_sig), sample_list[i], with=F]) 
                                                                              - log2(input_cons_prot_table[intersect(index_mouse, index_sig), sample_list[j], with=F]) - known_mouse_log2FC[i, j])) )))
        
        result <- rbind(result, data.frame(organism="ecoli", known_log2FC=known_ecoli_log2FC[i, j], sample1=i, sample2=j, size=length(intersect(index_ecoli, index_sig))
                                           , median_error=median_na( unlist(abs(log2(input_cons_prot_table[intersect(index_ecoli, index_sig), sample_list[i], with=F]) 
                                                                              - log2(input_cons_prot_table[intersect(index_ecoli, index_sig), sample_list[j], with=F]) - known_ecoli_log2FC[i, j])) )))
        
        cat("one_comparison: ", i, " vs ", j, " num_DE_Prot: ", length(intersect(index_human, index_sig))+length(intersect(index_mouse, index_sig))+length(intersect(index_ecoli, index_sig)), 
            " actual_FDR: ",  length(intersect(index_human, index_sig))/(length(intersect(index_human, index_sig))+length(intersect(index_mouse, index_sig))+length(intersect(index_ecoli, index_sig))), "\n")
        
        
      }
    }
  }
  
  result <- result[order(result$organism), ]
  
  return(result)
  
}








#' @export

validate_HEM_protein_anova <- function(input_cons_prot_table) {
  
  index_human <- which(grepl("HUMA", input_cons_prot_table$ProteinName))
  index_mouse <- which(grepl("MOUS", input_cons_prot_table$ProteinName))
  index_ecoli <- which(grepl("ECOL", input_cons_prot_table$ProteinName))
  
  aov_prot_list <- unique(input_cons_prot_table$ProteinName)
  
  v_pval <- rep(NA, dim(input_cons_prot_table)[1])
  
  test_prot_long <- melt(input_cons_prot_table, id.vars=c("ProteinName", "numPerProt", "protQuantProb"), measure.vars=anno$InjectionName, variable.name="run_id", value.name="Intensity")
  
  for (i in 1: length(aov_prot_list)) {
    
    case <- test_prot_long[ProteinName==aov_prot_list[i], ]
    case <- merge(case, anno, by.x="run_id", by.y="InjectionName")
    
    if(length(which(case[13:15, ]$Intensity > -999)) > 1 & length(which(case[10:12, ]$Intensity > -999)) > 1 ) {
      v_pval[i] <- unlist(summary(aov(Intensity~SampleName, case[7:15, ]))[[1]])[9]
    }
    
  }
  
  
  
  cat("HUMAN", length(intersect( which(v_pval < 0.01), index_human)), "\n")
  cat("MOUSE", length(intersect( which(v_pval < 0.01), index_mouse)), "\n")
  cat("ECOLI", length(intersect( which(v_pval < 0.01), index_ecoli)), "\n")
  
  
}










validate_HEM_data_pairwise <- function(input_dt) {
  
  result <- data.frame(organism=character(), known_log2FC=numeric(), sample1=numeric(), sample2=numeric(), size=numeric(), median_error=numeric(), stringsAsFactors = FALSE)
  
  known_ecoli <- c(2,3,4,6,8)
  known_mouse <- c(8,7,6,4,2)
  known_ecoli_log2FC <- matrix(0, 5, 5)
  known_mouse_log2FC <- matrix(0, 5, 5)
  
  for (i in 1:5) {
    for (j in 1:5) {
      known_ecoli_log2FC[i, j] <- log2(known_ecoli[i]) - log2(known_ecoli[j])
      known_mouse_log2FC[i, j] <- log2(known_mouse[i]) - log2(known_mouse[j])
    }
  } 
  

  
  sample_list <- c("mean_intensity_A", "mean_intensity_B", "mean_intensity_C", "mean_intensity_D", "mean_intensity_E")
  numNA_list <- c("numNA_intensity_A", "numNA_intensity_B", "numNA_intensity_C", "numNA_intensity_D", "numNA_intensity_E")
  
  for (i in 1:5) {
    
    for (j in 1:5) {
      
      if(i > j ) {
        
        filtered <- input_dt[get(numNA_list[i]) < 2 & get(numNA_list[j]) < 2, ]
        
        temp_prot <- merge_replicates(pept2prot(filtered, "prob", 99999), anno)
        
        human <- temp_prot[which(grepl("HUMA", temp_prot$ProteinName)), ]
        mouse <- temp_prot[which(grepl("MOUS", temp_prot$ProteinName)), ]
        ecoli <- temp_prot[which(grepl("ECOL", temp_prot$ProteinName)), ]
        
        #        cat("human", i, j, dim(human)[1], median_na( unlist(abs(log2(human[, sample_list[i], with=F]) -
        #                                                        log2(human[, sample_list[j], with=F])))), "\n")
        #        cat("mouse", known_mouse_log2FC[i, j], i, j, dim(mouse)[1], median_na(unlist(abs(log2(mouse[, sample_list[i], with=F]) -
        #                                                                                         log2(mouse[, sample_list[j], with=F]) - known_mouse_log2FC[i, j]))), "\n")
        #        cat("ecoli", known_ecoli_log2FC[i, j], i, j, dim(ecoli)[1], median_na(unlist(abs(log2(ecoli[, sample_list[i], with=F]) -
        #                                                                                         log2(ecoli[, sample_list[j], with=F]) - known_ecoli_log2FC[i, j]))), "\n")
        
        #        result <- rbind(result, c("human", 0, i, j, dim(human)[1], median_na( unlist(abs(log2(human[, sample_list[i], with=F]) -
        #                                                                                           log2(human[, sample_list[j], with=F])))) ))
        
        result <- rbind(result, data.frame(organism="human", known_log2FC=0, sample1=i, sample2=j, size=dim(human)[1]
                                           , median_error=median_na( unlist(abs(log2(human[, sample_list[i], with=F]) - log2(human[, sample_list[j], with=F])))) ))
        
        result <- rbind(result, data.frame(organism="mouse", known_log2FC=known_mouse_log2FC[i, j], sample1=i, sample2=j, size=dim(mouse)[1]
                                           , median_error=median_na( unlist(abs(log2(mouse[, sample_list[i], with=F]) - log2(mouse[, sample_list[j], with=F]) - known_mouse_log2FC[i, j]))) ))
        
        result <- rbind(result, data.frame(organism="ecoli", known_log2FC=known_ecoli_log2FC[i, j], sample1=i, sample2=j, size=dim(ecoli)[1]
                                           , median_error=median_na( unlist(abs(log2(ecoli[, sample_list[i], with=F]) - log2(ecoli[, sample_list[j], with=F]) - known_ecoli_log2FC[i, j]))) ))
        
        
      }
    }
  }
  
  result <- result[order(result$organism), ]
  
  return(result)
  
}


#' @export

validate_HEM_pairwise_comparison <- function(input_dt) {
  
  result <- data.frame(organism=character(), known_log2FC=numeric(), sample1=numeric(), sample2=numeric(), size=numeric(), median_error=numeric(), stringsAsFactors = FALSE)
  
  known_ecoli <- c(2,3,4,6,8)
  known_mouse <- c(8,7,6,4,2)
  known_ecoli_log2FC <- matrix(0, 5, 5)
  known_mouse_log2FC <- matrix(0, 5, 5)
  
  for (i in 1:5) {
    for (j in 1:5) {
      known_ecoli_log2FC[i, j] <- log2(known_ecoli[i]) - log2(known_ecoli[j])
      known_mouse_log2FC[i, j] <- log2(known_mouse[i]) - log2(known_mouse[j])
    }
  } 
  
  
  
  sample_list <- c("mean_intensity_A", "mean_intensity_B", "mean_intensity_C", "mean_intensity_D", "mean_intensity_E")
  numNA_list <- c("numNA_intensity_A", "numNA_intensity_B", "numNA_intensity_C", "numNA_intensity_D", "numNA_intensity_E")
  
  for (i in 1:5) {
    
    for (j in 1:5) {
      
      if(i > j ) {
        
        temp_prot <- input_dt[, .(median_log2fc = median_na( log2(get(sample_list[i])) - log2(get(sample_list[j])) )), by=ProteinName]
        
        human <- temp_prot[which(grepl("HUMA", temp_prot$ProteinName)), ]
        mouse <- temp_prot[which(grepl("MOUS", temp_prot$ProteinName)), ]
        ecoli <- temp_prot[which(grepl("ECOL", temp_prot$ProteinName)), ]
        
        result <- rbind(result, data.frame(organism="human", known_log2FC=0, sample1=i, sample2=j, size=length(which(human$median_log2fc > -9999999))
                                           , median_error=median_na( unlist(abs( human$median_log2fc - 0 ))) ))
        
        result <- rbind(result, data.frame(organism="mouse", known_log2FC=known_mouse_log2FC[i, j], sample1=i, sample2=j, size=length(which(mouse$median_log2fc > -9999999))
                                           , median_error=median_na( unlist(abs( mouse$median_log2fc - known_mouse_log2FC[i, j]))) ))
        
        result <- rbind(result, data.frame(organism="ecoli", known_log2FC=known_ecoli_log2FC[i, j], sample1=i, sample2=j, size=length(which(ecoli$median_log2fc > -9999999))
                                           , median_error=median_na( unlist(abs( ecoli$median_log2fc - known_ecoli_log2FC[i, j]))) ))
        
        
      }
    }
  }
  
  result <- result[order(result$organism), ]
  
  return(result)
  
}



#' @export

validate_HEM_pairwise_comparison_update <- function(input_dt) {
  
  result <- data.frame(organism=character(), known_log2FC=numeric(), sample1=numeric(), sample2=numeric(), size=numeric(), median_error=numeric(), stringsAsFactors = FALSE)
  
  known_ecoli <- c(2,3,4,6,8)
  known_mouse <- c(8,7,6,4,2)
  known_ecoli_log2FC <- matrix(0, 5, 5)
  known_mouse_log2FC <- matrix(0, 5, 5)
  
  for (i in 1:5) {
    for (j in 1:5) {
      known_ecoli_log2FC[i, j] <- log2(known_ecoli[i]) - log2(known_ecoli[j])
      known_mouse_log2FC[i, j] <- log2(known_mouse[i]) - log2(known_mouse[j])
    }
  } 
  
  
  log2_dt <- copy(input_dt)
  
  sample_list <- c("A", "B", "C", "D", "E")
  
  for(i in 1:length(sample_list)) {
    
    log2_dt[, paste0("mean_log2_intensity_", sample_list[i]) := apply(.SD, 1, function(x) mean_na(log2(x))), 
            .SDcols=anno[anno$SampleName %in% sample_list[i], ]$InjectionName] 

  }

  
  for (i in 1:5) {
    
    for (j in 1:5) {
      
      if(i > j ) {
        
        temp_prot <- log2_dt[, .(median_log2fc = median_na( get(paste0("mean_log2_intensity_", sample_list[i])) - get(paste0("mean_log2_intensity_", sample_list[j])) )), by=ProteinName]
        #temp_prot <- log2_dt[, .(median_log2fc = as.numeric(mean_na( get(paste0("mean_log2_intensity_", sample_list[i])) - get(paste0("mean_log2_intensity_", sample_list[j])) ))), by=ProteinName]
        
        human <- temp_prot[which(grepl("HUMA", temp_prot$ProteinName)), ]
        mouse <- temp_prot[which(grepl("MOUS", temp_prot$ProteinName)), ]
        ecoli <- temp_prot[which(grepl("ECOL", temp_prot$ProteinName)), ]
        
        result <- rbind(result, data.frame(organism="human", known_log2FC=0, sample1=i, sample2=j, size=length(which(human$median_log2fc > -9999999))
                                           , median_error=median_na( unlist(abs( human$median_log2fc - 0 ))) ))
        
        result <- rbind(result, data.frame(organism="mouse", known_log2FC=known_mouse_log2FC[i, j], sample1=i, sample2=j, size=length(which(mouse$median_log2fc > -9999999))
                                           , median_error=median_na( unlist(abs( mouse$median_log2fc - known_mouse_log2FC[i, j]))) ))
        
        result <- rbind(result, data.frame(organism="ecoli", known_log2FC=known_ecoli_log2FC[i, j], sample1=i, sample2=j, size=length(which(ecoli$median_log2fc > -9999999))
                                           , median_error=median_na( unlist(abs( ecoli$median_log2fc - known_ecoli_log2FC[i, j]))) ))
        
        
      }
    }
  }
  
  result <- result[order(result$organism), ]
  
  return(result)
  
}



#' @export

validate_HEM_pairwise_comparison_mapDIA <- function(input_dt) {
  
  result <- data.frame(organism=character(), known_log2FC=numeric(), sample1=numeric(), sample2=numeric(), size=numeric(), median_error=numeric(), stringsAsFactors = FALSE)
  
  known_ecoli <- c(2,3,4,6,8)
  known_mouse <- c(8,7,6,4,2)
  known_ecoli_log2FC <- matrix(0, 5, 5)
  known_mouse_log2FC <- matrix(0, 5, 5)
  
  for (i in 1:5) {
    for (j in 1:5) {
      known_ecoli_log2FC[i, j] <- log2(known_ecoli[i]) - log2(known_ecoli[j])
      known_mouse_log2FC[i, j] <- log2(known_mouse[i]) - log2(known_mouse[j])
    }
  } 
  
  
  log2_dt <- copy(input_dt)
  
  sample_list <- c("A", "B", "C", "D", "E")
  
  for(i in 1:length(anno$InjectionName)) {
    
    log2_dt[, paste0("mean_log2_intensity_", anno$InjectionName[i]) := log2( get(anno$InjectionName[i])  )] 
    
  }
  
  index_inj <- which(grepl("^mean_log2_intensity", names(log2_dt)))
  
  temp_log2 <- copy(log2_dt)
  for(i in 1:length(anno$InjectionName)) {
  
    temp_log2[, index_inj[i]] <- temp_log2[, index_inj[i], with=F]  - apply(log2_dt[, index_inj, with=F], 1, function(x) median(x, na.rm=T))
  
  }
  
  for (i in 1:5) {
    
    for (j in 1:5) {
      
      if(i > j ) {
        
        temp_prot <- log2_dt[, .(median_log2fc = median_na( get(paste0("mean_log2_intensity_", sample_list[i])) - get(paste0("mean_log2_intensity_", sample_list[j])) )), by=ProteinName]
        #temp_prot <- log2_dt[, .(median_log2fc = as.numeric(mean_na( get(paste0("mean_log2_intensity_", sample_list[i])) - get(paste0("mean_log2_intensity_", sample_list[j])) ))), by=ProteinName]
        
        human <- temp_prot[which(grepl("HUMA", temp_prot$ProteinName)), ]
        mouse <- temp_prot[which(grepl("MOUS", temp_prot$ProteinName)), ]
        ecoli <- temp_prot[which(grepl("ECOL", temp_prot$ProteinName)), ]
        
        result <- rbind(result, data.frame(organism="human", known_log2FC=0, sample1=i, sample2=j, size=length(which(human$median_log2fc > -9999999))
                                           , median_error=median_na( unlist(abs( human$median_log2fc - 0 ))) ))
        
        result <- rbind(result, data.frame(organism="mouse", known_log2FC=known_mouse_log2FC[i, j], sample1=i, sample2=j, size=length(which(mouse$median_log2fc > -9999999))
                                           , median_error=median_na( unlist(abs( mouse$median_log2fc - known_mouse_log2FC[i, j]))) ))
        
        result <- rbind(result, data.frame(organism="ecoli", known_log2FC=known_ecoli_log2FC[i, j], sample1=i, sample2=j, size=length(which(ecoli$median_log2fc > -9999999))
                                           , median_error=median_na( unlist(abs( ecoli$median_log2fc - known_ecoli_log2FC[i, j]))) ))
        
        
      }
    }
  }
  
  result <- result[order(result$organism), ]
  
  return(result)
  
}




#' @export
validate_HEM_pairwise_comparison_peca <- function(input_dt) {
  
  result <- data.frame(organism=character(), known_log2FC=numeric(), sample1=numeric(), sample2=numeric(), 
                             size=numeric(), median_error=numeric(), stringsAsFactors = FALSE)
  
  
  known_ecoli <- c(2,3,4,6,8)
  known_mouse <- c(8,7,6,4,2)
  known_ecoli_log2FC <- matrix(0, 5, 5)
  known_mouse_log2FC <- matrix(0, 5, 5)
  
  for (i in 1:5) {
    for (j in 1:5) {
      known_ecoli_log2FC[i, j] <- log2(known_ecoli[i]) - log2(known_ecoli[j])
      known_mouse_log2FC[i, j] <- log2(known_mouse[i]) - log2(known_mouse[j])
    }
  } 
  
  
  for (i in 1:5) {
    
    for (j in 1:5) {
      
      if(i > j ) {
        
        sample_list <- c("A", "B", "C", "D", "E")
        
        analysis <- input_dt[input_dt$Label==paste0(sample_list[i], "/", sample_list[j]), ]
        
        human <- analysis[grepl("HUMA", analysis$ProteinName), ]
        ecoli <- analysis[grepl("ECOL", analysis$ProteinName), ]
        mouse <- analysis[grepl("MOUS", analysis$ProteinName), ]
        
        
        result <- rbind(result, data.frame(organism="human", known_log2FC=0, sample1=i, sample2=j, size=dim(human)[1]
                                           , median_error=median_na( unlist(abs( human$slr - 0 ))) ))
        
        result <- rbind(result, data.frame(organism="mouse", known_log2FC=known_mouse_log2FC[i, j], sample1=i, sample2=j, size=dim(mouse)[1]
                                           , median_error=median_na( unlist(abs( mouse$slr - known_mouse_log2FC[i, j]))) ))
        
        result <- rbind(result, data.frame(organism="ecoli", known_log2FC=known_ecoli_log2FC[i, j], sample1=i, sample2=j, size=dim(ecoli)[1]
                                           , median_error=median_na( unlist(abs( ecoli$slr - known_ecoli_log2FC[i, j]))) ))
        
        
        
        cat("one_comparison: ", i, " vs ", j, " num_DE_Prot: ", dim(human)[1]+dim(mouse)[1]+dim(ecoli)[1], 
            " actual_FDR: ",  dim(human)[1]/(dim(human)[1]+dim(mouse)[1]+dim(ecoli)[1]), "\n")
        
        
        human <- analysis[grepl("HUMA", analysis$Protein), ]
        mouse <- analysis[grepl("MOUS", analysis$Protein), ]
        ecoli <- analysis[grepl("ECOL", analysis$Protein), ]

        
      }
      
    }
  }
  
  #result <- result[order(result$organism), ]
  
  return(result)
  
}




#' @export
validate_HYE_lumos_pairwise_comparison_peca <- function(input_dt) {
  
  result <- data.frame(organism=character(), known_log2FC=numeric(), sample1=numeric(), sample2=numeric(), 
                       size=numeric(), median_error=numeric(), stringsAsFactors = FALSE)
  
  
  known_ecoli <- c(1, 10)
  known_yeast <- c(10, 1)
  known_ecoli_log2FC <- matrix(0, 2, 2)
  known_yeast_log2FC <- matrix(0, 2, 2)
  
  for (i in 1:2) {
    for (j in 1:2) {
      known_ecoli_log2FC[i, j] <- log2(known_ecoli[i]) - log2(known_ecoli[j])
      known_yeast_log2FC[i, j] <- log2(known_yeast[i]) - log2(known_yeast[j])
    }
  } 
  
  
  for (i in 1:2) {
    
    for (j in 1:2) {
      
      if(i > j ) {
        
        sample_list <- c("A", "B")
        
        analysis <- input_dt[input_dt$Label==paste0(sample_list[i], "/", sample_list[j]), ]
        
        human <- analysis[grepl("HUMA", analysis$ProteinName), ]
        ecoli <- analysis[grepl("ECOL", analysis$ProteinName), ]
        yeast <- analysis[grepl("YEAS", analysis$ProteinName), ]
        
        
        result <- rbind(result, data.frame(organism="human", known_log2FC=0, sample1=i, sample2=j, size=dim(human)[1]
                                           , median_error=median_na( unlist(abs( human$slr - 0 ))) ))
        
        result <- rbind(result, data.frame(organism="yeast", known_log2FC=known_yeast_log2FC[i, j], sample1=i, sample2=j, size=dim(yeast)[1]
                                           , median_error=median_na( unlist(abs( yeast$slr - known_yeast_log2FC[i, j]))) ))
        
        result <- rbind(result, data.frame(organism="ecoli", known_log2FC=known_ecoli_log2FC[i, j], sample1=i, sample2=j, size=dim(ecoli)[1]
                                           , median_error=median_na( unlist(abs( ecoli$slr - known_ecoli_log2FC[i, j]))) ))
        
        
        
        cat("one_comparison: ", i, " vs ", j, " num_DE_Prot: ", dim(human)[1]+dim(yeast)[1]+dim(ecoli)[1], 
            " actual_FDR: ",  dim(human)[1]/(dim(human)[1]+dim(yeast)[1]+dim(ecoli)[1]), "\n")
        
        
        human <- analysis[grepl("HUMA", analysis$Protein), ]
        yeast <- analysis[grepl("MOUS", analysis$Protein), ]
        ecoli <- analysis[grepl("ECOL", analysis$Protein), ]
        
        
      }
      
    }
  }
  
  #result <- result[order(result$organism), ]
  
  return(result)
  
}


#' @export
check_a_column <- function(input_vector) {
 
  cat("length:", length(input_vector), "; max:", max(input_vector, na.rm=T), 
                                       "; min:", min(input_vector, na.rm=T), 
                                       "; median:", median(input_vector, na.rm=T),
                                       "; mean:", mean(input_vector, na.rm=T), "\n")
  cat("#number:", length(which(input_vector > -99999999999 & input_vector < 99999999999)), "\n")
  cat("#na:", length(which(is.na(input_vector))), "; #inf:", length(which(is.infinite(input_vector))), "; #nan:", length(which(is.nan(input_vector))), "\n")
  cat("if:                                                ", length(input_vector) == length(which(input_vector > -99999999999 & input_vector < 99999999999)) 
                                + length(which(is.na(input_vector))) + length(which(is.infinite(input_vector))) + length(which(is.nan(input_vector))), "\n")
   
}



my_test_cor <- function (x, y = NULL, use = "everything", method = c("pearson", 
                                                      "kendall", "spearman")) 
{
  
  C_cor=get("C_cor", asNamespace("stats"))
  
  na.method <- pmatch(use, c("all.obs", "complete.obs", "pairwise.complete.obs", 
                             "everything", "na.or.complete"))
  if (is.na(na.method)) 
    stop("invalid 'use' argument")
  method <- match.arg(method)
  if (is.data.frame(y)) 
    y <- as.matrix(y)
  if (is.data.frame(x)) 
    x <- as.matrix(x)
  if (!is.matrix(x) && is.null(y)) 
    stop("supply both 'x' and 'y' or a matrix-like 'x'")
  if (!(is.numeric(x) || is.logical(x))) 
    stop("'x' must be numeric")
  stopifnot(is.atomic(x))
  if (!is.null(y)) {
    if (!(is.numeric(y) || is.logical(y))) 
      stop("'y' must be numeric")
    stopifnot(is.atomic(y))
  }
  Rank <- function(u) {
    if (length(u) == 0L) 
      u
    else if (is.matrix(u)) {
      if (nrow(u) > 1L) 
        apply(u, 2L, rank, na.last = "keep")
      else row(u)
    }
    else rank(u, na.last = "keep")
  }
  if (method == "pearson") 
    .Call(C_cor, x, y, na.method, FALSE)
  else if (na.method %in% c(2L, 5L)) {
    if (is.null(y)) {
      .Call(C_cor, Rank(na.omit(x)), NULL, na.method, method == 
              "kendall")
    }
    else {
      nas <- attr(na.omit(cbind(x, y)), "na.action")
      dropNA <- function(x, nas) {
        if (length(nas)) {
          if (is.matrix(x)) 
            x[-nas, , drop = FALSE]
          else x[-nas]
        }
        else x
      }
      .Call(C_cor, Rank(dropNA(x, nas)), Rank(dropNA(y, 
                                                     nas)), na.method, method == "kendall")
    }
  }
  else if (na.method != 3L) {
    x <- Rank(x)
    if (!is.null(y)) 
      y <- Rank(y)
    .Call(C_cor, x, y, na.method, method == "kendall")
  }
  else {
    if (is.null(y)) {
      ncy <- ncx <- ncol(x)
      if (ncx == 0) 
        stop("'x' is empty")
      r <- matrix(0, nrow = ncx, ncol = ncy)
      for (i in seq_len(ncx)) {
        for (j in seq_len(i)) {
          x2 <- x[, i]
          y2 <- x[, j]
          ok <- complete.cases(x2, y2)
          x2 <- rank(x2[ok])
          y2 <- rank(y2[ok])
          r[i, j] <- if (any(ok)) 
            .Call(C_cor, x2, y2, 1L, method == "kendall")
          else NA
        }
      }
      r <- r + t(r) - diag(diag(r))
      rownames(r) <- colnames(x)
      colnames(r) <- colnames(x)
      r
    }
    else {
      if (length(x) == 0L || length(y) == 0L) 
        stop("both 'x' and 'y' must be non-empty")
      matrix_result <- is.matrix(x) || is.matrix(y)
      if (!is.matrix(x)) 
        x <- matrix(x, ncol = 1L)
      if (!is.matrix(y)) 
        y <- matrix(y, ncol = 1L)
      ncx <- ncol(x)
      ncy <- ncol(y)
      r <- matrix(0, nrow = ncx, ncol = ncy)
      for (i in seq_len(ncx)) {
        for (j in seq_len(ncy)) {
          x2 <- x[, i]
          y2 <- y[, j]
          ok <- complete.cases(x2, y2)
          x2 <- rank(x2[ok])
          y2 <- rank(y2[ok])
          r[i, j] <- if (any(ok)) 
            .Call(C_cor, x2, y2, 1L, method == "kendall")
          else NA
        }
      }
      rownames(r) <- colnames(x)
      colnames(r) <- colnames(y)
      if (matrix_result) 
        r
      else drop(r)
    }
  }
}



















.f <- function() {

`PECA_tsv` <- function(file=NULL, samplenames1=NULL, samplenames2=NULL, normalize=FALSE, test="t", type="median", paired=FALSE, progress=FALSE) {
  # Read tsv-file
  message("Reading data")
  flush.console()
  probeintensities <- read.csv(file, sep="\t")
  probenamesGene <- probeintensities[,1]
  probeintensities <- subset(probeintensities,select=c(samplenames1,samplenames2))
  probeintensities <- as.matrix(probeintensities)
  PECA(probenamesGene=probenamesGene, probeintensities=probeintensities, samplenames1=samplenames1, samplenames2=samplenames2, normalize=normalize, test=test, type=type, paired=paired, progress=progress)
}



# PECA function called by different wrappers
`PECA` <- function(probenamesGene=NULL, probeintensities=NULL, samplenames1=NULL, samplenames2=NULL, test="t", paired=FALSE) {
  
  
  # Log transformation
  message("Performing log-transformation")
  flush.console()
  probeintensities <- probeintensities + 1
  probeintensities <- log2(probeintensities)
  colnames(probeintensities) <- c(samplenames1,samplenames2)
  
  # PECA slr and t-statistic
  message("Calculating low-level statistics")
  flush.console()
  if (test == "t") {
    if (paired) {
      probeSLR <- matrix(nrow=nrow(probeintensities), ncol=length(samplenames1)) 
      for(i in 1:length(samplenames1)) probeSLR[,i] <- probeintensities[,samplenames1[i]] - probeintensities[,samplenames2[i]]
      t <- rowttests(probeSLR, tstatOnly=TRUE)
      df.total <- length(samplenames1)-1
    }
    else {
      labels <- factor(c(rep(1,length(samplenames1)), rep(2,length(samplenames2))))    
      t <- rowttests(as.matrix(cbind(probeintensities[,samplenames1], probeintensities[,samplenames2])), fac=labels, tstatOnly=TRUE)
      df.total <- length(samplenames1)+length(samplenames2)-2
    }
    probeSLR <- t$dm
    t <- t$statistic
  }
  if (test == "modt" | test == "rots") {
    if (paired) {
      probeSLR <- matrix(nrow=nrow(probeintensities), ncol=length(samplenames1)) 
      for(i in 1:length(samplenames1)) probeSLR[,i] <- probeintensities[,samplenames1[i]] - probeintensities[,samplenames2[i]]
      fit <- lmFit(probeSLR)
      fit <- eBayes(fit)
      probeSLR <- fit$coefficients
      t <- fit$t
    } else {
      design <- cbind(G1=1,G1vsG2=c(rep(1,length(samplenames1)), rep(0,length(samplenames2))))
      probeSLR <- as.matrix(cbind(probeintensities[,samplenames1], probeintensities[,samplenames2]))
      fit <- lmFit(probeSLR, design)
      fit <- eBayes(fit)
      probeSLR <- fit$coefficients[,2]
      t <- fit$t[,2]
    }
    df.total<- fit$df.residual[1] + fit$df.prior
    rm(fit)
    gc()
  }

  
  # Aggregating statistics
  message("Aggregating statistics")
  flush.console()
  gene.n <- tapply(t, probenamesGene, function(x) sum(!is.na(x)))
  if (type=="median") {
    geneSLR <- tapply(probeSLR, probenamesGene, median, na.rm=TRUE)
    t <- tapply(t, probenamesGene, median, na.rm=TRUE)
    if (test == "rots") {rots.p <- 1- abs(tapply(rots.p, probenamesGene, median, na.rm=TRUE))}
  }
  if (type=="tukey") {
    geneSLR <- tapply(probeSLR, probenamesGene, tukey)
    t <- tapply(t, probenamesGene, tukey)
    if (test == "rots") {rots.p <- 1 - abs(tapply(rots.p, probenamesGene, tukey))}
  }
  
  # P-values
  gene.p <- 2 * pt(abs(t), df=df.total, lower.tail=FALSE)
  if (test == "rots") {gene.p <- rots.p}
  gene.p2 <- gene.p
  if (type=="median") {
    gene.p2 <- pbeta(gene.p, gene.n/2 + 0.5, gene.n - (gene.n/2 + 0.5) + 1)
  }
  if (type=="tukey") {
    message("Simulating distributions")
    flush.console()
    distributions <- generateDistributions(max=max(gene.n), k=10000)
    gene.p2 <- mapply(psim, p=gene.p, list=distributions[gene.n])
  }
  gene.p.fdr <- p.adjust(gene.p2, method="BH")
  
  # Cleanup
  rm(probeintensities)
  rm(probeSLR)
  gc()
  
  # Return a table containing slr, t-statistic and p-value
  result <- data.frame(cbind(log2fc=geneSLR, t=t, score=gene.p, nunPerProt=gene.n, pval=gene.p2, padj=gene.p.fdr))
  message("Done")
  return(result)
  
}



}
