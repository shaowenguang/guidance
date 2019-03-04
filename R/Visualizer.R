#' Visualization of protein abundance 
#' 
#' 
#' @export
plot_protein_profile <- function(case, cutoff_prob=0.3) {
  
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
  p_heatmap <- ggplot(melt(a_reorder, value.name = "Correlation", 
                           varnames=c("peptideIon", "PeptideIon")), 
                      aes(x=peptideIon , y=PeptideIon, fill=Correlation)) + 
    geom_tile() + scale_fill_gradient2(limits=c(-1, 1), low="red", mid="black", 
                                       high="green", space = "Lab", midpoint=0) + 
    geom_text(aes(peptideIon , PeptideIon, label=Correlation), col="white", size=2.5) + 
    theme(axis.text.x=element_blank(), axis.text.y=element_blank())
  
  if(bool_isPeptideIon == TRUE) {
    
    case$PeptideIon <- factor(case$PeptideIon, levels=case[colnames(a_reorder), 
                                                           on="PeptideIon"]$PeptideIon)
    
    limits <- aes(ymax = case$feature_mean_intensity_all + 
                    case$feature_cv_intensity_all * case$feature_mean_intensity_all, 
                  ymin = case$feature_mean_intensity_all - 
                    case$feature_cv_intensity_all*case$feature_mean_intensity_all) 
    
    #p_bar_int <- ggplot(case, aes(x=PeptideIon, y=feature_mean_intensity_all, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) + guides(fill=FALSE)
    p_bar_int <- ggplot(case, aes(x=PeptideIon, y=feature_mean_intensity_all, fill=ProteinName)) + 
      geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) 
    
    #p_bar_prob <- ggplot(case, aes(x=PeptideIon, y=prob, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() 
    #wenguang: please note that coord_flip() and coord_cartesian() are exclusive.
    p_bar_prob <- ggplot(case, aes(x=PeptideIon, y=prob, fill=Selection)) + 
      geom_bar(stat="identity") + coord_flip(ylim=c(0,1)) 
    
    tmp_injections <- gsub("Intensity_", "", names(case)[grepl("^Intensity", names(case))])
    
    case_long <- melt(case, id.vars=c("PeptideIon", "Selection"), 
                      measure.vars = patterns("^Intensity_", "^Score_"), 
                      value.name=c("Intensity", "Score"), variable.name="Injections")
    
    setattr(case_long$Injections, "levels", tmp_injections)
    
    case_long[, Extraction := "Original"]
    case_long[case_long$Score == 2.0, ]$Extraction <- "Requant"
    
    
    p_matplot <- ggplot(case_long, aes(x=Injections, y=log2(Intensity))) + 
      geom_line(aes(colour=PeptideIon, group=PeptideIon, linetype=Selection)) + 
      geom_point(aes(shape=Extraction, solid=F), size=3) + 
      geom_text(data=case_long[case_long$Injections==tail(case_long$Injections, n=1), ], 
                aes(colour=PeptideIon, label=PeptideIon, vjust=-2, hjust=1), size=3) + 
      theme(axis.text.x = element_text(angle = 45))
    
    
  } else {
    
    case$aggr_Fragment_Annotation <- factor(case$aggr_Fragment_Annotation, 
                                            levels=case[colnames(a_reorder), 
                                                        on="aggr_Fragment_Annotation"]$aggr_Fragment_Annotation)
    
    limits <- aes(ymax = case$feature_mean_intensity_all + 
                    case$feature_cv_intensity_all * case$feature_mean_intensity_all, 
                  ymin = case$feature_mean_intensity_all - 
                    case$feature_cv_intensity_all*case$feature_mean_intensity_all) 
    
    #p_bar_int <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=feature_mean_intensity_all, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) + guides(fill=FALSE) 
    p_bar_int <- ggplot(case, aes(x=aggr_Fragment_Annotation, 
                                  y=feature_mean_intensity_all, fill=ProteinName)) +
      geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) 
    
    #p_bar_prob <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=prob, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip()
    p_bar_prob <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=prob, fill=Selection)) + 
      geom_bar(stat="identity") + coord_flip(ylim=c(0,1)) 
    
    tmp_injections <- gsub("aggr_Peak_Area_", "", 
                           names(case)[grepl("^aggr_Peak_Area", names(case))])
    
    case_long <- melt(case, id.vars=c("aggr_Fragment_Annotation", "PeptideIon", "Selection"), 
                      measure.vars = patterns("^aggr_Peak_Area_", "^Score_"), 
                      value.name=c("aggr_Peak_Area", "Score"), variable.name="Injections")
    
    setattr(case_long$Injections, "levels", tmp_injections)
    
    case_long[, Extraction := "Original"]
    case_long[case_long$Score == 2.0, ]$Extraction <- "Requant"
    
    
    p_matplot <- ggplot(case_long, aes(x=Injections, y=log2(aggr_Peak_Area))) + 
      geom_line(aes(colour=PeptideIon, group=aggr_Fragment_Annotation, 
                    linetype=Selection)) + 
      geom_point(aes(shape=Extraction, solid=F), size=3) + 
      geom_text(data=case_long[case_long$Injections==tail(case_long$Injections, n=1), ], 
                aes(colour=PeptideIon, label=aggr_Fragment_Annotation, 
                    vjust=-2, hjust=1), size=3) + 
      theme(axis.text.x = element_text(angle = 45)) 
    
  }
  
  
  p_plot <- grid.arrange( p_matplot, 
                          arrangeGrob(p_heatmap, p_bar_int, p_bar_prob, ncol=3), nrow=2 )
  
  
  return(p_plot)
  
}



#' @export
plot_peptide_intensity <- function(case, cutoff_prob=0.3) {
  
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
  
  if(bool_isPeptideIon == TRUE) {
    
    case$PeptideIon <- factor(case$PeptideIon, levels=case[colnames(a_reorder), 
                                                           on="PeptideIon"]$PeptideIon)
   
    tmp_injections <- gsub("Intensity_", "", names(case)[grepl("^Intensity", names(case))])
    
    case_long <- melt(case, id.vars=c("PeptideIon", "Selection"), 
                      measure.vars = patterns("^Intensity_", "^Score_"), 
                      value.name=c("Intensity", "Score"), variable.name="Injections")
    
    setattr(case_long$Injections, "levels", tmp_injections)
    
    case_long[, Extraction := "Original"]
    case_long[case_long$Score == 2.0, ]$Extraction <- "Requant"
    
    
    p_matplot <- ggplot(case_long, aes(x=Injections, y=log2(Intensity))) + 
      geom_line(aes(colour=PeptideIon, group=PeptideIon, linetype=Selection)) + 
      geom_point(aes(shape=Extraction, solid=F), size=3) + 
      geom_text(data=case_long[case_long$Injections==tail(case_long$Injections, n=1), ], 
                aes(colour=PeptideIon, label=PeptideIon, vjust=-2, hjust=1), size=3) + 
      theme(axis.text.x = element_text(angle = 45))
    
    
  } else {
    
    case$aggr_Fragment_Annotation <- factor(case$aggr_Fragment_Annotation, 
                                            levels=case[colnames(a_reorder), 
                                                        on="aggr_Fragment_Annotation"]$aggr_Fragment_Annotation)
    
    tmp_injections <- gsub("aggr_Peak_Area_", "", 
                           names(case)[grepl("^aggr_Peak_Area", names(case))])
    
    case_long <- melt(case, id.vars=c("aggr_Fragment_Annotation", "PeptideIon", "Selection"), 
                      measure.vars = patterns("^aggr_Peak_Area_", "^Score_"), 
                      value.name=c("aggr_Peak_Area", "Score"), variable.name="Injections")
    
    setattr(case_long$Injections, "levels", tmp_injections)
    
    case_long[, Extraction := "Original"]
    case_long[case_long$Score == 2.0, ]$Extraction <- "Requant"
    
    
    p_matplot <- ggplot(case_long, aes(x=Injections, y=log2(aggr_Peak_Area))) + 
      geom_line(aes(colour=PeptideIon, group=aggr_Fragment_Annotation, 
                    linetype=Selection)) + 
      geom_point(aes(shape=Extraction, solid=F), size=3) + 
      geom_text(data=case_long[case_long$Injections==tail(case_long$Injections, n=1), ], 
                aes(colour=PeptideIon, label=aggr_Fragment_Annotation, 
                    vjust=-2, hjust=1), size=3) + 
      theme(axis.text.x = element_text(angle = 45)) 
  } 
  
  return(p_matplot)
  
}

#' @export
plot_cor_heatmap <- function(case, cutoff_prob=0.3) {
  
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
  p_heatmap <- ggplot(melt(a_reorder, value.name = "Correlation", 
                           varnames=c("peptideIon", "PeptideIon")), 
                      aes(x=peptideIon , y=PeptideIon, fill=Correlation)) + 
    geom_tile() + scale_fill_gradient2(limits=c(-1, 1), low="red", mid="black", 
                                       high="green", space = "Lab", midpoint=0) + 
    geom_text(aes(peptideIon , PeptideIon, label=Correlation), col="white", size=2.5) + 
    theme(axis.text.x=element_blank())
  
 
  return(p_heatmap)
  
}

#' @export
plot_bar_intensity_n_probability <- function(case, cutoff_prob=0.3, plot_intensity = T, 
                                               plot_prob = T) {
  
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
  
  if(bool_isPeptideIon == TRUE) {
    
    case$PeptideIon <- factor(case$PeptideIon, levels=case[colnames(a_reorder), 
                                                           on="PeptideIon"]$PeptideIon)
    
    limits <- aes(ymax = case$feature_mean_intensity_all + 
                    case$feature_cv_intensity_all * case$feature_mean_intensity_all, 
                  ymin = case$feature_mean_intensity_all - 
                    case$feature_cv_intensity_all*case$feature_mean_intensity_all) 
    
    if(plot_intensity == TRUE){
      #p_bar_int <- ggplot(case, aes(x=PeptideIon, y=feature_mean_intensity_all, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) + guides(fill=FALSE)
      p_bar_int <- ggplot(case, aes(x=PeptideIon, y=feature_mean_intensity_all, fill=ProteinName)) + 
        geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) 
    }
    
    if(plot_prob == TRUE){
      #p_bar_prob <- ggplot(case, aes(x=PeptideIon, y=prob, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() 
      #wenguang: please note that coord_flip() and coord_cartesian() are exclusive.
      p_bar_prob <- ggplot(case, aes(x=PeptideIon, y=prob, fill=Selection)) + 
        geom_bar(stat="identity") + coord_flip(ylim=c(0,1)) 
    }

  } else {
    
    case$aggr_Fragment_Annotation <- factor(case$aggr_Fragment_Annotation, 
                                            levels=case[colnames(a_reorder), 
                                                        on="aggr_Fragment_Annotation"]$aggr_Fragment_Annotation)
    
    limits <- aes(ymax = case$feature_mean_intensity_all + 
                    case$feature_cv_intensity_all * case$feature_mean_intensity_all, 
                  ymin = case$feature_mean_intensity_all - 
                    case$feature_cv_intensity_all*case$feature_mean_intensity_all) 
    
    if(plot_intensity == TRUE){
      #p_bar_int <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=feature_mean_intensity_all, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) + guides(fill=FALSE) 
      p_bar_int <- ggplot(case, aes(x=aggr_Fragment_Annotation, 
                                    y=feature_mean_intensity_all, fill=ProteinName)) +
        geom_bar(stat="identity") + coord_flip() + geom_errorbar(limits, width=0.5) 
      p_plot <- p_bar_int
      }
    
    if(plot_prob == TRUE){
      #p_bar_prob <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=prob, fill=ProteinName)) + geom_bar(stat="identity") + coord_flip()
      p_bar_prob <- ggplot(case, aes(x=aggr_Fragment_Annotation, y=prob, fill=Selection)) + 
        geom_bar(stat="identity") + coord_flip(ylim=c(0,1))
      p_plot <- plot_prob
      }
  }
  
  if(plot_intensity == TRUE && plot_prob == TRUE){
    p_plot <- grid.arrange(p_bar_int, p_bar_prob)
  }
  
  
  return(p_plot)
  
}


#' Plot density distribution of feature
#' 
#' 
#' @export
plot_density <- function(data, feature = "feature_mean_intensity_all", fill = NULL,
                         font.size = 12, title = feature){
  
  if(feature == "feature_mean_intensity_all"){
    data[[feature]] <- log2(data[[feature]])
  }
  
  ggplot(data, aes_string(x=feature, fill=fill)) + 
    geom_density(alpha=0.8)  + theme(
      axis.text=element_text(size=font.size), 
      axis.title=element_text(size=font.size),    
      legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    ) + labs(fill="") + 
    theme(axis.line.x = element_line(color="black"), 
          axis.line.y = element_line(color="black")) +
    ggtitle(title)
}
