#' Visualize peptide profiles 
#' 
#' @param case data table or data frame in wide representation. The data typically 
#' contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
#' measurements of each peptide or precursor ions in rows. 
#' 
#' For a clear visualization of peptide profile, recommend to prepare a 
#' data table pre-filter 
#' peptides corresponding to a single protein  
#' 
#' @param cutoff_prob a numeric value denoting posterior probability threshold
#' to keep or remove peptides. 
#'
#' @export
#' 
#' @examples 
#' peptideIons_features <- calc_features(all_peptideIons)
#' test <- perform_selection(peptideIons_features)
#' 
#' prot_name <- c("1/O75976")
#' test_prot <- test[test$ProteinName==prot_name, ]
#' plot_protein_profile(test_prot)
#' 
plot_protein_profile <- function(case, cutoff_prob=0.3) {

  a_reorder <- compute_cor(case = case, cutoff_prob=cutoff_prob)
  
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



#' Visualize peptide intensity 
#' 
#' @param case data table or data frame in wide representation. The data typically 
#' contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
#' measurements of each peptide or precursor ions in rows. 
#' 
#' To visualize profile of peptides corresponding to the same protein, 
#' peptides corresponding to one protein are pre-filtered for clearer 
#' visualizations. 
#'
#' @param cutoff_prob a numeric value denoting posterior probability threshold
#' to keep or remove peptides. 
#'
#' @export
#' 
#' @examples 
#' peptideIons_features <- calc_features(all_peptideIons)
#' peptideIons_features_select <- perform_selection(peptideIons_features)
#' 
#' prot_name <- c("1/O75976")
#' test_prot <- test[test$ProteinName==prot_name, ]
#' plot_peptide_intensity(test_prot)
#' 
plot_peptide_intensity <- function(case, cutoff_prob=0.3) {
  
  a_reorder <- compute_cor(case = case, cutoff_prob=cutoff_prob)
  
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

#' Visualize correlation heatmap 
#' 
#' @param case data table or data frame in wide representation. The data typically 
#' contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
#' measurements of each peptide or precursor ions in rows. 
#' 
#' To visualize profile of peptides corresponding to the same protein, 
#' peptides corresponding to one protein are pre-filtered for clearer 
#' visualizations. 
#'
#' @param cutoff_prob a numeric value denoting posterior probability threshold
#' to keep or remove peptides. 
#'
#' @export
#' 
#' @examples 
#' peptideIons_features <- calc_features(all_peptideIons)
#' peptideIons_features_select <- perform_selection(peptideIons_features)
#' 
#' prot_name <- c("1/O75976")
#' test_prot <- test[test$ProteinName==prot_name, ]
#' plot_cor_heatmap(test_prot)
#' 
plot_cor_heatmap <- function(case, cutoff_prob=0.3) {
  
  a_reorder <- compute_cor(case = case, cutoff_prob=cutoff_prob)
  
  p_heatmap <- ggplot(melt(a_reorder, value.name = "Correlation", 
                           varnames=c("peptideIon", "PeptideIon")), 
                      aes(x=peptideIon , y=PeptideIon, fill=Correlation)) + 
    geom_tile() + scale_fill_gradient2(limits=c(-1, 1), low="red", mid="black", 
                                       high="green", space = "Lab", midpoint=0) + 
    geom_text(aes(peptideIon , PeptideIon, label=Correlation), col="white", size=2.5) + 
    theme(axis.text.x=element_blank())
  
 
  return(p_heatmap)
  
}

#' Visualize intensity and posterior probability in barplots 
#' 
#' @param case data table or data frame in wide representation. The data typically 
#' contains \code{"PeptideIon"}, \code{"ProteinName"} and sample names in columns and 
#' measurements of each peptide or precursor ions in rows. 
#' 
#' To visualize profile of peptides corresponding to the same protein, 
#' peptides corresponding to one protein are pre-filtered for clearer 
#' visualizations. 
#'
#' @param cutoff_prob a numeric value denoting posterior probability threshold
#' to keep or remove peptides. 
#' @param plot_intensity a logical indicating whether to plot intensity barplots
#' @param plot_prob a logical indicating whether to plot posterior probability
#' barplots 
#'
#' @export
#' 
#' @examples 
#' peptideIons_features <- calc_features(all_peptideIons)
#' peptideIons_features_select <- perform_selection(peptideIons_features)
#' 
#' prot_name <- c("1/O75976")
#' test_prot <- test[test$ProteinName==prot_name, ]
#' plot_bar_intensity_n_probability(test_prot)
#' 
plot_bar_intensity_n_probability <- function(case, cutoff_prob=0.3, 
                                             plot_intensity = T, 
                                               plot_prob = T) {
  
  a_reorder <- compute_cor(case = case, cutoff_prob=cutoff_prob)
  
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


#' Plot density distribution of peptide features
#' 
#' @param data data table or data frame in wide representation. The data typically 
#' contains \code{PeptideIon}, \code{ProteinName} and sample names in columns and 
#' measurements of each peptide or precursor ions in rows. The function is 
#' particularly useful after \code{calc_features()} step to visualize distribution of 
#' computed feature statistics. 
#' 
#' @param feature a character vector denoting a peptide feature for density 
#' plots. Examples include \code{"feature_mean_intensity_all"}, 
#' \code{"feature_cv_intensity_all"}, \code{"scaled_median_PCC"} and etc. 
#' @param fill a character vector denoting a column name to color by
#' @param font.size font size of x and y labels 
#' @param title title of density plot 
#'
#' @export
#' 
#' @examples 
#' peptideIons_features <- calc_features(all_peptideIons)
#' plot_density(peptideIons_features, feature = "feature_mean_intensity_all")
#' 
plot_density <- function(data, feature = "feature_mean_intensity_all", 
                         fill = NULL,
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


compute_cor <- function(case, cutoff_prob=0.3){
  
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
  
  return(a_reorder)
}
