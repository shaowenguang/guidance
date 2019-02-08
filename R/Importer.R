#' @importFrom data.table fread data.table


#' @export
import_openswath <- function(search_results, bool.removeDecoy = T, level = "PeptideIon", sample_annotation = NULL, remove_prefixInFileName = FALSE) {
  
  if (!level %in% c("PeptideIon","Transition", "Peptide", "PeptideWithMod")) {
    stop("Please select a valid type for imported data to be kept in Prom. Options:  \"Transition\", \"PeptideIon(default)\", \"PeptideWithMod\", \"Peptide\"")
  }
  
  global_level <<- level 
  
  message("start to read raw OpenSWATH search results...")
  
  raw <- fread(input=search_results, head=T)

  if(bool.removeDecoy==T) {
    raw <- subset(raw, raw$decoy==0)
  }
  
  names(raw)[1] <- "PeptideIon"
  names(raw)[which(names(raw)=="Sequence")] <- "Peptide"
  names(raw)[which(grepl("Full", names(raw)) & grepl("PeptideName", names(raw)))] <- "PeptideWithMod" #wenguang: to be compatible with both "FullUniModPeptideName" and "FullPeptideName"
  names(raw)[which(names(raw)=="m_score")] <- "Score"
  
  raw$PeptideIon <- paste0(raw$PeptideWithMod, "_", raw$Charge)
  raw[, Width := rightWidth - leftWidth]

  #wenguang: in the search results from euler portal, the filename was in the format as "/scratch/71239421.tmpdir/xuep_J180621_SW_3.mzXML.gz". this will remove this unnecessary prefix...
  if(remove_prefixInFileName == TRUE) {
    raw$filename <- sapply(strsplit(raw$filename, "/"), "[[", 4)
  }
  
  summarize_data(raw)
  
#  raw$PeptideIon <- paste(sapply(strsplit(raw$PeptideIon, "_"), "[[", 2), "_", sapply(strsplit(raw$PeptideIon, "_"), "[[", 3), sep="")
  
#  raw$filename <- sapply(strsplit(raw$filename, "/"), "[[", 4)
  
  if(global_level=="Transition") {
    
    numTrans <- length(unlist(strsplit(as.character(raw[1,]$aggr_Peak_Area), ";")))
    raw_expanded <- data.table( "PeptideIon"=rep(raw$PeptideIon, each=numTrans), 
                        "ProteinName"=rep(raw$ProteinName, each=numTrans), 
                        "filename"=rep(raw$filename, each=numTrans), 
                        "Intensity"=rep(raw$Intensity, each=numTrans), 
                        "RT"=rep(raw$RT, each=numTrans), 
                        "Score"=rep(raw$Score, each=numTrans), 
                        "aggr_Fragment_Annotation"=unlist(strsplit(as.character(raw$aggr_Fragment_Annotation), ";")), 
                        "aggr_Peak_Area"=as.numeric(as.character(unlist(strsplit(as.character(raw$aggr_Peak_Area), ";")))) )
    
    raw_expanded[raw_expanded$aggr_Peak_Area == 0 ,]$aggr_Peak_Area <- NA
    raw <- raw_expanded

  } else {
    
    keep.columns <- c(global_level, "ProteinName", "filename", "Intensity", "RT", "Score", "Width")
    raw <- subset(raw, select=keep.columns)
    
  }
  
  
  if(is.null(sample_annotation)) {
    
    message("no sample annotation table was provided. each injection will be treated independently")
    
    return(raw)
    
  } else {
    
    read_sample_annotation(input_file = sample_annotation)
    annotated <- annotate_sample(raw, anno)
    
    return(annotated)
    
  }
  
} 





#' @export
import_openswath_matrix_fromEulerPortal <- function(search_results, sample_annotation = NULL) {
  
  message("Message: It starts reading the search result matrix of OpenSWATH generated from Euler Portal...")
  
  global_level <<- "PeptideIon"
  
  raw <- fread(input=search_results, head=T)
  
  names(raw)[1:2] <- c("PeptideIon", "ProteinName")
  
  if(is.null(sample_annotation)) {
    
    message("no sample annotation table was provided. each injection will be treated independently")
    
  } else {
    
    read_sample_annotation(input_file = sample_annotation)

  }

  #wenguang: change "_Intensity" from suffix to prefix...
  names(raw)[grepl("Intensity$", names(raw))] <- paste0("Intensity_", gsub("_Intensity", "", names(raw)[grepl("Intensity$", names(raw))]))
  names(raw)[grepl("Score$", names(raw))] <- paste0("Score_", gsub("_Score", "", names(raw)[grepl("Score$", names(raw))]))
  names(raw)[grepl("RT$", names(raw))] <- paste0("RT_", gsub("_RT", "", names(raw)[grepl("RT$", names(raw))]))
  names(raw)[grepl("Width$", names(raw))] <- paste0("Width_", gsub("_Width", "", names(raw)[grepl("Width$", names(raw))]))
  
  return(raw)
    
}








#' @export
read_sample_annotation <- function(input_file="sample_annotation") {

#  sample_annotation <- as.data.frame(read.table(file=input_file, fill=T, header=T, stringsAsFactors=F))
  
#  anno <<- sample_annotation
  
#  return(sample_annotation)

  anno <<- as.data.frame(read.table(file=input_file, fill=T, header=T, stringsAsFactors=F))
    
}



#' @export
annotate_sample <- function(search_result, sample_annotation) {
  
  if( !all(unique(search_result$filename) %in% unique(sample_annotation$filename)) ) {
    message("******WARNING******: Some samples were not annotated in the annotation table.")
    message("******WARNING******: Those unannotated samples would be exclued in the further analysis.\n")
  }
  
  if( any(duplicated(sample_annotation$filename)) ) {
    stop("Error: some samples were duplicately annotated in the annotation table.")
  }

  merged <- merge(search_result, sample_annotation, all.y=T, by.x="filename", by.y="filename")
  
  if(global_level=="Transition") {
    merged <- subset(merged, select=c("aggr_Fragment_Annotation", "PeptideIon", "ProteinName", "aggr_Peak_Area", 
                                      "Intensity", "InjectionName", "SampleName", "RT", "Score"))
  } else {
    merged <- subset(merged, select=c(global_level, "ProteinName", "InjectionName", "SampleName", "Intensity", "RT", "Score", "Width"))
  }
  
  return(merged)
  
}


#' @export
remove_prefix <- function(input_dt) {
  
  input_dt$PeptideIon <- paste(sapply(strsplit(input_dt$PeptideIon, "_"), "[[", 2), "_", sapply(strsplit(input_dt$PeptideIon, "_"), "[[", 3), sep="")
  
  if(length(which(duplicated(input_dt$PeptideIon)))>0) {
    
    message("******WARNING******: The assay library contains multiple transition groups (",  length(which(duplicated(input_dt$PeptideIon))), ") with the same precursor. It may be due to Ludo phospho project or N-term Caramylation (not sure...).")
    message("******WARNING******: Now ALL these transitions are REMOVED!!!\n")
    
    n_occur <- data.frame(table(input_dt$PeptideIon))
    
    input_dt <- input_dt[input_dt$PeptideIon %in% n_occur$Var1[n_occur$Freq == 1], ]
    
  }
  
  return(input_dt)
  
}
