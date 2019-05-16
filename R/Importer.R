#' Import openSWATH data table 
#' 
#' @description A function to import openSWATH data and filter for relevant columns. 
#' The imported data table contains columns denoting name of peptide ion, 
#' corresponding protein name, injection name, sample name, intensity values, 
#' retention time, score and width of spectral peaks. In particular for transition-level 
#' data, the table contains aggregate fragment annotation and aggregated peak area 
#' in place of peak width. 
#' 
#' @param search_results A data frame containing the SWATH-MS data. This data typically
#' contains peptide precursors in each row with corresponding \code{ProteinName}, 
#' \code{Intensity}, \code{RT}, \code{Score} and etc in columns. The data can be loaded
#' from local directory.
#' @param bool.removeDecoy boolean value (\code{TRUE} or \code{FALSE}) determining 
#' if the decoy peptides should be removed from the imported data 
#' @param level the protein-level of imported data. The options include 
#' \code{"PeptideIon"} (by default), \code{"Transition"}, \code{"Peptide"} or 
#' \code{"PeptideWithMod"}
#' @param sample_annotation data matrix with \code{SampleName}, biological covariates 
#' (biological replicates) and technical covariates (technical replicates, batches, etc)
#' @param remove_prefixInFileName boolean value (\code{TRUE} or \code{FALSE}) whether to 
#' remove unnecessary prefix from euler portal file name. For example, a typical 
#' filename will look like \code{"/scratch/71239421.tmpdir/xuep_J180621_SW_3.mzXML.gz"} 
#' and \code{remove_prefixInFileName = TRUE} will result in \code{"xuep_J180621_SW_3.mzXML.gz"}
#' 
#' @return data.table data.frame containing SWATH-MS data and sample annotation
#' 
#' @export
#' 
#' @examples peptideIons <- import_openswath(search_results= "data/QGS_SWATH_data", 
#' sample_annotation="data/QGS_sample_annotation", 
#' level="PeptideIon")
#' 
import_openswath <- function(search_results, bool.removeDecoy = T, level = "PeptideIon", 
                             sample_annotation = NULL, remove_prefixInFileName = FALSE) {
  
  if (!level %in% c("PeptideIon","Transition", "Peptide", "PeptideWithMod")) {
    stop("Please select a valid type for imported data to be kept in Prom. Options:  \"Transition\", \"PeptideIon(default)\", \"PeptideWithMod\", \"Peptide\"")
  }
  
  global_level <<- level 
  
  message("It starts to read raw OpenSWATH search results...")
  
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
    
    message("  no sample annotation table was provided. each injection will be treated independently")
    message("Done with importing search results...")

    return(raw)
    
  } else {
    
    read_sample_annotation(input_file = sample_annotation)
    annotated <- annotate_sample(raw, anno)
    
    message("Done with importing search results...")

    return(annotated)
    
  }
  
} 


#' Import openSWATH data table from Eular portal 
#' 
#' @description A function to import a full openSWATH data table, 
#' typically pre-processed or pre-filtered for downstream analysis. 
#' The imported data table contains columns denoting peptide name, corresponding 
#' protein name, intensity, retention time and score of each peptide. 
#' 
#' @param search_results  A data frame containing the SWATH-MS data. This data typically
#' contains peptide ions in each row with corresponding \code{ProteinName}, 
#' \code{Intensity}, \code{RT}, \code{Score} and etc in columns. The data can be loaded
#' from local directory.
#'
#' @param sample_annotation data matrix with \code{SampleName}, biological covariates 
#' (biological replicates) and technical covariates (technical replicates, batches, etc)
#'
#' @return data.table data containing raw openSWATH matrix 
#' 
#' @export
#' 
#' @examples 
#' all_peptideIons <- import_openswath_matrix_fromEulerPortal(
#' search_results="data/YS_SWATH_data", 
#' sample_annotation="data/YS_sample_annotation")
#' 
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

#' Import Spectronaut data table
#' 
#' @description A function to import a full Spectronaut data table, typically 
#' pre-processed or pre-filtered for downstream analysis. The imported data 
#' table contains columns denoting peptide name, corresponding protein name, 
#' intensity, retention time and score of each peptide. From the Spectronaut 
#' output, this function takes PEP Quantity value as Intensity and EG q-value 
#' as Score. 
#' 
#' @param search_results A data frame containing the SWATH-MS data. This data typically
#' contains peptide ions in each row with corresponding \code{ProteinName}, 
#' \code{Intensity}, \code{RT}, \code{Score} and etc in columns. The data can be loaded
#' from local directory.
#'
#' @param sample_annotation data matrix with \code{SampleName}, biological covariates 
#' (biological replicates) and technical covariates (technical replicates, batches, etc)
#'
#' @export 
#' 
#' @examples 
#' peptideIons <- import_spectronaut_matrix(search_results= "data/QGS_SWATH_data", 
#' sample_annotation="data/QGS_sample_annotation")
#' 
import_spectronaut_matrix <- function(search_results, sample_annotation = NULL) {

  message("Message: It starts reading the search result matrix of Spectronaut...")
  
  global_level <<- "PeptideIon"
  
  raw <- fread(input=search_results, head=T)
  
  names(raw)[1:2] <- c("PeptideIon", "ProteinName")
  
  if(is.null(sample_annotation)) {
    
    message("no sample annotation table was provided. each injection will be treated independently")
    
  } else {
    
    read_sample_annotation(input_file = sample_annotation)

  }

# modify the column name accordingly...
if(length(which(grepl(".PEP.Quantity$", names(raw)))) > 0) {
  names(raw)[grepl(".PEP.Quantity$", names(raw))] <- paste0("Intensity_", gsub(".PEP.Quantity", "", names(raw)[grepl(".PEP.Quantity$", names(raw))]))
} 

if(length(which(grepl(".EG.Qvalue$", names(raw)))) > 0) {
  names(raw)[grepl(".EG.Qvalue$", names(raw))] <- paste0("Score_", gsub(".EG.Qvalue", "", names(raw)[grepl(".EG.Qvalue$", names(raw))]))
} 

  #wenguang: change "_Intensity" from suffix to prefix...
  names(raw)[grepl("Intensity$", names(raw))] <- paste0("Intensity_", gsub("_Intensity", "", names(raw)[grepl("Intensity$", names(raw))]))
  names(raw)[grepl("Score$", names(raw))] <- paste0("Score_", gsub("_Score", "", names(raw)[grepl("Score$", names(raw))]))
  names(raw)[grepl("RT$", names(raw))] <- paste0("RT_", gsub("_RT", "", names(raw)[grepl("RT$", names(raw))]))
  names(raw)[grepl("Width$", names(raw))] <- paste0("Width_", gsub("_Width", "", names(raw)[grepl("Width$", names(raw))]))
  
  return(raw)


}


read_sample_annotation <- function(input_file="sample_annotation") {

#  sample_annotation <- as.data.frame(read.table(file=input_file, fill=T, header=T, stringsAsFactors=F))
  
#  anno <<- sample_annotation
  
#  return(sample_annotation)

##  anno <<- as.data.frame(read.table(file=input_file, fill=T, header=T, stringsAsFactors=F))

  anno <<- as.data.frame(read.table(file = input_file, fill = T, header = T, stringsAsFactors = F))

    
}


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
