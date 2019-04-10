#' Export data table 
#' 
#' @description Export data table e.g. protein table to a tsv or csv file. 
#'
#' @param input_dt the object to be written, preferably data frame, data table or matrix
#' @param file_name name of exported file 
#' @param format format the object to be exported. Options include \code{"tsv"} and 
#' \code{"csv"}. 
#' @param ... other parameters of \code{write.table} function
#'  
#' @return csv or tsv output
#' @export 
#' 
#' @examples \dontrun{write_protein_table(protein_Filtered_top3_sum_ImputedWeighted, 
#' file_name = "protein_Filtered_top3_sum_ImputedWeighted", format = "tsv")}
#' 
write_protein_table <- function(input_dt, file_name = "protein_table", format = "tsv", ...){
  
  if (!format %in% c("tsv","csv")) {
    stop("Please select a valid format to export the table. Options:  \"tsv\", \"csv\"")
  }
  
  file <- paste(file_name, format, sep = ".")
  if(format == "tsv"){
    write.table(input_dt, file=file, quote=FALSE, sep='\t', col.names = NA, ...)
  }else{
    write.csv(input_dt, file, ...)
  }
}

