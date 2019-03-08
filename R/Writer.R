#' 
#' export protein table to tsv or csv file 
#'
#' @param table 
#' @param name 
#' @param format 
write_protein_table <- function(table, name = "protein_table", format = "tsv"){
  
  if (!format %in% c("tsv","csv")) {
    stop("Please select a valid format to export the table. Options:  \"tsv\", \"csv\"")
  }
  
  file <- paste(name, format, sep = ".")
  if(format == "tsv"){
    write.table(table, file=file, quote=FALSE, sep='\t', col.names = NA)
  }else{
    write.csv(table, file)
  }
}

