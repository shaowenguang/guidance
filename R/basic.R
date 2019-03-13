#' @importFrom data.table as.data.table
#' @importFrom data.table fread
#' @importFrom data.table dcast
#' @importFrom data.table copy
#' @importFrom MASS lda


#' 
#' @param input_dt a data table, data frame or matrix with numeric values 
#'
#' @return data.table 
#' @export

#replace infinity with NA
replace_inf <- function(input_dt) {
  
  temp <- lapply(input_dt, function(x) {
    if(any(is.infinite(x))) {
      x[is.infinite(x)] <- NA
    }
    return(x)
  })
  
  return(as.data.table(temp))
  
}


#' Count number of unique values 
#' 
#' @param x a vector (atomic or list) or an `expression` object. 
#' 
#' @export
u_count <- function(x) {
  return( length(unique(x)) )
}

#' Count number of values in a vector 
#' 
#' @param x a vector (atomic or list) or an `expression` object. 
#'
#' @export
l_count <- function(x) {
  return( length(which(x)))
}

#' Computes average excluding missing values 
#' 
#' @param x a vector (atomic or list) or an `expression` object. 
#'
#' @export
mean_na <- function(x) {
  
  if(is.nan(mean(x, na.rm=T))) {
    return(NA)
  } else {
    return(mean(x, na.rm=T))
  }
  
}

#' Computes standard deviation excluding missing values 
#' 
#' @param x a vector (atomic or list) or an `expression` object.
#'
#' @export
sd_na <- function(x) {
  
  if(is.nan(sd(x, na.rm=T))) {
    return(NA)
  } else {
    return(sd(x, na.rm=T))
  }

}

#' Computes coefficient of variation (CV) excluding missing values
#' 
#' @param x a vector (atomic or list) or an `expression` object.
#'
#' @export
cv_na <- function(x) {
  return(sd(x, na.rm = T)/mean(x, na.rm = T))
}

#' Computes median excluding missing values 
#' 
#' @param x a vector (atomic or list) or an `expression` object.
#'
#' @export
median_na <- function(x) {
  return(median(x, na.rm = T))
}

#' Computes sum excluding missing values 
#' 
#' @param x a vector (atomic or list) or an `expression` object.
#'
#' @export
sum_na <- function(x) {
  
  if(is.nan(sum(x, na.rm=T))) {
    return(NA)
  } else if( length(which(is.na(x))) == length(x)  ) { 
    #wenguang: basic functions mean and sum return different values for a NA vector. 
    # e.g. > sum(c(NA, NA), na.rm=T)
    #[1] 0
    #> mean(c(NA, NA), na.rm=T)
    #[1] NaN
    return(NA)
  } else {
    return(sum(x, na.rm = T))
  }
  
}

#' Number of intersection between two groups
#' 
#' @param input_a vectors, data frames, or ps objects
#'  containing a sequence of elements (conceptually).
#' @param input_b vectors, data frames, or ps objects
#'  containing a sequence of elements (conceptually).
#'
#' @export
pairwise_length_vector <- function(input_a, input_b) {
  
  return( length(intersect(which(!is.na(input_a)),  which(!is.na(input_b)) )) )
  
}


#' @export
pairwise_length_matrix <- function(input_matrix) {
  
  output_matrix <- matrix(0, dim(input_matrix)[1], dim(input_matrix)[1])
  
  for(i in 1:dim(output_matrix)[1]) {
    #wenguang: probably wrong... needs to be checked
    cat("probably wrong... needs to be checked", "\n")
    output_matrix[,i] <- apply(case[, 3:17, with=F], 1, function(x) pairwise_length_vector(x, case[i, 3:17, with=F]))
  }
  
  return(output_matrix)
}


#' Number of intersection betwen numeric groups 
#' 
#' @param input_vector_a vectors containing a sequence of numbers
#' @param input_vector_b vectors containing a sequence of numbers  
#'
#' @export
count_pairwise_number <- function(input_vector_a, input_vector_b) {
  
  return( length(intersect( which(input_vector_a > -99999999999 & input_vector_a < 99999999999), 
                            which(input_vector_b > -99999999999 & input_vector_b < 99999999999) )) )
  
}


#' @param input_matrix 
#'
#' @export
count_pairwise_number_matrix <- function(input_matrix) {
  
  out_m <- matrix(0, dim(input_matrix)[2], dim(input_matrix)[2])
  
  for(i in 1:dim(input_matrix)[2]) {
    out_m[i, ] <- apply(input_matrix, 2, function(x) count_pairwise_number(x, input_matrix[, i]))
  }
  
  return(out_m)
}
