#' @importFrom data.table as.data.table
#' @importFrom data.table fread
#' @importFrom data.table dcast
#' @importFrom data.table copy
#' @importFrom MASS lda


#' 
#' @description A function to replace infinity values by NA in a data table 
#' @param input_dt a data table, data frame or matrix with numeric values 
#' 
#' @return data.table 
#' @export
#' 
#' @examples 
#'
# Replace infinity values by NA
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
#' @description A function to count the number of unique values in a numeric 
#' or character vector
#' 
#' @param x a vector (atomic or list) or an `expression` object. 
#' 
#' @return a numeric value
#' @export 
#' 
#' @examples 
#' vector <- c("c", "d", "v", "d", "t", "c", "a")
#' u_count(vector)
u_count <- function(x) {
  return( length(unique(x)) )
}

#' Count TRUE values 
#' 
#' @description A function to count the number of TRUE values in a vector 
#' @param x a vector (atomic or list) in TRUE or FALSE 
#'
#' @export
#' @return a numeric value 
#' 
#' @examples
#' vector <- c(TRUE, TRUE, FALSE, TRUE)
#' l_count(vector)
#' 
l_count <- function(x) {
  return( length(which(x)))
}

#' Computes mean  
#' 
#' @description A function to compute average of a numeric vector 
#' after exclusion of missing values 
#' @param x a vector (atomic or list) or an `expression` object. 
#' 
#' @return a numeric value
#' @export
#' 
#' @examples 
#' vector <- c(1, 2, 4, NA, 2, 5, 6)
#' mean_na(vector)
mean_na <- function(x) {
  
  if(is.nan(mean(x, na.rm=T))) {
    return(NA)
  } else {
    return(mean(x, na.rm=T))
  }
  
}

#' Computes standard deviation 
#' 
#' @description A function to compute standard deviation of a numeric vector 
#' after exclusion of missing values 
#' 
#' @param x a vector (atomic or list) or an `expression` object.
#'
#' @return a numeric value
#' @export
#' 
#' @examples 
#' vector <- c(1, 2, 4, NA, 2, 5, 6)
#' sd_na(vector)
sd_na <- function(x) {
  
  if(is.nan(sd(x, na.rm=T))) {
    return(NA)
  } else {
    return(sd(x, na.rm=T))
  }

}

#' Computes coefficient of variation (CV) 
#' 
#' @description A function to compute coefficient of variance (CV) of a numeric vector 
#' after exclusion of missing values 
#' @param x a vector (atomic or list) or an `expression` object.
#'
#' @export
#' @return a numeric value
#' 
#' @examples 
#' vector <- c(1, 2, 4, NA, 2, 5, 6)
#' cv_na(vector)
cv_na <- function(x) {
  return(sd(x, na.rm = T)/mean(x, na.rm = T))
}

#' Computes median  
#' 
#' @description A function to compute median of a numeric vector 
#' after exclusion of missing values 
#' 
#' @param x a vector (atomic or list) or an `expression` object.
#'
#' @return a numeric value
#' @export
#' 
#' @examples 
#' vector <- c(1, 2, 4, NA, 2, 5, 6)
#' median_na(vector)
median_na <- function(x) {
  return(median(x, na.rm = T))
}

#' Computes sum 
#' 
#' @description A function to compute sum of a numeric vector 
#' after exclusion of missing values 
#' @param x a vector (atomic or list) or an `expression` object.
#'
#' @return a numeric value
#' @export
#' 
#' @examples 
#' vector <- c(1, 2, 4, NA, 2, 5, 6)
#' sum_na(vector)
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

#' Number of intersection
#' 
#' @description A function to calculate the number of intersecting 
#' values between two vectors, excluding missing values 
#' @param input_a vectors, data frames, or ps objects
#'  containing a sequence of elements (conceptually).
#' @param input_b vectors, data frames, or ps objects
#'  containing a sequence of elements (conceptually).
#'
#' @export
#' @examples 
#' vector1 <- c(1, 2, 4, NA, 2, 5, 6)
#' vector2 <- c(1, 2, 4, NA, 6, 7, 8)
#' pairwise_length_vector(vector1, vector2)
#' 
pairwise_length_vector <- function(input_a, input_b) {
  
  return( length(intersect(which(!is.na(input_a)),  which(!is.na(input_b)) )) )
  
}


#' Number of intersection within a matrix 
#'
#' @description A function to calculate the number of intersecting values 
#' between rows of input matrix 
#' @param input_matrix a matrix containing numeric values 
#'
#' @return
#' @export
#'
#' @examples
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
#' @description A function to calculate the number of intersecting values 
#' between two numeric vectors. 
#' 
#' @param input_vector_a vectors containing a sequence of numeric values 
#' @param input_vector_b vectors containing a sequence of numeric values
#' @export
#' @examples 
#' vector1 <- c(1, 2, 5, 100, 200, 50000)
#' vector2 <- c(1, 2, 4, 100, 500, 40000)
#' count_pairwise_number(vector1, vector2)
count_pairwise_number <- function(input_vector_a, input_vector_b) {
  
  return( length(intersect( which(input_vector_a > -99999999999 & input_vector_a < 99999999999), 
                            which(input_vector_b > -99999999999 & input_vector_b < 99999999999) )) )
  
}


#' Number of intersection between samples within a matrix 
#' 
#' @description A function to calculate the number of intersecting numeric values 
#' between columns of a matrix.
#' 
#' @param input_matrix matrix containing a numeric values 
#'
#' @export
count_pairwise_number_matrix <- function(input_matrix) {
  
  out_m <- matrix(0, dim(input_matrix)[2], dim(input_matrix)[2])
  
  for(i in 1:dim(input_matrix)[2]) {
    out_m[i, ] <- apply(input_matrix, 2, function(x) count_pairwise_number(x, input_matrix[, i]))
  }
  
  return(out_m)
}

