#' Check vectors length equality
#'
#' Check whether all vectors have the same length
#' @param ... real vectors
#' @return logical, return null if the vectors have not the same length
#' @keywords internal
#' @source \url{https://www.simonqueenborough.info/R/basic/lessons/lapply_and_sapply.html}
check_equal_length <- function(...) {
  lengths <- sapply(list(...), length)
  if (length(unique(lengths)) != 1){
    return(NULL)
  }else{
    return('ok')
  }
}

#' Check square matrix
#'
#' Check whether a matrix is square
#'
#' @param x matrix
#' @return logical, null if matrix is not square
#' @keywords internal
check_square_matrix <- function(x){
  if(ncol(x)==nrow(x)){
    return('ok')
  }else{
    return(NULL)
  }
}

#' Check integer
#'
#' Check is a number is an integer (in the mathematical sense, not the R one)
#'
#' @param x vector
#' @return logical, TRUE if all values in x are integer, FALSE otherwise
#' @keywords internal
is_integer <- function(x){
  out=all(x==as.integer(x))
  return(out)
}
