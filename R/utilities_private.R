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
