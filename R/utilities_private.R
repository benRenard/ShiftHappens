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

#' Numeric to time format
#'
#' Cast a numeric interpreted as the duration after an origin date into a time format.
#'
#' @param d real vector, duration (in days) after origin.date
#' @param origin.date value, origin date in character, POSIXct or Date format
#' @return the date nday days after origin.date, in the same format as the latter
#' @import lubridate
#' @keywords internal
numeric_to_time <- function(d,origin.date){
  if (lubridate::is.Date(origin.date)) {
    return(origin.date+d)
  }else if(lubridate::is.POSIXct(origin.date)){
    return(origin.date+d*86400)
  }else{ # character format
    # Attempt to parse the date using various formats
    parsed_date <- lubridate::parse_date_time(origin.date, tz = "UTC",
                                              orders  = c('ymd H:M:S', 'ymd', 'mdy',
                                                          'dmy', 'ymd HMS', 'y-m-d H:M:S',
                                                          'y/m/d H:M:S', 'y/m/d HMS' ))
    return(parsed_date+d*86400)
  }
}
