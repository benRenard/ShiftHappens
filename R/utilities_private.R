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
#' @return the date d days after origin.date, in the same format as the latter
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

#' Time to numeric format
#'
#' Cast a time or date into a numeric interpreted as the duration(in dayss) after an origin date.
#'
#' @param date vector, time in POSIXct, character or Date format
#'
#' @return List with the following components :
#' \enumerate{
#'   \item origin.date: time in POSIXct, character or Date format, corresponding to the
#'       oldest date from the input d.
#'   \item d: real vector, duration (in days) after origin.date.
#' }
#' @import lubridate
#' @keywords internal
time_to_numeric <- function(date){
  # Examples
  # time_to_numeric(c(as.POSIXct('2024-04-19 12:30:00',tz='UTC'),
  #                   as.POSIXct('2024-03-19 18:30:00',tz='UTC')))
  # time_to_numeric(c(as.Date('2024/02/19'),as.Date('2024/02/10')))
  # time_to_numeric(c(as.character('20240419'),as.character('20240119')))

  if(lubridate::is.Date(date)){
    origin=min(date)
    # Interval in days between origin and each date
    diff_days <- lubridate::time_length(lubridate::interval(origin, date),"day")
    return(list(origin.date=origin,d=diff_days))
  } else if (!lubridate::is.POSIXct(date)){
    date_string <- as.character(date)
    # Attempt to parse the date using various formats
    parsed_date <- lubridate::parse_date_time(date_string,tz = "UTC",
                                              orders  = c('ymd H:M:S', 'ymd', 'mdy',
                                                          'dmy', 'ymd HMS', 'y-m-d H:M:S',
                                                          'y/m/d H:M:S', 'y/m/d HMS' ))
    if(any(is.na(parsed_date)))stop('The format is not supported; please verify the input date or time format')
    # Transform date to numeric format: origin default is "1970-01-01 00:00:00 UTC"
    origin=min(parsed_date)
    diff_days <- lubridate::time_length(lubridate::interval(origin, parsed_date), "day")
    return(list(origin.date=origin,d=diff_days))
  } else {
    # Read time zone
    date.transf <- date
    origin=min(date.transf)
    diff_days <- lubridate::time_length(lubridate::interval(origin, date.transf), "day")
    return(list(origin.date=origin,d=diff_days))
  }
}

