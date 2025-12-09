#***************************************************************************----
# Constructors ----

#' simpleSegmentation object constructor.
#'
#' Creates a new instance of a 'simpleSegmentation' object, used to
#' store the results of a segmentation with an known number of segment
#' @param obs real vector, observations
#' @param time real vector, time
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param data data frame, data summary, should contain the following columns:
#'     time,obs,u,I95_lower,I95_upper,period.
#'     If NULL, will be automatically intialised from obs,time,u.
#'     If provided, will supersede obs, time, u.
#' @param shifts data frame, detected shifts (with columns tau, I95_lower and I95_upper)
#' @param mcmc data frame, MCMC simulations
#' @param DIC real, DIC estimation
#' @param origin.date dates origin (useful for date conversions)
#' @return An object of class [simpleSegmentation()], containing the following fields:
#' \enumerate{
#'   \item data: data frame, all data with their respective periods after segmentation
#'   \item shifts: data frame, all detected shift times in numeric or POSIXct format in UTC
#'   \item mcmc: data frame, MCMC simulations
#'   \item DIC: real, DIC estimation
#'   \item origin.date: positive real or date, date describing origin of the segmentation for a sample. Useful for recursive segmentation.
#' }
#' @examples
#' sg <- simpleSegmentation(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH)
#' @export
simpleSegmentation<-function(obs,time=1:length(obs),u=0*obs,data=NULL,
                             shifts=data.frame(tau=numeric(0),I95_lower=numeric(0),I95_upper=numeric(0)),
                             mcmc=data.frame(),DIC=NA,origin.date=min(time)){
  o<-new_simpleSegmentation(obs,time,u,data,shifts,mcmc,DIC,origin.date)
  return(validate_simpleSegmentation(o))
}

#' multipleSegmentation object constructor.
#'
#' Creates a new instance of a 'multipleSegmentation' object, used to
#' store the results of a segmentation with an unknown number of segment.
#' @param results vector of 'simpleSegmentation' objects, segmentation results for each tested number of segments
#'
#' @return An object of class [multipleSegmentation()], containing the following fields:
#' \enumerate{
#'   \item nS: integer, optimal number of segments (minimum DIC)
#'   \item DICs: real vector, DICs computed for each number of segment
#'   \item data: data frame, all data with their respective periods after segmentation (for optimal nS)
#'   \item shifts: data frame, all detected shift time in numeric or POSIXct format in UTC (for optimal nS)
#'   \item mcmc: data frame, MCMC simulations (for optimal nS)
#'   \item DIC: real, DIC estimation (for optimal nS)
#'   \item origin.date: positive real or date, date describing origin of the segmentation for a sample. Useful for recursive segmentation.
#'   \item results: list, intermediate results for all tested number of segments see ?Segmentation_Engine
#' }
#' @examples
#' results=list()
#' results[[1]]=Segmentation_Engine(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,
#'                                  u=RhoneRiverAMAX$uH,nS=1)
#' results[[2]]=Segmentation_Engine(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,
#'                                  u=RhoneRiverAMAX$uH,nS=2)
#' sg <- multipleSegmentation(results)
#' @export
multipleSegmentation<-function(results){
  o<-new_multipleSegmentation(results)
  return(validate_multipleSegmentation(o))
}

#***************************************************************************----
# is functions ----

#' simpleSegmentation tester
#'
#' Is an object of class 'simpleSegmentation'?
#'
#' @param o Object, an object.
#' @return A logical equal to TRUE if class(o)== 'simpleSegmentation', FALSE otherwise.
#' @keywords internal
is.simpleSegmentation<-function(o){
  return(class(o)=='simpleSegmentation')
}

#' multipleSegmentation tester
#'
#' Is an object of class 'multipleSegmentation'?
#'
#' @param o Object, an object.
#' @return A logical equal to TRUE if class(o)== 'multipleSegmentation', FALSE otherwise.
#' @keywords internal
is.multipleSegmentation<-function(o){
  return(class(o)=='multipleSegmentation')
}

#***************************************************************************----
# internal constructors ----

new_simpleSegmentation<-function(obs,time,u,data,shifts,mcmc,DIC,origin.date){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(is.numeric(obs))
  stopifnot(is.vector(obs))
  stopifnot(is.numeric(time))
  stopifnot(is.vector(time))
  stopifnot(is.numeric(u))
  stopifnot(is.vector(u))
  if(!is.null(data)){stopifnot(is.data.frame(data))}
  stopifnot(is.data.frame(shifts))
  stopifnot(is.data.frame(mcmc))
  stopifnot(is.na(DIC) | is.numeric(DIC))
  stopifnot(is.numeric(origin.date))
  if(is.null(check_equal_length(obs,time,u))){
    stop('The observations, time and uncertainty do not have the same length')
  }
  # assemble object
  o=list()
  if(is.null(data)){
    o$data=data.frame(time=time,obs=obs,u=u,
                      I95_lower=obs-1.96*u,I95_upper=obs+1.96*u,period=1)
  } else {
    o$data=data
  }
  o$shifts=shifts
  o$mcmc=mcmc
  o$DIC=DIC
  o$origin.date=origin.date
  class(o) <- 'simpleSegmentation'
  return(o)
}

new_multipleSegmentation<-function(results){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(is.list(results))
  if(length(results)<1){
    stop('Results list should be non-empty')
  }
  stopifnot(all(sapply(results,is.simpleSegmentation)))
  # assemble object
  DICs=sapply(results,function(x){x$DIC})
  nS=which.min(DICs)
  o=results[[nS]] # results for optimal number of segments
  o$nS=nS # Add optimal number of segments
  o$DICs=DICs # Add DICs for each number of segments
  o$results=results # Add all partial results for all tested number of segments
  class(o) <- 'multipleSegmentation'
  return(o)
}

#***************************************************************************----
# validators ----

validate_simpleSegmentation<-function(x){
  # nothing to do
  return(x)
}

validate_multipleSegmentation<-function(x){
  # nothing to do
  return(x)
}

