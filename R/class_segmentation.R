#***************************************************************************----
# Constructors ----

#' simpleSegmentation object constructor.
#'
#' Creates a new instance of a 'simpleSegmentation' object
#' @param obs real vector, observations
#' @param time real vector, time
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param mcmc data frame, MCMC simulations
#' @param DIC real, DIC estimation
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
simpleSegmentation<-function(obs,time=1:length(obs),u=0*obs,
                             mcmc=data.frame(),DIC=NA){
  o<-new_simpleSegmentation(obs,time,u,mcmc,DIC)
  return(validate_simpleSegmentation(o))
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

#***************************************************************************----
# internal constructors ----

new_simpleSegmentation<-function(obs,time,u,mcmc,DIC){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(is.numeric(obs))
  stopifnot(is.vector(obs))
  stopifnot(is.numeric(time))
  stopifnot(is.vector(time))
  stopifnot(is.numeric(u))
  stopifnot(is.vector(u))
  stopifnot(is.data.frame(mcmc))
  stopifnot(is.na(DIC) | is.numeric(DIC))
  if(is.null(check_equal_length(obs,time,u))){
    stop('The observations, time and uncertainty do not have the same length')
  }
  # assemble object
  o=list()
  o$data=data.frame(time=time,obs=obs,u=u,
                      I95_lower=obs-1.96*u,I95_upper=obs+1.96*u,period=1)
  o$shifts=data.frame(tau=numeric(0),I95_lower=numeric(0),I95_upper=numeric(0))
  o$mcmc=mcmc
  o$DIC=DIC
  o$origin.date=min(time)
  class(o) <- 'simpleSegmentation'
  return(o)
}

#***************************************************************************----
# validators ----

validate_simpleSegmentation<-function(x){
  # nothing to do
  return(x)
}

