#***************************************************************************----
# Constructor ----

#' recursiveModeling object constructor.
#'
#' Creates a new instance of a 'recursiveModeling' object, used to
#' store the results of a recursive modeling segmentation
#' @param segmentation [recursiveSegmentation()] object, segmentation results
#' @param fit list of [fittedModel()] objects, fitted model results
#' @return An object of class [recursiveModeling()], containing the following fields:
#' \enumerate{
#'   \item segmentation: an object of class [recursiveSegmentation()],
#'       containing the results of the segmentation of residuals
#'   \item fit: a list of objects of class [fittedModel()],
#'       containing the results of the model fit at each stage of the recursion
#' }
#' @examples
#' rec <- recursiveModeling()
#' @export
recursiveModeling<-function(segmentation=recursiveSegmentation(),fit=list(fittedModel())){
  o<-new_recursiveModeling(segmentation,fit)
  return(validate_recursiveModeling(o))
}

#***************************************************************************----
# is function ----

#' recursiveModeling object tester
#'
#' Is an object of class 'recursiveModeling'?
#'
#' @param o Object, an object.
#' @return A logical equal to TRUE if class(o)== 'recursiveModeling', FALSE otherwise.
#' @keywords internal
is.recursiveModeling<-function(o){
  return(class(o)=='recursiveModeling')
}

#***************************************************************************----
# internal constructor ----

new_recursiveModeling<-function(segmentation,fit){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(inherits(segmentation,'recursiveSegmentation'))
  stopifnot(is.list(fit))

  for (i in 1:length(fit)){
    stopifnot(inherits(fit[[i]],'fittedModel'))
  }
  # assemble object
  o=list(segmentation=segmentation,fit=fit)
  class(o) <- 'recursiveModeling'
  return(o)
}

#***************************************************************************----
# validator ----

validate_recursiveModeling<-function(x){
  # nothing to do
  return(x)
}
