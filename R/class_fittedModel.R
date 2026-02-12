#***************************************************************************----
# Constructor ----

#' fittedModel object constructor.
#'
#' Creates a new instance of a 'fittedModel' object, used to
#' store the results of a fitted model
#' @param data data.frame, observations, simulations and residuals.
#'     Should at least contain the following columns, which will be used
#'     for segmentation: 'time','res' and 'uRes'.
#' @param parameters data.frame, fitted parameters
#' @return An object of class [fittedModel()], containing the following fields:
#' \enumerate{
#'   \item data: data frame, observations, simulations and residuals
#'   \item parameters: data frame, fitted parameters
#' }
#' @examples
#' f <- fittedModel()
#' @export
fittedModel<-function(data=data.frame(time=numeric(0),res=numeric(0),uRes=numeric(0)),
                      parameters=numeric(0)){
  o<-new_fittedModel(data,parameters)
  return(validate_fittedModel(o))
}

#***************************************************************************----
# is function ----

#' fittedModel object tester
#'
#' Is an object of class 'fittedModel'?
#'
#' @param o Object, an object.
#' @return A logical equal to TRUE if class(o)== 'fittedModel', FALSE otherwise.
#' @keywords internal
is.fittedModel<-function(o){
  return(class(o)=='fittedModel')
}

#***************************************************************************----
# internal constructor ----

new_fittedModel<-function(data,parameters){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(is.data.frame(data))
  stopifnot(is.numeric(parameters))
  # assemble object
  o=list(data=data,parameters=parameters)
  class(o) <- 'fittedModel'
  return(o)
}

#***************************************************************************----
# validator ----

validate_fittedModel<-function(x){
  noms=names(x$data)
  if( ! ( ('time' %in% noms) & ('res' %in% noms) & ('uRes' %in% noms) ) ){
    mess=paste0("Invalid data: the data frame should at least contain columns named time, res and uRes")
    stop(mess,call.=FALSE)
  }
  return(x)
}
