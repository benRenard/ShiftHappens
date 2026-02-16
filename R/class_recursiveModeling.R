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

#***************************************************************************----
# plotting ----

plot_recursiveModeling_data <- function(x,type=c('xy','ty','tx')){
  typ=match.arg(type)
  if(typ=='xy'){
    xlab='x';ylab='y'
  } else if (typ=='ty'){
    xlab='Time';ylab='y'
  } else if (typ=='tx'){
    xlab='Time';ylab='x'
  } else {
    stop('Unknown plot type')
  }
  terminal=x$segmentation$tree$indx[x$segmentation$tree$nS==1]
  DF=c()
  for(i in 1:length(terminal)){
    node=terminal[i]
    foo=x$fit[[node]]$data
    nX=NCOL(foo)-5
    if(nX==1){
      names(foo)[2]='x'
      DF=rbind(DF,cbind(foo,period=i))
    } else if(nX >1) {

    } else {

    }
  }
  DF$period=as.factor(DF$period)
  # 2DO: review fitModel object to store original data in a standardized way, including uY/uX,
  # hence allowing plotting error bars
  # 2DO: handle multi-X case
  DF$uY=0
  colourCount_obs=length(unique(DF$period))
  getPalette_obs=scales::viridis_pal(option='D')
  if(typ=='xy'){
    g=ggplot(DF,aes(x=.data$x))
  } else {
    g=ggplot(DF,aes(x=.data$time))
  }
  if(typ=='tx'){
    g=g+geom_point(aes(y=.data$x,group=.data$period,color=.data$period))
  } else {
    g=g+geom_point(aes(y=.data$obs,group=.data$period,color=.data$period))
    if(any(DF$uY>0)){
      g=g+geom_errorbar(aes(y=.data$obs,
                            ymin=.data$obs-1.96*.data$uY,
                            ymax=.data$obs+1.96*.data$uY,
                            color=.data$period,group=.data$period))
    }
  }
  if(typ=='xy'){
    g=g+geom_line(aes(y=.data$sim,group=.data$period,color=.data$period))
  }
  g=g+scale_color_manual(values=getPalette_obs(colourCount_obs))
  g=g+labs(y=ylab,x=xlab,title='Segmented data')
  g=g+theme_bw()+theme(plot.title=element_text(hjust=0.5,face='bold',size=15),
                       legend.position="none")
  return(g)
}

#' recursiveModeling plotter
#'
#' Plot a recursiveModeling object.
#' @param x object of class recursiveModeling
#' @param type string, type of plot
#' @param dataPlotType string, type of data plot: 'xy' for a x-y scatterplot,
#'     'ty' for a time series plot of y, 'tx' for a time series plot of x
#' @param ... Optional arguments
#' @return a ggplot
#' @examples
#' x=Segmentation_RecursiveModeling(x=ArdecheRiverGaugings$H,y=ArdecheRiverGaugings$Q)
#' plot(x)
#' @export
plot.recursiveModeling <- function(x,type=c('both','data','shifts'),dataPlotType=c('xy','ty','tx'),...){
  typ=match.arg(type)
  dtyp=match.arg(dataPlotType)
  if(typ=='data'){
    g=plot_recursiveModeling_data(x,type=dtyp)
  } else if(typ=='shifts'){
    g=plot_shifts(x$segmentation)
  } else if(typ=='both'){
    g1=plot_recursiveModeling_data(x,type=dtyp)
    g2=plot_shifts(x$segmentation)
    g=wrap_plots(g1+labs(title='Segmented data and shifts'),
                 g2+labs(title=NULL),ncol=1)+
      plot_layout(axes='collect_x')
  } else {
    g=NULL
  }
  return(g)
}

