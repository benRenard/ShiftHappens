#***************************************************************************----
# Constructor ----

#' recursiveModeling object constructor.
#'
#' Creates a new instance of a 'recursiveModeling' object, used to
#' store the results of a recursive modeling segmentation
#' @param data data frame, initial dataset
#' @param segmentation [recursiveSegmentation()] object, segmentation results
#' @param fit list of [fittedModel()] objects, fitted model results
#' @return An object of class [recursiveModeling()], containing the following fields:
#' \enumerate{
#'   \item data: a data frame containing the original dataset used as input of the recursive modeling
#'   \item segmentation: an object of class [recursiveSegmentation()],
#'       containing the results of the segmentation of residuals
#'   \item fit: a list of objects of class [fittedModel()],
#'       containing the results of the model fit at each stage of the recursion
#' }
#' @examples
#' rec <- recursiveModeling()
#' @export
recursiveModeling<-function(data=data.frame(),segmentation=recursiveSegmentation(),fit=list(fittedModel())){
  o<-new_recursiveModeling(data,segmentation,fit)
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

new_recursiveModeling<-function(data,segmentation,fit){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(is.data.frame(data))
  stopifnot(inherits(segmentation,'recursiveSegmentation'))
  stopifnot(is.list(fit))

  for (i in 1:length(fit)){
    stopifnot(inherits(fit[[i]],'fittedModel'))
  }
  # assemble object
  o=list(data=data,segmentation=segmentation,fit=fit)
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
  terminal=x$segmentation$tree$index[x$segmentation$tree$nS==1]
  DF=c()
  nX=length(grep('x',names(x$data)))
  for(i in 1:length(terminal)){
    node=terminal[i]
    foo=x$fit[[node]]$data
    DF=rbind(DF,cbind(foo,period=i))
  }
  DF$period=as.factor(DF$period)
  # 2DO: handle multi-X case
  if(nX>1 & typ %in% c('xy','tx')){
    warning(paste('This plotting function does not handle multiple predictors yet.',
                  'Only the first predictor x1 is represented in the plot'))
  }
  colourCount_obs=length(unique(DF$period))
  getPalette_obs=scales::viridis_pal(option='D')
  if(typ=='xy'){
    g=ggplot(DF)+
      geom_point(aes(x=.data$x1,y=.data$y,group=.data$period,color=.data$period))+
      geom_line(aes(x=.data$x1,y=.data$ysim,group=.data$period,color=.data$period))
    if(any(x$data$uY>0)){
      g=g+geom_errorbar(data=x$data,aes(x=.data$x1,
                                        ymin=.data$y-1.96*.data$uY,
                                        ymax=.data$y+1.96*.data$uY,
                                        color=.data$period,group=.data$period))
    }
    if(any(x$data$uX1>0)){
      g=g+geom_errorbarh(data=x$data,aes(y=.data$y,
                                         xmin=.data$x1-1.96*.data$uX1,
                                         xmax=.data$x1+1.96*.data$uX1,
                                         color=.data$period,group=.data$period))
    }
  } else if(typ=='ty') {
    g=ggplot(x$data)+
      geom_point(aes(x=.data$time,y=.data$y,group=.data$period,color=.data$period))
    if(any(x$data$uY>0)){
      g=g+geom_errorbar(aes(x=.data$time,
                            ymin=.data$y-1.96*.data$uY,
                            ymax=.data$y+1.96*.data$uY,
                            color=.data$period,group=.data$period))
    }
    if(NROW(x$segmentation$shifts)>0){
      g=g+geom_vline(data=x$segmentation$shifts,aes(xintercept=.data$tau))
    }
  } else if(typ=='tx'){
    g=ggplot(x$data)+
      geom_point(aes(x=.data$time,y=.data$x1,group=.data$period,color=.data$period))
    if(any(x$data$uX1>0)){
      g=g+geom_errorbar(aes(x=.data$time,
                            ymin=.data$x1-1.96*.data$uX1,
                            ymax=.data$x1+1.96*.data$uX1,
                            color=.data$period,group=.data$period))
    }
    if(NROW(x$segmentation$shifts)>0){
      g=g+geom_vline(data=x$segmentation$shifts,aes(xintercept=.data$tau))
    }
  }
  g=g+scale_color_manual(values=getPalette_obs(colourCount_obs))
  g=g+labs(y=ylab,x=xlab,title='Segmented data')
  g=g+theme_bw()+theme(plot.title=element_text(hjust=0.5,face='bold',size=15),
                       legend.position="none")
  return(g)
}

#' recursiveModeling plotter
#'
#' Plot a [recursiveModeling()] object.
#' @param x object of class recursiveModeling
#' @param type string, type of plot: 'both' for a data panel on top and a shift panel on the bottom,
#'     'data' for the data panel only, 'shifts' for the shift panel only,
#'     'fits' for a plot of fitted models in each terminal mode.
#' @param dataPlotType string, type of data plot: 'xy' for a x-y scatterplot,
#'     'ty' for a time series plot of y, 'tx' for a time series plot of x
#' @param ... Optional arguments
#' @return a ggplot (or a list of ggplots if type='fits')
#' @examples
#' x=Segmentation_RecursiveModeling(x=ArdecheRiverGaugings$H,y=ArdecheRiverGaugings$Q)
#' plot(x)
#' @export
plot.recursiveModeling <- function(x,type=c('both','data','shifts','fits'),dataPlotType=c('xy','ty','tx'),...){
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
  } else if(typ=='fits'){
    terminal=x$segmentation$tree$index[x$segmentation$tree$nS==1]
    g=vector('list',length(terminal))
    for(i in 1:length(terminal)){
      g[[i]]=plot(x$fit[[terminal[i]]])+labs(title=paste('Node',terminal[i]))
    }
  } else {
    g=NULL
  }
  return(g)
}

