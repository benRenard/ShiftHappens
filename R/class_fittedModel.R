#***************************************************************************----
# Constructor ----

#' fittedModel object constructor.
#'
#' Creates a new instance of a 'fittedModel' object, used to
#' store the results of a fitted model
#'
#' @param time numeric vector, time
#' @param x numeric vector or matrix, observed inputs (predictors)
#' @param y numeric vector, observed output (predictand)
#' @param ysim numeric vector, simulated output
#' @param res numeric vector, residual
#' @param uX numeric vector or matrix, uncertainty in observed inputs (expressed as a standard deviation)
#' @param uY numeric vector, uncertainty in observed outputs (expressed as a standard deviation)
#' @param uYsim numeric vector, uncertainty in simulated outputs (expressed as a standard deviation)
#' @param uRes numeric vector, residual uncertainty (expressed as a standard deviation)
#' @param group factor, grouping variable
#' @param parameters numeric vector, fitted parameters
#' @param uParameters numeric vector, standard deviation or standard error of fitted parameters
#' @return An object of class [fittedModel()], containing the following fields:
#' \enumerate{
#'   \item data: data frame, column-binding all the arguments of the function except parameters
#'   \item parameters: data frame, fitted parameters + standard error
#' }
#' @examples
#' f <- fittedModel()
#' @export
fittedModel<-function(time=numeric(0),x=numeric(0),y=numeric(0),ysim=numeric(0),res=y-ysim,
                      uX=0*x,uY=0*y,uYsim=0*ysim,uRes=sqrt(uY^2+uYsim^2),
                      group=factor(rep(1,length(time))),
                      parameters=numeric(0),uParameters=NA*parameters){
  o<-new_fittedModel(time,x,y,ysim,res,uX,uY,uYsim,uRes,group,parameters,uParameters)
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

new_fittedModel<-function(time,x,y,ysim,res,uX,uY,uYsim,uRes,group,parameters,uParameters){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(y))
  stopifnot(is.numeric(ysim))
  stopifnot(is.numeric(res))
  stopifnot(is.numeric(uX))
  stopifnot(is.numeric(uY))
  stopifnot(is.numeric(uYsim))
  stopifnot(is.numeric(uRes))
  stopifnot(!is.null(check_equal_length(time,y,ysim,res,uY,uYsim,uRes,group)))
  stopifnot(NROW(x)==length(time))
  stopifnot(NROW(uX)==length(time))
  # assemble object
  DF=data.frame(time,x,y,ysim,res,uX,uY,uYsim,uRes,group)
  names(DF) <- c('time',paste0('x',1:NCOL(x)),'y','ysim','res',
                 paste0('uX',1:NCOL(uX)),'uY','uYsim','uRes','group')
  param=data.frame(value=as.numeric(parameters)[],u=as.numeric(uParameters))
  names(param) <- c('value','u')
  o=list(data=DF,parameters=param)
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


#' fittedModel Plotter.
#'
#' Plot a [fittedModel()] object.
#' @param x object of class fittedModel
#' @param type string, type of plot: 'xy' for a x-y scatterplot,
#'     'ty' for a time series plot of y and ysim,
#'     'tres' for a time series plot of residuals,
#'     'yysim' for a scatterplot y vs. ysim
#' @param ... Optional arguments
#' @return a ggplot
#' @examples
#' f=Fit_LinearRegression(x=RhoneRiverAMAX$H,y=RhoneRiverAMAX$Q,time=RhoneRiverAMAX$Date)
#' plot(f)
#' plot(f,type='ty')
#' plot(f,type='tres')
#' plot(f,type='yysim')
#' @export
#' @importFrom dplyr arrange %>%
#' @import ggplot2
plot.fittedModel <- function(x,type=c('xy','ty','tres','yysim'),...){
  typ=match.arg(type)
  if(typ=='xy'){
    xlab='x';ylab='y'
  } else if (typ=='ty'){
    xlab='Time';ylab='y'
  } else if (typ=='tres'){
    xlab='Time';ylab='Residual'
  } else if (typ=='yysim'){
    xlab='Observed';ylab='Simulated'
  } else {
    stop('Unknown plot type')
  }
  nX=length(grep('x',names(x$data)))
  # 2DO: handle multi-X case
  if(nX>1 & typ %in% c('xy')){
    warning(paste('This plotting function does not handle multiple predictors yet.',
                  'Only the first predictor x1 is represented in the plot'))
  }
  DF=x$data %>% arrange(.data$x1)
  DF$group=factor(DF$group)
  colourCount_obs=length(unique(DF$group))
  getPalette_obs=scales::viridis_pal(option='D')
  alf=0.2
  if(typ=='xy'){
    g=ggplot(DF)+
      geom_ribbon(aes(x=.data$x1,ymin=.data$ysim-1.96*.data$uYsim,
                      ymax=.data$ysim+1.96*.data$uYsim,
                      group=.data$group,fill=.data$group),alpha=alf)+
      geom_line(aes(x=.data$x1,y=.data$ysim,group=.data$group,color=.data$group))+
      geom_point(aes(x=.data$x1,y=.data$y,group=.data$group,color=.data$group))
    if(any(DF$uY>0)){
      g=g+geom_errorbar(aes(x=.data$x1,ymin=.data$y-1.96*.data$uY,
                            ymax=.data$y+1.96*.data$uY,
                            group=.data$group,color=.data$group),width=0)
    }
    if(any(DF$uX1>0)){
      g=g+geom_errorbarh(aes(y=.data$y,xmin=.data$x1-1.96*.data$uX1,
                             xmax=.data$x1+1.96*.data$uX1,
                             group=.data$group,color=.data$group))
    }
  } else if(typ=='ty'){
    g=ggplot(DF)+
      geom_ribbon(aes(x=.data$time,ymin=.data$ysim-1.96*.data$uYsim,
                      ymax=.data$ysim+1.96*.data$uYsim,
                      group=.data$group,fill=.data$group),alpha=alf)+
      geom_line(aes(x=.data$time,y=.data$ysim,group=.data$group,color=.data$group))+
      geom_point(aes(x=.data$time,y=.data$y,group=.data$group,color=.data$group))
    if(any(DF$uY>0)){
      g=g+geom_errorbar(aes(x=.data$time,ymin=.data$y-1.96*.data$uY,
                            ymax=.data$y+1.96*.data$uY,
                            group=.data$group,color=.data$group),width=0)
    }
  } else if(typ=='tres'){
    g=ggplot(DF)+
      geom_point(aes(x=.data$time,y=.data$res,group=.data$group,color=.data$group))+
      geom_errorbar(aes(x=.data$time,ymin=.data$res-1.96*.data$uRes,
                        ymax=.data$res+1.96*.data$uRes,group=.data$group,color=.data$group),width=0)
  } else if(typ=='yysim'){
    g=ggplot(DF)+
      geom_abline(slope=1,intercept=0)+
      geom_point(aes(x=.data$y,y=.data$ysim,group=.data$group,color=.data$group))
    if(any(DF$uYsim>0)){
      g=g+geom_errorbar(aes(x=.data$y,ymin=.data$ysim-1.96*.data$uYsim,
                            ymax=.data$ysim+1.96*.data$uYsim,
                            group=.data$group,color=.data$group),width=0)
    }
    if(any(DF$uY>0)){
      g=g+geom_errorbarh(aes(y=.data$ysim,xmin=.data$y-1.96*.data$uY,
                             xmax=.data$y+1.96*.data$uY,
                             group=.data$group,color=.data$group))
    }
    lim=range(c(DF$y,DF$ysim))
    g=g+coord_equal(xlim=lim,ylim=lim)
  }
  g=g+scale_color_manual(values=getPalette_obs(colourCount_obs))
  g=g+scale_fill_manual(values=getPalette_obs(colourCount_obs))
  g=g+labs(y=ylab,x=xlab)
  g=g+theme_bw()+theme(plot.title=element_text(hjust=0.5,face='bold',size=15),
                       legend.position="none")
  return(g)
}
