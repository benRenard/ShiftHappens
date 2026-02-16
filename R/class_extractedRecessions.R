#***************************************************************************----
# Constructor ----

#' extractedRecessions object constructor.
#'
#' Creates a new instance of a 'extractedRecessions' object, used to
#' store the results of a recession extraction algorithm.
#' @param date Date vector, date.
#' @param time numeric vector, within-recession time (in days), i.e. time is reset to zero at the beginning of each recession.
#' @param H numeric vector, stage.
#' @param uH numeric vector, stage uncertainty expressed as a standard deviation.
#' @param index integer vector, recession index.
#' @return An object of class [extractedRecessions()]: data frame containing the input vectors as columns.
#' @examples
#' rec <- extractedRecessions()
#' @export
extractedRecessions<-function(date=as.Date(numeric(0)),time=numeric(0),
                              H=numeric(0),uH=numeric(0),index=integer(0)){
  o<-new_extractedRecessions(date,time,H,uH,index)
  return(validate_extractedRecessions(o))
}

#***************************************************************************----
# is function ----

#' extractedRecessions object tester
#'
#' Is an object of class 'extractedRecessions'?
#'
#' @param o Object, an object.
#' @return A logical equal to TRUE if class(o)== 'extractedRecessions', FALSE otherwise.
#' @keywords internal
is.extractedRecessions<-function(o){
  return(any(class(o)=='extractedRecessions'))
}

#***************************************************************************----
# internal constructor ----

new_extractedRecessions<-function(date,time,H,uH,index){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(is.numeric(time))
  stopifnot(is.numeric(H))
  stopifnot(is.numeric(uH))
  stopifnot(is_integer(index))
  stopifnot(!is.null(check_equal_length(date,time,H,uH,index)))
  # assemble object
  o=data.frame(date=date,time=time,H=H,uH=uH,index=index)
  class(o) <- c('extractedRecessions','data.frame')
  return(o)
}

#***************************************************************************----
# validator ----

validate_extractedRecessions<-function(x){
  # nothing to do
  return(x)
}

#***************************************************************************----
# plotting ----

#' extractedRecessions plotter
#'
#' Plot a [extractedRecessions()] object.
#' @param x object of class extractedRecessions, resulting from a call to function [Extract_Recessions()]
#' @param type string, type of data plot: 'dh' for date (x) vs. stage (y),
#'     'th' for time (x) vs. stage (y),
#'     'dhmin' for date (x) vs. minimum stage of the recession (y).
#' @param showAnnotation logical, show annotations on plot?
#' @param ... Optional arguments
#' @return A ggplot.
#' @export
#' @examples
#' rec=Extract_Recessions(time=ArdecheRiverStage$Date,H=ArdecheRiverStage$H,
#'                        nSlim=10) # used to speed-up example, but nSlim=1 is recommended.
#' plot(rec)
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr group_by summarise %>%
#' @import ggplot2
plot.extractedRecessions <- function(x,type=c('dh','th','dhmin'),showAnnotation=FALSE,...){
  rec=x
  typ=match.arg(type)
  colors=RColorBrewer::brewer.pal(10,'Paired')
  if(typ %in% c('dh','th')){
    if(typ == 'dh'){
      g=ggplot(rec,aes(x=.data$date,y=.data$H,col=.data$index))
      xlab='Date'
      # position of annotations
      ax=min(rec$date)
      ay=max(rec$H)
    } else {
      g=ggplot(rec,aes(x=.data$time,y=.data$H,col=.data$index))
      xlab='Time [days]'
      # position of annotations
      ax=quantile(rec$time,probs = 0.95)
      ay=quantile(rec$H,probs = 0.99)
    }
    g=g+geom_point()
    if(any(rec$uH>0)){
      g=g+geom_errorbar(ymin=.data$H+stats::qnorm(0.025)*.data$uH,
                        ymax=.data$H+stats::qnorm(0.975)*.data$uH)
    }
    g=g+theme_bw()+
      labs(x = xlab,y = 'Stage [m]')+
      guides(color='none')+
      scale_color_gradientn(colors=colors)
    if(showAnnotation){
      g=g+annotate("text",x=ax,y=ay,label=paste0("Total points = ",NROW(rec)),hjust=0)+
        annotate("text",x=ax,y=ay*0.9,label=paste0("Total recessions = ",max(rec$index)),hjust=0)
    }
  }
  if(typ=='dhmin'){
    DF <- rec %>% group_by(.data$index) %>%
      summarise(H=min(.data$H),date=.data$date[which.min(.data$H)],uH=.data$uH[which.min(.data$H)])
    g=ggplot(DF,aes(x=.data$date,y=.data$H))+geom_point()
    if(any(DF$uH>0)){
      g=g+geom_errorbar(ymin=.data$H+stats::qnorm(0.025)*.data$uH,
                        ymax=.data$H+stats::qnorm(0.975)*.data$uH)
    }
    g=g +theme_bw()+
      labs(x = 'Time [date]', y = 'Minimum recession stage H [m]')
  }
  return(g)
}
