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
#'   \item shifts: data frame, all detected shift times and their uncertainty in numeric or POSIXct format in UTC
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
#' @param results list, vector of 'simpleSegmentation' objects, segmentation results for each tested number of segments
#'
#' @return An object of class [multipleSegmentation()], containing the following fields:
#' \enumerate{
#'   \item nS: integer, optimal number of segments (minimum DIC)
#'   \item DICs: real vector, DICs computed for each number of segment
#'   \item data: data frame, all data with their respective periods after segmentation (for optimal nS)
#'   \item shifts: data frame, all detected shift times and their uncertainty in numeric or POSIXct format in UTC (for optimal nS)
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

#' recursiveSegmentation object constructor.
#'
#' Creates a new instance of a 'recursiveSegmentation' object, used to
#' store the results of a recursive segmentation with an unknown number of segment.
#' @param data data frame, data summary.
#' @param shifts data frame, detected shifts (with columns tau, I95_lower and I95_upper)
#' @param tree [segmentationTree()] object, recursion tree.
#' @param results list, vector of 'multipleSegmentation' objects, segmentation results at each step of the recursion
#'
#' @return An object of class [recursiveSegmentation()], containing the following fields:
#' \enumerate{
#'   \item nS: integer, final number of segments
#'   \item data: data frame, all data with their respective periods after segmentation
#'   \item shifts: data frame, all detected shift times and their uncertainty in numeric or POSIXct format in UTC
#'   \item tree: [segmentationTree()] object, recursion tree.
#'   \item results: list, intermediate results at each stage of the recursion.
#' }
#' @examples
#' sg <- recursiveSegmentation()
#' @export
recursiveSegmentation<-function(data=data.frame(),shifts=data.frame(),
                                tree=segmentationTree(),results=list()){
  o<-new_recursiveSegmentation(data,shifts,tree,results)
  return(validate_recursiveSegmentation(o))
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

#' recursiveSegmentation tester
#'
#' Is an object of class 'recursiveSegmentation'?
#'
#' @param o Object, an object.
#' @return A logical equal to TRUE if class(o)== 'recursiveSegmentation', FALSE otherwise.
#' @keywords internal
is.recursiveSegmentation<-function(o){
  return(class(o)=='recursiveSegmentation')
}

#***************************************************************************----
# internal constructors ----

new_simpleSegmentation<-function(obs,time,u,data,shifts,mcmc,DIC,origin.date){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(is.numeric(obs))
  stopifnot(is.vector(obs))
  stopifnot(is.numeric(u))
  stopifnot(is.vector(u))
  if(!is.null(data)){stopifnot(is.data.frame(data))}
  stopifnot(is.data.frame(shifts))
  stopifnot(is.data.frame(mcmc))
  stopifnot(is.na(DIC) | is.numeric(DIC))
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
  o$nS=NROW(shifts)+1
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
  if(is.null(DICs) | any(is.na(DICs))){nS=1} else {nS=which.min(DICs)}
  o=results[[nS]] # results for optimal number of segments
  o$nS=nS # Add optimal number of segments
  o$DICs=DICs # Add DICs for each number of segments
  o$results=results # Add all partial results for all tested number of segments
  class(o) <- 'multipleSegmentation'
  return(o)
}

new_recursiveSegmentation<-function(data,shifts,tree,results){
  #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  # basic checks
  stopifnot(is.data.frame(data))
  if(!is.null(shifts)){stopifnot(is.data.frame(shifts))}
  stopifnot(is.segmentationTree(tree))
  stopifnot(is.list(results))
  # assemble object
  o=list(nS=NROW(shifts)+1,data=data,shifts=shifts,
         tree=tree,results=results)
  class(o) <- 'recursiveSegmentation'
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

validate_recursiveSegmentation<-function(x){
  # nothing to do
  return(x)
}

#***************************************************************************----
# plotting ----

#' simpleSegmentation plotter
#'
#' Plot a simpleSegmentation object.
#' @param x object of class simpleSegmentation
#' @inheritParams plot_anySegmentation
#' @inherit plot_anySegmentation return
#' @examples
#' x=Segmentation_Engine(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH)
#' plot(x)
#' @export
plot.simpleSegmentation <- function(x,type=c('both','data','shifts'),...){
  g=plot_anySegmentation(x,type,...)
  return(g)
}

#' multipleSegmentation plotter
#'
#' Plot a multipleSegmentation object.
#' @param x object of class multipleSegmentation
#' @inheritParams plot_anySegmentation
#' @inherit plot_anySegmentation return
#' @examples
#' x=Segmentation(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH)
#' plot(x)
#' @export
plot.multipleSegmentation <- function(x,type=c('both','data','shifts'),...){
  g=plot_anySegmentation(x,type,...)
  return(g)
}

#' recursiveSegmentation plotter
#'
#' Plot a recursiveSegmentation object.
#' @param x object of class recursiveSegmentation
#' @inheritParams plot_anySegmentation
#' @inherit plot_anySegmentation return
#' @examples
#' x=Segmentation_Recursive(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH)
#' plot(x)
#' @export
plot.recursiveSegmentation <- function(x,type=c('both','data','shifts'),...){
  g=plot_anySegmentation(x,type,...)
  return(g)
}

#' Segmentation plotter
#'
#' plot a segmentation object of class simpleSegmentation, multipleSegmentation or recursiveSegmentation.
#' @param x object, of class simpleSegmentation, multipleSegmentation or recursiveSegmentation.
#' @param type string, type of plot
#' @param ... Optional arguments
#' @return a ggplot
#' @keywords internal
#' @import ggplot2
#' @importFrom patchwork wrap_plots plot_layout
plot_anySegmentation <- function(x,type=c('both','data','shifts'),...){
  typ=match.arg(type)
  if(typ=='data'){
    g=plot_segmentedData(x)
  } else if(typ=='shifts'){
    g=plot_shifts(x)
  } else if(typ=='both'){
    g1=plot_segmentedData(x)
    g2=plot_shifts(x)
    g=wrap_plots(g1+labs(title='Segmented data and shifts'),
                 g2+labs(title=NULL),ncol=1)+
      plot_layout(axes='collect_x')
  } else {
    g=NULL
  }
  return(g)
}

#' Plot segmented data
#'
#' Plot segmented data after application of a segmentation procedure
#'
#' @param x object, of class simpleSegmentation, multipleSegmentation or recursiveSegmentation.
#' @return a ggplot.
#' @keywords internal
#' @import ggplot2
#' @importFrom scales viridis_pal
plot_segmentedData <- function(x){
  DF=x$data
  colourCount_obs=length(unique(DF$period))
  getPalette_obs=scales::viridis_pal(option='D')
  g=ggplot(data=DF)+
    geom_point(aes(x=.data$time,y=.data$obs,col=factor(.data$period)))
  if(any(DF$u>0)){
    g=g+geom_errorbar(aes(x=.data$time, y=.data$obs,
                          ymin=.data$I95_lower,ymax=.data$I95_upper,
                          col=factor(.data$period)),
                      width=3)
  }
  g=g+scale_color_manual(values=getPalette_obs(colourCount_obs))
  DF=x$shifts
  if(NROW(DF)>0){g=g+geom_vline(data=DF,aes(xintercept=.data$tau))}
  g=g+labs(y='Observation',x='Time',title='Segmented data')
  g=g+theme_bw()+theme(plot.title=element_text(hjust=0.5,face='bold',size=15),
                       legend.position="none")
  return(g)
}

#' Plot shifts
#'
#' Plot shifts detected after application of a segmentation procedure
#'
#' @param x object, of class simpleSegmentation, multipleSegmentation or recursiveSegmentation.
#' @return a ggplot.
#' @keywords internal
#' @import ggplot2
#' @importFrom dplyr left_join starts_with
#' @importFrom tidyr pivot_longer
plot_shifts <- function(x){
  tlim=range(x$data$time)
  shifts=getShifts(x)
  if(NROW(shifts)>0){
    DF=pivot_longer(shifts,cols=starts_with('tau'),names_to='shift')
    if(inherits(x,'recursiveSegmentation')){ #
      foo=data.frame(iteration=x$shifts$id_iteration,shift=paste0('tau',1:NROW(x$shifts)))
      DF=left_join(DF,foo,by='shift')
    } else {
      DF=cbind(DF,iteration=1)
    }
    g=ggplot(DF,aes(x=.data$value,group=.data$shift))+
      geom_density(fill='black',alpha=0.6,color=NA)+
      coord_cartesian(xlim=tlim)+
      facet_grid(rows=vars(.data$iteration),scales='free_y')
  } else {
    g=ggplot()
  }
  g=g+labs(y='Density',x='Time',title='Shifts')+theme_bw()+
    theme(plot.title=element_text(hjust=0.5,face='bold',size=15),legend.position="none",
          axis.text.y=element_blank(),axis.ticks.y=element_blank(),
          panel.grid=element_blank())
  return(g)
}

#***************************************************************************----
# utilities ----
#' Get Shifts
#'
#' Extract shifts (MCMC samples) from a of class simpleSegmentation, multipleSegmentation or recursiveSegmentation
#'
#' @param x object, of class simpleSegmentation, multipleSegmentation or recursiveSegmentation.
#' @param castAsDate boolean. The raw MCMC simulations use an internal numeric representation of time.
#'     If castAsDate is TRUE, they will be transformed into the time/date format that was used in the segmented dataset.
#' @return a dataframe containing the MCMC samples of detected shift times
#' @examples
#' x <- Segmentation(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Date,u=RhoneRiverAMAX$uH)
#' shifts=getShifts(x)
#' summary(shifts)
#' @export
getShifts <- function(x,castAsDate=TRUE){
  if(x$nS==1){return(data.frame())}
  if(inherits(x,'simpleSegmentation') | inherits(x,'multipleSegmentation')){
    nm=names(x$mcmc)
    ix=grep('tau',nm)
    col=nm[ix]
    out=x$mcmc[col]
    if(castAsDate){
      for(i in 1:NCOL(out)){
        out[,i]=numeric_to_time(out[,i],origin.date=x$origin.date)
      }
    }
  } else if(inherits(x,'recursiveSegmentation')){
    its=unique(x$shifts$id_iteration)
    for(i in 1:length(its)){
      mcmc=x$results[[its[i]]]$mcmc
      nm=names(mcmc)
      ix=grep('tau',nm)
      col=nm[ix]
      z=mcmc[col]
      if(castAsDate){
        for(j in 1:NCOL(z)){
          z[,j]=numeric_to_time(z[,j],origin.date=x$results[[its[i]]]$origin.date)
        }
      }
      if(i==1){out=z} else {out=cbind(out,z)}
    }
    names(out) <- paste0('tau',1:NCOL(out))
  } else {
    stop('this function only applies to objects of class simpleSegmentation, multipleSegmentation or recursiveSegmentation')
  }
  return(out)
}
