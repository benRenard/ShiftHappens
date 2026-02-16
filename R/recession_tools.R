#' Recessions extractor
#'
#' Extract recessions from a stage record. NAs are not allowed.
#'
#' @param time POSIXct vector, time. The POSIXct format is mandatory.
#' @param H real vector , stage (m).
#' @param uH real vector, stage uncertainty, expressed as a standard deviation (m).
#' @param deltahMax real value, maximum stage rise allowed within a recession (m).
#'     If the stage record rises by more than deltahMax m,
#'     the current recession is stopped and a new one is started.
#' @param deltatMax real value, maximum time between two stage values allowed within a recession (in days).
#'     When two successive values are separated by more than deltatMax days,
#'     the current recession is stopped and a new one is started.
#'     Useful for dealing with periods of missing data.
#' @param dMin real value, minimum duration (in days) for a recession to be retained in the output dataset.
#' @param nMin integer value, minimum number of data points for a recession to be retained in the output dataset.
#' @param burn numeric value >=0 and <1, burn factor. burn=0.4 means that the first 40 percents of the recession points are discarded.
#' @param deltatMin real value, minimum time between two successive values (in days).
#'     When successive values are too close (less than deltatMin apart), the algorithm will
#'     jump to the next value for speed-up purposes (and hence ignore a stage value).
#'     Useful for quick preliminary runs, but deltatMin=0 is recommended for definitive runs.
#' @param nSlim integer value, initial slimming of the stage record (for speed-up purposes).
#'     nSlim=10 means that only one stage value every 10 is kept before applying the extraction algorithm.
#'     Useful for quick preliminary runs, but nSlim=1 is recommended for definitive runs.
#'
#' @return A data frame with the following columns :
#' \enumerate{
#'    \item date: date
#'    \item time: within-recession time (in days), i.e. time is reset to zero at the beginning of each recession
#'    \item H: stage
#'    \item uH: stage uncertainty (expressed as a standard deviation)
#'    \item index: index of the recession event
#' }
#' @export
#'
#' @examples
#' rec=Extract_Recessions(time=ArdecheRiverStage$Date,H=ArdecheRiverStage$H,
#'                        nSlim=10) # used to speed-up example, but nSlim=1 is recommended.
#' plot(ArdecheRiverStage$Date,ArdecheRiverStage$H,type='l',col='lightgray')
#' points(rec$date,rec$H,col=rec$index)
#' plot(rec$time,rec$H,col=rec$index)
#' @importFrom dplyr between summarise group_by mutate first last slice n %>%
#' @importFrom lubridate is.POSIXct
#' @importFrom stats quantile
Extract_Recessions <- function(time,H,uH=0*H,
                               deltahMax=stats::quantile(H,probs = 0.95,na.rm = TRUE),
                               deltatMax=20,dMin=10,nMin=10,
                               burn=0,deltatMin=0,nSlim=1){
  if(any(is.na(H) | is.na(uH) | is.na(time)))stop('NAs are not alowed in input data (time, H and uH)')
  if(is.null(check_equal_length(H,time)))stop('Input data do not have the same length')
  if(is.null(check_equal_length(uH,time)))stop('Input data do not have the same length')
  if(nSlim  < 1 )stop('nSlim must be greater than or equal to 1')
  if(nSlim != round(nSlim))stop('nSlim must be a integer')
  # if(deltatMin != round(deltatMin))stop('deltatMin must be a integer')
  if(deltatMin < 0)stop('deltatMin must be strictly positive')
  # if(deltatMax != round(deltatMax))stop('deltatMax must be a integer')
  if(deltatMax < 0)stop('deltatMax must be strictly positive')
  # if(dMin != round(dMin))stop('dMin must be a integer')
  if(dMin < 0)stop('dMin must be strictly positive')
  if(nMin < 0)stop('nMin must be strictly positive')
  if(nMin != round(nMin))stop('nMin must be a integer')
  if(!lubridate::is.POSIXct(time))stop('time must be POSIXct format')
  if(!dplyr::between(burn,0,1))stop('burn must be in [0;1)')

  stage.record = data.frame(h=H, uH=uH,t=time)
  # Order by time
  stage.record=stage.record[order(stage.record$t),]
  # slim stage record
  stage.record = stage.record[seq(1, NROW(stage.record), nSlim), ]

  # Find all recession curves
  # initialize:
  hrec = uHrec = data_rec = indx = c()
  trec =  as.POSIXct(character(0))
  hrec[1]  = stage.record$h[1]
  trec[1]  = stage.record$t[1]
  uHrec[1] = stage.record$uH[1]
  indx[1]  = 1
  j = 2  # row index in the original input data
  m = 1  # row index in the output data (extracted recessions)
  k = 1 # recession index
  j_max = length(stage.record$t)
  while (j <= j_max) {
    # duration with previous
    diff_t=as.numeric(difftime(stage.record$t[j],trec[m],units='days'))
    # Check 1: minimum duration between successive values (deltatMin).
    # If successive values are too close, jump to the next for speed-up
    if(diff_t < deltatMin){
      j <- j + 1
      next
    }
    # Check 2: maximum distance between recession data (deltatMax).
    # If time gap is higher than deltatMax, a new recession is started.
    if(diff_t > deltatMax){
      # increase recession index (k) and output data index (m)
      k = k+1
      m=m+1
      # save data
      hrec[m]  = stage.record$h[j]
      trec[m]  = stage.record$t[j]
      uHrec[m] = stage.record$uH[j]
      indx[m]  = k
      j <- j + 1
      next
    }

    # Check 3: Delta h
    diff_h=stage.record$h[j]-hrec[m]
    # Delta h is negative ->  save data.
    if(diff_h < 0){
      m=m+1
      hrec[m]  = stage.record$h[j]
      trec[m]  = stage.record$t[j]
      uHrec[m] = stage.record$uH[j]
      indx[m]  = k
    } else { # Case when delta h is positive
      # if delta h is larger than deltahMax, start new recession
      if(diff_h >= deltahMax){
        m=m+1
        k = k+1
        hrec[m]  = stage.record$h[j]
        trec[m]  = stage.record$t[j]
        uHrec[m] = stage.record$uH[j]
        indx[m]  = k
      }
      # if delta h is smaller than deltahMax, do not take into account this data
      # and keep going. So j is increasing but m remains the same.
    }
    j <- j + 1
  }

  data_df=data.frame(H=hrec,
                     date=trec,
                     uH=uHrec,
                     index=indx)

  summary_temp <- data_df %>% group_by(.data$index) %>%
                  summarise(first_time = first(.data$date),  # First date for each group
                            last_time = last(.data$date),    # Last date for each group
                            count_values = n())  # Total number of values for each group
  # recession duration
  summary_temp <- summary_temp %>%
    mutate(diff_time=as.numeric(difftime(.data$last_time,.data$first_time,units="days")))

  # Check 4 : dMin
  # Select only recessions with duration >= dMin
  rm_indx_dMin = which(summary_temp$diff_time < dMin)
  # Check 5 : Nmin
  # Select only recessions with a number of point >= Nmin :
  rm_indx_Nmin = which(summary_temp$count_values < nMin)
  rm_indx=unique(c(rm_indx_dMin,rm_indx_Nmin))
  out <- data_df[!data_df$index %in% rm_indx, ]

  # Re-assign the recession index
  # Step 1: Get unique values of 'index' and sort them
  unique_indx <- unique(out$index)
  # Step 2: Create a mapping from old 'indx' values to new sequential values
  indx_map <- data.frame(old = unique_indx, new = seq_along(unique_indx))
  # Step 3: Replace 'index' values in the filtered data frame with new sequential values
  out$index <- indx_map$new[match(out$index, indx_map$old)]
  if(burn!=0){
    # Burn the first part of the recession to reduce data information
    out <- out %>% group_by(.data$index) %>% slice(-(1:floor(n() * burn)))
  }
  out <- out %>% group_by(.data$index) %>%
    mutate(time=as.numeric(difftime(.data$date,first(.data$date),units = "days")))

  return(out[c('date','time','H','uH','index')])
}

#' Plot recessions
#'
#' Plot recessions extracted using function [Extract_Recessions()]
#'
#' @param rec data frame, resulting from a call to function [Extract_Recessions()]
#' @param type string, type of data plot: 'dh' for date (x) vs. stage (y),
#'     'th' for time (x) vs. stage (y),
#'     'dhmin' for date (x) vs. minimum stage of the recession (y).
#' @param showAnnotation logical, show annotations on plot?#'
#' @return A ggplot.
#' @export
#' @examples
#' rec=Extract_Recessions(time=ArdecheRiverStage$Date,H=ArdecheRiverStage$H,
#'                        nSlim=10) # used to speed-up example, but nSlim=1 is recommended.
#' PlotRecessions(rec)
#' @importFrom RColorBrewer brewer.pal
#' @importFrom dplyr group_by summarise %>%
#' @import ggplot2
PlotRecessions <- function(rec,type=c('dh','th','dhmin'),showAnnotation=FALSE){
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

