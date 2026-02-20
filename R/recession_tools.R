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
#' @return An object of class [extractedRecessions()], which is a data frame with the following columns:
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
  if(is.null(check_equal_length(H,uH,time)))stop('Input data do not have the same length')
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
  return(extractedRecessions(date=out$date,time=out$time,H=out$H,uH=out$uH,index=out$index))
}

#' Fit recession model
#'
#' Fit a model to a set of recession events.
#'
#' @param rec [extractedRecessions()] object, resulting from a call to function [Extract_Recessions()]
#' @param equation string, equation used to model recessions, from M1 to M9.
#'     For more details see Chapter 3 of the PhD thesis of Matteo Darienzo
#'     (2021, \url{https://theses.hal.science/tel-03211343}) or call getRecessionEquations()
#' @inheritParams Fit_BaRatin
#' @inherit fittedModel return
#' @source \url{https://theses.hal.science/tel-03211343}
#' @examples
#' rec=Extract_Recessions(time=ArdecheRiverStage$Date,H=ArdecheRiverStage$H,
#'                        nSlim=10) # used to speed-up example, but nSlim=1 is recommended.
#' f=Fit_Recessions(rec=rec,equation='M7',
#'                  # MCMC options are modified to speed-up example.
#'                  # Using default MCMC options is safer and is recommended.
#'                  mcmc_options=mcmcOptions(nAdapt=20,nCycles=10))
#' if(!is.null(f)){
#'   plot(f)
#'   plot(f,type='ty')
#' }
#' @importFrom RBaM dataset xtraModelInfo model BaM mcmcOptions mcmcCooking
#' @importFrom utils read.table
#' @export
Fit_Recessions <- function(rec,equation=c('M6','M1','M2','M3','M4','M5','M7','M8','M9'),
                           mcmc_options=RBaM::mcmcOptions(),
                           mcmc_cooking=RBaM::mcmcCooking(),
                           remnant=list(getConstantRemnant(rec)),
                           temp.folder=file.path(tempdir(),'Recessions')){
  eq=match.arg(equation)
  # RBaM dataset object
  D=RBaM::dataset(X=data.frame(t=rec$time),Y=data.frame(H=rec$H),
                  Yu=data.frame(uH=rec$uH),VAR.indx=data.frame(index=rec$index),
                  data.dir=temp.folder)
  Ncurves=max(rec$index)
  # get parameter list
  param=getRecessionParameters(rec,eq,D)
  # Use xtraModelInfo to pass the names of the inputs and the formulas
  xtra=RBaM::xtraModelInfo(object=list(inputs=c('t'),formulas=getRecessionEquations()[[eq]]))
  # BaM model object
  M=RBaM::model(ID='TextFile',nX=1,nY=1,par=param,xtra=xtra)
  # Run BaM executable to calibrate RC
  ok=try(RBaM::BaM(mod=M,data=D,workspace=temp.folder,
                   mcmc=mcmc_options,cook=mcmc_cooking,remnant=remnant))
  if(inherits(ok, "try-error")){return()}
  # Read results
  res=utils::read.table(file=file.path(temp.folder,"Results_Residuals.txt"),header=TRUE)
  resume=utils::read.table(file=file.path(temp.folder,"Results_Summary.txt"),header=TRUE)
  # Assemble output
  stot=res$Y1_res/res$Y1_stdres # total uncertainty
  out=fittedModel(time=rec$date,x=rec$time,y=rec$H,ysim=res$Y1_sim,res=res$Y1_res,
                uY=rec$uH,uYsim=sqrt(stot^2-rec$uH^2),uRes=stot,group=rec$index,
                parameters=resume[16,],uParameters=resume[11,])
  return(out)
}

#' Recessions equations
#'
#' Get all available recession equations.
#' For details see Chapter 3 of the PhD thesis of Matteo Darienzo
#' (2021, \url{https://theses.hal.science/tel-03211343}).
#'
#' @return A list containing all available equations.
#' @source \url{https://theses.hal.science/tel-03211343}
#' @examples
#' getRecessionEquations()
#' @export
getRecessionEquations <- function(){
  out=list(
    M1='alpha_k*exp(-lambda*t)+beta_k',
    M2='alpha1_k*exp(-lambda1*t)+alpha2*exp(-lambda2*t)+beta_k',
    M3='alpha1_k*exp(-lambda1*t)+alpha2_k*exp(-lambda2*t)+beta_k',
    M4='alpha1_k*exp(-lambda1*t)+alpha2*exp(-lambda2*t)+alpha3*exp(-lambda3*t)+beta_k',
    M5='alpha1_k*exp(-lambda1*t)+alpha2_k*exp(-lambda2*t)+alpha3*exp(-lambda3*t)+beta_k',
    M6='alpha_k*exp(-lambda*t^eta)+beta_k',
    M7='alpha_k*exp(-lambda_k*t^eta)+beta_k',
    M8='alpha_k/((1+lambda*t)^eta)+beta_k',
    M9='alpha_k/((1+lambda_k*t)^eta)+beta_k'
  )
  return(out)
}

#' Recessions-based segmentation
#'
#' Segmentation based on a segmentation of the lowest points reached by stage recessions.
#' For more details see Chapter 3 of the PhD thesis of Matteo Darienzo
#' (2021, \url{https://theses.hal.science/tel-03211343}).
#'
#' @param equation string, equation used to model recessions. If 'none' (default), the minimum
#'     of each recession is extracted directly (hence no recession model is used) and seegmented.
#'     If a model between 'M1' and 'M9' is used, recessions are modeled and the 'asymptoti stage'
#'     parameter is segmented. For more details on the model equations, see Chapter 3 of the PhD thesis
#'     of Matteo Darienzo (2021, \url{https://theses.hal.science/tel-03211343}) or call getRecessionEquations().
#' @inheritParams Fit_Recessions
#' @inheritParams Segmentation_Recursive
#' @inherit fittedModel return
#' @source \url{https://theses.hal.science/tel-03211343}
#' @examples
#' rec=Extract_Recessions(time=ArdecheRiverStage$Date,H=ArdecheRiverStage$H,
#'                        nSlim=10) # used to speed-up example, but nSlim=1 is recommended.
#'  sg=Segmentation_Recessions(rec)
#'  plot(sg)
#' @export
Segmentation_Recessions <- function(rec,equation=c('none','M1','M2','M3','M4','M5','M6','M7','M8','M9'),
                                    nSmax=2,doQuickApprox=TRUE,nMin=ifelse(doQuickApprox,3,1),
                                    nSim=500,varShift=FALSE,alpha=0.1,
                                    mcmc_options=RBaM::mcmcOptions(),
                                    mcmc_cooking=RBaM::mcmcCooking(),
                                    mu_prior=list(),
                                    remnant=list(getConstantRemnant(rec)),
                                    temp.folder=file.path(tempdir(),'Recessions')){
  # Get recession lows to be segmented
  eq=match.arg(equation)
  recMin=getRecessionMin(rec)
  Ncurves=NROW(recMin)
  DF=data.frame(index=recMin$index,date=recMin$date)
  if(eq=='none'){
    DF=cbind(DF,value=recMin$H,u=recMin$uH)
  } else {
    f=Fit_Recessions(rec=rec,equation=eq,
                     mcmc_options=mcmc_options,mcmc_cooking=mcmc_cooking,
                     remnant=remnant,temp.folder=temp.folder)
    DF=cbind(DF,value=f$parameters$value[1:Ncurves],u=f$parameters$u[1:Ncurves])
  }
  # Segment
  out=Segmentation_Recursive(obs=DF$value,time=DF$date,u=DF$u,
                             nSmax=nSmax,doQuickApprox=doQuickApprox,nMin=nMin,
                             nSim=nSim,varShift=varShift,alpha=alpha,
                             mcmc_options=mcmc_options,mcmc_cooking=mcmc_cooking,
                             temp.folder=temp.folder,mu_prior=mu_prior)
  return(out)
}

#***************************************************************************----
# internal functions ----

#' Get constant remnant error model
#'
#' Get a remnant error model corresponding to a constant standard deviation.
#'
#' @param rec [extractedRecessions()] object, resulting from a call to function [Extract_Recessions()]
#' @return a remnantErrorModel object, see ?RBaM::remnantErrorModel
#' @keywords internal
#' @importFrom RBaM parameter remnantErrorModel
#' @importFrom stats sd
getConstantRemnant <- function(rec){
  param=RBaM::parameter(name='sigma',init=stats::sd(rec$H),prior.dist='FlatPrior+')
  out=RBaM::remnantErrorModel(funk='Constant',par=list(param))
  return(out)
}

#' Get recession parameters
#'
#' Get a list of parameters corresponding to a given recession equation.
#'
#' @param rec [extractedRecessions()] object, resulting from a call to function [Extract_Recessions()]
#' @param eq string, the requested equation
#' @param D RBaM dataset object, ?RBaM::dataset
#' @return a list of parameter objects, see ?RBaM::parameter
#' @keywords internal
#' @importFrom RBaM parameter
#' @importFrom stats median
#' @importFrom dplyr group_by summarise %>%
getRecessionParameters <- function(rec,eq,D){
  Ncurves=max(rec$index)
  foo=rec %>% group_by(.data$index) %>% summarise(minH=min(.data$H),maxH=max(.data$H))
  # beta_k is common to all models
  beta_k=RBaM::parameter_VAR(name='beta_k',index='index',d=D,
                             # Specify the parameter's initial guess and prior FOR EACH INDEX
                             init=foo$minH,prior.dist=rep('FlatPrior',Ncurves),prior.par =rep(list(NULL), Ncurves))
  # Model-specific parameters
  if(eq %in% c('M6','M7')){
    eta=RBaM::parameter(name='eta',init=1,prior.dist='FlatPrior+')
    alpha_k=RBaM::parameter_VAR(name='alpha_k',index='index',d=D,
                                # The next 3 lines specify the parameter's initial guess and prior FOR EACH INDEX
                                init=foo$maxH-foo$minH,
                                prior.dist=rep('FlatPrior+',Ncurves),
                                prior.par =rep(list(NULL), Ncurves))
    if(eq=='M6'){
      lambda=RBaM::parameter(name='lambda',
                             init=log(2)/stats::median(rec$time), # see https://en.wikipedia.org/wiki/Exponential_decay
                             prior.dist='FlatPrior+')
      out=list(beta_k,alpha_k,lambda,eta)
    } else {
      lambda_k=RBaM::parameter_VAR(name='lambda_k',index='index',d=D,
                             init=rep(log(2)/stats::median(rec$time),Ncurves),
                             prior.dist=rep('FlatPrior+',Ncurves),prior.par =rep(list(NULL), Ncurves))
      out=list(beta_k,alpha_k,lambda_k,eta)
    }
  } else {
    stop('The requested recession equation is unknown or is not yet implemented')
  }
  return(out)
}
