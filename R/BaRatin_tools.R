#' BaRatin
#'
#' Fit a rating curve with BaRatin
#'
#' @param x real vector, stage
#' @param y real vector, discharge
#' @param time vector (numeric or date), time
#' @param uX real vector, stage uncertainty (as a standard deviation)
#' @param uY real vector, discharge uncertainty (as a standard deviation)
#' @param flavor string, BaRatin flavor
#' @param controlMatrix square matrix, control matrix
#' @param priors list, list of RBaM::parameter objects containing information for each
#'     parameter of the rating curve. Length 3*NROW(controlMatrix)
#' @param mcmc_options mcmcOptions object, see ?RBaM::mcmcOptions
#' @param mcmc_cooking mcmcCooking object, see ?RBaM::mcmcCooking
#' @param remnant list of one remnantErrorModel object, see ?RBaM::remnantErrorModel
#' @param temp.folder string, directory where config and result files are stored.
#' @inherit fittedModel return
#' @examples
#' f=Fit_BaRatin(x=ArdecheRiverGaugings$H,y=ArdecheRiverGaugings$Q,
#'               uY=ArdecheRiverGaugings$uQ,time=ArdecheRiverGaugings$Date)
#' if(!is.null(f)){
#'   plot(f$data$x,f$data$obs);points(f$data$x,f$data$sim,col='red')
#'   plot(f$data$time,f$data$res)
#' }
#' @importFrom RBaM dataset xtraModelInfo model BaM prediction
#' @export
Fit_BaRatin <- function(x,y,time=1:NROW(y),uX=0*x,uY=0*y,
                        flavor=c('BaRatin','BaRatinBAC'),
                        controlMatrix=matrix(1),
                        priors=get3FlatPriors(x,y),
                        mcmc_options=RBaM::mcmcOptions(),
                        mcmc_cooking=RBaM::mcmcCooking(),
                        remnant=list(RBaM::remnantErrorModel()),
                        temp.folder=file.path(tempdir(),'BaRatin')){
  if(is.null(check_equal_length(x,y,time,uX,uY)))stop('Arguments x, y, time, uX, uY should all have the same length')
  if(!is.matrix(controlMatrix))stop('Control matrix must be a matrix')
  nControl=NROW(controlMatrix)
  if(NCOL(controlMatrix)!=nControl)stop('Hydraulic control must be square')
  if(any(controlMatrix!=0 & controlMatrix!=1))stop('Hydraulic control must be filled with 1 (active) and 0 (inactive) for describing hydraulic controls')
  if(length(priors)!=3*nControl)stop('The length of the priors list should be three times the number of controls')
  flav=match.arg(flavor)

  # Define the calibration dataset
  data=data.frame(time=time,H=x,Q=y,uQ=uY)
  names(data) <- c('time','H','Q','uQ')
  D=RBaM::dataset(X=data['H'],Y=data['Q'],Yu=data['uQ'],data.dir=temp.folder)
  # Config_xtra
  xtra=RBaM::xtraModelInfo(object=controlMatrix)
  # Stitch it all together into a model object
  M=RBaM::model(ID=flav,
                nX=1,nY=1, # number of input/output variables
                par=priors, # list of model parameters
                xtra=xtra) # use xtraModelInfo() to pass the control matrix
  # Run BaM executable to calibrate RC
  ok=try(RBaM::BaM(mod=M,data=D,workspace=temp.folder,
                   mcmc=mcmc_options,cook=mcmc_cooking,remnant=remnant))
  if(inherits(ok, "try-error")){return()}
  # Prepare predictions
  if(all(uX==0)){
    xpred=data['H']
  } else { # propagate stage uncertainty
    # Create 100 random values for x, with the uncertainty given by the user
    xrep=matrix(rnorm(n=length(uX)*100,mean=x,sd=uX),nrow=length(uX),ncol=100)
    xpred=list(xrep)
  }
  totalU=RBaM::prediction(X=xpred,doParametric=TRUE,doStructural=TRUE,
                          spagFiles='QRC_TotalU.spag',data.dir=temp.folder)
  # Run BaM executable in prediction mode
  ok=try(RBaM::BaM(mod=M,data=D,workspace=temp.folder,
                   mcmc=mcmc_options,cook=mcmc_cooking,remnant=remnant,
                   pred=totalU,doCalib=FALSE,doPred=TRUE))
  if(inherits(ok, "try-error")){return()}
  # Read and interpret result files
  ok=try({resid=utils::read.table(file=file.path(temp.folder,"Results_Residuals.txt"),header=TRUE);
          resum=utils::read.table(file=file.path(temp.folder,"Results_Summary.txt"),header=TRUE);
          env=utils::read.table(file.path(temp.folder,'QRC_TotalU.env'),header=TRUE)})
  if(inherits(ok, "try-error")){return()}
  # Assemble returned object
  DF=data.frame(time=time,x=x,obs=y,sim=resid$Y1_sim,
                res=resid$Y1_res,uRes=sqrt(as.matrix(uY)^2+env$Stdev^2))
  names(DF) <- c('time','x','obs','sim','res','uRes')
  pars=resum[NROW(resum),1:(3*nControl)]
  out=fittedModel(data=DF,parameters=as.numeric(pars))
  return(out)
}

#' Get BaRatin Flat Prior
#'
#' Get a default flat prior for the 3 parameters of a 1-control BaRatin Rating curve.
#' Reasonable starting values are computed from a log-linear regression.
#'
#' @param h real vector, stage values
#' @param Q real vector, discharge values
#' @return a list of length 3 containing the RBaM::parameter objects for k, a and c.
#' @keywords internal
#' @importFrom RBaM parameter
get3FlatPriors <- function(h,Q){
  if(length(h)>1){ # Do a log-linear regression to get reasonable starting points
    h0=min(h)-0.01*diff(range(h)) # a bit below the smallest stage
    logdepth=log(h-h0)
    w=lm(log(Q)~logdepth)
    pars=w$coefficients
    out=list(RBaM::parameter('k',init=h0),
             RBaM::parameter('a',init=exp(pars[1])),
             RBaM::parameter('c',init=pars[2]))
  } else { # order-of-magnitude starting points
    out=list(RBaM::parameter('k',init=0),
             RBaM::parameter('a',init=10),
             RBaM::parameter('c',init=1.5))
  }
  return(out)
}
