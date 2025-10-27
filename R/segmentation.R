#' No-shift likelihood
#'
#' No-shift likelihood function for the quick approximation algorithm.
#'
#' @param theta real vector, parameter vector (mu,sigma)
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @return a numeric value equal to the log-likelihood of the no-shift model
#' @keywords internal
#' @importFrom stats dnorm
quickApprox_llfunk0 <- function(theta,obs,u){
  out=sum(dnorm(obs,mean=theta[1],sd=sqrt(theta[2]^2+u^2),log=TRUE))
  return(out)
}

#' Single-shift likelihood
#'
#' Single-shift likelihood function used in the quick approximation algorithm.
#'
#' @param theta real vector, parameter vector: either (mu1,mu2,sigma) (treated as a fixed-var model) or
#'     (mu1,mu2,sigma1,sigma2) (treated as a shifting-var model)
#' @param tau real, shift time
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @return a numeric value equal to the log-likelihood of the single-shift model
#' @keywords internal
#' @importFrom stats dnorm
quickApprox_llfunk1 <- function(theta,tau,time,obs,u){
  mask=(time<=tau)
  if(length(theta)==4){ # shifting var
    s1=theta[3];s2=theta[4]
  } else {
    s1=s2=theta[3]
  }
  l1=sum(dnorm(obs[mask],mean=theta[1],sd=sqrt(s1^2+u[mask]^2),log=TRUE))
  l2=sum(dnorm(obs[!mask],mean=theta[2],sd=sqrt(s2^2+u[!mask]^2),log=TRUE))
  out=l1+l2
  return(out)
}

#' Compute theta given tau
#'
#' Compute theta (mu's and sigma's) given tau (shift time) by maximizing the log-likelihood,
#' or using explicit estimates if available (which is the case when u==0). Also returns the
#' associated log-likelihood.
#'
#' @param tau real, shift time
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param varShift logical, allow for a shifting variance?
#' @return A list with the following components:
#' \enumerate{
#'   \item par: numeric vector, parameter vector: either (mu1,mu2,sigma) (if varShift is FALSE) or
#'     (mu1,mu2,sigma1,sigma2) (if varShift is TRUE)
#'   \item value: real, corresponding log-likelihood
#'  }
#' @keywords internal
#' @importFrom stats optim sd
quickApprox_getTheta <- function(tau,time,obs,u,varShift){
  mask=(time<=tau)
  start=c()
  start[1]=mean(obs[mask])
  start[2]=mean(obs[!mask])
  if(varShift){
    start[3]=sd(obs[mask])
    start[4]=sd(obs[!mask])
  } else {
    z=obs
    z[mask]=obs[mask]-start[1]
    z[!mask]=obs[!mask]-start[2]
    start[3]=sd(z)
  }
  if(all(u==0)){ # no need to optimize, explicit estimates
    pars=start
    val=quickApprox_llfunk1(theta=start,tau=tau,time=time,obs=obs,u=u)
  } else { # maximize log-lkh
    opt=optim(start,quickApprox_llfunk1,tau=tau,time=time,obs=obs,u=u,control=list(fnscale=-1))
    pars=opt$par
    val=opt$value
  }
  out=list(par=pars,value=val)
  return(out)
}

#' Compute theta given tau
#'
#' Compute theta (mu's and sigma's) given tau (shift time) by maximizing the log-likelihood,
#' or using explicit estimates if available (which is the case when u==0).
#'
#' @inheritParams quickApprox_getTheta
#' @return A numeric vector: either (mu1,mu2,sigma) (if varShift is FALSE) or
#'     (mu1,mu2,sigma1,sigma2) (if varShift is TRUE)
#' @keywords internal
quickApprox_getThetaOnly <- function(tau,time,obs,u,varShift){
  out=quickApprox_getTheta(tau,time,obs,u,varShift)
  return(out$par)
}

#' Default output list
#'
#' Get the default output list returned by Segmentation_Engine when failing.
#'
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @return a list, see ?Recursive_Segmentation for details.
#' @keywords internal
getOutputList_Engine <- function(time,obs,u){
  out=list()
  out$summary=list(data=data.frame(time=time,obs=obs,u=u,
                                   I95_lower=obs-1.96*u,
                                   I95_upper=obs+1.96*u,
                                   period=1),
                   shift=data.frame(tau=numeric(0),I95_lower=numeric(0),I95_upper=numeric(0)))
  out$plot=list(density.tau=NULL,density.inc.tau=NULL)
  out$tau=numeric(0)
  out$segments=numeric(0)
  out$mcmc=data.frame()
  out$data.p=list()
  out$DIC=NA
  out$origin.date.p=min(time)
  return(out)
}

#' Default output list
#'
#' Get the default output list returned by Segmentation when failing.
#'
#' @param time real vector, time
#' @param obs real vector, observations
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @return a list, see ?Segmentation for details.
#' @keywords internal
getOutputList <- function(time,obs,u){
  res=getOutputList_Engine(time,obs,u)
  res$data.p=list(obs.p=res$summary$data$obs,time.p=res$summary$data$time,u.p=res$summary$data$u)
  out=list(summary=res$summary,plot=res$plot,results=list(res),nS=1,
           origin.date=res$origin.date.p)
  return(out)
}

#' Date transformer
#'
#' Transform the dates contained in a raw output list back into their original format.
#'
#' @param out list, output from Segmentation (see ?Segmentation for details)
#' @param origin.date date, origin used to transform date back
#' @return a list, see ?Segmentation for details. It is the same as the list `out`,
#'     but with all dates reformatted.
#' @keywords internal
transformDatesInOutput <- function(out,origin.date){
  # Data time (summary)
  out$summary$data$time = NumericFormatTransform(numeric.date = out$summary$data$time,
                                                 origin.date = origin.date)
  # time series indexed by number of segments identified
  if(is.list(out$data.p$time.p)==T){
    for(i in 1:length(out$data.p$time.p)){
      out$data.p$time.p[[i]] =  NumericFormatTransform(numeric.date = out$data.p$time.p[[i]],
                                                       origin.date = origin.date)
    }
  }else{
    out$data.p$time.p = NumericFormatTransform(numeric.date = out$data.p$time.p,
                                               origin.date = origin.date)
  }
  # Rating shift time summary
  if(all(out$summary$shift$tau!= 0)){
    # Transform all time in POSIXct format
    out$summary$shift <- data.frame(lapply(out$summary$shift, function(column) {
      NumericFormatTransform(numeric.date = column,
                             origin.date = origin.date)
    }))
  }
  # Tau in input date format
  out$tau=out$summary$shift$tau
  # Transform density data in POSIXct format if necessary
  if(NROW(out$summary$shift)>0){
    out$plot$density.tau$Value <- NumericFormatTransform(numeric.date = out$plot$density.tau$Value,
                                                         origin.date = origin.date)
    out$plot$density.inc.tau$tau_lower_inc <- NumericFormatTransform(numeric.date = out$plot$density.inc.tau$tau_lower_inc,
                                                                     origin.date = origin.date)
    out$plot$density.inc.tau$tau_upper_inc <- NumericFormatTransform(numeric.date = out$plot$density.inc.tau$tau_upper_inc,
                                                                     origin.date = origin.date)
    out$plot$density.inc.tau$taU_MAP <- NumericFormatTransform(numeric.date = out$plot$density.inc.tau$taU_MAP,
                                                               origin.date = origin.date)
  }
  out$origin.date.p=origin.date
  return(out)
}

#' Compute critical value
#'
#' Compute the critical value of a likelihood-ratio test for a single step change
#' at an unknown position. Based on the approximations suggested by Gombay and Horvath (1994, 1996a, 1996b, 1997),
#' as detailled in Renard (2006, chapter 3, section I.1.5, p. 101, https://hal.science/tel-02588353).
#' @param d integer>1, difference of dimension between M0 (no-change) and M1 (single-change)
#' @param alpha real in (0;1), type-I error level of the test.
#' @param n integer, sample size
#' @return a real number equal to the critical value of the test.
#' @keywords internal
#' @importFrom stats uniroot
getCriticalValue <- function(d,alpha,n){
  logn=log(n)
  if(d==1){
    h=(logn^1.5)/n
    S=log(((1-h)^2)/(h^2))
    foo=uniroot(BrownianFunk,interval=c(2,10),d=d,h=h,S=S,alpha=alpha)
    out=(foo$root)^2
  } else {
    crit=-log(-0.5*log(1-alpha))
    a=sqrt(2*log(logn))
    b=2*log(logn)+0.5*d*log(log(logn))-lgamma(0.5*d)
    out=((crit+b)/a)^2
  }
  return(out)
}

#' Gombay and Horvath's Brownian function
#'
#' Function to be nullified to compute the critical value when d=1.
#' See Gombay and Horvath (1994, 1996a, 1996b, 1997), and Renard
#' (2006, chapter 3, section I.1.5, p. 101, https://hal.science/tel-02588353).
#' @param z real, critical value
#' @param d integer>1, difference of dimension between M0 (no-change) and M1 (single-change)
#' @param h real, h variable (see suggested references)
#' @param S integer, S variable (noted T in suggested references)
#' @param alpha real in (0;1), type-I error level of the test.
#' @return a real number equal to the evaluated Brownian function.
#' @keywords internal
BrownianFunk <- function(z,d,h,S,alpha){
  out=(1/((2^(0.5*d))*gamma(0.5*d)))*(z^d)*exp(-0.5*z^2)*(S-(d/(z^2))*S+(4/(z^2)))-alpha
  return(out)
}
