#' Segmentation engine - quick approximation algorithm
#'
#' A quick approximation to the segmentation procedure for either one or two segments, and
#' no prior information. The approximation is based on a likelihood-ratio test for a single step change
#' at an unknown position. DIC0 is hence replaced with the deviance of a no-shift model M0,
#' while DIC1 is replaced by the maximum deviance of a single-shift model M1 + the critical value of the test.
#' This way, a change will be detected when DIC1<DIC0, i.e. when deviance1+critical value<deviance0,
#' i.e. when deviance0-deviance1>critical value, which precisely corresponds to the outcome of the likelihood ratio test.
#' \cr
#' Critical values are computed using the approximations suggested by Gombay and Horvath (1994, 1996a, 1996b, 1997),
#' as detailled in Renard (2006, chapter 3, section I.1.5, p. 101, https://hal.science/tel-02588353).
#' \cr
#' The uncertainty in tau is still estimated by plugging-in point-estimates of mu1, mu2 and sigma
#' for each possible value of tau, and computing the associated likelihood. This likelihood can then be normalized to
#' estimate the posterior pdf of tau, which can then be numerically integrated to give the posterior cdf of tau.
#' Samples of tau values can then be simulated from this cdf using uniform sampling then inverse-cdf transformation.
#' The resulting algorithm is MCMC-free and independent of RBaM.
#'
#' @param obs real vector, observations
#' @param time vector, time in POSIXct, string or numeric format
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param nS integer, number of segments, either 1 or 2
#' @param nMin integer, minimum number of observations by segment
#' @param nSim integer, number of simulated tau values
#' @param varShift logical, allow for a shifting variance?
#' @param alpha real in (0;1), type-I error level of the underlying step-change test.
#' @return list with the following components:
#' \enumerate{
#'   \item summary: list, summarize the information to present to the user
#'   \itemize{
#'        \item data: data frame, data augmented with a column denoting the period after segmentation
#'        \item shift: data frame, detected shift time in numeric or POSIXct format in UTC
#'   }
#'   \item plot : list, data formatted to use as input for some plot functions
#'   \itemize{
#'        \item density.tau: data frame, a table with three columns. The first column indicates the specific shift being analyzed.
#'        The second column contains the value of the shift time tau. The last columns shows the probability density associated
#'        with each tau value
#'        \item density.inc.tau: data frame, all information about the 95% credibility interval and the Maximum a posterior (MAP) estimation
#'        for shift times with their associated probability densities
#'   }
#'   \item tau: real, estimated shift time in numeric or POSIXct format in UTC
#'   \item segments: list, segment mean value indexed by the list number
#'   \item mcmc: data frame, Monte-Carlo simulations. Note that values for mu's and sigma's are fixed to the values corresponding
#'       to the maxpost tau.
#'   \item data.p: list, separate and assign information by identified stable period indexed by the list number
#'   \item DIC: real, pseudo-DIC (see description)
#'   \item origin.date.p: positive real or date, date describing origin of the segmentation for a sample. Useful for recursive segmentation.
#' }
#' @examples
#' # Segmentation into two segments for the RhoneRiverAMAX data set (details in ?RhoneRiverAMAX)
#' res=Segmentation_quickApprox(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH,nS=2)
#' res$summary$shift
#' @export
#' @importFrom stats quantile sd approxfun runif optim
Segmentation_quickApprox <- function(obs,time=1:length(obs),u=0*obs,nS=2,
                                     nMin=3,nSim=500,varShift=FALSE,alpha=0.1){
  n=NROW(obs)
  cll=rep(-Inf,n) # conditional log-likelihoods (for each tau)
  foo=sort.int(time,index.return=TRUE) # make sure everything is sorted chronologically
  obs=obs[foo$ix];time=time[foo$ix];u=u[foo$ix]
  ix=nMin:(n-nMin) # indices of candidate taus
  if(any(diff(ix)<0)){
    mess=paste0('Not enough data (n=',n,') to detect a shift with nMin=',nMin,
                ' values on each side of the shift')
    stop(mess)
  }
  if(nMin<3){stop('nMin cannot be smaller than 3')}
  if(nS>2){stop('Quick-approximation algorithm is only available for nS<3 (i.e. single-shift model at most)')}

  # Start preparing output list
  out=getOutputList_Engine(time,obs,u)
  # no-change model
  start=c(mean(obs),sd(obs))
  if(all(u==0)){ # no need to optimize, estimate is explicit
    theta0=start
    f0=quickApprox_llfunk0(start,obs,u)
  } else { # maximize log-lkh
    opt=optim(start,quickApprox_llfunk0,obs=obs,u=u,control=list(fnscale=-1))
    theta0=opt$par
    f0=opt$value
  }
  if(nS==1){ # stop here and return
    out$segments=rep(theta0[1],n)
    out$mcmc=data.frame(mu1=rep(theta0[1],nSim),
                        structural_sd=rep(theta0[2],nSim),
                        LogPost=rep(f0,nSim))
    out$data.p=list(obs.p=obs,time.p=time,u.p=u)
    out$DIC=-2*f0
    return(out)
  }

  # Single-shift model
  if(varShift){ # mu1,mu2,sig1,sig2
    start=c(theta0[1],theta0[1],theta0[2],theta0[2])
  } else {  # mu1,mu2,sig
    start=c(theta0[1],theta0[1],theta0[2])
  }
  taus=time[ix]
  pars=vector('list',length(time))
  for (i in 1:length(taus)){ # for each tau, find mus/sds by max. likelihood
    foo=quickApprox_getTheta(tau=taus[i],time=time,obs=obs,u=u,varShift=varShift)
    pars[[ix[i]]]=foo$par
    cll[ix[i]]=foo$value
  }
  # Approximate cdf by rectangle integration
  post=exp(cll-max(cll)) # avoid numerical zeros, remains proportional to pdf
  w=c(diff(time),0) # basis of each rectangle
  cdf0=c(0,cumsum(w*post))[1:n] # unnormalized cdf
  cdf=cdf0/max(cdf0) # normalized cdf
  post=post/max(cdf0) # normalized pdf
  # compute pseudo-DIC (see Description)
  d=ifelse(varShift,2,1)
  crit=getCriticalValue(d=d,alpha=alpha,n=n)
  DIC1=-2*max(cll)+crit
  # simulate taus from posterior
  unif=runif(nSim)
  foo=approxfun(x=cdf,y=time,ties='ordered')
  sim=sapply(unif,foo)
  # ready to return
  imax=which.max(post)
  tau=time[imax]
  out$summary$data$period[time>tau]=2
  n1=sum(out$summary$data$period==1)
  n2=sum(out$summary$data$period==2)
  out$summary$shift[1,]=c(tau,quantile(sim,probs=c(0.025,0.975)))
  out$plot$density.tau=data.frame(Shift='tau1',Value=time,Density=post)
  foo=approxfun(x=time,y=post)
  out$plot$density.inc.tau=data.frame(Shift='tau1',
                                      tau_lower_inc=out$summary$shift$I95_lower,
                                      density_tau_lower_inc=foo(out$summary$shift$I95_lower),
                                      tau_upper_inc=out$summary$shift$I95_upper,
                                      density_tau_upper_inc=foo(out$summary$shift$I95_upper),
                                      taU_MAP=tau,
                                      density_taU_MAP=foo(tau))
  out$tau=tau
  out$segments=list(rep(pars[[imax]][1],n1),rep(pars[[imax]][2],n2))
  out$mcmc=data.frame(mu1=pars[[imax]][1],mu2=pars[[imax]][2],tau1=sim)
  if(varShift){
    out$mcmc$sig1=pars[[imax]][3]
    out$mcmc$sig2=pars[[imax]][4]
  } else {
    out$mcmc$structural_sd=pars[[imax]][3]
  }
  foo=approxfun(x=time,y=post)
  out$mcmc$LogPost=foo(sim)
  mask=(time<=tau)
  out$data.p$obs.p=list(obs[mask],obs[!mask])
  out$data.p$time.p=list(time[mask],time[!mask])
  out$data.p$u.p=list(u[mask],u[!mask])
  out$DIC=DIC1
  return(out)
}


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
  out$summary$data$time = numeric_to_time(d=out$summary$data$time,origin.date = origin.date)
  # time series indexed by number of segments identified
  if(is.list(out$data.p$time.p)==T){
    for(i in 1:length(out$data.p$time.p)){
      out$data.p$time.p[[i]] =  numeric_to_time(d=out$data.p$time.p[[i]],origin.date = origin.date)
    }
  }else{
    out$data.p$time.p = numeric_to_time(d=out$data.p$time.p,origin.date=origin.date)
  }
  # Rating shift time summary
  if(all(out$summary$shift$tau!= 0)){
    # Transform all time in POSIXct format
    out$summary$shift <- data.frame(lapply(out$summary$shift, function(column) {
      numeric_to_time(d=column,origin.date=origin.date)
    }))
  }
  # Tau in input date format
  out$tau=out$summary$shift$tau
  # Transform density data in POSIXct format if necessary
  if(NROW(out$summary$shift)>0){
    out$plot$density.tau$Value <- numeric_to_time(d=out$plot$density.tau$Value,origin.date=origin.date)
    out$plot$density.inc.tau$tau_lower_inc <- numeric_to_time(d=out$plot$density.inc.tau$tau_lower_inc,origin.date = origin.date)
    out$plot$density.inc.tau$tau_upper_inc <- numeric_to_time(d=out$plot$density.inc.tau$tau_upper_inc,origin.date = origin.date)
    out$plot$density.inc.tau$taU_MAP <- numeric_to_time(d=out$plot$density.inc.tau$taU_MAP,origin.date=origin.date)
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
