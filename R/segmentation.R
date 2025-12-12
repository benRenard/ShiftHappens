#' Recursive Segmentation
#'
#' Recursive segmentation procedure for an \strong{unknown} number of segments
#'
#' @inheritParams Segmentation
#'
#' @inherit recursiveSegmentation return
#'
#' @examples
#' obs=RhoneRiverAMAX$H
#' res=Segmentation_Recursive(obs=obs)
#' res$shifts
#' res$tree
#' plot(res$data$obs,col=res$data$period)
#' # Create an artificial series with more changes to see recursion in action
#' obs=c(RhoneRiverAMAX$H,RhoneRiverAMAX$H)
#' res=Segmentation_Recursive(obs=obs)
#' res$shifts
#' res$tree
#' plot(res$data$obs,col=res$data$period)
#' @export
Segmentation_Recursive <- function(obs,
                                   time=1:length(obs),
                                   u=0*obs,
                                   nSmax=2,
                                   doQuickApprox=TRUE,
                                   nMin=ifelse(doQuickApprox,3,1),
                                   nSim=500,varShift=FALSE,alpha=0.1,
                                   mcmc_options=RBaM::mcmcOptions(),
                                   mcmc_cooking=RBaM::mcmcCooking(),
                                   temp.folder=file.path(tempdir(),'BaM'),
                                   mu_prior = list(),...){
  if(length(obs)<2){
    stop('At least 2 observations are required',call.=FALSE)
  }
  # Initialization
  allRes=list() # store segmentation results for all nodes in a sequential list
  k=0 # Main counter used to control indices in allRes
  tree=data.frame() # store tree structure (parents - children relationship)
  p=1 # Auxiliary counter needed to keep track of children / parents indices
  level=0 # Recursion level. The tree is created level-by-level rather than branch-by-branch
  X=list(obs) # List of all nodes (each corresponding to a subseries of x) to be segmented at this level. Start with a unique node corresponding to the whole series
  TIME=list(time) # List of corresponding times
  U=list(u) # List of corresponding uncertainties
  indices=c(1) # Vector containing the indices of each node - same size as X
  parents=c(0) # Vector containing the indices of the parents of each node - same size as X
  continue=TRUE # Logical determining whether recursion should continue

  while(continue){
    level=level+1 # Increment recursion level
    nX=length(X) # Number of nodes at this level
    keepgoing=rep(NA,nX) # Should recursion continue for each node?
    newX=newTIME=newU=newIndices=newParents=c() # Will be used to update subseries, indices and parents at the end of each recursion level
    m=0 # Local counter used to control indices in the 4 vectors above => reset to 0 at each new level of the recursion
    for(j in 1:nX){ # Loop on each node
      k=k+1 # Increment main counter
      if(NROW(X[[j]])<nSmax*nMin){ # Can't segment, return default result
        partial.segmentation=simpleSegmentation(time=TIME[[j]],obs=X[[j]],u=U[[j]])
      } else { # Apply segmentation to subseries stored in node X[[j]]
        partial.segmentation=Segmentation(obs=X[[j]],time=TIME[[j]],u=U[[j]],
                                          nSmax=nSmax,doQuickApprox=doQuickApprox,nMin=nMin,
                                          nSim=nSim,varShift=varShift,alpha=alpha,
                                          mcmc_options=mcmc_options,mcmc_cooking=mcmc_cooking,
                                          temp.folder,mu_prior=mu_prior) #,...)
      }
      # Save results for this node
      allRes[[k]]=partial.segmentation
      # Save optimal number of segments
      nSopt=partial.segmentation$nS
      # Update recursion tree
      tree=rbind(tree,data.frame(indx=k,level=level,parent=parents[j],nS=nSopt))
      # This was the trickiest part: keeping track of indices and parents
      keepgoing[j]=nSopt>1 # if nS=1, segmentation will not continue for this node which is hence terminal
      if(keepgoing[j]){ # Save results for segmentation at next level
        for(i in 1:nSopt){ # Loop on each segment detected for the current node
          p=p+1 # Increment auxiliary counter
          m=m+1 # Increment local counter
          mask=partial.segmentation$results[[nSopt]]$data$period==i
          newX[[m]]=partial.segmentation$results[[nSopt]]$data$obs[mask] # Save ith segment (on a total of nS)
          newTIME[[m]]=partial.segmentation$results[[nSopt]]$data$time[mask] # Save corresponding times
          newU[[m]]=partial.segmentation$results[[nSopt]]$data$u[mask] # Save corresponding uncertainty
          newParents[m]=indices[j] # At next level, the parent of this segment will be the index of current node
          newIndices[m]=p # At next level, the index of this segment will be p
        }
      }
    }
    # Check if recursion should continue at all, i.e. if at least one node is not terminal
    if(all(keepgoing==FALSE)) continue=FALSE
    # Update list of nodes to be further segmented at next level + parents and indices
    X=newX
    TIME=newTIME
    U=newU
    parents=newParents
    indices=newIndices
  }
  # Get terminal nodes
  terminal=which(tree$nS==1)
  # Assemble final dataset with period column by using information from terminal nodes
  data <- c()
  for(i in 1:length(terminal)){
    data.stable.p=allRes[[terminal[i]]]$results[[1]]$data #Save data from stable period
    node = data.frame(time=data.stable.p$time,obs=data.stable.p$obs,u=data.stable.p$u,
                      I95_lower=data.stable.p$obs+stats::qnorm(0.025)*data.stable.p$u,
                      I95_upper=data.stable.p$obs+stats::qnorm(0.975)*data.stable.p$u,
                      period = rep(i,NROW(data.stable.p)))
    data = rbind(data,node)
  }
  # Get info from nodes with shifts
  hasShift=which(tree$nS!=1)
  if(length(hasShift)>0){
    shift <- c()
    for(i in 1:length(hasShift)){
      nSopt.p = allRes[[hasShift[i]]]$nS
      results.p = allRes[[hasShift[i]]]$results
      shift.time.p=data.frame(c(results.p[[nSopt.p]]$shifts$tau))
      for(j in 1:(nSopt.p-1)){
        shift.time.p.unc=data.frame(tau=shift.time.p[j,],
                                    I95_lower=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],probs=c(0.025)),
                                    I95_upper=stats::quantile(results.p[[nSopt.p]]$mcmc[,nSopt.p+j],probs=c(0.975)),
                                    id_iteration=hasShift[[i]])
        shift <- rbind(shift,shift.time.p.unc)
      }
    }
    rownames(shift) <- NULL
    # Transform uncertainty on the shift in POSIXct format
    origin.date=allRes[[1]]$origin.date
    if(all(!is.numeric(time))){
      transformed.shift <- c()
      for(i in 1:NROW(shift)){
        transformed.shift.p <- data.frame(lapply(shift[i,1:3], function(column) {
          numeric_to_time(d=column,origin.date=origin.date)
        }))
        transformed.shift=rbind(transformed.shift,transformed.shift.p)
      }
      shift[,1:3]=transformed.shift
    }
    shift <- shift[order(shift$tau),]
  } else {
    shift=NULL
  }
  # Tidy up returned data
  data=data[order(data$time),] # sort
  data$period=c(1,1+cumsum(diff(data$period)!=0))
  # 2DO: review return object
  out=recursiveSegmentation(data=data,shifts=shift,tree=tree,
                            origin.date=origin.date,results=allRes)
  return(out)
}

#' Segmentation
#'
#' Segmentation procedure for an \strong{unknown} number of segments
#'
#' @inheritParams Segmentation_Engine
#' @param nSmax integer, maximum number of segments to assess
#'
#' @inherit multipleSegmentation return
#'
#' @examples
#' res=Segmentation(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Date,u=RhoneRiverAMAX$uH)
#' res$shifts
#' res$DICs
#' hist(res$mcmc$tau)
#' # Segmentation using BaM, allowing more than 2 segments
#' \dontrun{
#' # This line of code is wrapped in \dontrun{} since it relies
#' # on the installation of an executable with the package RBaM.
#' # See ?RBaM::downloadBaM for downloading and installing this executable.
#' res=Segmentation(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH,
#'      doQuickApprox=FALSE,nSmax=3)
#' }
#' @export
Segmentation <- function(obs,
                         time=1:length(obs),
                         u=0*obs,
                         nSmax=2,
                         doQuickApprox=TRUE,
                         nMin=ifelse(doQuickApprox,3,1),
                         nSim=500,varShift=FALSE,alpha=0.1,
                         mcmc_options=RBaM::mcmcOptions(),
                         mcmc_cooking=RBaM::mcmcCooking(),
                         temp.folder=file.path(tempdir(),'BaM'),
                         mu_prior = list(),...){

  if(nSmax<=0){
    stop('Maximum number of segments should be larger than 0',call.=FALSE)
  }
  if(length(obs)<2){
    stop('At least 2 observations are required',call.=FALSE)
  }

  res=vector(mode = 'list',length = nSmax)
  DICs <- rep(NA,nSmax)
  for(i in (1:nSmax)){
    nS <- i
    quick=ifelse(nS<3,doQuickApprox,FALSE)
    if(length(obs)<nS){
      warning(paste0('NA was returned because the number of observations (',length(obs),
                     ') is lower than the number of segments (',nS,')'))
      DICs [i] <- NA
    }else if(trunc(length(obs)/nS)<nMin){
      warning(paste0('The minimum number of observations per segment (',nMin,') cannot be matched with the number of observations (',length(obs),
                     ') and the number of segments (',nS,')'))
      DICs [i] <- NA
    }else{
      # Check if prior knowledge has been provided:
      if(length(mu_prior)==0){
        mu_args <- list(NULL)
      }else{
        mu_args <- mu_prior
      }
      res[[i]] <- Segmentation_Engine(obs=obs,time=time,u=u,nS=nS,nMin=nMin,
                                      mcmc_options=mcmc_options,mcmc_cooking=mcmc_cooking,
                                      temp.folder=temp.folder,mu_prior=mu_args,
                                      doQuickApprox=quick,varShift=varShift,alpha=alpha,...)
      DICs [i] <- res[[i]]$DIC
    }
  }
  out=multipleSegmentation(res)
  return(out)
}

#' Segmentation engine
#'
#' Segmentation procedure for a \strong{known} given number of segments
#'
#' @param obs real vector, observations
#' @param time vector, time in POSIXct, string or numeric format
#' @param u real vector, uncertainty in observations (as a standard deviation)
#' @param nS integer, number of segments
#' @param doQuickApprox logical, use quick approximation? see ?Segmentation_quickApprox
#' @param nMin integer, minimum number of observations by segment
#' @param nSim integer, number of Monte-Carlo simulated values for shift times. Only used when doQuickApprox=TRUE.
#' @param varShift logical, allow for a shifting variance? Only used when doQuickApprox=TRUE.
#' @param alpha real in (0;1), type-I error level of the underlying step-change test. Only used when doQuickApprox=TRUE.
#' @param mcmc_options mcmcOptions object, options used to generate MCMC samples (e.g. size and adaption tunings), see ?RBaM::mcmcOptions.
#'     Only used when doQuickApprox=FALSE.
#' @param mcmc_cooking mcmcCooking object, options used to post-process MCMC samples (e.g. burn and slice), see ?RBaM::mcmcCooking.
#'     Only used when doQuickApprox=FALSE.
#' @param temp.folder directory, temporary directory to write BaM computations.
#'     Only used when doQuickApprox=FALSE.
#' @param mu_prior list, object describing prior knowledge on the mean of residuals for each segment.
#'     Only used when doQuickApprox=FALSE.
#' @param ... other arguments passed to RBaM::BaM.
#'
#' @inherit simpleSegmentation return
#'
#' @source \url{https://theses.hal.science/tel-03211343}
#' @source \url{https://agupubs.onlinelibrary.wiley.com/doi/epdf/10.1029/2020WR028607}
#'
#' @examples
#' # Default segmentation using quick approximation
#' res=Segmentation_Engine(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH)
#' res$shifts
#' hist(res$mcmc$tau)
#' # Segmentation using BaM, allowing more than 2 segments and the use of priors
#' \dontrun{
#' # This line of code is wrapped in \dontrun{} since it relies
#' # on the installation of an executable with the package RBaM.
#' # See ?RBaM::downloadBaM for downloading and installing this executable.
#' prior=list(RBaM::parameter(name=paste0('mu1'),init=6,prior.dist='Gaussian',prior.par=c(6,5)),
#'            RBaM::parameter(name=paste0('mu2'),init=8,prior.dist='Gaussian',prior.par=c(8,2))))
#' res=Segmentation_Engine(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH,
#'      doQuickApprox=FALSE,nS=3,mu_prior=prior)
#' }
#' @export
#' @importFrom RBaM parameter xtraModelInfo model dataset mcmcOptions mcmcCooking remnantErrorModel BaM
#' @importFrom stats quantile sd
#' @importFrom utils read.table
Segmentation_Engine <- function(obs,
                                time=1:length(obs),
                                u=0*obs,
                                nS=2,
                                doQuickApprox=TRUE,
                                nMin=ifelse(doQuickApprox,3,1),
                                nSim=500,varShift=FALSE,alpha=0.1,
                                mcmc_options=RBaM::mcmcOptions(),
                                mcmc_cooking=RBaM::mcmcCooking(),
                                temp.folder=file.path(tempdir(),'BaM'),
                                mu_prior=list(NULL),...){

  if(length(obs)<2)stop('At least 2 observations are required',call.=FALSE)
  if(length(obs)<nS)stop('Number of observations is lower than the number of segments',call.=FALSE)
  if(any(is.na(obs)) | any(is.na(time)) | any(is.na(u)))stop('Missing values not allowed in observation, time and uncertainty')

  if(nS<=0)stop('Maximum number of segments should be larger than 0',call.=FALSE)
  if(trunc(length(obs)/nS)<nMin)stop(paste0('The minimum number of observations per segment (',nMin,') cannot be matched with the number of observations (',length(obs),
                                            ') and the number of segments (',nS,')'))
  if(is.null(check_equal_length(obs,time,u)))stop('The observations, time and uncertainty have not the same length')

  # Sort data frame case time not ascending
  DF.order <- data.frame(obs=obs,time=time,u=u)
  DF.order <- DF.order[order(DF.order$time),]

  # Check time format
  numeric.check=TRUE

  # Get origin date of the segmentation
  origin.date <- min(DF.order$time)

  # Date transformation function to passe to numeric format if necessary
  if(!is.numeric(DF.order$time)){
    DateTransformed <- time_to_numeric(date=DF.order$time)
    DF.order$time <- DateTransformed$d
    origin.date <- DateTransformed$origin.date
    numeric.check=FALSE
  }

  obs <- DF.order$obs
  time <- DF.order$time
  u <- DF.order$u

  if(doQuickApprox){
    hasPrior=!sapply(mu_prior,is.null)
    if(hasPrior){
      warning('Quick approximation does not handle prior information. The provided priors will be ignored.')
    }
    out=Segmentation_quickApprox(obs=obs,time=time,u=u,nS=nS,nMin=nMin,
                                 nSim=nSim,varShift=varShift,alpha=alpha)
  } else {
    npar = nS + nS - 1

    priors <- vector(mode = 'list',length = npar)
    # Extract mu parameter from if user-defined
    if(is.null(mu_prior[[1]])){
      # Set default for mu_prior if no mu_prior are provided
      for(i in 1:nS){
        priors [[i]] <- RBaM::parameter(name=paste0('mu',i),
                                        init=mean(obs),
                                        prior.dist = 'FlatPrior' ,
                                        prior.par = NULL)
      }
    }else{
      # Use the provided mu parameter from mu_prior
      for(i in 1:nS){
        mu_list <- mu_prior[[1]]
        mu_list$name <- paste0(mu_list$name,i)
        priors [[i]] <- mu_list
      }
    }

    if(i>1){
      # Set default for tau_prior
      prior_tau_init <- as.numeric(stats::quantile(time,probs = seq(1,nS-1)/nS))

      for(i in 1:(nS-1)){
        priors [[nS+i]] <- RBaM::parameter(name=paste0('tau',i),
                                           init= prior_tau_init[i],
                                           prior.dist = 'FlatPrior' ,
                                           prior.par = NULL)
      }
    }
    # Config_Xtra
    xtra=RBaM::xtraModelInfo(object=c(nS=nS,tmin_xtra=0,nmin_xtra=nMin,option_xtra=1))
    # Model
    mod=RBaM::model(
      ID='Segmentation',
      nX=1,
      nY=1,
      par=priors,
      xtra=xtra)

    # dataset object
    data=RBaM::dataset(X=data.frame(time),
                       Y=data.frame(obs),
                       Yu=data.frame(u),
                       data.dir=temp.folder)

    remnantInit=stats::sd(obs)
    if(is.na(remnantInit)){ # happens when nObs=1
      remnantInit=abs(mean(obs))
    }
    if(remnantInit==0){remnantInit=1}
    remnant_prior <- list(RBaM::remnantErrorModel(funk = "Constant",
                                                  par = list(RBaM::parameter(name="gamma1",
                                                                             init=remnantInit,
                                                                             prior.dist = "FlatPrior+"))))

    # Run BaM executable
    RBaM::BaM(mod=mod,
              data=data,
              workspace=temp.folder,
              mcmc=mcmc_options,
              cook = mcmc_cooking,
              remnant = remnant_prior,...)

    mcmc.segm    <- utils::read.table(file=file.path(temp.folder,mcmc_cooking$result.fname),header=TRUE)
    mcmc.DIC     <- utils::read.table(file=file.path(temp.folder,"Results_DIC.txt"),header=FALSE)
    resid.segm   <- utils::read.table(file=file.path(temp.folder,"Results_Residuals.txt"),header=TRUE)

    colnames(mcmc.segm)[ncol(mcmc.segm)-1] <- "structural_sd"

    simulation.MAP <- resid.segm$Y1_sim

    data = data.frame(time=time,obs=obs,u=u,
                      I95_lower=obs+stats::qnorm(0.025)*u,I95_upper=obs+stats::qnorm(0.975)*u,
                      period = 1)

    if(nS==1){
      shift=data.frame(tau=numeric(0), # no shift time
                       I95_lower=numeric(0),
                       I95_upper=numeric(0))
    } else {
      shift <- c()
      for(j in 1:(nS-1)){
        foo=data.frame(tau=mcmc.segm[which.max(mcmc.segm$LogPost),(nS+j)],
                                    I95_lower=stats::quantile(mcmc.segm[,(nS+j)],probs=c(0.025)),
                                    I95_upper=stats::quantile(mcmc.segm[,(nS+j)],probs=c(0.975)))
        shift <- rbind(shift,foo)
      }
      rownames(shift) <- NULL
      shift <- shift[order(shift$tau),]

      # Store sub series into a list
      obss=periods=vector(mode='list',length=nS)
      intervals.time.shift=c(-Inf,shift$tau,Inf) # intervals defined by time shifts
      for(i in 1:nS){
        position.ti.p <- which((time-intervals.time.shift[[i]])>=0)[1]
        position.tf.p <- rev(which((time-intervals.time.shift[[i+1]])<0))[1]
        obss[[i]]=obs[position.ti.p:position.tf.p]
        periods[[i]]=rep(i,length(obss[[i]]))
      }
      data$period = unlist(periods)
    }
    # Assemble output object
    out=simpleSegmentation(obs=obs,time=time,u=u,
                           data=data,shifts=shift,mcmc=mcmc.segm,
                           DIC=mcmc.DIC[1,2],origin.date=origin.date)
  }
  # Transform back all time into original units
  if(!numeric.check){out=transformDatesInOutput(out,origin.date)}
  return(out)
}

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
#'
#' @inherit simpleSegmentation return
#'
#' @examples
#' # Segmentation into two segments for the RhoneRiverAMAX data set (details in ?RhoneRiverAMAX)
#' res=Segmentation_quickApprox(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Year,u=RhoneRiverAMAX$uH,nS=2)
#' res$shifts
#' hist(res$mcmc$tau)
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
  out=simpleSegmentation(obs=obs,time=time,u=u)
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
    out$mcmc=data.frame(mu1=rep(theta0[1],nSim),
                        structural_sd=rep(theta0[2],nSim),
                        LogPost=rep(f0,nSim))
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
  out$data$period[time>tau]=2
  out$shifts[1,]=c(tau,quantile(sim,probs=c(0.025,0.975)))
  out$mcmc=data.frame(mu1=pars[[imax]][1],mu2=pars[[imax]][2],tau1=sim)
  if(varShift){
    out$mcmc$sig1=pars[[imax]][3]
    out$mcmc$sig2=pars[[imax]][4]
  } else {
    out$mcmc$structural_sd=pars[[imax]][3]
  }
  foo=approxfun(x=time,y=post)
  out$mcmc$LogPost=foo(sim)
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
  out$data$time = numeric_to_time(d=out$data$time,origin.date = origin.date)
  # Rating shift time summary
  if(all(out$shifts$tau!= 0)){
    # Transform all time in POSIXct format
    out$shifts <- data.frame(lapply(out$shifts, function(column) {
      numeric_to_time(d=column,origin.date=origin.date)
    }))
  }
  out$origin.date=origin.date
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
