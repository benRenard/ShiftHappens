library(ShiftHappens)
library(ggplot2)
library(patchwork)
library(RBaM)

# Recursive segmentation using default
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$H)
plot(sg)

# Use explicit dates
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Date)
plot(sg)

# Use unsorted dates
ix=sample(1:NROW(RhoneRiverAMAX))
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$H[ix],time=RhoneRiverAMAX$Date[ix])
plot(sg)

# Use uncertainties
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Date,u=RhoneRiverAMAX$uH)
plot(sg)

# No-shift case
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$Q,time=RhoneRiverAMAX$Date,u=RhoneRiverAMAX$uQ)
plot(sg)

# Many shifts (nSmax=7), non-recursive, using BaM
sg=Segmentation(obs=RhoneRiverAMAX$uQ,time=RhoneRiverAMAX$Date,nSmax=7)
if(!is.null(sg)){
  plot(sg$DICs)
  matplot(sg$mcmc[,grep('tau',names(sg$mcmc))])
  plot(sg)
}

# Many shifts, recursive (does not use BaM)
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$uQ,time=RhoneRiverAMAX$Date)
plot(sg)
plot(sg$tree)
shifts=getShifts(sg)
for(i in 1:NCOL(shifts)){
  if(i==1){
    plot(shifts[,i],type='l',ylim=range(sg$data$time),col=i)
  } else {
    lines(shifts[,i],col=i)
  }
}
# Unsorted dates
ix=sample(1:NROW(RhoneRiverAMAX))
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$uQ[ix],time=RhoneRiverAMAX$Date[ix])
plot(sg)
sg$shifts

# Recursive modeling with default linear regression
sg=Segmentation_RecursiveModeling(x=ArdecheRiverGaugings$H,
                                  y=ArdecheRiverGaugings$Q,
                                  uY=ArdecheRiverGaugings$uQ,
                                  time=ArdecheRiverGaugings$Date)
sg$segmentation$shifts
plot(sg) # default 'xy' scatterplot for recursiveModeling object
plot(sg,dataPlotType='ty') # 'time-y' alternative plot in the data panel
wrap_plots(plot(sg,type='fits')) # plot fitted model in each terminal node
plot(sg$segmentation) # plot segmentation object contained in sg
plot(sg$segmentation$tree) # plot segmentation tree contained in sg
# Unsorted dates
ix=sample(1:NROW(ArdecheRiverGaugings))
sg=Segmentation_RecursiveModeling(x=ArdecheRiverGaugings$H[ix],
                                  y=ArdecheRiverGaugings$Q[ix],
                                  uY=ArdecheRiverGaugings$uQ[ix],
                                  time=ArdecheRiverGaugings$Date[ix])
sg$segmentation$shifts
plot(sg) # default 'xy' scatterplot for recursiveModeling object

# Recursive modeling with BaRatin
controlMatrix=rbind(c(1,0,0),c(0,1,0),c(0,1,1))
k1=RBaM::parameter(name='k1',init=-0.6,prior.dist='Gaussian',prior.par=c(-0.6,0.5))
a1=RBaM::parameter(name='a1',init=exp(2.65),prior.dist='LogNormal',prior.par=c(2.65,0.35))
c1=RBaM::parameter(name='c1',init=1.5,prior.dist='Gaussian',prior.par=c(1.5,0.025))
k2=RBaM::parameter(name='k2',init=0,prior.dist='Gaussian',prior.par=c(0,1))
a2=RBaM::parameter(name='a2',init=exp(3.28),prior.dist='LogNormal',prior.par=c(3.28,0.33))
c2=RBaM::parameter(name='c2',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))
k3=RBaM::parameter(name='k3',init=1.2,prior.dist='Gaussian',prior.par=c(1.2,0.2))
a3=RBaM::parameter(name='a3',init=exp(3.46),prior.dist='LogNormal',prior.par=c(3.46,0.38))
c3=RBaM::parameter(name='c3',init=1.67,prior.dist='Gaussian',prior.par=c(1.67,0.025))
priors=list(k1,a1,c1,k2,a2,c2,k3,a3,c3)
# Note: possible to use wrapper function ?Segmentation_BaRatin instead of the function below
sg=Segmentation_RecursiveModeling(x=ArdecheRiverGaugings$H,
                                  y=ArdecheRiverGaugings$Q,
                                  uY=ArdecheRiverGaugings$uQ,
                                  time=ArdecheRiverGaugings$Date,
                                  Fit_funk=Fit_BaRatin,
                                  controlMatrix=controlMatrix,priors=priors)
if(!is.null(sg)){
  plot(sg)
  plot(sg,type='data',dataPlotType='xy')+scale_y_continuous(trans='log')
  plot(sg,dataPlotType='ty')
  plot(sg$segmentation$tree)
}

# Recession tools
rec=Extract_Recessions(time=ArdecheRiverStage$Date,
                       H=ArdecheRiverStage$H,
                       uH=rep(0.1,NROW(ArdecheRiverStage)))
plot(rec)
plot(rec,type='th')
plot(rec,type='dhmin')
# Segmentation based on recession minima
lows=getRecessionMin(rec)
sg=Segmentation_Recursive(obs=lows$H,time=lows$date)
plot(sg)
# Model recessions and segment asymptotic stages
f=Fit_Recessions(rec,equation='M7',
                 mcmc_options=mcmcOptions(nAdapt=20,nCycles=10)) # speed-up
betas=f$parameters[1:max(rec$index),]
sg=Segmentation_Recursive(obs=betas$value,u=betas$u)
plot(sg)
# Possible to do segmentation directly from rec using wrapper function Segmentation_Recession:
sg=Segmentation_Recessions(rec) # by default, uses recession minima with no modeling
plot(sg)
