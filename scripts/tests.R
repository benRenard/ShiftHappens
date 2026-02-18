library(ShiftHappens)
library(ggplot2)
library(RBaM)
# default
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$H)
plot(sg)

# use dates
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Date)
plot(sg)

# use uncertainties
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Date,u=RhoneRiverAMAX$uH)
plot(sg)

# no shift
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$Q,time=RhoneRiverAMAX$Date,u=RhoneRiverAMAX$uQ)
plot(sg)

# many shifts, non-recursive
sg=Segmentation(obs=RhoneRiverAMAX$uQ,time=RhoneRiverAMAX$Date,nSmax=7)
if(!is.null(sg)){
  plot(sg$DICs)
  matplot(sg$mcmc[,grep('tau',names(sg$mcmc))])
  plot(sg)
}

# many shifts, recursive
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

# Recursive modeling with default linear regression
sg=Segmentation_RecursiveModeling(x=ArdecheRiverGaugings$H,
                                  y=ArdecheRiverGaugings$Q,
                                  uY=ArdecheRiverGaugings$uQ,
                                  time=ArdecheRiverGaugings$Date)
sg$segmentation$shifts
plot(sg)
plot(sg,dataPlotType='ty')
plot(sg$segmentation)
plot(sg$segmentation$tree)

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
