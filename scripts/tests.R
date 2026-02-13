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
plot(sg$DICs)
matplot(sg$mcmc[,grep('tau',names(sg$mcmc))])
plot(sg)

# many shifts, recursive
sg=Segmentation_Recursive(obs=RhoneRiverAMAX$uQ,time=RhoneRiverAMAX$Date)
plot(sg)
PlotTree(sg$tree)
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
                                  uY=ArdecheRiverGaugings$uQ)
sg$segmentation$shifts
plot(sg$segmentation)
