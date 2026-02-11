# default
x=Segmentation_Recursive(obs=RhoneRiverAMAX$H)
plot(x)

# use dates
x=Segmentation_Recursive(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Date)
plot(x)

# use uncertainties
x=Segmentation_Recursive(obs=RhoneRiverAMAX$H,time=RhoneRiverAMAX$Date,u=RhoneRiverAMAX$uH)
plot(x)

# no shift
x=Segmentation_Recursive(obs=RhoneRiverAMAX$Q,time=RhoneRiverAMAX$Date,u=RhoneRiverAMAX$uQ)
plot(x)

# many shifts, non-recursive
x=Segmentation(obs=RhoneRiverAMAX$uQ,time=RhoneRiverAMAX$Date,nSmax=7)
plot(x$DICs)
matplot(x$mcmc[,grep('tau',names(x$mcmc))])
plot(x)

# many shifts, recursive
x=Segmentation_Recursive(obs=RhoneRiverAMAX$uQ,time=RhoneRiverAMAX$Date)
plot(x)
PlotTree(x$tree)
shifts=getShifts(x)
for(i in 1:NCOL(shifts)){
  if(i==1){
    plot(shifts[,i],type='l',ylim=range(x$data$time),col=i)
  } else {
    lines(shifts[,i],col=i)
  }
}
