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
