X=read.table('data-raw/ArdecheRiverGaugings.txt',header=T,sep='\t')
# Remove last two columns
X=X[,1:(NCOL(X)-2)]
DateString=paste(paste(X$Day,X$Month,X$Year,sep = "/"),paste(X$Hour,X$Minutes,'00',sep = ":"),sep=" ")

ArdecheRiverGaugings=data.frame(Date=as.POSIXct(DateString,format='%d/%m/%Y %H:%M:%S',tz='UTC'),
                                H=X$h,Q=X$Q,
                                uQ=(X$Q*0.01*X$Qsigma)/1.96)

save(ArdecheRiverGaugings,file='data/ArdecheRiverGaugings.RData')
