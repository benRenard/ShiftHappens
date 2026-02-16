X=read.table('data-raw/ArdecheRiverStage.txt',header=T,sep='\t')
colnames(X) <- c('NumericDate','H','Date','DateFromZero')

DateString <- paste(X$Date,rep('00',length(X$Date)),sep = ':')
ArdecheRiverStage <- data.frame(Date = as.POSIXct(DateString, format = '%d/%m/%Y %H:%M:%S', tz='UTC'),
                                H = X$H/100) # Stage record in meter
save(ArdecheRiverStage,file='data/ArdecheRiverStage.RData')
