load("~/Documents/Data/BiasCorrection/CollectedSPSResult_TrustSharpTrue_corrected")

Cols=c("#1b9e77","#d95f02","#7570b3","#e7298a")
par(mar=c(5,5,4,1)+.1) #Add this before calling plot to add margin for ylabel
plot(0,0,xlim=c(1,12),ylim=c(0,8),xlab="Leadtime (months)",ylab=expression(paste("SPS (10"^"6"," km"^"2",")")),pch='.',yaxs="i",main="Calibrated and raw SSIPS forecast")
for(i in 1:4) lines(Mon_calSPS[i,],col=Cols[i],lwd=2)
for(i in 1:4) lines(Mon_rawSPS[i,],col=Cols[i],lty=2,lwd=2)
legend("topleft",legend = month.abb[c(1,4,7,10)],col = Cols,lty=1,bty = "n",title = "Initialised")

load("~/Documents/Data/BiasCorrection/CollectedSPSResultPlusLongjiang",envir = (LongEn=new.env()))
for(i in 1:4) lines(LongEn$Mon_longSPS[i,],col=Cols[i],lty=3,lwd=2)
legend("bottomright",legend = c("Raw Forecast","Longjiang corrected","TAQM calibrated"),col = Cols[c(i,i,i)],lty=c(2,3,1),bty = "n")


load("~/Documents/Data/BiasCorrection/CollectedLongCalSPSResult",envir = (LongEn=new.env()))
for(i in 1:4) lines(LongEn$Mon_rawSPS[i,],col=Cols[i],lty=3,lwd=2)
legend("bottomright",legend = c("Raw Forecast","Longjiang corrected","TAQM calibrated"),col = Cols[c(i,i,i)],lty=c(2,3,1))



initlist=c(1,4,7,10)
i=1
par(mar=c(5,5,4,1)+.1) 
plot(Mon_rawSPS[i,],col=Cols[i],lty=2,lwd=2,type="l",ylim=c(0,8),yaxs="i",xlab="Leadtime (months)",ylab=expression(paste("SPS (10"^"6"," km"^"2",")")),main=paste0("Mean SPS for forecasts initialised in ",month.abb[initlist[i]]))
lines(Mon_calSPS[i,],col=Cols[i],lwd=2)
# legend("bottomright",legend = c("Raw Forecast","Calibrated Forecast"),col = Cols[c(i,i)],lty=c(2,1))

lines(LongEn$Mon_longSPS[i,],col=Cols[i],lty=3,lwd=2)



##### IIEE
load("~/Documents/Data/BiasCorrection/Collected_IIEE_Results")
Cols=c("#1b9e77","#d95f02","#7570b3","#e7298a")
par(mar=c(5,5,4,1)+.1) #Add this before calling plot to add margin for ylabel
plot(0,0,xlim=c(1,12),ylim=c(0,10),xlab="Leadtime (months)",ylab=expression(paste("IIEE (10"^"6"," km"^"2",")")),pch='.',yaxs="i",main="Calibrated and raw IIEE forecast")
for(i in 1:4) lines(Mon_calIIEE[i,],col=Cols[i],lwd=2)
for(i in 1:4) lines(Mon_rawIIEE[i,],col=Cols[i],lty=2,lwd=2)
legend("topleft",legend = c(month.abb[c(1,4,7,10)],"Raw forecast"),col = c(Cols,"black"),lty=c(1,1,1,1,2),bty = "n",title = "Initialised")
grid()



##
load("~/Documents/Data/BiasCorrection/EstMeanConc_IIEEResults")

Cols=c("#1b9e77","#d95f02","#7570b3","#e7298a")
par(mar=c(5,5,4,1)+.1) #Add this before calling plot to add margin for ylabel
plot(0,0,xlim=c(1,12),ylim=c(1,8),xlab="Leadtime (months)",ylab=expression(paste("IIEE (10"^"6"," km"^"2",")")),pch='.',yaxs="i",main="IIEE of ensemble mean")
for(i in 1:4) lines(Mon_calIIEE[i,],col=Cols[i],lwd=2)
for(i in 1:4) lines(Mon_rawIIEE[i,],col=Cols[i],lty=2,lwd=2)
legend("topleft",ncol = 2,legend = month.abb[c(1,4,7,10)],col = Cols,lty=1,bty = "n",title = "Initialised")
