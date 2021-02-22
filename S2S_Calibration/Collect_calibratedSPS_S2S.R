library(ncdf4);library(spheRlab)


data_path = '/work/ba1138/a270138/BiasCorrOutput/S2S_Results/'
modelname="KMA";flen=60
ylist=2005:2010
save_name=paste0(data_path,"CollectedCalibratedSPS_for_",modelname)

## Find the grid areas
Gridfile="/mnt/lustre01/work/ab0995/a270138/S2S_Output/GridandNApoints"
load(Gridfile,envir =(GridEnv=new.env()))
grd=GridEnv$grd.full
Elementareas=array(dim=dim(grd$elem)[[1]])
for (i in 1:length(Elementareas)){
  nodes=grd$elem[i,]   #What nodes does this element have?
  Elementareas[i]=sl.triag.area(grd$lon[nodes],grd$lat[nodes]) #And how much is the weight here?
}
Nodeareas=array(dim=length(grd$lon))
for (i in 1:length(Nodeareas)){
  elems=grd$neighelems[i,]  #Most nodes have upto 6 elements around it, 
  areaofelems=Elementareas[elems]
  areaofelems=areaofelems[!is.na(areaofelems)]  #but some have less than 6, so some of those points could be NA
  Nodeareas[i]=sum(areaofelems)*(1/3)  #Each element area is shared by 3 nodes, and nodes get shares from all the elements its part of.
}

## Where to save the SPS values?
calSPSarr=array(dim=c(length(ylist),12,flen))
rawSPSarr=calSPSarr
for(yy in 1:length(ylist)){
  for(mm in 1:12){
    file1=Sys.glob(sprintf("%s%s/TAQM_calibrated_*%d-%0.2d-??.nc",data_path,modelname,ylist[yy],mm))
    if(length(file1)==0) next()
    
    fl=nc_open(file1)
    # fl$var$SIP_
    CalSIP=ncvar_get(fl,"SIP_FCST_CORR")
    ObsSIP=ncvar_get(fl,"SIP_FCST_OBS")
    RawSIP=ncvar_get(fl,"SIP_FCST_RAW")
    nc_close(fl)
    
    ObsSIP[ObsSIP<0.15]=0
    ObsSIP[ObsSIP>=0.15]=1
    CalSIP[CalSIP<0]=NA
    RawSIP[RawSIP<0]=NA
    # CalSIP[is.na(RawSIP)]=NA
    # RawSIP[is.na(CalSIP)]=NA
    
    # Now calculate the SPS
    for (i in 1:flen) {
      temp=as.vector(CalSIP[,,i]-ObsSIP[,,i])
      temp2=(temp^2)*Nodeareas
      calSPSarr[yy,mm,i]=sum(temp2,na.rm = T)
      remove(temp,temp2)
      temp=as.vector(RawSIP[,,i]-ObsSIP[,,i])
      temp2=(temp^2)*Nodeareas
      rawSPSarr[yy,mm,i]=sum(temp2,na.rm = T)
      remove(temp,temp2)
    }
  }
}

calMean=colMeans(calSPSarr)
rawMean=colMeans(rawSPSarr)
Fctr=((6371^2)/1000000)



save(file = save_name,version = 2, calSPSarr,rawSPSarr,ylist,calMean,rawMean,Fctr)
print(paste0("File saved:",save_name))


cols=rainbow(12)
plot(-1,-1,xlim=c(1,46),ylim=c(0,2),ylab="SPS (mill sq km)",xlab="Leadtime (days)")

for(i in seq(1,12,2)){
  lines(rawMean[i,1:45]*Fctr,col=cols[i],lty=1)
  lines(calMean[i,1:45]*Fctr,col=cols[i],lty=2)
}
legend("topleft",legend = month.abb[seq(1,12,2)],col = cols[seq(1,12,2)],lty=3)

dev.off()


# load("~/Documents/Data/BiasCorrection/S2S_calibration/CollectedCalibratedSPS_for_ECMWF")
# 
# 
# plot(as.Date("2004-01-01"),0,xlim=as.Date(c("2008-01-01","2010-01-02")),ylim=c(0,2),ylab="SPS (mill sq km)",xlab="")
# yy=1;mm=2
# 
# for(yy in 1:6){
#   for(mm in seq(1,12,2)){
#     lines((as.Date(sprintf("%d-%0.2d-01",ylist[yy],mm))+(0:45)),rawSPSarr[yy,mm,1:46]*Fctr,col=cols[mm],lty=1)
#     lines((as.Date(sprintf("%d-%0.2d-01",ylist[yy],mm))+(0:45)),calSPSarr[yy,mm,1:46]*Fctr,col=cols[mm],lty=2)}
# }
# 
# 







