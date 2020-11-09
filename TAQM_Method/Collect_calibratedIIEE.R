### Script to open calibrated SSIPS forecast, compute IIEE, then save it into a table.

library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'
Allsavename=paste0("/work/ba1138/a270138/BiasCorrOutput/CollectedIIEE_Results")

binarise <-function (somearr,dlevel) {  #Function to binarise some given array  based on  this level
  ll=dim(somearr)
  if (is.null(ll)){ll=length(somearr)}  #because 1d arrays dont do
  ibinar=array(dim =ll)
  ibinar[somearr[]>=dlevel]=1
  ibinar[somearr[]<dlevel]=0
  return (ibinar)
}

##First we need the gridinfo to compute IIEE
grd = sl.grid.readNCDF("/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc")

inYR=2011:2018  #which year
# inMON=1:4   #And which initialisation 
rawOarr=array(dim=c(length(inYR),4,12))
rawUarr=rawOarr
calOarr=rawOarr
calUarr=rawUarr

for(yy in 1:length(inYR)){
  for(init in 1:4){
    for(mm in 1:12){
      loadname=sprintf("%s/Forecast_Calibration_BigHist_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
      if(!file.exists(loadname)) next()
      
      fl=nc_open(loadname)
      rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
      calSIP=ncvar_get(fl,"SIP_FCST_CORR")
      obsSIP=ncvar_get(fl,"SIP_FCST_OBS")
      nc_close(fl)
      
      rawSIPmed=binarise(rawSIP,0.5)
      temp=(obsSIP-rawSIPmed)
      tempU=array(0,dim=dim(temp));tempO=tempU
      
      tempU[temp>0]=1  #Underforecast area
      U=sum(tempU*grd$cell_area,na.rm=T)  #And sum
      tempO[temp<0]=1 #Overforecast area
      O=sum(tempO*grd$cell_area,na.rm=T)
      rawUarr[yy,init,mm]=U*(10^-12)
      rawOarr[yy,init,mm]=O*(10^-12)
      
      remove(O,U,temp,rawSIPmed,tempU,tempO)
      
      calSIPmed=binarise(calSIP,0.5)
      temp=(obsSIP-calSIPmed)
      tempU=array(0,dim=dim(temp));tempO=tempU
      
      tempU[temp>0]=1  #Underforecast area
      U=sum(tempU*grd$cell_area,na.rm=T)  #And sum
      tempO[temp<0]=1 #Overforecast area
      O=sum(tempO*grd$cell_area,na.rm=T)
      calUarr[yy,init,mm]=U*(10^-12)
      calOarr[yy,init,mm]=O*(10^-12)
      remove(O,U,temp,calSIPmed,tempU,tempO)
      
      remove(rawSIP,calSIP,obsSIP)
      
    }
  }
}

rawIIEEarr=rawOarr+rawUarr
calIIEEarr=calOarr+calUarr
Mon_rawIIEE=apply(rawIIEEarr, c(2,3),mean,na.rm=T)
Mon_calIIEE=apply(calIIEEarr, c(2,3),mean,na.rm=T)


# rawAEEarr=abs(rawOarr- rawUarr )
# rawMEarr=rawIIEEarr-rawAEEarr
# calAEEarr=abs(calOarr-calUarr)
# calMEarr=calIIEEarr-calAEEarr

save(file = Allsavename,version = 2,grd,calIIEEarr,rawIIEEarr,inYR,Mon_rawIIEE,Mon_calIIEE,rawOarr,calOarr)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)))
print("Done!")

