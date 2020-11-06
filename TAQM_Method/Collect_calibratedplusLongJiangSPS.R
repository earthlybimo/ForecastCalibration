### Script to open calibrated SSIPS forecast, compute SPS, then save it into a table.

library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'  #Where TAQMcalibrated forecasts were saved
data_path='/work/ba1138/a270112/awicm3/FCST_CLIM'

Allsavename=paste0("/work/ba1138/a270138/BiasCorrOutput/CollectedSPSResultPlusLongjiang")

##First we need the gridinfo to compute SPS
grd = sl.grid.readNCDF("/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc")

inYR=2011:2018  #which year
# inMON=1:4   #And which initialisation 
rawSPSarr=array(dim=c(length(inYR),4,12))
calSPSarr=rawSPSarr
longJcalSPSarr=rawSPSarr

for(yy in 1:length(inYR)){
  ## Section to collect SPS for forecast that is already bias corrected by longjian
  loadname=sprintf("%s/F%02d_MEM_ens_mon_mean_corr.nc",data_path,(inYR[yy]-2000))
  #Let's HOPE it all exists, because not sure how to do order of operation
  # if(!file.exists(loadname)) next()
  fl=nc_open(loadname)
  longcalSIC=ncvar_get(fl,"SIC_FCST_CORR") # num [1:126858, 1:12, 1:4, 1:30]
  nc_close(fl)
  longcalSIPo=array(dim=dim(longcalSIC))
  longcalSIPo[longcalSIC<0.15]=0
  longcalSIPo[longcalSIC>=0.15]=1
  longcalSIP=apply(longcalSIPo,c(1,2,3),mean)
  
  for(init in 1:4){
    for(mm in 1:12){
      loadname=sprintf("%s/Forecast_Calibration_BigHist_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
      if(!file.exists(loadname)) next()
      
      fl=nc_open(loadname)
      rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
      calSIP=ncvar_get(fl,"SIP_FCST_CORR")
      obsSIP=ncvar_get(fl,"SIP_FCST_OBS")
      nc_close(fl)
      
      preSPS=(rawSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPSraw=sum(preSPS2)*(10^-12)
      rawSPSarr[yy,init,mm]=SPSraw
      remove(preSPS2,preSPS)
      
      preSPS=(calSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPScal=sum(preSPS2)*(10^-12)
      calSPSarr[yy,init,mm]=SPScal
      remove(preSPS2,preSPS)
      
      preSPS=(longcalSIP[,mm,init]-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPScal=sum(preSPS2)*(10^-12)
      longJcalSPSarr[yy,init,mm]=SPScal
      
      remove(rawSIP,calSIP,obsSIP,SPSraw,SPScal,preSPS2,preSPS)
      
    }
  }
}

Mon_rawSPS=apply(rawSPSarr, c(2,3),mean,na.rm=T)
Mon_calSPS=apply(calSPSarr, c(2,3),mean,na.rm=T)
Mon_longSPS=apply(longJcalSPSarr, c(2,3),mean,na.rm=T)

save(file = Allsavename,version = 2,grd,calSPSarr,rawSPSarr,inYR,Mon_rawSPS,Mon_calSPS,longJcalSPSarr,Mon_longSPS)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)))
print("Done!")

