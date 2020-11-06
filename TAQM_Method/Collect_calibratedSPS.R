### Script to open calibrated SSIPS forecast, compute SPS, then save it into a table.

library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'
Allsavename=paste0(save_path,"/CollectedSPSResult")

##First we need the gridinfo to compute SPS
grd = sl.grid.readNCDF("/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc")

inYR=2011:2018  #which year
# inMON=1:4   #And which initialisation 
rawSPSarr=array(dim=c(length(inYR),4,12))
calSPSarr=rawSPSarr

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
      
      preSPS=(rawSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPSraw=sum(preSPS2)*(10^-12)
      rawSPSarr[yy,init,mm]=SPSraw
      remove(preSPS2,preSPS)
      preSPS=(calSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPScal=sum(preSPS2)*(10^-12)
      calSPSarr[yy,init,mm]=SPScal
      
      remove(rawSIP,calSIP,obsSIP,SPSraw,SPScal)
      
    }
  }
}

Mon_rawSPS=apply(rawSPSarr, c(2,3),mean,na.rm=T)
Mon_calSPS=apply(calSPSarr, c(2,3),mean,na.rm=T)

save(file = Allsavename,version = 2,grd,calSPSarr,rawSPSarr,inYR,Mon_rawSPS,Mon_calSPS)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)))
print("Done!")

