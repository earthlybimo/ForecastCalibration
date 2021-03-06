### Script to open calibrated SSIPS forecast, compute SPS, then save it into a table.

library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'
data_path='/work/ba1138/a270112/awicm3/FCST_CLIM'  #Where the obs are
Allsavename=paste0(save_path,"/CollectedSPSResult_TrustSharpFalse_corrected")

##First we need the gridinfo to compute SPS
grd = sl.grid.readNCDF("/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc")

inYR=2011:2018  #among which year
# inMON=1:4   #And which initialisation 
rawSPSarr=array(dim=c(length(inYR),4,12))
calSPSarr=rawSPSarr
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation

for(yy in 1:length(inYR)){
  ## Section to collect SPS for forecast that is already bias corrected by longjiang
  obsTyr=inYR[yy]
  for(init in 1:4){
    for(mm in 1:12){
      obsTyr=inYR[yy]
      obsTmnth = mm+strtm[init]-1
      if (obsTmnth>12){
        obsTyr=obsTyr+1
        obsTmnth=obsTmnth-12}
      
      # print(paste0('Targetyear = ',(inYR[yy]),',initialisation = ',(init),' which means from ', (strtm[init]) ,',leadtime ',(mm),' so target month is ',(obsTmnth),' of year ',(obsTyr)))  # Testing
      
      loadname=paste0("/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_",obsTyr,".nc")
      if(!file.exists(loadname)) next()
      fl=nc_open(loadname)
      obsVar1=ncvar_get(fl,"obs") # num [1:126858, 1:12]
      nc_close(fl)
      obsVar2=obsVar1[,obsTmnth]
      obsSIP=array(dim =length(grd$lat))
      obsSIP[obsVar2>=0.15]=1
      obsSIP[obsVar2<0.15]=0

      # loadname=sprintf("%s/Forecast_Calibration_TrustSharpFalse_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)      
      loadname=sprintf("%s/Forecast_Calibration_TrustSharpFalse_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
      if(!file.exists(loadname)) next()
      
      fl=nc_open(loadname)
      rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
      calSIP=ncvar_get(fl,"SIP_FCST_CORR")
      # obsSIP2=ncvar_get(fl,"SIP_FCST_OBS")  #This obs was found to have turned some NA points to 0, so slightly wrong. 
      nc_close(fl)
      
      preSPS=(rawSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPSraw=sum(preSPS2,na.rm = T)*(10^-12)
      rawSPSarr[yy,init,mm]=SPSraw
      remove(preSPS2,preSPS)
      preSPS=(calSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPScal=sum(preSPS2,na.rm = T)*(10^-12)
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

