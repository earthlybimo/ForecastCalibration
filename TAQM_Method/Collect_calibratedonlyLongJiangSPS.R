library(ncdf4);library(spheRlab)
### Script to compute and save SPS for forecasts calibrated by Longjiang's method, using the original files saved by Longjiang.

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'  #Where TAQMcalibrated forecasts were saved
data_path='/work/ba1138/a270112/awicm3/FCST_CLIM'

Allsavename=paste0("/work/ba1138/a270138/BiasCorrOutput/CollectedSPSResult_onlyLongjiang")

##First we need the gridinfo to compute SPS
grd = sl.grid.readNCDF("/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc")

inYR=2011:2018  #among which year
# inMON=1:4   #And which initialisation 
longJcalSPSarr=array(dim=c(length(inYR),4,12))
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation

for(yy in 1:length(inYR)){
  ## Section to collect SPS for forecast that is already bias corrected by longjiang
  loadname=sprintf("%s/F%02d_MEM_ens_mon_mean_corr.nc",data_path,(inYR[yy]-2000))
  #Let's assume it all exists
  # if(!file.exists(loadname)) next()
  fl=nc_open(loadname)
  longcalSIC=ncvar_get(fl,"SIC_FCST_CORR") # num [1:126858, 1:12, 1:4, 1:30]
  nc_close(fl)
  longcalSIPo=array(dim=dim(longcalSIC))
  longcalSIPo[longcalSIC<0.15]=0
  longcalSIPo[longcalSIC>=0.15]=1
  longcalSIP=apply(longcalSIPo,c(1,2,3),mean)
  obsTyr=inYR[yy]
  for(init in 1:4){
    for(mm in 1:12){
      obsTyr=inYR[yy]
      obsTmnth = mm+strtm[init]-1
      if (obsTmnth>12){
        obsTyr=obsTyr+1
        obsTmnth=obsTmnth-12}
      
      # print('Targetyear = '+str(targetyear)+',initialisation = '+str(init)+' which means from '+ str(strtm[init-1]) +',leadtime '+str(leadtimeMonth)+' so target month is '+str(obsTmnth)+' of year '+str(obsTyr))  # Testing
      
      loadname=paste0("/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_",obsTyr,".nc")
      if(!file.exists(loadname)) next()
      fl=nc_open(loadname)
      obsVar1=ncvar_get(fl,"obs") # num [1:126858, 1:12]
      nc_close(fl)
      
      obsVar2=obsVar1[,obsTmnth]
      obsSIP=array(dim =length(grd$lat))
      obsSIP[obsVar2>=0.15]=1
      obsSIP[obsVar2<0.15]=0
      
      preSPS=(longcalSIP[,mm,init]-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPScal=sum(preSPS2,na.rm = T)*(10^-12)
      longJcalSPSarr[yy,init,mm]=SPScal
      
      remove(obsSIP,obsVar1,obsVar2,SPScal,preSPS2,preSPS)
      
    }
  }
}
Mon_longSPS=apply(longJcalSPSarr, c(2,3),mean,na.rm=T)

save(file = Allsavename,version = 2,grd,inYR,longJcalSPSarr,Mon_longSPS)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)))
print("Done!")
