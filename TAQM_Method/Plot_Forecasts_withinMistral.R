### Script to plot maps from Longiang's outputs. or my output, if need be

library(ncdf4);library(spheRlab)
grd = sl.grid.readNCDF("/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc")

data_path='/work/ba1138/a270112/awicm3/FCST_CLIM'
Figpath="/pf/a/a270138/Data/BiasCorrOutput/Figs"
save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'

inYR=2011:2018  #among which year
# longJcalSPSarr=array(dim=c(length(inYR),4,12))
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation

#Which forecast, what time.
yy=2  #2012
init=3
mm=1 #lead time month
for (yy in 1:8) {
  
  
  ## First for Longjiang's output
  loadname=sprintf("%s/F%02d_MEM_ens_mon_mean_corr.nc",data_path,(inYR[yy]-2000))
  fl=nc_open(loadname)
  longcalSIC=ncvar_get(fl,"SIC_FCST_CORR") # num [1:126858, 1:12, 1:4, 1:30]
  nc_close(fl)
  longcalSIPo=array(dim=dim(longcalSIC))
  longcalSIPo[longcalSIC<0.15]=0
  longcalSIPo[longcalSIC>=0.15]=1
  longcalSIP=apply(longcalSIPo,c(1,2,3),mean)

  obsTyr=inYR[yy]
  obsTmnth = mm+strtm[init]-1
  if (obsTmnth>12){
    obsTyr=obsTyr+1
    obsTmnth=obsTmnth-12}

  ###print('Targetyear = '+str(targetyear)+',initialisation = '+str(init)+' which means from '+ str(strtm[init-1]) +',leadtime '+str(leadtimeMonth)+' so target month is '+str(obsTmnth)+' of year '+str(obsTyr))  # Testing, python code so need to correct

  loadname=paste0("/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_",obsTyr,".nc")
  if(!file.exists(loadname)) next()
  fl=nc_open(loadname)
  obsVar1=ncvar_get(fl,"obs") # num [1:126858, 1:12]
  nc_close(fl)

  obsVar2=obsVar1[,obsTmnth]
  obsSIP=array(dim =length(grd$lat))
  obsSIP[obsVar2>=0.15]=1;obsSIP[obsVar2<0.15]=0

  preSPS=(longcalSIP[,mm,init]-obsSIP)^2   #Diff between model and satelite
  preSPS2=preSPS*grd$cell_area
  SPScal=sum(preSPS2,na.rm = T)*(10^-12)

  ## Plot Forecast
  pir = sl.plot.init(projection="polar",polar.latbound = 50,file.name = sprintf("%s/LongJSIP_%d_init%d_leadtime%02d.pdf",Figpath,inYR[yy],init,mm))
  pcol= sl.plot.field.elem(pir,num=longcalSIP[,mm,init],lon=grd$lon,lat = grd$lat,elem = grd$elem)
  res = sl.plot.naturalearth(pir,fill.col = "grey",lines.col = "grey")
  res=sl.plot.text(pir,lon=120,lat=65,labels = paste0("SPS = ",SPScal))
  sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
  sl.plot.end(pir)

  ## Plot Obs
  pir = sl.plot.init(projection="polar",polar.latbound = 50,file.name = sprintf("%s/ObsSIP_%d_month%02d.pdf",Figpath,obsTyr,obsTmnth))
  pcol= sl.plot.field.elem(pir,num=obsSIP,lon=grd$lon,lat = grd$lat,elem = grd$elem)
  res = sl.plot.naturalearth(pir,fill.col = "grey",lines.col = "grey")
  sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
  sl.plot.end(pir)
  
  
  ### Let's add the TAQM calibrated forecasts to the lot:
  
  loadname=sprintf("%s/Forecast_Calibration_BigHist_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
  if(!file.exists(loadname)) next()
  
  fl=nc_open(loadname)
  rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
  calSIP=ncvar_get(fl,"SIP_FCST_CORR")
  # obsSIP=ncvar_get(fl,"SIP_FCST_OBS")
  nc_close(fl)
  
  preSPS=(rawSIP-obsSIP)^2   #Diff between model and satelite
  preSPS2=preSPS*grd$cell_area
  SPSraw=sum(preSPS2)*(10^-12)

  remove(preSPS2,preSPS)
  preSPS=(calSIP-obsSIP)^2   #Diff between model and satelite
  preSPS2=preSPS*grd$cell_area
  SPScal=sum(preSPS2)*(10^-12)

  ## Plot Forecast
  pir = sl.plot.init(projection="polar",polar.latbound = 50,file.name = sprintf("%s/TAQMcalSIP_%d_init%d_leadtime%02d.pdf",Figpath,inYR[yy],init,mm))
  pcol= sl.plot.field.elem(pir,num=calSIP,lon=grd$lon,lat = grd$lat,elem = grd$elem)
  res = sl.plot.naturalearth(pir,fill.col = "grey",lines.col = "grey")
  res=sl.plot.text(pir,lon=120,lat=65,labels = paste0("SPS = ",SPScal))
  sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
  sl.plot.end(pir)
  
  ## Plot 'Raw'  (quotes because of the way we went from ensemble SIC to SIP with TAQM method)
  pir = sl.plot.init(projection="polar",polar.latbound = 50,file.name = sprintf("%s/rawSIP_%d_init%d_leadtime%02d.pdf",Figpath,inYR[yy],init,mm))
  pcol= sl.plot.field.elem(pir,num=rawSIP,lon=grd$lon,lat = grd$lat,elem = grd$elem)
  res = sl.plot.naturalearth(pir,fill.col = "grey",lines.col = "grey")
  res=sl.plot.text(pir,lon=120,lat=65,labels = paste0("SPS = ",SPSraw))
  sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
  sl.plot.end(pir)
  
  
}