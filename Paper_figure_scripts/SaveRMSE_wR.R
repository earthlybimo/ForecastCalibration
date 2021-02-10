### Script to open calibrated SSIPS forecast, compute RMSE from estimated Mean, then save it into a table.

library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults/MeanConc'
Allsavename=paste0("/work/ba1138/a270138/BiasCorrOutput/EstMeanConc_RMSE_Results")

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
NHgrdpts=which(grd$lat>10)
SHgrdpts=which(grd$lat<(-10))
area=grd$cell_area
inYR=2011:2018  #which year
# inMON=1:4   #And which initialisation 
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation

rawNHarr=array(dim=c(length(inYR),4,12))
rawSHarr=rawNHarr
calNHarr=rawNHarr
calSHarr=rawNHarr

for(yy in 1:length(inYR)){
  for(init in 1:4){
    for(mm in 1:12){
      
      ## Get Obs
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
      obsSIC=obsVar1[,obsTmnth]
    
      
      loadname=sprintf("%s/Forecast_Calibration_TrustSharpFalse_wMeanConc_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
      if(!file.exists(loadname)) next()
      
      fl=nc_open(loadname)
      # rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
      # calSIP=ncvar_get(fl,"SIP_FCST_CORR")
      rawSIC=ncvar_get(fl,"SIC_RAW")
      calSIC=ncvar_get(fl,"SIC_Cal")
      # obsSIP2=ncvar_get(fl,"SIP_FCST_OBS")
      nc_close(fl)

      
      temp=(obsSIC-rawSIC)
      temp[abs(temp)>1]=NA
      rawNHarr[yy,init,mm]=sqrt(sum(temp[NHgrdpts],na.rm = T)/sum(area[NHgrdpts],na.rm = T))
      rawSHarr[yy,init,mm]=sqrt(sum(temp[SHgrdpts],na.rm = T)/sum(area[SHgrdpts],na.rm = T))
    
      temp=(obsSIC-calSIC);temp[abs(temp)>1]=NA  
      calNHarr[yy,init,mm]=sqrt(sum(temp[NHgrdpts],na.rm = T)/sum(area[NHgrdpts],na.rm = T))
      calSHarr[yy,init,mm]=sqrt(sum(temp[SHgrdpts],na.rm = T)/sum(area[SHgrdpts],na.rm = T))
    
      
    }
  }
}


save(file = Allsavename,version = 2,grd,rawSHarr,rawNHarr,inYR,calSHarr,calNHarr)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)))
print("Done!")

