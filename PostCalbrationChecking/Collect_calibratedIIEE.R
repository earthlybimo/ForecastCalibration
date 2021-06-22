### Script to open calibrated SSIPS forecast, compute IIEE, then save it into a table.

args = commandArgs(trailingOnly=TRUE)
ccc=as.integer(args[1]) #Choice of Hemisphere


library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'
Clima_path= '/work/ab0995/a270112/data_fesom2/sic'
SIMPdata_path='/pf/a/a270138/Data/LongjiangOutput/FCST_CLIM'

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
Hemfilter=array(0,dim=length(grd$lat))

if(ccc == 1)  {
  HEM="nh"
  Hemfilter[grd$lat>0]=1
}
if(ccc == 2)  {
  HEM="sh"
  Hemfilter[grd$lat<0]=1
}
Allsavename=paste0("/work/ba1138/a270138/BiasCorrOutput/CollectedIIEE_Results_",HEM)

inYR=2011:2018  #which year
# inMON=1:4   #And which initialisation 
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation

rawOarr=array(dim=c(length(inYR),4,12))
rawUarr=rawOarr
calOarr=rawOarr
calUarr=rawUarr
climaOarr=array(dim=c(length(inYR),4,12))
climaUarr=climaOarr
simpOarr=climaOarr
simpUarr=climaUarr

for(yy in 1:length(inYR)){ # yy=1;init=1;mm=1 #For tests
  loadname=sprintf("%s/F%0.2d_ens_mon_mean_SIP_corr.nc",Data_path,(inYR[yy]-2000))
  # if(!file.exists(loadname)) next()
  fl=nc_open(loadname)
  simpSIP_all=ncvar_get(fl,"SIP_FCST_CORR")
  nc_close(fl)
  
  for(init in 1:4){
    for(mm in 1:12){
      ### Get SIP:
      
      # Obs
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
      obsSIP[(obsVar2>=0.15)&(obsVar2<=1)]=1
      obsSIP[(obsVar2>=0)&(obsVar2<0.15)]=0
      
      # TAQM cal forecasts
      loadname=sprintf("%s/Forecast_Calibration_TrustSharpFalse_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
      if(!file.exists(loadname)) next()
      fl=nc_open(loadname)
      rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
      calSIP=ncvar_get(fl,"SIP_FCST_CORR")
      # obsSIP2=ncvar_get(fl,"SIP_FCST_OBS")
      nc_close(fl)
      
      ## Climatology:
      # Let's make our own climatlogy:
      ClimaSIParr=array(dim =c(9,length(grd$lat)))
      for(l in 9:1){
        loadname=paste0("/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_",(obsTyr-l),".nc")
        if(!file.exists(loadname)) {print(paste0("Not found: ",basename(loadname)));next()}
        fl=nc_open(loadname)
        obsVar1=ncvar_get(fl,"obs") # num [1:126858, 1:12]
        nc_close(fl)
        obsVar2=obsVar1[,obsTmnth]
        obsSIPtemp=array(dim =length(grd$lat))
        obsSIPtemp[(obsVar2>=0.15)&(obsVar2<=1)]=1
        obsSIPtemp[(obsVar2>=0)&(obsVar2<0.15)]=0
        ClimaSIParr[l,]=obsSIPtemp
      }
      climaSIP=colMeans(ClimaSIParr,na.rm = T)
      
      ## SIMP cal forecast:
      simpSIP=simpSIP_all[,mm,init]
      ## But simp calibration is just calibrated ensemble mean so this is potentially not correct
      
      
      ### Make medians (for IIEE)
      climaSIPmed=binarise(climaSIP,0.5)
      rawSIPmed=binarise(rawSIP,0.5)
      calSIPmed=binarise(calSIP,0.5)
      simpSIPmed=binarise(simpSIP,0.5)
      
      
      ### Measure IIEE (O + U)
    
      # raw
      uid=which((rawSIPmed==0)&(obsSIP==1)&(Hemfilter==1)) #identify underforecast grids
      U=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      oid=which((rawSIPmed==1)&(obsSIP==0)&(Hemfilter==1)) #identify overforecast grids
      O=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      rawUarr[yy,init,mm]=U*(10^-12)
      rawOarr[yy,init,mm]=O*(10^-12)
      remove(O,U,uid,oid,rawSIPmed,rawSIP)
      
      # Cal:
      uid=which((calSIPmed==0)&(obsSIP==1)&(Hemfilter==1)) #identify underforecast grids
      U=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      oid=which((calSIPmed==1)&(obsSIP==0)&(Hemfilter==1)) #identify overforecast grids
      O=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      calUarr[yy,init,mm]=U*(10^-12)
      calOarr[yy,init,mm]=O*(10^-12)
      
      remove(O,U,uid,oid,calSIPmed,calSIP)
      
      # Clima:
      uid=which((climaSIPmed==0)&(obsSIP==1)&(Hemfilter==1)) #identify underforecast grids
      U=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      oid=which((climaSIPmed==1)&(obsSIP==0)&(Hemfilter==1)) #identify overforecast grids
      O=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      climaUarr[yy,init,mm]=U*(10^-12)
      climaOarr[yy,init,mm]=O*(10^-12)
      remove(O,U,uid,oid,climaSIPmed,climaSIC)
      
      #SIMP:
      uid=which((simpSIPmed==0)&(obsSIP==1)&(Hemfilter==1)) #identify underforecast grids
      U=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      oid=which((simpSIPmed==1)&(obsSIP==0)&(Hemfilter==1)) #identify overforecast grids
      O=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      simpUarr[yy,init,mm]=U*(10^-12)
      simpOarr[yy,init,mm]=O*(10^-12)
      remove(O,U,uid,oid,simpSIPmed,simpSIP)
      
      remove(obsSIP)
    }
  }
}

rawIIEEarr=rawOarr+rawUarr
calIIEEarr=calOarr+calUarr
climaIIEEarr=climaOarr+climaUarr
simpIIEEarr=simpOarr+simpUarr
Mon_rawIIEE=apply(rawIIEEarr, c(2,3),mean,na.rm=T)
Mon_calIIEE=apply(calIIEEarr, c(2,3),mean,na.rm=T)
Mon_climaIIEE=apply(climaIIEEarr, c(2,3),mean,na.rm=T)
Mon_simpIIEE=apply(simpIIEEarr, c(2,3),mean,na.rm=T)

# rawAEEarr=abs(rawOarr- rawUarr )
# rawMEarr=rawIIEEarr-rawAEEarr
# calAEEarr=abs(calOarr-calUarr)
# calMEarr=calIIEEarr-calAEEarr

save(file = Allsavename,version = 2,grd,calIIEEarr,rawIIEEarr,climaIIEEarr,simpIIEEarr,inYR,Mon_rawIIEE,Mon_calIIEE,Mon_simpIIEE,Mon_climaIIEE,rawOarr,calOarr,simpOarr,climaOarr)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)),overwrite = T)
print(paste0("Done! saved file:",basename(Allsavename)))

