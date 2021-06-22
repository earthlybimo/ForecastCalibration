### Script to open calibrated SSIPS forecast, compute SPS, then save it into a table.
args = commandArgs(trailingOnly=TRUE)
ccc=as.integer(args[1]) #Choice of Hemisphere

library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'
Clima_path= '/work/ab0995/a270112/data_fesom2/sic'
SIMPdata_path='/pf/a/a270138/Data/LongjiangOutput/FCST_CLIM'

##First we need the gridinfo to compute SPS
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
Allsavename=paste0(save_path,"/CollectedSPSResult_AllinclTAQM_",HEM)



inYR=2011:2018  #which year
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation
rawSPSarr=array(dim=c(length(inYR),4,12))
calSPSarr=rawSPSarr
climaSPSarr=array(dim=c(length(inYR),4,12))
simpSPSarr=climaSPSarr

for(yy in 1:length(inYR)){
  
  loadname=sprintf("%s/F%0.2d_ens_mon_mean_SIP_corr.nc",SIMPdata_path,(inYR[yy]-2000))
  # if(!file.exists(loadname)) next()
  fl=nc_open(loadname)
  SIMPcalSIP_all=ncvar_get(fl,"SIP_FCST_CORR")
  # This assumes the simple forecasts are stored as SIP in  grid x mm x initialisation for each year
  nc_close(fl)
  
  for(init in 1:4){
    for(mm in 1:12){
      ### Load variables and get SIP:
      
      ## TAQM
      loadname=sprintf("%s/Forecast_Calibration_TrustSharpFalse_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
      if(!file.exists(loadname)) next()
      
      fl=nc_open(loadname)
      rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
      calSIP=ncvar_get(fl,"SIP_FCST_CORR")
      # obsSIP=ncvar_get(fl,"SIP_FCST_OBS")
      nc_close(fl)
      
      
      ###Obs! Let's compute it:
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
      
      ### Climatology:
      # ClimFile=paste0(Clima_path,"/OSISAF_MON_CLIM_",(obsTyr-9),"-",(obsTyr-1),".nc")
      # if(!file.exists(ClimFile)) {print(paste0("This Climafile doesnt exist! ",basename(ClimFile))); next()}
      # fl=nc_open(ClimFile)
      # climaSIC=ncvar_get(fl,"obs")
      # nc_close(fl)
      # climaSIP=binarise(climaSIC[,obsTmnth],0.15)  #This is not clear whether it is climatological SIC, or SIP..
      
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
      
      ## SIMP:
      SIMPcalSIP=SIMPcalSIP_all[,mm,init]
      ## But simp calibration is just calibrated ensemble mean so not correct
     
      
      ### Compute SPS ### ------------
      
      # Clima:
      preSPS=(climaSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPSclima=sum(preSPS2[Hemfilter==1],na.rm = T)*(10^-12)
      climaSPSarr[yy,init,mm]=SPSclima
      remove(preSPS2,preSPS)
      
      # Raw:
      preSPS=(rawSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPSraw=sum(preSPS2[Hemfilter==1],na.rm = T)*(10^-12)
      rawSPSarr[yy,init,mm]=SPSraw
      remove(preSPS2,preSPS)
      
      # TAQM:
      preSPS=(calSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPScal=sum(preSPS2[Hemfilter==1],na.rm = T)*(10^-12)
      calSPSarr[yy,init,mm]=SPScal
      
      # SIMP:
      preSPS=(SIMPcalSIP-obsSIP)^2   #Diff between model and satelite
      preSPS2=preSPS*grd$cell_area
      SPScal=sum(preSPS2[Hemfilter==1],na.rm = T)*(10^-12)
      simpSPSarr[yy,init,mm]=SPScal
      
      remove(rawSIP,calSIP,obsSIP,SPSraw,SPScal,climaSIP,SIMPcalSIP)
      
    }
  }
}

Mon_rawSPS=apply(rawSPSarr, c(2,3),mean,na.rm=T)
Mon_calSPS=apply(calSPSarr, c(2,3),mean,na.rm=T)
Mon_climaSPS=apply(climaSPSarr, c(2,3),mean,na.rm=T)
Mon_simpSPS=apply(simpSPSarr, c(2,3),mean,na.rm=T)

saved_Date=Sys.Date()
save(file = Allsavename,version = 2,grd,calSPSarr,rawSPSarr,climaSPSarr,simpSPSarr,Mon_rawSPS,Mon_calSPS,Mon_climaSPS,Mon_simpSPS,saved_Date,inYR)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)),overwrite = T)
print(paste0("Done! saved file:",basename(Allsavename)))

