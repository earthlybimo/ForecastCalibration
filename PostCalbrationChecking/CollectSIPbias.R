### Script to open Simple calibrated SSIPS forecast, TAQM calibrated forecasts, measure bias, then save it all in a table.

args = commandArgs(trailingOnly=TRUE)
ccc=as.integer(args[1]) #Choice of Hemisphere


library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'
Clima_path= '/work/ab0995/a270112/data_fesom2/sic'
Data_path='/pf/a/a270138/Data/LongjiangOutput/FCST_CLIM'

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
Allsavename=paste0("/work/ba1138/a270138/BiasCorrOutput/Collected_SIPbias_all_",HEM)

inYR=2011:2018  #which year
# inMON=1:4   #And which initialisation 
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation

climaBiasArr=array(dim=c(length(grd$lat),length(inYR),4,12))
SimpCalBiasArr=climaBiasArr
RawBiasArr=climaBiasArr
TAQMcalBiasArr=climaBiasArr

for(yy in 1:length(inYR)){ # yy=1;init=1;mm=1 #For tests
  print(paste0("Looking at forecasts of year: ",inYR[yy]))
  loadname=sprintf("%s/F%0.2d_ens_mon_mean_SIP_corr.nc",Data_path,(inYR[yy]-2000))
  # if(!file.exists(loadname)) next()
  fl=nc_open(loadname)
  SIMPcalSIP_all=ncvar_get(fl,"SIP_FCST_CORR")
  # This assumes the simple forecasts are stored as SIP in  grid x mm x initialisation for each year
  nc_close(fl)
  
  for(init in 1:4){
    for(mm in 1:12){
      
      loadname=sprintf("%s/Forecast_Calibration_TrustSharpFalse_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
      if(!file.exists(loadname)) {print(paste0("File doesn't exist!: ",basename(loadname)))
        next()}
      
      fl=nc_open(loadname)
      rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
      TAQMcalSIP=ncvar_get(fl,"SIP_FCST_CORR")
      # obsSIP=ncvar_get(fl,"SIP_FCST_OBS")
      nc_close(fl)
      
      ## Get Obs
      obsTyr=inYR[yy]
      obsTmnth = mm+strtm[init]-1
      if (obsTmnth>12){
        obsTyr=obsTyr+1
        obsTmnth=obsTmnth-12}
      
      # print(paste0('Targetyear = ',(inYR[yy]),',initialisation = ',(init),' which means from ', (strtm[init]) ,',leadtime ',(mm),' so target month is ',(obsTmnth),' of year ',(obsTyr)))  # Testing
      
      loadname=paste0("/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_",obsTyr,".nc")
      if(!file.exists(loadname)) {print(paste0("Not found: ",basename(loadname)));next()}
      fl=nc_open(loadname)
      obsVar1=ncvar_get(fl,"obs") # num [1:126858, 1:12]
      nc_close(fl)
      obsVar2=obsVar1[,obsTmnth]
      obsSIP=array(dim =length(grd$lat))
      obsSIP[(obsVar2>=0.15)&(obsVar2<=1)]=1
      obsSIP[(obsVar2>=0)&(obsVar2<0.15)]=0
      
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
    
      climaBiasArr[,yy,init,mm]=climaSIP-obsSIP
      SIMPcalSIP=SIMPcalSIP_all[,mm,init]
      SimpCalBiasArr[,yy,init,mm]=SIMPcalSIP-obsSIP
      RawBiasArr[,yy,init,mm]=rawSIP-obsSIP
      TAQMcalBiasArr[,yy,init,mm]=TAQMcalSIP-obsSIP

    }
  }
}

SpecialNote="Bias is for SIP not SIC! Climatology was derived using 15% presence cutoff for last 9 years that month. Not sure if SIMP calibrated files are SIP and how they are computed."
save(file = Allsavename,version = 2,grd,climaBiasArr,SimpCalBiasArr,inYR,RawBiasArr,TAQMcalBiasArr,SpecialNote)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)),overwrite = T)
print(paste0("Done! saved file:",basename(Allsavename)))

