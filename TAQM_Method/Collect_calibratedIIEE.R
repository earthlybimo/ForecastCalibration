### Script to open calibrated SSIPS forecast, compute IIEE, then save it into a table.

args = commandArgs(trailingOnly=TRUE)
ccc=as.integer(args[1]) #Choice of Hemisphere


library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'


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

for(yy in 1:length(inYR)){ # yy=1;init=1;mm=1 #For tests
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
      obsVar2=obsVar1[,obsTmnth]
      obsSIP=array(dim =length(grd$lat))
      obsSIP[(obsVar2>=0.15)&(obsVar2<=1)]=1
      obsSIP[(obsVar2>=0)&(obsVar2<0.15)]=1
      
      loadname=sprintf("%s/Forecast_Calibration_TrustSharpFalse_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
      if(!file.exists(loadname)) next()
      
      fl=nc_open(loadname)
      rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
      calSIP=ncvar_get(fl,"SIP_FCST_CORR")
      # obsSIP2=ncvar_get(fl,"SIP_FCST_OBS")
      nc_close(fl)
      
      rawSIPmed=binarise(rawSIP,0.5)
      uid=which((rawSIP==0)&(obsSIP==1)&(Hemfilter==1)) #identify underforecast grids
      U=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      oid=which((rawSIP==1)&(obsSIP==0)&(Hemfilter==1)) #identify overforecast grids
      O=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      rawUarr[yy,init,mm]=U*(10^-12)
      rawOarr[yy,init,mm]=O*(10^-12)
      
      remove(O,U,uid,oid,rawSIPmed)
      
      calSIPmed=binarise(calSIP,0.5)
      uid=which((calSIP==0)&(obsSIP==1)&(Hemfilter==1)) #identify underforecast grids
      U=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      oid=which((calSIP==1)&(obsSIP==0)&(Hemfilter==1)) #identify overforecast grids
      O=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      calUarr[yy,init,mm]=U*(10^-12)
      calOarr[yy,init,mm]=O*(10^-12)
      
      remove(O,U,uid,oid,calSIPmed)
      
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
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)),overwrite = T)
print(paste0("Done! saved file:",basename(Allsavename)))

