### Script to open Simple calibrated SSIPS forecast, compute IIEE, then save it into a table.

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
Allsavename=paste0("/work/ba1138/a270138/BiasCorrOutput/CollectedIIEE_LongCalandClima_",HEM)

inYR=2011:2018  #which year
# inMON=1:4   #And which initialisation 
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation

climaOarr=array(dim=c(length(inYR),4,12))
climaUarr=climaOarr
LcalOarr=climaOarr
LcalUarr=climaUarr

for(yy in 1:length(inYR)){ # yy=1;init=1;mm=1 #For tests
  loadname=sprintf("%s/F%0.2d_ens_mon_mean_SIP_corr.nc",Data_path,(inYR[yy]-2000))
  # if(!file.exists(loadname)) next()
  fl=nc_open(loadname)
  LcalSIP_all=ncvar_get(fl,"SIP_FCST_CORR")
  nc_close(fl)
  
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
      obsSIP[(obsVar2>=0)&(obsVar2<0.15)]=0
      
      ## Climatology:
      ClimFile=paste0(Clima_path,"/OSISAF_MON_CLIM_",(obsTyr-9),"-",(obsTyr-1),".nc")
      if(!file.exists(ClimFile)) {print(paste0("This Climafile doesnt exist! ",basename(ClimFile))); next()}
      fl=nc_open(ClimFile)
      climaSIC=ncvar_get(fl,"obs")
      nc_close(fl)
      climaSIPmed=binarise(climaSIC[,obsTmnth],0.15)
      
      uid=which((climaSIPmed==0)&(obsSIP==1)&(Hemfilter==1)) #identify underforecast grids
      U=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      oid=which((climaSIPmed==1)&(obsSIP==0)&(Hemfilter==1)) #identify overforecast grids
      O=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      climaUarr[yy,init,mm]=U*(10^-12)
      climaOarr[yy,init,mm]=O*(10^-12)
      
      remove(O,U,uid,oid,climaSIPmed,climaSIC)
      
      LcalSIP=LcalSIP_all[,mm,init]
      LcalSIPmed=binarise(LcalSIP,0.5)
      uid=which((LcalSIPmed==0)&(obsSIP==1)&(Hemfilter==1)) #identify underforecast grids
      U=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      oid=which((LcalSIPmed==1)&(obsSIP==0)&(Hemfilter==1)) #identify overforecast grids
      O=sum(grd$cell_area[uid],na.rm=T)  #And sum area
      LcalUarr[yy,init,mm]=U*(10^-12)
      LcalOarr[yy,init,mm]=O*(10^-12)
      
      remove(O,U,uid,oid,LcalSIPmed,LcalSIP)
      
    }
  }
}

climaIIEEarr=climaOarr+climaUarr
LcalIIEEarr=LcalOarr+LcalUarr
Mon_climaIIEE=apply(climaIIEEarr, c(2,3),mean,na.rm=T)
Mon_LcalIIEE=apply(LcalIIEEarr, c(2,3),mean,na.rm=T)


# climaAEEarr=abs(climaOarr- climaUarr )
# climaMEarr=climaIIEEarr-climaAEEarr
# LcalAEEarr=abs(LcalOarr-LcalUarr)
# LcalMEarr=LcalIIEEarr-LcalAEEarr

save(file = Allsavename,version = 2,grd,LcalIIEEarr,climaIIEEarr,inYR,Mon_climaIIEE,Mon_LcalIIEE,climaOarr,LcalOarr)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)),overwrite = T)
print(paste0("Done! saved file:",basename(Allsavename)))

