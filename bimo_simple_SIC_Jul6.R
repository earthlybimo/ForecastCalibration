### Script to go through each EM, find raw bias for each year, correct SIC directly, and measure SPS

args = commandArgs(trailingOnly=TRUE)
yc=as.integer(args[1]) # Target year from 2011 to 2018
Ylistt=2011:2018
targetyear=Ylistt[yc]

init=as.integer(args[2]) # Which initialisation? 1 to 4

library(ncdf4);library(spheRlab)

save_path = '/work/ba1138/a270138/BiasCorrOutput/SIMPB_SIC_results'
if(!dir.exists(save_path)) dir.create(save_path,recursive = T)
Clima_path= '/work/ab0995/a270112/data_fesom2/sic'
Data_path='/pf/a/a270138/Data/LongjiangOutput/FCST_CLIM'
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation


savename=paste0(save_path,"/BSIMPbiasCorrectSIC_J27_Yrs",targetyear,"_init_",init)
if(file.exists(savename)) stop("Already done?")
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
NH=which(grd$lat>0)
SH=which(grd$lat<0)

# targetyear=2014 #measurement years
histYrs=2003:(targetyear-1)
grlen=126858 # 


histbiasAll=array(dim=c(length(histYrs),grlen,12))
# Collect bias over all years in histYrs
for( y1 in 1:length(histYrs)){
  whichyear= histYrs[y1]
  
  ## Get historical obs fr this period:
  histObs=array(dim=c(grlen,12))
  for(mm in 1:12){
    obsTyr=whichyear
    obsTmnth = mm+strtm[init]-1
    if (obsTmnth>12){
      obsTyr=whichyear+1
      obsTmnth=obsTmnth-12}
    # print(paste0('Targetyear = ',(inYR[yy]),',initialisation = ',(init),' which means from ', (strtm[init]) ,',leadtime ',(mm),' so target month is ',(obsTmnth),' of year ',(obsTyr)))  # Testing
    
    loadname=paste0("/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_",obsTyr,".nc")
    if(!file.exists(loadname)) {print(paste0("Not found: ",basename(loadname)));next()}
    fl=nc_open(loadname)
    obsVar1=ncvar_get(fl,"obs") # num [1:126858, 1:12]
    nc_close(fl)
    obsVar2=obsVar1[,obsTmnth]
    histObs[,mm]=obsVar2
  }
  
  ## Get historical foorecast and bias
  histFcst=array(dim=c(30,grlen,12))
  histbias=array(dim=c(30,grlen,12))
  for(EM in 1:30){
    file_anom0=sprintf("%s/F%0.2d%0.1d_ens_mon_mean/SIC_mon_%0.2d.nc",Data_path,(whichyear-2000),init,EM)
    # print(file_anom0)
    if(!file.exists(file_anom0)) next()
    
    fl=nc_open(file_anom0)
    sic0=ncvar_get(fl,'a_ice') #grlen x 12 leadtimes
    nc_close(fl)
    histFcst[EM,,]=sic0
    histbias[EM,,]=sic0-histObs
  }
  
  histbiasAll[y1,,]=apply(histbias,c(2,3),mean,na.rm=T) #averaging over EM
}

meanHistbias=apply(histbiasAll,c(2,3),mean,na.rm=T)  #averaging over all years in histYrs

### Get raw forecast and also calibrated forecast:
rawFcst=array(dim=c(30,grlen,12))
bcalFcst=array(dim=c(30,grlen,12))
for(EM in 1:30){
  file_anom0=sprintf("%s/F%0.2d%0.1d_ens_mon_mean/SIC_mon_%0.2d.nc",Data_path,(targetyear-2000),init,EM)
  # print(file_anom0)
  if(!file.exists(file_anom0)) next()
  
  fl=nc_open(file_anom0)
  sic0=ncvar_get(fl,'a_ice') #grlen x 12 leadtimes
  nc_close(fl)
  rawFcst[EM,,]=sic0
  bcalFcst[EM,,]=sic0-meanHistbias
}
bcalFcst_prelimit=bcalFcst
bcalFcst[bcalFcst>1]=1
bcalFcst[bcalFcst<0]=0

## To find SPS, binarise, average 'cal fcst' over EM to get Probability
binbcal=binarise(bcalFcst,0.15)
binraw=binarise(rawFcst,0.15)
bcalSIP=apply(binbcal,c(2,3),mean,na.rm=T)
rawSIP=apply(binraw,c(2,3),mean,na.rm=T)

SPSbcalNH=array(dim=c(12))
SPSrawNH=array(dim=c(12))
SPSbcalSH=array(dim=c(12))
SPSrawSH=array(dim=c(12))
IIEEbcalNH=array(dim=c(12))
IIEErawNH=array(dim=c(12))
IIEEbcalSH=array(dim=c(12))
IIEErawSH=array(dim=c(12))

## Get target obs  and IIEE, SPS at the same time:
targetObs=array(dim=c(grlen,12))
for(mm in 1:12){
  obsTyr=targetyear
  obsTmnth = mm+strtm[init]-1
  if (obsTmnth>12){
    obsTyr=targetyear+1
    obsTmnth=obsTmnth-12}
  # print(paste0('Targetyear = ',(inYR[yy]),',initialisation = ',(init),' which means from ', (strtm[init]) ,',leadtime ',(mm),' so target month is ',(obsTmnth),' of year ',(obsTyr)))  # Testing
  
  loadname=paste0("/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_",obsTyr,".nc")
  if(!file.exists(loadname)) {print(paste0("Not found: ",basename(loadname)));next()}
  fl=nc_open(loadname)
  obsVar1=ncvar_get(fl,"obs") # num [1:126858, 1:12]
  nc_close(fl)
  obsVar2=obsVar1[,obsTmnth]
  obsSIP=array(dim =length(grlen))
  obsSIP[(obsVar2>=0.15)&(obsVar2<=1)]=1
  obsSIP[(obsVar2>=0)&(obsVar2<0.15)]=0
  targetObs[,mm]=obsSIP
  
  
  preSPS0=grd$cell_area*((bcalSIP[,mm]-obsSIP)^2)   #Diff between model and satelite
  SPSbcalNH[mm]=sum(preSPS0[NH],na.rm = T)
  SPSbcalSH[mm]=sum(preSPS0[SH],na.rm = T)
  preSPS0=grd$cell_area*((rawSIP[,mm]-obsSIP)^2)   #Diff between model and satelite
  SPSrawNH[mm]=sum(preSPS0[NH],na.rm = T)
  SPSrawSH[mm]=sum(preSPS0[SH],na.rm = T)
  
  preIIEE0=grd$cell_area*((binarise(bcalSIP[,mm],0.5)-obsSIP)^2)   #Diff between model and satelite
  IIEEbcalNH[mm]=sum(preIIEE0[NH],na.rm = T)
  IIEEbcalSH[mm]=sum(preIIEE0[SH],na.rm = T)
  preIIEE0=grd$cell_area*((binarise(rawSIP[,mm],0.5)-obsSIP)^2)   #Diff between model and satelite
  IIEErawNH[mm]=sum(preIIEE0[NH],na.rm = T)
  IIEErawSH[mm]=sum(preIIEE0[SH],na.rm = T)
}


save(file = savename,version = 2,SPSrawSH,SPSrawNH,SPSbcalSH,SPSbcalNH,IIEErawSH,IIEEbcalSH,IIEErawNH,IIEEbcalNH,bcalFcst,bcalFcst_prelimit,meanHistbias)
print(paste0("Done! Bias removed and saved SIC (and SPS/IIEE) here: ",basename(savename)))

