#! R
args = commandArgs(trailingOnly=TRUE)
# hmm=as.integer(args[1]) #Which HEM?  #For now, just NH  
yc=as.integer(args[2]) # Target year from 2011 to 2018
init=as.integer(args[3]) # Which initialisation? 1 to 4
mm=as.integer(args[4]) # lead time

Ylistt=2011:2018
targetyear=Ylistt[yc]
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation



library(ncdf4);library(spheRlab);library(RColorBrewer)
Figpath="/pf/a/a270138/Data/BiasCorrOutput/Figs/"

Clima_path= '/work/ab0995/a270112/data_fesom2/sic'
Data_path='/pf/a/a270138/Data/LongjiangOutput/FCST_CLIM' #Raw fcst
TAQMpath = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults'  #TAQM calibrated
SIMPbpath='/work/ba1138/a270138/BiasCorrOutput/bimoSimpResults'  #Simple calibd

grd = sl.grid.readNCDF("/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc")
grlen=length(grd$lat)
binarise <-function (somearr,dlevel) {  #Function to binarise some given array  based on  this level
  ll=dim(somearr)
  if (is.null(ll)){ll=length(somearr)}  #because 1d arrays dont do
  ibinar=array(dim =ll)
  ibinar[somearr[]>=dlevel]=1
  ibinar[somearr[]<dlevel]=0
  return (ibinar)
}

### Get variables

# 1. Obs
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
obsSIP=array(dim =length(obsVar2))
obsSIP[(obsVar2>=0.15)&(obsVar2<=1)]=1
obsSIP[(obsVar2>=0)&(obsVar2<0.15)]=0


# 2. Raw Fcst:
rawFcst=array(dim=c(30,grlen,12))
for(EM in 1:30){
  file_anom0=sprintf("%s/F%0.2d%0.1d_ens_mon_mean/SIC_mon_%0.2d.nc",Data_path,(targetyear-2000),init,EM)
  if(!file.exists(file_anom0)) next()
  
  fl=nc_open(file_anom0)
  sic0=ncvar_get(fl,'a_ice') #grlen x 12 leadtimes
  nc_close(fl)
  rawFcst[EM,,]=sic0
}
sip0=binarise(rawFcst[,,mm],0.15)
RawSIP=apply(sip0,2,mean,na.rm=T) #averaging over EM

#3. TAQM Fcst:
loadname=sprintf("%s/Forecast_Calibration_TrustSharpFalse_Yr%d_%02dMn_%02d.nc",TAQMpath,targetyear,init,mm)
if(!file.exists(loadname)) next()
fl=nc_open(loadname)
# rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
TAQMcalSIP=ncvar_get(fl,"SIP_FCST_CORR")
# obsSIP=ncvar_get(fl,"SIP_FCST_OBS")
nc_close(fl)

#4. SIMPb Fcst:
savename=paste0(SIMPbpath,"/SIMPbiasCorrect_bcalJ27_Yrs",targetyear,"_init_",init)
if(!file.exists(savename)) next()
load(savename,envir = (loadenv=new.env()))
temp=loadenv$bcalFcst[,,mm]
SIMPbcalSIP=apply(temp,2,mean,na.rm=T) #averaging over EM 



#### Let's plot
testcol=sl.colbar(brewer.pal(9,"YlGnBu"),11)
brks=seq(0,0.9,length.out = 10)

# Will make a 4 x 4 matrix for each instance
figname=paste0(Figpath,"SIP_4flavors_yr",targetyear,"_init",init,"_ltime_",mm,".pdf")
pdf(file = figname,width = 20,height = 20)
par(mfrow=c(2,2))
pir = sl.plot.init(projection="polar",polar.latbound = 55,do.init.device = F,main = "Obs SIP")
res = sl.plot.naturalearth(pir)
pcol= sl.plot.field.elem(pir,num=obsSIP,lon=grd$lon,lat = grd$lat,elem = grd$elem,colbar =testcol,colbar.breaks = brks)
sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
sl.plot.end(pir,do.close.device = F)

# par()
pir = sl.plot.init(projection="polar",polar.latbound = 55,do.init.device = F,main = "Raw SIP")
 
res = sl.plot.naturalearth(pir)
pcol= sl.plot.field.elem(pir,num=RawSIP,lon=grd$lon,lat = grd$lat,elem = grd$elem,colbar =testcol,colbar.breaks = brks)
sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
sl.plot.end(pir,do.close.device = F)

# par()
pir = sl.plot.init(projection="polar",polar.latbound = 55,do.init.device = F,main = "TQAM SIP")
res = sl.plot.naturalearth(pir)
pcol= sl.plot.field.elem(pir,num=TAQMcalSIP,lon=grd$lon,lat = grd$lat,elem = grd$elem,colbar =testcol,colbar.breaks = brks)
sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
sl.plot.end(pir,do.close.device = F)

# par()
pir = sl.plot.init(projection="polar",polar.latbound = 55,do.init.device = F,main = "SIMP-B SIP")
res = sl.plot.naturalearth(pir)
pcol= sl.plot.field.elem(pir,num=SIMPbcalSIP,lon=grd$lon,lat = grd$lat,elem = grd$elem,colbar =testcol,colbar.breaks = brks)
sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
sl.plot.end(pir,do.close.device = F)

dev.off()


file.copy(from = figname,to = paste0("~/Data/tomove/",basename(figname)),overwrite = T)
print(paste0("Done! saved figure:",basename(figname)))