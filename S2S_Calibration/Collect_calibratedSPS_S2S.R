library(ncdf4);library(spheRlab)


data_path = '/work/ba1138/a270138/BiasCorrOutput/S2S_Results/'
modelname="ECMWF";flen=46
ylist=2005:2010
save_name=paste0(data_path,"CollectedCalibratedSPS_for_",modelname,"_wMask")

HEM="nh"

if(HEM=="nh"){
  maskfile="/mnt/lustre01/work/ab0995/a270138/S2S_Output/Newmask"  
  Gridfile="/mnt/lustre01/work/ab0995/a270138/S2S_Output/GridandNApoints"
  OGsource="/mnt/lustre01/work/ab0995/a270099/S2S/"
}

if(HEM=="sh"){
  outputdir="/mnt/lustre01/work/ab0995/a270138/S2S_Output_south/Newmask" 
  Gridfile="/mnt/lustre01/work/ab0995/a270138/S2S_Output_south/GridandNApoints"
  OGsource="/mnt/lustre01/work/ab0995/a270099/S2S_south/"
}

load(maskfile);mask=as.vector(mask)
cellareafile=paste0(OGsource,"mask_weight/area.nc")
fl=nc_open(cellareafile)
cellarea=ncvar_get(fl,'cell_area')
nc_close(fl)
cellarea=as.vector(cellarea)
## Find the grid areas
# Gridfile="/mnt/lustre01/work/ab0995/a270138/S2S_Output/GridandNApoints"
load(Gridfile,envir =(GridEnv=new.env()))
grd=GridEnv$grd.full
# Elementareas=array(dim=dim(grd$elem)[[1]])
# for (i in 1:length(Elementareas)){
#   nodes=grd$elem[i,]   #What nodes does this element have?
#   Elementareas[i]=sl.triag.area(grd$lon[nodes],grd$lat[nodes]) #And how much is the weight here?
# }
# Nodeareas=array(dim=length(grd$lon))
# for (i in 1:length(Nodeareas)){
#   elems=grd$neighelems[i,]  #Most nodes have upto 6 elements around it,
#   areaofelems=Elementareas[elems]
#   areaofelems=areaofelems[!is.na(areaofelems)]  #but some have less than 6, so some of those points could be NA
#   Nodeareas[i]=sum(areaofelems)*(1/3)  #Each element area is shared by 3 nodes, and nodes get shares from all the elements its part of.
# }

## Where to save the SPS values?
calSPSarr=array(dim=c(length(ylist),12,flen))
rawSPSarr=calSPSarr
for(yy in 1:length(ylist)){
  for(mm in 1:12){
    file1=Sys.glob(sprintf("%s%s/TAQM_calibrated_*%d-%0.2d-??.nc",data_path,modelname,ylist[yy],mm))
    if(length(file1)==0) next()
    
    fl=nc_open(file1)
    # fl$var$SIP_
    CalSIP=ncvar_get(fl,"SIP_FCST_CORR")
    ObsSIP=ncvar_get(fl,"SIP_FCST_OBS")
    RawSIP=ncvar_get(fl,"SIP_FCST_RAW")
    nc_close(fl)
    
    ObsSIP[ObsSIP<0.15]=0
    ObsSIP[ObsSIP>=0.15]=1
    CalSIP[CalSIP<0]=NA
    RawSIP[RawSIP<0]=NA
    # CalSIP[is.na(RawSIP)]=NA
    # RawSIP[is.na(CalSIP)]=NA
    
    # Now calculate the SPS
    for (i in 1:flen) {
      temp=as.vector(CalSIP[,,i]-ObsSIP[,,i])
      temp2=(temp^2)*cellarea*(mask)
      calSPSarr[yy,mm,i]=sum(temp2,na.rm = T)
      remove(temp,temp2)
      temp=as.vector(RawSIP[,,i]-ObsSIP[,,i])
      temp2=(temp^2)*cellarea*(mask)
      rawSPSarr[yy,mm,i]=sum(temp2,na.rm = T)
      remove(temp,temp2)
    }
  }
}

calMean=colMeans(calSPSarr)
rawMean=colMeans(rawSPSarr)
Fctr=((6371^2)/1000000)



save(file = save_name,version = 2, calSPSarr,rawSPSarr,ylist,calMean,rawMean,Fctr,flen)
print(paste0("File saved:",save_name))






