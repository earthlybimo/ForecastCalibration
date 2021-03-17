### Script to collect IceExtent for all ECMWF forecasts

HEM="nh"
MODELS=c("UKMO","CMA","KMA","MFRANCE","ECMWF","ECMWF_PRES","NCEP")
# FLENGTHS=c(60,60,60,59,46,46,44) #This doesn't matter now, since we use forecast length
modelname=MODELS[5]
flen=46

library(ncdf4)
getncvar<-function(file1,var1){ #Function to make taking variables out of ncdf easy
  fl=nc_open(file1)
  v1=ncvar_get(fl,var1)
  nc_close(fl)
  return(v1)
}
binarise <-function (somearr,dlevel) {  #Function to binarise some given array  based on  this level, greater than this is 1, rest is 0.
  ll=dim(somearr)
  if (is.null(ll)){ll=length(somearr)}  #because 1d arrays dont have dim, they have length
  ibinar=array(dim =ll)
  ibinar[somearr[]>=dlevel]=1
  ibinar[somearr[]<dlevel]=0
  return (ibinar)
}

if(HEM=="nh"){
  datasource="/mnt/lustre01/work/ab0995/a270099/S2S/"  #In Lorenzos' side
  newmaskfile="/mnt/lustre01/work/ab0995/a270138/S2S_Output/Newmask"
}
if(HEM=="sh"){
  datasource="/mnt/lustre01/work/ab0995/a270099/S2S_south/"  #In Lorenzos' side
  newmaskfile="/mnt/lustre01/work/ab0995/a270138/S2S_Output_south/Newmask"
}
outputdir="/work/ba1138/a270138/BiasCorrOutput/S2S_Results"

# maskfile=paste0(datasource,"mask_weight/mask.nc")
# OGmask=getncvar(maskfile,'ci') # This is the mask that lorenzo used, in case we want to use that one

## Now we are using a new mask!
load(newmaskfile)

cellareafile=paste0(datasource,"mask_weight/area.nc")
cellarea=getncvar(cellareafile,'cell_area')

sat_folder=paste0(datasource,"satellite_interp")
Modelpath=paste0(datasource,"forecasts/",modelname)

modelArea_big=array(0,dim = c(12,12,flen))  #12 yrs, 12 mnths, 46 leadtime
satArea_big=modelArea_big
modelArea_small=modelArea_big
satArea_small=modelArea_big
dateArr=array(0,dim = c(12,12))

Ylist=1999:2010
for (yc in 1:length(Ylist)){
  yy=Ylist[yc];print(paste0('Saving for year: ',yy))
  for(mm in 1:12){
    donecount=0
    for(dd in 1:8){
      if(donecount>0) next() #If already did a file for this month, skip rest.
      start_date=format.Date(sprintf("%d-%02d-%02d",yy,mm,dd)) #Cycling through start date
      
      pert_forc=Sys.glob(paste0(Modelpath,"/",modelname,"_ref_pert_sic_",start_date,".nc"))
      cont_forc=Sys.glob(paste0(Modelpath,"/",modelname,"_ref_cont_sic_",start_date,".nc"))
      #Check file exists? If not skip?
      if(length(pert_forc)==0) next()
      start_timestep=as.Date(start_date)
      dateArr[yc,mm]=start_date
      if (((modelname == "MFRANCE") & (HEM=="nh")) | (modelname == "NCEP")) start_timestep= start_timestep+1
      # end_timestep=start_timestep+flength-1
      ## Because NCEP and MFrance (but only in NH) start from day +1, so we use start_timestep as the real timestep of comparison 
      
      ###### Part 1: merge model file, binarise 
      ci_p=getncvar(pert_forc,'ci')
      ci_c=getncvar(cont_forc,'ci')
      dm=dim(ci_p)
      ci=array(dim=c(dm[1],dm[2],(dm[3]+1),dm[4]))
      flength=dm[4]
      ci[,,1:dm[3],]=ci_p
      ci[,,(dm[3]+1),]=ci_c
      ci2=binarise(ci,0.15)
      
      Area1=apply(ci2, c(3,4),function(x) sum(x*cellarea,na.rm = T))
      modelArea_big[yc,mm,]=colMeans(Area1)  #If we later wanted to save the whole ensemble, we could do that here
      
      Area1=apply(ci2, c(3,4),function(x) sum(x*cellarea*mask,na.rm = T))
      modelArea_small[yc,mm,]=colMeans(Area1)  #If we later wanted to save the whole ensemble, we could do that here
      
      
      
      ###### Part 2: Find right Satelite Data and merge them ~~~~~~~~~~~
      
      #Gather all satelite files from this daterange
      sat_sip=array(dim=c(dm[1],dm[2],dm[4]))
      for(f in 0:(flength-1)) { #Along time, starting from initial day to final day
        satfile=Sys.glob(paste0(sat_folder,"/ice_conc_",HEM,"_ease-125_reproc_",format((start_timestep+f),"%Y%m%d"),"1200.nc4"))
        if(length(satfile)==0) next()
        sat_ci=getncvar(satfile,'ice_conc')
        sat_sip[,,(f+1)]=binarise(sat_ci,15)
      }
      Area1=apply(sat_sip, 3,function(x) sum(x*cellarea,na.rm = T))
      satArea_big[yc,mm,]=(Area1) 
      
      Area1=apply(sat_sip, 3,function(x) sum(x*cellarea*mask,na.rm = T))
      satArea_small[yc,mm,]=(Area1) 
      
      remove(sat_sip,ci_c,ci_p,ci2,ci,pert_forc,cont_forc,Area1)
      invisible(gc(verbose = FALSE))
      donecount=donecount+1 
    }
  }}

savename=paste0(outputdir,"/IceExtent_",modelname,"_",HEM)
save(file = savename,version = 2,cellarea,mask,satArea_small,satArea_big,modelArea_small,modelArea_big,flen,dm,Ylist,dateArr)
file.copy(savename,to = paste0("~/Data/tomove/",basename(savename)))
print(paste0("File ",basename(savename) ," saved"))