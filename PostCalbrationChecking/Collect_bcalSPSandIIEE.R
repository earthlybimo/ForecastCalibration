
save_path = '/work/ba1138/a270138/BiasCorrOutput/bimoSimpResults'

Ylistt=2011:2018
IIEErawNH=array(dim=c(8,4,12))
IIEErawSH=array(dim=c(8,4,12))
IIEEbcalNH=array(dim=c(8,4,12))
IIEEbcalSH=array(dim=c(8,4,12))
SPSrawNH=array(dim=c(8,4,12))
SPSrawSH=array(dim=c(8,4,12))
SPSbcalNH=array(dim=c(8,4,12))
SPSbcalSH=array(dim=c(8,4,12))

for (yc in 1:8){ # yc=2;init=3
  
  for(init in 1:4){
    targetyear=Ylistt[yc]
    savename=paste0(save_path,"/SIMPbiasCorrect_bcalJ27_Yrs",targetyear,"_init_",init)
    if(!file.exists(savename)) {print(paste0("This doesn't exist?: ",basename(savename)))
    next()}
    load(savename,envir = (loadenv=new.env()))
    IIEErawNH[yc,init,]=loadenv$IIEErawNH
    IIEErawSH[yc,init,]=loadenv$IIEErawSH
    IIEEbcalNH[yc,init,]=loadenv$IIEEbcalNH
    IIEEbcalSH[yc,init,]=loadenv$IIEEbcalSH
    SPSrawNH[yc,init,]=loadenv$SPSrawNH
    SPSrawSH[yc,init,]=loadenv$SPSrawSH
    SPSbcalNH[yc,init,]=loadenv$SPSbcalNH
    SPSbcalSH[yc,init,]=loadenv$SPSbcalSH
    
  }}

Avg=list()
Avg$IIEErawNH=apply(IIEErawNH,3,mean,na.rm=T)
Avg$IIEErawSH=apply(IIEErawSH,3,mean,na.rm=T)
Avg$IIEEbcalNH=apply(IIEEbcalNH,3,mean,na.rm=T)
Avg$IIEEbcalSH=apply(IIEEbcalSH,3,mean,na.rm=T)

Avg$SPSrawNH=apply(SPSrawNH,3,mean,na.rm=T)
Avg$SPSrawSH=apply(SPSrawSH,3,mean,na.rm=T)
Avg$SPSbcalNH=apply(SPSbcalNH,3,mean,na.rm=T)
Avg$SPSbcalSH=apply(SPSbcalSH,3,mean,na.rm=T)

savedate=Sys.Date()
savename=paste0(save_path,"/SIMPbiasCorrect_bcalJ27_CollectedAll")
save(file = savename,version = 2,Avg,IIEEbcalSH,IIEEbcalNH,IIEErawSH,IIEErawNH,SPSbcalSH,SPSbcalNH,SPSrawSH,SPSrawNH,savedate,Ylistt)
print(paste0("Done! SPS IIEE collected file here: ",basename(savename)))

file.copy(from = savename,to = paste0("~/Data/tomove/",basename(savename)),overwrite = T)

