### Script to open calibrated SSIPS forecast, compute mean bias from estimated Mean, then save it into a table.

library(ncdf4);library(spheRlab);library(RColorBrewer)

save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults/MeanConc'
Allsavename=paste0("/work/ba1138/a270138/BiasCorrOutput/EstMeanConc_MeanBias")
Figpath="/pf/a/a270138/Data/BiasCorrOutput/Figs"

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
colBar=sl.colbar(brewer.pal(7,name = "BrBG"),10)
brks=seq(-0.8,0.8,length.out = 9)
inYR=2011:2018  #which year
# inMON=1:4   #And which initialisation 
strtm = c(1, 4, 7, 10)  # which is the starting month for each initialisation

grdlen=length(grd$lon)

rawBiasAvg=array(dim=c(4,12,grdlen)) # Bias averaged over the years
calBiasAvg=array(dim=c(4,12,grdlen)) # Bias averaged over the years

for(init in 1:4){
  for(mm in 1:12){ #leadtime in months
    rawBias_arr1=array(dim=c(grdlen,length(inYR)))
    calBias_arr1=array(dim=c(grdlen,length(inYR)))
    for(yy in 1:length(inYR)){      
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
      obsSIC=obsVar1[,obsTmnth]
      
      loadname=sprintf("%s/Forecast_Calibration_TrustSharpFalse_wMeanConc_Yr%d_%02dMn_%02d.nc",save_path,inYR[yy],init,mm)
      if(!file.exists(loadname)) next()
      
      fl=nc_open(loadname)
      rawSIC=ncvar_get(fl,"SIC_RAW") #Est mean 
      calSIC=ncvar_get(fl,"SIC_Cal")
      nc_close(fl)
      
      rawBias_arr1[,yy]=rawSIC-obsSIC
      calBias_arr1[,yy]=calSIC-obsSIC
      
      remove(rawSIC,calSIC,obsSIC)
    }
    
    ## Here we did not do na.rm=T in rowMeans, assuming that NA positions should not have changed. Not 100% sure about this.
    rawYavgBias=rowMeans(rawBias_arr1)
    calYavgBias=rowMeans(calBias_arr1)
    rawBiasAvg[init,mm,]= rawYavgBias 
    calBiasAvg[init,mm,]= calYavgBias
    
    ## Plot the average? maybe for only lead times of 3, 6, 9, 12
    if(mm %in% c(3,6,9,12)){
      #Raw average
      pir = sl.plot.init(projection="polar",polar.latbound = 50,file.name = sprintf("%s/meanBias_estRawMeanSIC_init_%s_lead%02d_months.pdf",Figpath,month.abb[strtm[init]],(mm-1)))
      # pir=sl.plot.init(projection = "polar",polar.latbound = 50,do.init.device = F)
      pcol= sl.plot.field.elem(pir,num=rawYavgBias,lon=grd$lon,lat = grd$lat,elem = grd$elem,colbar = colBar,colbar.breaks =brks)
      res = sl.plot.naturalearth(pir,fill.col = "grey",lines.col = "grey")
      # res=sl.plot.text(pir,lon=120,lat=65,labels = paste0("SPS = ",SPScal))
      sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
      sl.plot.end(pir)
      
      #And cal average
      pir = sl.plot.init(projection="polar",polar.latbound = 50,file.name = sprintf("%s/meanBias_estCalMeanSIC_init_%s_lead%02d_months.pdf",Figpath,month.abb[strtm[init]],(mm-1)))
      # pir=sl.plot.init(projection = "polar",polar.latbound = 50,do.init.device = F)
      pcol= sl.plot.field.elem(pir,num=calYavgBias,lon=grd$lon,lat = grd$lat,elem = grd$elem,colbar = colBar,colbar.breaks =brks)
      res = sl.plot.naturalearth(pir,fill.col = "grey",lines.col = "grey")
      # res=sl.plot.text(pir,lon=120,lat=65,labels = paste0("SPS = ",SPScal))
      sl.plot.lonlatgrid(pir, pole.hole=TRUE,col = "grey",labels = TRUE)
      sl.plot.end(pir)
    }
    
  } # leadtime loop
} # initialisation loop

sl.plot.colbar(file.name =paste0(Figpath,"/meanBias_colorscale.pdf"),breaks = brks,colbar = colBar) # labels.at = c(rep(c(T,F),4),T))
save(file = Allsavename,version = 2,grd,rawBiasAvg,calBiasAvg,inYR)
file.copy(from = Allsavename,to = paste0("~/Data/tomove/",basename(Allsavename)))
print("Done!")

