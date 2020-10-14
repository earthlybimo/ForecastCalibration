library(ncdf4);library(spheRlab)

# fl=nc_open("~/Documents/Data/Others/LonjiangOutput/Fesom2_griddes.nc")
# lon=ncvar_get(fl,"lon")
# lat=ncvar_get(fl,"lon")
# lonBnds=ncvar_get(fl,"lon_bnds")
# latBnds=ncvar_get(fl,"lon_bnds")
# cellarea=ncvar_get(fl,"cell_area")
# nc_close(fl)
#Last time I did this: but now we are missing the elem (checked inside the grid desc, it's not there)
# Gridinfo=list(lon=lon,lat=lat,elem=elem,neighmat=neighmat,cell_area=cell_area)
# grd = sl.grid.curvilin2unstr(lon = lon, lat = lat, close.sides = TRUE)
## This above line gave a vector memory exhausted error. But I can read the grid directly from the file

grd = sl.grid.readNCDF("~/Documents/Data/Others/LonjiangOutput/Fesom2_griddes.nc")


fl=nc_open("~/Documents/Scripts/Others/Arlan_SICcalibration/SIC-probability-master/code/TEST_forecasts.nc")
rawSIP=ncvar_get(fl,"SIP_FCST_RAW")
calSIP=ncvar_get(fl,"SIP_FCST_CORR")
obsSIP=ncvar_get(fl,"SIP_FCST_OBS")
nc_close(fl)


preSPS=(obsSIP-rawSIP)^2   #Diff between model and satelite
preSPS2=preSPS*grd$cell_area
SPSraw=sum(preSPS2)*(10^-12)

preSPS=(calSIP-rawSIP)^2   #Diff between model and satelite
preSPS2=preSPS*grd$cell_area
SPScal=sum(preSPS2)*(10^-12)

# pir = sl.plot.init(projection="polar",polar.latbound = 50,do.init.device = F)
pir = sl.plot.init(projection="polar",polar.latbound = 55,file.name ="~/Documents/Figures/Oct/SSIPS_RawFcst_2013_01_leadtime_2.pdf")
pcol= sl.plot.field.elem(pir,num=rawSIP,lon=grd$lon,lat = grd$lat,elem = grd$elem)
res = sl.plot.naturalearth(pir)
cnt=sl.contours(var=obsSIP,lat = grd$lat,lon = grd$lon,elem=grd$elem,levels = 0.15)
res=sl.plot.contours(plot.init.res = pir,contours.res = cnt,col = "green")
sl.plot.text(pir,lon = 120,lat = 70,labels = paste0("SPS = ",round(SPSraw,4)," mill sq km"),cex=2)
sl.plot.end(pir)

pir = sl.plot.init(projection="polar",polar.latbound = 55,file.name ="~/Documents/Figures/Oct/SSIPS_CalFcst_2013_01_leadtime_2.pdf")
pcol= sl.plot.field.elem(pir,num=calSIP,lon=grd$lon,lat = grd$lat,elem = grd$elem)
res = sl.plot.naturalearth(pir)
# cnt=sl.contours(var=obsSIP,lat = grd$lat,lon = grd$lon,elem=grd$elem,levels = 0.15)
res=sl.plot.contours(plot.init.res = pir,contours.res = cnt,col = "green")
sl.plot.text(pir,lon = 120,lat = 70,labels = paste0("SPS = ",round(SPScal,4)," mill sq km"),cex=2)
sl.plot.end(pir)


