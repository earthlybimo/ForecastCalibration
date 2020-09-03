## Script for step 2 of Hannah Director remapping, stacking pred and obs into arrays and generating mapping
library(ncdf4)
library(IceCast)


Preddir="/pf/a/a270138/Data/BiasCorrOutput/Remapped_fcst"
Obsdir="/pf/a/a270138/Data/BiasCorrOutput/NSIDC_bootstrap_iceconc"
save_mapping_name="/pf/a/a270138/Data/BiasCorrOutput/Mapping_2003to2012_try1"
pred_ci=array(dim=c(10,304,448,12))
for(y in 3:12){
  infile=sprintf("%s/F%02d1_MON_CLIM.nc",Preddir,y)
  if(!file.exists(infile)){print(paste("File doesn't exist:",basename(infile)));next()}
  fl=nc_open(infile)
  ci=ncvar_get(fl,"a_ice")
  pred_ci[y-2,,,]=ci
  nc_close(fl)
} 
obs_ci=array(dim=c(10,304,448,12))
mm=2
for(y in 3:12){
  # seaice_conc_monthly_nh_f13_200402_v03r01.nc or seaice_conc_monthly_nh_f17_200902_v03r01.nc
  infile=Sys.glob(sprintf("%s/seaice_conc_monthly_nh*_20%02d%02d_v03r01.nc",Obsdir,y,mm))
  if(length(infile)==0){print(paste("File doesn't exist for year",y));next()}
  fl=nc_open(infile)
  ci=ncvar_get(fl,"seaice_conc_monthly_cdr")
  obs_ci[y-2,,,mm]=ci
  nc_close(fl)
}

pred_ci=aperm(pred_ci,c(1,4,2,3))
obs_ci=aperm(obs_ci,c(1,4,2,3))

obs <- get_region(dat = obs_ci[1,mm, ,],
                  dat_type = "simple", level = 0.15)
obs_map <- get_map(ice = obs, plotting = TRUE, reg_info,
                   main = "Observed Mapping \n September 2007")
month=mm
discrep <- create_mapping(start_year = 2003, end_year = 2012,
                          obs_start_year = 2003, pred_start_year = 2003,
                          observed = obs_ci[,mm,,],predicted = pred_ci[,month,,],
                          reg_info, month=mm,level = 0.15, dat_type_obs = "simple",
                          dat_type_pred = "simple", plotting = F)
save_mapping_name="/pf/a/a270138/Data/BiasCorrOutput/Mapping_2003to2012_try1"
save(file = save_mapping_name,version = 2,discrep,pred_ci,obs_ci,mm)