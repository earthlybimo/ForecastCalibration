## Script to interpolate Longjiang's model outputs to polar stereographic grid and start the bias correction process

system("module load cdo",wait = FALSE) #In case cdo is not loaded


Indir="/pf/a/a270138/Data/LongjiangOutput"
# For ensemble means, the output is:
Indir="/pf/a/a270138/Data/LongjiangOutput/FCST_CLIM"

Outdir="/pf/a/a270138/Data/BiasCorrOutput"

remapGrdfile=paste0(Outdir,"/grd_descp_polstr.txt")
FeSGrdDesc=paste0(Outdir,"/Fesom2_griddes.nc")
if(!file.exists(remapGrdfile)) {stop("Missing Polar stereo grid descp")}
if(!file.exists(FeSGrdDesc)) {stop("Missing Fesom2 grid descp")}

for(y in 3:19){
  infile=sprintf("%s/F%02d1_ens_mon_mean/F%02d1_MON_CLIM.nc",Indir,y,y)
  if(!file.exists(infile)) next()
  outfile=paste0(Outdir,"/Remapped_fcst/",basename(infile))
  if (file.exists(outfile)) next()
  system(paste0("cdo remapcon,",remapGrdfile," -setgrid,",FeSGrdDesc," ",infile," ",outfile),wait =T)
  print(paste0("Succesfully remapped: ",basename(infile)))
}

