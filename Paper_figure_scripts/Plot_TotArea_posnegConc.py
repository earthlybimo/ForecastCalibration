import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from netCDF4 import Dataset, num2date

def getFcst(File1,init,mm):
    temp=Dataset(File1)
    temp2=temp.variables['SIC_FCST_CORR'][:]
    temp.close()
    temp3=temp2[init,mm,:]
    return temp3
def getObs(File1,init,mm):
    temp=Dataset(File1)
    temp2=temp.variables['SIC_FCST_CORR'][:]
    temp.close()
    temp3=temp2[init,mm,:]
    return temp3


# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults/MeanConc/LongStyle/'
Fig_path = '/work/ba1138/a270138/BiasCorrOutput/Figs/'
# CLIM_BIAS=readbin('CLIM_BIAS_04-10_4x12x126858f_method2.bin',[4,12,126858])
## Problem. where is readbin, and where is this bin file?

mesh_area_file=Dataset('/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc')
mesh_area=mesh_area_file.variables['cell_area'][:]
mesh_lat=mesh_area_file.variables['lat'][:]
mesh_lon=mesh_area_file.variables['lon'][:]
mesh_area_file.close()

# ?? where are these mesh.y2 coming from?
# nod_NH = (mesh.y2>40.)
# nod_SH = (mesh.y2<0.)

# From another file, but couldn't make it work
# sys.path.append("/pf/a/a270112/pyfesom2/")
# import pyfesom2 as pf
# meshpath='/work/ba1138/a270112/awicm3_input/fesom2/core2_meanz/'
# mesh=pf.load_mesh(meshpath, abg=[50, 15, -90], usepickle=False)

#Could be something like this?
# nod_NH= (mesh_lat>30.0)
# nod_SH= (mesh_lat<-30.0)

# osisaf climatology


fcst_anom_sum_pos_NH = np.zeros((4,108))
fcst_anom_sum_neg_NH = np.zeros((4,108))
fcst_anom_sum_pos_SH = np.zeros((4,108))
fcst_anom_sum_neg_SH = np.zeros((4,108))
osisaf_anom_sum_pos_NH = np.zeros(108)
osisaf_anom_sum_neg_NH = np.zeros(108)
osisaf_anom_sum_pos_SH = np.zeros(108)
osisaf_anom_sum_neg_SH = np.zeros(108)

months = 0
# osisaf climatology
data_path='/work/ba1138/a270112/awicm3/FCST_CLIM/'
file_osisaf_clim=Dataset(data_path+'OSISAF_MON_CLIM.nc')

for year in np.arange(2011,2018+1):
    yr=str(year-2000)       # current year
    lyr=str(year-1-2000)    # last year
    file_osisaf = Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(year)+'.nc')
    # file_anom = Dataset('F'+str(yr)+'_ens_mon_mean_corr_anom.nc')
    for mon in np.arange(1,12+1):
        print(year,mon)
        osisaf = file_osisaf.variables['obs'][mon-1,:]
        osisaf_clim=file_osisaf_clim.variables['obs'][mon-1,:]
        osisaf_anom = osisaf - osisaf_clim
        osisaf_anom = osisaf_anom * mesh_area
        nod_positive = np.where((osisaf_anom > 0) & (mesh_lat>40.) & (osisaf > .15)) #mesh.y2 replaced with mesh_lat
        nod_negative = np.where((osisaf_anom < 0) & (mesh_lat>40.) & (osisaf > .15))
        osisaf_anom_sum_pos_NH[months] = np.sum(osisaf_anom[nod_positive])
        osisaf_anom_sum_neg_NH[months] = np.sum(osisaf_anom[nod_negative])

        nod_positive = np.where((osisaf_anom > 0) & (mesh_lat<0.) & (osisaf > .15))
        nod_negative = np.where((osisaf_anom < 0) & (mesh_lat<0.) & (osisaf > .15))
        osisaf_anom_sum_pos_SH[months] = np.sum(osisaf_anom[nod_positive])
        osisaf_anom_sum_neg_SH[months] = np.sum(osisaf_anom[nod_negative])

        months = months + 1

    file_osisaf.close()
file_osisaf_clim.close()


strtm = (1, 4, 7, 10)  # which is the starting month for each initialisation
months = -1
file_osisaf_clim=Dataset(data_path+'OSISAF_MON_CLIM.nc')
for year in np.arange(2011,2019+1):
    yr=str(year-2000)       # current year

    file_fcst = Dataset(save_path+'F'+yr+'_ens_mon_mean_corr.nc')
    fcst = file_fcst.variables['SIC_FCST_CORR'][:]
    file_fcst.close()

    for mon in np.arange(1,12+1):
        months = months + 1
        for lead in np.arange(0,4):

            obsTmnth = mon+strtm[lead]-1
            obsTyr=year
            if (obsTmnth>12):
                obsTyr=obsTyr+1
                obsTmnth=obsTmnth-12

            ## Not sure what anomaly means, is it anomaly against obs or anomaly against clima?

            # truobsfile=('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(obsTyr)+'.nc')
            # ## Could add a file exists check here.
            # file_osisaf = Dataset(truobsfile)
            # truobs=file_osisaf.variables['obs'][obsTmnth-1,:]
            # file_osisaf.close()
            # diff = (fcst[lead,mon-1,:] - truobs)


            osisaf_clim=file_osisaf_clim.variables['obs'][obsTmnth-1,:]
            diff = (fcst[lead,mon-1,:] - osisaf_clim)

            fcst_anomAREA = diff * mesh_area

            nod_positive = np.where((fcst_anomAREA > 0) & (mesh_lat>40.))
            nod_negative = np.where((fcst_anomAREA < 0) & (mesh_lat>40.))
            fcst_anom_sum_pos_NH[lead,months] = np.sum(fcst_anomAREA[nod_positive])
            fcst_anom_sum_neg_NH[lead,months] = np.sum(fcst_anomAREA[nod_negative])

            nod_positive = np.where((fcst_anomAREA > 0) & (mesh_lat<0.))
            nod_negative = np.where((fcst_anomAREA < 0) & (mesh_lat<0.))
            fcst_anom_sum_pos_SH[lead,months] = np.sum(fcst_anomAREA[nod_positive])
            fcst_anom_sum_neg_SH[lead,months] = np.sum(fcst_anomAREA[nod_negative])





##

# Now we plot
marker = 'o'
markersize = [3,4,6]
plt.figure(figsize=(8,8))
plt.subplot(211)
RMSE_fcst_year = fcst_anom_sum_pos_NH[0,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='blue',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='blue',label='L0-2')
plt.text(-4,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='blue')

RMSE_fcst_year = fcst_anom_sum_pos_NH[1,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='orange',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='orange',label='L3-5')
plt.text(-4,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='orange')

RMSE_fcst_year = fcst_anom_sum_pos_NH[2,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='green',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='green',label='L6-8')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='green')

RMSE_fcst_year = fcst_anom_sum_pos_NH[3,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='red',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='red',label='L9-11')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='red')

RMSE_fcst_year = osisaf_anom_sum_pos_NH/1.e12
plt.plot(RMSE_fcst_year,color='k',label='OSI SAF')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='k')

RMSE_fcst_year = fcst_anom_sum_neg_NH[0,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='blue',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='blue')
plt.text(-4,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='blue')

RMSE_fcst_year = fcst_anom_sum_neg_NH[1,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='orange',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='orange')
plt.text(-4,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='orange')

RMSE_fcst_year = fcst_anom_sum_neg_NH[2,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='green',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='green')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='green')

RMSE_fcst_year = fcst_anom_sum_neg_NH[3,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='red',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='red')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='red')

RMSE_fcst_year = osisaf_anom_sum_neg_NH/1.e12
plt.plot(RMSE_fcst_year,color='k')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='k')

plt.legend(frameon=False,loc=3)
plt.xticks(np.arange(0,108,12),np.arange(2011,2020,1))
plt.xlim([-6,114])
plt.ylabel(r'Sea ice concentration anomaly $(10^6 km^2)$')
plt.title('Arctic')
#-------
plt.subplot(212)
RMSE_fcst_year = fcst_anom_sum_pos_SH[0,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='blue',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='blue',label='L0-2')
plt.text(-1,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='blue')

RMSE_fcst_year = fcst_anom_sum_pos_SH[1,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='orange',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='orange',label='L3-5')
plt.text(-4,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='orange')

RMSE_fcst_year = fcst_anom_sum_pos_SH[2,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='green',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='green',label='L6-8')
plt.text(111,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='green')

RMSE_fcst_year = fcst_anom_sum_pos_SH[3,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='red',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='red',label='L9-11')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='red')

RMSE_fcst_year = osisaf_anom_sum_pos_SH/1.e12
plt.plot(RMSE_fcst_year,color='k',label='OSI SAF')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='k')

RMSE_fcst_year = fcst_anom_sum_neg_SH[0,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='blue',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='blue')
plt.text(-1,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='blue')

RMSE_fcst_year = fcst_anom_sum_neg_SH[1,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='orange',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='orange')
plt.text(-4,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='orange')

RMSE_fcst_year = fcst_anom_sum_neg_SH[2,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='green',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='green')
plt.text(111,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='green')

RMSE_fcst_year = fcst_anom_sum_neg_SH[3,:]/1.e12
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='red',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='red')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='red')

RMSE_fcst_year = osisaf_anom_sum_neg_SH/1.e12
plt.plot(RMSE_fcst_year,color='k')
plt.text(108,np.round(np.mean(RMSE_fcst_year),1),np.round(np.mean(RMSE_fcst_year),1),color='k')

plt.xticks(np.arange(0,108,12),np.arange(2011,2020,1))
plt.ylabel(r'Sea ice concentration anomaly $(10^6 km^2)$')

plt.title('Antarctic')
plt.xlim([-6,114])
plt.tight_layout()

plt.savefig((Fig_path+'SIC_AnomalyPosandNeg.png'),dpi=300)
