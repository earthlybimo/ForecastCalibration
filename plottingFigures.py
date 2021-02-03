#!/usr/bin/env python
# coding: utf-8
### Script for plotting figures to match Longjiang figures for the calibrated SIC forecasts

import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from netCDF4 import Dataset, num2date
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
#Could be something like this?
nod_NH= (mesh_lat>0.0)
nod_SH= (mesh_lat<0.0)

# osisaf climatology
RMSE_osisaf_cim_NH=np.zeros((2019-2011+1,12))
RMSE_osisaf_cim_SH=np.zeros((2019-2011+1,12))
RMSE_fcst_NH=np.zeros((2019-2011+1,4,12))
RMSE_fcst_SH=np.zeros((2019-2011+1,4,12))

year=2015
for year in np.arange(2011,2019+1):
    yr=str(year-2000)       # current year
    lyr=str(year-1-2000)    # last year

    file_osisaf = Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(year)+'.nc')
    file_osisaf_clim=Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_MON_CLIM_'+str(year-9)+'-'+str(year-1)+'.nc')
    osisaf = file_osisaf.variables['obs'][:]
    osisaf_clim=file_osisaf_clim.variables['obs'][:]
    file_osisaf.close()
    file_osisaf_clim.close()

    file_fcst = Dataset(save_path+'F'+yr+'_ens_mon_mean_corr.nc')
    fcst = file_fcst.variables['SIC_FCST_CORR'][:]
    file_fcst.close()

    for mon in np.arange(1,12+1):
        print(year,mon)

        # clim forecast
        diff = (osisaf_clim[mon-1,:] - osisaf[mon-1,:])**2
        diff[diff==0.0] = np.nan
        diff[np.abs(diff)>1.0] = np.nan
        areamask=mesh_area.copy()
        areamask[np.isnan(diff)]=np.nan
        diff=diff*areamask

        RMSE_osisaf_cim_NH[year-2011,mon-1] = np.sqrt(np.nansum(diff[nod_NH])/np.nansum(areamask[nod_NH]))
        RMSE_osisaf_cim_SH[year-2011,mon-1] = np.sqrt(np.nansum(diff[nod_SH])/np.nansum(areamask[nod_SH]))

            # forecast
        for lead in np.arange(4):
            diff = (fcst[lead,mon-1,:] - osisaf[mon-1,:])**2
            diff[diff==0.0] = np.nan
            diff[np.abs(diff)>1.0] = np.nan
            areamask=mesh_area.copy()
            areamask[np.isnan(diff)]=np.nan
            diff=diff*areamask

            RMSE_fcst_NH[year-2011,lead,mon-1] = np.sqrt(np.nansum(diff[nod_NH])/np.nansum(areamask[nod_NH]))
            RMSE_fcst_SH[year-2011,lead,mon-1] = np.sqrt(np.nansum(diff[nod_SH])/np.nansum(areamask[nod_SH]))



## Now we start plotting
marker = 'o'
markersize = [3,4,6]
plt.figure(figsize=(10,7))

plt.subplot(211)
RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_fcst_NH[year-2011,0,:])
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='blue',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='blue',label='L0-2')
plt.text(-4,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='blue')


RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_fcst_NH[year-2011,1,:])
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='orange',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='orange',label='L3-5')
plt.text(-4,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='orange')
#print(np.mean(RMSE_fcst_year))

RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_fcst_NH[year-2011,2,:])
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='green',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='green',label='L6-8')
plt.text(108,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='green')
#print(np.mean(RMSE_fcst_year))

RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_fcst_NH[year-2011,3,:])
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='red',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='red',label='L9-11')
plt.text(108,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='red')
#print(np.mean(RMSE_fcst_year))

RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_osisaf_cim_NH[year-2011,:])
plt.plot(RMSE_fcst_year,color='k',label='OSI SAF')
plt.text(108,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='k')
#print(np.mean(RMSE_fcst_year))

plt.legend(frameon=False,loc=2)
plt.xticks(np.arange(0,108,12),np.arange(2011,2020,1))
plt.ylabel('RMSE of sea ice concentration')
plt.title('Arctic')

plt.subplot(212)
RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_fcst_SH[year-2011,0,:])
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='C0',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='C0',label='L0-2')
plt.text(108,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='C0')
#print(np.mean(RMSE_fcst_year))

RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_fcst_SH[year-2011,1,:])
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='C1',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='C1',label='L3-5')
plt.text(108,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='C1')
#print(np.mean(RMSE_fcst_year))

RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_fcst_SH[year-2011,2,:])
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='C2',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='C2',label='L6-8')
plt.text(108,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='C2')
#print(np.mean(RMSE_fcst_year))

RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_fcst_SH[year-2011,3,:])
for leading in np.arange(0,3):
    plt.plot(np.arange(leading,108,3),RMSE_fcst_year[leading:108:3],color='C3',marker=marker,markersize=markersize[leading],alpha=0.5,markeredgecolor='None',linestyle='')
plt.plot(RMSE_fcst_year,color='C3',label='L9-11')
plt.text(108,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='C3')
#print(np.mean(RMSE_fcst_year))

RMSE_fcst_year = []
for year in np.arange(2011,2019+1):
    RMSE_fcst_year = np.append(RMSE_fcst_year, RMSE_osisaf_cim_SH[year-2011,:])
plt.plot(RMSE_fcst_year,color='k',label='OSI SAF')
plt.text(108,np.round(np.mean(RMSE_fcst_year),2),np.round(np.mean(RMSE_fcst_year),2),color='k')
#print(np.mean(RMSE_fcst_year))

plt.xticks(np.arange(0,108,12),np.arange(2011,2020,1))
plt.yticks(np.round(np.arange(0.08,0.3,0.05),2),np.round(np.arange(0.08,0.3,0.05),2))
plt.ylabel('RMSE of sea ice concentration')
plt.title('Antarctic')

plt.tight_layout()
plt.savefig(('/RMSE_Corrected.png'),dpi=300)
