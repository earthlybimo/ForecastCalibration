#!/usr/bin/env python
# coding: utf-8
### Script for plotting IIEE figures to match Longjiang figures for the calibrated SIC forecasts

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

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
data_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults/MeanConc/LongStyle/'
Fig_path = '/work/ba1138/a270138/BiasCorrOutput/Figs/'
# CLIM_BIAS=readbin('CLIM_BIAS_04-10_4x12x126858f_method2.bin',[4,12,126858])
## Problem. where is readbin, and where is this bin file?

mesh_area_file=Dataset('/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc')
mesh_area=mesh_area_file.variables['cell_area'][:]
mesh_lat=mesh_area_file.variables['lat'][:]
mesh_lon=mesh_area_file.variables['lon'][:]
mesh_area_file.close()


IIEE_osisaf_cim_NH=np.zeros((2019-2011+1,12))
IIEE_osisaf_cim_SH=np.zeros((2019-2011+1,12))
IIEE_fcst_NH=np.zeros((2019-2011+1,4,12))
IIEE_fcst_SH=np.zeros((2019-2011+1,4,12))
year=2015;mon=4
for year in np.arange(2011,2019+1):
    yr=str(year-2000)       # current year
    lyr=str(year-1-2000)    # last year
    file_osisaf = Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(year)+'.nc')
    file_osisaf_clim=Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_MON_CLIM_'+str(year-9)+'-'+str(year-1)+'.nc')
    for mon in np.arange(1,12+1):

        print(year,mon)
        osisaf = file_osisaf.variables['obs'][mon-1,:]
        osisaf_clim=file_osisaf_clim.variables['obs'][mon-1,:]

        nod_mismatch = ((mesh_lat<-50.) & (osisaf > 0.15) & (osisaf <= 1.0) & (osisaf_clim < 0.15) & (osisaf_clim >= 0.)) | \
                ((mesh_lat<-50.) & (osisaf>=0.0) & (osisaf < 0.15) & (osisaf_clim > 0.15) & (osisaf_clim <= 1.))
        IIEE_osisaf_cim_SH[year-2011,mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12    # million km^2
        nod_mismatch = ((mesh_lat>20.) & (osisaf > 0.15) & (osisaf <= 1.0) & (osisaf_clim < 0.15) & (osisaf_clim >= 0.)) | \
                ((mesh_lat>20.) & (osisaf>=0.0) & (osisaf < 0.15) & (osisaf_clim > 0.15) & (osisaf_clim <= 1.))
        IIEE_osisaf_cim_NH[year-2011,mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12    # million km^2

        if (mon>=1 and mon<=3):
                    Fcst1=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 0, mon-1)
                    Fcst2=getFcst((data_path+'F'+lyr+'_ens_mon_mean_corr.nc'), 3, 3+mon-1)
                    Fcst3=getFcst((data_path+'F'+lyr+'_ens_mon_mean_corr.nc'), 2, 5+mon-1)
                    Fcst4=getFcst((data_path+'F'+lyr+'_ens_mon_mean_corr.nc'), 1, 8+mon-1)
        elif (mon>=4 and mon<=6):
                    Fcst1=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 1, mon-4)
                    Fcst2=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 0, 3+mon-4)
                    Fcst3=getFcst((data_path+'F'+lyr+'_ens_mon_mean_corr.nc'), 3, 5+mon-4)
                    Fcst4=getFcst((data_path+'F'+lyr+'_ens_mon_mean_corr.nc'), 2, 8+mon-4)

        elif (mon>=7 and mon<=9):
                    Fcst1=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 2, mon-7)
                    Fcst2=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 1, 3+mon-7)
                    Fcst3=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 0, 5+mon-7)
                    Fcst4=getFcst((data_path+'F'+lyr+'_ens_mon_mean_corr.nc'), 3, 8+mon-7)
        elif (mon>=10 and mon<=12):
                    Fcst1=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 3, mon-10)
                    Fcst2=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 2, 3+mon-10)
                    Fcst3=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 1, 5+mon-10)
                    Fcst4=getFcst((data_path+'F'+yr+'_ens_mon_mean_corr.nc'), 0, 8+mon-10)


        # tmp = fcst_bias_remove(get_fcst(year,mon,file1),CLIM_BIAS[0,mon-1,:])
        tmp = Fcst1
        nod_mismatch = ((mesh_lat<20.) & (osisaf > 0.15) & (osisaf <= 1.0) & (tmp < 0.15)) |  \
            ((mesh_lat<-50.) & (osisaf>=0.0) & (osisaf < 0.15) & (tmp > 0.15))
        IIEE_fcst_SH[year-2011, 0, mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12
        nod_mismatch = ((mesh_lat>20.) & (osisaf > 0.15) & (osisaf <= 1.0) & (tmp < 0.15)) |  \
            ((mesh_lat>20.) & (osisaf>=0.0) & (osisaf < 0.15) & (tmp > 0.15))
        IIEE_fcst_NH[year-2011, 0, mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12

        tmp = Fcst2
        nod_mismatch = ((mesh_lat<20.) & (osisaf > 0.15) & (osisaf <= 1.0) & (tmp < 0.15)) |  \
            ((mesh_lat<-50.) & (osisaf>=0.0) & (osisaf < 0.15) & (tmp > 0.15))
        IIEE_fcst_SH[year-2011, 1, mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12
        nod_mismatch = ((mesh_lat>20.) & (osisaf > 0.15) & (osisaf <= 1.0) & (tmp < 0.15)) |  \
            ((mesh_lat>20.) & (osisaf>=0.0) & (osisaf < 0.15) & (tmp > 0.15))
        IIEE_fcst_NH[year-2011, 1, mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12

        tmp = Fcst3
        nod_mismatch = ((mesh_lat<20.) & (osisaf > 0.15) & (osisaf <= 1.0) & (tmp < 0.15)) |  \
            ((mesh_lat<-50.) & (osisaf>=0.0) & (osisaf < 0.15) & (tmp > 0.15))
        IIEE_fcst_SH[year-2011, 2, mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12
        nod_mismatch = ((mesh_lat>20.) & (osisaf > 0.15) & (osisaf <= 1.0) & (tmp < 0.15)) |  \
            ((mesh_lat>20.) & (osisaf>=0.0) & (osisaf < 0.15) & (tmp > 0.15))
        IIEE_fcst_NH[year-2011, 2, mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12


        tmp = Fcst4
        nod_mismatch = ((mesh_lat<20.) & (osisaf > 0.15) & (osisaf <= 1.0) & (tmp < 0.15)) |  \
            ((mesh_lat<-50.) & (osisaf>=0.0) & (osisaf < 0.15) & (tmp > 0.15))
        IIEE_fcst_SH[year-2011, 3, mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12
        nod_mismatch = ((mesh_lat>20.) & (osisaf > 0.15) & (osisaf <= 1.0) & (tmp < 0.15)) |  \
            ((mesh_lat>20.) & (osisaf>=0.0) & (osisaf < 0.15) & (tmp > 0.15))
        IIEE_fcst_NH[year-2011, 3, mon-1] = np.sum(mesh_area[nod_mismatch])/1.e12

    file_osisaf.close()
file_osisaf_clim.close()


IIEE_fcst_yearmonmean_NH = np.mean(IIEE_fcst_NH,axis=0)
IIEE_fcst_yearmonmean_SH = np.mean(IIEE_fcst_SH,axis=0)
IIEE_osisaf_yearmonmean_NH = np.mean(IIEE_osisaf_cim_NH,axis=0)
IIEE_osisaf_yearmonmean_SH = np.mean(IIEE_osisaf_cim_SH,axis=0)


plt.figure(figsize=(8,6))
plt.subplot(211)
#plt.plot(IIEE_osisaf_yearmonmean,'s-k',label='OSI SAF')
for i in np.arange(4):
    if i==1:
        plt.plot(np.arange(3*i,3*i+3),IIEE_osisaf_yearmonmean_NH[3*i:3*i+3],'s-k',label='OSI SAF')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_NH[0,3*i:3*i+3],'.-',color='C0',label='L0-2')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_NH[1,3*i:3*i+3],'.-',color='C1',label='L3-5')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_NH[2,3*i:3*i+3],'.-',color='C2',label='L6-8')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_NH[3,3*i:3*i+3],'.-',color='C3',label='L9-11')
    else:
        plt.plot(np.arange(3*i,3*i+3),IIEE_osisaf_yearmonmean_NH[3*i:3*i+3],'s-k')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_NH[0,3*i:3*i+3],'.-',color='C0')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_NH[1,3*i:3*i+3],'.-',color='C1')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_NH[2,3*i:3*i+3],'.-',color='C2')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_NH[3,3*i:3*i+3],'.-',color='C3')
plt.legend(frameon=False,loc=2)
plt.xticks(np.arange(12),('J', 'F', 'M', 'A', 'M','J', 'J', 'A', 'S', 'O','N','D'),fontsize=12)
plt.ylabel(r'IIEE $(10^6 km^2)$',fontsize=12)
plt.title('Arctic')

plt.subplot(212)
#plt.plot(IIEE_osisaf_yearmonmean,'s-k',label='OSI SAF')
for i in np.arange(4):
    if i==1:
        plt.plot(np.arange(3*i,3*i+3),IIEE_osisaf_yearmonmean_SH[3*i:3*i+3],'s-k',label='OSI SAF')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_SH[0,3*i:3*i+3],'.-',color='C0',label='L0-2')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_SH[1,3*i:3*i+3],'.-',color='C1',label='L3-5')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_SH[2,3*i:3*i+3],'.-',color='C2',label='L6-8')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_SH[3,3*i:3*i+3],'.-',color='C3',label='L9-11')
    else:
        plt.plot(np.arange(3*i,3*i+3),IIEE_osisaf_yearmonmean_SH[3*i:3*i+3],'s-k')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_SH[0,3*i:3*i+3],'.-',color='C0')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_SH[1,3*i:3*i+3],'.-',color='C1')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_SH[2,3*i:3*i+3],'.-',color='C2')
        plt.plot(np.arange(3*i,3*i+3),IIEE_fcst_yearmonmean_SH[3,3*i:3*i+3],'.-',color='C3')
#plt.legend()
plt.xticks(np.arange(12),('J', 'F', 'M', 'A', 'M','J', 'J', 'A', 'S', 'O','N','D'),fontsize=12)
plt.ylabel(r'IIEE $(10^6 km^2)$',fontsize=12)
plt.title('Antarctic')
plt.tight_layout()

plt.savefig((Fig_path+'/IIEE_seasonalCycle_OBSvarCLIMBiasCorr_method2.png'),dpi=300)
