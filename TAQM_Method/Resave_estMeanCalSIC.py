#!/usr/bin/env python
# coding: utf-8
### Script for taking calibrated mean SIC forecasts and saving them in arrangement of Fcst[leadtime, mm, grid] to use with other scripts from Longjiang.

import os
import sys
import numpy as np
from netCDF4 import Dataset, num2date
# data_path ='/work/ba1138/a270112/awicm3/FCST_CLIM/'
data_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults/MeanConc/'
Grdlen = 126858  # For now we are pre-setting these. Dims[3] = 126858

for targetyear in np.arange(2011,(2019+1)):
# targetyear=2014
    SIC_cal_arr=np.empty((4,12,Grdlen))
    SIC_raw_arr=np.empty((4,12,Grdlen))
    for whichinit in np.arange(1,5):
        for leadtimeMonth in np.arange(1,13):
            filename = data_path+'Forecast_Calibration_TrustSharpFalse_wMeanConc_Yr'+str(targetyear)+'_'+str(whichinit).zfill(2)+'Mn_'+str(leadtimeMonth).zfill(2)+'.nc'
            fileconcMean=Dataset(filename)
            calMeanSIC=fileconcMean.variables['SIC_Cal'][:]
            rawMeanSIC=fileconcMean.variables['SIC_RAW'][:]
            fileconcMean.close()
            SIC_cal_arr[(whichinit-1),(leadtimeMonth-1),:]=calMeanSIC
            SIC_raw_arr[(whichinit-1),(leadtimeMonth-1),:]=rawMeanSIC



    # Name for the variable later = SIC_FCST_CORRs

    # = Dataset('FCST_CLIM/F'+yr+'_ens_mon_mean_corr.nc')
    save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults/MeanConc/LongStyle/'
    yr=str(targetyear-2000)

    filename = (save_path+'F'+str(yr)+'_ens_mon_mean_corr.nc')
    ncfile = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
    # create dimensions
    ncfile.createDimension('n2d', 126858)
    ncfile.createDimension('fcst',4)
    ncfile.createDimension('time',12)
    # ncfile.createDimension('ens',30)
    # define variables
    Fcst_calsic_to_write = ncfile.createVariable('SIC_FCST_CORR', 'd', ('fcst','time','n2d'))
    Fcst_rawsic_to_write = ncfile.createVariable('SIC_RAW', 'd', ('fcst','time','n2d'))

    # attributes
    # Fcst_corr_to_write.units = 'change in hours'
    # Fcst_corr_to_write.units = '%'
    # Fcst_corr_to_write.description = 'SIP with bias correction with 4 leadtimes from JAN to DEC for each member.'
    # populate the data

    Fcst_calsic_to_write[:] = SIC_cal_arr
    Fcst_rawsic_to_write[:] = SIC_raw_arr
    # close ncfile
    ncfile.close()
    print('Saved file: ' + filename)
