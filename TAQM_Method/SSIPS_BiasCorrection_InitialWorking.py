#!/usr/bin/env python
# coding: utf-8

import numpy as np
from taqm import taqm
from scipy.stats import linregress
from beinf import beinf
import matplotlib.pyplot as plt
import os
from netCDF4 import Dataset, num2date
data_path='/work/ba1138/a270112/awicm3/FCST_CLIM/'


## Which year are we targetting? And what month?
histYrs=np.arange(2011,2015)

targetyear=2013
targetmonth=2  #For 4 fcsts of each year, it will be different month. For now, we are only doing Fcst 1


## First let's save targetFcst
yr=str(targetyear-2000)       # current year
file_anom0 = Dataset(data_path+'F'+str(yr).zfill(2)+'_MEM_ens_mon_mean_corr.nc')
fcst_target=file_anom0.variables['SIC_FCST_CORR'][:]
file_anom0.close()
#plt.plot(fcst_anom[1,1,1,])
#plt.show()

Dims=fcst_target.shape
Grdlen=Dims[3]  #roughly 126858
EMlen=Dims[0]  #30 Ensemble Members

## Observation for the date already exists?
file_osisaf = Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(targetyear).zfill(2)+'.nc')
truobs=file_osisaf.variables['obs'][targetmonth-1,:]
file_osisaf.close()

rawFcstCRPSS=np.empty(Grdlen)
calFcstCRPSS=np.empty(Grdlen)

obsSIP=np.empty(Grdlen)
rawSIP=np.empty(Grdlen)
calSIP=np.empty(Grdlen)

## Let's save Historical Forecast, for only month 2 of forecast 1 of each year in histYRs for now
histFcst=np.empty((len(histYrs),EMlen,Grdlen))
c=0
for year in histYrs:
    #print(year)
    yr=str(year-2000)       # current year
    lyr=str(year-1-2000)    # last year
    file_anom = Dataset(data_path+'F'+str(yr).zfill(2)+'_MEM_ens_mon_mean_corr.nc')
    fcst_anom=file_anom.variables['SIC_FCST_CORR'][:]
    file_anom.close()
    histFcst[c,:,:]=fcst_anom[:,1,targetmonth-1,:]
    c=c+1


## Similarly, let's save Historical observation, for only Feb of each year in histYRs for now
histObs=np.empty((len(histYrs),Grdlen))
c=0
for year in histYrs:
    #print(year)
    yr=str(year-2000)       # current year
    lyr=str(year-1-2000)    # last year
    file_osisaf = Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(year).zfill(2)+'.nc')
    obs=file_osisaf.variables['obs'][targetmonth-1,:]
    file_osisaf.close()
    histObs[c,:]=obs
    c=c+1


# Now the bias correction steps are run on each grid point. Let's start by trying just one, then later run a loop?

g=2000
for g in np.arange(Grdlen):
    X=histFcst[:,:,g]
    Y=histObs[:,g]
    X_t=fcst_target[:,1,targetmonth-1,g]
    taqminst = taqm()
    tau_t=histYrs
    t=targetyear

    # if(np.std(Y)<0.025):  #This should probably not be run the same way
        # print("Should skip this grid point. Somehow.")
        # continue
        # print("low stDev for Y")
    # Get TAMH from MH
    pval_x = linregress(tau_t,X.mean(axis=1))[3]  #check the p-value for MH trend over tau_t
    if pval_x<0.05:
        # if significant, then adjust MH for the trend to create TAMH
        X_ta = taqminst.trend_adjust_1p(X,tau_t,t)
    else:
        # else, set TAMH equal to MH (i.e. don't perform the trend adjustment)
        X_ta = np.copy(X)


    pval_y = linregress(tau_t,Y)[3]     #check p-value for OH trend over tau_t
    if pval_y<0.05:
        # if significant, then adjust OH for the trend to create TAOH
        Y_ta = taqminst.trend_adjust_1p(Y,tau_t,t)
    else:
        # else, set TAOH equal to OH (i.e. don't perform the trend adjustment)
        Y_ta = np.copy(Y)

    X_t=fcst_target[:,1,targetmonth-1,g]
    X_ta_params, Y_ta_params, X_t_params = taqminst.fit_params(X_ta,Y_ta,X_t)

    trust_sharp_fcst = True

    # Now calibrate the forecast ensemble using the calibrate() method:

    X_t_cal_params, X_t_cal = taqminst.calibrate(X_ta_params, Y_ta_params, X_t_params,X_ta, Y_ta, X_t,trust_sharp_fcst)

    ## Why are some points problematic!!!!

    # Next, we’re going to compute the SIP quantity for the raw and calibrated forecast, plot all cumulative distributions, and calculate the continuous rank probability score (CRPS) for the raw and calibrated forecast.

    # First, evaluate the cdf for each of these using the cdf_eval() method in the beinf class. This method handles instances when a and b aren’t known (and given the value np.inf), in which case the cdf over (0,1) is computed using the ecdf() method. When a and b are known (as is the case in this example), cdf_eval() evaluates the cdf using the cdf() method. We can also use the cdf_eval() method to compute SIP.


    x = np.linspace(0, 1, 1000)
    x_c = 0.15

    # Evaluate cdf for the TAMH distribution at x
    cdf_x_ta = beinf.cdf_eval(x, X_ta_params, X_ta)

    # Evaluate cdf for the TAOH distribution at x
    cdf_y_ta = beinf.cdf_eval(x, Y_ta_params, Y_ta)

    # Evaluate cdf for the forecast distribution at x and calculate SIP
    cdf_x_t = beinf.cdf_eval(x, X_t_params, X_t)
    sip_x_t = 1.0 - beinf.cdf_eval(x_c, X_t_params, X_t)


    p_x_t = X_t_params[2] # raw forecast
    p_x_ta = X_ta_params[2] # TAMH climatology
    p_y_ta = Y_ta_params[2] # TAOH climatology

    # Evaluate cdf for the calibrated forecast distribution at x and calculate SIP
    if trust_sharp_fcst==True and p_x_t==1.0:
        # go with the original forecast data/distribution when any of the p parameters are one
        # for the three distributions used in calibration
        cdf_x_t_cal = beinf.cdf_eval(x, X_t_params, X_t)
        sip_x_t_cal = 1.0 - beinf.cdf_eval(x_c, X_t_params, X_t)
    else:
        if p_x_t==1.0 or p_x_ta==1.0 or p_y_ta==1.0:
            # go with the TAOH data/distribution when any of the p parameters are
            # one for the three distributions used in calibration
            cdf_x_t_cal = beinf.cdf_eval(x, Y_ta_params, Y_ta)
            sip_x_t_cal = 1.0 - beinf.cdf_eval(x_c, Y_ta_params, Y_ta)
        else:
            # go with the calibrated forecast data/distribution
            cdf_x_t_cal = beinf.cdf_eval(x, X_t_cal_params, X_t_cal)
            sip_x_t_cal = 1.0 - beinf.cdf_eval(x_c, X_t_cal_params, X_t_cal)


    # In[30]: Plot CDF

    # fig = plt.figure()
    # ax1 = fig.add_subplot(1,2,1)
    # ax1.set_ylim((0.0,1.01))
    # ax1.plot(x, cdf_y_ta, 'g-',label='TAMH', lw=1.5)
    # ax1.plot(x, cdf_x_ta, 'y-',label='TAOH', lw=1.5)
    # ax1.legend(loc='lower right')

    # # fig = plt.figure()
    # ax2 = fig.add_subplot(1,2,2)
    # ax2.set_ylim((0.0,1.01))
    # ax2.plot(x, cdf_x_t, 'b-',label='Raw Fcst', lw=1.5)
    # ax2.plot(x, cdf_x_t_cal, 'r-',label='TAQM Fcst', lw=1.5)
    # ax2.legend(loc='lower right')

    # plt.show()

    Y_t=truobs[g]

    if(Y_t>=0):  #To avoid the mask
        obsSIP[g]=np.int(Y_t>=0.15)
    rawSIP[g]=sip_x_t
    calSIP[g]=sip_x_t_cal
    cdf_obs = np.zeros(len(x))
    cdf_obs[Y_t*np.ones(len(x))<=x] = 1.0

    # CRPS for the raw forecast
    crps_x_t = np.trapz((cdf_x_t - cdf_obs)**2.,x)
    # print (crps_x_t)
    # >>> 0.0277481871254

    # CRPS for the calibrated forecast
    crps_x_t_cal = np.trapz((cdf_x_t_cal - cdf_obs)**2.,x)
    # print (crps_x_t_cal)
    rawFcstCRPSS[g]=crps_x_t
    calFcstCRPSS[g]=crps_x_t_cal
    del taqminst


print("Calibratting done! Now saving a file")


### the part where I try to save the calibrated forecasts:

file = 'Forecast_Calibration_Yr',targetyear,'_01Mn_',targetmonth,'.nc'
ncfile = Dataset(file,'w',format='NETCDF4_CLASSIC')
#create dimensions
ncfile.createDimension('n2d',126858)
# ncfile.createDimension('fcst',4)
# ncfile.createDimension('time',12)
# ncfile.createDimension('ens',30)
#define variables
Fcst_corr_to_write = ncfile.createVariable('SIP_FCST_CORR','d',('n2d'))
Fcst_raw_to_write = ncfile.createVariable('SIP_FCST_RAW','d',('n2d'))
Fcst_obs_to_write = ncfile.createVariable('SIP_FCST_OBS','d',('n2d'))


# attributes
# Fcst_corr_to_write.units = 'change in hours'
# Fcst_corr_to_write.units = '%'
# Fcst_corr_to_write.description = 'SIP with bias correction with 4 leadtimes from JAN to DEC for each member.'
#populate the data
Fcst_corr_to_write[:] =calSIP
Fcst_raw_to_write[:] =rawSIP
Fcst_obs_to_write[:] =obsSIP


#close ncfile
ncfile.close()

print("DONE!!!")
