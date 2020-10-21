#!/usr/bin/env python
# coding: utf-8
### Script for doing bias correction. Target year, which initialisation (1 to 4) and target month will be system/cmndline argument, so we can run this using sbatch. The 'historic' period will always be years 2003 to 2019.


import numpy as np
from taqm import taqm
from scipy.stats import linregress
from beinf import beinf
import matplotlib.pyplot as plt
import os
import sys
from netCDF4 import Dataset, num2date
data_path ='/work/ba1138/a270112/awicm3/FCST_CLIM/'
save_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults/'

# Which year are we targetting? And what month?
# python filenam.py 2015 1 6
targetyear = int(sys.argv[1])
leadtimeMonth = int(sys.argv[3])  #Leadtime, in month
# whichinit=1#
whichinit = int(sys.argv[2])  #Must be within 1 to 4
# histYrs=np.arange(2003,2011)
histYrs=np.arange(2003,targetyear)  # Now let's include all years until target within hist
obsTyr=targetyear

## The target of the leadtime should change according to which initialisation.
strtm = (1, 4, 7, 10)  # which is the starting month for each initialisation
obsTmnth = leadtimeMonth+strtm[whichinit-1]-1
if (obsTmnth>12):
    obsTyr=obsTyr+1
    obsTmnth=obsTmnth-12

print('Targetyear = '+str(targetyear)+',initialisation = '+str(whichinit)+' which means from '+ str(strtm[whichinit-1]) +',leadtime '+str(leadtimeMonth)+' so target month is '+str(obsTmnth)+' of year '+str(obsTyr))  # Testing

Grdlen = 126858  # For now we are pre-setting these. Dims[3] = 126858
EMlen = 10  # Dims[0]  #30 Ensemble Members
fcst_target = np.empty((EMlen,12,Grdlen))
## First let's save targetFcst
yr=str(targetyear-2000)       # current year
# file_anom0 = Dataset(data_path+'F'+str(yr)+'_MEM_ens_mon_mean_corr.nc')
for EM in np.arange(EMlen):
    file_anom0 = Dataset(data_path+'F'+str(yr).zfill(2)+str(whichinit)+'_ens_mon_mean/SIC_mon_'+str(EM+1).zfill(2)+'.nc')
    temp=file_anom0.variables['a_ice'][:]
    fcst_target[EM,:,:]=temp
    file_anom0.close()
    del temp

#plt.plot(fcst_anom[1,1,1,])
#plt.show()


## Observation for the date already exists?
file_osisaf = Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(obsTyr)+'.nc')
truobs=file_osisaf.variables['obs'][obsTmnth-1,:]
file_osisaf.close()

rawFcstCRPSS=np.empty(Grdlen)
calFcstCRPSS=np.empty(Grdlen)

obsSIP=np.empty(Grdlen)
rawSIP=np.empty(Grdlen)
calSIP=np.empty(Grdlen)

## Let's save Historical Forecast, for only month 2 of forecast 1 of each year in histYRs for now
histFcst=np.empty((len(histYrs),EMlen,Grdlen))

for c,year in  enumerate(histYrs):
    #print(year)
    yr=str(year-2000)       # current year
    lyr=str(year-1-2000)    # last year
    for EM in np.arange(EMlen):
        file_anom0 = Dataset(data_path+'F'+str(yr).zfill(2)+str(whichinit)+'_ens_mon_mean/SIC_mon_'+str(EM+1).zfill(2)+'.nc')
        temp=file_anom0.variables['a_ice'][:]
        histFcst[c,EM,:]=temp[leadtimeMonth-1,:]
        file_anom0.close()
        del temp



    # file_anom = Dataset(data_path+'F'+str(yr)+'_MEM_ens_mon_mean_corr.nc')
    # fcst_anom=file_anom.variables['SIC_FCST_CORR'][:]
    # histFcst[c,:,:]=fcst_anom[:,1,leadtimeMonth-1,:]
    # file_anom.close()



## Similarly, let's save Historical observation, for only Feb of each year in histYRs for now
histYrsforObs=histYrs
if (obsTyr>targetyear):
    histYrsforObs=np.arange(2003,targetyear)  # Now let's include all years until target within hist
histObs=np.empty((len(histYrsforObs),Grdlen))

for c,year in  enumerate(histYrsforObs):
    #print(year)
    file_osisaf = Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(year)+'.nc')
    obs=file_osisaf.variables['obs'][obsTmnth-1,:]
    histObs[c,:]=obs
    file_osisaf.close()


# Now the bias correction steps are run on each grid point. Let's start by trying just one, then later run a loop?


print("Input done, now calibrating")

g=2000
for g in np.arange(Grdlen):
    X=histFcst[:,:,g]
    Y=histObs[:,g]
    X_t=fcst_target[:,leadtimeMonth-1,g]
    taqminst = taqm()
    tau_t=histYrs
    t=targetyear

    # if(np.std(Y)<0.025):  #Maybe this should not be run the same way
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

    # X_t=fcst_target[:,1,leadtimeMonth-1,g] # Do we need this?
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

    Y_t = truobs[g]

    if(Y_t >= 0):  # To avoid the mask
        obsSIP[g] = np.int(Y_t >= 0.15)
    rawSIP[g] = sip_x_t
    calSIP[g] = sip_x_t_cal
    cdf_obs = np.zeros(len(x))
    cdf_obs[Y_t*np.ones(len(x)) <= x] = 1.0

    # CRPS for the raw forecast
    crps_x_t = np.trapz((cdf_x_t - cdf_obs)**2., x)
    # print (crps_x_t)
    # >>> 0.0277481871254

    # CRPS for the calibrated forecast
    crps_x_t_cal = np.trapz((cdf_x_t_cal - cdf_obs)**2., x)
    # print (crps_x_t_cal)
    rawFcstCRPSS[g] = crps_x_t
    calFcstCRPSS[g] = crps_x_t_cal
    del taqminst


print("Calibratting done! Now saving result file")


### the part where I try to save the calibrated forecasts:

filename = save_path+'Forecast_Calibration_BigHist_Yr'+str(targetyear)+'_'+str(whichinit).zfill(2)+'Mn_'+str(leadtimeMonth).zfill(2)+'.nc'
ncfile = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
# create dimensions
ncfile.createDimension('n2d', 126858)
# ncfile.createDimension('fcst',4)
# ncfile.createDimension('time',12)
# ncfile.createDimension('ens',30)
# define variables
Fcst_corr_to_write = ncfile.createVariable('SIP_FCST_CORR', 'd', ('n2d'))
Fcst_raw_to_write = ncfile.createVariable('SIP_FCST_RAW', 'd', ('n2d'))
Fcst_obs_to_write = ncfile.createVariable('SIP_FCST_OBS', 'd', ('n2d'))


# attributes
# Fcst_corr_to_write.units = 'change in hours'
# Fcst_corr_to_write.units = '%'
# Fcst_corr_to_write.description = 'SIP with bias correction with 4 leadtimes from JAN to DEC for each member.'
# populate the data
Fcst_corr_to_write[:] = calSIP
Fcst_raw_to_write[:] = rawSIP
Fcst_obs_to_write[:] = obsSIP


# close ncfile
ncfile.close()
print('Saved file: ' + filename)
