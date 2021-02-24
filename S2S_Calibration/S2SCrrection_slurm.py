#!/usr/bin/env python
# coding: utf-8
### Script for doing bias correction. Target year, which initialisation (1 to 4) and target month will be system/cmndline argument, so we can run this using sbatch. The 'historic' period will always be years 2003 to 2019.

import os
import sys
from netCDF4 import Dataset, num2date
import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
from datetime import datetime, timedelta, date
import glob
sys.path.append("/mnt/lustre01/pf/a/a270138/Scripts/ForecastCalibration/TAQM_Method")
from taqm import taqm
from beinf import beinf
trust_sharp_fcst = False

#
HEM="nh"
datasource = "/mnt/lustre01/work/ab0995/a270099/S2S/"  # In Lorenzos' side
sat_folder=datasource+"satellite_interp"

model_name="ECMWF"
model_path=datasource+"forecasts/"+model_name+"/"
save_path = '/work/ba1138/a270138/BiasCorrOutput/S2S_Results/'+model_name+"/"
if not(os.path.isdir(save_path)):
    os.makedirs(save_path)
#

# # Which year are we targetting? And what month?
# python filenam.py 2015 1
# targetyear = 2005;initMonth=5
## Actually now we can make this into a sys arg so we can run a slurm loop
targetyear = int(sys.argv[1])
initMonth=int(sys.argv[2])

histYrs=np.arange(1999,targetyear)  # Now let's include all years until target within hist

fname1=(model_path+model_name+"_ref_cont_sic_"+str(targetyear)+"-"+str(initMonth).zfill(2)+"*.nc")
flist=glob.glob(fname1);flist.sort();file1=flist[0]
concNC=Dataset(file1)
ci_cont=concNC.variables['ci'][:]
concNC.close()

file2=file1.replace("cont", "pert")
if not(os.path.isfile(file2)):
    sys.exit((" Pert file: "+ file2+ " does not exist!"))
concNC=Dataset(file2)
ci_pert=concNC.variables['ci'][:]
concNC.close()
Dims=ci_pert.shape
rawFcst=np.ma.empty((Dims[0],(Dims[1]+1),Dims[2],Dims[3]))  #  Leadtime, EMno, i, j
rawFcst[:,0:Dims[1],:,:]=ci_pert.copy()
rawFcst[:,Dims[1],:,:]=ci_cont.copy()

calFcst=np.ma.empty((Dims[0],Dims[2],Dims[3]))
rawSIP=calFcst.ma.copy()
obsSIP=calFcst.ma.copy()

filename = save_path+'TAQM_calibrated_'+os.path.basename(file2)
if (os.path.isfile(filename)):
    sys.exit((" File: "+ filename+ " already exists!"))

#
histFcst=np.ma.empty((len(histYrs),Dims[0],(Dims[1]),Dims[2],Dims[3]))  # Fcstyear, Leadtime, EMno, i, j
## Or If we wish to include the control runs as well
histFcst=np.ma.empty((len(histYrs),Dims[0],(Dims[1]+1),Dims[2],Dims[3]))  # Fcstyear, Leadtime, EMno, i, j

for c,year in  enumerate(histYrs):
    file3=file2.replace(str(targetyear), str(year))
    # os.path.isfile(file2)  DO THIS CHECK LATER!
    concNC=Dataset(file3)
    ci_pert=concNC.variables['ci'][:]
    concNC.close()
    # histFcst[c,:,:,:]=ci_pert
    #
    ## Or If we wish to include the control runs as well
    file4=file3.replace("pert", "cont")
    # os.path.isfile(file1)  DO THIS CHECK LATER!
    concNC=Dataset(file4)
    ci_cont=concNC.variables['ci'][:]
    concNC.close()
    histFcst[c,:,0:Dims[1],:,:]=ci_pert
    histFcst[c,:,Dims[1],:,:]=ci_cont;


temp=file1.find(".nc")
datetemp=(file1[(temp-2):temp])
initdate=date(targetyear,initMonth,int(datetemp))
targetObs=np.ma.empty((Dims[0],Dims[2],Dims[3]))  # Fcstyear, Leadtime
for f in np.arange((Dims[0])): #Along time, starting from initial day to final day
    trgtdate=initdate+timedelta(int(f))
    satfile=(sat_folder+"/ice_conc_"+HEM+"_ease-125_reproc_"+trgtdate.strftime("%Y%m%d")+"1200.nc4")
    if not(os.path.isfile(satfile)):
        continue
    concNC=Dataset(satfile)
    sat_ci=concNC.variables['ice_conc'][:]
    concNC.close()
    targetObs[f,:,:]=sat_ci

targetObs=targetObs/100
targetObs[targetObs<0]=np.nan

histObs=np.empty((len(histYrs),Dims[0],Dims[2],Dims[3]))  # Fcstyear, Leadtime
for c,year in  enumerate(histYrs):
    initdate=date(year,initMonth,int(datetemp))
    for f in np.arange((Dims[0])): #Along time, starting from initial day to final day
        trgtdate=initdate+timedelta(int(f))
        satfile=(sat_folder+"/ice_conc_"+HEM+"_ease-125_reproc_"+trgtdate.strftime("%Y%m%d")+"1200.nc4")
        if not(os.path.isfile(satfile)):
            continue
        concNC=Dataset(satfile)
        sat_ci=concNC.variables['ice_conc'][:]
        concNC.close()
        histObs[c,f,:,:]=sat_ci


histObs=histObs/100
histObs[histObs<0]=np.nan
print("Input done, now calibrating")

## Grid loops should be here:
i=1;j=1;lt=1
for lt in np.arange(2):#(Dims[0])):
    print("Leadtime: "+str(lt))
    for i in np.arange((Dims[2])):
        for j in np.arange((Dims[3])):

            X=histFcst[:,lt,:,i,j]
            Y=histObs[:,lt,i,j]
            X_t=rawFcst[lt,:,i,j]
            x2=np.ma.empty_like(X_t)
            x2.mask=X_t.mask.copy
            x2.data[X_t.data>=0.15]=1
            if all(X_t.mask== True):
                rawSIP[lt,i,j].mask=True
                calFcst[lt,i,j].mask=True
                continue
            if all(np.isnan(Y)):
                rawSIP[lt,i,j].mask=True
                calFcst[lt,i,j].mask=True
                continue

            taqminst = taqm()
            tau_t=histYrs
            t=targetyear
            rawSIP[lt,i,j]=np.ma.mean(x2)
            pval_x = linregress(tau_t,X.mean(axis=1))[3]  #check the p-value for MH trend over tau_t
            if pval_x<0.05:
                # if significant, then adjust MH for the trend to create TAMH
                X_ta = taqminst.trend_adjust_1p(X,tau_t,t)
                X_ta=X_ta.compressed()
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
            #
            X_ta_params, Y_ta_params, X_t_params = taqminst.fit_params(X_ta,Y_ta,X_t)

            # Now calibrate the forecast ensemble using the calibrate() method:

            X_t_cal_params, X_t_cal = taqminst.calibrate(X_ta_params, Y_ta_params, X_t_params,X_ta, Y_ta, X_t,trust_sharp_fcst)

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

            # fig = plt.figure()
            # ax1 = fig.add_subplot(1,2,1)
            # ax1.set_ylim((0.0,1.01))
            # ax1.plot(x, cdf_y_ta, 'g-',label='TAMH', lw=1.5)
            # ax1.plot(x, cdf_x_ta, 'y-',label='TAOH', lw=1.5)
            # ax1.legend(loc='lower right')
            #
            # # fig = plt.figure()
            # ax2 = fig.add_subplot(1,2,2)
            # ax2.set_ylim((0.0,1.01))
            # ax2.plot(x, cdf_x_t, 'b-',label='Raw Fcst', lw=1.5)
            # ax2.plot(x, cdf_x_t_cal, 'r-',label='TAQM Fcst', lw=1.5)
            # ax2.legend(loc='lower right')

            # plt.show()

            Y_t = targetObs[lt,i,j]
            #
            # if(Y_t >= 0):  # To avoid the mask
            #     obsSIP[g] = np.int(Y_t >= 0.15)  # Condition, so gives 1 or 0

            cdf_obs = np.zeros(len(x))
            cdf_obs[Y_t*np.ones(len(x)) <= x] = 1.0
            # CRPS for the raw forecast
            crps_x_t = np.trapz((cdf_x_t - cdf_obs)**2., x)
            # print (crps_x_t)
            # >>> 0.0277481871254
            # CRPS for the calibrated forecast
            crps_x_t_cal = np.trapz((cdf_x_t_cal - cdf_obs)**2., x)
            # print (crps_x_t_cal)
            del taqminst

            calFcst[lt,i,j]=sip_x_t_cal
            # rawSIP[lt,i,j]=sip_x_t

### the part where I try to save the calibrated forecasts:

filename = save_path+'TAQM_calibrated_'+os.path.basename(file2)
ncfile = Dataset(filename, 'w', format='NETCDF4_CLASSIC')
# create dimensions
ncfile.createDimension('time',Dims[0])
ncfile.createDimension('i', Dims[2])
ncfile.createDimension('j',Dims[3])
# ncfile.createDimension('ens',30)
# define variables
Fcst_corr_to_write = ncfile.createVariable('SIP_FCST_CORR', 'd', ('time','i','j'))
Fcst_raw_to_write = ncfile.createVariable('SIP_FCST_RAW', 'd',('time','i','j'))
Fcst_obs_to_write = ncfile.createVariable('SIP_FCST_OBS', 'd', ('time','i','j'))


# attributes
# Fcst_corr_to_write.units = 'change in hours'
# Fcst_corr_to_write.units = '%'
# Fcst_corr_to_write.description = 'SIP with bias correction with 4 leadtimes from JAN to DEC for each member.'
# populate the data
Fcst_corr_to_write[:] = calFcst
Fcst_raw_to_write[:] = rawSIP
Fcst_obs_to_write[:] = targetObs


# close ncfile
ncfile.close()
print('Saved file: ' + filename)
