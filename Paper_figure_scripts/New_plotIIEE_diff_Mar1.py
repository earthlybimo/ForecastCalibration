#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import os
import sys
from netCDF4 import Dataset, num2date

# from mpl_toolkits.axes_grid1.inset_locator import inset_axes
data_path = '/work/ba1138/a270138/BiasCorrOutput/TAQMResults/MeanConc/LongStyle/'
Fig_path = '/work/ba1138/a270138/BiasCorrOutput/Figs/'
# CLIM_BIAS=readbin('CLIM_BIAS_04-10_4x12x126858f_method2.bin',[4,12,126858])
## Problem. where is readbin, and where is this bin file?


# In[2]:



mesh_area_file=Dataset('/mnt/lustre02/work/ab0995/a270112/data_fesom2/griddes.nc')
mesh_area=mesh_area_file.variables['cell_area'][:]
mesh_lat=mesh_area_file.variables['lat'][:]
mesh_lon=mesh_area_file.variables['lon'][:]
mesh_area_file.close()

IIEE_osisaf_cim_NH=np.zeros((2019-2011+1,12))
IIEE_osisaf_cim_SH=np.zeros((2019-2011+1,12))
IIEE_fcst_NH=np.zeros((2018-2011+1,4,12))
IIEE_fcst_SH=np.zeros((2018-2011+1,4,12))
IIEE_fcst_diff_NH=np.zeros((2018-2011+1,4,12))
IIEE_fcst_diff_SH=np.zeros((2018-2011+1,4,12))


# In[3]:


## Let's calculate IIEE for Climatology

for yr in np.arange(2011,2019+1):
    file_osisaf_clim=Dataset('/work/ab0995/a270112/data_fesom2/sic/OSISAF_MON_CLIM_'+str(yr-9)+'-'+str(yr-1)+'.nc')
    osisaf_climALL=file_osisaf_clim.variables['obs'][:]
    file_osisaf_clim.close()
    truobsfile=('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(yr)+'.nc')
    if not os.path.isfile(truobsfile):
        truobs=np.NaN
    else:
        file_osisaf = Dataset(truobsfile)
        obsAll=file_osisaf.variables['obs'][:]
        file_osisaf.close()
    for mm in np.arange(1,12+1):
        tmp=osisaf_climALL[(mm-1),:]
        obs=obsAll[(mm-1),:]
        nod_mismatch = ((mesh_lat>40.) & (obs > 0.15) & (obs <= 1.0) & (tmp < 0.15)) |                      ((mesh_lat>40.) & (obs>=0.0) & (obs < 0.15) & (tmp > 0.15))

        IIEE_osisaf_cim_NH[yr-2011, mm-1] =np.sum(mesh_area[nod_mismatch])/1.e12

        nod_mismatch = ((mesh_lat<-30.) & (obs > 0.15) & (obs <= 1.0) & (tmp < 0.15)) |                      ((mesh_lat<-30.) & (obs>=0.0) & (obs < 0.15) & (tmp > 0.15))

        IIEE_osisaf_cim_SH[yr-2011, mm-1] =np.sum(mesh_area[nod_mismatch])/1.e12


# In[19]:


## Let's calculate IIEE for forecasts

strtm=(1,4,7,10)
for yr in np.arange(2011,2018+1):
    File1=data_path+'F'+str((yr-2000))+'_ens_mon_mean_corr.nc'
    print(File1)
    temp=Dataset(File1)
    meanSIC=temp.variables['SIC_FCST_CORR'][:]
    temp.close()
    for init in np.arange(1,4+1):
        for mm in np.arange(1,12+1):
            obsTyr=yr
            obsTmnth = mm+strtm[init-1]-1
            if (obsTmnth>12):
                obsTyr=obsTyr+1
                obsTmnth=obsTmnth-12
            print('Targetyear = '+str(yr)+',initialisation = '+str(init)+' which means from '+ str(strtm[init-1]) +',leadtime '+str(mm)+' so target month is '+str(obsTmnth)+' of year '+str(obsTyr))  # Testing

            truobsfile=('/work/ab0995/a270112/data_fesom2/sic/OSISAF_monthly_'+str(obsTyr)+'.nc')
            if not os.path.isfile(truobsfile):
                obs=np.NaN
            else:
                file_osisaf = Dataset(truobsfile)
                obs=file_osisaf.variables['obs'][obsTmnth-1,:]
                file_osisaf.close()
            tmp=meanSIC[(init-1),(mm-1),:]
            nod_mismatch = ((mesh_lat>40.) & (obs > 0.15) & (obs <= 1.0) & (tmp < 0.15)) |                          ((mesh_lat>40.) & (obs>=0.0) & (obs < 0.15) & (tmp > 0.15))

            IIEE_fcst_NH[yr-2011, init-1, mm-1] =np.sum(mesh_area[nod_mismatch])/1.e12
            IIEE_fcst_diff_NH[yr-2011, init-1, mm-1]=IIEE_fcst_NH[yr-2011, init-1, mm-1]-IIEE_osisaf_cim_NH[obsTyr-2011, obsTmnth-1]
            
            
            
            nod_mismatch = ((mesh_lat<-30.) & (obs > 0.15) & (obs <= 1.0) & (tmp < 0.15)) |                          ((mesh_lat<-30.) & (obs>=0.0) & (obs < 0.15) & (tmp > 0.15))

            IIEE_fcst_SH[yr-2011, init-1, mm-1] =np.sum(mesh_area[nod_mismatch])/1.e12
            IIEE_fcst_diff_SH[yr-2011, init-1, mm-1]=IIEE_fcst_SH[yr-2011, init-1, mm-1]-IIEE_osisaf_cim_SH[obsTyr-2011, obsTmnth-1]


# In[32]:


## IIEE line-plot for OSI SAF climatology

plt.close()
plt.figure(figsize=(9,4))
plt.plot(np.mean(IIEE_osisaf_cim_NH,0),marker='D',color='k',label = 'Arctic ')
plt.plot(np.mean(IIEE_osisaf_cim_SH,0),linestyle="dashed",marker='v',color='k',label = 'Antarctic ')
plt.xlim([-1, 12])
plt.xticks(np.arange(12),('J', 'F', 'M', 'A', 'M','J', 'J', 'A', 'S', 'O','N','D'),fontsize=12)
plt.ylabel(r'IIEE $(10^6 km^2)$',fontsize=12)
plt.legend(frameon=False,loc=9)
plt.title("OSI SAF climatology")
plt.savefig((Fig_path+'/IIEE_seasonalCycle_OSISAF_climatologyonly.png'),dpi=300)

plt.show()


# In[16]:


IIEE_mean_displacedNH=np.ma.zeros((4,21))
IIEE_mean_displacedNH[0,0:12]=np.squeeze(np.mean(IIEE_fcst_diff_NH[:,0,:],0))
IIEE_mean_displacedNH[1,3:(12+3)]=np.squeeze(np.mean(IIEE_fcst_diff_NH[:,1,:],0))
IIEE_mean_displacedNH[2,6:(12+6)]=np.squeeze(np.mean(IIEE_fcst_diff_NH[:,1,:],0))
IIEE_mean_displacedNH[3,9:(12+9)]=np.squeeze(np.mean(IIEE_fcst_diff_NH[:,1,:],0))
IIEE_mean_displacedNH[IIEE_mean_displacedNH==0]=np.ma.masked

IIEE_mean_displacedSH=np.ma.zeros((4,21))
IIEE_mean_displacedSH[0,0:12]=np.squeeze(np.mean(IIEE_fcst_diff_SH[:,0,:],0))
IIEE_mean_displacedSH[1,3:(12+3)]=np.squeeze(np.mean(IIEE_fcst_diff_SH[:,1,:],0))
IIEE_mean_displacedSH[2,6:(12+6)]=np.squeeze(np.mean(IIEE_fcst_diff_SH[:,1,:],0))
IIEE_mean_displacedSH[3,9:(12+9)]=np.squeeze(np.mean(IIEE_fcst_diff_SH[:,1,:],0))
IIEE_mean_displacedSH[IIEE_mean_displacedSH==0]=np.ma.masked


# plt.pcolormesh(IIEE_mean_displacedNH)
# plt.show()
# plt.close()


# In[34]:



plt.close()
plt.figure(figsize=(9.5,6))
plt.subplot(211)
ax = plt.gca()
plt.pcolormesh(np.flipud(IIEE_mean_displacedNH),cmap=plt.cm.RdBu_r);cbar=plt.colorbar(extend='both')
cbar.set_label('IIEE ('+r'$x10^6 {km}^2$)')
ax.set_yticks(np.arange(0.5,4.0,1))
ax.set_yticklabels(["Oct","Jul","Apr","Jan"])
ax.set_xticks(np.arange(0.5,21.0,1))
ax.set_xticklabels(["J","F","M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A","S"])
plt.xlim([0, 21])
plt.title('Arctic')
# plt.show()

plt.subplot(212)
ax = plt.gca()
plt.pcolormesh(np.flipud(IIEE_mean_displacedSH),cmap=plt.cm.RdBu_r);cbar=plt.colorbar(extend='both')
cbar.set_label('IIEE ('+r'$x10^6 {km}^2$)')
ax.set_yticks(np.arange(0.5,4.0,1))
ax.set_yticklabels(["Oct","Jul","Apr","Jan"])
ax.set_xticks(np.arange(0.5,21.0,1))
ax.set_xticklabels(["J","F","M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A","S"])
plt.xlim([0, 21])
plt.title('Antarctic')
plt.tight_layout()
plt.savefig((Fig_path+'/IIEE_calibratedMeanSICminusClimatology.png'),dpi=300)
plt.show()


# In[ ]:




