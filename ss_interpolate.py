# -*- coding: utf-8 -*-
"""
Created on Mon May 30 15:00:28 2016

@author: siirias
"""

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf

col_map='cool'
vmin=0.0
vmax=20.0
file_n="IM_6902014_20130814_20140821.nc"

kuva='temp_b'

fmk=netcdf.netcdf_file(file_n,'r')
temp=fmk.variables['TEMP'][:].copy()
press=fmk.variables['PRES'][:].copy()

reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
jultime = fmk.variables['JULD'][:].tolist()
apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])


mask=press>9000
press_m=press[:].copy()
press_m[mask]=np.nan

temp_m=temp[:].copy()
temp_m[mask]=np.nan

apetime_a=apetime[::2].copy()
apetime_b=apetime[1::2].copy()
press_a=press_m[::2][:].copy()
press_b=press_m[1::2][:].copy()

temp_a=temp_m[::2][:].copy()
temp_b=temp_m[1::2][:].copy()

z_source=temp_a;y_source=press_a;tt=apetime_a;title_txt='Temperature [$^\circ$C]'
x_source=np.tile(tt,(z_source.shape[1],1))
x,y=np.mgrid[0:z_source.shape[0],0:z_source.shape[1]]
x=x.T
y=-1*y.T
dat=np.ma.masked_invalid(z_source.T) # z_source.T
pre=np.ma.masked_invalid(y_source.T)
time=x_source
plt.pcolor( np.ma.masked_invalid(time), \
                np.ma.masked_invalid(pre),  \
                np.ma.masked_invalid(dat), \
                cmap=col_map,vmin=vmin,vmax=vmax)
    
ax=plt.gca()
ax.set_axis_bgcolor('gray')
plt.ylabel('Pressure [dbar]')
plt.xlabel('Time [month-year]')
plt.gca().invert_yaxis()
plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%m-%y'))
plt.colorbar(label=title_txt)
locs,labels = plt.xticks()
plt.setp(labels,rotation=45)
