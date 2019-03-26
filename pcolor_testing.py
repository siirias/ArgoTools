# -*- coding: utf-8 -*-
"""
Created on Thu Apr 07 16:13:02 2016

@author: siirias
"""

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
files=["IM_6902014_20130814_20140821.nc", \
       "IM_6902019_20140821_20150805.nc", \
       "IM_6902020_20150805_20160331_active.nc"]
       
       
for kuva,col_map in zip(['temp_b', 'salt_b', 'oxyg_a','scat_a' ],\
                        ['hot',    'jet',    'viridis','BrBG_r']):
    plt.close()
    fig=plt.figure(figsize=(16.0, 10.0))
    for file_n in files:
        fmk=netcdf.netcdf_file(file_n,'r')
        temp=fmk.variables['TEMP'][:].copy()
        salt=fmk.variables['PSAL'][:].copy()
        press=fmk.variables['PRES'][:].copy()
        oxyg=fmk.variables['DOXY'][:].copy()
        scat=fmk.variables['SCATTERING'][:].copy()
        
        reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
        jultime = fmk.variables['JULD'][:].tolist()
        apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
        
        
        mask=press>9000
        press_m=press[:].copy()
        press_m[mask]=np.nan
        
        temp_m=temp[:].copy()
        temp_m[mask]=np.nan
        
        salt_m=salt[:].copy()
        salt_m[mask]=np.nan
        
        oxyg_m=oxyg[:].copy()
        oxyg_m[mask]=np.nan
        
        scat_m=scat[:].copy()
        scat_m[mask]=np.nan
        
        apetime_a=apetime[::2].copy()
        apetime_b=apetime[1::2].copy()
        press_a=press_m[::2][:].copy()
        press_b=press_m[1::2][:].copy()
    
        salt_a=salt_m[::2][:].copy()
        salt_b=salt_m[1::2][:].copy()
        temp_a=temp_m[::2][:].copy()
        temp_b=temp_m[1::2][:].copy()
        oxyg_a=oxyg_m[::2][:].copy()
        oxyg_b=oxyg_m[1::2][:].copy()
        scat_a=scat_m[::2][:].copy()
        scat_b=scat_m[1::2][:].copy()
        
        if kuva=='salt_a':
            z_source=salt_a;y_source=press_a;tt=apetime_a;title_txt='PSU'
        if kuva=='salt_b':
            z_source=salt_b;y_source=press_b;tt=apetime_b;title_txt='PSU'
    
        if kuva=='temp_a':
            z_source=temp_a;y_source=press_a;tt=apetime_a;title_txt='Temperature [C$^\circ$]'
        if kuva=='temp_b':
            z_source=temp_b;y_source=press_b;tt=apetime_b;title_txt='Temperature [C$^\circ$]'
    
        if kuva=='oxyg_a':
            z_source=oxyg_a;y_source=press_a;tt=apetime_a;title_txt='Oxygen [$\mu mol/kg$]'
        if kuva=='oxyg_b':
            z_source=oxyg_b;y_source=press_b;tt=apetime_b;title_txt='Oxygen [$\mu mol/kg$]'
        
        if kuva=='scat_a':
            z_source=scat_a;y_source=press_a;tt=apetime_a;title_txt='Scattering [$M^{-1} sr^{-1}$]'
        if kuva=='scat_b':
            z_source=scat_b;y_source=press_b;tt=apetime_b;title_txt='Scattering [$M^{-1} sr^{-1}$]'
        
        
        x_source=np.tile(tt,(z_source.shape[1],1))
        
        x,y=np.mgrid[0:z_source.shape[0],0:z_source.shape[1]]
        x=x.T
        y=-1*y.T
        dat=np.ma.masked_invalid(z_source.T) # z_source.T
        pre=-1*np.ma.masked_invalid(y_source.T)
        time=x_source
        #plt.pcolor(x,y,np.ma.masked_invalid(np.rot90(temp_m)))
    #    plt.figure()
    #    plt.pcolor(x,y,dat)
    #    plt.gca().set_axis_bgcolor('gray')
    #   plt.colorbar()
        
    #    plt.figure()

        #ACTUAL PLOT 
        plt.pcolor(time,pre,dat,cmap=col_map)
        
        
        ax=plt.gca()
        ax.set_axis_bgcolor('gray')
        plt.ylabel('Pressure [dbar]')
        plt.xlabel('Time [month-year]')
        plt.title(title_txt)
        plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
        plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%m-%y'))
        #plt.set_cmap('gist_stern')
        
    plt.colorbar()
    locs,labels = plt.xticks()
    plt.setp(labels,rotation=45)
    plt.savefig('%s.png' % (kuva),dpi=300)

