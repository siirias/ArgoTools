# -*- coding: utf-8 -*-
"""
Created on Fri Apr 22 12:14:18 2016

@author: siirias
"""

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap



def plot_full_diff_data( value="temp",col_map='cool', time_highlight=None, \
                         set_no=None, high_no=None, vmin=None, vmax=None, \
                         new_fig=None, save_file=None, distances=None):
    def nearest(data,value):
        near=np.nan
        for i in range(data.size):
            if(not np.isnan(data[i])):
                if(np.isnan(near) or np.abs(data[i]-value)<np.abs(data[near]-value)):
                    near=i
        return near
    
    
    def prof_diff(a_dat, a_depth, b_dat, b_depth):
        diff=a_dat[:].copy()
        max_depth=0.0
        for d in a_depth:
            if(not np.isnan(d) and np.abs(d)>max_depth):
                max_depth=np.abs(d)
        max_depth_temp=0.0
        for d in b_depth:
            if(not np.isnan(d) and np.abs(d)>max_depth_temp):
                max_depth_temp=np.abs(d)
        max_depth=np.min([max_depth,max_depth_temp]) #max depth, which both have
                
        for d in range(a_dat.size):
            if(not np.isnan(a_depth[d])):
                diff[d]=b_dat[nearest(b_depth,a_depth[d])]-a_dat[d]
    
        length=nearest(a_depth,max_depth)
        diff[length:]=np.nan
        return [diff, a_depth]

#Actual plot_full_diff_data starts here        

    files=["IM_6902014_20130814_20140821.nc", \
           "IM_6902019_20140821_20150805.nc", \
           "IM_6902020_20150805_20160331_active.nc"]



    if(new_fig!=None):
        try:
            fig=plt.figure(figsize=new_fig)   
        except:
            fig=plt.figure()   
    plt.clf()
                

    if(time_highlight!=None):
        time_highlight=mp.dates.date2num(datetime.datetime.strptime(time_highlight,"%Y%m%d"))
    if(value=='temp'):
        kuva='temp_b'
    if(value=='salt'):
        kuva='salt_b'
    if(value=='temp_coarse'):
        kuva='temp_a'
    if(value=='salt_coarse'):
        kuva='salt_a'
    if(value=='oxygen'):
        kuva='oxyg_a'
    if(value=='scatter'):
        kuva='scat_a'
    current_set=-1
    for file_n in files:
        current_set+=1
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
            z_source=temp_a;y_source=press_a;tt=apetime_a;title_txt='Temperature [$^\circ$C]'
        if kuva=='temp_b':
            z_source=temp_b;y_source=press_b;tt=apetime_b;title_txt='Temperature [$^\circ$C]'
    
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
        pre=np.ma.masked_invalid(y_source.T)
        time=x_source
        
        if(set_no is not None and high_no is not None):
            if(set_no==current_set):
                time_highlight=time[0][high_no]


        diff_data=dat[:][:].copy()
        for i in range(dat.shape[1]-1):
            p1_dat=dat[:,i].copy()
            p1_dep=pre[:,i].copy()
            p2_dat=dat[:,i+1].copy()
            p2_dep=pre[:,i+1].copy()
            [d_dat,d_dep]=prof_diff(p1_dat,p1_dep,p2_dat,p2_dep)
            diff_data[:,i]=d_dat.copy()

        #diff_data=dat[:][:].copy()
        #ACTUAL PLOT 
        plt.pcolor( np.ma.masked_invalid(time), \
                    np.ma.masked_invalid(pre),  \
                    np.ma.masked_invalid(diff_data), \
                    cmap=col_map,vmin=vmin,vmax=vmax)
        
        ax=plt.gca()
        ax.set_axis_bgcolor('#101010')
        plt.ylabel('Pressure [dbar]')
        plt.xlabel('Time [month-year]')
        plt.title('Difference: ' + title_txt)
        plt.gca().invert_yaxis()
        plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
        plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%m-%y'))
        #plt.set_cmap('gist_stern')
        

    if time_highlight is not None:
        plt.plot([time_highlight,time_highlight],[plt.ylim()[0],plt.ylim()[1]] \
        ,color="#ff0000",linewidth=3)
    plt.colorbar(label=title_txt)
    locs,labels = plt.xticks()
    plt.setp(labels,rotation=45)
    
    if(distances is not None):
        target_range=0.2*(plt.ylim()[0]-plt.ylim()[1])
        multip=target_range/distances[1].max()
        plt.plot(distances[0],plt.ylim()[0]-(distances[1]*multip),'w.-')
        l1=plt.ylim()[0]-30*multip
        l2=plt.ylim()[0]-60*multip
        l3=plt.ylim()[0]-90*multip
        plt.plot([plt.xlim()[0],plt.xlim()[1]],[l1,l1],'g--',alpha=0.7)
        plt.plot([plt.xlim()[0],plt.xlim()[1]],[l2,l2],'y--',alpha=0.4)
        plt.plot([plt.xlim()[0],plt.xlim()[1]],[l3,l3],'r--',alpha=0.4)
        

    if(save_file!=None):
        plt.savefig(save_file,dpi=300)

    