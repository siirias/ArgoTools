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
import math 
import argohelper as ah

def plot_full_distance(new_fig=None, save_file=None, ref_point=None, no_plot=False, ref_is_no=None, vmin=None, vmax=None):


    def distance(origin, destination): 
        lat1, lon1 = origin 
        lat2, lon2 = destination 
        radius = 6371 # km 
        dlat = math.radians(lat2-lat1) 
        dlon = math.radians(lon2-lon1) 
        a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
            * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
        d = radius * c 
        return d

#Actual plot_full_distance starts here        

    files=["IM_6902014_20130814_20140821.nc", \
           "IM_6902019_20140821_20150805.nc", \
           "IM_6902020_20150805_20160331_active.nc"]
    
    files=ah.file_names_cleaned
    distance_aggregate=np.array(());
    distance_r_aggregate=np.array(());
    lats_aggregate=np.array(());
    lons_aggregate=np.array(());
    time_aggregate=np.array(());
    if(new_fig!=None):
        try:
            fig=plt.figure(figsize=new_fig)   
        except:
            fig=plt.figure()   
    if(not no_plot):
        plt.clf()

    current_set=-1
    tmin=None
    tmax=None
    prev_lat=None
    prev_lon=None
    if(ref_is_no is not None):
        for file_n in files:
            fmk=netcdf.netcdf_file(file_n,'r')
            lats=fmk.variables['LATITUDE'][::2].copy()
            lons=fmk.variables['LONGITUDE'][::2].copy()
            if(np.size(lats)>ref_is_no):
                ref_point=(lats[ref_is_no],lons[ref_is_no])
                break
            else:
                ref_is_no=ref_is_no-np.size(lats)
    for file_n in files:
        current_set+=1
        fmk=netcdf.netcdf_file(file_n,'r')
        
        reftime = datetime.datetime.strptime(fmk.variables['REFERENCE_DATE_TIME'][:].tostring(), '%Y%m%d%H%M%S') #"YYYYMMDDHHMISS"
        jultime = fmk.variables['JULD'][:].tolist()
        apetime = np.array([mp.dates.date2num(reftime + datetime.timedelta(days=x)) for x in jultime])
        
        lats=fmk.variables['LATITUDE'][::2].copy()
        lons=fmk.variables['LONGITUDE'][::2].copy()
        time=apetime[::2]
        if(tmin is None or time[0]<tmin):
            tmin=time[0]
        if(tmax is None or time[-1]>tmax):
            tmax=time[-1]
        distances=lats[:].copy() #Just to get right size
        #first:
        if prev_lat is not None:
            distances[0]=distance((lats[0],lons[0]),(prev_lat,prev_lon))
        else:
            distances[0]=0.0 #start place
        for i in range(1,distances.shape[0]):
            distances[i]=distance((lats[i],lons[i]),(lats[i-1],lons[i-1]))
        
        distance_from_ref=distances[:]*0.0
        if(ref_point is not None):
            for i in range(0,distances.shape[0]):
                distance_from_ref[i]=distance((lats[i],lons[i]),ref_point)
            
        prev_lat=lats[-1]
        prev_lon=lons[-1]
        distance_aggregate=np.concatenate([distance_aggregate,distances])
        distance_r_aggregate=np.concatenate([distance_r_aggregate,distance_from_ref])
        time_aggregate=np.concatenate([time_aggregate,time])
        lats_aggregate=np.concatenate([lats_aggregate,lats])
        lons_aggregate=np.concatenate([lons_aggregate,lons])
        #diff_data=dat[:][:].copy()
        #ACTUAL PLOT 
        if(not no_plot):
            if(ref_point is not None):
                plt.plot( time,distance_from_ref,'k')
            else:
                plt.plot( time,distances,'k')
            

    if(not no_plot):
        ax=plt.gca()
        ax.set_facecolor('#ffffff')
        plt.xlim((tmin,tmax))
        plt.ylim((vmin,vmax))
        plt.ylabel('Distance (km)')
        plt.xlabel('Time [month-year]')
        if(ref_point is not None):
            plt.title('Distance From (%.2f, %.2f)' % (ref_point[0],ref_point[1]))
            plt.plot([plt.xlim()[0],plt.xlim()[1]],[30,30],'g--',alpha=0.4)
#            plt.plot([plt.xlim()[0],plt.xlim()[1]],[60,60],'y--',alpha=0.4)
#            plt.plot([plt.xlim()[0],plt.xlim()[1]],[90,90],'r--',alpha=0.4)
        else:
            plt.title('Distance between profiles')
        plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
        plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%d-%m-%y'))
        plt.grid(axis='y')
        plt.yticks(range(0,int(plt.ylim()[1]),15))
        #plt.set_cmap('gist_stern')
        locs,labels = plt.xticks()
        plt.ylim([0,plt.ylim()[1]]) #Ensure that the lower limit is zero.
        plt.setp(labels,rotation=45)
    if(save_file!=None):
        plt.savefig(save_file,dpi=300)
    if(ref_point is not None):
        return (time_aggregate,distance_r_aggregate,lats_aggregate,lons_aggregate)
    else:
        return (time_aggregate,distance_aggregate,lats_aggregate,lons_aggregate)

#plt.clf()
#plot_full_distance(ref_point=(57.315,20.04))
    