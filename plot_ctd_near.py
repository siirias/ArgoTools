# -*- coding: utf-8 -*-
"""
Created on Thu Dec 21 15:33:07 2017

@author: siirias
"""

import sys
sys.path.insert(0,'D:\\svnfmi_merimallit\\qa\\nemo')
import datetime as dt
import calendar
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import ModelQATools as qa
import math

#runfile('D:/ArgoData/plot_full_data.py', wdir='D:/ArgoData')
draw_images=True

#km/day range 
day_in_km=3
how_many_shown=10

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

def split_csv_profiles(pressure, other_vars, invalid_val=-10000000000.000):
    profiles=0
    depths=0
    max_depths=0
    parameter_num=len(other_vars)+1 #as pressure is one val
    for i in range(1,len(pressure)):
        depths+=1
        if(pressure[i]>pressure[i-1]): #Hypattiin seuraavaan profiiliin
            profiles+=1
            if(max_depths<depths):
                max_depths=depths
            depths=0
    result=np.ones((profiles+1,max_depths+10,parameter_num))*invalid_val
    #Fill the data:
    depth=0
    profile=0
    for i in range(1,len(pressure)):
        if(pressure[i]>pressure[i-1]): #Hypattiin seuraavaan profiiliin
            profile+=1
            depth=0
        if(result[profile][depth][0]!=invalid_val):
            if(abs(result[profile][depth][0]-pressure[i])>0.1):
                print "voi kräpylä", result[profile][depth][0]-pressure[i]
        result[profile][depth][0]=pressure[i]
        for j in range(len(other_vars)):
            result[profile][depth][j+1]=other_vars[j][i]
            
        depth+=1

    result=np.ma.masked_where(result==invalid_val,result)
    return result
        


"""
nc Variable	Variable	Units	Description

metavar1	Cruise		
metavar2	Station		
metavar3	Type		
longitude	Longitude	degrees_east	
latitude	Latitude	degrees_north	
metavar4	Bot. Depth	m	
metavar5	Secchi Depth	m	
date_time	Decimal Gregorian Days of the station	days since 2013-01-01 00:00:00 UTC	Relative Gregorian Days with decimal part
			
var1	PRES	db	
var2	TEMP	deg C	
var3	PSAL	psu	
var4	DOXY	ml/l	
var5	PHOS	umol/l	
var6	TPHS	umol/l	
var7	SLCA	umol/l	
var8	NTRA	umol/l	
var9	NTRI	umol/l	
var10	AMON	umol/l	
var11	NTOT	umol/l	
var12	H2SX	umol/l	
var13	PHPH		
var14	ALKY	meq/l	
var15	CPHL	ug/l	
var16	Year (station date)	
"""
#Gotlands deep
#lon_min=17;lat_min=56;lon_max=22;lat_max=59;
lon_min=16.5;lat_min=55.5;lon_max=22.5;lat_max=59.5;
target_lat=57.3; target_lon=20; target_rad=1.8*6 #rad in km

filetype='csv' # 'nc' tai 'csv'

if filetype=='nc':
    invalid_val=-10000000000.000
    file_n='d:/ArgoData/Siiriaetal2017/2013-2016_GotlDeep_data_from_helcom.nc'
    fmk=netcdf.netcdf_file(file_n,'r')
    press=-1.0*fmk.variables['var1'][:]
    longitude=fmk.variables['longitude'][:]
    latitude=fmk.variables['latitude'][:]
    start_epoch=dt.datetime(2013,1,1,0,0)
    start_secs=(start_epoch-dt.datetime.utcfromtimestamp(0.0)).total_seconds() #this should be amount of seconds to add to actual timestamps
    times_s=fmk.variables['date_time'][:]
    times=[]
    for i in range(len(times_s)):
        times.append(dt.datetime.utcfromtimestamp(times_s[i]*24.0*60.0*60.0+start_secs))
else:
    if filetype=='csv':
        invalid_val=-10000000000.000
        file_n='d:/ArgoData/Siiriaetal2017/Uudet_CTDt_16102017.csv'
        fmk=pd.read_csv(file_n)
        press=-1.0*fmk[u'PRES [db]'][:]
        longitude=fmk[u'Longitude [degrees_east]'][:]
        latitude=fmk[ u'Latitude [degrees_north]'][:]
        temperature=fmk[u'TEMP [deg C]'][:]
        salinity=fmk[u'PSAL [psu]'][:]

        times_s=fmk[u'yyyy-mm-ddThh:mm'][:]
        times=[]
        for i in range(len(times_s)):
            times.append(dt.datetime.strptime(times_s[i],'%Y-%m-%dT%H:%M'))
        ts=[]
        for i in range(len(times)):
            ts.append(calendar.timegm(times[i].timetuple()))
        ctd_data=split_csv_profiles(press,[longitude,latitude,temperature,salinity,ts])
        
        #tehdän se perus aikajana
        timesx=ctd_data[:,0,5]
        times=[]
        for i in range(len(timesx)):
            times.append(dt.datetime.utcfromtimestamp(timesx[i]))  
#        d=np.ma.masked_where(d==invalid_val,d)
        #create mask, for values close to main point:
        distance_mask=[False]*len(times_s)
        for i in range(len(times_s)):
            if(target_rad>distance((target_lat,target_lon),(latitude[i],longitude[i]))):
                distance_mask[i]=True
        #Ja perus paikkajanatkin
        pressure=ctd_data[:,:,0]
        latitude=ctd_data[:,0,2]
        longitude=ctd_data[:,0,1]
        temperature=ctd_data[:,:,3]
        salinity=ctd_data[:,:,4]
    else:
        print "wrong filetype,",filetype," aborting!"
        filetype='exitnow'


if(filetype!='exitnow'):
    #
    #
    #create mask, for values close to main point:
    distance_mask=[False]*len(times)
    print target_rad, target_lat,target_lon
    for i in range(len(times)):
        if(target_rad>distance((target_lat,target_lon),(latitude[i],longitude[i]))):
            distance_mask[i]=True


    #Let's plot times and measurements.
    if draw_images==True:
        plt.figure()
#        plt.plot(times,latitude,'*')
        plt.plot(longitude,latitude,'*')
#clf();plt.hist2d(latitude,longitude,bins=20);plt.colorbar()    
    
    
print "FINISHED!"