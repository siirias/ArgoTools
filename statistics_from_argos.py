# -*- coding: utf-8 -*-
"""
Created on Mon May 14 11:27:31 2018

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
import pandas as pd
import datetime as dt
import calendar
import time
import pandas as pd


last_time=time.time()
invalid_val=-10000000000.000
#        file_n='./Siiriaetal2017/uudet_CTDt_16102017.csv'
#file_n='./Siiriaetal2017/0316544c.csv'
target_lat=57.32; target_lon=20.05; target_rad=30.0 #rad in km

file_n=ah.ctd_data_file
fmk=pd.read_csv(file_n)
press=-1.0*fmk[u'PRES [db]'][:]
longitude=fmk[u'Longitude [degrees_east]'][:]
latitude=fmk[ u'Latitude [degrees_north]'][:]
temperature=fmk[u'TEMP [deg C]'][:]
salinity=fmk[u'PSAL [psu]'][:]
oxygen=fmk[u'DOXY [ml/l]'][:]
print "CTD's read from file, time T+{}".format(time.time()-last_time)
last_time=time.time()
#Muutetaan saliniteetit g/kg:ksi
for point in range(salinity.shape[0]):
    salinity[point]=ah.abs_suolaisuus(salinity[point],longitude[point],latitude[point])[0]
for i in range(len(oxygen)):
    if(type(oxygen[i])==str and oxygen[i][0]=='<'):
       oxygen[i]=None 
    else:
       oxygen[i]=float(oxygen[i])
print "CTD's converted, time T+{}".format(time.time()-last_time)
last_time=time.time()

times_s=fmk[u'yyyy-mm-ddThh:mm'][:]
times=[]
for i in range(len(times_s)):
    times.append(dt.datetime.strptime(times_s[i],'%Y-%m-%dT%H:%M'))
ts=[]
for i in range(len(times)):
    ts.append(calendar.timegm(times[i].timetuple()))
ctd_data=ah.split_csv_profiles(press,[longitude,latitude,temperature,salinity,ts,oxygen])
   #tehdän se perus aikajana
timesx=ctd_data[:,0,5]
times=[]
for i in range(len(timesx)):
    times.append(dt.datetime.utcfromtimestamp(timesx[i]))  
#        d=np.ma.masked_where(d==invalid_val,d)

print "CTD's split, time-line set, time T+{}".format(time.time()-last_time)
last_time=time.time()


#create mask, for values close to main point:
#distance_mask=[False]*len(times_s)
#for i in range(len(times_s)):
#    if(target_rad>ah.distance((target_lat,target_lon),(latitude[i],longitude[i]))):
#        distance_mask[i]=True

#Ja perus paikkajanatkin
pressure=ctd_data[:,:,0]
latitude=ctd_data[:,0,2]
longitude=ctd_data[:,0,1]
temperature=ctd_data[:,:,3]
salinity=ctd_data[:,:,4]
oxygen=ctd_data[:,:,6]
print "CTD's masked, time T+{}".format(time.time()-last_time)
last_time=time.time()

#create mask, for values close to main point:
distance_mask=[False]*len(times)
for i in range(len(times)):
    if(target_rad>ah.distance((target_lat,target_lon),(latitude[i],longitude[i]))):
        distance_mask[i]=True


files=ah.file_names_converted
timesteps=[]
last_accepted=True
a_time_all=[]
a_distance_mask_all=[]
for file_n in files:
    fmk=netcdf.netcdf_file(file_n,'r')
    temp=fmk.variables['TEMP_ADJUSTED'][:].copy()
    salt=fmk.variables['PSAL_ADJUSTED'][:].copy()
    press=fmk.variables['PRES_ADJUSTED'][:].copy()
    oxyg=fmk.variables['DOXY'][:].copy()
    scat=fmk.variables['SCATTERING'][:].copy()
    apetime = np.array([datetime.datetime.fromordinal(x) for x in fmk.variables['TIME'][:]])
    a_lats=fmk.variables['LATITUDE'][:].copy()
    a_lons=fmk.variables['LONGITUDE'][:].copy()
    
    unique_times=[]
    for i in apetime:
        if len(unique_times)==0 or i-unique_times[-1] > datetime.timedelta(0) or last_accepted==False:
            if(len(unique_times)>0):
                timesteps.append((i-unique_times[-1]).days)
            unique_times.append(i)
            last_accepted=True
        else:
            last_accepted=False

    #calculate oxygen on 20 m
    a_ox=[]
    a_time=[]
    for i in range(len(apetime)):
        dat=ah.get_closest(press[i,:],oxyg[i,:],20.0)
        if(not np.isnan(dat[1]) and target_rad>ah.distance((target_lat,target_lon),(a_lats[i],a_lons[i]))):
            a_ox.append(dat[1])
            a_time.append(apetime[i])
    plt.plot(a_time,a_ox,'.-')



    print "Unique measurements in float: {}".format(len(unique_times))
    a_time_all+=a_time
print "Argo's read, time T+{}".format(time.time()-last_time)
last_time=time.time()


#check also the oxygen levels from ctd's
f_depths=[];
f_oxygen=[];
f_time=[];
for i in range(len(pressure[:,0])):
    stats=ah.get_closest(np.array(pressure[i,:]),np.array(oxygen[i,:]),-20.)
    if(stats[1]>-1000.0):
        f_oxygen.append(stats[1])
        f_depths.append(stats[0])    
        f_time.append(times[i])
        
#figure()
plt.plot(f_time,f_oxygen,'.')

print "CTD's plotted T+{}".format(time.time()-last_time)
last_time=time.time()

#make the histogram style image:
a_time_all=pd.DatetimeIndex(a_time_all)
times=pd.DatetimeIndex(times)
plt.figure(figsize=[15,6])
dates_per_month=pd.date_range(min(a_time_all)-datetime.timedelta(days=30),max(a_time_all+datetime.timedelta(days=30)),freq='M')
argo_profiles=np.zeros(len(dates_per_month))
ctd_profiles=np.zeros(len(dates_per_month))
for i in range(len(dates_per_month)-1):
     argo_profiles[i]=sum(np.array(a_time_all>=dates_per_month[i]) * np.array(a_time_all<=dates_per_month[i+1]))
     ctd_profiles[i]=sum(np.array(times>=dates_per_month[i]) * np.array(times<=dates_per_month[i+1]) * np.array(distance_mask))

plt.gca().grid(color='lightgray',linewidth=0.5,axis='y',zorder=0)
p1=plt.bar(dates_per_month, argo_profiles,width=26,zorder=10)
p2=plt.bar(dates_per_month, ctd_profiles,width=16,zorder=20)
plt.ylabel('Profiles')
plt.yticks(np.arange(0, 15, 1))
plt.legend((p1[0], p2[0]), ('Argo', 'Ship'))
locs,labels = plt.xticks()
plt.setp(labels,rotation=45)
plt.gca().xaxis.set_major_locator(mp.dates.MonthLocator(range(1,12,2)))
plt.gca().xaxis.set_major_formatter(mp.dates.DateFormatter('%m-%y'))
