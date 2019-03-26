# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 18:11:00 2018

@author: siirias
"""

# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 17:02:26 2018

@author: siirias
"""
import sys
import os
import re

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from netCDF4 import Dataset
import pandas as pd
from datetime import datetime

data_path="C:\\Data\\BarentsinMeri\\"
save_path=data_path
file_name='gps_log.txt'
all_file_names=os.listdir(data_path)
data_file_names=[]
for i in all_file_names:
    if bool(re.search('\.csv$',i)):
        data_file_names.append(i)
        

#Plot route
fig=plt.figure(figsize=(10,10))
plt.clf()
bmap = Basemap(width=300*1000,height=300*1000,projection='stere',lat_0=78.,lon_0=25., resolution = 'i',fix_aspect=False)

bmap.etopo()
bmap.drawcoastlines()
bmap.fillcontinents()
bmap.drawparallels(np.arange(70.,98.,2.),labels=[1,1,1,1],linewidth=0.5)
bmap.drawmeridians(np.arange(-180.0,180.0,5.),labels=[1,1,1,1],linewidth=0.5)


data_path="C:\\Data\\BarentsinMeri\\"
file_name='gps_log.txt'
float_route=pd.read_csv(data_path+file_name,names=['number','time','lat','lon','unknown'])
col='r'
x,y=bmap(np.array(float_route['lon']),np.array(float_route['lat']))
bmap.plot(x,y,'.-',color=col,linewidth=2, alpha=0.8)
bmap.plot(x[-1],y[-1],'x',color='k',markersize=8,alpha=1.0)
bmap.plot(x[0],y[0],'o',color='k',markersize=8,alpha=1.0)
plt.savefig(save_path+'barents_route.png')

#GET THE DATA
all_pts_data=[]
all_profile_data=[]
for i in data_file_names:
    tmp_data=open(data_path+i,'r').readlines()
    pts_data=[]
    profile_data=[]
    for data_line in tmp_data:
        d=(data_line.rstrip()).split(',')
        if(d[0] ==  'CTD_PTS' or d[0] == 'CTD_CP'):   #either PTS or profile data. format is the same.
            d[1]=datetime.strptime(d[1],"%Y%m%dT%H%M%S")
            d[2]=-1.0*float(d[2])
            d[3]=float(d[3])
            d[4]=float(d[4])
            if(d[0] ==  'CTD_PTS'):  #Now separate two types of data.
                pts_data.append(d[1:])
            else:
                profile_data.append(d[1:])
    pts_data=pd.DataFrame(pts_data,columns=['time','P','T','S'])
    profile_data=pd.DataFrame(profile_data,columns=['time','P','T','S','samples'])
    if(len(pts_data)>1):
        all_pts_data.append(pts_data)
    if(len(profile_data)>1):
        all_profile_data.append(profile_data)



plot_mark='.'
plot_size=2
fig_size=(14,8)
#plot Temperature  PTS
plt.figure(figsize=fig_size)
plt.clf()
order_no=0.0
for i in all_pts_data:
    rate=order_no/len(all_pts_data)
    order_no+=1.
    plt.plot(i['T'],i['P'],plot_mark,color=(rate,rate*0.6,rate*0.6),markersize=plot_size)
plt.xlabel('Temperature')
plt.ylabel('Pressure')
plt.title('Temperature (PTS points)')
plt.savefig(save_path+'barents_T_PTS.png')


#Plot Salinity PTS
plt.figure(figsize=fig_size)
plt.clf()
order_no=0.0
for i in all_pts_data:
    rate=order_no/len(all_pts_data)
    order_no+=1.
    plt.plot(i['S'],i['P'],plot_mark,color=(rate*0.6,rate*0.6,rate),markersize=plot_size)
plt.xlabel('Salinity')
plt.ylabel('Pressure')
plt.title('Salinity (PTS points)')
plt.savefig(save_path+'barents_S_PTS.png')


plot_mark='-o'
plot_size=2
plot_linewidth=0.5

#plot Temperature  profiles

plt.figure(figsize=fig_size)
plt.clf()
order_no=0.0
for i in all_profile_data:
    rate=order_no/len(all_pts_data)
    order_no+=1.
    plt.plot(i['T'],i['P'],plot_mark,color=(rate,rate*0.6,rate*0.6),markersize=plot_size,linewidth=plot_linewidth)
plt.xlabel('Temperature')
plt.ylabel('Pressure')
plt.title('Temperature (Profiles)')
plt.savefig(save_path+'barents_T_profile.png')


#plot Salinity  profiles
plt.figure(figsize=fig_size)
plt.clf()
order_no=0.0
for i in all_profile_data:
    rate=order_no/len(all_pts_data)
    order_no+=1.
    plt.plot(i['S'],i['P'],plot_mark,color=(rate*0.6,rate*0.6,rate),markersize=plot_size,linewidth=plot_linewidth)
plt.xlabel('Salinity')
plt.ylabel('Pressure')
plt.title('Salinity (Profiles)')
plt.savefig(save_path+'barents_S_profile.png')

