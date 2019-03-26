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
import ModelQATools as qa
import ModelPltTools
import argohelper as ah

data_path="C:\\Data\\BarentsinMeri\\HaetutDatat\\"
save_path=data_path
alue='barents' #'barents' or 'bob'
#data_file_names=['6903695_20190104102220073.nc']  #barents
#data_file_names=['6902026_20190104103025617.nc']  #bob

shift_with_time=True
shift_step=1.
shift_type='index'  #'date','index'
single_color=True
start=mp.dates.datetime.datetime(1000,5,5)
end=mp.dates.datetime.datetime(3030,5,5)



#Plot route
fig=plt.figure(figsize=(10,10))
plt.clf()

if(alue=='barents'):
    data_file_names=['6903695_20190104102220073.nc']  #barents
    bmap = Basemap(width=300*1000,height=300*1000,projection='stere',lat_0=78.,lon_0=25., resolution = 'i',fix_aspect=False)
    
    bmap.etopo()
    bmap.drawcoastlines()
    bmap.fillcontinents()
    bmap.drawparallels(np.arange(70.,98.,2.),labels=[1,1,1,1],linewidth=0.5)
    bmap.drawmeridians(np.arange(-180.0,180.0,5.),labels=[1,1,1,1],linewidth=0.5)
else:
    #Bothnian Sea
    data_file_names=['6902026_20190104103025617.nc']  #bob
    lon_min=20;lat_min=63;lon_max=25;lat_max=66;
    
    bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
    resolution = 'i',fix_aspect=False)
    
    bmap.drawcoastlines()
    bmap.fillcontinents()
    bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
    bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])
    data_path="C:\\Data\\Pape1\\"




a=qa.PointData(data_path+data_file_names[0],1,start,end,"argonc");
lon_dat=a.obs['ape']['lon'][~a.obs['ape']['lon'].mask]
lat_dat=a.obs['ape']['lat'][~a.obs['ape']['lat'].mask]
date_axis=a.obs['ape']['date']
col='r'
x,y=bmap(np.array(lon_dat),np.array(lat_dat))
bmap.plot(x,y,'-',color=col,linewidth=2, alpha=0.4)
bmap.plot(x,y,'.',color=col,linewidth=2, alpha=0.8)
bmap.plot(x[-1],y[-1],'x',color='k',markersize=8,alpha=1.0)
bmap.plot(x[0],y[0],'o',color='k',markersize=8,alpha=1.0)
plt.savefig(save_path+'{}_route_new.png'.format(alue))

print("some statistics")
ah.file_names_converted=[data_path+data_file_names[0]]
ah.give_statistics()


profile_num=a.obs['ape']['tem'].shape[0]
#plot Temperature  profiles
plot_mark='-o'
plot_size=2
plot_linewidth=0.5
fig_size=(14,8)
the_colormap=mp.cm.get_cmap('Reds')

plt.figure(figsize=fig_size)
plt.clf()
order_no=0.0
shifting=-shift_step

for i in range(profile_num):
    rate=order_no/profile_num
    if(single_color):
        rate=1.0
    order_no+=1.
    used_color=the_colormap(rate)
    if(shift_with_time):
        if(shift_type=='index'):
            shifting+=shift_step
        if(shift_type=='date'):
            shifting=date_axis[i]
    plt.plot(a.obs['ape']['tem'][i,:]+shifting,-1.*a.obs['ape']['depth'][i,:],plot_mark,color=used_color,markersize=plot_size,linewidth=plot_linewidth)
plt.xlabel('Temperature')
plt.ylabel('Depth (m)')
plt.title('Temperature')
plt.savefig(save_path+'{}_T_profile.png'.format(alue))

#plot Salinity  profiles
shifting=-shift_step
the_colormap=mp.cm.get_cmap('Blues')

plt.figure(figsize=fig_size)
plt.clf()
order_no=0.0
for i in range(profile_num):
    rate=order_no/profile_num
    if(single_color):
        rate=1.0
    order_no+=1.
    used_color=the_colormap(rate)
    if(shift_with_time):
        if(shift_type=='index'):
            shifting+=shift_step
        if(shift_type=='date'):
            shifting=date_axis[i]

    plt.plot(a.obs['ape']['sal'][i,:]+shifting,-1.*a.obs['ape']['depth'][i,:],plot_mark,color=used_color,markersize=plot_size,linewidth=plot_linewidth)
plt.xlabel('Salinity')
plt.ylabel('Depth (m)')
plt.title('Salinity')
plt.savefig(save_path+'barents_S_profile.png')
