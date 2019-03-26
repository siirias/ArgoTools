# -*- coding: utf-8 -*-
"""
Created on Fri Nov 16 17:02:26 2018

@author: siirias
"""
import sys
sys.path.insert(0,'C:\\svnfmi_merimallit\\qa\\nemo')

import datetime
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import ModelQATools as qa
import ModelPltTools
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from netCDF4 import Dataset
import pandas as pd

#Barent's
lon_min=20;lat_min=70;lon_max=70;lat_max=80;
#Full baltic
#lon_min=12;lat_min=53;lon_max=30;lat_max=66;
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;
#Bothnian Sea
#lon_min=17;lat_min=60;lon_max=22;lat_max=63;

plot_bathymetry=True
plot_legends=True
plot_routes=True

fig=plt.figure(figsize=(10,10))

#fig=plt.figure(figsize=(10,13))
#fig=plt.figure(figsize=(10,10))


plt.clf()
#bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
#resolution = 'l',fix_aspect=False)
bmap = Basemap(width=300*1000,height=300*1000,projection='stere',lat_0=78.,lon_0=25., resolution = 'f',fix_aspect=False)

bmap.etopo()
bmap.drawcoastlines()
bmap.fillcontinents()
#bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
#bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])
bmap.drawparallels(np.arange(70.,98.,2.),labels=[1,1,1,1],linewidth=0.5)
bmap.drawmeridians(np.arange(-180.0,180.0,5.),labels=[1,1,1,1],linewidth=0.5)

filest_to_plot=[\
                "6902020_prof.nc", "6902021_prof.nc", "6902022_prof.nc", \
                "6902023_prof.nc", "6902024_prof.nc"
                ]
                
colors=["#ff0000","#00ff00","#0000ff",\
        "#000000","#d00000" \
        ]

labels=[\
        "6902020", "6902021", "6902022", \
        "6902023", "6902024" \
        ]

               
start=mp.dates.datetime.datetime(1000,5,5)
end=mp.dates.datetime.datetime(3030,5,5)

"""
#TOPOGRAPHY EXPERIMENT
if plot_bathymetry:
    topodata = Dataset('iowtopo2_rev03.nc')
    
    topoin = topodata.variables['Z_WATER'][:]
    lons = topodata.variables['XT_I'][:]
    lats = topodata.variables['YT_J'][:]
    x=np.tile(lons,(lats.shape[0],1))
    y=np.tile(lats,(lons.shape[0],1)).T
    bmap.pcolor(x,y,-1*topoin,cmap='bone_r',vmin=0,vmax=300)
    cb=plt.colorbar()
    cb.ax.invert_yaxis()
"""

"""
if plot_routes:
    for f,col,lab in zip(filest_to_plot,colors,labels):
        a=qa.PointData(f,1,start,end,"argonc");
        x,y=bmap(a.obs['ape']['lon'][:],a.obs['ape']['lat'][:])
    #    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
        bmap.plot(x,y,color=col,linewidth=2, alpha=0.8)
        bmap.plot(x[-1],y[-1],'o',color=col,markersize=10,alpha=1.0,label=lab)
    #    print lab, mp.dates.num2date(a.obs['ape']['date'][0]).date() \
    #             , mp.dates.num2date(a.obs['ape']['date'][-1]).date()
"""
data_path="C:\\Data\\BarentsinMeri\\"
file_name='gps_log.txt'
float_route=pd.read_csv(data_path+file_name,names=['number','time','lat','lon','unknown'])
col='r'
x,y=bmap(np.array(float_route['lon']),np.array(float_route['lat']))
bmap.plot(x,y,'.-',color=col,linewidth=2, alpha=0.8)
bmap.plot(x[-1],y[-1],'x',color='k',markersize=8,alpha=1.0)
bmap.plot(x[0],y[0],'o',color='k',markersize=8,alpha=1.0)

#if plot_legends:
#    plt.legend(bbox_to_anchor=(1.0,0.5),numpoints=1)

#plt.savefig('ArgoRoutes2015-2016.png' ,facecolor='w',dpi=300)
#plt.savefig('ArgoRoutes2015-2016.eps' ,facecolor='w',dpi=300)
