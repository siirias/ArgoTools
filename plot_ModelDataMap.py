# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:11:08 2016

@author: siirias
"""
import sys
sys.path.insert(0,'D:\\svnfmi_merimallit\\qa\\nemo')

import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import ModelQATools as qa
import ModelPltTools
from scipy.io import netcdf
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from netCDF4 import Dataset

#Full baltic
lon_min=9;lat_min=53.5;lon_max=30.3;lat_max=66; #26 GoB
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;

plot_bathymetry=False
plot_salt=True
plot_legends=True
plot_routes=True
DataDirectory="F:\\ArgoData\\"
fig=plt.figure(figsize=(13,13))
#fig=plt.figure(figsize=(10,10))


plt.clf()
bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = 'i',fix_aspect=False)

bmap.drawcoastlines(linewidth=0.2)
bmap.fillcontinents()
bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])
                

#TOPOGRAPHY EXPERIMENT
if plot_bathymetry:
    topodata = Dataset('iowtopo2_rev03.nc')
    
    topoin = topodata.variables['Z_WATER'][:].T
    lons = topodata.variables['XT_I'][:]
    lats = topodata.variables['YT_J'][:]
    x=np.tile(lons,(lats.shape[0],1))
    y=np.tile(lats,(lons.shape[0],1)).T
    bmap.pcolor(x,y,topoin,cmap='bone_r',vmin=0,vmax=200)
#    cb=plt.colorbar()
#    cb.ax.invert_yaxis()
#Model data experiment
if plot_salt:
    salt_min=0
    salt_max=12
#    dcolormap='viridis_r'
    dcolormap='jet'
    levs=[0,1,2,3,4,5,6,7,8,9,10,11,12]

    topodata = Dataset(DataDirectory+'fmi_hirlam_forecastv4_sd_20151112_12_D4.nc')
    
    topoin4 = topodata.variables['salt'][0,0,:,:]
    lons = topodata.variables['lon'][:]
    lats = topodata.variables['lat'][:]
    lons=lons+(lons[1]-lons[0])
    lats=lats+(lats[1]-lats[0])
    x=np.tile(lons,(lats.shape[0],1)).T
    y=np.tile(lats,(lons.shape[0],1))
    bmap.pcolor(x,y,topoin4,cmap=dcolormap,vmin=salt_min, vmax=salt_max)
 
    cb=plt.colorbar()
    cb.set_ticks(levs)

    CS=bmap.contour(x,y,topoin4,levels=levs)
    plt.clabel(CS, inline=1, fontsize=10,fmt='%1.1f')


    topodata = Dataset(DataDirectory+'fmi_hirlam_forecastv4_sd_20151112_12_D2.nc')
    
    topoin2 = topodata.variables['salt'][0,0,:,:]
    lons = topodata.variables['lon'][:]
    lats = topodata.variables['lat'][:]
    lons=lons+(lons[1]-lons[0])
    lats=lats+(lats[1]-lats[0])
    x=np.tile(lons,(lats.shape[0],1)).T
    y=np.tile(lats,(lons.shape[0],1))
    bmap.pcolor(x,y,topoin2,cmap=dcolormap,vmin=salt_min, vmax=salt_max)
    CS=bmap.contour(x,y,topoin2,levels=levs)
    plt.clabel(CS, inline=1, fontsize=10,fmt='%1.1f')



    topodata = Dataset(DataDirectory+'fmi_hirlam_forecastv4_sd_20151112_12_D1.nc')
    
    topoin1 = topodata.variables['salt'][0,0,:,:]
    lons = topodata.variables['lon'][:]
    lats = topodata.variables['lat'][:]
    lons=lons+(lons[1]-lons[0])
    lats=lats+(lats[1]-lats[0])
    x=np.tile(lons,(lats.shape[0],1)).T
    y=np.tile(lats,(lons.shape[0],1))
    bmap.pcolor(x,y,topoin1,cmap=dcolormap,vmin=salt_min, vmax=salt_max)

    CS=bmap.contour(x,y,topoin1,levels=levs)
    plt.clabel(CS, inline=1, fontsize=10,fmt='%1.1f')


#    cb=plt.colorbar()
#    cb.ax.invert_yaxis()

plt.savefig('BS_SaltMap.png' ,facecolor='w',dpi=400)
plt.savefig('BS_SaltMap.pdf' ,facecolor='w',dpi=400)
plt.savefig('BS_SaltMap.jpg' ,facecolor='w',dpi=400)
