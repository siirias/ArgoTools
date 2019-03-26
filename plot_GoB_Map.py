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
lon_min=17;lat_min=60;lon_max=30;lat_max=66; #26 GoB
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;

plot_bathymetry=True
plot_legends=True
plot_routes=True

fig=plt.figure(figsize=(13,13))
#fig=plt.figure(figsize=(10,10))


plt.clf()
bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = 'f',fix_aspect=False)

bmap.drawcoastlines(linewidth=0.2)
bmap.fillcontinents()
bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])
                

#TOPOGRAPHY EXPERIMENT
if plot_bathymetry:
    topodata = Dataset('iowtopo2_rev03.nc')
    
    topoin = topodata.variables['Z_WATER'][:]
    lons = topodata.variables['XT_I'][:]
    lats = topodata.variables['YT_J'][:]
    x=np.tile(lons,(lats.shape[0],1))
    y=np.tile(lats,(lons.shape[0],1)).T
    bmap.pcolor(x,y,-1*topoin,cmap='bone_r',vmin=0,vmax=200)
#    cb=plt.colorbar()
#    cb.ax.invert_yaxis()

plt.savefig('GoB_Map.png' ,facecolor='w',dpi=400)
