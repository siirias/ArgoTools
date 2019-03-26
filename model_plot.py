# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 12:17:58 2016

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

data=qa.GriddedData('fmi_hirlam_forecastv4_sd_20151216_12_D4.nc','hbm',varlist=['temp','salt'])

salt=data.get_var('salt')[0,0,:,:].copy().T

"""
#Oma yritys lataamiseen
ncf=netcdf.netcdf_file('fmi_hirlam_forecastv4_sd_20151216_12_D4.nc','r')
salt=ncf.variables['salt'][0][0][:][:].copy()
salt=salt
mask=salt<-9000
salt_m=salt[:].copy()
salt_m[mask]=np.nan

lat=ncf.variables['lat'][:].copy()
lat=np.tile(lat,(salt_m.shape[0],1))
lon=ncf.variables['lon'][:].copy()
lon=np.tile(lon,(salt_m.shape[1],1)).T
"""
lat=data.get_axis('lat')
#lat=np.tile(x,(1,salt.shape[0]))
lon=data.get_axis('lon')
#lon=np.tile(y,(1,salt.shape[1])).T

#BASEMAP TESTS
#Full baltic
lon_min=15;lat_min=53.5;lon_max=30;lat_max=66;




fig=plt.figure(figsize=(12.0, 16.0))
plt.clf()

bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = 'l',fix_aspect=False)
bmap.drawcoastlines()
bmap.fillcontinents()
bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])

"""
x=np.tile(lons,(lats.shape[0],1))
y=np.tile(lats,(lons.shape[0],1)).T
bmap.pcolor(x,y,-1*topoin,cmap='bone_r',vmin=0,vmax=300)
cb=plt.colorbar()
cb.ax.invert_yaxis()
"""
#plt.pcolor(lon,lat,np.ma.masked_invalid(salt_m),cmap='cool')
plt.pcolor(lon,lat,salt,cmap='cool')
plt.title("Salinity")
ax.set_axis_bgcolor('gray')

plt.colorbar()
plt.savefig('Salinity_map.png',dpi=300)
