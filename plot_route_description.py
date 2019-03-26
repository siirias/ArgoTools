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
from matplotlib import gridspec
#Full baltic
#lon_min=12;lat_min=53;lon_max=30;lat_max=66;
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;

#Gulf of Bothnia
lon_min=-5;lat_min=45;lon_max=31;lat_max=66;
#Gotlands deep
glon_min=18.0;glat_min=55.0;glon_max=21.0;glat_max=59.0;


lon_mean=0.5*(lon_min+lon_max)
lat_mean=0.5*(lat_min+lat_max)

plot_bathymetry=False

fig=plt.figure(figsize=(14,14))
#fig=plt.figure(figsize=(10,10))


plt.clf()
bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = 'i',fix_aspect=True)

#bmap = Basemap(projection='nsper',lon_0=lon_mean,lat_0=lat_mean,
#        satellite_height=500000.,resolution='h')

#bmap.drawcoastlines()
#bmap.fillcontinents()

#gludge to get rid of the lakes
for i,cp in enumerate(bmap.coastpolygons):
     if bmap.coastpolygontypes[i]<2:
         bmap.plot(cp[0],cp[1],'k-')
bmap.fillcontinents([0.9,0.9,0.9],lake_color=[0.85,0.85,0.85]) 


#bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
#bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])

#bmap.drawparallels(np.arange(50.,69,1.),linewidth=1,dashes=[3,10])
#bmap.drawmeridians(np.arange(12.,30,1.),linewidth=1,dashes=[3,10])

#TOPOGRAPHY EXPERIMENT
if plot_bathymetry:
    topodata = Dataset('iowtopo2_rev03.nc')
    
    topoin = topodata.variables['Z_WATER'][:]
    lons = topodata.variables['XT_I'][:]
    lats = topodata.variables['YT_J'][:]
#    x=np.tile(lons,(lats.shape[0],1))
#    y=np.tile(lats,(lons.shape[0],1)).T
    tmp_lon,tmp_lat =bmap(*np.meshgrid(lons,lats))
    print "LON", tmp_lon.shape
    print "LAT", tmp_lat.shape
#    bmap.pcolor(tmp_lon,tmp_lat,-1*topoin,cmap='bone_r',vmin=0,vmax=150)
    bmap.contourf(tmp_lon,tmp_lat,-1*topoin,20,cmap='bone_r',vmin=0,vmax=200)
    #cb=plt.colorbar()
    #cb.ax.invert_yaxis()

#PLOT HIGHLIGHT
x,y=bmap([18.0,18.0,21.0,21.0,18.0],[55.0,59.0,59.0,55.0,55.0])
x,y=bmap([glon_min,glon_min,glon_max,glon_max,glon_min],[glat_min,glat_max,glat_max,glat_min,glat_min])
#    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
bmap.plot(x,y,color='red',linewidth=2, alpha=0.8)
        
#PLOT INSET FIGURE:
#ax2=fig.add_subplot(121)
#plt.subplots_adjust(top=0.6,bottom=0.4)
gs=gridspec.GridSpec(3,3, width_ratios=[1,5,4], height_ratios=[1,5,1])
ax2=plt.subplot(gs[4])
bmap2=Basemap(llcrnrlon=glon_min,llcrnrlat=glat_min,urcrnrlon=glon_max,urcrnrlat=glat_max, \
resolution = 'i',fix_aspect=True)
for i,cp in enumerate(bmap2.coastpolygons):
     if bmap2.coastpolygontypes[i]<2:
         bmap2.plot(cp[0],cp[1],'k-')
bmap2.fillcontinents([0.9,0.9,0.9],lake_color=[0.85,0.85,0.85]) 


plt.savefig('GulfOfBothnia_far.png' ,facecolor='w',dpi=600)
plt.savefig('GulfOfBothnia_far.pdf' ,facecolor='w',dpi=600)
plt.savefig('GulfOfBothnia_far.svg' ,facecolor='w',dpi=600)
plt.savefig('GulfOfBothnia_far.eps' ,facecolor='w',dpi=600)



