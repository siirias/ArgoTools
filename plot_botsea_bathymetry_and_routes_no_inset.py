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
#lon_min=12;lat_min=53;lon_max=30;lat_max=66;
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;
def plot_SR5(bmap,size_multiplier=1.0):
    x,y=bmap([19.3555],[61.0520])
    bmap.plot(x,y,'X',color="#ff2020",markersize=8*size_multiplier,alpha=1.0,zorder=20,label="SR5")    

#Gulf of Bothnia
#lon_min=-5;lat_min=53;lon_max=31;lat_max=66;
#lon_min=14.00;lat_min=60.00;lon_max=22.3;lat_max=63.5;
lon_min=17.00;lat_min=60.50;lon_max=22.0;lat_max=63.0; #as reviewer suggester
lon_mean=0.5*(lon_min+lon_max)
lat_mean=0.5*(lat_min+lat_max)
filest_to_plot=["petranc//6901901_prof.nc","petranc//6902013_prof.nc","petranc//6902017_prof.nc", \
                "petranc//6902018_prof.nc","petranc//6902021_prof.nc", \
                "petranc//6902022_prof.nc"]
                
colors=["#fa8072","#ef9aef","#fff550","#4682b4","#00ffff","#00ff00"]
labels=["6901901","6902017","6902013", \
        "6902018","6902021","6902022"]
cmapdef="binary"
plot_bathymetry=True

fig=plt.figure(figsize=(10,10))
#fig=plt.figure(figsize=(10,10))


plt.clf()
bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = 'i',fix_aspect=False)

#bmap = Basemap(projection='nsper',lon_0=lon_mean,lat_0=lat_mean,
#        satellite_height=500000.,resolution='h')

#bmap.drawcoastlines()
#bmap.fillcontinents()

#gludge to get rid of the lakes
for i,cp in enumerate(bmap.coastpolygons):
     if bmap.coastpolygontypes[i]<2:
         bmap.plot(cp[0],cp[1],'k-')
bmap.fillcontinents([0.9,0.9,0.9],lake_color=[0.85,0.85,0.85]) 
#bmap.drawmapboundary(fill_color='aqua')
#batym_levs=[0,10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,1000]
batym_levs=[0,20,40,60,80,100,120,140,160,180]
baty_min=0
baty_max=150
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
    topoin.mask[400:,1:200]=True #Remove the atlantic points from the data.
#    x=np.tile(lons,(lats.shape[0],1))
#    y=np.tile(lats,(lons.shape[0],1)).T
    tmp_lon,tmp_lat =bmap(*np.meshgrid(lons,lats))
    print "LON", tmp_lon.shape
    print "LAT", tmp_lat.shape
#    bmap.pcolor(tmp_lon,tmp_lat,-1*topoin,cmap='bone_r',vmin=0,vmax=150)
    bmap.contourf(tmp_lon,tmp_lat,-1*topoin,levels=batym_levs,cmap=cmapdef,vmin=baty_min,vmax=baty_max,extend='max')
    cb=plt.colorbar(fraction=0.027, pad=0.04)
#    cb.set_ticks(range(0,200,20))
#    plt.clim(0,200)
    cb.ax.invert_yaxis()
#ACTUAL DATA
    start=mp.dates.datetime.datetime(1000,5,5)
    end=mp.dates.datetime.datetime(3030,5,5)    
    for f,col,lab in zip(filest_to_plot,colors,labels):
        a=qa.PointData(f,1,start,end,"argonc");
        x,y=bmap(a.obs['ape']['lon'][:],a.obs['ape']['lat'][:])
    #    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
        bmap.plot(x[~x.mask],y[~y.mask],color=col,linewidth=2, alpha=0.8) #x,y masks are there to crop non values.
        bmap.plot(x[0],y[0],'*',color=col,markersize=12,alpha=1.0,zorder=10)    
        bmap.plot(x[-1],y[-1],'o',color=col,markersize=8,alpha=1.0,label=lab,zorder=10)    
#    plt.gca().set_facecolor('#000000')
    #plot the SR5 on the map
    plot_SR5(bmap)
    

    bmap.drawparallels(np.arange(45.,69,1.),labels=[1,0,0,0],linewidth=0)
    bmap.drawmeridians(np.arange(-5.,30,2.),labels=[0,0,0,1],linewidth=0)
#    plt.legend(bbox_to_anchor=(1.9,0.2),numpoints=1)
    plt.legend(numpoints=1)


plt.savefig('petranc//ArgoRoutesSmallArea2016.png' ,facecolor='w',dpi=600)
plt.savefig('petranc//ArgoRoutesSmallArea2016.pdf' ,facecolor='w',dpi=600)
plt.savefig('petranc//ArgoRoutesSmallArea2016.svg' ,facecolor='w',dpi=600)
plt.savefig('petranc//ArgoRoutesSmallArea2016.eps' ,facecolor='w',dpi=600)