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
import argohelper as ah

#Full baltic
#lon_min=12;lat_min=53;lon_max=30;lat_max=66;
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;

cmapdef="binary"
levdefs=[0,20,40,80,100,120,140,160,180,200,220,240,1000]
#levdefs=[0,20,40,60,80,100,120,140,160,180,200,220,240,1000]
#levdefs=[0,20,40,60,80,100,120,140,160,1000]
#Gulf of Bothnia
lon_min=-5;lat_min=45;lon_max=31;lat_max=66;

lon_mean=0.5*(lon_min+lon_max)
lat_mean=0.5*(lat_min+lat_max)

plot_bathymetry=True

fig=plt.figure(figsize=(14,7.5))
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
    topoin.mask[400:,1:200]=True #Remove the atlantic points from the data.
#    x=np.tile(lons,(lats.shape[0],1))
#    y=np.tile(lats,(lons.shape[0],1)).T
    tmp_lon,tmp_lat =bmap(*np.meshgrid(lons,lats))
    print "LON", tmp_lon.shape
    print "LAT", tmp_lat.shape
#    bmap.pcolor(tmp_lon,tmp_lat,-1*topoin,cmap='bone_r',vmin=0,vmax=150)
    bmap.contourf(tmp_lon,tmp_lat,-1*topoin,levels=levdefs,cmap=cmapdef,vmin=0,vmax=200)
    cb=plt.colorbar(fraction=0.027, pad=0.04)
    bmap.contour(tmp_lon,tmp_lat,-1*topoin,levels=levdefs,cmap=cmapdef+'_r',vmin=0,vmax=200,alpha=0.3)

#    cb.set_ticks(range(0,200,20))
#    plt.clim(0,200)
    cb.ax.invert_yaxis()
    
    
#HIGHLIGHT
    #Gotlands deep
    lon_min=18;lat_min=55;lon_max=21;lat_max=59;
    x,y=bmap([lon_min, lon_min, lon_max, lon_max, lon_min],[lat_min, lat_max, lat_max, lat_min,lat_min])
    bmap.plot(x,y,color='r',linewidth=2, alpha=0.8)

    zlon_max=8.7;zlon_min=-3.76;zlat_max=63.75;zlat_min=47.25    
    
    x,y=bmap([lon_max, zlon_max],[lat_max, zlat_max])
    bmap.plot(x,y,color='r',linewidth=2, alpha=0.4)
    x,y=bmap([lon_max, zlon_max],[lat_min, zlat_min])
    bmap.plot(x,y,color='r',linewidth=2, alpha=0.4)
    x,y=bmap([lon_min, zlon_min],[lat_max, zlat_max])
    bmap.plot(x,y,color='r',linewidth=2, alpha=0.2)
    x,y=bmap([lon_min, zlon_min],[lat_min, zlat_min])
    bmap.plot(x,y,color='r',linewidth=2, alpha=0.2)

#    x,y=bmap([zlon_min, zlon_min, zlon_max, zlon_max, zlon_min],[zlat_min, zlat_max, zlat_max, zlat_min,zlat_min])
#    bmap.plot(x,y,color='k',linewidth=8, alpha=0.5)

    bmap.drawparallels(np.arange(40.,69,2.),labels=[1,0,0,0],linewidth=0)
    bmap.drawmeridians(np.arange(-12.,30,2.),labels=[0,0,0,1],linewidth=0)

#ACTUAL DATA
    filest_to_plot=["6902014_prof.nc","6902019_prof.nc", \
                    "6902020_prof.nc"]
    ah.file_names_converted
    start=mp.dates.datetime.datetime(1000,5,5)
    end=mp.dates.datetime.datetime(3030,5,5)
                    
    colors=["#fa8072","#ef9aef","#fff5ee"]
    labels=["6902014","6902019", \
            "6902020"]
    
    for f,col,lab in zip(filest_to_plot,colors,labels):
        a=qa.PointData(f,1,start,end,"argonc");
        x,y=bmap(a.obs['ape']['lon'][:],a.obs['ape']['lat'][:])
    #    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
        bmap.plot(x,y,color=col,linewidth=2, alpha=0.8)
        bmap.plot(x[-1],y[-1],'o',color=col,markersize=10,alpha=1.0,label=lab)    
    plt.legend(bbox_to_anchor=(0.75,0.3),numpoints=1)
    #plt.legend(numpoints=1)

#ANOTHER Plot inside this one
    sub_a=plt.axes([.15,.15,.25,.7])
    plt.xticks([])
    plt.yticks([])
    #Gotlands deep
    lon_min=18;lat_min=55;lon_max=21;lat_max=59;
    
    bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
    resolution = 'i',fix_aspect=True)
    bmap.drawcoastlines()
    tmp_lon,tmp_lat =bmap(*np.meshgrid(lons,lats))
    print "LON", tmp_lon.shape
    print "LAT", tmp_lat.shape
    bmap.contourf(tmp_lon,tmp_lat,-1*topoin,levels=levdefs,cmap=cmapdef,vmin=0,vmax=200)
    bmap.contour(tmp_lon,tmp_lat,-1*topoin,levels=levdefs,cmap=cmapdef+'_r',vmin=0,vmax=200,alpha=0.3)
    start=mp.dates.datetime.datetime(1000,5,5)
    end=mp.dates.datetime.datetime(3030,5,5)    
    for f,col,lab in zip(filest_to_plot,colors,labels):
        a=qa.PointData(f,1,start,end,"argonc");
        x,y=bmap(a.obs['ape']['lon'][:],a.obs['ape']['lat'][:])
    #    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
        bmap.plot(x,y,color=col,linewidth=2, alpha=0.8)
        bmap.plot(x[-1],y[-1],'o',color=col,markersize=10,alpha=1.0,label=lab)    


plt.savefig('GulfOfBothnia_routes.png' ,facecolor='w',dpi=600)
plt.savefig('GulfOfBothnia_routes.pdf' ,facecolor='w',dpi=600)
plt.savefig('GulfOfBothnia_routes.svg' ,facecolor='w',dpi=600)
plt.savefig('GulfOfBothnia_routes.eps' ,facecolor='w',dpi=600)

