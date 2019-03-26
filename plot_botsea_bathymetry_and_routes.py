# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:11:08 2016

@author: siirias
"""
from __future__ import unicode_literals

#Full baltic

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
import pandas as pd
import string
import re
#lon_min=12;lat_min=53;lon_max=30;lat_max=66;
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;
def plot_SR5(bmap,size_multiplier=1.0):
    x,y=bmap([19.3555],[61.0520])
    bmap.plot(x,y,'X',color="#ff2020",markersize=8*size_multiplier,alpha=1.0,zorder=20,label="SR5")    

make_inset=True
#Gulf of Bothnia
lon_min=-5;lat_min=53;lon_max=31;lat_max=66;

lon_mean=0.5*(lon_min+lon_max)
lat_mean=0.5*(lat_min+lat_max)
filest_to_plot=["petranc//6901901_prof.nc","petranc//6902013_prof.nc","petranc//6902017_prof.nc", \
                "petranc//6902018_prof.nc","petranc//6902021_prof.nc", \
                "petranc//6902022_prof.nc"]
                
colors=["#fa8072","#ef9aef","#fff550","#4682b4","#00ffff","#00ff00"]
labels=["6901901","6902013","6902017", \
        "6902018","6902021","6902022"]
cmapdef="binary"
plot_bathymetry=True

cities=[{'name':'Pori', 'lat':61.4851,'lon':21.7974},
        {'name':'Rauma', 'lat':61.1309,'lon':21.5059},
        {'name':'Sundvall', 'lat':62.39008,'lon':17.3069},
        {'name':'Gävle', 'lat':60.6749,'lon':17.1413}
        ]

Figure2plot=1
if(Figure2plot==1):
    fig=plt.figure(figsize=(14,10))
if(Figure2plot  in [2,3]):
    make_inset=False
    lon_min=19.00;lat_min=60.50;lon_max=22.0;lat_max=62.0; #close zoom
    if(Figure2plot==3):
        lon_min=17.00;lat_min=60.50;lon_max=22.0;lat_max=63.0; #as reviewer suggester
    lon_mean=0.5*(lon_min+lon_max)
    lat_mean=0.5*(lat_min+lat_max)
    fig=plt.figure(figsize=(10,10))

res='h' #'h','c'
plt.clf()
bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = res,fix_aspect=False) #h

#bmap = Basemap(projection='nsper',lon_0=lon_mean,lat_0=lat_mean,
#        satellite_height=500000.,resolution='h')

#bmap.drawcoastlines()
#bmap.fillcontinents()

#gludge to get rid of the lakes
for i,cp in enumerate(bmap.coastpolygons):
     if bmap.coastpolygontypes[i]<2:
         bmap.plot(cp[0],cp[1],color="#505050",linewidth=0.5)
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
    bmap.drawcountries()
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
        if(hasattr(x,'mask')):
            x=x[~x.mask]
            y=y[~y.mask]
        bmap.plot(x,y,color=col,linewidth=2, alpha=0.8) #x,y masks are there to crop non values.
        bmap.plot(x[0],y[0],'*',color=col,markersize=12,alpha=1.0,zorder=10)    
        bmap.plot(x[-1],y[-1],'o',color=col,markersize=8,alpha=1.0,label=lab,zorder=10)    
#    plt.gca().set_facecolor('#000000')
    #plot the SR5 on the map
    plot_SR5(bmap)
    if(Figure2plot==1):
        bmap.drawparallels(np.arange(45.,69,2.),labels=[1,0,0,0],linewidth=0)
        bmap.drawmeridians(np.arange(-5.,30,2.),labels=[0,0,0,1],linewidth=0)
    if(Figure2plot in [2,3]):
        bmap.drawparallels(np.arange(45.,69,0.5),labels=[1,0,0,0],linewidth=0)
        bmap.drawmeridians(np.arange(-5.,30,0.5),labels=[0,0,0,1],linewidth=0)

    #Let's plot some country names:
    countries=pd.read_csv('countries.txt','\t')
    for count in range(len(countries)):
        i=countries.iloc[count]
        if(i['lat']>lat_min and i['lat']<lat_max and i['lon']>lon_min and i['lon']<lon_max):
            x,y=bmap([i['lon']],[i['lat']])
            name=re.sub('[^'+string.printable+']',' ',i['name'])
            plt.text(x[0],y[0],name)
    #Oh and why not some cities too while we are at it:
    if Figure2plot in [2,3]:
        for i in cities:
            if(i['lat']>lat_min and i['lat']<lat_max and i['lon']>lon_min and i['lon']<lon_max):
                x,y=bmap([i['lon']],[i['lat']])
                name=i['name']
                plt.text(x[0],y[0]-0.02,'\n   '+name,va='center',ha='center')
                plt.plot([x],[y],'.k',markersize=13)
    if(Figure2plot==3):
        lon_min=19.97;lat_min=61.331;lon_max=20.523;lat_max=61.862;
        x,y=bmap([lon_min, lon_min, lon_max, lon_max, lon_min],[lat_min, lat_max, lat_max, lat_min,lat_min])
        bmap.plot(x,y,color='r',linewidth=2, alpha=0.8)
    if(Figure2plot in [2]):
        plt.legend(numpoints=1)
    if(Figure2plot in [3]):
        plt.legend(loc="upper right", numpoints=1)
    
    if make_inset:
    #ANOTHER Plot inside this one
#HIGHLIGHT
        #Baltic Sea deep
        lon_min=17.00;lat_min=60.55;lon_max=21.85;lat_max=62.85;
        x,y=bmap([lon_min, lon_min, lon_max, lon_max, lon_min],[lat_min, lat_max, lat_max, lat_min,lat_min])
        bmap.plot(x,y,color='r',linewidth=2, alpha=0.8)
    
        zlon_max=13.9;zlon_min=-3.76;zlat_max=64.35;zlat_min=56.2    
        
        x,y=bmap([lon_max, zlon_max],[lat_max, zlat_max])
        bmap.plot(x,y,color='r',linewidth=2, alpha=0.4)
        x,y=bmap([lon_max, zlon_max],[lat_min, zlat_min])
        bmap.plot(x,y,color='r',linewidth=2, alpha=0.4)
        x,y=bmap([lon_min, zlon_min],[lat_max, zlat_max])
        bmap.plot(x,y,color='r',linewidth=2, alpha=0.2)
        x,y=bmap([lon_min, zlon_min],[lat_min, zlat_min])
        bmap.plot(x,y,color='r',linewidth=2, alpha=0.2)
    
        x,y=bmap([zlon_min, zlon_min, zlon_max, zlon_max, zlon_min],[zlat_min, zlat_max, zlat_max, zlat_min,zlat_min])
        bmap.plot(x,y,color='k',linewidth=8, alpha=0.5)
    
        sub_a=plt.axes([.15,.3,.35,.5])
        plt.xticks([])
        plt.yticks([])
        #Gotlands deep
    #    lon_min=18;lat_min=55;lon_max=21;lat_max=59;
        
        bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
        resolution = res,fix_aspect=False)
        bmap.drawcoastlines(color="#505050",linewidth=0.5)
        bmap.drawcountries()
        bmap.fillcontinents([0.7,0.7,0.7],lake_color=[0.85,0.85,0.85]) 
        tmp_lon,tmp_lat =bmap(*np.meshgrid(lons,lats))
        print "LON", tmp_lon.shape
        print "LAT", tmp_lat.shape
        bmap.contourf(tmp_lon,tmp_lat,-1*topoin,levels=batym_levs,cmap=cmapdef,vmin=baty_min,vmax=baty_max,extend='max')
        bmap.contour(tmp_lon,tmp_lat,-1*topoin,levels=batym_levs,cmap=cmapdef+'_r',vmin=baty_max,vmax=baty_max,alpha=0.3)
    #ACTUAL DATA
        #plot the SR5 on the map
        start=mp.dates.datetime.datetime(1000,5,5)
        end=mp.dates.datetime.datetime(3030,5,5)    
        for f,col,lab in zip(filest_to_plot,colors,labels):
            a=qa.PointData(f,1,start,end,"argonc");
            x,y=bmap(a.obs['ape']['lon'][:],a.obs['ape']['lat'][:])
        #    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
            if(hasattr(x,'mask')):
                x=x[~x.mask]
                y=y[~y.mask]
            bmap.plot(x,y,color=col,linewidth=2, alpha=0.8)#x,y masks are there to crop non values.
            bmap.plot(x[0],y[0],'*',color=col,markersize=12,alpha=1.0,zorder=10)    
            bmap.plot(x[-1],y[-1],'o',color=col,markersize=8,alpha=1.0,label=lab,zorder=10)    
        plot_SR5(bmap,1.3)
        plt.legend(bbox_to_anchor=(1.9,0.2),numpoints=1)


if(Figure2plot==1):
    plt.savefig('petranc//ArgoRoutes2016.png' ,facecolor='w',dpi=600)
    plt.savefig('petranc//ArgoRoutes2016.pdf' ,facecolor='w',dpi=600)
    plt.savefig('petranc//ArgoRoutes2016.svg' ,facecolor='w',dpi=600)
    plt.savefig('petranc//ArgoRoutes2016.eps' ,facecolor='w',dpi=600)
if(Figure2plot==2):
    plt.savefig('petranc//ArgoRoutesSmallArea2016.png' ,facecolor='w',dpi=600)
    plt.savefig('petranc//ArgoRoutesSmallArea2016.pdf' ,facecolor='w',dpi=600)
    plt.savefig('petranc//ArgoRoutesSmallArea2016.svg' ,facecolor='w',dpi=600)
    plt.savefig('petranc//ArgoRoutesSmallArea2016.eps' ,facecolor='w',dpi=600)
if(Figure2plot==3):
    plt.savefig('petranc//ArgoRoutesSmallArea2016b.png' ,facecolor='w',dpi=600)
    plt.savefig('petranc//ArgoRoutesSmallArea2016b.pdf' ,facecolor='w',dpi=600)
    plt.savefig('petranc//ArgoRoutesSmallArea2016b.svg' ,facecolor='w',dpi=600)
    plt.savefig('petranc//ArgoRoutesSmallArea2016b.eps' ,facecolor='w',dpi=600)
