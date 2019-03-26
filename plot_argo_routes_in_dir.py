# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:11:08 2016

@author: siirias
"""
import sys
sys.path.insert(0,'D:\\svnfmi_merimallit\\qa\\nemo')
import os
import re
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import ModelQATools as qa
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset


dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\Riikan\\"
figure_name="argos_in_baltic"
line_width=0.7
line_alpha=0.8
marker_end_size=5
marker_start_size=5
legend_size=10
#Full baltic
lon_min=14;lat_min=53;lon_max=30;lat_max=66;
figure_size=(9,10)

#Full baltic (Finnish areas)
#lon_min=16;lat_min=53;lon_max=30;lat_max=66;
#figure_size=(8,10)
#Gotlands deep
#lon_min=18;lat_min=55;lon_max=21;lat_max=59;
#Bothnian Sea
#lon_min=17;lat_min=60;lon_max=22;lat_max=63;
all_colors= ["#ff0000","#00ff00","#0000ff",\
             "#000000","#d00000","#d000d0",\
             "#005050","#888888","#ffff00",\
             "#ff00ff", "#00ffff","#600000",\
             "#aa0055", "#50ff50", "#ff5050",\
             "#5050ff", "#505000", "#500050",\
             "#005050", "#50ff00", "#ff5000"]



plot_bathymetry=True
plot_legends=True
plot_routes=True
start=mp.dates.datetime.datetime(1000,5,5)
end=mp.dates.datetime.datetime(3030,5,5)

fig=plt.figure(figsize=(8,10))
plt.clf()
bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = 'i',fix_aspect=False)

bmap.drawcoastlines()
bmap.fillcontinents()
bmap.drawparallels(np.arange(50.,69,2.),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
bmap.drawmeridians(np.arange(12.,30,2.),labels=[0,0,0,1],linewidth=0,dashes=[5,10])



files_to_plot=os.listdir(dir_to_plot)
labels= map(lambda x: re.search('\d{7}',files_to_plot[x]).group(0),range(len(files_to_plot)))                
colors=all_colors[0:len(files_to_plot)]
               

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

if plot_routes:
    for f,col,lab in zip(files_to_plot,colors,labels):
        a=qa.PointData(dir_to_plot+f,1,start,end,"argonc");
        lon_dat=a.obs['ape']['lon'][~a.obs['ape']['lon'].mask]
        lat_dat=a.obs['ape']['lat'][~a.obs['ape']['lat'].mask]
        x,y=bmap(lon_dat,lat_dat)
    #    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
        bmap.plot(x,y,color=col,linewidth=line_width, alpha=line_alpha)
        bmap.plot(x[-1],y[-1],'x',color=col,markersize=marker_end_size,alpha=1.0)
        bmap.plot(x[0],y[0],'o',color=col,markersize=marker_start_size,alpha=1.0,label=lab)
    #    print lab, mp.dates.num2date(a.obs['ape']['date'][0]).date() \
    #             , mp.dates.num2date(a.obs['ape']['date'][-1]).date()

if plot_legends:
    #plt.legend(bbox_to_anchor=(1.0,0.5),numpoints=1)
    plt.legend(loc='lower right',numpoints=1,prop={'size': legend_size})
plt.savefig(figure_name+'.png' ,facecolor='w',dpi=300)
plt.savefig(figure_name+'.eps' ,facecolor='w',dpi=300)
