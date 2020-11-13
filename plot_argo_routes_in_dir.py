# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:11:08 2016

@author: siirias
"""
import sys
import os
import re
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import xarray as xr
import argohelper as ah
import cmocean as cmo
from itertools import cycle

dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\" #default value
output_dir = "D:\\Data\\ArgoData\\Figures\\"
data_dir = "D:\\Data\\ArgoData\\"  # mainly for topography data
figure_setup = "GoB"#"Bothnian Sea Aranda" # "Bothnian Sea Aranda" # "GotlandD"#May change dir_to_plot
#figure_setup ="Bothnian Sea"  #"EARISE_BP" #May change dir_to_plot
figure_name="ArgoPlot"
plot_contours = False  # default. specific etups may change this
fig_dpi = 300
line_width = 1.2  #0.7
line_alpha = 0.8
marker_end_size = 5
marker_start_size = 5
marker_size = 5
legend_size = 10
label_step = 2.0
bathy_max = 300 # meters
replace_labels = {}
all_colors= ["#ff0000","#000000","#0000ff",\
             "#00ff00","#007060","#d000d0",\
             "#d00000","#888888","#ffff00",\
             "#ff00ff", "#00ffff","#600000",\
             "#aa0055", "#50ff50", "#ff5050",\
             "#5050ff", "#505000", "#500050",\
             "#005050", "#50ff00", "#ff5000"]
plot_bathymetry=True
plot_legends=True
plot_routes=True
plot_points = True
start=mp.dates.datetime.datetime(2000,3,1)
end=mp.dates.datetime.datetime(2230,5,5)
figure_size=(10,5)  #default value!

if( figure_setup == "GoB"):
    dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\GotlandDeep\\"
    figure_name = "Gulf of Bothnia"
    lon_min=17;lat_min=60;lon_max=26;lat_max=66;
    figure_size=(10,10)
    
if(figure_setup == "FullBaltic"):
    lon_min=16;lat_min=53;lon_max=30.5;lat_max=66;
    figure_size=(8,10)
    
if(figure_setup == "FullBalticEAR"):
    dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\EARise\\" 
    line_alpha = 0.4
    plot_points = False
    plot_legends = False
    bathy_max = 400 # meters
    lon_min=10;lat_min=53;lon_max=30.5;lat_max=66;
    figure_size=(12,10)
    all_colors = ['#000000']
                  
if(figure_setup == "AllFinnish"):
    lon_min=10;lat_min=53;lon_max=30.5;lat_max=66;
    figure_name = 'AllFinnishFloats'
    figure_size=(12,10)
    dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\AllFinnish\\"
    
if(figure_setup == "GotlandD"):
    lon_min=18;lat_min=56;lon_max=21;lat_max=59;
    start=mp.dates.datetime.datetime(2018,3,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_size=(5,5)  #default value!
    figure_name="FMIGotlandDeep"
    dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\GotlandDeep\\"
    
if(figure_setup == "Bothnian Sea"):
    start=mp.dates.datetime.datetime(2019,3,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_name="FMIBothnianSea"
    dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\BothnianSea\\"
    lon_min=17;lat_min=60;lon_max=22;lat_max=63;
    
if(figure_setup == "Bothnian Sea Aranda"):
    start=mp.dates.datetime.datetime(2019,3,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_name="FMIBothnianSeaAranda"
    dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\BothnianSea\\"
    lon_min=19;lat_min=61;lon_max=22;lat_max=63;
    plot_contours = True
    figure_size=(10,10)
    line_alpha=0.5
    label_step = 0.4
    bathy_max = 200.0
    marker_end_size = 7
    marker_start_size = 5
    marker_size = 5

if(figure_setup == "Bay of Bothnia"):
    start=mp.dates.datetime.datetime(2000,3,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_name="FMIBothnianBay"
    dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\BayOfBothnia\\"
    lon_min=20;lat_min=64;lon_max=26;lat_max=66;
    
if(figure_setup == "EARISE_BP"):
    figure_name="EuroArgoRISE"
    dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\" 
    lon_min=19.0;lat_min=58.3;lon_max=21.2;lat_max=59.2;
    plot_contours = True
    replace_labels = {'6903703':'ARVOR-I(6903703)',\
                      '6903704':'APEX(6903704)'}
    figure_size=(10,5)
    line_alpha=0.3
    label_step = 0.2
    bathy_max = 200.0
    marker_end_size = 10
    marker_start_size = 5
    marker_size = 5




fig=plt.figure(figsize=figure_size)
plt.clf()
bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
resolution = 'i',fix_aspect=False)

bmap.drawcoastlines(linewidth=0.5)
bmap.fillcontinents()
bmap.drawparallels(np.arange(50.,69,label_step),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
bmap.drawmeridians(np.arange(12.,30,label_step),labels=[0,0,0,1],linewidth=0,dashes=[5,10])



files_to_plot=[i for i in os.listdir(dir_to_plot) if i.endswith('.nc')]
labels= map(lambda x: re.search('\d{7}',files_to_plot[x]).group(0),range(len(files_to_plot)))                
colors=all_colors[0:len(files_to_plot)]
               

#TOPOGRAPHY EXPERIMENT
if plot_bathymetry:
    topodata = Dataset(data_dir+'iowtopo2_rev03.nc')
    
    topoin = topodata.variables['Z_WATER'][:]
    lons = topodata.variables['XT_I'][:]
    lats = topodata.variables['YT_J'][:]
    x=np.tile(lons,(lats.shape[0],1))
    y=np.tile(lats,(lons.shape[0],1)).T
    if plot_contours:
        cn = bmap.contour(x,y,-1*topoin,colors='k',vmin=0,vmax=bathy_max, alpha=0.3)
        plt.clabel(cn,fmt='%1.0f')
    bmap.pcolor(x,y,-1*topoin,cmap=cmo.cm.deep,vmin=0,vmax=bathy_max)
    cb=plt.colorbar()
    cb.ax.invert_yaxis()
    cb.set_label('Depth (m)')
color_stack = colors[:]
if plot_routes:
    for f,label in zip(files_to_plot,labels):
        if(len(color_stack)==0):
            color_stack = colors[:]
        lab = label
        if(lab in replace_labels.keys()):
            lab = replace_labels[lab]  #some image setups want specific labels
        d=xr.open_dataset(dir_to_plot+f)
        primaries = ah.get_primary_indices(d)
        primaries = np.asarray(primaries) & \
                    np.asarray(d['JULD']>np.datetime64(start)) &\
                    np.asarray(d['JULD']<np.datetime64(end))
        lat_dat = np.array(d['LATITUDE'])[primaries]
        lon_dat = np.array(d['LONGITUDE'])[primaries]
            
        x,y=bmap(lon_dat,lat_dat)
    #    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
        if(len(x)>0):
            col = color_stack.pop(0)
            bmap.plot(x,y,color=col,linewidth=line_width, alpha=line_alpha)
            bmap.plot(x[-1],y[-1],'x',color=col,markersize=marker_end_size,alpha=1.0)
            bmap.plot(x[0],y[0],'o',color=col,markersize=marker_start_size,alpha=1.0,label=lab)
            if(plot_points):
                bmap.plot(x,y,'.',color=col,markersize=marker_size,alpha=line_alpha)
            print(lon_dat[-1],lat_dat[-1], lab)
    #    print lab, mp.dates.num2date(a.obs['ape']['date'][0]).date() \
    #             , mp.dates.num2date(a.obs['ape']['date'][-1]).date()

if plot_legends:
    #plt.legend(bbox_to_anchor=(1.0,0.5),numpoints=1)
    plt.legend(loc='lower right',numpoints=1,prop={'size': legend_size})
plt.savefig(output_dir+figure_name+'.png' ,\
            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
plt.savefig(output_dir+figure_name+'.eps' ,\
            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
