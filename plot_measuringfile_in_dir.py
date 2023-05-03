# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:11:08 2016
Plans to become way to plot point where measurement data has been taken for given area.
@author: siirias
"""
import sys
import os
import re
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from netCDF4 import Dataset
import xarray as xr
import argohelper as ah
import cmocean as cmo
from itertools import cycle
import pandas as pd

file_to_plot="C:\\Data\\EARiseQC\\FTP\\FMI\\ICESCTD00-20.csv" #default value
lat_i = 'Latitude [degrees_north]'
lon_i = 'Longitude [degrees_east]'

output_dir = "C:\\Data\\ArgoData\\Figures\\"
data_dir = "C:\\Data\\ArgoData\\"  # mainly for topography data
figure_setup = "DMQCAreas" #"EAR_UseCase" #"EARISE_deployment"#"Bothnian Sea Aranda" # "Bothnian Sea Aranda" # "GotlandD"#May change dir_to_plot
figure_name="DMQCAreas"  #default value
do_plot = True
plot_contours = False  # default. specific etups may change this
draw_labels = False
draw_EEZ = True
draw_measuring_points = False
contour_levels = [50,100,150,200,250,300]
shore_resolution = "10m"  # "10m" "50m"
fig_dpi = 300
line_width = 1.2  #0.7
line_alpha = 0.8
marker_end_size = 5
marker_start_size = 5
marker_size = 5
legend_size = 10
label_step = 2.0
bathy_max = 300 # meters
the_proj = ccrs.PlateCarree()
requested_proj = ccrs.PlateCarree()
requested_aspect = 'equal'

replace_labels = {}
all_colors= ["#ff0000","#000000","#0000ff",\
             "#00ff00","#007060","#d000d0",\
             "#d00000","#888888","#ffff00",\
             "#ff00ff", "#00ffff","#600000",\
             "#aa0055", "#50ff50", "#ff5050",\
             "#5050ff", "#505000", "#500050",\
             "#005050", "#50ff00", "#ff5000"]
bathy_colormap = cmo.cm.deep
plot_bathymetry=True
plot_legends=False
plot_routes=True
plot_points = True
start=mp.dates.datetime.datetime(2000,3,1)
end=mp.dates.datetime.datetime(2230,5,5)
figure_size=(10,5)  #default value!

    
if( figure_setup == "DMQCAreas"):
    figure_name = "DMQCAreas"
    lon_min=12;lat_min=53;lon_max=30.5;lat_max=66;
    figure_size=(10,10)
    C_LAT = (lat_min+lat_max)/2.0
    C_LON = (lon_min+lon_max)/2.0
    requested_proj = ccrs.TransverseMercator(\
               central_latitude = C_LAT,\
               central_longitude = C_LON,\
           scale_factor = 0.001, approx=True) #kilometers, center at radar
    #requested_proj = ccrs.PlateCarree()
if draw_measuring_points:
    data = pd.read_csv(file_to_plot, delimiter='\t')


if do_plot:
    fig=plt.figure(figsize=figure_size)
    plt.clf()
    
    ax = plt.axes(projection=requested_proj)
    ax.set_extent([lon_min,lon_max,lat_min,lat_max])
    ax.set_aspect(requested_aspect)
    ax.coastlines(shore_resolution)
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', shore_resolution,\
                                            edgecolor='#000000', linewidth = 0.1,\
                                            facecolor='#ccccdd', alpha = 1.0))
    ax.add_feature(cfeature.NaturalEarthFeature('physical', 'lakes', shore_resolution,\
                                            edgecolor='#000000', linewidth=0.5, \
                                            facecolor='#ffffff', alpha = 1.0))
    ax.add_feature(cfeature.NaturalEarthFeature('cultural', 'admin_0_boundary_lines_land', shore_resolution,\
                                            edgecolor='#600000', linewidth=0.7, \
                                            facecolor = 'none', alpha = 0.4, zorder = 5))
    grid_proj = ccrs.PlateCarree()
    gl = ax.gridlines(crs=grid_proj, draw_labels=draw_labels,
              linewidth=2, color='gray', alpha=0.1, linestyle='-')
    gl.xlabels_top = True
    gl.ylabels_right = True
    gl.top_labels = False
    gl.right_labels = False
    
    
        
    
    #TOPOGRAPHY
    if plot_bathymetry:
        topodata = Dataset(data_dir+'iowtopo2_rev03.nc')
        topoin = topodata.variables['Z_WATER'][:]
        lons = topodata.variables['XT_I'][:]
        lats = topodata.variables['YT_J'][:]
        x=np.tile(lons,(lats.shape[0],1))
        y=np.tile(lats,(lons.shape[0],1)).T
    
        if plot_contours:
    #        cn = bmap.contour(x,y,-1*topoin,colors='k',vmin=0,vmax=bathy_max, alpha=0.3)
            cn = plt.contour(x,y,-1*topoin,levels = contour_levels,\
                             colors='k',vmin=0,vmax=bathy_max,\
                             alpha=0.3,transform = ccrs.PlateCarree())
            plt.clabel(cn,fmt='%1.0f')
    #    bmap.pcolor(x,y,-1*topoin,cmap=cmo.cm.deep,vmin=0,vmax=bathy_max)
        plt.pcolor(x,y,-1*topoin,cmap=bathy_colormap,vmin=0,\
                   vmax=bathy_max, transform = ccrs.PlateCarree(), zorder = -1)
        cb=plt.colorbar()
        cb.ax.invert_yaxis()
        cb.set_label('Depth (m)')
        # gl.xlabels_top = False
        # gl.ylabels_right = False
       
    if(draw_EEZ):
        ax.add_wms('http://geo.vliz.be/geoserver/MarineRegions/wms?',\
                   layers='eez_boundaries', alpha = 1.0)
            
    if draw_measuring_points:
        plt.plot(data[lon_i], data[lat_i], '.', transform = ccrs.PlateCarree(), zorder = 4, color = 'b', markersize= 3.0,alpha=0.002)
    
        
    if plot_legends:
        plt.legend(loc='lower right',numpoints=1,prop={'size': legend_size})
    plt.savefig(output_dir+figure_name+'.png' ,\
                facecolor='w',dpi=fig_dpi,bbox_inches='tight')
    print("saved: {}".format(output_dir+figure_name+'.png'))
    plt.savefig(output_dir+figure_name+'.eps' ,\
                facecolor='w',dpi=fig_dpi,bbox_inches='tight')
