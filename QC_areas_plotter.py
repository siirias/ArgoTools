# -*- coding: utf-8 -*-
"""
Created on Mon May  9 16:46:30 2022

@author: siirias
"""

import cartopy.crs as ccrs
# import cartopy.feature as cfeature
# from cartopy.feature import ShapelyFeature
# import cartopy.io.shapereader as shpreader

import matplotlib.pyplot as plt          # just for plots
import matplotlib.ticker as mticker      # just for plots
import numpy as np
# from PIL import Image

output_dir = "C:\\Data\\ArgoData\\Figures\\"
C_LAT = 59.5  
C_LON = 21.0  
Proj = ccrs.TransverseMercator(\
           central_latitude = C_LAT,\
           central_longitude = C_LON,\
           scale_factor = 0.001, approx=True) #kilometers, center at radar
Crs_LatLon = ccrs.PlateCarree()


MIN_LAT = C_LAT - 7.0
MAX_LAT = C_LAT + 7.0
MIN_LON = C_LON - 10.0
MAX_LON = C_LON + 10.0

lines = [
    [np.arange(19.0,22.3,0.1), 63.5*np.ones(np.arange(19.0,22.3,0.1).size)], #BB/BS
    [np.arange(19.0,22.9,0.1), 59.8*np.ones(np.arange(19.0,22.9,0.1).size)], #BS/NBP
    [np.arange(17.0,24.0,0.1), 58.4*np.ones(np.arange(17.0,24.0,0.1).size)], #NBP/BP
#    [18.6*np.ones(np.arange(56.2,58.4,0.1).size), np.arange(56.2,58.4,0.1)],#WBP/EBP
    [np.arange(16.0,21.0,0.1), 56.2*np.ones(np.arange(16.0,21.0,0.1).size)],#BP/BORNHOLM
    [17.5*np.ones(np.arange(56.2,54.7,-0.1).size), np.arange(56.2,54.7,-0.1)],
    [14.4*np.ones(np.arange(55.5,53.9,-0.1).size), np.arange(55.5,53.9,-0.1)],#Arkona/Bornholm
    ]

texts = [
    ["Bothnian Bay\n4.0",22.9,64.9],
    ["Bothnian Sea\n6.0",19.2,61.9],
    ["Northern\nBaltic Proper\n7.5",20.4,59.0],
    ["Baltic Proper\n8.0",20.3, 57.6],
    ["Bornholm Basin\n8.1",15.9, 55.2],
    ["Gdansk Basin\n8.1",19.2, 55.2],
    
    ]
plt.figure(figsize = (8,10))
ax = plt.axes(projection=Proj)
ax.coastlines(resolution = '10m')
ax.add_wms('https://geoserver2.ymparisto.fi/geoserver/eo/wms',\
           layers = ['HELCOM offshore sub-basins (2016)'],\
            alpha = 0.6)
#ax.add_feature(cfeature.GSHHSFeature())
ax.set_extent((MIN_LON,MAX_LON,MIN_LAT, MAX_LAT))
ax.set_aspect('auto')
for i in lines:
    plt.plot(i[0],i[1], transform = ccrs.PlateCarree(), linewidth = 2.0, color = 'blue')
for i in texts:
    plt.text(i[1], i[2], i[0], transform = ccrs.PlateCarree(), ha = 'center', color = 'blue')
gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True,\
                  linewidth=1, color='gray', alpha=0.5)
gl.xlocator = mticker.FixedLocator(np.arange(MIN_LON,MAX_LON,1.0))
gl.ylocator = mticker.FixedLocator(np.arange(MIN_LAT,MAX_LAT,0.5))
gl.xlabels_top = False
gl.ylabels_left = False        

figure_name = "QC-areas"
plt.savefig(output_dir+figure_name+'.png' ,\
            facecolor='w',dpi=300,bbox_inches='tight')
print("saved: {}".format(output_dir+figure_name+'.png'))

