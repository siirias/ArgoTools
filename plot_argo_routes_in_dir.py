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
#from mpl_toolkits.basemap import Basemap
import cartopy.crs as ccrs
import cartopy.feature as cfeature

from netCDF4 import Dataset
import xarray as xr
import argohelper as ah
import cmocean as cmo
from itertools import cycle

import folium # only needed for interactive leaflet maps
from PIL import Image # only needed for interactive leaflet maps

dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\" #default value
output_dir = "C:\\Data\\ArgoData\\Figures\\"
data_dir = "C:\\Data\\ArgoData\\"  # mainly for topography data
figure_setup = "NBalticProper" #"EAR_UseCase" #"EARISE_deployment"#"Bothnian Sea Aranda" # "Bothnian Sea Aranda" # "GotlandD"#May change dir_to_plot
#figure_setup ="Bothnian Sea"  #"EARISE_BP" #May change dir_to_plot

make_leaflet = True

figure_name="ArgoPlot"  #default value
plot_contours = False  # default. specific etups may change this
draw_labels = True
draw_EEZ = False
contour_levels = [50,100,150,200,250,300]
shore_resolution = "50m"  # "10m" "50m"
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

# Proj = ccrs.TransverseMercator(\
#            central_latitude = 78.0,\
#            central_longitude = 30.0,\
#            approx=True)
requested_aspect = 'auto'

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
plot_legends=True
plot_routes=True
plot_points = True
start=mp.dates.datetime.datetime(2000,3,1)
end=mp.dates.datetime.datetime(2230,5,5)
figure_size=(10,5)  #default value!
if(figure_setup == "ArvorC"):
    start=mp.dates.datetime.datetime(2019,3,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_name="FMI_ArvorC"
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\arvorc\\"
    lon_min=17;lat_min=60;lon_max=22;lat_max=63;
#    lon_min=20.1;lat_min=61.395;lon_max=20.25;lat_max=61.425;
    figure_size=(10,11)
    shore_resolution = "10m"  # "10m" "50m"
    

if( figure_setup == "BGC_BS"):
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BGC_BS\\"
    figure_name = "BGC_BothnianSea"
    lon_min=17;lat_min=60;lon_max=24;lat_max=64;
    figure_size=(10,10)

if( figure_setup == "BGC_BP"):
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BGC_BP\\"
    figure_name = "BGC_BalticProper"
    lon_min=16;lat_min=55;lon_max=23;lat_max=60;
    figure_size=(10,10)

if( figure_setup == "RBR"):
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\RBR\\"
    figure_name = "RBR_BalticProper"
    lon_min=16;lat_min=53;lon_max=30.5;lat_max=66;
    figure_size=(10,10)

if(figure_setup == "NBalticProper"):
    figure_name="NorthernBalticProper"
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\NBalticProper\\" 
#    lon_min=19.0-4.0;lat_min=58.0-2.0;lon_max=21.0+4.0;lat_max=59.8+2.0;
#    lon_min=19.15;lat_min=58.5;lon_max=20.8;lat_max=59.5;
    lon_min=19.6;lat_min=57.0;lon_max=22.2;lat_max=59.5;
    plot_contours = True
    contour_levels = list(range(0,250,25))
    replace_labels = {'6903703':'ARVOR-I(6903703)',\
                      '6903704':'APEX(6903704)'}
    figure_size=(10,5)
    line_alpha=0.3
    label_step = 0.2
    bathy_max = 200.0
    marker_end_size = 10
    marker_start_size = 5
    marker_size = 5
    line_width = 1.5  #0.7
    line_alpha = 0.9
    all_colors = ['#ff0000','#ffffff']
    requested_aspect = 1.0

if( figure_setup == "GoB"):
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BothnianSea\\"
    figure_name = "Gulf of Bothnia"
    lon_min=17;lat_min=60;lon_max=26;lat_max=66;
    figure_size=(10,10)
    
if(figure_setup == "FullBaltic"):
    lon_min=16;lat_min=53;lon_max=30.5;lat_max=66;
    figure_size=(8,10)
    
if(figure_setup == "FullBalticEAR"):
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\EARise\\" 
    line_alpha = 0.4
    plot_points = False
    plot_legends = False
    bathy_max = 400 # meters
    lon_min=10;lat_min=53;lon_max=30.5;lat_max=66;
    figure_size=(12,10)
    all_colors = ['#000000']

if(figure_setup == "EAR_UseCase"):
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\EARise_UseCase\\" 
    line_alpha = 0.4
    plot_points = False
    plot_legends = False
    bathy_max = 400 # meters
    lon_min=10;lat_min=53;lon_max=30.5;lat_max=66;
    figure_size=(12,10)
    all_colors = ['#000000']
    shore_resolution = "10m"  # "10m" "50m"
    fig_dpi  =300
    line_width = 1.0  #0.7


if(figure_setup == "Barents Sea"):
    figure_name = 'FMI_Barents_Sea_Argos'
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BarentsSea\\" 
    line_alpha = 0.5
    plot_points = True
    plot_legends = True
    plot_bathymetry = False
    bathy_max = 400 # meters
    lon_min=10;lat_min=75;lon_max=50.0;lat_max=80.0;
    center = [(lon_min+lon_max)*0.5, (lat_min+lat_max)*0.5]
    requested_proj = ccrs.LambertAzimuthalEqualArea(center[0],center[1])
    figure_size=np.array((12,10))*0.5
    marker_end_size = 5
    marker_start_size = 5
    marker_size = 5
    all_colors= ["#ff0000","#000000","#0000ff"]
    bathy_colormap = 'gist_gray_r'


if(figure_setup == "Barents Sea all"):
    figure_name = 'EARise_Barents_Sea_Argos'
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BarentsSeaAll\\" 
    line_alpha = 0.5
    plot_points = True
    plot_legends = True
    plot_bathymetry = False
    bathy_max = 400 # meters
    lon_min=10;lat_min=68;lon_max=50.0;lat_max=80.0;
    center = [(lon_min+lon_max)*0.5, (lat_min+lat_max)*0.5]
    requested_proj = ccrs.LambertAzimuthalEqualArea(center[0],center[1])
    figure_size=np.array((12,10))*0.75
    marker_end_size = 5
    marker_start_size = 5
    marker_size = 5
    bathy_colormap = 'gist_gray_r'


if(figure_setup == "EARISE_deployment"):
    figure_name = 'EARise_deployment'
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\EARise_deployment\\" 
    line_alpha = 0.5
    plot_points = False
    plot_legends = False
    bathy_max = 400 # meters
    lon_min=15;lat_min=55;lon_max=30.5;lat_max=66;
    figure_size=(12,10)
    marker_end_size = 0
    marker_start_size = 0
    marker_size = 0
    all_colors = ['#303040']
    bathy_colormap = 'gist_gray_r'
                  
if(figure_setup == "BS_ISA"):
    lon_min=16;lat_min=60;lon_max=25.5;lat_max=66;
    figure_name = 'BalticSeaISA'
    figure_size=np.array((12,10))*0.5
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BS_ICE\\"
    
if(figure_setup == "AllFinnish"):
    lon_min=10;lat_min=53;lon_max=30.5;lat_max=66;
    figure_name = 'AllFinnishFloats'
    figure_size=(12,10)
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\AllFinnish\\"
    
if(figure_setup == "NationalReport2020"):
    lon_min=10;lat_min=53;lon_max=30.5;lat_max=66;
    start=mp.dates.datetime.datetime(2010,1,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_name = 'NationalReport2020'
    figure_size=(12,9)
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\NationalReport2020\\"
    bathy_colormap = 'gist_gray_r'
    
if(figure_setup == "NationalReport2021"):
    lon_min=10;lat_min=53;lon_max=30.5;lat_max=66;
    start=mp.dates.datetime.datetime(2010,1,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_name = 'NationalReport2021'
    figure_size=(12,9)
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\NationalReport2021\\"
    bathy_colormap = 'gist_gray_r'
    
if(figure_setup == "GotlandD"):
    lon_min=18;lat_min=56;lon_max=21;lat_max=59;
    start=mp.dates.datetime.datetime(2013,3,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_size=(10,15)  #default value!
    figure_name="FMIGotlandDeep"
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\GotlandDeep\\"
    
if(figure_setup == "Bothnian Sea"):
    start=mp.dates.datetime.datetime(2019,3,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_name="FMIBothnianSea"
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BothnianSea\\"
    lon_min=17;lat_min=60;lon_max=22;lat_max=63;
    
if(figure_setup == "Bothnian Sea Aranda"):
    start=mp.dates.datetime.datetime(2019,3,1)
    end=mp.dates.datetime.datetime(2230,5,5)
    figure_name="FMIBothnianSeaAranda"
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BothnianSea\\"
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
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BayOfBothnia\\"
    lon_min=20;lat_min=64;lon_max=26;lat_max=66;
    
if(figure_setup == "EARISE_BP"):
    figure_name="EuroArgoRISE"
    dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\" 
#    lon_min=18.2;lat_min=57.7;lon_max=22.4;lat_max=59.3;
#    lon_min=19.0-4.0;lat_min=58.0-2.0;lon_max=21.0+4.0;lat_max=59.8+2.0;
    lon_min=19.0-1.0;lat_min=58.0-1.0;lon_max=21.0+1.0;lat_max=59.8+1.0;
    plot_contours = True
    replace_labels = {'6903703':'ARVOR-I(6903703)',\
                      '6903704':'APEX(6903704)'}
    figure_size=(10,5)
    line_alpha=0.8
    label_step = 0.2
    bathy_max = 200.0
    marker_end_size = 10
    marker_start_size = 5
    marker_size = 5
    all_colors = ['#ff0000','#ffffff']
    requested_aspect = 1.0




fig=plt.figure(figsize=figure_size)
plt.clf()
#bmap = Basemap(llcrnrlon=lon_min,llcrnrlat=lat_min,urcrnrlon=lon_max,urcrnrlat=lat_max, \
#resolution = 'i',fix_aspect=False)


#ax = plt.axes(projection=the_proj)
ax = plt.axes(projection=requested_proj)
ax.set_extent([lon_min,lon_max,lat_min,lat_max])
ax.set_aspect(requested_aspect)
ax.coastlines(shore_resolution)
#ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', '10m',\
#                                        edgecolor='face', facecolor=(0.6,0.6,0.65)))
#gl = ax.gridlines(crs=the_proj, draw_labels=True,
#          linewidth=2, color='gray', alpha=0.1, linestyle='-')



#ax.coastlines(shore_resolution,zorder=4, alpha = 0.5)
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'land', shore_resolution,\
                                        edgecolor='#000000', linewidth = 0.1,\
                                        facecolor='#ccccdd', alpha = 1.0))
ax.add_feature(cfeature.NaturalEarthFeature('physical', 'lakes', shore_resolution,\
                                        edgecolor='#000000', linewidth=0.5, \
                                        facecolor='#ffffff', alpha = 1.0))
grid_proj = ccrs.PlateCarree()
gl = ax.gridlines(crs=grid_proj, draw_labels=draw_labels,
          linewidth=2, color='gray', alpha=0.1, linestyle='-')
gl.xlabels_top = False
gl.ylabels_right = False
gl.top_labels = False
gl.right_labels = False
#bmap.drawcoastlines(linewidth=0.5)
#bmap.fillcontinents()
#bmap.drawparallels(np.arange(50.,69,label_step),labels=[1,0,0,0],linewidth=0,dashes=[5,10])
#bmap.drawmeridians(np.arange(12.,30,label_step),labels=[0,0,0,1],linewidth=0,dashes=[5,10])



files_to_plot=[i for i in os.listdir(dir_to_plot) if i.endswith('.nc')]
labels= list(map(lambda x: re.search('\d{7}',files_to_plot[x]).group(0),range(len(files_to_plot))))
colors=all_colors[0:len(files_to_plot)]
print(files_to_plot,labels)

#TOPOGRAPHY EXPERIMENT
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
    gl.xlabels_top = False
    gl.ylabels_right = False
   
if(draw_EEZ):
    ax.add_wms('http://geo.vliz.be/geoserver/MarineRegions/wms?',\
               layers='eez_boundaries')
        
color_stack = colors[:]

if make_leaflet:
    map_center = [(lat_min+lat_max)/2.0, (lon_min+lon_max)/2.0]
    map_radius = ah.distance(map_center, [lat_min,lon_min])
    zoom_level = np.log2(40000.0/(map_radius*0.75))
    leaflet_map = folium.Map(location = map_center, zoom_start = zoom_level)
    if plot_bathymetry:
        tmp_topo = topoin.data.copy()
        tmp_topo[topoin.mask] = 0.0
        tmp_topo[tmp_topo>bathy_max] = bathy_max
        tmp_topo = -1.0*tmp_topo/bathy_max
        tmp_topo = cmo.cm.ice(tmp_topo)
        tmp_topo[:,:,3] = ~topoin.mask # transparency
        tmp_topo = tmp_topo[-1::-1,:]
        leaflet_bathy = folium.raster_layers.ImageOverlay(\
                        tmp_topo, [[lats[0],lons[0]],[lats[-1],lons[-1]]]
#                        ,colomap = cmo.cm.topo
                         ,opacity = 0.5
                         ,mercator_project = True
                        )
        leaflet_map.add_child(leaflet_bathy)
        
    
if plot_routes:
    for f,label in zip(files_to_plot,labels):
        if(len(color_stack)==0):
            color_stack = colors[:]
        lab = label
        if(lab in replace_labels.keys()):
            lab = replace_labels[lab]  #some image setups want specific labels
        d=xr.open_dataset(dir_to_plot+f)
        time_var = 'JULD' 
        if('JULD' not in list(d.keys())):
            time_var = 'TIME'        
        primaries = ah.get_primary_indices(d)
        primaries = np.asarray(primaries) & \
                    np.asarray(d[time_var]>np.datetime64(start)) &\
                    np.asarray(d[time_var]<np.datetime64(end))
        lat_dat = np.array(d['LATITUDE'])[primaries]
        lon_dat = np.array(d['LONGITUDE'])[primaries]
        d.close()
#        x,y=bmap(lon_dat,lat_dat)
    #    bmap.plot(x,y,color=col,linewidth=2,alpha=0.5)
        print(len(lon_dat))
        print(label,d[time_var].min(),d[time_var].max() )
        if(len(lon_dat)>0):
            col = color_stack.pop(0)
#            bmap.plot(x,y,color=col,linewidth=line_width, alpha=line_alpha)
#            bmap.plot(x[-1],y[-1],'x',color=col,markersize=marker_end_size,alpha=1.0)
#            bmap.plot(x[0],y[0],'o',color=col,markersize=marker_start_size,alpha=1.0,label=lab)
#            if(plot_points):
#                bmap.plot(x,y,'.',color=col,markersize=marker_size,alpha=line_alpha)
            plt.plot(lon_dat,lat_dat,color=col,linewidth=line_width, alpha=line_alpha,\
                     transform = ccrs.PlateCarree())
            plt.plot(lon_dat[-1],lat_dat[-1],'x',color=col,markersize=marker_end_size,\
                     alpha=1.0, transform = ccrs.PlateCarree())
            plt.plot(lon_dat[0],lat_dat[0],'o',color=col,markersize=marker_start_size,\
                     alpha=1.0,label=lab, transform = ccrs.PlateCarree())
            if(plot_points):
                plt.plot(lon_dat,lat_dat,'.',color=col,markersize=marker_size,\
                         alpha=line_alpha, transform = ccrs.PlateCarree())
            print(lon_dat[-1],lat_dat[-1], lab)
    #    print lab, mp.dates.num2date(a.obs['ape']['date'][0]).date() \
    #             , mp.dates.num2date(a.obs['ape']['date'][-1]).date()
        if make_leaflet:
            locations = list(zip(lat_dat[~np.isnan(lat_dat)],\
                                 lon_dat[~np.isnan(lon_dat)]))
            folium_line = folium.PolyLine(locations = locations, \
                                weight = 2, color = col, opacity = line_alpha)
            leaflet_map.add_child(folium_line)
            leaflet_map.add_child(\
                    folium.CircleMarker(locations[0], marker_start_size, color = col))
            leaflet_map.add_child(\
                    folium.Marker(locations[-1],
                        icon = folium.Icon(color='#009900'),
                        riseOnHover = True,
                        title = lab))

if plot_legends:
    #plt.legend(bbox_to_anchor=(1.0,0.5),numpoints=1)
    plt.legend(loc='lower right',numpoints=1,prop={'size': legend_size})
plt.savefig(output_dir+figure_name+'.png' ,\
            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
print("saved: {}".format(output_dir+figure_name+'.png'))
#plt.savefig(output_dir+figure_name+'.eps' ,\
#            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
if make_leaflet:
    leaflet_map.save(output_dir+figure_name+'.html')
    print("saved: {}{}.html".format(output_dir,figure_name))
