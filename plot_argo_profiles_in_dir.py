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
import xarray as xr
import pandas as pd
import argohelper as ah
import cmocean as cmo
import plotly.express as px #Only needed for plotly output
import plotly.io as pio     #Only needed for plotly output
import gsw

make_plotly = False
close_figures = False
file_format = "new_server"  # "old_server" "new_server"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\arvorc\\"
#dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\"
dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\ice_examples\\"
#dir_to_plot="C:\\Data\\ArgoData\\BSSC2025\\set2\\coriolis\\2903899\\profiles\\"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\Cape\\"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BGC_BP\\"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\AllFinnish\\"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BarentsSea\\"
output_dir = "C:\\Data\\ArgoData\\Figures\\"
figure_size_timeline =(6,2.5) #(10,4)
figure_size_profile = (3.5,5)#(7,10)

timeline_xtics_rotation = 45.0  # 0.0
timeline_xtics_fontsize = 8.0 #or none

max_depth = 250.0
tl_min = None #4.5 #None #4.5   # fixes axes for each variable, 
tl_max = None #7.0 #None #23.0  # so usualy work for just one at a time.
fig_dpi = 300
c_map = 'viridis'
interp_depths = np.array(np.arange(0,max_depth,0.1))
plot_profile_timelines = True
plot_profile_clusters = True
cluster_grid = True
profile_cloud_alpha = 0.2
enhance_temperature_min = -100.0 # -100.0 would ignore this
#variables = ['TEMP','PSAL_ADJUSTED', 'DENSITY', 'DOX2']
variables = ['TEMP']
#variables = ['TEMP','PSAL','DOX2', 'BBP700', 'CPHL_ADJUSTED', 'CDOM', \
#             'DOWN_IRRADIANCE380', 'DOWN_IRRADIANCE412', 'DOWN_IRRADIANCE490']
#start=mp.dates.datetime.datetime(1000,5,5)
#end=mp.dates.datetime.datetime(3030,5,5)

start=np.datetime64('2000-11-01T00:00:00.0')
end=np.datetime64('2121-12-20T00:00:00.0')
def calculate_density(data, temperature_field = 'TEMP', 
                      salinity_field = 'PSAL',
                      pressure_field = 'PRES'):
    #returns a field that can be added into loaded
    #argo data. 
    temperature = data[temperature_field]
    salinity = data[salinity_field]
    pressure = data[pressure_field]
    density = gsw.rho(salinity, temperature, pressure)
    density.attrs.update({'standard_name':"sea_water_density",
                          'long_name':"sea water density",
                          'units':"kgm-3"})
    return density
    
    
if make_plotly:
    pio.renderers.default='browser'
    out_file = "{}.html".format(re.search("\\\([^\\\]*)\\\$",\
                                        dir_to_plot).groups()[0])
    plotly_file = open(output_dir + out_file,'w')

time_var = 'JULD' 
if(file_format == 'new_server'):
    time_var = 'TIME'

files_to_plot=[i for i in os.listdir(dir_to_plot) if i.endswith('.nc')]
if plot_profile_timelines:
    for f in files_to_plot:
        d=xr.open_dataset(dir_to_plot+f)
        if(not time_var in d.keys()): #gludge for different type of files
            if 'JULD' in d.keys():
                time_var = 'JULD'
            elif 'TIME' in d.keys():
                time_var = 'TIME'
            else:
                print("No suitable time variable!")
#        print("Availabe variables: {}".format(list(d.keys())))
        f_latitude = float(d['LATITUDE'][0])
        f_longitude = float(d['LONGITUDE'][0])
        f_area = ah.give_area(f_latitude, f_longitude)
        for var in variables:
            if var == 'DENSITY': #special case, let's calculate
                d[var] = calculate_density(d)
            if(var in d.keys()):
                variable_print_name = d[var].attrs['long_name']

                variable_print_unit = d[var].attrs['units']

                if(variable_print_name == 'SEA TEMPERATURE IN SITU ITS-90 SCALE'):
                    variable_print_name = 'Temperature'
                if(variable_print_unit == 'degree_Celsius'):
                    variable_print_unit = '°C'

                fig=plt.figure(figsize=figure_size_timeline)
                float_name = re.search("(\d{7})",f).groups()[0]
                plt.clf()
                cmap = ah.axes_label_from_variable_name(var, give_colormap=True)[1]            
                primaries = ah.get_primary_indices(d)
                primaries = np.asarray(primaries) & \
                            np.asarray(d[time_var]>start) &\
                            np.asarray(d[time_var]<end)
                if(np.sum(primaries)>0):
                    print("PLOTTING!")
                    if(var in ['DOXY', 'DOX2', 'BBT700', 'CPHL_ADJUSTED']): #these must be taken from secondary profiles
                        primaries = []
                        for i in d[var]:
                            if np.isnan(np.max(i)):
                                primaries.append(True)
                            else:
                                primaries.append(False)
                        primaries = list(map(lambda x: not x,primaries))
                        
                        
                    interp_data = ah.interpolate_data_to_depths(\
                                    np.array(d[var])[primaries,:],\
                                    np.array(d['PRES'])[primaries,:],\
                                    interp_depths)
            
                    plt.pcolormesh(\
                            np.array(d[time_var])[primaries],\
                            interp_depths,\
                            np.transpose(interp_data),\
                            cmap = cmap,\
                            shading = 'auto',
                            vmin = tl_min,
                            vmax = tl_max)
                    plt.gca().invert_yaxis()
                    cbar = plt.colorbar(pad = 0.01, fraction = 0.05)
    #                cbar.set_label(ah.axes_label_from_variable_name(var))
                    cbar.set_label("{}/{}".format(variable_print_name,
                                              variable_print_unit))
                    if(var == 'TEMP' and enhance_temperature_min>-50.0):
                        tmp_data = interp_data.copy()
                        tmp_data[tmp_data>enhance_temperature_min] = np.nan
                        # plt.pcolormesh(\
                        #         np.array(d[time_var])[primaries],\
                        #         interp_depths,\
                        #         np.transpose(tmp_data),\
                        #         cmap = cmo.cm.gray,\
                        #         shading = 'auto',
                        #         zorder = 10)
                        plt.axhline(35,color = 'b',zorder = 11)
                        plt.axhline(18,color = 'r',zorder = 11)
                        # cbar = plt.colorbar(pad = 0.01, fraction = 0.05)
                    plt.title(f"Float {float_name} ({f_area})\n")
                    plt.ylabel(ah.axes_label_from_variable_name('PRES'))
                    plt.xlabel(ah.axes_label_from_variable_name(time_var))
                    plt.xticks(rotation=timeline_xtics_rotation)
                    xticks_labels = plt.gca().get_xticklabels()
                    if(timeline_xtics_fontsize):
                        # Set the font size for the tick labels
                        for tick in xticks_labels:
                            tick.set_fontsize(timeline_xtics_fontsize)
                    filename = "{}_{}_tl".format(float_name,var)
                    plt.savefig(output_dir+filename+'.png' ,\
                                facecolor='w',dpi=fig_dpi,bbox_inches='tight')
                    # plt.savefig(output_dir+filename+'.eps' ,\
                    #             facecolor='w',dpi=fig_dpi,bbox_inches='tight')
                    print("Saved: {}".format(output_dir+filename+'.png'))
                    if(close_figures):
                        plt.close()
                else:
                    print('nothing to plot')
if plot_profile_clusters:
    for f in files_to_plot:
        for var in variables:
            d=xr.open_dataset(dir_to_plot+f)
            if(not time_var in d.keys()): #gludge for different type of files
                if 'JULD' in d.keys():
                    time_var = 'JULD'
                elif 'TIME' in d.keys():
                    time_var = 'TIME'
                else:
                    print("No suitable time variable!")
            f_latitude = float(d['LATITUDE'][0])
            f_longitude = float(d['LONGITUDE'][0])
            f_area = ah.give_area(f_latitude, f_longitude)
            
            if var == 'DENSITY': #special case, let's calculate
                d[var] = calculate_density(d)
            if(var in d.keys()):
                variable_print_name = d[var].attrs['long_name']
                variable_print_unit = d[var].attrs['units']
                float_name = re.search("(\d{7})",f).groups()[0]
                fig=plt.figure(figsize=figure_size_profile)
                plt.clf()
                primaries = ah.get_primary_indices(d)
                # primaries = np.asarray(primaries) & \
                #             np.asarray(d[time_var]>np.datetime64(start)) &\
                #             np.asarray(d[time_var]<np.datetime64(end))
                d_sel = d[var][primaries]
                d_sel_time = d[time_var][primaries]
                d_sel_pres = d['PRES'][primaries]
                #create colorbar for time-indices
                sm = plt.cm.ScalarMappable(cmap = c_map, \
                            norm=plt.Normalize(vmin = d[time_var].min(),\
                                               vmax = d[time_var].max()))
                sm._A=[] # this is needed as scalar mappable needs someting to map.
                clims = sm.get_clim()
                for i in range(d_sel.shape[0]):
                    color = float(d_sel_time[i]-clims[0])/float(clims[1]-clims[0])
                    color = sm.get_cmap()(color)
                    plt.plot(d_sel[i,:], d_sel_pres[i,:],\
                             color = color, alpha = profile_cloud_alpha)
                if tl_min and tl_max:
                    plt.gca().set_xlim(tl_min, tl_max)
                plt.title(f"Float {float_name} ({f_area})\n")
                plt.ylabel(ah.axes_label_from_variable_name('PRES'))
#                plt.xlabel(ah.axes_label_from_variable_name(var))
                plt.xlabel("{}/{}".format(variable_print_name,
                                          variable_print_unit))
                plt.gca().invert_yaxis()
                cbar = plt.colorbar(sm)            
                cbar.ax.set_yticklabels(\
                        pd.to_datetime(\
                            cbar.get_ticks()).strftime(date_format='%d %b %Y'))
                if cluster_grid:
                    plt.grid()
                filename = "{}_{}_cl".format(float_name,var)
                plt.savefig(output_dir+filename+'.png' ,\
                            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
                # plt.savefig(output_dir+filename+'.eps' ,\
                #             facecolor='w',dpi=fig_dpi,bbox_inches='tight')
                print("Saved: {}".format(output_dir+filename+'.png'))
                if(close_figures):
                    plt.close()
                if make_plotly:
                    pfig = px.line()
if make_plotly:
    plotly_file.close()