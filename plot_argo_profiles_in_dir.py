# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:11:08 2016

@author: siirias
"""
import sys
sys.path.insert(0,'D:\\Data\\svnfmi_merimallit\\qa\\nemo')
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

make_plotly = True
close_figures = False
file_format = "new_server"  # "old_server" "new_server"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\arvorc\\"
#dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\"
dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\RBR_BalticProper\\"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\Cape\\"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BGC_BP\\"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\RBR\\"
#dir_to_plot="C:\\Data\\ArgoData\\ArgosForPlot\\BarentsSea\\"
output_dir = "C:\\Data\\ArgoData\\Figures\\"
figure_size_timeline = (10,4)
figure_size_profile = (7,10)
max_depth = 250.0
tl_min = None #4.5   # fixes axes for each variable, 
tl_max = None #23.0  # so usualy work for just one at a time.
fig_dpi = 300
c_map = 'viridis'
interp_depths = np.array(np.arange(0,max_depth,0.1))
plot_profile_timelines = True
plot_profile_clusters = True
cluster_grid = True
profile_cloud_alpha = 0.2
enhance_temperature_min = -100.0 # -100.0 would ignore this
variables = ['TEMP','PSAL']
#variables = ['TEMP','PSAL','DOX2', 'BBP700', 'CPHL_ADJUSTED', 'CDOM', \
#             'DOWN_IRRADIANCE380', 'DOWN_IRRADIANCE412', 'DOWN_IRRADIANCE490']
#start=mp.dates.datetime.datetime(1000,5,5)
#end=mp.dates.datetime.datetime(3030,5,5)

start=mp.dates.datetime.datetime(1021,6,29)
end=mp.dates.datetime.datetime(3021,12,20)

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
        print("Availabe variables: {}".format(list(d.keys())))
        for var in variables:
            if(var in d.keys()):
                variable_print_name = d[var].attrs['long_name']
                variable_print_unit = d[var].attrs['units']
                fig=plt.figure(figsize=figure_size_timeline)
                float_name = re.search("(\d{7})",f).groups()[0]
                plt.clf()
                cmap = ah.axes_label_from_variable_name(var, give_colormap=True)[1]            
                primaries = ah.get_primary_indices(d)
                # primaries = np.asarray(primaries) & \
                #             np.asarray(d[time_var]>np.datetime64(start)) &\
                #             np.asarray(d[time_var]<np.datetime64(end))
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
                plt.title("Float "+float_name)
                plt.ylabel(ah.axes_label_from_variable_name('PRES'))
                plt.xlabel(ah.axes_label_from_variable_name(time_var))
                filename = "{}_{}_tl".format(float_name,var)
                plt.savefig(output_dir+filename+'.png' ,\
                            facecolor='w',dpi=fig_dpi,bbox_inches='tight')
                # plt.savefig(output_dir+filename+'.eps' ,\
                #             facecolor='w',dpi=fig_dpi,bbox_inches='tight')
                print("Saved: {}".format(output_dir+filename+'.png'))
                if(close_figures):
                    plt.close()
if plot_profile_clusters:
    for f in files_to_plot:
        for var in variables:
            d=xr.open_dataset(dir_to_plot+f)
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
                            norm=plt.Normalize(vmin = d_sel[time_var].min(),\
                                               vmax=d_sel[time_var].max()))
                sm._A=[] # this is needed as scalar mappable needs someting to map.
                clims = sm.get_clim()
                for i in range(d_sel.shape[0]):
                    color = float(d_sel_time[i]-clims[0])/float(clims[1]-clims[0])
                    color = sm.get_cmap()(color)
                    plt.plot(d_sel[i,:], d_sel_pres[i,:],\
                             color = color, alpha = profile_cloud_alpha)
                plt.title("Float " + float_name)
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