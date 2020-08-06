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


dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\"
#dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\Cape\\"
output_dir = "D:\\Data\\ArgoData\\Figures\\"
figure_setup = "EARISE_BP"
figure_name = "ArgoPlot_profile"
figure_size = (8,10)
fig_dpi = 300
c_map = 'viridis'
interp_depths = np.array(np.arange(0,210,0.1))
plot_profile_timelines = True
plot_profile_clusters = True

variables = ['TEMP','PSAL']
start=mp.dates.datetime.datetime(1000,5,5)
end=mp.dates.datetime.datetime(3030,5,5)


files_to_plot=[i for i in os.listdir(dir_to_plot) if i.endswith('.nc')]
if plot_profile_timelines:
    for f in files_to_plot:
        for var in variables:
            fig=plt.figure(figsize=figure_size)
            float_name = re.search("(\d*)_",f).groups()[0]
            plt.clf()
            cmap = ah.axes_label_from_variable_name(var, give_colormap=True)[1]            
            d=xr.open_dataset(dir_to_plot+f)
            primaries = ah.get_primary_indices(d)
            interp_data = ah.interpolate_data_to_depths(\
                            np.array(d[var])[primaries,:],\
                            np.array(d['PRES'])[primaries,:],\
                            interp_depths)
    
            plt.pcolormesh(\
                    np.array(d['JULD'])[primaries],\
                    interp_depths,\
                    np.transpose(interp_data),\
                    cmap = cmap,\
                    shading = 'auto')
            plt.gca().invert_yaxis()
            cbar = plt.colorbar()
            cbar.set_label(ah.axes_label_from_variable_name(var))
            plt.title("Float "+float_name)
            plt.ylabel(ah.axes_label_from_variable_name('PRES'))
            plt.xlabel(ah.axes_label_from_variable_name('JULD'))
            filename = "{}_{}_tl".format(float_name,var)
            plt.savefig(output_dir+filename+'.png' ,\
                        facecolor='w',dpi=fig_dpi,bbox_inches='tight')
            plt.savefig(output_dir+filename+'.eps' ,\
                        facecolor='w',dpi=fig_dpi,bbox_inches='tight')

if plot_profile_clusters:
    for f in files_to_plot:
        for var in variables:
            float_name = re.search("(\d*)_",f).groups()[0]
            fig=plt.figure(figsize=figure_size)
            plt.clf()
            d=xr.open_dataset(dir_to_plot+f)
            primaries = ah.get_primary_indices(d)
            #create colorbar for time-indices
            sm = plt.cm.ScalarMappable(cmap = c_map, \
                        norm=plt.Normalize(vmin = d['JULD'].min(),\
                                           vmax=d['JULD'].max()))
            sm._A=[] # this is needed as scalar mappable needs someting to map.
            clims = sm.get_clim()
            for i in range(d[var].shape[0]):
                color = float(d['JULD'][i]-clims[0])/float(clims[1]-clims[0])
                color = sm.get_cmap()(color)
                plt.plot(d[var][i,:], d['PRES'][i,:], color = color)
            plt.title("Float " + float_name)
            plt.ylabel(ah.axes_label_from_variable_name('PRES'))
            plt.xlabel(ah.axes_label_from_variable_name(var))
            plt.gca().invert_yaxis()
            cbar = plt.colorbar(sm)            
            cbar.ax.set_yticklabels(\
                    pd.to_datetime(\
                        cbar.get_ticks()).strftime(date_format='%d %b %Y'))
            filename = "{}_{}_cl".format(float_name,var)
            plt.savefig(output_dir+filename+'.png' ,\
                        facecolor='w',dpi=fig_dpi,bbox_inches='tight')
            plt.savefig(output_dir+filename+'.eps' ,\
                        facecolor='w',dpi=fig_dpi,bbox_inches='tight')
