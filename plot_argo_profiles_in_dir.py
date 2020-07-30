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
from mpl_toolkits.basemap import Basemap
from netCDF4 import Dataset
import xarray as xr
import argohelper as ah

dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\Cape\\"
figure_setup = "EARISE_BP"
figure_name = "ArgoPlot"
fig_dpi = 300
interp_depths = np.array(range(210))

figure_name+="_profile"
plot_profiles = True
variables = ['TEMP','PSAL']
start=mp.dates.datetime.datetime(1000,5,5)
end=mp.dates.datetime.datetime(3030,5,5)


files_to_plot=[i for i in os.listdir(dir_to_plot) if i.endswith('.nc')]
labels= map(lambda x: re.search('\d{7}',\
                    files_to_plot[x]).group(0),range(len(files_to_plot)))                
colors=all_colors[0:len(files_to_plot)]
if plot_profiles:
    for f,col,lab in zip(files_to_plot,colors,labels):
        for var in variables:
            fig=plt.figure(figsize=figure_size)
            plt.clf()
            d=xr.open_dataset(dir_to_plot+f)
            primaries = ah.get_primary_indices(d)
            interp_data = ah.interpolate_data_to_depths(\
                            np.array(d[var])[primaries,:],\
                            np.array(d['PRES'])[primaries,:],\
                            interp_depths)
    
            plt.pcolormesh(\
                    np.array(d['JULD'])[primaries],\
                    interp_depths,\
                    np.transpose(interp_data))
            ah.give_statistics(dir_to_plot+f)
            plt.gca().invert_yaxis()
            plt.colorbar()
            plt.title(var)

plt.savefig(dir_to_plot+figure_name+'.png' ,facecolor='w',dpi=fig_dpi)
plt.savefig(dir_to_plot+figure_name+'.eps' ,facecolor='w',dpi=fig_dpi)
