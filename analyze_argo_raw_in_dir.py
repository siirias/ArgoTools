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
import datetime as dt
import pandas as pd

dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\rawdata\\"
output_dir = "D:\\Data\\ArgoData\\Figures\\"
figure_name = 'park_pressures'
fig_dpi = 300
start=mp.dates.datetime.datetime(1000,5,5)
end=mp.dates.datetime.datetime(3030,5,5)
plot_park_pressures = True
profiles = {}
float_name = ""
#first things from system log
files_to_plot=[i for i in os.listdir(dir_to_plot) if i.endswith('system_log.txt')]
for f in files_to_plot:
    txt = open(dir_to_plot+f).readlines()
    dat ={}
    prof_num = int(re.search('[^\.]*\.([0-9]*)\.',f).groups()[0])
    for l in txt:
        # find the target pressure for this profile
        if('DeepDescentPressure' in l):
            dat['DeepDescentPressure'] = float(\
            re.search('DeepDescentPressure\s*([0-9\.]*)',l).groups()[0])
            print(f,dat['DeepDescentPressure'])
    profiles[prof_num] = dat

#next things from science logs
files_to_plot=[i for i in os.listdir(dir_to_plot) if i.endswith('science_log.csv')]
float_name = re.search('([^\.]*)\.',files_to_plot[0]).groups()[0]
for f in files_to_plot:
    current_mission = "none"
    txt = open(dir_to_plot+f).readlines()
    prof_num = int(re.search('[^\.]*\.([0-9]*)\.',f).groups()[0])
    dat = profiles[prof_num] # pointer to the created object
    print(f)
    park_pressure = []
    park_times = []
    for l in txt:
        # switch to current mission
        if('Mission' in l):
            current_mission = (re.search(',([^,]* Mission)',l).groups()[0])
        if('Park Mission' == current_mission):
            if('CTD_P,' in l):
                tmp = re.search('CTD_P,([0-9T]*),([-+.0-9]*)',l).groups()
                time_stamp = dt.datetime.strptime(tmp[0],'%Y%m%dT%H%M%S')
                pressure = float(tmp[1])
                park_pressure.append(pressure)
                park_times.append(time_stamp)
            
    dat['park_pressures']=pd.DataFrame({'pressure':park_pressure,\
                                     'time':park_times})

if( plot_park_pressures):
    plt.figure(figsize=(10,5))
    for prof in profiles:
        dat = profiles[prof]['park_pressures']
        tgt_pres = profiles[prof]['DeepDescentPressure']
        if(len(dat)>0):
            plt.plot(dat['time'], dat['pressure'])
            plt.plot([dat['time'].min(),dat['time'].max()],[tgt_pres,tgt_pres],'k')
            plt.fill_between(dat['time'].values, dat['pressure'],\
                             [tgt_pres]*len(dat['time'].values),
                             alpha = 0.2)
    plt.gca().invert_yaxis()
    plt.title(float_name)
    plt.savefig(output_dir+figure_name+'.png' ,facecolor='w',dpi=fig_dpi)
    plt.savefig(output_dir+figure_name+'.eps' ,facecolor='w',dpi=fig_dpi)
