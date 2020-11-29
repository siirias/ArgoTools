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
plot_park_pressure_deviations = False
profiles = {}
float_name = ""
stuck_time = 5.0 # hours from which being stuck is calculated
stuck_limit = 0.05  # dbar in stuck_time hours 
def get_pressure_deviations(dat,remove_too_short = True):
        #calculate pressure difference per minute:
        time_d = np.array(\
            map(lambda x: x.total_seconds(),\
            dat['time'].diff()))
        press_d =  dat['pressure'].diff()/time_d
        if(remove_too_short):
            press_d[time_d<60.0] = np.nan
        press_d*=60.0  # to get number in minutes
        return press_d

def get_bottom_contacts(dat):
        windowsize = dt.timedelta(hours=stuck_time/2.0)
        bottom_contacts = [False]*dat['time'].shape[0]
        for i in range(dat['time'].shape[0]):
            t_window = abs(dat['time']-dat['time'][i])<windowsize
            max_diff = abs(dat['pressure'][t_window].min() - \
                       dat['pressure'][t_window].max())
            if(max_diff)<stuck_limit:
                bottom_contacts[i] = True
        return bottom_contacts
    
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

   #PLOTTTING STARTS 
    
if( plot_park_pressures):
    plt.figure(figsize=(10,5))
    for prof in profiles:
        dat = profiles[prof]['park_pressures']
        tgt_pres = profiles[prof]['DeepDescentPressure']
        if(len(dat)>0):
            contacts = get_bottom_contacts(dat)
            plt.plot(dat['time'], dat['pressure'])
            plt.plot([dat['time'].min(),dat['time'].max()],[tgt_pres,tgt_pres],'k')
            plt.fill_between(dat['time'].values, dat['pressure'],\
                             [tgt_pres]*len(dat['time'].values),
                             alpha = 0.2)
            plt.plot(dat['time'][contacts],\
                     dat['pressure'][contacts],\
                        'r.',alpha = 0.05)
    plt.gca().invert_yaxis()
    plt.title(float_name)
    plt.savefig(output_dir+figure_name+'.png' ,facecolor='w',dpi=fig_dpi)
    plt.savefig(output_dir+figure_name+'.eps' ,facecolor='w',dpi=fig_dpi)

if( plot_park_pressure_deviations):
    plt.figure(figsize=(10,5))
    for prof in profiles:
        dat = profiles[prof]['park_pressures']
        if(len(dat)>0):
            #calculate pressure difference per minute:
            press_d = get_pressure_deviations(dat, remove_too_short=False)      
            contacts = get_bottom_contacts(dat)
            plt.plot(dat['time'], press_d)
            plt.plot(dat['time'][contacts], \
                     press_d[contacts],\
                     'r.')
            
    plt.gca().invert_yaxis()
    plt.title(float_name)
