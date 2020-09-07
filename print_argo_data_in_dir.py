# -*- coding: utf-8 -*-
"""
Created on Mon Apr 11 18:11:08 2016

@author: siirias
"""
import sys
sys.path.insert(0,'D:\\Data\\svnfmi_merimallit\\qa\\nemo')
import os
import re
import numpy as np
import xarray as xr
import pandas as pd
import argohelper as ah
import datetime as dt

#dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\GotlandDeep\\"
dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\"
#dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\AllFinnish\\"
#dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\Cape\\"
output_dir = "D:\\Data\\ArgoData\\Figures\\"
figure_name = "FMI_Argos"
start=dt.datetime(1000,5,5)
end=dt.datetime(3030,5,5)
print_extras = True
print_extras_to_file = True
statistics = {'wmo':'WMO',
              'type':'Type',
              'area':'Area',
              'time_deployed':'Mission start',
              'time_last_profile':'Last profile',
              'serial':'Serial number',
              'nickname':'Nickname',
              'sensors':'Sensors',
              'deployment_lat': 'Deployment latitude',
              'deployment_lon': 'Deployment longitude'}
stats_shown = ['wmo','serial','nickname','type',
               'area','sensors','time_deployed','time_last_profile',
               'deployment_lat','deployment_lon']

files_to_print=[i for i in os.listdir(dir_to_plot) if i.endswith('.nc')]
file_out = open(output_dir+figure_name+'.csv','w')
if(print_extras_to_file):
    extras_file_out = open(output_dir+figure_name+'extras.txt','w')
#print header
file_out.write("\t".join(stats_shown)+"\n")
    
for f in files_to_print:
    float_name = re.search("[^\d]*(\d*)[^\d]?",f).groups()[0]
    d=xr.open_dataset(dir_to_plot+f)
    primaries = ah.get_primary_indices(d)
    stats = ah.gather_statistics(d,primaries)
    txt = ""
    for i in stats_shown:
        tmp = stats[i]
        if(isinstance(tmp,dt.datetime)):
            tmp=tmp.strftime("%Y-%m-%d")
        if txt == "":
            txt+=tmp
        else:
            txt+="\t"+str(tmp)
    txt += "\n"
    print(txt)
    file_out.write(txt)
    if (print_extras):
        txt=[]
        txt.append("Max distance: {:.1f} km, max jump between profiles: {:.1f} km".\
              format(np.nanmax(stats['distance_from_origin']),\
                     np.nanmax(stats['distance_since_last'])))
        txt.append("Median distance: {:.1f} km".\
              format(np.nanmedian(stats['distance_since_last'])))
        txt.append("Diving depth, max: {:.1f} m, mean: {:.1f} m, min: {:.1f} m".\
              format(stats['depth_max'], stats['depth_avg'], stats['depth_min'], ))
        for i,j,cn, d_tot, d_last in zip(d['JULD'][primaries],\
                          stats['times_between'],\
                          d['CYCLE_NUMBER'][primaries],\
                          stats['distance_from_origin'],\
                          stats['distance_since_last']):
            txt.append(\
             "cycle: {},\t{}\tdiff:{:.1f} hours ({:.1f} days) Travelled:{:.1f}(d:{:.1f})".format(\
                      cn.data,\
                      pd.to_datetime(i.data).strftime("%Y-%m-%d %H.%M"),\
                      j,j/24.0,\
                      d_tot,d_last))
        for t in txt:
            print(t)
            if(print_extras_to_file):
                extras_file_out.write(t+"\n")
            
file_out.close()
if(print_extras_to_file):
    extras_file_out.close()        
