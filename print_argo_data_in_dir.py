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

dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\EARise_BP\\"
#dir_to_plot="D:\\Data\\ArgoData\\ArgosForPlot\\Cape\\"
output_dir = "D:\\Data\\ArgoData\\Figures\\"
figure_setup = "EARISE_BP"
figure_name = "ArgoPlot_profile"
variables = ['TEMP']
start=dt.datetime(1000,5,5)
end=dt.datetime(3030,5,5)


files_to_print=[i for i in os.listdir(dir_to_plot) if i.endswith('.nc')]
file_out = open(output_dir+figure_name+'.csv','w')
file_out.write("WMO\ttype\tArea\tMission start\tLast profile\n")
for f in files_to_print:
    for var in variables:
        float_name = re.search("[^\d]*(\d*)[^\d]?",f).groups()[0]
        d=xr.open_dataset(dir_to_plot+f)
        primaries = ah.get_primary_indices(d)
        stats = ah.gather_statistics(d,primaries)
        txt = "{}\t{}\t{}\t{}\t{}\n".format(\
              str(stats['wmo']),\
              str(stats['type']),\
              stats['area'],\
              stats['time_deployed'].strftime("%Y-%m-%d"), \
              stats['time_last_profile'].strftime("%Y-%m-%d"))
        print(txt)
        file_out.write(txt)
        for i,j,cn in zip(d['JULD'][primaries],\
                          stats['times_between'],\
                          d['CYCLE_NUMBER'][primaries]):
            print("{},\t{}\tdiff:{:.1f} hours ({:.1f} days)".format(\
                      cn.data,\
                      pd.to_datetime(i.data).strftime("%Y-%m-%d %H.%M"),\
                      j,j/24.0))
            
file_out.close()
        
