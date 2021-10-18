# -*- coding: utf-8 -*-
"""
Created on Tue Oct 12 12:03:54 2021

@author: siirias
"""
import xarray as xr
import pandas as pd
import cmocean as cmo
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt

in_dir = "C:\\Data\\ArgoData\\ForComparison\\"
out_dir = "C:\\Data\\ArgoData\\ForComparison\\output\\"
files = ["GL_PR_PF_6903703.nc","GL_PR_PF_6903704.nc"]
#files = ["GL_PR_PF_6903704.nc"]
#files = ["6903704_20211012093549117.nc","6903703_20211012093527695.nc"]
min_time = dt.datetime(2020,2,1)
max_time = dt.datetime(2021,9,30)
time_name = 'time'
lat_name = 'lat'
lon_name = 'lon'
pres_name = 'pressure (db)'
var_name = {'TEMP':'Temperature (C)',
            'PSAL':'salinity (g/kg)',
            'DOX2_ADJUSTED':'Oxygen (µmol kg-1)'}
qc_name = 'QC'
var_qc = {'TEMP':'TEMP_QC',
            'PSAL':'PSAL_QC',
            'DOX2_ADJUSTED':'DOX2_ADJUSTED_QC'}

"""
QC flags:
0 = no_qc_performed
1 = good_data
2 = probably_good_data
3 = bad_data_that_are_potentially_correctable
4 = bad_data
5 = value_changed
6 = value_below_detection
7 = nominal_value
8 = interpolated_value
9 = missing_value'
"""
variables = ['TEMP', 'PSAL', 'DOX2_ADJUSTED']
f = files[0]
for f in files:
    data = xr.open_dataset(in_dir+f)
    #check the dataset type:
    d_type = "CTD"
    if "DOX2_ADJUSTED" in data.keys():
        d_type = "CTD+O"
    
    wmo = data.platform_code
    
    for var in variables:
        if( var in data.keys()):
            profile_no = 0
            for i in range(len(data.TIME)):
                collected_data = pd.DataFrame(columns = \
                        [time_name, lat_name, lon_name, pres_name, var_name[var], qc_name])
                the_time = pd.to_datetime(data.TIME[i].to_pandas())
                depths = data.PRES[i]
                variable = data[var][i]
                length = len(np.array(variable[~np.isnan(variable)]))
            
            
                profile_pass = True
                if(the_time > max_time):
                    print("too late:{}".format(the_time))
                    Profile_pass = False
                if(the_time < min_time):
                    print("too early:{}".format(the_time))
                    Profile_pass = False
                if(d_type == "CTD+O"):
                    if(var in ["TEMP", "PSAL"] and length < 40):
                        profile_pass = False
                    if(var in ["DOX2_ADJUSTED"] and length < 1):
                        profile_pass = False
                if(d_type == "CTD"):
                    if(var in ["TEMP", "PSAL"] and length < 40):
                        profile_pass = False
                if(profile_pass):
                    profile_no += 1
                    print("{:05} \t{:03.01f} \t{:03.02f} \t{}".format(\
                            profile_no,float(depths[0]),\
                            float(variable[0]), length))
                    for j in range(len(data.DEPTH)):
                        qc_result = data[var_qc[var]][i][j]
                        if(not np.isnan(qc_result)):
                            qc_result = int(qc_result)
                        tmp = {time_name:data.TIME[i].to_pandas(),
                               lat_name:float(data.LATITUDE[i]),
                               lon_name:float(data.LONGITUDE[i]),
                               pres_name:float(data.PRES[i][j]),
                               var_name[var]:float(data[var][i][j]),
                               qc_name:qc_result
                               }
                        if(not np.isnan(tmp[var_name[var]])):
                            collected_data = \
                                collected_data.append(tmp, ignore_index = True)
                    out_file_name = "{}_{:03}_{}.csv".format(wmo, profile_no,var)
                    collected_data.to_csv(out_dir + out_file_name, float_format ='%0.5f')
                    print("saved: {}".format(out_dir + out_file_name))
