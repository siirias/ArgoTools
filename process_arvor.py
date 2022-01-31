# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 16:49:04 2020

@author: siirias
"""

import os
import shutil
import re
import datetime as dt
process_dir = "C:\\ArvorGDrive\\ToProcess\\"
processed_dir = "C:\\ArvorGDrive\\Archived\\"
extract_dir = "C:\\ArvorGDrive\\Processed\\"
exec_dir = "C:\Data\ArgoData\ArvorSoftware\\"
tmp_dir = "C:\\ArvorGDrive\\tmp\\"
exec_name = "nkeInstrumentationParser.exe"
convert_style = " ARVOR_PROVOR_5900A02_and_higher "
#convert_style = " ARVOR_C_5603L12_and_higher "

def parse_location(in_text):
    if(type(in_text) == str):
        lines = in_text.split('\n')
    else:
        lines = in_text
    lat_d = None
    lat_m = None
    lat_f = None
    lon_d = None
    lon_m = None
    lon_f = None
    NS = 'N'
    EW = 'E'
    f_hour = None
    f_minute = None
    f_day = None
    f_month = None
    f_year = None
    batteryt_voltage = None
    for l in lines:
        if('"Batteries voltage at Pmax (V)"' in l):
            batteryt_voltage = float(re.search(";([\d\.]+)",l).groups()[0])
        if('"Float time : Hour"' in l):
            f_hour = int(re.search(";(\d+)",l).groups()[0])
        if('"Float time : Minute"' in l):
            f_minute = int(re.search(";(\d+)",l).groups()[0])
        if('"Float time : Day' in l):
            f_day = int(re.search(";(\d+)",l).groups()[0])
        if('"Float time : Month"' in l):
            f_month = int(re.search(";(\d+)",l).groups()[0])
        if('"Float time : Year"' in l):
            f_year = int(re.search(";(\d+)",l).groups()[0])
        
        if('GPS latitude (°)' in l):
            lat_d = int(re.search(";(\d+)",l).groups()[0])
        if('GPS latitude (minutes)' in l):
            lat_m = int(re.search(";(\d+)",l).groups()[0])
        if('GPS latitude (minutes fractions' in l):
            lat_f = int(re.search(";(\d+)",l).groups()[0])
        if('GPS longitude (°)' in l):
            lon_d = int(re.search(";(\d+)",l).groups()[0])
        if('GPS longitude (minutes)' in l):
            lon_m = int(re.search(";(\d+)",l).groups()[0])
        if('GPS longitude (minutes fractions' in l):
            lon_f = int(re.search(";(\d+)",l).groups()[0])
        if('GPS latitude orientation (0=North 1=South)' in l):
            if(int(re.search(";(\d+)",l).groups()[0]) == 0):
                NS = 'N'
            else:
                NS = 'S'
        if('GPS longitude orientation (0=East 1=West)' in l):
            if(int(re.search(";(\d+)",l).groups()[0]) == 0):
                EW = 'E'
            else:
                EW = 'W'
    if(f_year<1000):
        f_year+=2000  #should work for newxt few hundred years.
    print("Location: {}° {:06.3f}' {} , {}° {:06.3f}' {} @ {:02}.{:02} {:02}.{:02}.{} (UTC)".format(\
            lat_d, float("{}.{}".format(lat_m, lat_f)), NS,\
            lon_d, float("{}.{}".format(lon_m, lon_f)), EW,\
            f_hour, f_minute, f_day, f_month, f_year))
    print("Battery: {} V".format(batteryt_voltage))
    lat = float(lat_d) + float("{}.{}".format(lat_m,lat_f))/60.0
    lon = float(lon_d) + float("{}.{}".format(lon_m,lon_f))/60.0
    if(NS == 'S'):
        lat *= -1.0
    if(EW == 'W'):
        lon *= -1.0
    return (lat,lon)

files_to_process = os.listdir(process_dir)
files_to_process = [x for x in files_to_process if x.endswith('.sbd')]
# we want to process each file separately to avoid automtical joinings
# and since the nke program only accepts directories, those have to be handled
# in their own directory

#find how many floats in the set
float_names = {}
for f in files_to_process:
#    float_names[re.search("Unit_\s(\d*)[^d]",f).groups()[0]] = True
    float_names[re.search("(\d*)_\d*\.sbd",f).groups()[0]] = True
print(float_names)
    

process_stamp = dt.datetime.now().strftime("%Y%m%d_%H%M")  # date and time 

for current_float in float_names:
    for f in files_to_process:
        if(current_float in f): #matches the current float name
            shutil.move(os.path.join(process_dir,f), \
                        os.path.join(tmp_dir,f))
#            float_name = re.search("Unit_\s(\d*)[^d]",f).groups()[0]
            float_name = re.search("(\d*)_\d*\.sbd",f).groups()[0]
            series_number = re.search("_(\d*)\.sbd",f)
            if(type(series_number) == re.Match):
                series_number = series_number.groups()[0]
            else:
                series_number = 'xxx'
    command = os.path.join(exec_dir,exec_name)
    command += convert_style
    command += tmp_dir
    os.system(command)
    #shutil.move(os.path.join(tmp_dir,f), \
    #            os.path.join(process_dir,f))
    created_files = [x for x in os.listdir(tmp_dir) if x.endswith('.csv')]
    #rename_created_files
    for c_f in created_files:
        new_name = float_name+"_"+process_stamp+"_"+c_f
        shutil.move(os.path.join(tmp_dir,c_f), \
                    os.path.join(extract_dir,new_name))
        if('Technical Message 1' in new_name):
            print(parse_location(\
                    open(os.path.join(extract_dir,new_name),'r')\
                        .readlines()))
            
    for f in files_to_process:
        shutil.move(os.path.join(tmp_dir,f),\
                    os.path.join(processed_dir,f))

