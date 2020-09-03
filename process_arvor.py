# -*- coding: utf-8 -*-
"""
Created on Mon Aug 17 16:49:04 2020

@author: siirias
"""

import os
import shutil
import re
import datetime as dt
process_dir = "D:\\ArvorGDrive\\ToProcess\\"
processed_dir = "D:\\ArvorGDrive\\Archived\\"
extract_dir = "D:\\ArvorGDrive\\Processed\\"
exec_dir = "D:\Data\ArgoData\ArvorSoftware\\"
tmp_dir = "D:\\ArvorGDrive\\tmp\\"
exec_name = "nkeInstrumentationParser.exe"
convert_style = " ARVOR_PROVOR_5900A02_and_higher "
#convert_style = " ARVOR_C_5603L12_and_higher "

files_to_process = os.listdir(process_dir)
files_to_process = [x for x in files_to_process if x.endswith('.sbd')]
# we want to process each file separately to avoid automtical joinings
# and since the nke program only accepts directories, those have to be handled
# in their own directory
process_stamp = dt.datetime.now().strftime("%Y%m%d_%H%M")  # date and time 
for f in files_to_process:
    shutil.move(os.path.join(process_dir,f), \
                os.path.join(tmp_dir,f))
    float_name = re.search("Unit_\s(\d*)[^d]",f).groups()[0]
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
        
for f in files_to_process:
    shutil.move(os.path.join(tmp_dir,f),\
                os.path.join(processed_dir,f))

