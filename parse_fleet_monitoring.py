# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 17:20:48 2022

@author: siirias
Meant for parce out webpage data from Argo Fleet Monitor,
which is copy&pasted here.

There is actually a way to download that as csv, so this would became 
just a simple way to read it in pandas and plot some statistics.
THat seems to have some issues atm tho.

"""


import pandas as pd
import datetime as dt
import re

in_dir = r'C:\Users\siirias\Downloads\\'
in_file = 'argo_fleet.csv'
# in_file = 'dashboard_20-04-2022_5-21-13.csv'
# def afm_to_datetime(in_string):
#     return dt.datetime.strptime(in_string, '%Y-%d-%m %H:%M:%S')

# data = pd.read_csv(in_dir+in_file,\
#                    delimiter=';',\
#                    parse_dates=['Last Tx','Launch Date'],
#                    date_parser=afm_to_datetime)
    
    
raw_dat = open(in_dir+in_file,'r').readlines()

mode ='done'
floats = []
curr_float = {}
for l in raw_dat:
    l = l.strip()
    if mode == 'done':
        if(len(l)>6 and bool(re.match('^\d+$',l))):
            curr_float['wmo'] = l
            mode = 'end_date'
    elif mode == 'end_date':
        if(bool(re.match('^\d\d/\d\d/\d\d\d\d',l))):
            curr_float['end_date'] = dt.datetime.strptime(l,'%d/%m/%Y')
            mode = 'start_date'
    elif mode == 'start_date':
        if(bool(re.match('^\d\d/\d\d/\d\d\d\d',l))):
            curr_float['start_date'] = dt.datetime.strptime(l,'%d/%m/%Y')
            curr_float['age'] = curr_float['end_date'] -curr_float['start_date']
            floats.append(curr_float)
            curr_float = {}
            mode = 'done'
            
for i in floats:
    print("{}\t{}\t{}".format(i['wmo'],i['start_date'],i['age']))
