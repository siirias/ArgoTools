# -*- coding: utf-8 -*-
"""
Created on Fri May 22 12:04:54 2020

@author: siirias
"""
import re
import numpy as np
import pandas as pd
import datetime as dt
in_file = "D:\\Data\\ArgoData\\PoijuData9234\\9234.007.log"

whole_file = open(in_file).readlines()
values = 11
valuesh = values-4  # koska osa headereistä on yhdessä solussa
headers_cell = 16
full_data = []
for current_line in whole_file:
    if 'Profile()' in current_line:
        time = dt.datetime.strptime(current_line.split(',')[0],'(%b %d %Y %H:%M:%S')
        line = re.sub('/',' ',current_line)
        line = re.sub(',',' ',line)
        line = re.sub(':','',line)
        dat = line.split()
        headers_tmp = ['Time', 'P','T','S','O2']+dat[headers_cell+2:headers_cell+values-2]
        data = dat[headers_cell+valuesh+2:headers_cell+2+valuesh+values]
        units = ['']+map(lambda x: re.sub('[0-9\.\-\s]*','',x),data)
        headers = []
        for i,j in zip(headers_tmp, units):
            if(j is not ''):
                headers.append("{}({})".format(i,j))
            else:
                headers.append(i)
        data = map(lambda x: float(re.search('-?[0-9\.]*',x).group()),data)
        data = [time] + data
        full_data.insert(0,data)
full_data = np.array(full_data)
full_data = pd.DataFrame(full_data,columns = headers)

#to print the dataframe:
with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', 200):
   print(full_data)