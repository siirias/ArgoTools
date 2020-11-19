# -*- coding: utf-8 -*-
"""
Created on Wed Nov 18 16:49:08 2020

@author: siirias
"""

import re
import os
import pandas as pd
import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import matplotlib.dates as mdates
import math


input_dir = "D:\\Data\\ArgoData\\ArgoRawData\\"

float_sets = [
        {'n':"f9234", 'loc':"f9234\\"},
        {'n':"f9234_a", 'loc':"f9234\\Old_data\\"},
        {'n':"f9234_b", 'loc':"f9234\\Old_data\\2013_2015\\"},
        {'n':"f9089", 'loc':"f9089\\"},
        {'n':"f9089", 'loc':"f9089\\BAPE2_2014\\"},
        {'n':"f9089", 'loc':"f9089\\BAPE2_2016\\"},
        {'n':"f8907", 'loc':"f8907\\"},
        {'n':"f8543", 'loc':"f8543\\"},
        {'n':"f8540", 'loc':"f8540\\"},
        {'n':"f8348", 'loc':"f8348\\"},
#        {'n':"f8348_b", 'loc':"f8348\\iceseason18-19\\"},
        {'n':"f7126", 'loc':"f7126\\"},
        {'n':"f7126_b", 'loc':"f7126\\APE1_2012\\"},
        {'n':"f7126_c", 'loc':"f7126\\APE1_2014-2015\\"},
        {'n':"f7087", 'loc':"f7087\\"},
        {'n':"f7087_b", 'loc':"f7087\\APE2_2013\\"},
        {'n':"f7087_c", 'loc':"f7087\\APE2_2016\\"},
        {'n':"f7087_d", 'loc':"f7087\\APE2_2017\\"}
        
        ]
for f_s in float_sets:
    current_input_dir = input_dir+f_s['loc']
    f_s['type'] = "apex9"
    files_to_handle = [i for i in os.listdir(current_input_dir) if re.match(".*\.\d\d\d\.log$",i)]
    # Newer floats have different file setup.
    # if these files are not found, let's search the others
    if(len(files_to_handle)==0):
        f_s['type'] = "apex11"
        files_to_handle = [i for i in os.listdir(current_input_dir) if re.match(".*vitals_log.csv$",i)]
    print("FLoat: {}, type: {}".format(f_s['n'],f_s['type']))
     
    voltages = []
    cycles = []
    for file in files_to_handle:
        with open(current_input_dir+file,'r') as f:
            lines = f.readlines()
            #Then let's split on which kind of files we handle:
            voltage = None
            if(f_s['type'] == "apex9"):
                for l in lines:
                    if(re.match(".*CtdPower.*Volts",l)):
                        voltage = float(re.search(".*CtdPower.* ([\d\.]+)Volts",l).groups()[0])
            if(f_s['type'] == "apex11"):
                for l in lines:
                    if(re.match("VITALS_CORE.*",l)):
                        voltage = float(re.search(\
                        "VITALS_CORE,[^,]*,[^,]*,[^,]*,([^,]*),"\
                        ,l).groups()[0])
                        break # there are several of these lines, but first is most recent.
            if(voltage):
                voltages.append(voltage)
                cycles.append(int(re.search(".*\.(\d\d\d)\.",file).groups()[0]))
        f_s['voltages'] = voltages.copy()
        f_s['cycles'] = cycles.copy()
plt.figure()        
for i in float_sets:
    if(i['type'] == 'apex11'):
        lw = 2.0
    else:
        lw = 1.0
    if('voltages' in i.keys()):
        plt.plot(i['cycles'][1:],i['voltages'][1:], label = i['n'], linewidth = lw)
    else:
        print("Different format: {}".format(i['n']))
    plt.legend()
    plt.ylabel('Voltage(V)')
    plt.xlabel('Cycle')
    