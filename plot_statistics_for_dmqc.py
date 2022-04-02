# -*- coding: utf-8 -*-
"""
Created on Sat Apr  2 14:17:25 2022

Plotting the statistical analysis for DMQC 
@author: siirias
"""
import numpy as np
import pandas as pd
import re
import matplotlib.pyplot as plt

data_dir="C:\\Data\\EARiseQC\\FTP\\FMI\\" #default value
output_dir = "C:\\Data\\ArgoData\\Figures\\"

files =['ICES_Statistics_Practical_Salinity_dmnless_by_depth_in_BothSea.csv',
        'ICES_Statistics_Practical_Salinity_dmnless_by_depth_in_BP.csv',
        'ICES_Statistics_Temperature_degC_by_depth_in_BothSea.csv', 
        'ICES_Statistics_Temperature_degC_by_depth_in_BP.csv']


fig_size = [8,10]
fig_dpi = 300

for f in files:
    data = pd.read_csv(data_dir+f)
    #split the name and get same variables
    variable_name = f.split('_')[2]
    unit = f.split('_')[3]
    area = f.split('_')[-1].split('.')[0]
    #then some specifics to get them right:
    if variable_name in ['Practical']:
        variable_name = 'Practical Salinity'
        unit = 'dmless'
    if variable_name in ['Temperature']:
        unit = '°C'
    
    #bit of gludge to get the pressure window into a mean value:
    depths = \
        np.array(list(\
        map(lambda y: (float(y[0])+float(y[1]))/2.0,\
        map(lambda x: re.search('([\d\.]+)[^\d]*([\d\.]+)',x)\
        .groups(), data['PRES_window']))))
    
    plt.figure(figsize = fig_size)
    plt.plot(data['mean'], depths,'k.-')
    plt.fill_betweenx(depths, data['min'],data['max'], color = 'b', alpha = 0.25)
    plt.fill_betweenx(depths, data['25%'],data['75%'], color = 'b', alpha = 0.5)
    plt.grid('on')
    plt.gca().invert_yaxis()
    plt.title("Mean {} in {}".format(variable_name, area))
    plt.ylabel('Pressure({})'.format(unit))
    plt.xlabel("{}({})".format(variable_name, unit))
    filename = "DMQC{}_{}".format(variable_name, area)
    plt.savefig(output_dir+filename+'.png' ,\
                facecolor='w',dpi=fig_dpi,bbox_inches='tight')


