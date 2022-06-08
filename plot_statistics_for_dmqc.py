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

ylimits = [(0.0,120.0), (0.0,230.0), (0.0,120.0), (0.0,230.0)]
fig_size = [8,10]
fig_dpi = 300
add_std = True
add_minmax = True
add_50p = False
for f,ylims in zip(files,ylimits):
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
    name_extra = ''
    plt.figure(figsize = fig_size)
    plt.plot(data['mean'], depths,'k.-')
    if add_minmax:
        plt.fill_betweenx(depths, data['min'],data['max'], color = 'b', alpha = 0.25)
        name_extra+="_MinMax"
    if add_50p:
        plt.fill_betweenx(depths, data['25%'],data['75%'], color = 'b', alpha = 0.5)
        name_extra+="_50p"
    if add_std:
        plt.fill_betweenx(depths, data['mean'] - data['std'], data['mean'] + data['std'], color = 'g', alpha = 0.5)
        name_extra += "STD"
    plt.grid('on')
    plt.ylim(ylims)
    plt.gca().invert_yaxis()
    plt.title("Mean {} in {}".format(variable_name, area))
    plt.ylabel('Pressure({})'.format(unit))
    plt.xlabel("{}({})".format(variable_name, unit))
    filename = "DMQC{}_{}{}".format(variable_name, area,name_extra)
    plt.savefig(output_dir+filename+'.png' ,\
                facecolor='w',dpi=fig_dpi,bbox_inches='tight')


