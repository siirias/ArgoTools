# -*- coding: utf-8 -*-
"""
Created on Mon Jun 13 11:55:51 2022

@author: siirias
"""

import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import xarray as xr
import pandas as pd

import arandapy as apy  # this is from: https://github.com/siirias/ArandaTools
import argohelper as ah
import re


# indir1 = "C:/Data/ArandaVEMIVE2020/"
# indir2 = "C:/Data/ArgoData/ArgosForPlot/NBalticProper/"

# infile1 = "a200257a.cnv"
# infile2 = "GL_PR_PF_6903704.nc"
# profile_no2 = [129] # [226,230] for Euro-Argo cruise

# indir1 = "C:/Data/ArandaVEMIVE2020/"
# indir2 = "C:/Data/ArgoData/ArgosForPlot/NBalticProper/"

# infile1 = "a200259a.cnv"
# infile2 = "GL_PR_PF_6903703.nc"
# profile_no2 = [80] # [226,230] for Euro-Argo cruise


# indir1 = "C:/Data/ArandaVEMIVE2020/"
# indir2 = "C:/Data/ArgoData/ArgosForPlot/ArvorC/"

indir1 = "C:/Data/EuroArgoCruise/sbetulos/"
indir2 = "C:/Data/ArgoData/ArgosForPlot/RBR/"
infile1 = "a210127a.cnv"
infile2 = "GL_PR_PF_6903710.nc"
profile_no2 = [0] # [226,230] for Euro-Argo cruise
figure_handles = []

output_dir = "C:/Data/ArgoData/Figures/"
fig_dpi = 300
float_no = re.search("\d+",infile2)[0]
values1 = ['t090C','sal00'] #,'sbeox0ML/L']
values2 = ['TEMP','PSAL'] #,'DOX2_ADJUSTED']
labels = [{'t':'Temperature','x':'T (°C)', 'y':'Pressure (db)', 'col':'r', 'sty':'.r'},\
          {'t':'Salinity','x':'Practical salinity', 'y':'Pressure (db)', 'col':'b', 'sty':'.b'}]


for val1,val2,lab in zip(values1, values2,labels):
    gludge=0
    if(val2 in ['DOX2_ADJUSTED']):
        gludge = 1
    figure_handles.append(plt.figure(figsize=[5,10]))
    dat1 = apy.read_aranda_file(indir1+infile1)[0]
    
    plt.plot(dat1[val1],dat1['prDM'], 'k', alpha = 0.7)
    
    
    dat2 = xr.open_dataset(indir2+infile2)
    for i_tmp in profile_no2:
        i = i_tmp+gludge
        plt.plot(dat2[val2][i],dat2['PRES'][i],lab['col'])
        print(dat2.TIME[i].data, i)
    plt.title("WMO {},{}".format(float_no, lab['t']))
    plt.xlabel(lab['x'])
    plt.ylabel(lab['y'])
    plt.gca().invert_yaxis()
    plt.grid(True)
    filename = "{}_{}_{}".format("Comparisons",float_no,val2)
    plt.savefig(output_dir+filename+'.png' ,\
                facecolor='w',dpi=fig_dpi,bbox_inches='tight')

# Difference profile
    gludge=0
    if(val2 in ['DOX2_ADJUSTED']):
        gludge = 1

    figure_handles.append(plt.figure(figsize=[5,10]))
    dat1 = apy.read_aranda_file(indir1+infile1)[0]
    dat2 = xr.open_dataset(indir2+infile2)
    interp_dat1 = ah.interpolate_data_to_depths(dat1[val1], dat1['prDM'], dat2['PRES'][i].data)[0,:]
    for i_tmp in profile_no2:
        i = i_tmp+gludge
        plt.plot(dat2[val2][i]- interp_dat1, dat2['PRES'][i] , lab['sty'])
        plt.plot([0,0],plt.ylim(),'k', alpha = 0.5)
        print(dat2.TIME[i].data, dat1['Time'][0], i)
    plt.title("Diff, WMO {},{}".format(float_no, lab['t']))
    plt.xlabel(lab['x'])
    plt.ylabel(lab['y'])
#    plt.xlim((-0.2,0.2))
    plt.gca().invert_yaxis()
    plt.grid(True)
    filename = "{}_{}_{}".format("Diff", float_no,val2)
    plt.savefig(output_dir+filename+'.png' ,\
                facecolor='w',dpi=fig_dpi,bbox_inches='tight')
    
