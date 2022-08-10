# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 12:05:22 2022

@author: siirias
"""


import matplotlib as mp
import matplotlib.pyplot as plt
import numpy as np
from netCDF4 import Dataset
import xarray as xr
import pandas as pd
import cmocean as cmo
import arandapy as apy  # this is from: https://github.com/siirias/ArandaTools
import argohelper as ah
import re


indir1 = "C:/Data/ArgoData/ArgosForPlot/RBR_BothnianSea/"
indir2 = "C:/Data/ArgoData/ArgosForPlot/RBR_BothnianSea/"
infile1 = "GL_PR_PF_6903710.nc"
infile2 = "GL_PR_PF_6903711.nc"

output_dir = "C:/Data/ArgoData/Figures/"

values1 = ['t090C','sal00'] #,'sbeox0ML/L']
values2 = ['TEMP','PSAL'] #,'DOX2_ADJUSTED']
labels = [{'t':'Temperature','x':'T (°C)', 'y':'Pressure (db)', 'col':'r', 'sty':'.r'},\
          {'t':'Salinity','x':'Practical salinity', 'y':'Pressure (db)', 'col':'b', 'sty':'.b'}]

t_s_coeff = 3.0  # how many km one day corresponds, while calculating distance
#first read the main profile:

def plot_comparison_for_profiles(filename1, filename2, prof_no1, prof_no2,\
                ax = None, parameters = ['TEMP', 'PSAL'], figsize = [5,5]):
    for param in parameters:
        prof1 = xr.open_dataset(filename1)
        prof2 = xr.open_dataset(filename2)
        if(not ax):
            plt.figure(figsize = figsize)
        else:
            plt.sca(ax)
        plt.plot(prof1[param][prof_no1],prof1['PRES'][prof_no1])
        plt.plot(prof2[param][prof_no2],prof2['PRES'][prof_no2])
        plt.gca().invert_yaxis()
        plt.grid(True)
        # filename = "{}_{}_{}".format("Comparisons",float_no,val2)
        # plt.savefig(output_dir+filename+'.png' ,\
        #             facecolor='w',dpi=fig_dpi,bbox_inches='tight')    
        return plt.gca()

if(True):
    main_profile = xr.open_dataset(indir1+infile1)
    compare_profile = xr.open_dataset(indir2+infile2)
    distances_t = np.zeros((len(main_profile['TIME']), len(compare_profile['TIME'])))
    distances_s = distances_t.copy()
    for i in range(distances_t.shape[0]): #the original profiles
        for j in range(distances_t.shape[1]): #the comparison profiles
            #distances in kilometers
            distances_s[i,j] = ah.distance( 
                                (main_profile['LATITUDE'][i], 
                                main_profile['LONGITUDE'][i]), 
                                (compare_profile['LATITUDE'][j],
                                 compare_profile['LONGITUDE'][j])
                                ) 
            #distances in time  (microseconds initially)
            distances_t[i,j] = compare_profile['TIME'][j] - main_profile['TIME'][i]
            distances_t[i,j] = distances_t[i,j] / (1000000000.0*60.0*60.0*24.0) #days
    distances_tot = distances_s + t_s_coeff*np.abs(distances_t)
    closest_ones = np.argmin(distances_tot, axis = 1)
    closest_values = np.min(distances_tot, axis = 1)
    #let's gather the best pairs,
    #First gather both indices, and the distance in a table
    matches = np.concatenate((
        np.reshape(range(len(main_profile['TIME'])),(-1,1)),
        np.reshape(closest_ones,(-1,1)), 
        np.reshape(closest_values,(-1,1)),
        ), 
        axis=1)
    #Then sort the table by distances:
    matches = matches[matches[:,2].argsort(),:]
    #And last, drop pairs, where
    # a closer original profile has been found
    tmp_pair = matches[:,1].astype(int)
    tmp_bool = tmp_pair.copy().astype(bool)
    for i in range(len(tmp_bool)):
        if(tmp_pair[i] in tmp_pair[:i]):
            tmp_bool[i] = False
        else:
            tmp_bool[i] = True
    best_matches = matches[tmp_bool,:]
    
    
fig, axes = plt.subplots(nrows = 1, ncols =3, figsize = (15,5))
fn = 0
plot1 = axes[fn].imshow(distances_t, cmap = cmo.cm.diff)
axes[fn].invert_yaxis()
axes[fn].title.set_text('time distance')
plt.colorbar(plot1, ax = axes[fn])

fn = 1
plot2 = axes[fn].imshow(distances_s, cmap = cmo.cm.matter)
axes[fn].invert_yaxis()
axes[fn].title.set_text('space distance')
plt.colorbar(plot2, ax = axes[fn])

fn = 2
plot3 = axes[fn].imshow(distances_tot, cmap = cmo.cm.algae)
axes[fn].plot(closest_ones, range(len(closest_ones)))
axes[fn].invert_yaxis()
axes[fn].title.set_text('total distance')
plt.colorbar(plot3, ax = axes[fn])

fig2 = plt.figure(figsize = (15,5))
plt.plot(main_profile['TIME'], closest_values)

comparisons_to_plot = min(5,best_matches.shape[0])
fig3, axes3 = plt.subplots(nrows = 1, ncols =comparisons_to_plot, 
                           figsize = (3*comparisons_to_plot,5))
curr_pair = 0
for ax in axes3:
    plot_comparison_for_profiles(indir1 + infile1,
                                 indir2 + infile2,
                                 int(np.round(best_matches[curr_pair,1])),
                                 int(np.round(best_matches[curr_pair,0])),
                                 ax = ax)
    curr_pair += 1



if(False):
    for i_tmp in range(len(main_profile['TIME'])):
        i = i_tmp
        plt.plot(main_profile['TEMP'][i],main_profile['PRES'][i],'r', alpha = 0.3)
        print(main_profile.TIME[i].data, i)
    
    for i_tmp in range(len(compare_profile['TIME'])):
        i = i_tmp
        plt.plot(compare_profile['TEMP'][i],compare_profile['PRES'][i],'b', alpha = 0.3)
        print(compare_profile.TIME[i].data, i)
