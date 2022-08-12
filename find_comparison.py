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


#Bothnian Sea
indir1 = "C:/Data/ArgoData/ArgosForPlot/RBR_BothnianSea/"
indir2 = "C:/Data/ArgoData/ArgosForPlot/RBR_BothnianSea/"
infile1 = "GL_PR_PF_6903710.nc"
infile2 = "GL_PR_PF_6903711.nc"

# #Gotland Deep
# indir1 = "C:/Data/ArgoData/ArgosForPlot/RBR_BalticProper/"
# indir2 = "C:/Data/ArgoData/ArgosForPlot/RBR_BalticProper/"
# infile1 = "GL_PR_PF_6903709.nc"
# #infile2 = "GL_PR_PF_6903708.nc"
# infile2 = "GL_PR_PF_6903706.nc"
#params_list = [['TEMP_ADJUSTED'], ['PSAL_ADJUSTED']]

params_list = [['TEMP'], ['PSAL']]

output_dir = "C:/Data/ArgoData/Figures/"



t_s_coeff = 3.0  # how many km one day corresponds, while calculating distance
#first read the main profile:
WMO_original = re.search('[^\d](\d*)\.nc',infile1).groups()[0]
WMO_compare = re.search('[^\d](\d*)\.nc',infile2).groups()[0]
def plot_comparison_for_profiles(filename1, filename2, prof_no1, prof_no2,\
                ax = None, parameters = ['TEMP'], figsize = [5,5], \
                col1 = 'r', col2 = 'b'):
    for param in parameters:
        prof1 = xr.open_dataset(filename1)
        prof2 = xr.open_dataset(filename2)
        param1 = param
        param2 = param
        if param1 not in prof1:
            if(param1 == 'TEMP'):
                param1 = 'TEMP_ADJUSTED'
            elif(param1 == 'TEMP_ADJUSTED'):
                param1 = 'TEMP'
            elif(param1 == 'PSAL'):
                param1 = 'PSAL_ADJUSTED'
            elif(param1 == 'PSAL_ADJUSTED'):
                param1 = 'PSAL'
        if param2 not in prof2:
            if(param2 == 'TEMP'):
                param2 = 'TEMP_ADJUSTED'
            elif(param2 == 'TEMP_ADJUSTED'):
                param2 = 'TEMP'
            elif(param2 == 'PSAL'):
                param2 = 'PSAL_ADJUSTED'
            elif(param2 == 'PSAL_ADJUSTED'):
                param2 = 'PSAL'

        if(not ax):
            plt.figure(figsize = figsize)
        else:
            plt.sca(ax)
        plt.plot(prof1[param1][prof_no1],prof1['PRES'][prof_no1], col1)
        plt.plot(prof2[param2][prof_no2],prof2['PRES'][prof_no2], col2)
        plt.gca().invert_yaxis()
        plt.grid(True)
        # filename = "{}_{}_{}".format("Comparisons",float_no,val2)
        # plt.savefig(output_dir+filename+'.png' ,\
        #             facecolor='w',dpi=fig_dpi,bbox_inches='tight')    
        return plt.gca()

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
axes[fn].set_ylabel('Profile index (WMO {})'.format(WMO_original))
axes[fn].set_xlabel('Profile index (WMO {})'.format(WMO_compare))
plt.colorbar(plot1, ax = axes[fn])

fn = 1
plot2 = axes[fn].imshow(distances_s, cmap = cmo.cm.matter)
axes[fn].invert_yaxis()
axes[fn].title.set_text('space distance')
axes[fn].set_ylabel('Profile index (WMO {})'.format(WMO_original))
axes[fn].set_xlabel('Profile index (WMO {})'.format(WMO_compare))
plt.colorbar(plot2, ax = axes[fn])

fn = 2
plot3 = axes[fn].imshow(distances_tot, cmap = cmo.cm.algae)
axes[fn].plot(closest_ones, range(len(closest_ones)))
axes[fn].invert_yaxis()
axes[fn].title.set_text('total distance')
axes[fn].set_ylabel('Profile index (WMO {})'.format(WMO_original))
axes[fn].set_xlabel('Profile index (WMO {})'.format(WMO_compare))
plt.colorbar(plot3, ax = axes[fn])

fig2 = plt.figure(figsize = (15,5))
plt.plot(main_profile['TIME'], closest_values)

for values, col1, col2 in zip(params_list,['r','c'],['b','m']):
    comparisons_to_plot = min(15,best_matches.shape[0])
    fig3, axes3 = plt.subplots(nrows = 1, ncols =comparisons_to_plot, 
                               figsize = (3*comparisons_to_plot,5))
    curr_pair = 0
    for ax in axes3:
        plot_comparison_for_profiles(indir1 + infile1,
                                     indir2 + infile2,
                                     int(np.round(best_matches[curr_pair,0])),
                                     int(np.round(best_matches[curr_pair,1])),
                                     ax = ax,
                                     parameters = values,
                                     col1 = col1,
                                     col2 = col2)
        ax.title.set_text('{}\n{:.1f} km, {:.1f} h'.format(
            values[0],
            distances_s[int(best_matches[curr_pair,0]), int(best_matches[curr_pair,1])],
            distances_t[int(best_matches[curr_pair,0]), int(best_matches[curr_pair,1])]*24,
            ))
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
