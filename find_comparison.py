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
from itertools import cycle  # used to cycle a short list to match others
import time #Just to track how long the script runs

#Bothnian Sea
# indir1 = "C:/Data/ArgoData/ArgosForPlot/RBR_BothnianSea/"
# indir2 = "C:/Data/ArgoData/ArgosForPlot/RBR_BothnianSea/"
# infile1 = "GL_PR_PF_6903710.nc"
# infiles2 = ["GL_PR_PF_6903711.nc"]

# #Gotland Deep
indir1 = "C:/Data/ArgoData/ArgosForPlot/RBR_BalticProper/"
indir2 = "C:/Data/ArgoData/ArgosForPlot/RBR_BalticProper/"
infile1 = "GL_PR_PF_6903709.nc"
#infile2 = "GL_PR_PF_6903708.nc"
infiles2 = ["GL_PR_PF_6903706.nc", "GL_PR_PF_6903708.nc"]


params_list = [['TEMP'], ['PSAL']]

output_dir = "C:/Data/ArgoData/Figures/"

plot_range = 20.0 
t_s_coeff = 3.0  # how many km one day corresponds, while calculating distance
example_number = 8
start_time = time.time()

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

best_matches = {} 
for infile2 in infiles2:
    #first read the main profile:
    WMO_original = re.search('[^\d](\d*)\.nc',infile1).groups()[0]
    #then teh comparison profile
    WMO_compare = re.search('[^\d](\d*)\.nc',infile2).groups()[0]
    
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
    matches = pd.DataFrame(zip(
                range(len(main_profile['TIME'])),
                closest_ones, 
                closest_values,
                [distances_t[x,closest_ones[x]] for x in range(len(closest_ones))],
                [distances_s[x,closest_ones[x]] for x in range(len(closest_ones))],
                cycle([infile2]),
                cycle([WMO_original]),
                cycle([WMO_compare])
                ),
                columns = ["prof1", "prof2", "distance", "distance_t", 
                           "distance_s", "infile2", "WMO1", "WMO2"])
    matches = matches.sort_values("distance")
    matches = matches.reset_index(drop = True)

    ok_profs = []
    for i in range(matches.shape[0]):
        if(matches.iloc[i]['prof2'] in ok_profs): # there is already a better match
           matches.iloc[i,matches.columns.get_loc('distance')] = -1.0
        else:
            ok_profs.append(matches.iloc[i]['prof2'])

    #best_matches[WMO_compare] = matches[tmp_bool,:]
    best_matches[WMO_compare] = matches[matches['distance']> 0.0] #drop the dublicates detected        
    best_matches[WMO_compare] = best_matches[WMO_compare].reset_index(drop = True)    
    fig, axes = plt.subplots(nrows = 1, ncols =3, figsize = (15,5))
    fn = 0
    plot1 = axes[fn].imshow(distances_t, cmap = cmo.cm.diff)
    axes[fn].invert_yaxis()
    axes[fn].title.set_text('time distance (hours)')
    axes[fn].set_ylabel('Profile index (WMO {})'.format(WMO_original))
    axes[fn].set_xlabel('Profile index (WMO {})'.format(WMO_compare))
    plt.colorbar(plot1, ax = axes[fn])
    plot1.set_clim(-plot_range/t_s_coeff,plot_range/t_s_coeff)  #distance accepted in hours
    fn = 1
    plot2 = axes[fn].imshow(distances_s, cmap = cmo.cm.matter)
    axes[fn].invert_yaxis()
    axes[fn].title.set_text('space distance (km)')
    axes[fn].set_ylabel('Profile index (WMO {})'.format(WMO_original))
    axes[fn].set_xlabel('Profile index (WMO {})'.format(WMO_compare))
    plt.colorbar(plot2, ax = axes[fn])
    plot2.set_clim(0.0,plot_range)  #distance accepted in km
    
    fn = 2
    plot3 = axes[fn].imshow(distances_tot, cmap = cmo.cm.algae)
    axes[fn].plot(closest_ones, range(len(closest_ones)))
    axes[fn].invert_yaxis()
    axes[fn].title.set_text('total distance')
    axes[fn].set_ylabel('Profile index (WMO {})'.format(WMO_original))
    axes[fn].set_xlabel('Profile index (WMO {})'.format(WMO_compare))
    plt.colorbar(plot3, ax = axes[fn])
    plot3.set_clim(0.0,plot_range*10.0)  #distance accepted in km
    
    # fig2 = plt.figure(figsize = (15,5))
    # plt.plot(main_profile['TIME'], closest_values)

best_of_best = pd.concat(best_matches).sort_values('distance')
best_of_best = best_of_best.reset_index(drop = True)

    
for values, col1, col2 in zip(params_list,['r','c'],['b','m']):
    comparisons_to_plot = min(example_number,best_of_best.shape[0])
    fig3, axes3 = plt.subplots(nrows = 1, ncols =comparisons_to_plot, 
                               figsize = (3*comparisons_to_plot,5))
    curr_pair = 0
    for ax in axes3:
        the_match = best_of_best.iloc[curr_pair]
        plot_comparison_for_profiles(indir1 + infile1,
                                     indir2 + the_match['infile2'],
                                     the_match['prof1'],
                                     the_match['prof2'],
                                     ax = ax,
                                     parameters = values,
                                     col1 = col1,
                                     col2 = col2)
        ax.title.set_text('{} vs {}-{}\n{}\n{:.1f} km, {:.1f} d'.format(
            the_match['prof1'],the_match['prof2'], the_match['WMO2'],
            values[0],
            the_match['distance_s'],
            the_match['distance_t'],
            ))
        curr_pair += 1
            

completion_time = time.time() - start_time