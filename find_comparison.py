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

# #Bothnian Sea
# indir1 = "C:/Data/ArgoData/ArgosForPlot/RBR_BothnianSea/"
# indir2 = "C:/Data/ArgoData/ArgosForPlot/RBR_BothnianSea/"
# infile1 = "GL_PR_PF_6903710.nc"
# infiles2 = ["GL_PR_PF_6903711.nc", "GL_PR_PF_6903698.nc", "GL_PR_PF_6903699.nc",
#             "GL_PR_PF_6903702.nc"]

#Gotland Deep
indir1 = "C:/Data/ArgoData/ArgosForPlot/RBR_BalticProper/"
indir2 = "C:/Data/ArgoData/ArgosForPlot/RBR_BalticProper/"
infile1 = "GL_PR_PF_6903706.nc"
#infile2 = "GL_PR_PF_6903708.nc"
infiles2 = ["GL_PR_PF_6903708.nc", "GL_PR_PF_3902110.nc",
            "GL_PR_PF_6904116.nc", "GL_PR_PF_6904117.nc", 
            "GL_PR_PF_6904226.nc", "GL_PR_PF_7900587.nc"]


#params_list = [['TEMP'], ['PSAL']]
params_list = [ ['PSAL']]

output_dir = "C:/Data/ArgoData/Figures/D412/"
csv_decimal = ','
csv_sep = ';'
fig_dpi = 300
plot_range = 20.0 
t_s_coeff = 3.0  # how many km one day corresponds, while calculating distance
example_number = 8
start_time = time.time()

def get_prof_param(filename, param, prof_no):
    #small tool to get from file one parameter of given profile number
    with xr.open_dataset(filename) as f:
        if param not in f:
            if(param == 'TEMP'):
                param = 'TEMP_ADJUSTED'
            elif(param == 'TEMP_ADJUSTED'):
                param = 'TEMP'
            elif(param == 'PSAL'):
                param = 'PSAL_ADJUSTED'
            elif(param == 'PSAL_ADJUSTED'):
                param = 'PSAL'
        return f[param][prof_no]
        
def calculate_match_fitness(filename1, filename2, 
                            prof_no1, prof_no2, 
                            parameters = ['TEMP', 'PSAL'] ):
    fitnesses = {}
    for param in parameters:
        prof1_x = get_prof_param(filename1, param, prof_no1)
        prof1_y = get_prof_param(filename1, 'PRES', prof_no1)
        prof2_x = get_prof_param(filename2, param, prof_no2)
        prof2_y = get_prof_param(filename2, 'PRES', prof_no2)
        print(param)
        print((prof1_y, prof1_x, prof2_y, prof2_x))
        fitnesses[param] = ah.compare_profiles(prof1_y, prof1_x, prof2_y, prof2_x)
    return fitnesses

def calculate_match_bias(filename1, filename2, 
                            prof_no1, prof_no2, 
                            parameters = ['TEMP', 'PSAL'] ):
    biases = {}
    for param in parameters:
        prof1_x = get_prof_param(filename1, param, prof_no1)
        prof1_y = get_prof_param(filename1, 'PRES', prof_no1)
        prof2_x = get_prof_param(filename2, param, prof_no2)
        prof2_y = get_prof_param(filename2, 'PRES', prof_no2)
        biases[param] = ah.profile_bias(prof1_y, prof1_x, prof2_y, prof2_x)
    return biases
    

def plot_comparison_for_profiles(filename1, filename2, prof_no1, prof_no2,\
                ax = None, parameters = ['TEMP'], figsize = [5,5], \
                col1 = 'r', col2 = 'b'):
    for param in parameters:
        prof1_x = get_prof_param(filename1, param, prof_no1)
        prof1_y = get_prof_param(filename1, 'PRES', prof_no1)
        prof2_x = get_prof_param(filename2, param, prof_no2)
        prof2_y = get_prof_param(filename2, 'PRES', prof_no2)
        if(not ax):
            plt.figure(figsize = figsize)
        else:
            plt.sca(ax)
        plt.plot(prof1_x,prof1_y, col1)
        plt.plot(prof2_x,prof2_y, col2)
        plt.gca().invert_yaxis()
        plt.grid(True)
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
                           "distance_s", "infile2", "WMO1", "WMO2",
                           
                          ])
    
    
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
    filename = "{}_{}_{}".format("distance",WMO_original, WMO_compare)
    plt.savefig(output_dir+filename+'.png' ,\
                facecolor='w',dpi=fig_dpi,bbox_inches='tight')
    print("saved {}".format(output_dir+filename+'.png'))

best_of_best = pd.concat(best_matches).sort_values('distance')
best_of_best = best_of_best.reset_index(drop = True)

    
for values, col1, col2 in zip(params_list,['r','m'],['b','c']):
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
        ax.title.set_text('{} vs {} WMO{}\n{}\n{:.1f} km, {:.1f} d'.format(
            the_match['prof1'],the_match['prof2'], the_match['WMO2'],
            values[0],
            the_match['distance_s'],
            np.abs(the_match['distance_t']),
            ))
        curr_pair += 1
    filename = "{}_{}_{}".format("profiles",the_match['WMO1'], values[0])
    plt.savefig(output_dir+filename+'.png' ,\
                facecolor='w',dpi=fig_dpi,bbox_inches='tight')
    print("saved {}".format(output_dir+filename+'.png'))

    fitnesses = []
    biases = []
    for i in range(15):        
            the_match = best_of_best.iloc[i]
            fitnesses.append(calculate_match_fitness(
                                        indir1 + infile1,
                                        indir2 + the_match['infile2'],
                                        the_match['prof1'],
                                        the_match['prof2']))
            biases.append(calculate_match_bias(
                                        indir1 + infile1,
                                        indir2 + the_match['infile2'],
                                        the_match['prof1'],
                                        the_match['prof2']))



analysis_data = pd.DataFrame(zip([x['TEMP'] for x in biases], [x['PSAL'] for x in biases],
                                 [x['TEMP'] for x in fitnesses], [x['PSAL'] for x in fitnesses]), 
                             columns=['TEMP_bias','PSAL_bias', 'TEMP_fit','PSAL_fit'])

print(analysis_data)
print(analysis_data.describe())
filename = "{}_{}.csv".format("analysis",the_match['WMO1'])
tmp = analysis_data.to_csv(
    line_terminator='\n', decimal = csv_decimal, sep = csv_sep)
tmp += analysis_data.describe().to_csv(
    line_terminator='\n', decimal = csv_decimal, sep = csv_sep)
open(output_dir+filename,'w').writelines(tmp)
