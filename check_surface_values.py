# -*- coding: utf-8 -*-
"""
This script contains functions to check the surface salinity of profiles 
in Argo floats, and generate plots to visualize the results. 

The main purpose of this script is to determine if the surface salinity 
of profiles in a given float exceeds pre-defined limits for the area 
in which the float is located.

Functions:
- get_profiles(WMO): 
    Returns a list of profiles for a given Argo float (identified by its WMO ID).
- determine_dmqc_area(lat, lon): 
    Determines the DMQC area for a given set of latitude and longitude coordinates.
- plot_salinity_profiles(dataset, result, WMO): 
    Generates a plot of salinity profiles for a given Argo float, 
    with faulty profiles highlighted.
- check_surface_salinity(dataset, WMO = None): 
    Checks the surface salinity of profiles in a given dataset, 
    and returns a list of profile numbers and a boolean indicating 
    whether the profile is faulty.
- check_surface_salinity_for_float(WMO, output_dir): 
    Calls the above functions to check the surface salinity of profiles 
    for a given Argo float, and generates a plot of the results.

Created on Sun Apr 16 15:46:43 2023

@author: siirias
"""

import matplotlib.pyplot as plt
import numpy as np

import argopy
from argopy import DataFetcher as ArgoDataFetcher


# Define a function to retrieve profiles for a given WMO number using argopy
def get_profiles(WMO):
    ArgoSet = ArgoDataFetcher().float([WMO])
    ds = ArgoSet.to_xarray()
    dsg = ds.groupby("CYCLE_NUMBER")
    return dsg

# Define a function to plot salinity profiles
def plot_salinity_profiles(dataset, result = None, WMO = None):
    colors = []
    if not result:
        result = [True]*len(dataset)
    for c in result:
            if c[1]:
                colors.append((0.0,0.4,0.4,0.1))  # Greenish color with transparency
            else:
                colors.append((1.0,0.3,0.3,0.8))  # Reddish color with transparency
            
    plt.figure()
    for data, color in zip(dataset,colors):
        prof_no = data[0]
        data = data[1]
        plt.plot(data.PSAL, data.PRES, color = color)
    plt.gca().invert_yaxis()
    plt.title(f"WMO {WMO}")
    plt.ylabel('Pressure (dB)')
    plt.xlabel("Salinity (PSU)")
    plt.grid()

# Define a function to determine which DMQC area a 
# given latitude and longitude coordinates belong to
def determine_dmqc_area(lat, lon):
    areas = [
        {'name': 'BothnianBay',
         'lat_min': 63.5, 'lat_max':66.0,
         'lon_min': 15.0, 'lon_max': 26.0},
        {'name': 'BothnianSea',
         'lat_min': 59.8, 'lat_max': 63.5,
         'lon_min': 15.0, 'lon_max': 26.0},
        {'name': 'NorthernBalticProper',
         'lat_min': 58.4, 'lat_max': 59.8,
         'lon_min': 15.0, 'lon_max': 24.0},
        {'name': 'BalticProper',
         'lat_min': 56.2, 'lat_max': 58.4,
         'lon_min': 15.0, 'lon_max': 26.0},
        {'name': 'BornholmBasin',
         'lat_min': 53.5, 'lat_max': 56.2,
         'lon_min': 14.4, 'lon_max': 17.5},
        {'name': 'GdanskBasin',
         'lat_min': 53.5, 'lat_max': 56.2,
         'lon_min': 17.5, 'lon_max': 22.0}
        ]
    area = 'Unknown'
    for a in areas:
        if lat>=a['lat_min'] and lat<=a['lat_max'] and\
           lon>=a['lon_min'] and lon<=a['lon_max'] :
               area = a['name']
               return area

# A function to check surface salinity of Argo floats against 
# predefined limits for each area in the Baltic Sea

def check_surface_salinity(dataset, WMO = None):
    sss_limits ={'BothnianBay':4.0,
                 'BothnianSea':6.0,
                 'NorthernBalticProper':7.5,
                 'BalticProper':8.0,
                 'BornholmBasin':8.5,
                 'GdanskBasin':8.2}
    max_pres = 10.0  #deciBar
    result = []
    for num, profile in dataset:
        lat = float(profile.LATITUDE[0])
        lon = float(profile.LONGITUDE[0])
        area = determine_dmqc_area(lat, lon)
        sss_lim = sss_limits[area]
        # Check if the surface salinity exceeds the limit 
        # for pressures less than the maximum pressure determined as surface
        if(np.sum(profile.where(profile.PRES<max_pres).PSAL.values>sss_lim)>0):
            result.append((num,False))
            print(f"Profile {num} in WMO: {WMO} Faulty")
        else:
            result.append((num,True))
    return result

def check_surface_salinity_for_float(WMO, output_dir = "C:/Data/ArgoData/Figures/DMQC/"):
    dataset = get_profiles(WMO)
    result = check_surface_salinity(dataset,WMO)
    plot_salinity_profiles(dataset,result,WMO)
    filename = f"Surface_test_WMO{WMO}"
    plt.savefig(output_dir+filename+'.png' ,\
                facecolor='w', dpi=300, bbox_inches='tight')
#Example usage
WMOS = [6903697, 6902020, 6902024]
for WMO in WMOS:
    check_surface_salinity_for_float(WMO)